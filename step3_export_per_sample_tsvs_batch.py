import asyncio
import argparse
import logging
import os
from typing import Dict, List

import hailtop.batch as hb
import hail as hl

from gnomad.utils.file_utils import parallel_file_exists

from .utils import get_sample_ids

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("run_per_sample_tsv_export")
logger.setLevel(logging.INFO)


def export_tsv(ht_path: str, sample_id: str, tsv_path: str, success_path: str) -> None:
    """
    Read in hail Table, filter to specified sample ID, and export TSV.

    Also write an empty success file to ensure no Batch jobs finish with only partially exported TSV.

    :param str ht_path: Path to input hail Table.
    :param str sample_id: Sample for which to export a TSV.
    :param str tsv_path: Path to output TSV.
    :param str success_path: Path to output success file.
    :return: None
    """
    ht = hl.read_table(ht_path)
    ht = ht.filter(ht.S == sample_id)
    ht = ht.naive_coalesce(1)
    ht.export(tsv_path)

    with hl.hadoop_open(success_path, "w") as o:
        o.write("")


def main(args):

    backend = hb.ServiceBackend(
        billing_project=args.billing_project,
        bucket=args.tmp_bucket,
        google_project=args.google_project,
    )
    b = hb.Batch(
        backend=backend, default_cpu=1, default_python_image=args.python_image,
    )

    logger.info("Extracting sample IDs...")
    sample_ids = get_sample_ids(args.ids_file, args.header)

    logger.info("Preparing to start batch job...")
    output_bucket = args.output_bucket
    files = [f"{output_bucket}/{sample}_success.txt" for sample in sample_ids]
    file_exists = asyncio.get_event_loop().run_until_complete(
        parallel_file_exists(files)
    )

    for sample in sample_ids:
        logger.info("Working on %s", sample)
        if file_exists[f"{output_bucket}/{sample}_success.txt"]:
            logger.info(
                "Output success txt file already exists, skipping %s...", sample
            )
            continue
        j = b.new_python_job(name=sample)
        j.call(
            export_tsv,
            args.ht_path,
            sample,
            f"{output_bucket}/{sample}.tsv.bgz",
            f"{output_bucket}/{sample}_success.txt",
        )

    b.run(wait=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--billing-project", help="Billing project to use with hail Batch."
    )
    parser.add_argument(
        "--tmp-bucket",
        help="Bucket to use for batch temporary data. Use name of bucket only (do not specify 'gs://' prefix.)",
    )
    parser.add_argument(
        "--google-project",
        help="Google project to use with hail Batch when authenticating storage.",
    )
    parser.add_argument(
        "--python-image",
        help="Docker image to use. Must have python and dill installed.",
        default="gcr.io/broad-mpg-gnomad/tgg-methods-vm:20210623",
    )
    parser.add_argument("--ids-file", help="File with sample IDs")
    parser.add_argument(
        "--header",
        help="Whether file with sample IDs has a header line",
        action="store_true",
    )
    parser.add_argument("--output-bucket", help="Full path to output bucket.")
    parser.add_argument(
        "--ht-path",
        help="Full path to hail Table with samples to be extracted for readviz. Should be keyed by sample ID.",
    )
    args = parser.parse_args()

    main(args)
