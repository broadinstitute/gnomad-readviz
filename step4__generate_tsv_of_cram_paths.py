import argparse
import logging

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import call_parallel_file_exists

from .utils import get_sample_ids

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("get_tsv_for_haplotype_caller")
logger.setLevel(logging.INFO)


def main(args):
    """
    Generate single TSV with sample IDs, cram paths, and variant TSV file paths.

    NOTE: Run locally (`parallel_file_exists` only works locally).
    """
    output_bucket = args.output_bucket

    logger.info("Getting sample IDs...")
    sample_ids = get_sample_ids(args.ids_file)

    tsvs = [f"{output_bucket}/{sample}.tsv.bgz" for sample in sample_ids]
    tsv_files_exist = call_parallel_file_exists(tsvs)

    logger.info("Starting cram existence checks...")
    cram_map = {}
    with hl.hadoop_open(args.cram_map) as c:
        for line in c:
            sample, cram = line.strip().split("\t")
            cram_map[sample] = cram
    cram_files_exist = call_parallel_file_exists(list(cram_map.values()))

    logger.info("Starting to write to output TSV...")
    with hl.hadoop_open(args.cram_map) as s, hl.hadoop_open(
        f"{output_bucket}/inputs/step4_output_cram_and_tsv_paths_table.tsv", "w",
    ) as o:
        o.write("sample_id\tcram\tcrai\tvariants_tsv_bgz\n")

        for line in s:
            sample = line.strip()
            if sample not in cram_map:
                raise DataException(
                    f"{sample} is missing a cram path. Please double check and restart!"
                )

            cram = cram_map[sample]
            tsv = f"{output_bucket}/{sample}.tsv.bgz"
            if not cram_files_exist[cram]:
                raise DataException(
                    f"{sample}'s cram does not exist. Please double check and restart!"
                )
            if not tsv_files_exist[tsv]:
                raise DataException(
                    f"{sample} is missing their variants TSV file. Please double check and restart!"
                )
            o.write(f"{sample}\t{cram}\t{cram}.crai\t{tsv}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ids-file", help="Path to local file with sample IDs.")
    parser.add_argument(
        "--cram-map",
        help="Path to file in cloud storage with cram to sample mapping. File should have two columns: sample ID and cram path.",
    )
    parser.add_argument(
        "--output-bucket", help="Path to output bucket (where to store output TSVs)."
    )
    args = parser.parse_args()

    main(args)
