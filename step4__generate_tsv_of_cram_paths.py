import argparse
import logging
import hail as hl
import os
import pandas as pd
import re

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("get_tsv_for_haplotype_caller")
logger.setLevel(logging.INFO)

hl.init(log="/dev/null", idempotent=True)

#args = argparse.Namespace()
#args.gnomad_metadata_tsv = "gnomad.exomes.v4.0.metadata.tsv.gz"
#args.gnomad_readviz_tsvs_directory = "gs://gnomad-readviz/v4.0/readviz_tsvs"
#args.output_bucket = "gs://gnomad-readviz/v4.0/readviz_tsvs"

#%%
def main(args):
    logger.info(f"Listing {args.gnomad_readviz_tsvs_directory}")

    tsv_paths = [x["path"] for x in hl.hadoop_ls(args.gnomad_readviz_tsvs_directory)]
    tsv_paths = {re.sub(".tsv.bgz$", "", os.path.basename(path)): path for path in tsv_paths}

    #%%
    logger.info(f"Found {len(tsv_paths):,d} .tsv.bgz files in {args.gnomad_readviz_tsvs_directory}")

    #%%
    df = pd.read_table(args.gnomad_metadata_tsv)
    df["s_join"] = df["s"].str.replace("/", "_")

    mismatched_sample_ids = set(tsv_paths.keys()) - set(df.s_join)
    if mismatched_sample_ids:
        raise ValueError(f"Found {len(mismatched_sample_ids)} .tsv.bgz files in {args.gnomad_readviz_tsvs_directory} "
                         f"with filenames that mismatch the sample ids in the gnomAD metadata from "
                         f"{args.gnomad_metadata_tsv}: {mismatched_sample_ids}")

    df["variants_tsv_bgz"] = df["s_join"].map(tsv_paths)
    df = df[~df["variants_tsv_bgz"].isna()]
    df = df.rename(columns={
        "s": "sample_id",
        "cram_path": "cram",
        "crai_path": "crai",
    })

    df = df[["sample_id", "cram", "crai", "variants_tsv_bgz", "gatk_version"]]

    #%%

    logger.info(f"Writing {len(df):,d} rows to {args.output_tsv}")
    with hl.hadoop_open(args.output_tsv, "w") as f:
        df.to_csv(f, header=True, sep="\t", index=False)

    logger.info("Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gnomad-readviz-tsvs-directory",
        help="Cloud storage directory where step3 wrote all per-sample TSVs",
        default="gs://gnomad-readviz/v4.0/eadviz_tsvs",
    )
    parser.add_argument(
        "--gnomad-metadata-tsv",
        help="Path of local gnomAD metadata table with columns: 's', 'cram_path', 'crai_path', and 'gatk version'",
        default="gnomad.exomes.v4.0.metadata.tsv.gz",
    )
    parser.add_argument(
        "--output-tsv",
        help="Path where to write the output tsv",
        default="gs://gnomad-readviz/v4.0/step4_output_cram_and_tsv_paths.tsv.gz",
    )
    args = parser.parse_args()

    main(args)
