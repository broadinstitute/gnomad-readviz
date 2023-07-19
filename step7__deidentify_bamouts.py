from datetime import datetime
import hail as hl   # used for hadoop file utils
import logging
import math
import os
import pandas as pd
import subprocess
import sys
import tqdm

# TODO create & keep a mapping of variant to sample ids?

from batch import batch_utils

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


GCLOUD_PROJECT = "broad-mpg-gnomad"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"

DOCKER_IMAGE = "weisburd/gnomad-readviz@sha256:555a77391da1ce4b7a77615f830e9a566d7c3d018902b5b8af2f50ecf071f1c7"

OUTPUT_BUCKET = "gs://gnomad-bw2/gnomad_all_readviz_bamout_deidentified_v3_and_v31_fixed__20210101"


def parse_args():
    """Parse command line args."""

    p = batch_utils.init_arg_parser(default_cpu=0.25, default_billing_project="gnomAD-readviz", gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    p.add_argument("-n", "--num-samples-to-process", help="For testing, process only the given sample id(s).", type=int)
    p.add_argument("--random", action="store_true", help="Select random sample")
    p.add_argument("-s", "--sample-to-process", help="For testing, process only the given sample id(s).", action="append")
    p.add_argument("--output-dir", help="Where to write combined bams.", default=OUTPUT_BUCKET)
    p.add_argument("cram_and_tsv_paths_table", help="A text file containing at least these columns: sample_id, cram_path")
    args = p.parse_args()

    return p, args


def main():

    p, args = parse_args()

    df = pd.read_table(args.cram_and_tsv_paths_table)
    if {"sample_id", "output_bamout_bam", "output_bamout_bai", "variants_tsv_bgz"} - set(df.columns):
        p.error(f"{args.tsv_path} must contain 'sample_id', 'output_bamout_bam', 'variants_tsv_bgz' columns")

    if args.num_samples_to_process:
        if args.random:
            df = df.sample(n=args.num_samples_to_process)
        else:
            df = df.iloc[:args.num_samples_to_process]

    if args.sample_to_process:
        df = df[df.sample_id.isin(set(args.sample_to_process))]

    logging.info(f"Processing {len(df)} samples")

    # check that all buckets are in "US-CENTRAL1" or are multi-regional to avoid egress charges to the Batch cluster
    batch_utils.set_gcloud_project(GCLOUD_PROJECT)
    with open("deidentify_bamout.py", "rt") as f:
        deidentify_bamouts_script = f.read()

    # process sample(s)
    if not args.sample_to_process and not args.num_samples_to_process:
        # if processing entire table, listing all files up front ends up being faster
        existing_deidentify_output_bams = subprocess.check_output(f"gsutil -m ls {args.output_dir}/*.deidentify_output.bam", shell=True, encoding="UTF-8").strip().split("\n")
        existing_deidentify_output_sorted_bams = subprocess.check_output(f"gsutil -m ls {args.output_dir}/*.deidentify_output.sorted.bam", shell=True, encoding="UTF-8").strip().split("\n")

    hl.init(log="/dev/null")
    with batch_utils.run_batch(args, batch_name=f"deidentify bamouts: {len(df)} samples") as batch:
        for _, row in tqdm.tqdm(df.iterrows(), unit=" samples"):
            output_bam_path = os.path.join(args.output_dir, f"{row.sample_id}.deidentify_output.bam")
            output_sorted_bam_path = os.path.join(args.output_dir, f"{row.sample_id}.deidentify_output.sorted.bam")

            if args.sample_to_process or args.num_samples_to_process:
                run_deidentify = args.force or not hl.hadoop_is_file(output_bam_path)
                run_sort = run_deidentify or not hl.hadoop_is_file(output_sorted_bam_path)
            else:
                run_deidentify = args.force or output_bam_path not in existing_deidentify_output_bams
                run_sort = run_deidentify or output_sorted_bam_path not in existing_deidentify_output_sorted_bams

            if run_deidentify or run_sort:
                bamout_stat = hl.hadoop_stat(row.output_bamout_bam)
                cpu = 0.25
                if bamout_stat['size_bytes'] > 0.25 * 20_000_000_000:
                    cpu = 0.5
                if bamout_stat['size_bytes'] > 0.5 * 20_000_000_000:
                    cpu = 1
                if bamout_stat['size_bytes'] > 1 * 20_000_000_000:
                    cpu = 2

            if run_deidentify:
                j = batch_utils.init_job(batch, f"{row.sample_id} - deidentify - cpu:{cpu}", DOCKER_IMAGE if not args.raw else None, cpu=cpu, disk_size=21*cpu)
                batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT)

                local_tsv_path = batch_utils.localize_file(j, row.variants_tsv_bgz, use_gcsfuse=True)
                local_exclude_tsv_path = batch_utils.localize_file(j, row.exclude_variants_tsv_bgz, use_gcsfuse=True)
                local_bamout_path = batch_utils.localize_file(j, row.output_bamout_bam, use_gcsfuse=True)

                batch_utils.localize_file(j, row.output_bamout_bai, use_gcsfuse=True)

                j.command(f"""echo --------------

echo "Start - time: $(date)"
df -kh

cat <<EOF > deidentify_bamout.py
{deidentify_bamouts_script}
EOF

time python3 deidentify_bamout.py -x "{local_exclude_tsv_path}" "{row.sample_id}" "{local_bamout_path}" "{local_tsv_path}"

ls -lh

gsutil -m cp "{row.sample_id}.deidentify_output.bam" {args.output_dir}/
gsutil -m cp "{row.sample_id}.deidentify_output.db"  {args.output_dir}/

echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------

""")
            else:
                logger.info(f"Skipping deidentify {row.sample_id}...")

            if run_sort:
                j2 = batch_utils.init_job(batch, f"{row.sample_id} - sort - cpu:{cpu}", DOCKER_IMAGE if not args.raw else None, cpu=cpu)
                batch_utils.switch_gcloud_auth_to_user_account(j2, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT)

                if run_deidentify:
                    j2.depends_on(j)

                local_bamout_path = batch_utils.localize_file(j2, output_bam_path, use_gcsfuse=True)

                j2.command(f"""echo --------------

echo "Start - time: $(date)"
df -kh

samtools sort -o "{row.sample_id}.deidentify_output.sorted.bam" "{local_bamout_path}"
samtools index "{row.sample_id}.deidentify_output.sorted.bam"

ls -lh

gsutil -m cp "{row.sample_id}.deidentify_output.sorted.bam"      {args.output_dir}/
gsutil -m cp "{row.sample_id}.deidentify_output.sorted.bam.bai"  {args.output_dir}/

echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------

""")
            elif run_sort:
                logger.info(f"Sorted output files exist (eg. {output_sorted_bam_path}). Skipping sort for {row.sample_id}...")


if __name__ == "__main__":
    main()


