import collections
from datetime import datetime
import hail as hl   # used for hadoop file utils
import hashlib
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

DOCKER_IMAGE = "gcr.io/broad-mpg-gnomad/gnomad-readviz@sha256:7013fc57e3471617a314b08e2bcefe4711d401f83500c5c57e9a3e79ee8efebd"

#INPUT_BAM_BUCKET = "gs://gnomad-bw2/gnomad_all_readviz_bamout_deidentified"
INPUT_BAM_BUCKET = "gs://gnomad-bw2/gnomad_all_readviz_bamout_deidentified_v3_and_v31_fixed__20210101"

#OUTPUT_BUCKET = "gs://gnomad-bw2/gnomad_all_combined_bamout"
OUTPUT_BUCKET =  "gs://gnomad-bw2/gnomad_all_combined_bamout_deidentified_v3_and_v31_fixed__20210101"

DEFAULT_GROUP_SIZE = 50
ALL_CHROMOSOMES = [str(c) for c in range(1, 23)] + ["X", "Y", "M"]


def parse_args():
    """Parse command line args."""

    p = batch_utils.init_arg_parser(default_cpu=0.25, default_billing_project="gnomAD-readviz", gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    p.add_argument("-p", "--output-dir", help="Where to write combined bams.", default=OUTPUT_BUCKET)
    p.add_argument("-g", "--group-size", help="How many samples to include in each group.", default=DEFAULT_GROUP_SIZE, type=int)
    p.add_argument("-n", "--num-groups-to-process", help="For testing, process only the given sample id(s).", type=int)
    p.add_argument("-d", "--db-names-to-process", help="Process only the given dbs", action="append")
    p.add_argument("--skip-step1", action="store_true", help="Skip the step that combines bams in a group")
    p.add_argument("--skip-step2", action="store_true", help="Skip the step that combines dbs in a group")
    p.add_argument("--skip-step3", action="store_true", help="Skip the step that combines all dbs from step2 into a single db for each chromosome")
    p.add_argument("cram_and_tsv_paths_table", help="A text file containing at least these columns: sample_id, cram_path", default="step4_output__cram_and_tsv_paths_table.tsv")
    args = p.parse_args()

    return p, args


def add_command_to_combine_dbs(j, output_db_filename, input_db_paths, select_chrom=None, set_combined_bamout_id=None, create_index=False):

    sqlite_queries = []
    sqlite_queries.append('CREATE TABLE "variants" ('
        '"id" INTEGER NOT NULL PRIMARY KEY, '
        '"chrom" VARCHAR(2) NOT NULL, '
        '"pos" INTEGER NOT NULL, '
        '"ref" TEXT NOT NULL, '
        '"alt" TEXT NOT NULL, '
        '"zygosity" INTEGER NOT NULL, '
        '"qual" INTEGER NOT NULL, '
        '"combined_bamout_id" TEXT, '
        '"read_group_id" INTEGER NOT NULL);')

    column_names_string = "chrom, pos, ref, alt, zygosity, qual, combined_bamout_id, read_group_id"
    where_clause = f'WHERE chrom="{select_chrom}"' if select_chrom else ""
    for input_db_path in input_db_paths:
        sqlite_queries.append(
            f'ATTACH "{input_db_path}" as toMerge; '
            f'BEGIN; '
            f'INSERT INTO variants ({column_names_string}) SELECT {column_names_string} FROM toMerge.variants {where_clause}; '
            f'COMMIT; '
            f'DETACH toMerge;')

    if set_combined_bamout_id:
        sqlite_queries.append(
            f'UPDATE variants SET combined_bamout_id="{set_combined_bamout_id}";')

    if create_index:
        sqlite_queries.append(
            'CREATE INDEX variant_index ON "variants" ("chrom", "pos", "ref", "alt", "zygosity", "qual");')

    sqlite_queries = "\n".join(sqlite_queries)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    sqlite_queries_filename = f"sqlite_queries__{output_db_filename}__merge_{len(input_db_paths)}_dbs__{timestamp}.sql"
    sqlite_queries_temp_google_bucket_path = f"gs://gnomad-bw2/temp/{sqlite_queries_filename}"
    with open(sqlite_queries_filename, "wt") as f:
        f.write(sqlite_queries)

    os.system(f"gsutil -m cp {sqlite_queries_filename} {sqlite_queries_temp_google_bucket_path}")
    local_sqlite_queries_file_path = batch_utils.localize_file(j, sqlite_queries_temp_google_bucket_path)

    j.command(f"""echo --------------
echo "Start - time: $(date)"
df -kh
ls -lh

wc  -l {local_sqlite_queries_file_path}

time sqlite3 {output_db_filename} < {local_sqlite_queries_file_path}
""")


def combine_bam_files_in_group(args, batch, combined_bamout_id, group, input_bam_size_dict, existing_combined_bamout_bams):
    output_bam_path = os.path.join(args.output_dir, f"{combined_bamout_id}.bam")
    if not args.force and output_bam_path in existing_combined_bamout_bams:
        logger.info(f"Combined bam already exists: {output_bam_path}. Skipping {combined_bamout_id}...")
        return 0

    # check how much disk will be needed
    total_bam_size = 0
    try:
        for _, row in group.iterrows():
            total_bam_size += input_bam_size_dict[f"{INPUT_BAM_BUCKET}/{row.sample_id}.deidentify_output.sorted.bam"]

    except Exception as e:
        logger.error(f"ERROR in group {combined_bamout_id}: {e}. Unable to combine bams for group {combined_bamout_id}. Skipping...")
        return 1

    cpu = 0.25
    if total_bam_size > 0.25 * 20_000_000_000:
        cpu = 0.5
    if total_bam_size > 0.5 * 20_000_000_000:
        cpu = 1
    if total_bam_size > 1 * 20_000_000_000:
        cpu = 2
    if total_bam_size > 2 * 20_000_000_000:
        cpu = 4


    j = batch_utils.init_job(batch, f"combine bams (cpu: {cpu}): {combined_bamout_id}", DOCKER_IMAGE if not args.raw else None, cpu)
    batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT)

    picard_merge_bam_inputs = ""
    for _, row in group.iterrows():
        local_input_bam_path = batch_utils.localize_file(j, f"{INPUT_BAM_BUCKET}/{row.sample_id}.deidentify_output.sorted.bam", use_gcsfuse=True)
        local_input_bai_path = batch_utils.localize_file(j, f"{INPUT_BAM_BUCKET}/{row.sample_id}.deidentify_output.sorted.bam.bai", use_gcsfuse=True)
        picard_merge_bam_inputs += f" -I {local_input_bam_path} "

    j.command(f"""echo --------------

echo "Start - time: $(date)"
df -kh

ls -lh

java -jar /gatk/gatk.jar MergeSamFiles --VALIDATION_STRINGENCY SILENT --ASSUME_SORTED --CREATE_INDEX {picard_merge_bam_inputs} -O {combined_bamout_id}.bam
gsutil -m cp {combined_bamout_id}.bam {combined_bamout_id}.bai {args.output_dir}/

echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------

""")
    return 0


def combine_db_files_in_group_for_chrom(
    args,
    batch,
    combined_bamout_id,
    group,
    chrom_to_combine_db_jobs,
    input_db_size_dict,
    existing_combined_dbs):

    # check how much disk will be needed
    try:
        chr1_db_size_estimate = 0
        for _, row in group.iterrows():
            chr1_db_size_estimate += input_db_size_dict[f"{INPUT_BAM_BUCKET}/{row.sample_id}.deidentify_output.db"] * 0.1  # multipy by 0.1 because chr1 is < 10% of the genome

    except Exception as e:
        logger.error(f"ERROR in group {combined_bamout_id}: {e}. Unable to combine dbs for group {combined_bamout_id}. Skipping...")
        return 1

    for chrom in ALL_CHROMOSOMES:
        cpu = 0.25
        if chr1_db_size_estimate > 0.25 * 20_000_000_000:
            cpu = 0.5
        if chr1_db_size_estimate > 0.5 * 20_000_000_000:
            cpu = 1

        combined_db_filename = f"{combined_bamout_id}.chr{chrom}.db"
        output_db_path = os.path.join(args.output_dir, combined_db_filename)
        if args.db_names_to_process and combined_db_filename not in args.db_names_to_process:
            continue

        if not args.force and output_db_path in existing_combined_dbs:
            logger.info(f"Combined db already exists: {output_db_path}. Skipping combine db step for {combined_bamout_id}...")
            continue

        j2 = batch_utils.init_job(batch, f"combine dbs (cpu: {cpu}): {combined_db_filename}", DOCKER_IMAGE if not args.raw else None, cpu)
        batch_utils.switch_gcloud_auth_to_user_account(j2, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT)
        chrom_to_combine_db_jobs[chrom].append(j2)

        local_input_db_paths = []
        for _, row in group.iterrows():
            local_input_db_path = batch_utils.localize_file(j2, f"{INPUT_BAM_BUCKET}/{row.sample_id}.deidentify_output.db", use_gcsfuse=False)
            local_input_db_paths.append(local_input_db_path)

        add_command_to_combine_dbs(
            j2,
            combined_db_filename,
            local_input_db_paths,
            select_chrom=chrom,
            set_combined_bamout_id=combined_bamout_id,
            create_index=False)
        j2.command(f"gsutil -m cp {combined_db_filename} {args.output_dir}/")

    return 0


def combine_all_dbs_for_chrom(args, batch, output_filename_prefix, chrom_to_combined_db_paths, chrom_to_combine_db_jobs):
    for chrom, combined_db_paths in chrom_to_combined_db_paths.items():
        output_filename = f"all_variants_{output_filename_prefix}.chr{chrom}.db"
        combine_db_jobs = chrom_to_combine_db_jobs[chrom]

        if not args.force and hl.hadoop_is_file(f"{args.output_dir}/{output_filename}"):
            logger.info(f"{output_filename} already exists. Skipping...")
            continue

        cpu = 2
        j3 = batch_utils.init_job(batch, f"combine all dbs (cpu: {cpu}): {output_filename}", DOCKER_IMAGE if not args.raw else None, cpu=cpu, disk_size=cpu*21)
        batch_utils.switch_gcloud_auth_to_user_account(j3, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT)

        # don't use batch_utils.localize here because the command becomes too large
        j3.command("gsutil -m cp " + " ".join(combined_db_paths) + " .")
        local_input_db_paths = [os.path.basename(p) for p in combined_db_paths]

        add_command_to_combine_dbs(
            j3,
            output_filename,
            local_input_db_paths,
            select_chrom=None,
            set_combined_bamout_id=None,
            create_index=True)
        j3.command(f"gsutil -m cp {output_filename} {args.output_dir}/")

        for j2 in combine_db_jobs:
            j3.depends_on(j2)


def main():
    p, args = parse_args()

    df = pd.read_table(args.cram_and_tsv_paths_table)
    if {"sample_id", "output_bamout_bam", "output_bamout_bai", "variants_tsv_bgz"} - set(df.columns):
        p.error(f"{args.tsv_path} must contain 'sample_id', 'output_bamout_bam', 'variants_tsv_bgz' columns")

    num_groups = int(math.ceil(len(df)/args.group_size))
    logging.info(f"Creating {num_groups} group(s) with {args.group_size} samples in each")

    groups = []
    for i in range(num_groups):
        if args.num_groups_to_process and i >= args.num_groups_to_process:
            break
        group = df.iloc[i::num_groups]
        groups.append(group)

        logging.info(f"--- group #{i}:")
        logging.info(group)

    # check that all buckets are in "US-CENTRAL1" or are multi-regional to avoid egress charges to the Batch cluster
    batch_utils.set_gcloud_project(GCLOUD_PROJECT)

    if not args.skip_step1:
        existing_combined_bamout_bams = subprocess.check_output(f"gsutil ls {OUTPUT_BUCKET}/*.bam", shell=True, encoding="UTF-8").strip().split("\n")
        input_bam_size_dict = batch_utils.generate_path_to_file_size_dict(f"{INPUT_BAM_BUCKET}/*.deidentify_output.sorted.bam")

    if not args.skip_step2:
        existing_combined_dbs = subprocess.check_output(f"gsutil ls {OUTPUT_BUCKET}/*.chr*.db", shell=True, encoding="UTF-8").strip().split("\n")
        input_db_size_dict = batch_utils.generate_path_to_file_size_dict(f"{INPUT_BAM_BUCKET}/*.deidentify_output.db")

    # process groups
    hl.init(log="/dev/null")
    with batch_utils.run_batch(args, batch_name=f"combine readviz bams: {len(groups)} group(s) (gs{args.group_size}_gn{num_groups}__s{len(df)})") as batch:
        chrom_to_combine_db_jobs = collections.defaultdict(list)
        chrom_to_combined_db_paths = collections.defaultdict(list)
        errors = 0
        for i, group in enumerate(tqdm.tqdm(groups, unit=" groups")):
            md5_hash = hashlib.md5(", ".join(sorted(list(group.sample_id))).encode('utf-8')).hexdigest()
            combined_bamout_id = f"s{len(df)}_gs{args.group_size}_gn{num_groups}_gi{i:04d}_h{md5_hash[-9:]}"

            for chrom in ALL_CHROMOSOMES:
                chrom_to_combined_db_paths[chrom].append(f"{args.output_dir}/{combined_bamout_id}.chr{chrom}.db")

            if not args.skip_step1 and not args.db_names_to_process:
                errors += combine_bam_files_in_group(args, batch, combined_bamout_id, group, input_bam_size_dict, existing_combined_bamout_bams)

            if not args.skip_step2:
                errors += combine_db_files_in_group_for_chrom(args, batch, combined_bamout_id, group, chrom_to_combine_db_jobs, input_db_size_dict, existing_combined_dbs)

        if not args.skip_step3 and not args.num_groups_to_process and not errors:
            # only do this after processing all groups
            combine_all_dbs_for_chrom(args, batch, f"s{len(df)}_gs{args.group_size}_gn{num_groups}", chrom_to_combined_db_paths, chrom_to_combine_db_jobs)


if __name__ == "__main__":
    main()


