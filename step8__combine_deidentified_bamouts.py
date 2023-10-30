import collections
from datetime import datetime
import gzip
import hashlib
import json
import logging
import math
import os
import re
import tqdm

from step_pipeline import pipeline, Backend, Localize, Delocalize


# TODO create & keep a mapping of variant to sample ids?


logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

DOCKER_IMAGE = "weisburd/gnomad-readviz@sha256:aaad6b9f76badabe41b2f74e7445fe9cd7b61b41770d7a7c1b8e360e0f97b66a"

INPUT_DEIDENTIFIED_BAMS_DIR = "gs://gnomad-readviz/v4.0/deidentified_bamout"
LOCAL_TEMP_DIR = "sqlite_queries"
TEMP_BUCKET = "gs://bw2-delete-after-15-days"
OUTPUT_DIR = "gs://gnomad-readviz/v4.0/combined_deidentified_bamout"

DEFAULT_GROUP_SIZE = 800  # as a rule of thumb, total samples / group size should approximately be between 1000 and 1500
ALL_CHROMOSOMES = [str(c) for c in range(1, 23)] + ["X", "Y", "M"]

def parse_args(batch_pipeline):
    """Parse command line args."""

    p = batch_pipeline.get_config_arg_parser()
    p.add_argument("--ignore-input-paths-cache", help="Regenerate the input paths cache even if it already exists", action="store_true")
    p.add_argument("--input-bams-dir", help="The bucket containing all deidentified bam paths", default=INPUT_DEIDENTIFIED_BAMS_DIR)
    p.add_argument("--output-dir", help="Where to write combined bams.", default=OUTPUT_DIR)
    p.add_argument("--temp-bucket", help="Temp bucket path", default=TEMP_BUCKET)
    p.add_argument("-g", "--group-size", type=int, help="How many samples to include in each group.", default=DEFAULT_GROUP_SIZE)
    p.add_argument("-n", "--num-groups-to-process", type=int, help="For testing, process only the given sample id(s).")
    p.add_argument("-d", "--db-names-to-process", action="append", help="Process only the given dbs")
    skip_steps_group = p.add_argument_group("skip steps")
    skip_steps_group.add_argument("--skip-creating-sql-files", action="store_true", help="Skip the step writes .sql.gz files and uploads them to a temp bucket")
    skip_steps_group.add_argument("--skip-step1", action="store_true", help="Skip the step that combines bams in a group")
    skip_steps_group.add_argument("--skip-step2", action="store_true", help="Skip the step that combines dbs in a group")
    skip_steps_group.add_argument("--skip-step3", action="store_true", help="Skip the step that combines all dbs from step2 into a single db for each chromosome")

    args = batch_pipeline.parse_known_args()

    return p, args


def add_command_to_combine_dbs(
        step, output_db_filename, input_db_paths, select_chrom=None, set_combined_bamout_id=None, create_index=False, skip_creating_sql_files=False, remote_temp_dir=TEMP_BUCKET):

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
    sqlite_queries_filename = f"sqlite_queries__{output_db_filename}__merge_{len(input_db_paths)}_dbs__{timestamp}.sql.gz"

    if not skip_creating_sql_files:
        local_temp_dir = "sqlite_queries"
        if not os.path.isdir(local_temp_dir):
            os.mkdir(local_temp_dir)
        with gzip.open(os.path.join(LOCAL_TEMP_DIR, sqlite_queries_filename), "wt") as f:
            f.write(sqlite_queries)

    local_sqlite_queries_file_path = step.input(os.path.join(remote_temp_dir, LOCAL_TEMP_DIR, sqlite_queries_filename))
    step.command(f"""echo --------------
set -x
echo "Start - time: $(date)"
ls -lh
cp {local_sqlite_queries_file_path} . 
ls -lh
gunzip {local_sqlite_queries_file_path.filename}

time sqlite3 {output_db_filename} < {re.sub('.gz$', '', str(local_sqlite_queries_file_path.filename))}
ls -lh
echo Done - time: $(date)
""")


def combine_bam_files_in_group(bp, args, combined_bamout_id, group, input_bam_size_dict, existing_combined_bamout_bams):
    output_bam_path = os.path.join(args.output_dir, f"{combined_bamout_id}.bam")
    if not args.force and output_bam_path in existing_combined_bamout_bams:
        logger.info(f"Combined bam already exists: {output_bam_path}. Skipping {combined_bamout_id}...")
        return 0

    # check how much disk will be needed
    total_bam_size = 0
    try:
        for sample_id in group:
            total_bam_size += input_bam_size_dict[f"{args.input_bams_dir}/{sample_id}.deidentified.bam"]

    except Exception as e:
        logger.error(f"ERROR in group {combined_bamout_id}: {e}. Unable to combine bams for group {combined_bamout_id}. Skipping...")
        return 1

    cpu = 1
    if total_bam_size > 1 * 20_000_000_000:
        cpu = 2
    if total_bam_size > 2 * 20_000_000_000:
        cpu = 4

    s1 = bp.new_step(
        f"combine bams (cpu: {cpu}): {combined_bamout_id}",
        arg_suffix=f"step1",
        image=DOCKER_IMAGE,
        step_number=1,
        cpu=cpu,
        localize_by=Localize.HAIL_BATCH_CLOUDFUSE,
        delocalize_by=Delocalize.COPY,
        timeout=2*60*60,  # 2 hours
        output_dir=args.output_dir,
        add_skip_command_line_args=False,
        all_outputs_precached=True)

    for_loop_bam_list = ""
    picard_merge_bam_inputs = ""
    for sample_id in group:
        local_input_bam_path = s1.input(f"{args.input_bams_dir}/{sample_id}.deidentified.bam")
        for_loop_bam_list += f" '{local_input_bam_path}'"
        picard_merge_bam_inputs += f" -I '{os.path.basename(str(local_input_bam_path)).replace(' ', '_').replace(':', '_')}' "

    s1.command(f"""echo --------------

echo "Start - time: $(date)"
df -kh

# create symlinks to make the filenames shorter, so the merge command doesn't get too long 
for p in {for_loop_bam_list};
do
    ln -s "${{p}}" $(basename "${{p}}" | sed "s/ /_/g" |  sed "s/:/_/g")
done

ls -lh

# run the merge command
java -jar /gatk/gatk.jar MergeSamFiles --VALIDATION_STRINGENCY SILENT --ASSUME_SORTED --CREATE_INDEX \
    {picard_merge_bam_inputs} \
    -O {combined_bamout_id}.bam

ls

echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------

""")
    s1.output(f"{combined_bamout_id}.bam")
    s1.output(f"{combined_bamout_id}.bai")
    return 0


def combine_db_files_in_group_for_chrom(
        bp,
        args,
        combined_bamout_id,
        group,
        chrom_to_combine_db_steps,
        existing_combined_dbs,
        skip_creating_sql_files=False,
        remote_temp_dir=TEMP_BUCKET):

    for chrom in ALL_CHROMOSOMES:
        combined_db_filename = f"{combined_bamout_id}.chr{chrom}.db"
        output_db_path = os.path.join(args.output_dir, combined_db_filename)
        if args.db_names_to_process and combined_db_filename not in args.db_names_to_process:
            continue

        if not args.force and output_db_path in existing_combined_dbs:
            logger.info(f"Combined db already exists: {output_db_path}. Skipping combine db step for {combined_bamout_id}...")
            continue

        cpu = 0.25
        s2 = bp.new_step(
            f"combine dbs (cpu: {cpu}): {combined_db_filename}",
            arg_suffix=f"step2",
            image=DOCKER_IMAGE,
            step_number=2,
            cpu=cpu,
            timeout=60*60,  # 60 minutes
            localize_by=Localize.HAIL_BATCH_CLOUDFUSE,
            delocalize_by=Delocalize.COPY,
            output_dir=args.output_dir,
            add_skip_command_line_args=False,
            all_outputs_precached=True)

        chrom_to_combine_db_steps[chrom].append(s2)

        local_input_db_paths = []
        for sample_id in group:
            local_input_db_path = s2.input(f"{args.input_bams_dir}/{sample_id}.deidentified.db")
            local_input_db_paths.append(local_input_db_path)

        add_command_to_combine_dbs(
            s2,
            combined_db_filename,
            local_input_db_paths,
            select_chrom=chrom,
            set_combined_bamout_id=combined_bamout_id,
            create_index=False,
            skip_creating_sql_files=skip_creating_sql_files,
            remote_temp_dir=remote_temp_dir)
        s2.output(combined_db_filename)

    return 0


def combine_all_dbs_for_chrom(bp, args, output_filename_prefix, chrom_to_combined_db_paths, chrom_to_combine_db_steps, skip_creating_sql_files=False, remote_temp_dir=TEMP_BUCKET):
    for chrom, combined_db_paths in chrom_to_combined_db_paths.items():
        output_filename = f"all_variants_{output_filename_prefix}.chr{chrom}.db"
        combine_db_steps = chrom_to_combine_db_steps[chrom]

        s3 = bp.new_step(
            f"combine all dbs: {output_filename}",
            arg_suffix=f"step3",
            image=DOCKER_IMAGE,
            step_number=3,
            cpu=2,
            timeout=60*60,  # 60 minutes
            localize_by=Localize.HAIL_BATCH_CLOUDFUSE,
            delocalize_by=Delocalize.COPY,
            output_dir=args.output_dir,
            add_skip_command_line_args=False,
            all_outputs_precached=True)

        # don't use batch_utils.localize here because the command becomes too large
        local_input_db_paths = s3.inputs(*combined_db_paths)

        add_command_to_combine_dbs(
            s3,
            output_filename,
            local_input_db_paths,
            select_chrom=None,
            set_combined_bamout_id=None,
            create_index=True,
            skip_creating_sql_files=skip_creating_sql_files,
            remote_temp_dir=remote_temp_dir)
        s3.output(output_filename)

        for s2 in combine_db_steps:
            s3.depends_on(s2)


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline_gnomad")

    p, args = parse_args(bp)

    cache_filename = "step8_input_paths_cache.json"
    if not os.path.isfile(cache_filename) or args.ignore_input_paths_cache:
        input_bam_and_db_size_dict = {
            x["path"]: x["size_bytes"] for x in bp.precache_file_paths(f"{args.input_bams_dir}/*.deidentified.*")
        }
        with open(cache_filename, "wt") as f:
            json.dump(input_bam_and_db_size_dict, f, indent=4)

    with open(cache_filename, "rt") as f:
        input_bam_and_db_size_dict = json.load(f)

    print(f"Read {len(input_bam_and_db_size_dict)} input bam/db paths from cache file: {cache_filename}")
    all_sample_ids = {
        re.sub(".deidentified(.bam|.db|.bam.bai)$", "", os.path.basename(p)) for p in input_bam_and_db_size_dict.keys()
    }
    all_sample_ids = list(sorted(all_sample_ids))

    num_groups = int(math.floor(len(all_sample_ids)/args.group_size))
    logging.info(f"Creating {num_groups} group(s) with {args.group_size} samples in each")

    groups = []
    for i in range(num_groups):
        if args.num_groups_to_process and i >= args.num_groups_to_process:
            break
        groups.append(all_sample_ids[i::num_groups])

    if not args.skip_step1:
        existing_combined_bamout_bams = {
            x["path"]: x["size_bytes"] for x in bp.precache_file_paths(f"{args.output_dir}/*.bam")
        }

    if not args.skip_step2:
        existing_combined_dbs = {
            x["path"]: x["size_bytes"] for x in bp.precache_file_paths(f"{args.output_dir}/*.chr*.db")
        }

    bp.name = f"combine readviz bams: {len(groups)} group(s) (gs{args.group_size}_gn{num_groups}__s{len(all_sample_ids)})"

    # process each group
    chrom_to_combine_db_steps = collections.defaultdict(list)
    chrom_to_combined_db_paths = collections.defaultdict(list)
    errors = 0
    print("Processing sample id groups:")
    for i, sample_ids in enumerate(tqdm.tqdm(groups, unit=" groups")):
        md5_hash = hashlib.md5(", ".join(sorted(sample_ids)).encode('utf-8')).hexdigest()
        combined_bamout_id = f"s{len(sample_ids)}_gs{args.group_size}_gn{num_groups}_gi{i:04d}_h{md5_hash[-9:]}"

        for chrom in ALL_CHROMOSOMES:
            chrom_to_combined_db_paths[chrom].append(f"{args.output_dir}/{combined_bamout_id}.chr{chrom}.db")

        if not args.skip_step1 and not args.db_names_to_process:
            errors += combine_bam_files_in_group(
                bp, args, combined_bamout_id, sample_ids, input_bam_and_db_size_dict, existing_combined_bamout_bams)

        if not args.skip_step2:
            errors += combine_db_files_in_group_for_chrom(
                bp, args, combined_bamout_id, sample_ids, chrom_to_combine_db_steps, existing_combined_dbs,
                skip_creating_sql_files=args.skip_creating_sql_files, remote_temp_dir=TEMP_BUCKET)

    if not args.skip_step3 and not errors:
        # only do this after processing all groups
        combine_all_dbs_for_chrom(
            bp, args, f"s{len(sample_ids)}_gs{args.group_size}_gn{num_groups}", chrom_to_combined_db_paths,
            chrom_to_combine_db_steps, skip_creating_sql_files=args.skip_creating_sql_files, remote_temp_dir=TEMP_BUCKET)

    if (not args.skip_step2 or not args.skip_step3) and not args.skip_creating_sql_files:
        os.system(f"gsutil -m cp -r {LOCAL_TEMP_DIR} {TEMP_BUCKET}")

    bp.run()


if __name__ == "__main__":
    main()


