import hail as hl
import logging
import os
import pandas as pd
import re
from tqdm import tqdm

from step_pipeline import pipeline, Backend, Localize, Delocalize

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("run HaplotypeCaller")
logger.setLevel(logging.INFO)

hl.init(log="/dev/null", idempotent=True)


# Exclude variants in these regions when running HaplotypeCaller.
EXCLUDE_INTERVALS = "gs://gnomad-readviz/v4.0/exclude_intervals_with_non_ACGT_bases_in_GRCh38__150bp_window.bed"

# Amount of padding to add around each variant when running HaplotypeCaller.
PADDING_AROUND_VARIANT = 200

DOCKER_IMAGE = "weisburd/gnomad-readviz@sha256:c8f02e79221d643dcdbfb0d73e7e28bb6ba9ee7e6d108382b1cf0f35d3dab86c"

HG38_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
HG38_FASTA_INDEX_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
HG38_FASTA_DICT_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"


def parse_args(batch_pipeline):
    """Parse command line args."""

    p = batch_pipeline.get_config_arg_parser()
    p.add_argument(
        "--output-dir",
        help="Where to write HaplotypeCaller output bams.",
        default="gs://gnomad-readviz/v4.0/bamout",
    )
    p.add_argument(
        "--cram-and-tsv-table-path",
        help="A text file containing at least these columns: sample_id, cram_path",
        default=f"gs://gnomad-readviz/v4.0/step4_output_cram_and_tsv_paths.tsv.gz",
    )
    debugging_group = p.add_mutually_exclusive_group()
    debugging_group.add_argument(
        "-s", "--sample-id",
        help="Only process this sample id",
        action="append",
    )
    debugging_group.add_argument(
        "-n",
        help="For testing, process only the first N samples",
        type=int)
    p.add_argument(
        "--offset",
        help="Skip the first this many sampmles",
        default=0,
        type=int)

    args = p.parse_args()

    return p, args


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE,
                  config_file_path="~/.step_pipeline_gnomad")

    parser, args = parse_args(bp)

    logger.info(f"Loading metadata table: {args.cram_and_tsv_table_path}")
    with hl.hadoop_open(args.cram_and_tsv_table_path) as f:
        df = pd.read_table(f)

    missing_columns = {"sample_id", "cram", "crai", "variants_tsv_bgz", "gatk_version"} - set(df.columns)
    if missing_columns:
        parser.error(f"{args.cram_and_tsv_table_path} is missing these columns: {missing_columns}")

    #df = df[~df["cram"].str.contains("gs://timi-cram")]

    if args.sample_id:
        df = df[df["sample_id"].isin(args.sample_id)]

    if not args.force:
        paths = bp.precache_file_paths(os.path.join(args.output_dir, f"*.*"))
        logger.info(f"Found {len(paths)} exising .bam and .g.vcf.gz files")

    bp.name = f"step5: run HaplotypeCaller ({min(args.n or 10**9, (len(df) - args.offset))} samples)"

    logging.info(f"Processing {len(df)} samples")
    for i, (_, row) in tqdm(enumerate(df.iterrows()), unit=" samples"):
        if args.offset and i < args.offset:
            continue
        if args.n and i >= args.offset + args.n:
            break

        gatk_version = row["gatk_version"]  # example: "4.0.10.1"

        output_filename_prefix = re.sub(".tsv.bgz$", "", os.path.basename(row["variants_tsv_bgz"]))

        s1 = bp.new_step(
            f"run HaplotypeCaller: {output_filename_prefix}",
            arg_suffix=f"step1",
            image=DOCKER_IMAGE,
            step_number=1,
            cpu=0.5,
            storage="20Gi",
            delocalize_by=Delocalize.COPY,
            output_dir=args.output_dir,
            all_outputs_precached=True)

        s1.switch_gcloud_auth_to_user_account()
        local_fasta = s1.input(HG38_FASTA_PATH, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
        local_fasta_fai = s1.input(HG38_FASTA_INDEX_PATH, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
        local_fasta_dict = s1.input(HG38_FASTA_DICT_PATH, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

        local_tsv_bgz = s1.input(row["variants_tsv_bgz"], localize_by=Localize.GSUTIL_COPY)
        local_exclude_intervals = s1.input(EXCLUDE_INTERVALS, localize_by=Localize.GSUTIL_COPY)

        local_cram_path = s1.input(row["cram"], localize_by=Localize.GSUTIL_COPY)
        local_crai_path = s1.input(row["crai"], localize_by=Localize.GSUTIL_COPY)
        s1.command(f"ln -s {local_cram_path} {local_cram_path.filename}")
        s1.command(f"ln -s {local_crai_path} {local_crai_path.filename}")
        local_cram_path = local_cram_path.filename
        local_crai_path = local_crai_path.filename

        s1.command("set +exuo pipefail")
        s1.command(f"samtools view -b {local_cram_path} > " + local_cram_path.replace(".cram", ".bam"))
        local_cram_path = local_cram_path.replace(".cram", ".bam")
        s1.command(f"samtools index {local_cram_path}")

        if "gs://timi-cram" in row["cram"]:
            s1.regions("us-east1")
        else:
            s1.regions("us-central1")
        s1.command(
            f"""unset GOOGLE_APPLICATION_CREDENTIALS
env
echo --------------
echo "Start - time: $(date)"
set -exuo pipefail
df -kh

# 1) Convert variants_tsv_bgz to sorted interval list

gunzip -c "{local_tsv_bgz}" | awk '{{ OFS="\t" }} {{ print( "chr"$1, $2, $2 ) }}' | bedtools slop -b {PADDING_AROUND_VARIANT} -g {local_fasta_fai} > variant_windows.bed

# Sort the .bed file so that chromosomes are in the same order as in the input_cram file.
# Without this, if the input_cram has a different chromosome ordering (eg. chr1, chr10, .. vs. chr1, chr2, ..)
# than the interval list passed to GATK tools' -L arg, then GATK may silently skip some of regions in the -L intervals.
# The sort is done by first retrieving the input_cram header and passing it to GATK BedToIntervalList.

java -Xms2g -jar /gatk/gatk-v4.1.8.0.jar PrintReadsHeader \
    --gcs-project-for-requester-pays {args.gcloud_project} \
    -R {local_fasta} \
    -I "{local_cram_path}" \
    -O header.bam

java -Xms2g -jar /gatk/gatk-v4.1.8.0.jar BedToIntervalList \
    --SORT true \
    --SEQUENCE_DICTIONARY header.bam \
    --INPUT variant_windows.bed \
    --OUTPUT variant_windows.interval_list

# 2) Get reads from the input_cram for the intervals in variant_windows.interval_list

time java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+DisableAttachMechanism -XX:MaxHeapSize=2000m -Xmx30000m \
    -jar /gatk/gatk-v{gatk_version}.jar HaplotypeCaller \
    -R {local_fasta} \
    -I "{local_cram_path}" \
    -L variant_windows.interval_list \
    -XL {local_exclude_intervals} \
    -ERC GVCF \
    -bamout "{output_filename_prefix}.bamout.bam" \
    -O "{output_filename_prefix}.g.vcf.gz"  |& grep -v "^DEBUG"

ls -lh
echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------
""")

        s1.output(f"{output_filename_prefix}.bamout.bam")
        s1.output(f"{output_filename_prefix}.bamout.bai")
        s1.output(f"{output_filename_prefix}.g.vcf.gz")
        s1.output(f"{output_filename_prefix}.g.vcf.gz.tbi")

    bp.run()


if __name__ == "__main__":
    main()
