import hail as hl   # used for hadoop file utils
import logging
import os
import pandas as pd
import tqdm

from batch import batch_utils

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


GCLOUD_PROJECT = "broad-mpg-gnomad"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"

DOCKER_IMAGE = "weisburd/gnomad-readviz@sha256:933128bf5219e2a219d77682425c1facd8f4f68e1985cdc874fd4aa1b145aa55"

EXCLUDE_INTERVALS = "gs://gnomad-bw2/exclude_intervals_with_non_ACGT_bases_in_GRCh38.bed"


PADDING_AROUND_VARIANT = 200


def parse_args():
	"""Parse command line args."""

	p = batch_utils.init_arg_parser(default_cpu=1, default_billing_project="gnomAD-readviz", gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
	p.add_argument("-p", "--output-dir", help="Where to write haplotype caller output.", default="gs://gnomad-bw2/gnomad_v3_1_readviz_bamout")
	p.add_argument("-n", "--num-samples-to-process", help="For testing, process only the first N samples.", type=int)
	p.add_argument("-s", "--sample-to-process", help="For testing, process only the given sample id(s).", nargs="+")
	p.add_argument("cram_and_tsv_paths_table", help="A text file containing at least these columns: sample_id, cram_path", default="step4_output__cram_and_tsv_paths_table.tsv")
	args = p.parse_args()

	return p, args


def main():
	p, args = parse_args()

	df = pd.read_table(args.cram_and_tsv_paths_table)
	if {"sample_id", "cram_path", "crai_path", "variants_tsv_bgz"} - set(df.columns):
		p.error(f"{args.tsv_path} must contain 'sample_id', 'cram_path' columns")

	# check that all buckets are in "US-CENTRAL1" or are multi-regional to avoid egress charges to the Batch cluster
	batch_utils.set_gcloud_project(GCLOUD_PROJECT)
	if args.cluster:
		batch_utils.check_storage_bucket_region(df.cram_path)

	if not args.force:
		hl.init(log="/dev/null", quiet=True)

	# process samples
	with batch_utils.run_batch(args, batch_name=f"HaplotypeCaller -bamout") as batch:
		counter = 0
		for _, row in tqdm.tqdm(df.iterrows(), unit=" rows", total=len(df)):
			if args.sample_to_process and row.sample_id not in set(args.sample_to_process):
				continue
			counter += 1
			if args.num_samples_to_process and counter > args.num_samples_to_process:
				break

			input_filename = os.path.basename(row.cram_path)
			output_prefix = input_filename.replace(".bam", "").replace(".cram", "")

			output_bam_path = os.path.join(args.output_dir, f"{output_prefix}.bamout.bam")
			output_bai_path = os.path.join(args.output_dir, f"{output_prefix}.bamout.bai")

			if not args.force and hl.hadoop_is_file(output_bam_path) and hl.hadoop_is_file(output_bai_path):
				logger.info(f"Output files exist (eg. {output_bam_path}). Skipping {input_filename}...")
				continue

			j = batch_utils.init_job(batch, f"readviz: {row.sample_id}", DOCKER_IMAGE if not args.raw else None, args.cpu, args.memory)
			batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT)

			local_fasta = batch_utils.localize_file(j, batch_utils.HG38_REF_PATHS.fasta, use_gcsfuse=True)
			local_fasta_fai = batch_utils.localize_file(j, batch_utils.HG38_REF_PATHS.fai, use_gcsfuse=True)
			batch_utils.localize_file(j, batch_utils.HG38_REF_PATHS.dict, use_gcsfuse=True)
			local_tsv_bgz = batch_utils.localize_file(j, row.variants_tsv_bgz)
			local_cram_path = batch_utils.localize_file(j, row.cram_path)
			local_crai_path = batch_utils.localize_file(j, row.crai_path)


			j.command(f"""echo --------------

echo "Start - time: $(date)"
df -kh


# 1) Convert variants_tsv_bgz to sorted interval list

gunzip -c "{local_tsv_bgz}" | awk '{{ OFS="\t" }} {{ print( "chr"$1, $2, $2 ) }}' | bedtools slop -b {PADDING_AROUND_VARIANT} -g {local_fasta_fai} > variant_windows.bed

# Sort the .bed file so that chromosomes are in the same order as in the input_cram file.
# Without this, if the input_cram has a different chromosome ordering (eg. chr1, chr10, .. vs. chr1, chr2, ..)
# than the interval list passed to GATK tools' -L arg, then GATK may silently skip some of regions in the -L intervals.
# The sort is done by first retrieving the input_cram header and passing it to GATK BedToIntervalList.

java -Xms2g -jar /gatk/gatk.jar PrintReadsHeader \
	--gcs-project-for-requester-pays {GCLOUD_PROJECT} \
	-R {local_fasta} \
	-I "{local_cram_path}" \
	-O header.bam

java -Xms2g -jar /gatk/gatk.jar BedToIntervalList \
	--SORT true \
	--SEQUENCE_DICTIONARY header.bam \
	--INPUT variant_windows.bed \
	--OUTPUT variant_windows.interval_list

# 2) Get reads from the input_cram for the intervals in variant_windows.interval_list

time java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+DisableAttachMechanism -XX:MaxHeapSize=2000m -Xmx30000m \
	-jar /gatk/GATK35.jar \
	-T HaplotypeCaller \
	-R {local_fasta} \
	-I "{local_cram_path}" \
	-L variant_windows.interval_list \
	--disable_auto_index_creation_and_locking_when_reading_rods \
	-ERC GVCF \
	--max_alternate_alleles 3 \
	-variant_index_parameter 128000 \
	-variant_index_type LINEAR \
	--read_filter OverclippedRead \
	-bamout "{output_prefix}.bamout.bam" \
	-o "{output_prefix}.gvcf"  |& grep -v "^DEBUG"

bgzip "{output_prefix}.gvcf"
tabix "{output_prefix}.gvcf.gz"

gsutil -m cp "{output_prefix}.bamout.bam" {args.output_dir}
gsutil -m cp "{output_prefix}.bamout.bai" {args.output_dir}
gsutil -m cp "{output_prefix}.gvcf.gz" {args.output_dir}
gsutil -m cp "{output_prefix}.gvcf.gz.tbi" {args.output_dir}

ls -lh
echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------

""")

if __name__ == "__main__":
	main()


