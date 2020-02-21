import argparse
import hail as hl
import logging
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


def parse_args():
	"""Parse command line args."""

	p = argparse.ArgumentParser()
	p.add_argument(
		"-o", "--output-bucket-path",
		help="Path where the .tsvs for each sample will be written",
		default="gs://gnomad-bw2/gnomad_readviz_tsvs",
	)
	p.add_argument(
		"-p", "--output-partitions",
		help="Split data for each sample into this many tsv's. Set to 1 to output a single .tsv.bgz file. "
			 "Setting this >> 1 will produce faster runtimes on large clusters.",
		default=1000,
	)
	p.add_argument(
		"-f", "--overwrite-checkpoints",
		action="store_true",
		help="Overwrite/update all hail checkpoints.",
	)
	p.add_argument(
		"sample_ids_path",
		help="A text file containing one sample id per line",
	)
	args = p.parse_args()

	return args


def read_sample_ids(sample_ids_path, n_sample_ids_to_print = 10):
	"""Read sample ids file.

	Args:
		sample_ids_path (str): sample ids path
		n_sample_ids_to_print (int): log no more than this many sample ids to stdout.
	Return:
		list: sample id strings
	"""
	sample_ids = []
	with hl.hadoop_open(sample_ids_path) if sample_ids_path.startswith("gs://" ) else open(sample_ids_path, "rt") as f:
		for i, line in enumerate(f):
			sample_id = line.rstrip("\n")
			sample_ids.append(sample_id)

			if i <= n_sample_ids_to_print:
				logging.info(sample_id)
				if i == n_sample_ids_to_print and n_sample_ids_to_print > 0:
					logging.info("...")

	logging.info(f"Parsed {len(sample_ids)} sample ids from {sample_ids_path}")

	return sample_ids


def explode_table_by_sample(ht):
	"""Explode variant-level ht by sample."""

	logging.info("Input schema:")
	ht.describe()

	ht = ht.annotate(
		samples_w_het_var=ht.samples_w_het_var.map(lambda x: hl.struct(S=x.S, GQ=x.GQ, het_or_hom_or_hemi=1)),
		samples_w_hom_var=ht.samples_w_hom_var.map(lambda x: hl.struct(S=x.S, GQ=x.GQ, het_or_hom_or_hemi=2)),
		samples_w_hemi_var=ht.samples_w_hemi_var.map(lambda x: hl.struct(S=x.S, GQ=x.GQ, het_or_hom_or_hemi=3)),
	)

	ht = ht.select(samples=ht.samples_w_het_var.extend(ht.samples_w_hom_var.extend(ht.samples_w_hemi_var)))

	ht = ht.explode(ht.samples)
	ht = ht.key_by(ht.locus)  # keep locus key to work around bug in naive_coalesce(..) when key is empty (https://github.com/hail-is/hail/issues/8138)
	ht = ht.transmute(
		#chrom=ht.locus.contig.replace("chr", ""),
		#pos=ht.locus.position,
		ref=ht.alleles[0],
		alt=ht.alleles[1],
		het_or_hom_or_hemi=ht.samples.het_or_hom_or_hemi,
		GQ=ht.samples.GQ,
		S=ht.samples.S,
	)

	return ht


def export_per_sample_tsvs(ht, sample_ids, output_bucket_path, n_partitions_per_sample):
	"""Iterate over sample_ids and export all records with that sample id to separate tsv(s).

	Args:
		ht (hail table): output of explode_table_by_sample(..)
		sample_ids (list): collection of sample id strings
		output_bucket_path (str):
		n_partitions_per_sample (int):
	"""

	logging.info("Output schema:")
	ht.describe()

	for s in sorted(sample_ids):
		per_sample_ht = ht.filter(ht.S==s, keep=True)
		per_sample_ht = per_sample_ht.drop(per_sample_ht.S)  # drop redundant column. Sample id is in the .tsv filename.
		if n_partitions_per_sample > 1:
			tsv_output_path = os.path.join(output_bucket_path, s)
			per_sample_ht = per_sample_ht.naive_coalesce(n_partitions_per_sample)
			per_sample_ht.export(tsv_output_path, parallel="separate_header")
		else:
			tsv_output_path = os.path.join(output_bucket_path, f"{s}.tsv.bgz")
			per_sample_ht.export(tsv_output_path, header=True)

		logging.info(f"Done exporting {tsv_output_path} to {n_partitions_per_sample} file(s)")


def main():
	args = parse_args()

	sample_ids = read_sample_ids(args.sample_ids_path)

	ht = hl.read_table("gs://gnomad/readviz/genomes_v3/gnomad_v3_readviz_crams.ht")

	ht = explode_table_by_sample(ht)

	ht = ht.checkpoint(
		"gs://gnomad/readviz/genomes_v3/gnomad_v3_readviz_crams_exploded_with_key.ht",
		overwrite=args.overwrite_checkpoints,
		_read_if_exists=not args.overwrite_checkpoints,
	)

	export_per_sample_tsvs(ht, sample_ids, args.output_bucket_path, args.output_partitions)


if __name__ == "__main__":
	main()
