import argparse
import hail as hl
import os

p = argparse.ArgumentParser()
p.add_argument("-o", "--output-bucket-path", help="Path where the .tsvs for each sample will be written.", default="gs://gnomad-bw2/gnomad_readviz_tsvs")
p.add_argument("sample_ids_file_path", help="A text file containing one sample id per line")
args = p.parse_args()

# parse sample ids
sample_ids = []
with hl.hadoop_open(args.sample_ids_file_path) if args.sample_ids_file_path.startswith("gs://" ) else open(args.sample_ids_file_path, "rt") as f:
    for i, line in enumerate(f):
        sample_id = line.rstrip("\n")
        sample_ids.append(sample_id)

        if i < 20:
            print(sample_id)
        elif i == 20:
            print("...")

print(f"Parsed {len(sample_ids)} sample ids from {args.sample_ids_file_path}\n")

# read gnomad ht
ht = hl.read_table("gs://gnomad/readviz/genomes_v3/gnomad_v3_readviz_crams.ht")
print("Input schema:")
ht.describe()

# add het_or_hom_or_hemi field to each struct
ht = ht.annotate(
    samples_w_het_var=ht.samples_w_het_var.map(lambda x: hl.struct(S=x.S, GQ=x.GQ, het_or_hom_or_hemi=1)),
    samples_w_hom_var=ht.samples_w_hom_var.map(lambda x: hl.struct(S=x.S, GQ=x.GQ, het_or_hom_or_hemi=2)),
    samples_w_hemi_var=ht.samples_w_hemi_var.map(lambda x: hl.struct(S=x.S, GQ=x.GQ, het_or_hom_or_hemi=3)),
)

# combine ht.samples_w_het_var, ht.samples_w_hom_var, ht.samples_w_hemi_var into a single ht.samples array
ht = ht.annotate(samples=ht.samples_w_het_var.extend(ht.samples_w_hom_var.extend(ht.samples_w_hemi_var)))
ht = ht.drop(ht.samples_w_het_var, ht.samples_w_hom_var, ht.samples_w_hemi_var)

# explode
ht = ht.explode(ht.samples).key_by()
ht = ht.transmute(
    chrom=ht.locus.contig.replace("chr", ""),
    pos=ht.locus.position,
    ref=ht.alleles[0],
    alt=ht.alleles[1],
    het_or_hom_or_hemi=ht.samples.het_or_hom_or_hemi,
    GQ=ht.samples.GQ,
    S=ht.samples.S,
)

ht = ht.checkpoint("gs://gnomad/readviz/genomes_v3/gnomad_v3_readviz_crams_exploded.ht", _read_if_exists=True)

print("Output schema:")
ht.describe()

# write a .tsv for each sample
for s in sorted(sample_ids):
    ht = ht.filter(ht.S==s, keep=True)
    tsv_output_path = os.path.join(args.output_bucket_path, s.replace(' ', '__').replace(":", "_"))  # + ".tsv.bgz"
    ht.export(tsv_output_path, parallel="separate_header")
