import argparse
import hail as hl
import os

p = argparse.ArgumentParser()
p.add_argument("-o", "--output-bucket-path", help="Path where the .tsvs for each sample will be written.", default="gs://gnomad-bw2/gnomad_readviz_tsvs")
p.add_argument("sample_ids_file_path", help="A text file containing one sample id per line")
args = p.parse_args()

sample_ids = []
with open(args.sample_ids_file_path, "rt") as f:
    for line in f:
        sample_ids.append(line.rstrip("\n"))


ht = hl.read_table("gs://gnomad/readviz/genomes_v3/gnomad_v3_readviz_crams.ht")
print(ht.describe())

# add het_or_hom_or_hemi field to each struct
ht = ht.annotate(
    samples_w_het_var=ht.samples_w_het_var.map(lambda x: hl.struct(S=x.S, GQ=x.GQ, het_or_hom_or_hemi=1)),
    samples_w_hom_var=ht.samples_w_hom_var.map(lambda x: hl.struct(S=x.S, GQ=x.GQ, het_or_hom_or_hemi=2)),
    samples_w_hemi_var=ht.samples_w_hemi_var.map(lambda x: hl.struct(S=x.S, GQ=x.GQ, het_or_hom_or_hemi=3))
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

print(ht.describe())

ht = ht.checkpoint("gs://gnomad/readviz/genomes_v3/gnomad_v3_readviz_crams_exploded.ht", _read_if_exists=True)

# write a .tsv for each sample
for s in sample_ids:
    ht = ht.filter(ht.samples.S==s, keep=True)
    s = s.replace(' ', '__').replace(":", "_")
    ht.export(os.path.join(args.sample_ids_file_path, f"{s}.tsv.bgz"), header=True)
