import argparse
import glob
import hail as hl
import os
import subprocess

v3_table = "gs://gnomad/readviz/genomes_v3/gnomad_v3_readviz_crams.ht"
v3_1_table = "gs://gnomad-bw2/gnomad_v3_1_readviz_crams.ht"

p = argparse.ArgumentParser()
p.add_argument("variant", nargs="*", help="1 or more variants like 1-123456-T-G. Prints paths of files containing read data for each variant.")
args = p.parse_args()

def run(command):
    print(command)
    return subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True, encoding="UTF-8")
    
def print_bam_paths(match):
    for sample_id in [m.S for m in match.samples_w_het_var][:3]:
        print(f"het sample: {sample_id}")
        os.system(f"gsutil ls -l 'gs://gnomad-bw2/gnomad_all_readviz_bamout/{sample_id}*'")
        os.system(f"gsutil ls -l 'gs://gnomad-bw2/gnomad_all_readviz_bamout_deidentified/{sample_id}*'")
    for sample_id in [m.S for m in match.samples_w_hom_var][:3]:
        print(f"hom sample: {sample_id}")
        os.system(f"gsutil ls -l 'gs://gnomad-bw2/gnomad_all_readviz_bamout/{sample_id}*'")
        os.system(f"gsutil ls -l 'gs://gnomad-bw2/gnomad_all_readviz_bamout_deidentified/{sample_id}*'")
    for sample_id in [m.S for m in match.samples_w_hemi_var][:3]:
        print(f"hemi sample: {sample_id}")
        os.system(f"gsutil ls -l 'gs://gnomad-bw2/gnomad_all_readviz_bamout/{sample_id}*'")
        os.system(f"gsutil ls -l 'gs://gnomad-bw2/gnomad_all_readviz_bamout_deidentified/{sample_id}*'")

hl.init(log="/dev/null")

for variant in args.variant:
    print("------------")
    try:
        chrom, pos, ref, alt = variant.split("-")
        chrom_without_prefix = chrom.replace("chr", "")
        chrom = "chr" + chrom_without_prefix

        pos = int(pos)
    except:
        p.error(f"Unable to parse variant: {variant}")
        break

    print(f"locus: {chrom}:{pos-200}-{pos+200}")
    
    locus = hl.parse_locus(f"{chrom}:{pos}", reference_genome="GRCh38")

    print(f"checking v3: {chrom}-{pos}")
    ht_v3 = hl.read_table(v3_table)
    matches = ht_v3.filter(ht_v3.locus==locus, keep=True).collect()
    print()
    for match in matches:
        print("----")
        print(f"   {match}")
        print_bam_paths(match)
        
    print(f"checking v3.1: {chrom}-{pos}")
    ht_v3_1 = hl.read_table(v3_1_table)
    matches = ht_v3_1.filter(ht_v3_1.locus==locus, keep=True).collect()
    print()
    for match in matches:
        print("----")
        print(f"   {match}")
        print_bam_paths(match)

    print(f"locus: {chrom}:{pos-200}-{pos+200}")

    db_filename = f"all_variants_s45189_gs50_gn904.{chrom}.db"
    if not os.path.isfile(db_filename):
        print(f"{db_filename} not available locally. Download it from gs://gnomad-bw2/gnomad_all_combined_bamout/{db_filename}")
        continue

    result = run(f"sqlite3 {db_filename} 'select * from variants where chrom=\"{chrom_without_prefix}\" and pos={pos}'")
    print(result)
    if result.strip():
        combined_bamout_ids = run(f"sqlite3 {db_filename} 'select combined_bamout_id from variants where chrom=\"{chrom_without_prefix}\" and pos={pos}'")
        combined_bamout_ids = combined_bamout_ids.strip()
        for combined_bamout_id in combined_bamout_ids.split("\n"):
            os.system(f"gsutil ls -l 'gs://gnomad-bw2/gnomad_all_combined_bamout/*{combined_bamout_id}*.ba*'")
    else:
        print(f"ERROR: No db records found.")
