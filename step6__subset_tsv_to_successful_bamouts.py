import hail as hl
import logging
import os
import pandas as pd
import re
import subprocess

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

hl.init(log="/dev/null")

#%%
ht = hl.read_table("gs://gnomad/metadata/genomes_v3.1/gnomad_v3.1_sample_qc_metadata.ht")
ht = ht.filter(ht.release)
release_samples = ht.s.collect()


#%%

sample_ids_gnomad_v3 = hl.hadoop_open("gs://gnomad-bw2/sample_ids_gnomad_v3__20210131.txt").read().split("\n")
release_sample_ids_gnomad_v3 = list(set(sample_ids_gnomad_v3) & set(release_samples))  # 39285 samples

sample_ids_gnomad_v3_1 = hl.hadoop_open("gs://gnomad-bw2/sample_ids_gnomad_v3_1__20210131.txt").read().split("\n")
release_sample_ids_gnomad_v3_1 = list(set(sample_ids_gnomad_v3_1) & set(release_samples))  # 3526 samples

#%%

all_bamouts = subprocess.check_output("gsutil ls gs://gnomad-bw2/gnomad_all_readviz_bamout/*.bam", shell=True, encoding="UTF-8")
all_bamouts = all_bamouts.strip().split("\n")

#%%

all_tsvs = subprocess.check_output("gsutil ls gs://gnomad-bw2/gnomad_all_readviz_tsvs/*.tsv.bgz", shell=True, encoding="UTF-8")
all_tsvs = all_tsvs.strip().split("\n")

#%%

exclude_tsvs_gnomad_v3 = subprocess.check_output("gsutil ls gs://gnomad-bw2/gnomad_v3__readviz_tsvs__variants_that_failed_AB_filter/*.tsv.bgz", shell=True, encoding="UTF-8")
exclude_tsvs_gnomad_v3 = exclude_tsvs_gnomad_v3.strip().split("\n")

#%%

exclude_tsvs_gnomad_v3_1 = subprocess.check_output("gsutil ls gs://gnomad-bw2/gnomad_v3_1_readviz_tsvs__variants_that_failed_AB_filter/*.tsv.bgz", shell=True, encoding="UTF-8")
exclude_tsvs_gnomad_v3_1 = exclude_tsvs_gnomad_v3_1.strip().split("\n")



#%%

def get_v3_sample_id(file_path):
    sample_id = os.path.basename(file_path)
    sample_id = sample_id.replace(".bamout.bam", "")
    sample_id = sample_id.replace(".tsv.bgz", "")

    sample_id = sample_id.replace(".final.bamout", "")
    sample_id = sample_id.split(".alt_bwamem")[0]
    sample_id = sample_id.split(".srt.aln")[0]
    sample_id = sample_id.split(".final")[0]
    sample_id = sample_id.replace("v3.1::", "")

    if sample_id == "NA18874":
        sample_id = "NA18874A"
    if sample_id == "NA12830":
        sample_id = "NA12830A"
    elif (sample_id.startswith("I-PAL") and " " not in sample_id) or sample_id in {
        'CS-0041-01_C',
        'CS-0055-01_A',
        'CS-0060-01_C',
        'CS-0071-01_C',
        'CS-0073-01_A',
        'MCL_021_A',
        'MEI_CON_133_B',
        'ML-0693-01_B',
        'ML-0879-01_C',
        'ML-0907-01_B',
        'ML-0930-01_A',
        'ML-0974-01_A',
        'ML-0988-01_A',
        'ML-1182-01_D',
    }:
        i = sample_id.rfind("_")
        sample_id = sample_id[:i] + " " + sample_id[i+1:]
    elif sample_id == '11129':
        sample_id = "CCDG::11129"

    return sample_id

sample_id_to_bamout_map = {get_v3_sample_id(p): p for p in all_bamouts}
sample_id_to_tsv_map = {get_v3_sample_id(p): p for p in all_tsvs}
sample_id_to_exclude_tsv_map = {get_v3_sample_id(p): p for p in exclude_tsvs_gnomad_v3 + exclude_tsvs_gnomad_v3_1}


#%%

output_table = []
for sample_id in set(release_sample_ids_gnomad_v3) | set(release_sample_ids_gnomad_v3_1):
    sample_id = sample_id.replace("v3.1::", "")
    if sample_id not in sample_id_to_bamout_map:
        print(f"sample_id_to_bamout_map is missing '{sample_id}'")

    if sample_id not in sample_id_to_tsv_map:
        print(f"sample_id_to_tsv_map is missing '{sample_id}'")

    if sample_id not in sample_id_to_exclude_tsv_map:
        print(f"sample_id_to_exclude_tsv_map is missing '{sample_id}'")

    output_table.append({
        "sample_id": sample_id,
        "output_bamout_bam": sample_id_to_bamout_map[sample_id],
        "output_bamout_bai": re.sub(".bam$", ".bai", sample_id_to_bamout_map[sample_id]),
        "variants_tsv_bgz": sample_id_to_tsv_map[sample_id],
        "exclude_variants_tsv_bgz": sample_id_to_exclude_tsv_map[sample_id],
    })

#%%

[x for x in sample_id_to_exclude_tsv_map.keys() if "HG01776" in x]

#%%

df = pd.DataFrame(output_table)
df.to_csv("./cram_and_tsv_paths_table_v3_and_v31_with_exclude_variants.tsv", header=True, index=False, sep="\t")

#%%
