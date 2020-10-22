import hail as hl
import logging
import os
import pandas as pd
import re
import subprocess

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

#hl.init(log="/dev/null")

#%%

bams_v3 = subprocess.check_output("gsutil ls -l gs://gnomad-bw2/gnomad_readviz_bamout/*.bam", shell=True, encoding="UTF-8")
bams_v3 = bams_v3.strip().split("\n")

#%%

def get_v3_sample_id(bam_path):
    sample_id = os.path.basename(bam_path).replace(".bamout.bam", "")
    if sample_id.startswith("I-PAL") or sample_id in {
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
        sample_id = sample_id[:i] + "-" + sample_id[i+1:]
    elif sample_id == '11129':
        sample_id = "CCDG--11129"

    return sample_id

bams_v3_sample_id_map = {get_v3_sample_id(p): p for p in bams_v3}


#%%

df_v3 = pd.read_table("./cram_and_tsv_paths_table_v3_original.tsv")
df_v3.variants_tsv_bgz = df_v3.variants_tsv_bgz.fillna(df_v3.variants_tsv_path)
df_v3 = df_v3[['s', "entity:participant_id", 'variants_tsv_bgz', 'output_bamout_bam', "output_bamout_bai"]]

print(df_v3.columns)

assert len(set(bams_v3_sample_id_map.keys()) - set(df_v3["entity:participant_id"])) == 0

#%%

df_v3 = pd.DataFrame(bams_v3_sample_id_map.items(), columns=["sample_id", "output_bamout_bam"]).set_index("sample_id").join(
    df_v3[["entity:participant_id", "variants_tsv_bgz"]].set_index("entity:participant_id"),
    how="inner")

#%%
df_v3["output_bamout_bai"] = df_v3.output_bamout_bam.apply(lambda s: re.sub(".bam$", ".bai", s))

df_v3 = df_v3.reset_index().rename(columns={"index": "sample_id"})

df_v3[["sample_id", "variants_tsv_bgz", "output_bamout_bam", "output_bamout_bai"]].to_csv("./cram_and_tsv_paths_table_v3.tsv", header=True, index=False, sep="\t")
df_v3.columns




#%%

bams_v3_1 = subprocess.check_output("gsutil ls gs://gnomad-bw2/gnomad_v3_1_readviz_bamout/*.bam", shell=True, encoding="UTF-8")
bams_v3_1 = bams_v3_1.strip().split("\n")

#%%


def get_v3_1_sample_id(bam_path):
    sample_id = os.path.basename(bam_path).replace(".bamout.bam", "")
    sample_id = sample_id.replace(".final", "")
    sample_id = sample_id.replace(".srt.aln", "")
    sample_id = sample_id.split(".alt_bwamem_GRCh38D")[0]
    if sample_id == "NA12546": sample_id = "NA12546B"
    if sample_id == "NA12830": sample_id = "NA12830A"
    if sample_id == "NA18874": sample_id = "NA18874A"
    return sample_id

bams_v3_1_sample_id_map = {get_v3_1_sample_id(p): p for p in bams_v3_1}


#%%
df_v3_1 = pd.read_table("./step4_output__cram_and_tsv_paths_table_v3_1.tsv")


df_v3_1.columns

#%%

df_v3_1 = pd.DataFrame(bams_v3_1_sample_id_map.items(), columns=["sample_id", "output_bamout_bam"]).set_index("sample_id").join(
    df_v3_1[["sample_id", "variants_tsv_bgz"]].set_index("sample_id"),
    how="inner")

df_v3_1 = df_v3_1.reset_index().rename(columns={"index": "sample_id"})

df_v3_1["output_bamout_bai"] = df_v3_1.output_bamout_bam.apply(lambda s: re.sub(".bam$", ".bai", s))

df_v3_1 = df_v3_1[["sample_id", "variants_tsv_bgz", "output_bamout_bam", "output_bamout_bai"]]

df_v3_1.columns


#%%


df = pd.concat([df_v3, df_v3_1])

df

#%%

df[["sample_id", "variants_tsv_bgz", "output_bamout_bam", "output_bamout_bai"]].to_csv(
    "./cram_and_tsv_paths_table__v3_and_v3_1.tsv", header=True, index=False, sep="\t")


#%%
