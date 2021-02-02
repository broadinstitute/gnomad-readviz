#%%
import datetime
import pandas as pd
import re
import hail as hl
import os
#%%

"""
ht_meta = hl.read_table("gs://gnomad/metadata/genomes_v3.1/gnomad_v3.1_sample_qc_metadata.ht")
v31_all_samples = ht_meta.s.collect()      # 153,030 samples
ht_meta =  ht_meta.filter(ht_meta.release)
v31_release_samples = ht_meta.s.collect()  # 76,156 samples


ht = hl.read_table("gs://gnomad-bw2/gnomad_v3_1_readviz_crams__that_failed_AB_filter_exploded_keyed_by_sample.ht")
current_samples_v31 = ht.distinct().S.collect()  # 4,445 samples
len(set(current_samples_v31) - set(v31_release_samples))   # 919 samples  

ht = hl.read_table("gs://gnomad-bw2/gnomad_v3_readviz_crams__that_failed_AB_filter_exploded_keyed_by_sample.ht")
current_samples_v3 = ht.distinct().S.collect()  # 68,639 samples

len(set(current_samples_v3) - set(v31_release_samples)) # 1,655

"""

#%%

tsv_paths = hl.hadoop_ls("gs://gnomad-bw2/gnomad_v3_1_readviz_tsvs")

#%%

path_tuples = [(os.path.basename(t['path']).replace(".tsv.bgz", ""), t['path']) for t in tsv_paths]
df = pd.DataFrame(path_tuples, columns=['entity:participant_id', 'variants_tsv_bgz'])
df = df.set_index('entity:participant_id')

#%%

df2 = pd.read_table("./metadata/v3_1_new_releasable_cram_paths_with_sex.txt").rename(columns={'CRAM': 'cram_path', 'CRAI': 'crai_path'})
df2 = df2[['sample_id', 'cram_path', 'crai_path']]
df2 = df2.set_index('sample_id')

#%%

df_final = df.join(df2, how="inner")
assert sum(df_final.variants_tsv_bgz.isna()) == 0
assert sum(df_final.cram_path.isna()) == 0

df_final.columns

#%%

df_final.reset_index().rename(
    columns={'index': 'sample_id'}).to_csv(
    "step4_output__cram_and_tsv_paths_table.tsv", index=False, header=True, sep="\t")


#%%

#%%
df_final = df_final.reset_index().rename(columns={'index': 'entity:participant_id'})
df_final.to_csv("step4_output__terra_tsv_paricipants.tsv", index=False, header=True, sep="\t")

#%%

df_final['membership:participant_set_id'] = "sample_set_" + datetime.datetime.now().strftime("%Y_%m_%d")
df_final['participant'] = df_final["entity:participant_id"]
df_final = df_final[['membership:participant_set_id', "participant"]]
df_final.to_csv("step4_output__terra_tsv_participant_set.tsv", index=False, header=True, sep="\t")

#%%
