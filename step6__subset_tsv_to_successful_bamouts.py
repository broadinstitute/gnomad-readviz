import logging
import os
import pandas as pd
import re
import subprocess

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

#%%
all_bamouts = subprocess.check_output("gsutil -m ls gs://gnomad-readviz/v4.0/bamout/*.bam", shell=True, encoding="UTF-8")
all_bamouts = all_bamouts.strip().split("\n")
all_bamouts_by_sample_id = {os.path.basename(p).replace(".bamout.bam", ""): p for p in all_bamouts}
#%%

all_tsvs = subprocess.check_output("gsutil ls gs://gnomad-readviz/v4.0/readviz_tsvs/*.tsv.bgz", shell=True, encoding="UTF-8")
all_tsvs = all_tsvs.strip().split("\n")
all_tsvs_by_sample_id = {os.path.basename(p).replace(".tsv.bgz", ""): p for p in all_tsvs}

#%%
sample_ids_without_tsvs = set(all_tsvs_by_sample_id) - set(all_bamouts_by_sample_id)
if len(sample_ids_without_tsvs) > 0:
    print(f"{len(sample_ids_without_tsvs)} samples are missing tsvs: {sample_ids_without_tsvs}")

sample_ids_without_bamout_bams = set(all_bamouts_by_sample_id) - set(all_tsvs_by_sample_id)
if len(sample_ids_without_bamout_bams) > 0:
    print(f"{len(sample_ids_without_bamout_bams)} samples are missing bamouts: {sample_ids_without_bamout_bams}")

sample_ids = set(all_tsvs_by_sample_id) & set(all_bamouts_by_sample_id)
print(f"Found {len(sample_ids):,d} sample ids that have both bamouts and tsvs")

#%%

output_table = []
for sample_id in sample_ids:
    output_table.append({
        "sample_id": sample_id,
        "output_bamout_bam": all_bamouts_by_sample_id[sample_id],
        "output_bamout_bai": re.sub(".bam$", ".bai", all_bamouts_by_sample_id[sample_id]),
        "variants_tsv_bgz": all_tsvs_by_sample_id[sample_id],
    })

df = pd.DataFrame(output_table)
df.to_csv("cram_and_tsv_paths_table_for_step7.tsv.gz", header=True, index=False, sep="\t")
subprocess.check_call("gsutil -m cp cram_and_tsv_paths_table_for_step7.tsv.gz gs://gnomad-readviz/v4.0/", shell=True)

print("Wrote gs://gnomad-readviz/v4.0/cram_and_tsv_paths_table_for_step7.tsv.gz")
#%%
