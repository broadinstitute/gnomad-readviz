import logging
import os
import pandas as pd
import re
import subprocess

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

#all_tsvs_path = "gs://gnomad-readviz/v4.0/readviz_tsvs/*.tsv.bgz"
#all_bamouts_path = "gs://gnomad-readviz/v4.0/bamout/*.bam"
#output_filename = "cram_and_tsv_paths_table_for_step7.tsv.gz"

all_tsvs_path = "gs://gnomad-readviz/ukbb/per_sample_tsv/*.tsv.bgz"
all_bamouts_path = "gs://gnomad-readviz/ukbb/haplotypecaller/outputs/*.bamout.bam"
output_filename = "cram_and_tsv_paths_table_for_step7_from_UKBB.tsv.gz"

#%%
all_bamouts = subprocess.check_output(f"gsutil -m ls {all_bamouts_path}", shell=True, encoding="UTF-8")
all_bamouts = all_bamouts.strip().split("\n")
all_bamouts_by_sample_id = {os.path.basename(p).replace(".bamout.bam", ""): p for p in all_bamouts}
#%%

all_tsvs = subprocess.check_output(f"gsutil ls {all_tsvs_path}", shell=True, encoding="UTF-8")
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
df.to_csv(output_filename, header=True, index=False, sep="\t")
subprocess.check_call(f"gsutil -m cp {output_filename} gs://gnomad-readviz/v4.0/", shell=True)

print(f"Wrote gs://gnomad-readviz/v4.0/{output_filename}")
#%%
