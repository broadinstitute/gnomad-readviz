INPUT_HT="gs://gnomad-readviz/v4.0/gnomad.exomes.v4.0.readviz_crams_exploded_keyed_by_sample.ht"
OUTPUT_DIR="gs://gnomad-readviz/v4.0/readviz_tsvs"

# path of text file containing all sample ids (created via hl.read_table($INPUT_HT); ht.S.collect())
SAMPLE_IDS_PATH="gs://gnomad-readviz/v4.0/gnomad_v4_sample_ids_2023_07_16.txt"

python3 step3__export_per_sample_tsvs.py \
  --input-ht ${INPUT_HT} \
  --output-bucket-path ${OUTPUT_DIR} \
  ${SAMPLE_IDS_PATH}