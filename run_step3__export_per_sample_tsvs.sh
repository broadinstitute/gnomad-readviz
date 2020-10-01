# path of file containing all sample ids

SAMPLE_IDS_PATH="gs://gnomad-bw2/v3_1_sample_ids.txt"

if [ -f "${SAMPLE_IDS_PATH}" ]; then
    TOTAL_SAMPLE_IDS=$(cat ${SAMPLE_IDS_PATH} | wc -l)
else
    TOTAL_SAMPLE_IDS=$(gsutil cat ${SAMPLE_IDS_PATH} | wc -l)
fi

echo "Found ${TOTAL_SAMPLE_IDS} sample ids in: ${SAMPLE_IDS_PATH}"

INPUT_HT="gs://gnomad-bw2/gnomad_v3_1_readviz_crams_exploded_keyed_by_sample.ht"

# number of sample ids to process per dataproc cluster
# N_SAMPLES_PER_CLUSTER was set to 3432 for gnomADv3 since each sample id takes ~15 seconds, each cluster should take ~14.5 hours, and 20 clusters are needed to process all 68339 sample ids
N_SAMPLES_PER_CLUSTER=900  # there are 4445 samples in gnomAD v3 and each sample id requires ~15 seconds to generate the tsv, so this should take ~ 4 hours using 5 clusters

# hailctl dataproc args
N_WORKERS=2
N_PREEMPTIBLES=0
WORKER_MACHINE_TYPE="n1-standard-1"
MASTER_MACHINE_TYPE="n1-standard-4"
PROJECT=broad-mpg-gnomad
PACKAGES_ARG="slackclient==2.0.0,websocket-client,sklearn,statsmodels,scikit-learn,hdbscan"
INIT_ARG="gs://gnomad-public/tools/inits/master-init.sh"

for i in $(seq 0 ${N_SAMPLES_PER_CLUSTER} ${TOTAL_SAMPLE_IDS} )
do
    let k=${i}/${N_SAMPLES_PER_CLUSTER}
    CLUSTER_NAME="readviz-tsv-cluster${k}"

    echo "Starting ${CLUSTER_NAME}. Log file: ${CLUSTER_NAME}.log"
    set -x
    (
	hailctl dataproc start --num-workers ${N_WORKERS} --num-preemptible ${N_PREEMPTIBLES} --max-idle 10m --project ${PROJECT} --init ${INIT_ARG} --packages ${PACKAGES_ARG} ${CLUSTER_NAME} \
	    --master-machine-type ${MASTER_MACHINE_TYPE} --worker-machine-type ${WORKER_MACHINE_TYPE} 
    	hailctl dataproc submit ${CLUSTER_NAME} step3__export_per_sample_tsvs.py --start-with-sample-i ${i} --n-samples-to-process ${N_SAMPLES_PER_CLUSTER} -i ${INPUT_HT} ${SAMPLE_IDS_PATH} \
    	&& gcloud --project ${PROJECT} dataproc clusters delete --quiet ${CLUSTER_NAME} \
    ) >& ${CLUSTER_NAME}.log &
    set +x
done
