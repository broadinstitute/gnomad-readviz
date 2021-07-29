set -x

hailctl dataproc start \
	--num-workers 2 \
	--num-preemptible 0 \
	--region us-central1 \
	--autoscaling-policy=dataproc-autoscale \
	--pkgs luigi,google-api-python-client,gnomad \
	--max-idle 8h \
	readviz-select-samples

hailctl dataproc submit \
	readviz-select-samples \
	step1__select_samples.py \
	--meta-table-ht gs://gnomad/metadata/genomes_v3.1/gnomad_v3.1_sample_qc_metadata.ht \
	--output-ht-path gs://gnomad-readviz/v3_and_v3.1/gnomad_v3_1_readviz_crams.ht
