set -x

	#   --off-heap-memory-fraction 0.5 \
	#   --autoscaling-policy=dataproc-autoscale \
hailctl dataproc start \
	--requester-pays-allow-all \
	--no-off-heap-memory \
	--num-workers 2 \
	--num-preemptible-workers 0 \
	--region us-central1 \
  --project broad-mpg-gnomad \
	--autoscaling-policy=lfran_v3_1000 \
	--packages google-api-python-client,gnomad,"git+https://github.com/broadinstitute/gnomad_methods.git@main","git+https://github.com/broadinstitute/gnomad_qc.git@main" \
	--max-idle 2h \
	--master-machine-type n1-standard-16 \
	--worker-machine-type n1-standard-8 \
  readviz-select-samples

#	--master-machine-type n1-highmem-16 \
#	--worker-machine-type n1-highmem-8 \


hailctl dataproc submit  readviz-select-samples  step1__select_samples.py  --overwrite \
	--sample-metadata-tsv gs://gnomad-readviz/v4.0/gnomad.exomes.v4.0.metadata.tsv.gz \
	--output-ht-path gs://gnomad-readviz/v4.0/gnomad.exomes.v4.0.readviz_crams.ht

	#--output-ht-path gs://bw2-delete-after-60-days/v4.0/gnomad.exomes.v4.0.readviz_crams.ht

