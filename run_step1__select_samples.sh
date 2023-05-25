set -x

cluster_name=readviz-select-samples

hailctl dataproc start \
	--requester-pays-allow-all \
	--num-workers 2 \
	--num-preemptible 0 \
	--region us-central1 \
	--autoscaling-policy=dataproc-autoscale \
	--pkgs luigi,google-api-python-client,gnomad,"git+https://github.com/broadinstitute/gnomad_methods.git@main","git+https://github.com/broadinstitute/gnomad_qc.git@main" \
	--max-idle 3h \
	--worker-machine-type n1-highmem-8 \
	$cluster_name

hailctl dataproc submit  $cluster_name step1__select_samples.py  --overwrite \
  --sample-metadata-tsv gs://gnomad-readviz/v4.0/gnomad.exomes.v4.0.metadata.tsv.gz \
  --output-ht-path gs://gnomad-readviz/v4.0/gnomad.exomes.v4.0.readviz_crams.ht
