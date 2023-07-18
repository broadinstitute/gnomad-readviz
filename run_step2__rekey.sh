
hailctl dataproc start \
	--requester-pays-allow-all \
	--no-off-heap-memory \
  --num-workers 2 \
  --num-preemptible-workers 30 \
  --region us-central1 \
  --project broad-mpg-gnomad \
	--packages google-api-python-client,gnomad,"git+https://github.com/broadinstitute/gnomad_methods.git@main","git+https://github.com/broadinstitute/gnomad_qc.git@main" \
  --max-idle 2h \
  --master-machine-type n1-highmem-8 \
  --worker-machine-type n1-highmem-4 \
  readviz-cluster-step2-rekey

hailctl dataproc submit readviz-cluster-step2-rekey step2__rekey.py --input-ht gs://gnomad-readviz/v4.0/gnomad.exomes.v4.0.readviz_crams.ht

# runs for 18 minutes

