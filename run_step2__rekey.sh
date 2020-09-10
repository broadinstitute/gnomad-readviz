hailctl dataproc start --num-workers 2 \
  --num-preemptible 30 \
  --max-idle 30m \
  --project broad-mpg-gnomad \
  --init gs://gnomad-public/tools/inits/master-init.sh \
  --packages slackclient==2.0.0,websocket-client,sklearn,statsmodels,scikit-learn,hdbscan readviz-cluster
  --master-machine-type n1-highmem-8 \
  --worker-machine-type n1-highmem-4

hailctl dataproc submit readviz-cluster step2__rekey_by_sample.py --input-ht gs://gnomad-bw2/gnomad_v3_1_readviz_crams.ht

# runs for 18 minutes
