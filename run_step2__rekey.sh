hailctl dataproc start --num-workers 2 \
  --num-preemptible 30 \
  --max-idle 30m \
  --project broad-mpg-gnomad \
  --init gs://gnomad-public/tools/inits/master-init.sh \
  --packages slackclient==2.0.0,websocket-client,sklearn,statsmodels,scikit-learn,hdbscan \
  --master-machine-type n1-highmem-8 \
  --worker-machine-type n1-highmem-4 \
  readviz-cluster-step2-rekey

hailctl dataproc submit readviz-cluster-step2-rekey step2__rekey.py --input-ht gs://gnomad-readviz/v4.0/gnomad.exomes.v4.0.readviz_crams.ht

# runs for 18 minutes

