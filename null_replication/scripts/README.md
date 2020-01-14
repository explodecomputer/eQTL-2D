
Begin with simulation setup

```
Rscript get_variants.r
```

Create baseline dataset

```
./create_base_dataset.sh
```

Use Snakemake to run simulations

```
module add languages/anaconda3/5.2.0-tflow-1.11
snakemake -prk \
-j 1000 \
--cluster-config bc4-cluster.json \
--cluster "sbatch \
  --job-name={cluster.name} \
  --partition={cluster.partition} \
  --nodes={cluster.nodes} \
  --ntasks-per-node={cluster.ntask} \
  --cpus-per-task={cluster.ncpu} \
  --time={cluster.time} \
  --mem={cluster.mem} \
  --output={cluster.output}"
```
