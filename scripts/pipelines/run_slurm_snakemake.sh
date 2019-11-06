snakemake -p --configfile config_test.yaml -j 1 --cluster-config cluster_test.json --cluster "sbatch --ntasks 1 --cpus-per-task {cluster.n}  -t {cluster.time} --mem {cluster.mem}"
