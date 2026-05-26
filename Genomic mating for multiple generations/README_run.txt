Clean package for multi-generation expected-progeny-mean simulation.

Files:
1. simulate_100rep_four_strategies_expected_mean_large.R
   Main R script. It compares PS mating, GS mating, GS+GM, and GM over 10 generations.
   The duplicate seed/n_restart problem has been fixed.

2. submit_multigen_expected_mean_large.sbatch
   Slurm submission script. It runs in the directory where sbatch is submitted, so the package can be copied to any new folder.

Default simulation design:
- 100 independent replicates
- 10 generations
- Founder population: 100 females and 100 males
- 1,000 SNPs, 200 effective loci/QTL
- Strong dominance scenario by default: SCENARIO=S3_strong
- H2 = 0.50, mu = 50
- Each generation selects 12 females and 6 males
- Each selected cross produces 20 offspring
- PS and GS use n_restart = 1
- GS+GM and GM use N_RESTART = 2 by default

How to run on skylark:

cd /rhome/jzou043/GM_run
mkdir -p multigen_expected_mean_large_clean
cd multigen_expected_mean_large_clean

# upload/copy the two files in this package into this folder, then run:
Rscript -e 'install.packages("lpSolve", repos="https://cloud.r-project.org")'

# quick test:
N_REP=1 N_WORKERS=1 OUT_DIR=test_results Rscript simulate_100rep_four_strategies_expected_mean_large.R

# full run:
sbatch submit_multigen_expected_mean_large.sbatch

# check progress:
squeue -u jzou043
tail -f mg_mean_clean_*.out

Main output folder:
multigen_expected_mean_large_results/

Main output files:
- multigen_expected_mean_large_design.csv
- multigen_expected_mean_large_history.csv
- multigen_expected_mean_large_pair_history.csv
- multigen_expected_mean_large_truth_summary.csv
- multigen_expected_mean_large_generation_summary.csv
- multigen_expected_mean_large_final_summary.csv
- multigen_expected_mean_large_gain_plot.png
- multigen_expected_mean_large_results.tar.gz
