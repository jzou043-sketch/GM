Current optimizer-comparison benchmark

This folder contains the modified version of the previous optimizer-only benchmark:

Original reference:
C:\Users\15280\OneDrive\Desktop\GM再分析\旧版\Comparative performance of optimization algorithms in a simulation\benchmark_1000rep_optimizer_only_parallel.R

Current requested scale:
- 100 female parents
- 100 male parents
- 1,000 SNPs
- 200 effective SNPs
- Select 12 female parents
- Select 6 male parents
- 1,000 simulation replicates
- 2 starts/runs per algorithm

Algorithms:
- Particle swarm optimization
- Tabu search
- Proposed method
- Hill climbing
- Simulated annealing
- Genetic algorithm
- Differential evolution

Hyperparameters retained from the previous version:
- Tabu search: max_iter = 80, tenure = 7
- Simulated annealing: max_iter = 140, temp0 = 1, cooling = 0.97
- Differential evolution: itermax = 15, NP = 70
- Particle swarm optimization: maxit = 45, swarm size = 50
- Genetic algorithm: maxiter = 35, popSize = 50, run = 12

Important implementation note:
The previous full-neighborhood local search was designed for a much smaller parent set. With 100 females and 100 males, a full one-swap neighborhood at every step becomes very large. Therefore, this version keeps the old core algorithm hyperparameters but uses MAX_NEIGHBORS_PER_ITER = 300 by default to make the 1000-replicate benchmark computationally feasible. This value is recorded in the design CSV and can be changed through an environment variable.

Main files:
- benchmark_1000rep_optimizer_current_parallel.R
- submit_benchmark_optimizer_current_1000rep.sbatch

Server usage:
1. Upload this folder to:
   /rhome/jzou043/GM_run/optimizer_comparison_current_1000rep

2. Run:
   cd /rhome/jzou043/GM_run/optimizer_comparison_current_1000rep
   sbatch submit_benchmark_optimizer_current_1000rep.sbatch

3. Check:
   squeue -u $USER
   sacct -u $USER --starttime today --format=JobID,JobName%30,State,ExitCode,Elapsed,MaxRSS
   tail -n 100 opt_cur_JOBID.out
   tail -n 100 opt_cur_JOBID.err

Expected outputs:
- optimizer_comparison_current_1000rep_results/benchmark_1000rep_optimizer_current_replicate_best.csv
- optimizer_comparison_current_1000rep_results/benchmark_1000rep_optimizer_current_start_level.csv
- optimizer_comparison_current_1000rep_results/benchmark_1000rep_optimizer_current_summary.csv
- optimizer_comparison_current_1000rep_results/benchmark_1000rep_optimizer_current_design.csv
- optimizer_comparison_current_1000rep_results.tar.gz
