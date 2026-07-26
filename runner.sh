CURSO=/mnt/lustre/scratch/nlsas/home/ulc/cursos/curso385
cd $CURSO/SingleCheck

# 1. baseline = the commit right before the optimizations (they landed in 9730eb6)
git worktree add $CURSO/SingleCheckBaseline 35ec840

# 2. the env file must be exported: benchmark_singlecheck.sh has no --env option,
#    it inherits this and every SingleCheck it launches sources it
export SINGLECHECK_ENV=$CURSO/SingleCheck/singlecheck_env.sh

# 3. submit
./submit --cpus 32 --mem 128G --name bench_ab benchmark_singlecheck.sh \
    -i $CURSO/Wang5.bam \
    -r 3 \
    -V $CURSO/SingleCheckBaseline,$CURSO/SingleCheck \
    -b SingleCheckBaseline \
    -T 1e-6 \
    -- -t 32
