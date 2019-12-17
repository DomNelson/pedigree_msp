set -e

N_CHROMS=4
N_SAMPLES=100
POPSIZE=100
LEN=1e6
MU=1e-8
THREADS=2

PEDFILE=~/project/pedigree_msp/data/pedEx.txt
# PEDFILE=~/project/pedigree_msp/data/BALasc_probands1930.txt

OUT_DIR=~/temp

PREFIX=$OUT_DIR/${N_SAMPLES}_samples_${POPSIZE}_Ne_${LEN}_len_${MU}_mu

cd ~/project/pedigree_msp/scripts
CHROM_NUMS=$( seq 1 $N_CHROMS )

for I in $CHROM_NUMS
do
    TS_FILES[$I]=${PREFIX}_chr$I.h5
done

## Build commands to run in parallel
for I in $CHROM_NUMS
do
    F=${TS_FILES[$I]};
    echo Simulating $F
    CMDS[$I]="python pedigree_simulate.py "
    CMDS[$I]+="-p $PEDFILE -o $F -n $N_SAMPLES -N $POPSIZE -m $MU -l $LEN;"
    # CMDS[$I]="echo testing $I;"
    # CMDS[$I]+=" sleep 1"
done;

## Simulate tree sequences
for C in "${CMDS[@]}"
do
    echo Running $C
done
parallel -j $THREADS ::: "${CMDS[@]}"

## Combine to single BCF
CMD="python ts_to_bcf.py -o ${PREFIX}_combined.bcf"
for F in "${TS_FILES[@]}"
do
    CMD+=" -t $F"
done
eval $CMD -v

## Run UMAP on combined BCF
python bcf_umap.py --bcf-file ${PREFIX}_combined.bcf -o $OUT_DIR
