#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2000
#SBATCH --account=clas12
#SBATCH --job-name=mc_rgm_gcf
#SBATCH --partition=production
#SBATCH --time=20:00:00
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --array=1-100 #Number of files 1-N

#EXACUTABLES
GCF_EX=/w/hallb-scshelf2102/clas12/users/PATH/TO/GCF/EX
LUND_EX=/w/hallb-scshelf2102/clas12/users/PATH/TO/LUND/EX
INPATH=/w/hallb-scshelf2102/clas12/users/PATH/TO/GCARD/AND/YAML/FILES
OUTPATH=/volatile/clas12/rg-m/PATH/TO/OUTPUT/FILES

#FILE NAMES
FILE_PREFIX=6gev_LD2_epn
PADDED_ID=$(printf "%05d" ${SLURM_ARRAY_TASK_ID})
#PADDED_ID=00001 #For testing

#PATHS TO OUTPUTS
ROOTFILE=$OUTPATH/rootfiles/gen_${FILE_PREFIX}_${PADDED_ID}.root
LUNDFILE=${OUTPATH}/lundfiles/lund_${FILE_PREFIX}_${PADDED_ID}.txt
MCFILE=${OUTPATH}/mchipo/mc_${FILE_PREFIX}_${PADDED_ID}.hipo
BKMERGEFILE=${OUTPATH}/bkgdmergehipo/bkmerge_${FILE_PREFIX}_${PADDED_ID}.hipo
RECONFILE=${OUTPATH}/reconhipo/recon_${FILE_PREFIX}_${PADDED_ID}.hipo

#PATHS TO INPUTS
BEAM_E=5.98636  #5.98636, 4.02962, 2.07052
NEVENTS=100000 #10 times the next number
NEVENTS_BKG=10000 #Max 10000 because this is how many background merging events there are per file
TARGET=liquid 
GCARD=${INPATH}/rgm_fall2021_D_511.gcard
BKGDFILE=/cache/clas12/rg-m/production/bkgfiles/tor-1.00_sol-1.00/D_5986MeV/d_${PADDED_ID}.hipo
YAML=${INPATH}/rgm_fall2021-ai_6Gev.yaml


##################
#RUN CODE
##################

#GENERATE EVENTS WITH GCF
$GCF_EX 1 1 $BEAM_E $ROOTFILE $NEVENTS -v -C -O

#CONVERT TO LUND
root -b -q "${LUND_EX}(\"${ROOTFILE}\",\"${LUNDFILE}\",\"${TARGET}\")"

#SIMULATE WITH GEMC
gemc -USE_GUI=0  -SCALE_FIELD="binary_torus, -1.0" -SCALE_FIELD="binary_solenoid, -1.0" -N=$NEVENTS -INPUT_GEN_FILE="lund, ${LUNDFILE}" -OUTPUT="hipo, ${MCFILE}" $GCARD

#BACKGROUND MERGING
bg-merger -d "ALL" -n $NEVENTS_BKG -b $BKGDFILE -i $MCFILE -o $BKMERGEFILE

#RECONSTRUCTION
recon-util -y $YAML -n $NEVENTS_BKG -i $BKMERGEFILE -o $RECONFILE
