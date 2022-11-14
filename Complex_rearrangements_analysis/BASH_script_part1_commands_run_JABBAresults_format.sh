#!/bin/sh

####################################################
#### RETRIEVE AND FORMAT THE RESULTS FROM JABBA ####
####################################################

#example run: sh BASH_script_part1_commands_run_JABBAresults_format.sh FILTEREDATLEAST2CALLERS normal
#arg1: select the run between FILTEREDATLEAST2CALLERS, PASS, or nofiltered
#arg2: select between AA classification forced, or normal

#### ACTIVATE ENVIRONMENT ####

DIR1="/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results"
cd $DIR1
source /fast/users/rodrigue_c/work/miniconda/etc/profile.d/conda.sh
conda activate jabba

#### SET UP THE WORKING DIR #### 

#variables jabba run and AA run
Jrun=$1
AArun=$2

#working directories
DIR="/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_$Jrun+AA$AArun"
mkdir $DIR
mkdir $DIR/complex_rearrang_counts
mkdir $DIR/complex_rearrang_info

#### RUN THE FORMAT
#set the launch file
sed -i "s/script_classify_rearrangements_jabba_allpatients_selectrun_NEW.R.*/script_classify_rearrangements_jabba_allpatients_selectrun_v2.0_NEW.R $Jrun $AArun/g" $DIR1/scripts_forallruns/run_jabba_classify.sh 
sbatch $DIR1/scripts_forallruns/run_jabba_classify.sh
