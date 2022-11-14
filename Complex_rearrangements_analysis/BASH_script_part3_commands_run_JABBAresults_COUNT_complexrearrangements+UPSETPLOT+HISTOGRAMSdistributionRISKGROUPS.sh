#!/bin/sh

##################################################
#### COUNT OF COMPLEX REARRANGEMENTS + PLOTS  ####
##################################################

#example run: sh BASH_script_part3_commands_run_JABBAresults_COUNT_complexrearrangements+UPSETPLOT+HISTOGRAMSdistributionRISKGROUPS.sh FILTEREDATLEAST2CALLERS normal
#arg1: select the run between FILTEREDATLEAST2CALLERS, PASS, or nofiltered
#arg2: select between AA classification forced, or normal

#### ACTIVATE ENVIRONMENT ####

DIR1="/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results"
cd $DIR1
source /fast/users/rodrigue_c/work/miniconda/etc/profile.d/conda.sh
conda activate R4


#### SET UP THE WORKING DIR ####

#variables jabba run and AA run
Jrun=$1
AArun=$2

#working directories
DIR="/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_$Jrun+AA$AArun"
mkdir $DIR/CR_allplots

#### ANALYSIS COUNT COMPLEX REARRANGEMENTS ####

### Count complex rearrangements
#Execute
Rscript $DIR1/scripts_forallruns/script_COUNT_complexrearrangements_jabba+AA_frominfo_results_allruns_NEW_v2.R $Jrun $AArun

### Generate count and freq matrix for all patients
#Execute
Rscript $DIR1/scripts_forallruns/script_generateMATRIX_forallpatients_counts+freq_perpatient_riskgroupINFO.R $Jrun $AArun

### Do UPSET plot
#Execute
Rscript $DIR1/scripts_forallruns/script_plot_UPSET_complexrearrangements_JABBA+AA+MC_allpatients_3.0_results_allruns_NEW.R $Jrun $AArun

### DO HISTOGRAMS distribution complex rearrangements across risk groups (mut signature style)
#Execute
Rscript $DIR1/scripts_forallruns/script_plot_HISTOGRAM_distribution_complexrearrangements_JABBA+AA+MC_allpatients_perRISKGROUPS_allruns_v2.R $Jrun $AArun

### DO VIOLIN PLOTS comparing the number of complex rearrangements per risk group (1 plot per complex rearrangement type)
#Execute
Rscript $DIR1/scripts_forallruns/script_plot_VIOLINPLOTS_distribution_complexrearrangements_JABBA+AA+MC_allpatients_compareRISKGROUPS_perCRtype_allruns.R $Jrun $AArun

### DO PIECHARTS ploting the distribution of complex rearrangements across the whole cohort - CR counts and counts of patients having CR
#Execute
Rscript $DIR1/scripts_forallruns/script_plot_PIECHART_distribution_complexrearrangements_WHOLECOHORT_allruns.R $Jrun $AArun
