#!/bin/sh

#############################################################
#### ANALYSIS OF SVS INVOLVED IN COMPLEX REARRANGEMENTS  ####
#############################################################

#example run: sh BASH_script_part4_commands_run_JABBAresults_SVS_in_complexrearrangements+PLOTS+generaltableALLregionsinCR.sh FILTEREDATLEAST2CALLERS normal
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
mkdir $DIR/INTERSECT_SVS_complexrearrangements


#### ANALYSIS SVS IN COMPLEX REARRANGEMENTS ####

### Get general table with all the regions in complex rearrangements for downstream analysis (with patient, risk group, and CR type info)
### Get table with proportion of SVs involved in complex rearrangements per patient
#Execute
Rscript $DIR1/scripts_forallruns/script_INTERSECT_getinfo_SVSinvolvedinCOMPLXREARRANG+generaltableregionsCR.R $Jrun $AArun

### Get table with proportion of SVs involved in complex rearrangements per complex rearrangement type
#Execute
Rscript $DIR1/scripts_forallruns/script_INTERSECT_getinfo_SVSinvolvedinCOMPLXREARRANG_perCRtype_WHOLECOHORT.R $Jrun $AArun


### PLOTS
#With the same script we plot the histogram to see % of SVs involved in CR type
# We also plot the boxplot showing the % of SVs involved in CR per risk group (per patient)
Rscript $DIR1/scripts_forallruns/script_plot_HISTOGRAM+BOXPLOTS_SVsinvolvedinCMPLXREARRANG_perpatient_perCRtype_perriskgroups_v2.R $Jrun $AArun

#Genome plot - plot the density of regions involved in CR rearrangements across the human genome per CR type
Rscript $DIR1/scripts_forallruns/script_plot_GENOMEPLOT_allCRtypes_CRdistributioninthegenome.R  $Jrun $AArun

#Genome plot - plot the FREQUENCY of complex rearrangements in each of the CHROMOSOMES (patient based)
Rscript $DIR1/scripts_forallruns/script_plot_HISTOGRAM_frequency_CR_inchromosomes.R  $Jrun $AArun

