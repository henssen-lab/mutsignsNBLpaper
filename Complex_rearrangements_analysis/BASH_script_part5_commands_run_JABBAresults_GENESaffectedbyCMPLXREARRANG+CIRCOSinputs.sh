#!/bin/sh

###############################################################
#### ANALYSIS OF GENES AFFECTED BY COMPLEX REARRANGEMENTS  ####
###############################################################

#example run: sh BASH_script_part5_commands_run_JABBAresults_GENESaffectedbyCMPLXREARRANG+CIRCOSinputs.sh FILTEREDATLEAST2CALLERS normal
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
mkdir $DIR/GENES_analysis_CR
mkdir $DIR/GENES_analysis_CR/CIRCOSinputs


#### ANALYSIS GENES AFFECTED BY COMPLEX REARRANGEMENTS ####

### Get a list of genes overlaping CR regions/bkps
#Execute
intersectBed -a $DIR/table_ALLregions_ALLcomplexrearrangements_114patients_infopatient_inforiskgroup.txt -b /fast/users/rodrigue_c/work/refs/gencode.v19.annotation_onlygenes.gtf -wo | grep gene_name | grep protein_coding | grep KNOWN | cut -f4-6,15 | cut -d ';' -f1,5 | sed 's/ /\t/g' | awk '{ print $7 "\t" $3 "\t" $2 "\t" $1}' | sed 's/\"//g' | sort -n | uniq > $DIR/GENES_analysis_CR/table_GENESoverlapingcomplexrearrangements_INTERSECT_114patients_infopatient_inforiskgroup_proteincoding.txt


#### PLOTS ####

#Plot histogram with the frequency of genes affected by complex rearrangements (info per CR type too)
Rscript $DIR1/scripts_forallruns/script_plot_multipleHISTOGRAM_freqCRaffectingGENES_allCRtypes_v3.0_NEW.R $Jrun $AArun

#Generate inputs for CIRCOS plot with regions and Freq of mutated genes by CR
Rscript $DIR1/scripts_forallruns/script_generate_inputs_CIRCOS_hist+tiles_genebased.R $Jrun $AArun


