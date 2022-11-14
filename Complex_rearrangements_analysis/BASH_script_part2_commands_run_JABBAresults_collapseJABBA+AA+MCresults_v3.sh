#!/bin/sh

#################################################################
#### COLLAPSE THE RESULTS FROM JABBA + AA + MANUAL DETECTION ####
#################################################################

#example run: sh BASH_script_part2_commands_run_JABBAresults_collapseJABBA+AA+MCresults_v2.sh FILTEREDATLEAST2CALLERS normal
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
mkdir $DIR/COMPLX_REARRANG_bypatient_infocoords
#per patient
for i in `cat /fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_risk_groups_simplified_Peifer+Berlin_infodataavail_v2.2.txt | cut -f1`; do echo mkdir $DIR/COMPLX_REARRANG_bypatient_infocoords/$i; done | sh


#### COLLAPSE RESULTS FROM JABBA + AA + MANUAL DETECTION ####

#### Re-format and COLLAPSE info from AA and Jabba(Complex rearrangements per file per patient) ####
#Execute
for i in `cat /fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_risk_groups_simplified_Peifer+Berlin_infodataavail_v2.2.txt | cut -f1`; do echo Rscript $DIR1/scripts_forallruns/COLLAPSE_jabba+AA_overlap_complexrearrangements_coords_classifier_10.0_jabba_results_allruns_NEW_filternumSVs.R $Jrun $AArun $i; done | sh

#### COLLAPSE info from AA+Jabba and Shatterseek (chromothripsis) ####
#Execute
for i in `cat /fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_risk_groups_simplified_Peifer+Berlin_infodataavail_v2.2.txt | cut -f1`; do echo Rscript $DIR1/scripts_forallruns/COLLAPSE_collapsedjabba+AA_+_SHATTERSEEK_clusteredrearrang_overlap_coords_classifier_results_allruns.R $Jrun $AArun $i; done | sh

#### COLLAPSE info from AA+Jabba and Manual curation (Clustered rearrangements) ####
#Execute
for i in `cat /fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_risk_groups_simplified_Peifer+Berlin_infodataavail_v2.2.txt | cut -f1`; do echo Rscript $DIR1/scripts_forallruns/COLLAPSE_collapsedjabba+AA_+_MANUALLYCURATED_clusteredrearrang_overlap_coords_classifier_results_allruns.R $Jrun $AArun $i; done | sh

#Remove patient directories without complex rearrang.
find $DIR/COMPLX_REARRANG_bypatient_infocoords/ -type d -empty -delete
