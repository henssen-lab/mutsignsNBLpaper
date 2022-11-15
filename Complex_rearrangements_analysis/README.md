# Complex rearrangements analysis

#### 1. Format JABBA results
> ###### BASH_script_part1_commands_run_JABBAresults_format.sh JABBArun AArun
> ###### Rscript script_classify_rearrangements_jabba_allpatients_selectrun_v2.0_NEW.R JABBArun AArun
#### 2. Collapse results from JABBA, Amplicon Architect, Shutterseek, and Clustered rearrangements
> ###### BASH_script_part2_commands_run_JABBAresults_collapseJABBA+AA+MCresults_v3.sh JABBArun AArun 
> ###### Rscript COLLAPSE_jabba+AA_overlap_complexrearrangements_coords_classifier_10.0_jabba_results_allruns_NEW_filternumSVs.R JABBArun AArun sampleID
> ###### Rscript COLLAPSE_collapsedjabba+AA_+_SHATTERSEEK_clusteredrearrang_overlap_coords_classifier_results_allruns.R JABBArun AArun sampleID
> ###### Rscript COLLAPSE_collapsedjabba+AA_+_MANUALLYCURATED_clusteredrearrang_overlap_coords_classifier_results_allruns.R JABBArun AArun sampleID
#### 3. Count complex rearrangements and generate plots (including UPSET plot)
> ###### BASH_script_part3_commands_run_JABBAresults_COUNT_complexrearrangements+UPSETPLOT+HISTOGRAMSdistributionRISKGROUPS.sh JABBArun AArun
> ###### Rscript script_COUNT_complexrearrangements_jabba+AA_frominfo_results_allruns_NEW_v2.R JABBArun AArun
> ###### Rscript script_generateMATRIX_forallpatients_counts+freq_perpatient_riskgroupINFO.R JABBArun AArun
##### Plots
> ###### Rscript script_plot_UPSET_complexrearrangements_JABBA+AA+MC_allpatients_3.0_results_allruns_NEW.R JABBArun AArun
> ###### Rscript script_plot_HISTOGRAM_distribution_complexrearrangements_JABBA+AA+MC_allpatients_perRISKGROUPS_allruns_v2.R JABBArun AArun
> ###### Rscript script_plot_VIOLINPLOTS_distribution_complexrearrangements_JABBA+AA+MC_allpatients_compareRISKGROUPS_perCRtype_allruns.R JABBArun AArun
> ###### Rscript script_plot_PIECHART_distribution_complexrearrangements_WHOLECOHORT_allruns.R JABBArun AArun
#### 4. Analysis of SVs involved in complex rearrangements
> ###### BASH_script_part4_commands_run_JABBAresults_SVS_in_complexrearrangements+PLOTS+generaltableALLregionsinCR.sh JABBArun AArun
> ###### Rscript script_INTERSECT_getinfo_SVSinvolvedinCOMPLXREARRANG+generaltableregionsCR.R JABBArun AArun
> ###### Rscript script_INTERSECT_getinfo_SVSinvolvedinCOMPLXREARRANG_perCRtype_WHOLECOHORT.R JABBArun AArun
##### Plots
> ###### Rscript script_plot_HISTOGRAM+BOXPLOTS_SVsinvolvedinCMPLXREARRANG_perpatient_perCRtype_perriskgroups_v2.R JABBArun AArun
> ###### Rscript script_plot_GENOMEPLOT_allCRtypes_CRdistributioninthegenome.R JABBArun AArun
> ###### Rscript script_plot_HISTOGRAM_frequency_CR_inchromosomes.R JABBArun AArun
#### 5. Analysis of genes affected by complex rearrangements
> ###### BASH_script_part5_commands_run_JABBAresults_GENESaffectedbyCMPLXREARRANG+CIRCOSinputs.sh JABBArun AArun
##### Plots
> ###### Rscript script_plot_multipleHISTOGRAM_freqCRaffectingGENES_allCRtypes_v3.0_NEW.R $Jrun $AArun<br/>
<br/>

#### Script to run shatterseek
> ###### Rscript script_identify_chromothripsis_perpatient_shatterseek_ATLEAST2CALLERS.R
 
