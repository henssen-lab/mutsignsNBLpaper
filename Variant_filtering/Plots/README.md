## Commands to generate the plots about mutations (Fig. 1, Suppl. Fig. 2)

### Oncoplot

#### Generate inputs
> ###### Rscript script_generate_oncoprint_inputs_from_intersect_genesandchromosomearms_allmuttypes_newt11-17_newchrarms_6.0.R
#### Plot
> ###### Rscript script_plot_ONCOPRINT_NBLessentialgenes+chrarms_ALLmuts_ALL_PATIENTS_frominputsoncoprint_allpatientsinplot_HISTOGRAMSMUTATIONS_newt11-17_newchrarms_8.0.R

### Violin plots

#### SNVs + indels
> ###### Rscript script_SNV+indels_plots_histogram_distribution_muts_percent_mutstype_landscapeNBL+violinplots_riskgroups_COHORTS_comparison_indelsmutect_3.0.R
#### SVs (per variant type)
> ###### Rscript script_SVS_ALLtypes_plots_violinplots_riskgroups_comparison_v3.0.R
#### CNAs (per variant type)
> ###### Rscript script_CNVS_ALLtypes_plots_violinplots_riskgroups_comparison_v3.0.R

