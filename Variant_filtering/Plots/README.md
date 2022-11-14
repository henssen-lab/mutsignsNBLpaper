# Commands to generate the plots about mutations (Fig. 1, Suppl. Fig. 2)

### Oncoplot

#### Generate inputs
Rscript script_generate_oncoprint_inputs_from_intersect_genesandchromosomearms_allmuttypes_newt11-17_newchrarms_6.0.R
#### Plot
Rscript script_plot_ONCOPRINT_NBLessentialgenes+chrarms_ALLmuts_ALL_PATIENTS_frominputsoncoprint_allpatientsinplot_HISTOGRAMSMUTATIONS_newt11-17_newchrarms_8.0.R

### Violin plots
#### SNVs + indels
Rscript script_indels_format_4.0.py INDELS_mutect_indels_patient_sampleID_nofilters_newinfostrands.txt sampleID
#### SVs (per variant type)
Rscript script_indels_format_4.0.py INDELS_mutect_indels_patient_sampleID_nofilters_newinfostrands.txt sampleID
#### CNAs (per variant type)
Rscript script_indels_format_4.0.py INDELS_mutect_indels_patient_sampleID_nofilters_newinfostrands.txt sampleID

