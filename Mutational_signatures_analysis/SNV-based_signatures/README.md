# SNV-based mutational signatures analysis

#### Generate tri-nucleotide context from SNS calling
###### Rscript script_SIGNATUREANALYSIS_1stpart_generateTRInucleotidecontext_ALLSNVS.R
#### Extract de-novo signatures
###### Rscript script_SIGNATUREANALYSIS_2ndpart_extractdenovoSIGNATURES_ALLSNVS_v4signs.R
#### Compare de-novo signatures with COSMIC signatures, get exposures COSMIC signatures
###### Rscript script_SIGNATUREANALYSIS_3rdpart_getCOSMICsignatures_EXPOSURE_ABS+FREQ_ALLSNVS.R
#### Plot exposure per patient and per risk group
###### Rscript script_SIGNATUREANALYSIS_4thpart_PLOTexposure_perpatient+perriskgroup.R
