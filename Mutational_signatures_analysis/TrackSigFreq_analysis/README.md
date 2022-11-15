# TrackSigFreq analysis

#### 1. Generate the tracksig information per patient
> ###### Rscript template_tracksigfreq.R -v output.muse_call.sampleID.txt.MuSE.txt.vcf.wo_tier5.vcf.tracksig_format -c tumor-neural-tissue-01_sampleID_merged_subclones.txt.sel_colmns.tracksig -p "purity" --file1 Plot_trajectory_linear_sampleID_th0.00_sig2_22.07.28.pdf --file2 Plot_trajectory_non-linear_sampleID_th0.00_sig2_22.05.16.pdf --title sampleID-th0 --sampleid sampleID
#### 2. Plot the signatures trajectories per risk group with median trajectory
> ###### Rscript PLOT_tracksig_with_median.R
