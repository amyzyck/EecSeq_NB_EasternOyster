# Population and Seascape Genomics Analyses 

Various population genomics and seascape genomics analyses were completed for all three SNP datasets (All, Neutral, Outlier) for both the Haplotig Masked SNP datasets and Non-Haplotig Masked SNP datasets.

### PopGen_SeaGen_Analysis Repository Content 

1. `HaplotigMasked_Genome` contains RMarkdown files and knitted PDFs of the full population and seascape genomics analyses for the Haplotig Masked datasets, organized into two directories: 
- `NeutralHap_SNPs`
- `OutlierHap_SNPs`

2. `NonHaplotigMasked_Genome` contains RMarkdown files and knitted PDFs of the full population and seascape genomics analyses for the Non-Haplotig Masked datasets, organized into three directories:
- `All_SNPs`
- `Neutral_SNPs`
- `Outlier_SNPs`

**`strata` and `strata_pop` are input files used for making genind objects for population and seascape genomics analyses for all SNP datasets. `strata ` contains population, environmental, and library information for each sample. `strata_pop` contains environmental information for each population.** 

3. `EEMS` contains input files, Markdown files, and RMarkdown files for running the Estimating Effective Migration Surfaces analysis and plotting the results. This directory only documents the analysis for the Haplotig Masked Outlier SNP dataset, however the code is duplicated across all of the SNP datasets. Directory content:
- `input_files`: all files needed to run `runeems_snps`, organized by each SNP dataset. Files were either created in R or created manually.
- `NB_EEMS_Input_Output.Rmd`: RMarkdown documenting steps taken to create input files and plot results after `runeems_snps` is run.
- `NB_EEMS_OutlierHap.md`: Markdown file documenting steps taken to run `runeems_snps` analysis. 

4. `RedundancyAnalysis` contains input files, Markdown files, and RMarkdown files for running the Seascape Redundancy Analysis and plotting the results. This directory only documents the analysis for the Haplotig Masked Outlier SNP dataset, however the code is duplicated across all of the SNP datasets. Directory content:
- `input_files`: all files needed to run the redundancy analysis, organized by each SNP dataset. Files were either created in R or created manually.
- `PrepSpatialData.Rmd`: RMarkdown documenting steps taken to convert geographic coordinates for each sample site into distance-based Moran's eigenmaps. `coords.csv` was created manually.
- `Prep_Genetic_EnvironmentData.Rmd` and `Prep_Genetic_EnvironmentData.pdf`: RMarkdown and knitted PDF documenting steps taken to prepare genetic (SNP) and environmental data input files for the redundancy analysis. 
- `RDA_Outlier_Hap.Rmd`: RMarkdown documenting steps taken to run a redundancy analysis and partial redundancy analysis, and identify SNPs putatively under selection. 

Outputs from the EEMS and Redundancy Analysis can be accessed [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/tree/master/Output).

