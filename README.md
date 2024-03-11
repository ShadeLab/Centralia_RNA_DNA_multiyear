## Github Repository for
# Arrive and wait: inactive bacterial taxa contribute to perceived soil microbiome resilience after a multidecadal press disturbance
## by Samuel Barnett and Ashley Shade
<i>This work is published.</i>


### Data
The raw data for this study are available in the NCBI SRA under bioproject [PRJNA973689](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA973689/)

The processed data for this study are available on figshare under DOI [10.6084/m9.figshare.23060354](https://figshare.com/articles/dataset/Centralia_multiyear_DNA_RNA/23060354)

### To cite this work or code

Barnett, S.E., Shade, A., 2024. Arrive and wait: Inactive bacterial taxa contribute to perceived soil microbiome resilience after a multidecadal press disturbance. Ecology Letters 27, e14393. doi:https://doi.org/10.1111/ele.14393

### Abstract

Long-term (press) disturbances like the climate crisis and other anthropogenic pressures are fundamentally altering ecosystems and their functions. Many critical ecosystem functions, such as biogeochemical cycling, are facilitated by microbial communities. Understanding the functional consequences of microbiome responses to press disturbances requires ongoing observations of the active populations that contribute functions. This study leverages a 7-year time series of a 60-year-old coal seam fire (Centralia, Pennsylvania, USA) to examine the resilience of soil bacterial microbiomes to a press disturbance. Using 16S rRNA and 16S rRNA gene amplicon sequencing, we assessed the interannual dynamics of the active subset and the "whole" bacterial community. Contrary to our hypothesis, the whole communities demonstrated greater resilience than active subsets, suggesting that inactive members contributed to overall resilience. Thus, in addition to selection mechanisms of active populations, perceived microbiome resilience is also supported by mechanisms of dispersal, persistence, and revival from the local dormant pool.

### Contents

Code is split up into two directories: [Sequence_processing](https://github.com/ShadeLab/Centralia_RNA_DNA_multiyear/tree/main/Sequence_processing) and [Analysis](https://github.com/ShadeLab/Centralia_RNA_DNA_multiyear/tree/main/Analysis).

#### Sequence processing
Code used for sequence processing including read QC, OTU clustering, taxonomy assignment, and tree building can be found under [Sequence_processing](https://github.com/ShadeLab/Centralia_RNA_DNA_multiyear/tree/main/Sequence_processing). Scripts were run using SLURM on the MSU HPCC using slurm batch files with suffix .sb and are numbered by their order in the processing workflow. Outputs such as logs, warnings, or errors if any, are designated by the suffix .out and named in accordence with the library, run number, and slurm batch file. 

#### Analysis
Formal analysis can be found under [Analysis](https://github.com/ShadeLab/Centralia_RNA_DNA_multiyear/tree/main/Analysis). All analysis was run with R and code was run in Rmarkdown. In the analysis directory you'll find the raw Rmarkdown files (.Rmd), a github friendly markdown rendering (.md) and the associated figure files from the rendering in separate sub-directories. The analysis was broken down into multiple chunks in separate Rmarkdown files:

  *  [build_phyloseq](https://github.com/ShadeLab/Centralia_RNA_DNA_multiyear/blob/main/Analysis/build_phyloseq.Rmd): Used to combine data relevant to community analysis (OTU tables, taxonomy, metadata, phylogenetic tree) into a master phyloseq object for easy import into multiple analyses.
  *  [Centralia_map](https://github.com/ShadeLab/Centralia_RNA_DNA_multiyear/blob/main/Analysis/Centralia_map.md): Used to generate a map of Centralia.
  *  [DNA_diversity_analysis](https://github.com/ShadeLab/Centralia_RNA_DNA_multiyear/blob/main/Analysis/DNA_diversity_analysis.md): Primary analysis of *whole* bacterial communities (DNA only analyses) including alpha and beta diversity and trajectory analyses.
  *  [RNA_DNA_ratio_analysis](https://github.com/ShadeLab/Centralia_RNA_DNA_multiyear/blob/main/Analysis/RNA_DNA_ratio_analysis.md): Primary analuysis of *active* bacterial community subsets (RNA:DNA ratio based) including alpha and beta diversity and trajectory analysis. Also included is the analysis looking at percent of Bray-Curtis dissimilarity attributed to active taxa.
  *  [Community_assembly](https://github.com/ShadeLab/Centralia_RNA_DNA_multiyear/blob/main/Analysis/Community_assembly.md): Analysis of the results from bNTI and RCbray community assembly analysis. These indexes were calculated separately using R scripts for parallel processing but here we are examining the data in relation to disturbance.
  *  [qPCR_analysis](https://github.com/ShadeLab/Centralia_RNA_DNA_multiyear/blob/main/Analysis/qPCR_analysis.md): Analysis of 16S rRNA gene copy counts from qPCR.
  *  [Centralia_HOBO_temps](https://github.com/ShadeLab/Centralia_RNA_DNA_multiyear/blob/main/Analysis/Centralia_HOBO_temps.md): Examination of daily temperature over 7 years at 5 sites from the in situ HOBO temperature trackers.
  *  [Combined_figures](https://github.com/ShadeLab/Centralia_RNA_DNA_multiyear/blob/main/Analysis/Combined_figures.md): Code for combining various figures from the above analyses into publication ready figures.
  *  [Make_SRA_Files](https://github.com/ShadeLab/Centralia_RNA_DNA_multiyear/blob/main/Analysis/Make_SRA_Files.Rmd): code for preparing data for easier upload to the SRA.

### Funding
This work was supported by the U.S. National Science Foundation CAREER award 1749544. This work was supported in part by Michigan State University through computational resources provided by the [Institute for Cyber-Enabled Research](https://icer.msu.edu/).

### More info
[ShadeLab](http://ashley17061.wixsite.com/shadelab/home)
