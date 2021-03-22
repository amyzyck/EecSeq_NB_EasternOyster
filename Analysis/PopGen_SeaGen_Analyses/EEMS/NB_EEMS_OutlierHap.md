# Estimating effective migration surfaces [(EEMS)](https://github.com/dipetkov/eems)

Here is the code used to run `runeems_snps` in command-line. Detailed steps for running this program can be accessed [here](https://github.com/dipetkov/eems/tree/master/runeems_snps).

The program `runeems_sats` implements the EEMS method for analyzing spatial population structure. This version uses raw microsatellite data.

### Input data format (microsatellites)
`runeems_snps` requires three data input files that have the same file name but different extension. The description below assumes that datapath is the full path + the file name (but without the extension). Here: `datapath = /home/azyck/NB_capture/NB_ddhaplo/NB_PopGenHap/OutlierHap/NB_EEMS_OutlierHap/`

__1. datapath.diffs__ `(outlierdata-filt.v1.diffs)`

`datapath.diffs` is the matrix of average pairwise genetic dissimilarities. This can be computed with bed2diffs from genetic data in plink binary format. 

The dissimilarity matrix is nonnegative, symmetric, with 0s on the main diagonal. These conditions are necessary but not sufficient for diffs to be a valid dissimilarity matrix. Mathematically, diffs should be conditionally negative definite.

__2. datapath.coord__ `(outlierdata-filt.v1.coord)`

`datapath.coord` are the sample coordinates, two coordinates per sample, one sample per line. The sampling locations should be given in the same order as the rows and columns of the dissimilarity matrix.

**`datapath.diffs` and `datapath.coord` are created in RStudio following steps documented [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/tree/master/Analysis/PopGen_SeaGen_Analyses/EEMS/NB_EEMS_Input_Output.Rmd)**

__3. datapath.outer__ `(outlierdata-filt.v1.outer)`

`datapath.outer` are the habitat coordinates, as a sequence of vertices that form a closed polygon. The habitat vertices should be listed counterclockwise and the first vertex should also be the last vertex, so that the outline is a closed ring. Otherwise, EEMS attempts to "correct" the polygon and prints a warning message.

**`datapath.outer` is created manually in Excel,based on site coordinates gathered from Google Maps, copied into terminal using `nano`, and saved as the file with the appropriate extension.**


#### Create separate directory for EEMS files

`PATH: /home/azyck/NB_capture/NB_ddhaplo/NB_PopGenHap/OutlierHap/`

```javascript
$ mkdir NB_EEMS_OutlierHap
$ cd NB_EEMS_OutlierHap
```

All [input files](https://github.com/amyzyck/EecSeq_NB_EasternOyster/tree/master/Analysis/PopGen_SeaGen_Analyses/EEMS/input_files/HapMasked_Outlier) should be saved in this directory.

#### Create a parameter file to set EEMS program arguments
There are a number of program parameters that can be set by the user. Here is the parameter file I created, named `params-outlierdata-chain1.ini`. Read the [instruction manual](https://github.com/dipetkov/eems/blob/master/Documentation/EEMS-doc.pdf) for descriptions of each parameter.
```javascript
datapath = /home/azyck/NB_capture/NB_ddhaplo/NB_PopGenHap/OutlierHap/NB_EEMS_OutlierHap/outlierdata-filt.v1
mcmcpath = /home/azyck/NB_capture/NB_ddhaplo/NB_PopGenHap/OutlierHap/NB_EEMS_OutlierHap/outlierdata-D200-chain1
nIndiv = 40
nSites = 841
nDemes = 200
diploid = true
numMCMCIter = 5000000
numBurnIter = 1000000
numThinIter = 9999
qEffctProposalS2 = 0.05
mEffctProposalS2 = 6.0
mrateMuProposalS2 = 0.05
mSeedsProposalS2 = 0.04
qSeedsProposalS2 = 0.07
```
- The first 9 parameter have to be specified (they have no default parameters).
- The last 5 parameters are *variances for the proposal distribution*. (I had to play around with these values, and changed most from the default values)
- As `runeems_snps` is running, it outputs information about the frequency of accepting proposals of different types.
  - The goal is to choose the parameters so that proposals are accepted about 20% − 30% of the time.
  - In practice, it seems to be sufficient to choose the variances so that proposals are not accepted too rarely (less than 10% of the time) or
too often (more than 40% of the time).

I created two additional parameter files, one with 300 demes `params-outlierdata-chain2.ini`, the other with 600 demes `params-outlierdata-chain3.ini`. All other parameters were kept the same.

#### Run the EEMS program
```javascript
$ runeems_snps --params params-alldata-chain1.ini
```
Repeat with the other two parameter files.

**Plot in RStudio following steps documented [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/tree/master/Analysis/PopGen_SeaGen_Analyses/EEMS/NB_EEMS_Input_Output.Rmd).**

**Copy from KITT to computer to view figures:**

```
$ scp -r -P xxxx azyck@kitt.uri.edu:/home/azyck/NB_capture/NB_ddhaplo/NB_PopGenHap/OutlierHap/NB_EEMS_OutlierHap/outlierdata-All-plots* /Users/azyck/Documents/Research/NB_AdultOyster/SeaGen_Analysis/Haplotig_Masked/OutlierSNPS/EEMS/
```

Figures can be viewed [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/tree/master/Output/EEMS/HapMasked_Outlier).
