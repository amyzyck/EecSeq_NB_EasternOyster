# *Crassostrea virginica* EecSeq Outlier Detection

Author: Amy Zyck

Last Edited: March 12, 2021

[PCAdapt](https://bcm-uga.github.io/pcadapt/articles/pcadapt.html), [OutFLANK](https://github.com/whitlock/OutFLANK), [BayeScan](http://cmpg.unibe.ch/software/BayeScan/), and [BayEnv2](https://gcbias.org/bayenv/) were used to identify any outliers. Outlier detection programs were run on `SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A.recode.vcf`.

All 4 programs were also performed on `SNP.TRSdp5g5mafMIp9g9dnDNAmaf052Ahap.recode.vcf` a VCF file containing SNPs called after reads were aligned to the oyster genome with masked haplotigs. The final number of outlier and neutral SNPs for this SNP dataset will be reported at the end of the markdown file.

VCF file is currently located: `PATH: /home/azyck/NB_capture/NB_ddocent/NB_SNPFiltering`

**First, move VCF file to new directory for outlier detection**

```
# Move to NB_ddocent Directory
$ cp ../

# Make new directory for Outlier detection
$ mkdir NB_OutlierDetect
$ cd NB_OutlierDetect

# Link VCF file to new directory
$ ln -s ../NB_SNPFiltering/SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A.recode.vcf .

# Link popmap file to new directory
$ ln -s ../NB_SNPFiltering/popmap .
```

**Remember to activate conda environment**

Steps for creating conda environment can be found in [EecSeq_Crassostrea_virginica.md](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Scripts/EecSeq_Crassostrea_virginica.md).

```
$ conda activate nb_capture
```

# 1. PCAdapt

> [pcadapt](https://bcm-uga.github.io/pcadapt/articles/pcadapt.html) has been developed to detect genetic markers involved in biological adaptation. pcadapt provides statistical tools for outlier detection based on Principal Component Analysis (PCA).

A majority of this is completed in R. I will specify which steps are completed in R or in terminal.

Make sure to save the R script or R markdown file to the same directory as `SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A.recode.vcf` and `popmap` or else you will get errors.

##### _In RStudio_

**First load the pcadapt package**

```
$ install.packages("pcadapt")
$ library(pcadapt)
```

**Convert VCF file to pcadapt format**

```
$ vcf2pcadapt("SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A.recode.vcf", output = "tmp.pcadapt", allele.sep = c("/", "|"))

$ pcadapt_file <- read.pcadapt("tmp.pcadapt", type = "pcadapt")

output:
75402 lines detected. #number of SNPs
50 columns detected. #number of individuals
```

**Choosing the number `K` of Principal Components**

We start off with a large number of PCs

```
$ x <- pcadapt(input = pcadapt_file, K = 20)
```

**Plot the likelihoods using a screeplot**

> The scree plot displays in decreasing order the percentage of variance explained by each PC. Up to a constant, it corresponds to the eigenvalues in decreasing order. The ideal pattern in a scree plot is a steep curve followed by a bend and a straight line. The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve.

```
$ plot(x, option = "screeplot")
```
![ScreeplotK20](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/ScreeplotK20.png)

I would consider 3 to be the "elbow" in this screeplot.

**Plotting likelihoods again, but with 10 Principal Components**

```
$ plot(x, option = "screeplot", K = 10)
```

![ScreeplotK10](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/ScreeplotK10.png)

**Create population designations**

```
$ poplist.names <- c(rep("BIS", 10),rep("GB", 10),rep("NAR", 10), rep("NIN", 10), rep("PVD", 10))
```

**Score Plot**

> Another option to choose the number of PCs is based on the ‘score plot’ that displays population structure. The score plot displays the projections of the individuals onto the specified principal components. Using the score plot, the choice of K can be limited to the values of K that correspond to a relevant level of population structure.

```
# Plot the first two PCs
$ plot(x, option = "scores", pop = poplist.names)
```

![PCAplot](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/PCAplot.png)

Ninigret Pond individuals are clearly separated from the other populations. These individuals are farmed oysters (unable to find wild adult oysters in Ninigret Pond).

#### I'm going to remove the NIN individuals before continuing on.

##### _In Terminal_

**Create a list of the 10 NIN Individuals**

```
$ cat NIN_Ind

output:
NIN_1
NIN_2
NIN_3
NIN_4
NIN_5
NIN_6
NIN_7
NIN_8
NIN_9
NIN_10
```

**Remove those individuals from the VCF file**

```
$ vcftools --vcf SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A.recode.vcf --remove NIN_Ind --recode --recode-INFO-all --out SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4

output:
VCFtools - 0.1.17
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A.recode.vcf
	--remove NIN_Ind
	--recode-INFO-all
	--out SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4
	--recode
Excluding individuals in 'exclude' list
After filtering, kept 40 out of 50 Individuals
Outputting VCF file...
After filtering, kept 75402 out of a possible 75402 Sites
Run Time = 18.00 seconds
```

### Repeating steps above on VCF file with NIN individuals removed

##### _In RStudio_

**Convert VCF file to pcadapt format**

```
vcf2pcadapt("SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4.recode.vcf", output = "tmp.pcadapt1", allele.sep = c("/", "|"))

filename <- read.pcadapt("tmp.pcadapt1", type = "pcadapt")

output:
75402 lines detected.
40 columns detected.
```

**Choosing the number `K` of Principal Components**

We start off with a large number of PCs

```
x1 <- pcadapt(input = filename, K = 20)
```

**Plot the likelihoods using a screeplot**

```
$ plot(x1, option = "screeplot")
```
![Screeplot2K20](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/ScreePlot2K20.png)

I would pick 2 PCs as they correspond to population structure (lie on the steep curve)

**Plotting likelihoods again, but with 10 Principal Components**

```
$ plot(x1, option = "screeplot", K = 10)
```

![Screeplot2K10](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/ScreePlot2K10.png)

**Create population designations**

```
$ poplist.names1 <- c(rep("BIS", 10),rep("GB", 10),rep("NAR", 10), rep("PVD", 10))
```

**Score Plot**

```
# Plot the first two PCs
$ plot(x1, option = "scores", pop = poplist.names1)
```

![PCA1&2](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/PCA1%262.png)

```
# PCA plot with PCs 3 and 4
$ plot(x1, option = "scores", i = 3, j = 4, pop = poplist.names1)
```

![PCA3&4](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/PCA3%264.png)

PCs 3 and 4 do not ascertain population structure, so I will stick with 2 PCs.

```
$ y <- pcadapt(input = filename, K = 2)
```

##### Start looking for outliers

**Manhattan Plot**

A Manhattan plot displays −log10 of the p-values.

```
$ plot(y , option = "manhattan")
```

![ManhattanPlot](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/ManhattanPlot.png)

**Q-Q Plots**

Used to check the expected uniform distribution of the p-values.

```
$ plot(y, option = "qqplot", threshold = 0.05)
```

![QQPlot](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/QQPlot.png)

**Look at P-value distribution**

The presence of outliers is also visible when plotting a histogram of the test statistic 𝐷𝑗.

```
$ plot(y, option = "stat.distribution")
```

![Pdist](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/Pdist.png)


**Assigning an alpha-value**

For a given 𝛼 (real valued number between 0 and 1), SNPs with q-values less than 𝛼 will be considered as outliers with an expected false discovery rate bounded by 𝛼. The false discovery rate is defined as the percentage of false discoveries among the list of candidate SNPs.

```
$ library(qvalue) # transforms p-values into q-values.
```

```
$ qval1 <- qvalue(y$pvalues)$qvalues
$ alpha <- 0.05 # expected false discovery rate lower than 5%
```

**Save Outliers**

```
$ outliers <- which(qval1 < alpha)
$ outliers
```

**Save outliers to .txt file**

```
$ invisible(lapply(outliers1, write, "outliers.txt", append=TRUE))
```

##### _In Terminal_

**Check out `outliers.txt`**

```
$ head outliers.txt

output:
359
595
961
3787
3788
3789
5913
6890
7219
7910
```

How many outliers are there?

```
$ cat outliers.txt | wc -w
```

### 865

**Save outliers to new file**

```
$ mawk '!/#/' SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4.recode.vcf | cut -f1,2 > totalloci
$ NUM=(`cat totalloci | wc -l`)
$ paste <(seq 1 $NUM) totalloci > loci.plus.index
$ cat outliers.txt | parallel "grep -w ^{} loci.plus.index" | cut -f2,3> outlier.loci.txt

$ head outlier.loci.txt

output:
NC_035780.1	2492392
NC_035780.1	5086100
NC_035780.1	7365541
NC_035780.1	26660140
NC_035780.1	26660142
NC_035780.1	26660152
NC_035780.1	39311951
NC_035780.1	52941236
NC_035780.1	55044993
NC_035780.1	58161593
```

After all outlier detection programs are completed, all the outlier SNPs will be merged into a single file.

# 2. OutFLANK

Steps for performing OutFLANK follow those documented in K. Lotterhos's [OutFLANK Vignette](https://htmlpreview.github.io/?https://github.com/whitlock/OutFLANK/blob/master/inst/doc/OutFLANKAnalysis.html).

##### _In RStudio_

**Load several necessary packages**

```
$ library(OutFLANK)  # outflank package
$ library(vcfR)
$ library(bigsnpr)   # package for Linkage Disequilibrium pruning
```

**Converting VCF file to OutFLANK format**

Using functions from the R package `vcfR`.

```
$ my_vcf <- read.vcfR("SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4.recode.vcf")

$ geno <- extract.gt(my_vcf) # Character matrix containing the genotypes
$ position <- getPOS(my_vcf) # Positions in bp
$ chromosome <- getCHROM(my_vcf) # Chromosome information

$ G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))

$ G[geno %in% c("0/0", "0|0")] <- 0
$ G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
$ G[geno %in% c("1/1", "1|1")] <- 2

# NA should be replaced with “9” to work with the functions in the OutFLANK package
$ G[is.na(G)] <- 9

$ head(G[,1:10])

output:
Scanning file to determine attributes.
File attributes:
  meta lines: 63
  header_line: 64
  variant count: 75402
  column count: 59
Meta line 63 read in.
All meta lines processed.
gt matrix initialized.
Character matrix gt created.
  Character matrix gt rows: 75402
  Character matrix gt cols: 59
  skip: 0
  nrows: 75402
  row_num: 0
Processed variant: 75402
All variants processed
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
[1,]    0    0    0    0    0    0    0    0    0     1
[2,]    0    2    2    2    2    1    2    2    2     2
[3,]    0    0    0    0    0    0    0    0    0     0
[4,]    0    0    1    0    0    0    2    0    0     1
[5,]    0    2    2    2    2    1    2    2    2     2
[6,]    1    0    0    0    0    1    0    0    0     0
```

The object "G" is now in OutFLANK format.

**Convert pop designations to a vector**

```
$ pop <- as.vector(poplist.names1)

$ pop

output:
[1] "BIS" "BIS" "BIS" "BIS" "BIS" "BIS" "BIS" "BIS" "BIS"
[10] "BIS" "GB"  "GB"  "GB"  "GB"  "GB"  "GB"  "GB"  "GB"
[19] "GB"  "GB"  "NAR" "NAR" "NAR" "NAR" "NAR" "NAR" "NAR"
[28] "NAR" "NAR" "NAR" "NIN" "NIN" "NIN" "NIN" "NIN" "NIN"
[37] "NIN" "NIN" "NIN" "NIN" "PVD" "PVD" "PVD" "PVD" "PVD"
[46] "PVD" "PVD" "PVD" "PVD" "PVD"
```

**Calculate FST on the data**

Calculate FST on all the loci in our dataset.

```
$ my_fst <- MakeDiploidFSTMat(t(G), locusNames = paste0(chromosome,"_", position), popNames = pop)

$ head(my_fst)
# the output will be a table showing Fst statistics for each locus
```
**Data checks: Heterozygosity vs. FST**

Looking for low H loci with high FST. These are all neutral loci, and it is important to exclude them from the OutFLANK algorithm.

```{r}
$ plot(my_fst$He, my_fst$FST)
```

![FstHE](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/FstHE.png)

**Data checks: FST vs. FSTNoCorr**

> To fit the FST distribution to chi-square, OutFLANK requires the FST uncorrected for sample size (`FSTNoCorr`). This is a valid approach as long as all loci have equal sample sizes within populations. The effect of correcting for sample size will make the corrected FST estimate (`FST`) lower than the uncorrected FST estimate (`FSTNoCorr`). Note that all loci deviate between `FST` and  `FSTNoCorr`, but OutFLANK assumes that these deviations are the same for each locus. If a locus has a much lower sample size compared to the rest, it could have a broader error distribution (and therefore incorrectly inferred as an outlier).

Look for loci that deviate from the linear relationship in this plot, and remove those loci.

```{r}
$ lot(my_fst$FST, my_fst$FSTNoCorr)
$ abline(0,1)
```

![FstNoCorr](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/FstNoCorr.png)


**_We need to give OUTFlank a set of quasi-independent SNPs to estimate the neutral FST distribution._**

> "quasi-independent" - SNPs located within the same genome can never be truly independent due to shared evolutionary history, but we mean a set of SNPs that are not in linkage disequilbrium due to physical linkage in the genome


##### _In RStudio_

**Following the _Bonus: code used to obtain trimmed SNPs_ in the [OutFLANK Vignette](https://htmlpreview.github.io/?https://github.com/whitlock/OutFLANK/blob/master/inst/doc/OutFLANKAnalysis.html).**

> This filtering process uses truncated [singular value decomposition (SVD)](http://faculty.washington.edu/sbrunton/me565/pdf/CHAPTER1.pdf) while limiting [linkage disequilibrium (LD)](https://en.wikipedia.org/wiki/Linkage_disequilibrium).

**Load additional packages**

```
$ library("dplyr")
$ library("stringr")
```

**Chromosomes need to be of class integer for this to work.**

```
# Visualizing chromosomes
$ chrom_unique <- unique(chromosome)
$ print(chrom_unique)

output:
[1] "NC_035780.1" "NC_035781.1" "NC_035782.1" "NC_035783.1" "NC_035784.1"
[6] "NC_035785.1" "NC_035786.1" "NC_035787.1" "NC_035788.1" "NC_035789.1"
```

```
# Removing "NC_" from the chromosome name so it can be converted to an integer
$ chrom_new <- chromosome%>%str_replace("NC_","")

# Converting the character string into an integer
$ chrom1 <- as.integer(chrom_new)
```

**`chromosome` and `position` need to be sorted for imputation to work.**

```
$ chrom_sort <- sort(chrom1)
$ pos_sort <- sort(position)
```

**_Note: This filtering program does not allow for missing genotype values._**

**Using `snp_fastImputeSimple` to impute missing genotypes**

This requires the genotpe matrix to be converted to class [FBM.code256](https://privefl.github.io/bigstatsr/reference/FBM.code256-class.html).

> A reference class for storing and accessing up to 256 arbitrary different values using a Filebacked Big Matrix of type `unsigned char`.

```
# Converting the genotype matrix to class FBM.code256
$ G1 <- add_code256(big_copy(t(G),type="raw"),code=bigsnpr:::CODE_012)

# Imputing missing genotypes
$ G_missSimple <- snp_fastImputeSimple(G1,method = c("mode", "mean0", "mean2", "random", "zero"))
```

**Truncated SVD while limiting LD.**

```
$ newpc<-snp_autoSVD(G_missSimple,infos.chr =chrom_sort,infos.pos = pos_sort)
$ which_pruned <- attr(newpc, which="subset") # Indexes of remaining SNPS after pruning
$ length(which_pruned)

output:
Phase of clumping (on MAF) at r^2 > 0.2.. keep 18362 SNPs.
Discarding 5798 variants with MAC < 10.

Iteration 1:
Computing SVD..
0 outlier variant detected..

Converged!
[1] 12564
```

**Run the `OutFLANK()` function to estimate the parameters on the neutral FST distribution.**

```
$ out_trim <- OutFLANK(my_fst[which_pruned,], NumberOfSamples=40, qthreshold = 0.05, Hmin = 0.1)
$ str(out_trim)

output:
List of 6
 $ FSTbar               : num 0.00401
 $ FSTNoCorrbar         : num 0.0516
 $ dfInferred           : num 2.76
 $ numberLowFstOutliers : int 0
 $ numberHighFstOutliers: int 0
 $ results              :'data.frame':	12564 obs. of  15 variables:
  ..$ LocusName       : Factor w/ 75402 levels "NC_035780.1_1000435",..: 6860 6861 6864 6866 6867 6868 6871 6878 6881 6885 ...
  ..$ He              : num [1:12564] 0.393 0.242 0.497 0.242 0.497 ...
  ..$ FST             : num [1:12564] -0.0324 0.0163 -0.0239 0.0612 -0.0195 ...
  ..$ T1              : num [1:12564] -0.00639 0.00201 -0.00597 0.0076 -0.00488 ...
  ..$ T2              : num [1:12564] 0.198 0.123 0.25 0.124 0.251 ...
  ..$ FSTNoCorr       : num [1:12564] 0.017 0.0692 0.0159 0.0993 0.038 ...
  ..$ T1NoCorr        : num [1:12564] 0.00336 0.00853 0.00396 0.01233 0.00952 ...
  ..$ T2NoCorr        : num [1:12564] 0.198 0.123 0.25 0.124 0.251 ...
  ..$ meanAlleleFreq  : num [1:12564] 0.731 0.141 0.463 0.141 0.538 ...
  ..$ indexOrder      : int [1:12564] 1 2 3 4 5 6 7 8 9 10 ...
  ..$ GoodH           : Factor w/ 2 levels "goodH","lowH": 1 1 1 1 1 1 1 1 1 1 ...
  ..$ qvalues         : num [1:12564] 0.914 0.895 0.914 0.895 0.91 ...
  ..$ pvalues         : num [1:12564] 0.425 0.521 0.393 0.257 0.96 ...
  ..$ pvaluesRightTail: num [1:12564] 0.788 0.26 0.804 0.129 0.52 ...
  ..$ OutlierFlag     : logi [1:12564] FALSE FALSE FALSE FALSE FALSE FALSE ...
```

**Check the fit and make sure it looks good, especially in the right tail:**

```
$ OutFLANKResultsPlotter(out_trim, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1,binwidth = 0.001, Zoom = FALSE, RightZoomFraction = 0.05, titletext = NULL)
```

![ResultsPlotter](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/OutFLANKResultsPlotter.png)

**Check the P-value histogram:**

Here, we plot the “right-tailed” P-values, which means that outliers in the right tail of the FST distribution will have a P-value near zero. Because we ran the algorithm on a trimmed set of SNPs, this will remove some of the signal around selected sites. So we expect this histogram to be flat and maybe have a bump near 0 for selected sites.

```
$ hist(out_trim$results$pvaluesRightTail)
```

![Histo](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/HistoofRightTail.png)

**Using estimated neutral mean FST and df to calculate P-values for all loci.**

Now that we’ve estimated neutral mean FST and df to a quasi-independent set of SNPs, we can go back and calculate P-values and q-values for all the loci in our dataset.

```
$ P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)
```

**Identify outliers SNPs**

```
$ my_out <- P1$OutlierFlag==TRUE
plot(P1$He, P1$FST, pch=19, col=rgb(0,0,0,0.1))
points(P1$He[my_out], P1$FST[my_out], col="blue")
```

![P1Outlier](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/P1Outlier.png)

```
$ P1[which(P1$OutlierFlag == TRUE),]
```

|      |       LocusName      |    He    |    FST   |     T1    |    T2    | FSTNoCorr |  T1NoCorr |  T2NoCorr | meanAlleleFreq |    pvalues   | pvaluesRightTail |   qvalues   | OutlierFlag |
|:----:|:--------------------:|:--------:|:--------:|:---------:|:--------:|:---------:|:---------:|:---------:|:--------------:|:------------:|:----------------:|:-----------:|:-----------:|
| 7983 | NC_035780.1_58431616 | 0.432133 | 0.622561 | 0.1602553 | 0.257413 |  0.640775 | 0.1649466 | 0.2574174 |    0.6842105   | 2.454826e-07 |    1.227413e-07   | 0.006736656 |     TRUE    |

**Save Outlier**

##### _In Terminal_

```
$ mawk '!/#/' SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4.recode.vcf | cut -f1,2 > totalloci1
$ NUM=(`cat totalloci1 | wc -l`)
$ paste <(seq 1 $NUM) totalloci1 > loci1.plus.index
$ echo -e "7983" | parallel "grep -w ^{} loci1.plus.index" | cut -f2,3> outlier1.loci.txt

$ head outlier1.loci.txt

output:
NC_035780.1	58431616
```

# 3. BayeScan

> BayeScan aims at identifying candidate loci under natural selection from genetic data, using differences in allele frequencies between populations. BayeScan is based on the multinomial-Dirichlet model. One of the simplest possible scenarios covered consists of an island model in which subpopulation allele frequencies are correlated through a common migrant gene pool from which they differ in varying degrees. The difference in allele frequency between this common gene pool and each subpopulation is measured by a subpopulation specific FST coefficient. Therefore, this formulation can consider realistic ecological scenarios where the effective size and the immigration rate may differ among subpopulations.

**_Note_**: PGDSpider2 does not run on Python version 3.7.6, which is the version in the `nb_capture` conda environment. Therefore, I am creating and activating a new conda environment with a previous version of Python that I know PGDSpider2 works with.

##### _In Terminal_

```
$ conda create -n nb_capture2 python=3.7.3
$ conda activate nb_capture2
```

**Convert VCF file to other outputs**

```
# load in configuration file that will convert VCF to BayeScan format
$ curl -L -O https://raw.githubusercontent.com/amyzyck/EecSeq_NB_EasternOyster/master/Scripts/BSsnp.spid
$ chmod +x BSsnp.spid

$ java -jar /usr/local/bin/PGDSpider2-cli.jar -spid BSsnp.spid -inputfile SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4.recode.vcf -outputfile SNP1.BS

output:
WARN  09:47:54 - PGDSpider configuration file not found! Loading default configuration.
initialize convert process...
read input file...
read input file done.
write output file...
write output file done.
```

**Run BayeScan**

```
$ BayeScan2.1_linux64bits SNP1.BS -thin 50 -nbp 30
```

**Copy R source file**

This file is used to plot figures for the software Bayescan in R.

```
$ curl -L -O https://raw.githubusercontent.com/z0on/2bRAD_denovo/master/plot_R.r
$ chmod +x plot_R.r
```

##### _In RStudio_

**Plotting outliers**

```
$ source("plot_R.r")
$ plot_bayescan("SNP_fst.txt")

output:
$outliers
[1]  7983 24006 24007 51866 68683

$nb_outliers
[1] 5
```

![SNPfst_outlier](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/SNPfst_outliers.png).

**Identifying outliers using a False discovery rate of 10%.**

Choosing a q-value threshold of 10% means that 10% of the corresponding outlier markers (those having a q-value lower than 10%) are expected to be false positives.

```
$ bs <- plot_bayescan("SNP_fst.txt", FDR = 0.1)
```

![FDR10](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/SNPFstFDR10.png)

```
$ bs$outliers

output:
[1]   350  7983 24003 24006 24007 51866 68683
```

**Save outliers to .txt file**

```
$ invisible(lapply(bs, write, "outliers2.txt", append=TRUE))
```

##### _In Terminal_

**Check out `outliers2.txt`**

```
$ head outliers2.txt

output:
350
7983
24003
24006
24007
51866
68683
```

**Save outliers to new file**

```
$ mawk '!/#/' SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4.recode.vcf | cut -f1,2 > totalloci2
$ NUM=(`cat totalloci2 | wc -l`)
$ paste <(seq 1 $NUM) totalloci2 > loci2.plus.index
$ cat outliers2.txt | parallel "grep -w ^{} loci2.plus.index" | cut -f2,3> outlier2.loci.txt

$ head outlier2.loci.txt

output:
NC_035780.1	58431616
NC_035780.1	2492005
NC_035782.1	39657559
NC_035782.1	39657589
NC_035782.1	39657590
NC_035784.1	65251317
NC_035788.1	7098770
```

# 4. BayEnv2

> BayEnv2 is a Bayesian method that estimates the empirical pattern of covariance in allele frequencies between populations from a set of markers, and then uses this as a null model for a test at individual SNPs.

Steps to performing BayEnv2 follow the [bayenv2_manual.pdf](https://bitbucket.org/tguenther/bayenv2_public/src/default/bayenv2_manual.pdf).

##### _In Terminal_

_`bayenv2` needs to be downloaded to your working repository. Information for downloading `bayenv2` and other scripts can be accessed [here](https://bitbucket.org/tguenther/bayenv2_public/src/default/)._

**Convert VCF to BayEnv format**

```
# load in configuration file that will convert VCF to BayEnv format
$ curl -L -O https://raw.githubusercontent.com/amyzyck/EecSeq_NB_EasternOyster/master/Scripts/SNPBayEnv.spid
$ chmod +x SNPBayEnv.spid

$ java -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4.recode.vcf -outputfile SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4BayEnv.txt -spid SNPBayEnv.spid

output:
WARN  19:56:18 - PGDSpider configuration file not found! Loading default configuration.
initialize convert process...
read input file...
read input file done.
write output file...
WARN  19:56:59 - Locus SNP_1638 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_2031 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_2254 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_5194 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_6391 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_6868 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_6873 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_7272 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_8080 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_8082 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_8085 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_8189 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_8398 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_8442 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_8549 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_8564 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_8577 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_8578 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_9292 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_9294 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_9663 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_10145 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_10153 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_10210 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_10347 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_10455 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_10458 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_10607 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_10642 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_10643 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_11082 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_11083 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_13622 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_13744 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_13777 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_13958 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_14003 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15285 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15318 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15320 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15332 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15714 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15715 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15716 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15717 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15718 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15719 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15724 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15740 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15744 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15755 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15756 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15812 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15813 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15871 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15911 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15919 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_15998 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_16106 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_16176 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_16267 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_16674 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_17061 was removed, as it is not polymorphic
WARN  19:56:59 - Locus SNP_17119 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_18458 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_18846 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_19073 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_19250 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_19527 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_19528 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_19552 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_19554 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_19556 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_20056 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_20625 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_21503 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_23482 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_24618 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_25654 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_26037 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_26942 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_28640 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_28643 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_28644 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_28646 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_28649 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_28651 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_28652 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_28653 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_28654 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_29409 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_29923 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_31150 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_36351 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_37565 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_38409 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_38838 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_38908 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_39621 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_40156 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_40538 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_40550 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_40853 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_41048 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_41072 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_41117 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_41126 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_41162 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_41169 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_41273 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_41280 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_41299 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_41482 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_41733 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_41737 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_41739 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_41754 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_43886 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_44773 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_44813 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_46343 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_46358 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_46524 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_46615 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_47115 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_47549 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_47590 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_47933 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_48354 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_48443 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_48524 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_48793 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_48897 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_48910 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_48911 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_48966 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_48978 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_49182 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_49255 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_49256 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_49257 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_49262 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_49307 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_49460 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_49461 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_49864 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_49932 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_50002 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_50198 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_50310 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_50313 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_50321 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_50325 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_50466 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_50469 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_50656 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_50772 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_50844 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_50867 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_50875 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_51109 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_51461 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_51647 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_51680 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_51729 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_51888 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_51957 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_52311 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_53631 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_53771 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_54290 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_54835 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_54836 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_55237 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_55630 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_55631 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_55868 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_56233 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_56785 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_56926 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_56931 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_57103 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_58280 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_58320 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_58447 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_58455 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_58512 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_58571 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_58720 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_58740 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_59338 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_59390 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_59421 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_59431 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_59514 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_59619 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_59649 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_59692 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_59746 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_59766 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_60111 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_60166 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_60167 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_60315 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_60376 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_60377 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_60395 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_60397 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_60660 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_61049 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_62008 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_62752 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_63362 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_63367 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_63967 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_64012 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_64509 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_64697 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_64823 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_65607 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_66062 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_66916 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_66992 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_68294 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_68438 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_68447 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_68459 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_68522 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_69893 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_70012 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_70013 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_70066 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_70069 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_70123 was removed, as it is not polymorphic
WARN  19:57:00 - Locus SNP_70336 was removed, as it is not polymorphic
WARN  19:57:01 - Locus SNP_72890 was removed, as it is not polymorphic
WARN  19:57:01 - Locus SNP_72892 was removed, as it is not polymorphic
WARN  19:57:01 - Locus SNP_72991 was removed, as it is not polymorphic
WARN  19:57:01 - Locus SNP_73165 was removed, as it is not polymorphic
WARN  19:57:01 - Locus SNP_73282 was removed, as it is not polymorphic
WARN  19:57:01 - Locus SNP_73359 was removed, as it is not polymorphic
WARN  19:57:01 - Locus SNP_73562 was removed, as it is not polymorphic
WARN  19:57:01 - Locus SNP_74056 was removed, as it is not polymorphic
WARN  19:57:01 - Locus SNP_74863 was removed, as it is not polymorphic
WARN  19:57:01 - Locus SNP_74901 was removed, as it is not polymorphic
WARN  19:57:01 - Locus SNP_74965 was removed, as it is not polymorphic
write output file done.
```

**Obtaining a covariance matrix**

```
$ bayenv2 -i SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4BayEnv.txt -p 4 -k 100000 -r 63479 > matrix.out
```

**This code generates 100,000 iterations. We only need the last one.**

```
$ tail -5 matrix.out | head -4 > matrix
```

**With the matrix we will use our environmental factor file:**

```
$ cat environ

output:
-0.6111998032	0.1725740621	-0.89882324	1.337448981
-0.03630567357	-0.5445851035	-0.8350304921	1.415921269
1.470881402	-0.2569358707	-0.4772423125	-0.7367032193
-0.7833494518	0.2611164839	1.30558242	-0.7833494518
0.9044201628	0.5236116732	-1.380430775	-0.0476010612
1.375423951	-0.7406128967	0.1058018424	-0.7406128967
-0.6308364127	1.456931715	-0.6758961565	-0.1501991459
1.204749684	-0.5163212931	0.378635615	-1.067064006
```

The environmental file contains standardized environmental data with each line representing an environmental factor with the value for each population tab delimited.

* Row 1 - Latitude
* Row 2 - Longitude
* Row 3 - Sewage Effluent
* Row 4 - Temperature (°C)
* Row 5 - Salinity
* Row 6 - pH
* Row 7 - Chlorophyll- a (μg/L)
* Row 8 - Dissolved oxygen (mg/L)

**Calculate the Bayes Factor for each SNP for each environmental variable:**

```
$ curl -L -O https://bitbucket.org/tguenther/bayenv2_public/raw/2b2b7f20bb62fedaf5ea10b57e30bcb2807994b9/calc_bf.sh
$ chmod +x calc_bf.sh

$ calc_bf.sh SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4BayEnv.txt environ matrix 4 10000 8
```

The bayes factor for each SNP for each of the seven environmental variables is saved in `bf_environ.environ`.

**Convert `bf_environ.environ` into something suitable to input into R.**

```
# In column 1 replace the SNP names with numbers: 1-75152
# There were originally 75402 SNPs but 250 were removed when converting from VCF to BayEnv format because they were not polymorphic
$ paste <(seq 1 75156) <(cut -f2,3,4,5,6,7,8,9 bf_environ.environ ) > bayenv.out

# Add a header for each column
$ cat <(echo -e "Locus\tBF1\tBF2\tBF3\tBF4\tBF5\tBF6\tBF7\tBF8") bayenv.out > bayenv.final
```

##### _In RStudio_

**Idenitfy outlier SNPs (BayesFactor 1 value greater than 100).**
```
# BayesFactor 1
table_bay <- read.table("bayenv.final",header=TRUE)
plot(table_bay$BF1)

outliers1 <- table_bay[which(table_bay$BF1 > 100),]

outliers1
```

No outliers identified:

![BayF1](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/Outlier_Detection/BayF1.png)

This was repeated for BayesFactors 2-8.

**Compile outliers into one variable `outliers_total`.**

```
$ outliers <- rbind(outliers1,outliers2,outliers3,outliers4,outliers5,outliers6,outliers7)
$ outliers_total <- unique.data.frame(outliers)

$ outliers_total
```

| Locus |    BF1   |    BF2   |    BF3   |    BF4   |     BF5    |    BF6   |     BF7    |
|:-----:|:--------:|:--------:|:--------:|:--------:|:----------:|:--------:|:----------:|
| 14964 |  0.42959 | 0.045923 | 0.045564 | 0.043329 | 1.6724e+02 | 7.737900 | 4.4629e+01 |
| 17572 |   2.9279 | 0.044583 | 0.041009 | 0.045147 | 3.1868e+03 | 3.617000 | 1.0896e+03 |
| 19987 |  0.15004 |  0.13188 |   4.5686 |   334.63 | 1.4725e-01 |  0.10059 | 1.1171e-01 |
| 20899 |  0.21218 |  0.73971 |   110.56 |   9.6912 | 7.2967e-02 | 0.054825 | 6.7746e-02 |
| 33279 | 0.065375 |   0.7632 |   133.70 |   5.5882 | 3.8977e-02 | 0.037459 | 3.2558e-02 |
| 38968 |  0.10745 |  0.16364 |   98.426 |   179.51 | 4.8184e-02 | 0.052549 | 5.4444e-02 |
| 42549 |   2.2624 | 0.037646 | 0.021242 | 0.022797 | 6.0928e+02 | 0.049653 | 5.3130e+02 |
| 50789 |   0.1222 |   1.1601 |   495.04 |   19.272 | 6.6858e-02 | 0.058311 | 4.7168e-02 |
| 54194 |  0.05102 |  0.13102 | 0.082396 | 0.048275 | 3.5243e-01 |   167.97 | 7.6152e-02 |
| 58383 |   6.8223 |   21.862 |   155.34 |   1.5474 | 4.1303e-02 | 0.046958 | 5.8122e-02 |
| 69248 | 0.025009 | 0.022183 |   1.8039 |   345.83 | 3.3897e-02 | 0.033261 | 2.8623e-02 |

**Save outliers to .txt file**

```
# Save only the values from column 1 and exclude the header_line
$ total_outliers <- outliers_total[,1, drop=TRUE]

$ invisible(lapply(total_outliers, write, "outliers3.txt", append=TRUE))
```

##### _In Terminal_

**Check out outliers3.txt**

```
$ head outliers3.txt

output:
20899
33279
50789
58383
19987
38968
69248
14964
17572
42549
```

**Save outliers to new file**

```
$ mawk '!/#/' SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4.recode.vcf | cut -f1,2 > totalloci3
$ NUM=(`cat totalloci3 | wc -l`)
$ paste <(seq 1 $NUM) totalloci3 > loci3.plus.index
$ cat outliers3.txt | parallel "grep -w ^{} loci3.plus.index" | cut -f2,3> outlier3.loci.txt

$ head outlier3.loci.txt

output:
NC_035782.1	3811652
NC_035783.1	18475147
NC_035784.1	59862728
NC_035785.1	32466402
NC_035781.1	57449356
NC_035783.1	56550670
NC_035788.1	19035433
NC_035781.1	26985674
NC_035781.1	41554587
NC_035784.1	16551538
```

# 5. Combine all outlier loci into one file

```
$ cat outlier*.loci.txt > all.outliers
$ cut -f1 all.outliers | sort | uniq | wc -l

output:
9
```

# 6. Separating outlier and neutral loci into 2 VCF files

**Create VCF file with just _neutral_ loci**

```
$ vcftools --vcf SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4.recode.vcf --exclude-positions all.outliers --recode-INFO-all --out neutralloci --recode

output:
VCFtools - 0.1.17
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4.recode.vcf
	--exclude-positions all.outliers
	--recode-INFO-all
	--out neutralloci
	--recode
After filtering, kept 40 out of 40 Individuals
Outputting VCF file...
After filtering, kept 74520 out of a possible 75402 Sites
Run Time = 15.00 seconds
```

_74520 neutral loci in `neutralloci.recode.vcf`._

**Create VCF file with just _outlier_ loci**

```
$ vcftools --vcf SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4.recode.vcf --recode --recode-INFO-all --positions all.outliers --out outlierloci

output:
VCFtools - 0.1.17
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4.recode.vcf
	--recode-INFO-all
	--out outlierloci
	--positions all.outliers
	--recode
After filtering, kept 40 out of 40 Individuals
Outputting VCF file...
After filtering, kept 882 out of a possible 75402 Sites
Run Time = 1.00 seconds
```

_882 outlier loci in `outlierloci.recode.vcf`._

#### Haplotig masked genome

_90861 neutral loci in `neutralloci.recode.vcf`._

_841 outlier loci in `outlierloci.recode.vcf`._

***

**Three VCF files will be used to perform seascape genomics analyses:**

* All SNPs (oultier and neutral combined): `SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A4.recode.vcf`
* Neutral SNPs only: `neutralloci.recode.vcf`
* Outlier SNPs only: `outlierloci.recode.vcf`

**Haplotig Masked files:**

* All SNPs (oultier and neutral combined): `SNP.TRSdp5g5mafMIp9g9dnDNAmaf052Ahap.recode.vcf`
* Neutral SNPs only: `neutrallocihap.recode.vcf`
* Outlier SNPs only: `outlierlocihap.recode.vcf`

**_Note_**: I chose to include all outliers identified across all 4 outlier programs. I am also planning to use Latent Factors Mixed Models (LFMM) as an additional genotype-environment association analysis for identifying SNPs putatively under selection based on associations with environmental variables. The documentation for this analysis will be added upon completion.

Population genomic and seascape genomic analyses were performed on all 3 SNP datasets (+ Haplotig Masked SNP datasets) and documented [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/tree/master/Analysis/PopGen_SeaGen_Analyses). 
