# *Crassostrea virginica* EecSeq Seascape Genomics Re-Analysis - dDocent pipeline

Author: Amy Zyck

Last Edited: February 19, 2023

Initial analysis was performed in March 2020 on sequence files for 4 different oyster populations in Narragansett Bay. Following results of this initial analysis, additional sampling of adult oysters was completed Fall 2020 for 4 additional populations (and one repeated population - NAR). Mantle tissue was extracted from the 10 adult oysters per population and then processed using EecSeq.

The analysis steps performed in March 2020 are documented [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/tree/master/Analysis/Analysis_Part1).

Data uploaded and analyzed on KITT. User logged in before following steps are completed.

**Data location: `PATH /home/azyck/NB_capture_both/NB_ddhaplo`**

## In Terminal:

### 1. Setup: Downloading software programs, creating environments anaconda folders

**Downloading Bioconda (skip this step if Bioconda has already been downloaded).**

```javascript
# downloading Miniconda software
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ chmod +x Miniconda3-latest-Linux-x86_64.sh
$ ./Miniconda3-latest-Linux-x86_64.sh

# restarting window with source command
$ source ~/.bashrc

# adding different channels
$ conda config --add channels defaults
	# should get: "warning: 'defaults' already in 'channels' list, moving to the top"
$ conda config --add channels bioconda
$ conda config --add channels conda-forge

# to see if this worked with cat command
$ cat .condarc
```
> Bioconda is a bioinformatics software package manager. See more at [https://bioconda.github.io](https:///bioconda.github.io).

**Create and activate a dDocent conda environment:**

```
$ conda create -n nb_capture ddocent
$ conda activate nb_capture

# the beginning of the line should then look like: (nb_capture) [username@KITT ~]$
```

**Creating and entering a directory for this project**

```
$ mkdir NB_capture_both # creates a directory
$ cd NB_capture_both # navigates into this new directory
# the beginning of the line should now look like: (nb_capture) [username@KITT NB_capture_both]$
```

**Link raw data to `NB_capture_both` directory**

The raw data was uploaded to KITT by [J. Puritz](https://github.com/jpuritz) and linked to the working directory.

```
$ ln -s /RAID_STORAGE2/azyck/NB_capture /home/azyck/NB_capture_both

# NB_capture directory (containing all files) within RAID_STORAGE2 will be linked to the azyck home directory in KITT
```

### 2-3. FastQC Analysis & Demultiplexing Data Files

[J. Puritz](https://github.com/jpuritz) performed the FastQC analysis and demultiplexing steps for this analysis. Depending on how barcodes are used in your sequencing protocol, you may have barcodes on both the forward and paired-end reads. Make sure that both the forward (.F.fq.gz) and reverse (.R.fq.gz) raw sequence files are fully demultiplexed. Otherwise, "...it will potentially introduce some level of spurious SNPs that would make samples with the same barcode appear artificially similar." - J. Puritz

You can check the demultiplexed files for lingering barcodes:
```
zcat PVD_9.R.fq.gz | mawk 'NR%4==2' | head -1000 | cut -c1-6|sort | uniq -c |sort -k1 -r | head
```

If the paired end read was not demultiplexed properly, there will be a large number of the same sequence like this:
```
   980 ATCAGT
     5 ACCAGT
     3 TCAGTT
     3 ATCAGA
     2 TCAGTC
     2 GTCAGT
     1 TTCAGT
     1 ATCGGT
     1 ATCCGT
     1 ATCAGG
```

All steps to complete the FastQC analysis and demultiplexing can be accessed [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/Analysis_Part1/EecSeq_Cvirginica_dDocent.md).

**Checking the quality of data post-demultiplexing**

```
$ mkdir demultiplexed_fastqc_results
$ cd demultiplexed_fastqc_results
$ fastqc ../*fq.gz
$ mv *fastqc.* demultiplexed_fastqc_results
```

**Multiqc Analysis**

```
$ multiqc .

$ mv multiqc_report.html demultiplexed_multiqc_report.html

# copy multiqc report to directory outside of KITT
# run the following commands outside of KITT
$ scp -r -P XXXX x@kitt.uri.edu:/home/azyck/NB_capture_both/NB_ddhaplo/demultiplexed_fastqc_results/demultiplexed_multiqc_report.html /Users/azyck/Documents/Research/NB_AdultOyster/SeaGen_Analysis/fastqc
# XXXX above indicates the password specific for our lab's KITT.
# scp = secure copy, -r = indicates copying a folder/repository

# open multiqc_report to view the quality of raw, multiplexed data files
```

Multiqc report can be viewed [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/multiqc/demultiplexed_multiqc_report_Analysis2.html).

_Some personal notes: There appears to be a large difference in the % duplicate reads between the forward and paired end reads of the newer five populations. The % duplicate reads of the paired end reads is about 1/2 the value of the forward reads. This is not the case for reads of the original four populations. After looking into this, it's possible this is an artifact of caused by sequencing errors. The reverse reads are still high quality (although slightly lower quality than the forward reads). I will move forward with read trimming, but will keep an eye on these samples._

_Additional note: Samples MCD_2 and BAR_6 are duplicated, resulting in an "original" version and a second version with _501 appended to the sample name (both versions appear in the multiqc report). In Part 1 of the analysis with the original four populations, samples with index primer set 501/701 had much lower total sequences than the other samples. To test any issues with this primer set, these two samples in the second sample set with the newest five populations were duplicated, so that one version was prepared with index primer set 501/701 and the other version was prepared with  a different index primer set. From the demultiplexed multiqc report, BAR_6 had a much higher total sequence count compared to BAR_6_501, whereas MCD_2_501 had a higher total sequence count compared to MCD_2. There does not appear to be an issue with the 501/701 index primer set. Therefore, I am going to separate BAR_6_501 and MCD_2 from the total dataset and move forward without them._

```
# Create a new directory to save duplicated sample files

$ mkdir duplicated_samples
$ cd duplicated_samples

# Move unwanted sample versions from NB_capture_both to duplicated_samples directory

$ mv ../BAR_6_501.* .
$ mv ../MCD_2.* .
```

### 4. Read Trimming, Mapping, and SNP Calling

Using [**dDocent**](http://www.ddocent.com/)

> dDocent is a bioinformatics program created by Dr. Jon Purtiz that is specifically designed for different types of RAD sequencing.

**Jon downloaded a version of dDocent on my KITT account `dDocent_ngs` that can be used for Expressed Exome Capture Sequencing (EecSeq). It is accessible [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Scripts/dDocent_ngs).

If this is your first time running dDocent, I recommend going through the [Quick Start Guide](http://www.ddocent.com/quick/). I also recommend:

1. Reading through the [User Guide](http://www.ddocent.com/UserGuide/).
2. Completing the [Assembly Tutorial](http://www.ddocent.com/assembly/), using the simulated dataset.

**Create and activate a dDocent conda environment (if you did not do so previously):**

```
$ conda create -n nb_capture ddocent
$ conda activate nb_capture

# the beginning of the line should then look like: (nb_capture) [azyck@KITT ~]$
```

**Make directory for dDocent and link files and `dDocent_ngs` into directory**

Starting in `NB_capture_both` directory:

```
$ mkdir NB_ddhaplo
$ cd NB_ddhaplo/

$ ln -s ../*.fq.gz .
$ ln -s ../dDocent_ngs .
```

**First run `dDocent_ngs` with just Read Trimming**
```
$ bash dDocent_ngs
```

```
dDocent 3.0.0

Contact jpuritz@uri.edu with any problems


Checking for required software

All required software is installed!

dDocent version 3.0.0 started Fri Feb 10 09:41:53 EST 2023

90 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
$ yes
Proceeding with 50 individuals
dDocent detects 80 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
$ 20

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
$ yes

Do you want to perform an assembly?
Type yes or no and press [ENTER].
$ no

Reference contigs need to be in a file named reference.fasta

Do you want to map reads?  Type yes or no and press [ENTER]
$ no
Mapping will not be performed
Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
$ no

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.


At this point, all configuration information has been entered and dDocent may take several hours to run.
It is recommended that you move this script to a background operation and disable terminal input and output.
All data and logfiles will still be recorded.
To do this:
Press control and Z simultaneously
Type bg and press enter
Type disown -h and press enter
Now sit back, relax, and wait for your analysis to finish

Trimming reads

dDocent has finished with an analysis in /home/azyck/NB_capture_both/NB_ddhaplo

dDocent started Fri Feb 10 09:41:53 EST 2023

dDocent finished Fri Feb 10 09:58:09 EST 2023
```

**Before continuing with dDocent, I am checking the quality of the trimmed reads with FastQC**

Back in `NB_capture_both/NB_ddhaplo` Directory

```
$ mkdir cleanedreads_fastqc_results
$ cd cleanedreads_fastqc_results/

$ ln -s ../*.fq.gz .
```

```
# fastqc takes awhile to run on larger data sets, so I recommend running it overnight
# you may be able to run fastqc on the .R1.fq.gz and .R2.fq.gz files simultaneously using two different terminal windows

$ fastqc *.R1.fq.gz
$ fastqc *.R2.fq.gz
```

```
$ multiqc .

$ mv multiqc_report.html cleanedreads_multiqc_report_Analysis2.html

$ scp -r -P XXXX x@kitt.uri.edu:/home/azyck/NB_capture_both/NB_ddhaplo/cleanedreads_fastqc_results/cleanedreads_multiqc_report.html /Users/azyck/Documents/Research/NB_AdultOyster/SeaGen_Analysis/fastqc
```

Multiqc report can be viewed [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/multiqc/cleanedreads_multiqc_report_Analysis2.html).

The difference in % duplicate reads between the .R1 and .R2 files for the four newer populations (BAR, KIC, MCD, and NAR2) is smaller after read trimming. The sequence quality is still good and adapter content is much better (lower percentage of sequences).

**Download _Crassostrea virginica_ genome**

Jon added the Eastern Oyster genome to the `NB_ddhaplo` Directory for me. Resources for the Eastern Oyster genome are on [NCBI](https://www.ncbi.nlm.nih.gov/Traces/wgs/MWPT03?display=contigs&page=1).

```
# In NB_ddhaplo directory
$ wget #paste link here

# Reference contigs need to be in a file named reference.fasta
$ mv #file_name_here reference.fasta
```

**Read Mapping**

Read mapping has to be performed at least once before SNP calling. If trimming and assembly remains the same when using dDocent repeatedly, then read mapping is not necessary. A (match score), B (mismatch score), and 0 (gap penalty) parameters are used.

*Match Score*: score for a matching base.

*Mismatch Score*: score for a mismatching base.

*Gap Penalty*: Mutations, caused by either insertions or deletions, are annotated as gaps in the sequence and refered to as indels. These gaps are represented as dashes in the alignment. Optimization of gap penalty allows for the most accurate alignment possible. This will avoid low scores in alignments. Too small of a gap penalty will prevent a high level of alignment, and too large of a gap penalty will prevent accurate alignments.

I'm going to start with A = 1, B = 3, and O = 5 for each parameter. These were the parameter values used in my first go of the analysis with the initial 5 populations.

```
$ bash dDocent_ngs
```

```
dDocent 3.0.0

Contact jpuritz@uri.edu with any problems


Checking for required software

All required software is installed!

dDocent version 3.0.0 started Fri Feb 10 17:26:58 EST 2023

90 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes
Proceeding with 90 individuals
dDocent detects 80 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
20

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
no

Do you want to perform an assembly?
Type yes or no and press [ENTER].
no

Reference contigs need to be in a file named reference.fasta

Do you want to map reads?  Type yes or no and press [ENTER]
yes
BWA will be used to map reads.  You may need to adjust -A -B and -O parameters for your taxa.
Would you like to enter a new parameters now? Type yes or no and press [ENTER]
yes
Please enter new value for A (match score).  It should be an integer.  Default is 1.
1
Please enter new value for B (mismatch score).  It should be an integer.  Default is 4.
3
Please enter new value for O (gap penalty).  It should be an integer.  Default is 6.
5
Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
no

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.

At this point, all configuration information has been entered and dDocent may take several hours to run.
It is recommended that you move this script to a background operation and disable terminal input and output.
All data and logfiles will still be recorded.
To do this:
Press control and Z simultaneously
Type bg and press enter
Type disown -h and press enter
Now sit back, relax, and wait for your analysis to finish

Using BWA to map reads

Creating alignment intervals

dDocent has finished with an analysis in /home/azyck/NB_capture_both/NB_ddhaplo

dDocent started Fri Feb 10 17:26:58 EST 2023

dDocent finished Sat Feb 11 07:01:10 EST 2023
```

All individuals bam files are merged into a single bam file for easier parallelized variant calling, under name `filter.merged.bam`.

Use `samtools flagstat` to investigate `bam` files.

https://github.com/bahlolab/bioinfotools/blob/master/SAMtools/flagstat.md

```
$ samtools flagstat KIC_9.F.bam

output:
5184643 + 0 in total (QC-passed reads + QC-failed reads)
5184643 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
5184643 + 0 mapped (100.00% : N/A)
5184643 + 0 primary mapped (100.00% : N/A)
5184643 + 0 paired in sequencing
2598389 + 0 read1
2586254 + 0 read2
5047683 + 0 properly paired (97.36% : N/A)
5175950 + 0 with itself and mate mapped
8693 + 0 singletons (0.17% : N/A)
84955 + 0 with mate mapped to a different chr
84955 + 0 with mate mapped to a different chr (mapQ>=5)
```

```
$ samtools flagstat GB_6.F.bam

output:
14054285 + 0 in total (QC-passed reads + QC-failed reads)
14054285 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
14054285 + 0 mapped (100.00% : N/A)
14054285 + 0 primary mapped (100.00% : N/A)
14054285 + 0 paired in sequencing
7028938 + 0 read1
7025347 + 0 read2
13732780 + 0 properly paired (97.71% : N/A)
14037815 + 0 with itself and mate mapped
16470 + 0 singletons (0.12% : N/A)
201853 + 0 with mate mapped to a different chr
201853 + 0 with mate mapped to a different chr (mapQ>=5)
```

I looked at a few other `bam` files. All had 100% mapped. A low percentage of the mappings were singletons (only one read from a pair). Percentage of proper pairings was between 97 and 99%.

I want to test out mapping parameter values before moving forward. Rather than running the dDocent mapping on all the samples, I'm going to save a subset of individuals into a new directory and perform the mapping test on that subset.

```
$ mkdir mapping_test

$ cd mapping_test

$ ln -s ../BAR_2.* .
# repeat the linking command for one or two individuals from each population

# link reference.fasta file
$ ln -s ../reference.fasta .
```

List of samples used in this test:

```
BAR_2 BIS_5
GB_7  GHP_8
KIC_5 MCD_3
NAR_1 NAR2_8
PVD_10
```

Test 1: Parameter values A = 1, B = 3, and O = 5
```
Output:

GB_7
16263343 + 0 in total (QC-passed reads + QC-failed reads)
16263343 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
16263343 + 0 mapped (100.00% : N/A)
16263343 + 0 primary mapped (100.00% : N/A)
16263343 + 0 paired in sequencing
8134962 + 0 read1
8128381 + 0 read2
16019655 + 0 properly paired (98.50% : N/A)
16253396 + 0 with itself and mate mapped
9947 + 0 singletons (0.06% : N/A)
146384 + 0 with mate mapped to a different chr
146384 + 0 with mate mapped to a different chr (mapQ>=5)

NAR2_8
15236395 + 0 in total (QC-passed reads + QC-failed reads)
15236395 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
15236395 + 0 mapped (100.00% : N/A)
15236395 + 0 primary mapped (100.00% : N/A)
15236395 + 0 paired in sequencing
7634776 + 0 read1
7601619 + 0 read2
14901247 + 0 properly paired (97.80% : N/A)
15213764 + 0 with itself and mate mapped
22631 + 0 singletons (0.15% : N/A)
182583 + 0 with mate mapped to a different chr
182583 + 0 with mate mapped to a different chr (mapQ>=5)
```

Test 2: Parameter values A = 2, B = 4, and O = 6

```
Output:

GB_7
16294659 + 0 in total (QC-passed reads + QC-failed reads)
16294659 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
16294659 + 0 mapped (100.00% : N/A)
16294659 + 0 primary mapped (100.00% : N/A)
16294659 + 0 paired in sequencing
8150048 + 0 read1
8144611 + 0 read2
16064791 + 0 properly paired (98.59% : N/A)
16286784 + 0 with itself and mate mapped
7875 + 0 singletons (0.05% : N/A)
141407 + 0 with mate mapped to a different chr
141407 + 0 with mate mapped to a different chr (mapQ>=5)

NAR2_8
15251787 + 0 in total (QC-passed reads + QC-failed reads)
15251787 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
15251787 + 0 mapped (100.00% : N/A)
15251787 + 0 primary mapped (100.00% : N/A)
15251787 + 0 paired in sequencing
7635293 + 0 read1
7616494 + 0 read2
14941204 + 0 properly paired (97.96% : N/A)
15234746 + 0 with itself and mate mapped
17041 + 0 singletons (0.11% : N/A)
173844 + 0 with mate mapped to a different chr
173844 + 0 with mate mapped to a different chr (mapQ>=5)
```

These parameter values (A = 2, B = 4, and O = 6) resulted in more passed reads, more mapped reads, higher perceentage of properly paired reads, and a lower percentage of singletons.

Test 3: Default parameter values (A = 1, B = 4, and O = 6)

```
Output:

GB_7
16192461 + 0 in total (QC-passed reads + QC-failed reads)
16192461 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
16192461 + 0 mapped (100.00% : N/A)
16192461 + 0 primary mapped (100.00% : N/A)
16192461 + 0 paired in sequencing
8100099 + 0 read1
8092362 + 0 read2
15938546 + 0 properly paired (98.43% : N/A)
16180887 + 0 with itself and mate mapped
11574 + 0 singletons (0.07% : N/A)
150431 + 0 with mate mapped to a different chr
150431 + 0 with mate mapped to a different chr (mapQ>=5)

NAR2_8
15143828 + 0 in total (QC-passed reads + QC-failed reads)
15143828 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
15143828 + 0 mapped (100.00% : N/A)
15143828 + 0 primary mapped (100.00% : N/A)
15143828 + 0 paired in sequencing
7600775 + 0 read1
7543053 + 0 read2
14790527 + 0 properly paired (97.67% : N/A)
15114876 + 0 with itself and mate mapped
28952 + 0 singletons (0.19% : N/A)
188784 + 0 with mate mapped to a different chr
188784 + 0 with mate mapped to a different chr (mapQ>=5)
```

The outputs from these parameter values (A = 1, B = 4, and O = 6) are slightly worse than the previous tests. I am going to stick with a higher A value (match score). I'm going to test higher values than the default for all three parameters next.

Test 4: Parameter values A = 3, B = 5, and O = 7

```
Output:

GB_7
16282725 + 0 in total (QC-passed reads + QC-failed reads)
16282725 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
16282725 + 0 mapped (100.00% : N/A)
16282725 + 0 primary mapped (100.00% : N/A)
16282725 + 0 paired in sequencing
8143714 + 0 read1
8139011 + 0 read2
16061927 + 0 properly paired (98.64% : N/A)
16275461 + 0 with itself and mate mapped
7264 + 0 singletons (0.04% : N/A)
136750 + 0 with mate mapped to a different chr
136750 + 0 with mate mapped to a different chr (mapQ>=5)

NAR2_8
15229350 + 0 in total (QC-passed reads + QC-failed reads)
15229350 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
15229350 + 0 mapped (100.00% : N/A)
15229350 + 0 primary mapped (100.00% : N/A)
15229350 + 0 paired in sequencing
7622615 + 0 read1
7606735 + 0 read2
14934542 + 0 properly paired (98.06% : N/A)
15214161 + 0 with itself and mate mapped
15189 + 0 singletons (0.10% : N/A)
165681 + 0 with mate mapped to a different chr
165681 + 0 with mate mapped to a different chr (mapQ>=5)
```

Compared to parameter values for test 2 (A = 2, B = 4, and O = 6), this set of parameter values (A = 3, B = 5, and O = 7) resulted in a lower number of reads, but a higher percentage of properly paired reads and lower percentage of singletons.

Test 5: Parameter values A = 2, B = 3, and O = 5

```
Output:

GB_7
16245770 + 0 in total (QC-passed reads + QC-failed reads)
16245770 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
16245770 + 0 mapped (100.00% : N/A)
16245770 + 0 primary mapped (100.00% : N/A)
16245770 + 0 paired in sequencing
8125213 + 0 read1
8120557 + 0 read2
16027463 + 0 properly paired (98.66% : N/A)
16238843 + 0 with itself and mate mapped
6927 + 0 singletons (0.04% : N/A)
135402 + 0 with mate mapped to a different chr
135402 + 0 with mate mapped to a different chr (mapQ>=5)

NAR2_8
15176941 + 0 in total (QC-passed reads + QC-failed reads)
15176941 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
15176941 + 0 mapped (100.00% : N/A)
15176941 + 0 primary mapped (100.00% : N/A)
15176941 + 0 paired in sequencing
7596249 + 0 read1
7580692 + 0 read2
14885742 + 0 properly paired (98.08% : N/A)
15162587 + 0 with itself and mate mapped
14354 + 0 singletons (0.09% : N/A)
163891 + 0 with mate mapped to a different chr
163891 + 0 with mate mapped to a different chr (mapQ>=5)
```

**Based on these five tests, I am going to redo the read mapping on the full dataset using parameter values A = 2, B = 4, and O = 6. With these values, passed reads were maximized with a high percentage of properly paired reads and a lower percentage of singletons.**

Return back to directory with all samples (`NB_ddhaplo`) and rerun `dDocent_ngs`.

```
$ cd ../

$ bash dDocent_ngs
```

```
dDocent 3.0.0

Contact jpuritz@uri.edu with any problems


Checking for required software

All required software is installed!

dDocent version 3.0.0 started Tue Feb 14 19:02:55 EST 2023

90 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes
Proceeding with 90 individuals
dDocent detects 80 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
20

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
no

Do you want to perform an assembly?
Type yes or no and press [ENTER].
no

Reference contigs need to be in a file named reference.fasta

Do you want to map reads?  Type yes or no and press [ENTER]
yes
BWA will be used to map reads.  You may need to adjust -A -B and -O parameters for your taxa.
Would you like to enter a new parameters now? Type yes or no and press [ENTER]
yes
Please enter new value for A (match score).  It should be an integer.  Default is 1.
2
Please enter new value for B (mismatch score).  It should be an integer.  Default is 4.
4
Please enter new value for O (gap penalty).  It should be an integer.  Default is 6.
6
Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
no

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.

At this point, all configuration information has been entered and dDocent may take several hours to run.
It is recommended that you move this script to a background operation and disable terminal input and output.
All data and logfiles will still be recorded.
To do this:
Press control and Z simultaneously
Type bg and press enter
Type disown -h and press enter
Now sit back, relax, and wait for your analysis to finish

Using BWA to map reads

Creating alignment intervals

dDocent has finished with an analysis in /home/azyck/NB_capture_both/NB_ddhaplo

dDocent started Tue Feb 14 19:02:55 EST 2023

dDocent finished Wed Feb 15 09:33:49 EST 2023



dDocent 3.0.0
The 'd' is silent, hillbilly.
```

Use `samtools flagstat` to investigate `bam` files.

```
$ samtools flagstat BAR_3.F.bam

14219450 + 0 in total (QC-passed reads + QC-failed reads)
14219450 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
14219450 + 0 mapped (100.00% : N/A)
14219450 + 0 primary mapped (100.00% : N/A)
14219450 + 0 paired in sequencing
7116509 + 0 read1
7102941 + 0 read2
13977879 + 0 properly paired (98.30% : N/A)
14208692 + 0 with itself and mate mapped
10758 + 0 singletons (0.08% : N/A)
144313 + 0 with mate mapped to a different chr
144313 + 0 with mate mapped to a different chr (mapQ>=5)
```

```
$ samtools flagstat KIC_9.F.bam

5191343 + 0 in total (QC-passed reads + QC-failed reads)
5191343 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
5191343 + 0 mapped (100.00% : N/A)
5191343 + 0 primary mapped (100.00% : N/A)
5191343 + 0 paired in sequencing
2599310 + 0 read1
2592033 + 0 read2
5062121 + 0 properly paired (97.51% : N/A)
5184920 + 0 with itself and mate mapped
6423 + 0 singletons (0.12% : N/A)
82344 + 0 with mate mapped to a different chr
82344 + 0 with mate mapped to a different chr (mapQ>=5)
```

Overall, a larger number of passed reads, a higher percentage of properly paired reads, and a lower percentage of singletons compared to the first mapping attempt on the full dataset with parameters (A = 1, B = 3, O = 5).

**Variant Calling**

Run dDocent again to generate the variant calls and then zip the huge VCF file!

```
$ bash dDocent_ngs
```

```
dDocent 3.0.0

Contact jpuritz@uri.edu with any problems


Checking for required software

All required software is installed!

dDocent version 3.0.0 started Wed Feb 15 12:35:51 EST 2023

90 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes
Proceeding with 90 individuals
dDocent detects 80 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
20

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
no

Do you want to perform an assembly?
Type yes or no and press [ENTER].
no

Reference contigs need to be in a file named reference.fasta

Do you want to map reads?  Type yes or no and press [ENTER]
no
Mapping will not be performed
Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
yes

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.

At this point, all configuration information has been entered and dDocent may take several hours to run.
It is recommended that you move this script to a background operation and disable terminal input and output.
All data and logfiles will still be recorded.
To do this:
Press control and Z simultaneously
Type bg and press enter
Type disown -h and press enter
Now sit back, relax, and wait for your analysis to finish

Using FreeBayes to call SNPs

Checking the headers of 7 files.
Done, the headers are compatible.
Concatenating Collated.raw.0.bcf[E::hts_open_format] Failed to open file "TotalRawSNPs.vcf.gz" : Permission denied
TotalRawSNPs.vcf.gz: Permission denied

Using VCFtools to parse TotalRawSNPS.vcf for SNPs that are called in at least 90% of individuals

After filtering, kept 90 out of 90 Individuals
After filtering, kept 6585214 out of a possible 55253967 Sites
```

> FreeBayes detects variants based on differences in haplotypes, including single nucleotide polymorphisms (SNPs), and indels (insertions and deletions). [FreeBayes](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionVI/lessons/02_variant-calling.html).


__`dDocent_ngs` variant calling results in a g-zipped VCF file.__


Filtering of the VCF file was then performed and documented [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/Analysis_Part2/EecSeq_Cvirginica_Filtering_Analysis_Part2.md).
