# Narragansett Bay *Crassostrea virginica* EecSeq Seascape Genomics Final Analysis - Full Bioinformatic Analysis 

Author: Amy Zyck

Last Edited: September 12, 2024

This markdown includes documentation of the finalized analysis for the Narragansett Bay *C.virginica* seascape genomics project. 

Initial analyses was performed in March 2020 on sequence files for 4 different oyster populations in Narragansett Bay. Following results of this initial analysis, additional sampling of adult oysters was completed Fall 2020 for 4 additional populations (and one repeated population - NAR). Mantle tissue was extracted from the 10 adult oysters per population and then processed using EecSeq. Sequence data was processed and analyzed in February of 2023. The analysis steps performed in March 2020 are documented [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/tree/master/Analysis/Analysis_Part1) and the steps performed in February 2023 are documented [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/tree/master/Analysis/Analysis_Part2). 

__During the second iteration of the analysis, it was discovered that two of the sample sites likely had non-wild oysters introduced at some point in time (NAR and GHP). Because they are not truly wild oysters, they were removed from the analysis.The discovery of the non-wild oysters is documented [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/Analysis_Part2/NB_OutlierDetection_6pops.Rmd).__ 

__In May 2024, 40 individuals with the lowest amount of sequence data were resequenced to increase the number of reads in those samples. This markdown documents the demultiplexing of the resequenced samples followed by the complete bioinformatic analysis of the full sample set (with NAR and GHP excluded).__ 

## Sequence data prep

**Data location: PATH: `/RAID_STORAGE2/Raw_Data/NB_and_EAGER_2024/01.RawData/`**

### In Terminal:

Data uploaded and analyzed on KITT. User logged in before following steps are completed.

Within `NB_capture_both/NB_ddhaplo` I made new directory `NB_ReSeq_Files`  and linked the multiplexed file folders from the `RAID_STORAGE2` directory to this new directory: 

```
$ ln -s /RAID_STORAGE2/Raw_Data/NB_and_EAGER_2024/01.RawData/D* .
```

Activate conda environment ([steps for setting up a conda environment](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/Analysis_Part2/EecSeq_Cvirginica_dDocent_Part2.md))

```
$ conda activate nb_capture4
```

Make new directory for multiplexed multiqc analysis
```
$ mkdir multiplexed_fastqc_results
$ cd multiplexed_fastqc_results
```

Copied all .fq.gz files to this directory 

```
$ cp ../D501_cap1/*.fq.gz . 
```

Run fastqc 

```
$ fastqc *.fq.gz
```

Create multiqc report
```
$ multiqc .
```

Renaming multiqc report for organization
```
$ mv multiqc_report.html reseq_multiplexed_multiqc_report.html 
```

Multiqc report can be viewed [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/final_analysis_outputs/multiqc_reports/reseq_multiplexed_multiqc_report.html).

## Demultiplexing 

Make a .txt file for each set of barcodes 

```
$ nano Library501cap13_barcodes.txt 

PVD_1	ATCACG	ATCACG
PVD_5	ACAGTG	ACAGTG
PVD_3	TTAGGC	TTAGGC
```

```
# remove first column and just save barcodes
$ cut -f2,3 Library501cap13_barcodes.txt > Library501cap13
```

Create new directory for each library group 

```
$ mkdir Library501cap13
```

Move barcodes .txt file to directory for each library 

```
$ mv Library501cap13_barcodes.txt Library501cap13/
```

In main directory, save custom script for renaming files for dDocent

```
$ curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/Rename_for_dDocent.sh
$ chmod +x Rename_for_dDocent.sh
```

Run process_shortreads for each library 

```
$ process_shortreads -1 D501_cap1_CKDL240018752-1A_H3VJJDSXC_L3_1.fq.gz -2 D501_cap1_CKDL240018752-1A_H3VJJDSXC_L3_2.fq.gz -b barcodes501cap13 -r -i gzfastq -o /home/azyck/NB_capture_both/NB_ddhaplo/NB_ReSeq_Files/Library501cap13/ --inline-inline -c -q -D

output: 
Using Phred+33 encoding for quality scores.
Reads trimmed shorter than 31 nucleotides will be discarded.
Found 1 paired input file(s).
Searching for single and paired-end, inlined barcodes.
Loaded 3 barcodes (6bp / 6bp).
Will attempt to recover barcodes with at most 1 / 1 mismatches.
Processing file 1 of 1 [D501_cap1_CKDL240018752-1A_H3VJJDSXC_L3_1.fq.gz]
  Reading data from:
  D501_cap1_CKDL240018752-1A_H3VJJDSXC_L3_1.fq.gz and
  D501_cap1_CKDL240018752-1A_H3VJJDSXC_L3_2.fq.gz
  6326904 total reads; -1056184 ambiguous barcodes; +62542 recovered; -1169 low quality reads; 5269551 retained reads.
    0 trimmed reads; 0 orphaned paired-ends.
Closing files, flushing buffers...
Outputing details to log: '/home/azyck/NB_capture_both/NB_ddhaplo/NB_ReSeq_Files/Library501cap13/process_shortreads.log'

6326904 total sequences;
  1056184 ambiguous barcode drops;
  1169 low quality read drops;
  0 trimmed reads;
  0 orphaned paired-end reads;
5269551 retained reads.
```

You can get the `process_shortreads.log` file - saved in each library-specific directory

```
$ head -n 50 process_shortreads.log

output: 
process_shortreads v2.62, executed 2024-06-25 12:00:24 (zlib-1.2.7)
process_shortreads -1 D501_cap1_CKDL240018752-1A_H3VJJDSXC_L3_1.fq.gz -2 D501_cap1_CKDL240018752-1A_H3VJJDSXC_L3_2.fq.gz -b barcodes501cap13 -r -i gzfastq -o /home/azyck/NB_capture_both/NB_ddhaplo/NB_ReSeq_Files/Library501cap13/ --inline-inline -c -q -D
File	Retained Reads	Low Quality	Ambiguous Barcodes	Trimmed Reads	Orphaned paired-end reads	Total
D501_cap1_CKDL240018752-1A_H3VJJDSXC_L3_1.fq.gz	5269551	1169	1056184	0	0	6326904

Total Sequences	6326904
Ambiguous Barcodes	1056184
Low Quality	1169
Trimmed Reads	0
Orphaned Paired-ends	0
Retained Reads	5269551

Barcode	Total	Retained
ATCACG-ATCACG	1574884	1574550
ACAGTG-ACAGTG	1638722	1638334
TTAGGC-TTAGGC	2057114	2056667

Sequences not recorded
Barcode	Total
CCATGG-CCATGG	148970
TTAGGC-GGGGGG	15672
ATCACG-GGGGGG	10376
TTAGGC-TAGGCT	4002
CCATGG-CATGGA	3868
ACAGTG-GGGGGG	3188
CCATGG-GGGGGG	2776
ATCACG-AGGGGG	2534
TAGGCT-TTAGGC	2496
TTAGGC-TTAGCT	1520
AAGTGT-ACAGTG	1484
ACAGTG-AAGTGT	1188
TTAGCT-TTAGGC	1140
CATGGA-CCATGG	978
CCATGG-CATTGA	878
TTGGCT-TTAGGC	856
TTAGGC-TTGGCT	738
CCATGG-CAATGG	678
ATACGT-ATCACG	636
CCATGG-ACATGG	632
CCATTG-CCATGG	620
ATCACG-ATACGT	508
CCATGG-CGTTGA	502
ATCACG-AGGGCG	460
CCATGG-CGTGGA	458
CAGTGT-ACAGTG	414
ATCACG-TCACGT	342
CCATGG-CCTGGA	336
CCATGT-CCATGG	334
ACGTGT-ACAGTG	314
CCATGG-CCTTGA	312
```

**The top sequences not recorded belong to the NCOI adapter used for the probes.**

Because two inline adapter are use, the demultiplexed files are named as `sample-{barcode}-{barcode}.1.fq.gz`. The bash script is looking for a sample name with only one barcode, so I modified the `Rename_for_dDocent.sh` script so it can recognize the double barcode in the name.  

In `Library501cap13` Directory

```
# rename for dDocent
$ bash ../Rename_for_dDocent.sh Library501cap13_barcodes.txt
```

Create directory for files that failed demultiplexing

```
$ mkdir Lib501cap13_failed_demultiplexed
$ mv sample* Lib501cap13_failed_demultiplexed
```

Moved all demultiplexed fq.gz files from their library directories to the main directory for the resequenced data 

In `NB_ReSeq_Files`

```
$ cp Library501cap13/*.fq.gz .
```

Then made a new directory for demultiplexed qc results 

```
$ mkdir demultiplexed_fastqc_results
$ cd demultiplexed_fastqc_results
$ fastqc ../*fq.gz
$ mv *fastqc.* demultiplexed_fastqc_results
```

Create multiqc report

```
$ multiqc .

# Rename 
$ mv multiqc_report.html reseq_demultiplexed_multiqc_report.html
```

Multiqc report can be viewed [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/final_analysis_outputs/multiqc_reports/reseq_demultiplexed_multiqc_report.html).

We didn't get a ton of data back for a lot of the samples. I can still combine this new data with the old data for these samples to hopefully increase the number of reads per sample and reduce missing data later on. 

To combine, first I will rename the new files to differentiate them from the old files.

In the `NB_ReSeq_Files` directory

```
declare -a StringArray=("BAR_10" "BAR_1" "BAR_2" "BAR_3" "BAR_4" "BAR_5" "BAR_6" "BAR_7" "BAR_8" "BAR_9" "BIS_3" "BIS_6" "BIS_7" "BIS_8" "GB_8" "KIC_10" "KIC_1" "KIC_2" "KIC_3" "KIC_4" "KIC_5" "KIC_6" "KIC_7" "KIC_8" "KIC_9" "MCD_10" "MCD_1" "MCD_2" "MCD_3" "MCD_4" "MCD_5" "MCD_6" "MCD_7" "MCD_8" "MCD_9" "PVD_1" "PVD_2" "PVD_3" "PVD_4" "PVD_5")

for i in "${StringArray[@]}"
do
mv ${i}.F.fq.gz ${i}_v2.F.fq.gz
mv ${i}.R.fq.gz ${i}_v2.R.fq.gz
done
```

In the `NB_ddhaplo` directory, I'm going to make a new directory where I will compile all old and new files for the samples that were resequenced. 

First, I will link the old version of the samples. 

```
$ mkdir NB_combined_samples
$ cd NB_combined_samples

$ ln -s ~/NB_capture_Set2/BAR*.fq.gz .
# repeated for KIC and MCD

$ ln -s ../GB_8*.fq.gz
# Repeat for the specific BIS and PVD samples 
```

Then, I'll rename these older versions to include `_v1` using a for loop. 

```
declare -a StringArray=("BAR_10" "BAR_1" "BAR_2" "BAR_3" "BAR_4" "BAR_5" "BAR_6" "BAR_7" "BAR_8" "BAR_9" "BIS_3" "BIS_6" "BIS_7" "BIS_8" "GB_8" "KIC_10" "KIC_1" "KIC_2" "KIC_3" "KIC_4" "KIC_5" "KIC_6" "KIC_7" "KIC_8" "KIC_9" "MCD_10" "MCD_1" "MCD_2" "MCD_3" "MCD_4" "MCD_5" "MCD_6" "MCD_7" "MCD_8" "MCD_9" "PVD_1" "PVD_2" "PVD_3" "PVD_4" "PVD_5")

for i in "${StringArray[@]}"
do
mv ${i}.F.fq.gz ${i}_v1.F.fq.gz
mv ${i}.R.fq.gz ${i}_v1.R.fq.gz
done
```

Then, I can link the version 2 files to this directory.

```
$ ln -s ../NB_ReSeq_Files/*.fq.gz . 
```

To combine all samples, I'm going to use a for-loop to merge the two versions for each sample. Repeat for the reverse reads files. 

```
declare -a StringArray=("BAR_10" "BAR_1" "BAR_2" "BAR_3" "BAR_4" "BAR_5" "BAR_6" "BAR_7" "BAR_8" "BAR_9" "BIS_3" "BIS_6" "BIS_7" "BIS_8" "GB_8" "KIC_10" "KIC_1" "KIC_2" "KIC_3" "KIC_4" "KIC_5" "KIC_6" "KIC_7" "KIC_8" "KIC_9" "MCD_10" "MCD_1" "MCD_2" "MCD_3" "MCD_4" "MCD_5" "MCD_6" "MCD_7" "MCD_8" "MCD_9" "PVD_1" "PVD_2" "PVD_3" "PVD_4" "PVD_5")

for i in "${StringArray[@]}"
do
    file1="${i}_v1.R.fq.gz"
    file2="${i}_v2.R.fq.gz"

    output="${i}_combined.R.fq.gz"

    if [[ -f "$file1" && -f "$file2" ]]; then

    cat "$file1" "$file2" > "$output"

    fi

done
``` 

Performing raw counts to make sure v1 and v2 files combineed properly:

```
$ for fq in *.fq.gz
> do
> echo $fq
> zcat $fq | echo $((`wc -l`/4))
> done
```

I'm glad I did this because I realized the forward read file for BAR_10_v1 had been overwritten. Resolved the issue 

Making new directory in `NB_ddhaplo` to add all the combined files and remaining samples for BIS, GB, and PVD. 

```
$ mkdir NB_ddhaplo_working
$ cd NB_ddhaplo_working 

$ ln -s ../NB_combined_samples/*combined*.fq.gz . 
$ ln -s ../BIS_1*.fq.gz . 
```

## Read Trimming, Mapping, and SNP Calling

Using [**dDocent**](http://www.ddocent.com/)

> dDocent is a bioinformatics program created by Dr. Jon Purtiz that is specifically designed for different types of RAD sequencing.

**Jon downloaded a version of dDocent on my KITT account `dDocent_ngs` that can be used for Expressed Exome Capture Sequencing (EecSeq). It is accessible [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Scripts/dDocent_ngs).**

If this is your first time running dDocent, I recommend going through the [Quick Start Guide](http://www.ddocent.com/quick/). I also recommend:

1. Reading through the [User Guide](http://www.ddocent.com/UserGuide/).
2. Completing the [Assembly Tutorial](http://www.ddocent.com/assembly/), using the simulated dataset.

**Create and activate a dDocent conda environment (if you did not do so previously):**

```
$ conda create -n nb_capture ddocent
$ conda activate nb_capture

# the beginning of the line should then look like: (nb_capture) [azyck@KITT ~]$
```

In `NB_ddhaplo_working` directory:

```
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

dDocent version 3.0.0 started Wed Jun 26 15:45:13 EDT 2024 

60 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes
Proceeding with 60 individuals
dDocent detects 80 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
20

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
yes

Do you want to perform an assembly?
Type yes or no and press [ENTER].
no

Reference contigs need to be in a file named reference.fasta

Do you want to map reads?  Type yes or no and press [ENTER]
no
Mapping will not be performed
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

Trimming reads
\
dDocent has finished with an analysis in /home/azyck/NB_capture_both/NB_ddhaplo/NB_ddhaplo_working 

dDocent started Wed Jun 26 15:45:13 EDT 2024 

dDocent finished Wed Jun 26 16:02:02 EDT 2024 

 

dDocent 3.0.0 
The 'd' is silent, hillbilly.
```

Running fastqc again to check quality of cleaned reads:

```
$ mkdir cleanedreads_fastqc_results
$ cd cleanedreads_fastqc_results/

$ ln -s ../*.fq.gz .

$ fastqc *.R1.fq.gz
$ fastqc *.R2.fq.gz
```

```
$ multiqc .

$ mv multiqc_report.html 6pops_cleanedreads_multiqc_report.html
```

Multiqc report can be viewed [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/final_analysis_outputs/multiqc_reports/6pops_cleanedreads_multiqc_report.html). Report looks good, essentially a mix of past QC reports for these samples but with more data.

__Read Mapping__

Copy `reference.fasta` file into this directory so I can do alignment

```
$ cp ../reference.fasta .
```

Run dDocent just with read mapping. I'm going to use the same parameters I determined to best most optimal previously: A = 2, B = 4, and O = 6. I'll compare to previously runs to see how the output looks before moving forward. 

```
$ bash dDocent_ngs
```

```
dDocent 3.0.0 

Contact jpuritz@uri.edu with any problems 

 
Checking for required software

All required software is installed!

dDocent version 3.0.0 started Wed Jun 26 20:42:39 EDT 2024 

60 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes
Proceeding with 60 individuals
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

dDocent has finished with an analysis in /home/azyck/NB_capture_both/NB_ddhaplo/NB_ddhaplo_working 

dDocent started Wed Jun 26 20:42:39 EDT 2024 

dDocent finished Thu Jun 27 08:01:27 EDT 2024
```

Use `samtools flagstat` to investigate `bam` files.

```
$ samtools flagstat BAR_3_combined.F.bam

output:
17851434 + 0 in total (QC-passed reads + QC-failed reads)
17851434 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
17851434 + 0 mapped (100.00% : N/A)
17851434 + 0 primary mapped (100.00% : N/A)
17851434 + 0 paired in sequencing
8932972 + 0 read1
8918462 + 0 read2
17541734 + 0 properly paired (98.27% : N/A)
17837952 + 0 with itself and mate mapped
13482 + 0 singletons (0.08% : N/A)
186070 + 0 with mate mapped to a different chr
186070 + 0 with mate mapped to a different chr (mapQ>=5)
```

```
$ samtools flagstat KIC_9_combined.F.bam

output:
9684983 + 0 in total (QC-passed reads + QC-failed reads)
9684983 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
9684983 + 0 mapped (100.00% : N/A)
9684983 + 0 primary mapped (100.00% : N/A)
9684983 + 0 paired in sequencing
4847161 + 0 read1
4837822 + 0 read2
9440965 + 0 properly paired (97.48% : N/A)
9673950 + 0 with itself and mate mapped
11033 + 0 singletons (0.11% : N/A)
155993 + 0 with mate mapped to a different chr
155993 + 0 with mate mapped to a different chr (mapQ>=5)
```

```
$ samtools flagstat GB_6.F.bam

14061614 + 0 in total (QC-passed reads + QC-failed reads)
14061614 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
14061614 + 0 mapped (100.00% : N/A)
14061614 + 0 primary mapped (100.00% : N/A)
14061614 + 0 paired in sequencing
7032151 + 0 read1
7029463 + 0 read2
13758781 + 0 properly paired (97.85% : N/A)
14048216 + 0 with itself and mate mapped
13398 + 0 singletons (0.10% : N/A)
195111 + 0 with mate mapped to a different chr
195111 + 0 with mate mapped to a different chr (mapQ>=5)
```

```
$ samtools flagstat PVD_4_combined.F.bam

3926336 + 0 in total (QC-passed reads + QC-failed reads)
3926336 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
3926336 + 0 mapped (100.00% : N/A)
3926336 + 0 primary mapped (100.00% : N/A)
3926336 + 0 paired in sequencing
1963638 + 0 read1
1962698 + 0 read2
3853324 + 0 properly paired (98.14% : N/A)
3923901 + 0 with itself and mate mapped
2435 + 0 singletons (0.06% : N/A)
49513 + 0 with mate mapped to a different chr
49513 + 0 with mate mapped to a different chr (mapQ>=5)
```

Since I previously optimized the alignment metrics and the output looks very similar to the previous set, I'm going to move forward. 

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

dDocent version 3.0.0 started Thu Jun 27 09:21:43 EDT 2024 

60 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes
Proceeding with 60 individuals
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

Using VCFtools to parse TotalRawSNPS.vcf for SNPs that are called in at least 90% of individuals

dDocent has finished with an analysis in /home/azyck/NB_capture_both/NB_ddhaplo/NB_ddhaplo_working 

dDocent started Thu Jun 27 09:21:43 EDT 2024 

dDocent finished Thu Jun 27 21:18:59 EDT 2024 

After filtering, kept 10921812 out of a possible 52027016 Sites
```

## Variant Filtering

I'm going to make a new directory within `NB_ddhaplo_working`

```
$ mkdir NB_SNPFiltering_working
$ cd NB_SNPFiltering_working

# Link Raw VCF file to new directory
$ ln -s ../TotalRawSNPs.vcf.gz .
```

**First step for filtering is to remove sites with a QUALITY score below 20 and those sites that have more than 50% missing data.** I talked to Jon about this and it makes the most sense to just remove 50% of the low call genotypes now. That would happen eventually and this will make the filtering much faster. 

```
$ vcftools --gzvcf TotalRawSNPs.vcf.gz --minQ 20 --max-missing 0.50 --recode --recode-INFO-all --out TRS6new

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 33237806 out of a possible 52027016 Sites
Run Time = 7797.00 seconds
```

**Next step is to mark any genotypes with less than 10 reads as missing.**

```
$ vcftools --vcf TRS6new.recode.vcf --minDP 10 --recode --recode-INFO-all --out TRS6newdp10

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 33237806 out of a possible 33237806 Sites
Run Time = 7479.00 seconds
```

**MAF filtering for FreeBayes Output.**

```
$ curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/untested/multi.maf.sh
$ chmod +x multi.maf.sh
$ multi.maf.sh TRS6newdp10.recode.vcf 0.001 TRS6newdp10maf

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 32756214 out of a possible 33237806 Sites
Run Time = 6565.00 seconds
```

**Use a custom script called `filter_missing_ind.sh` to filter out bad individuals.**

Here I mostly just want to see % of missing data within individuals. We are working with a small number of individuals, so I don't want to remove any.

_I gzipped all vcf files in this directory before this step to free up space in KITT. I modified the filter_missing_ind script to recognize the gzipped vcf file input._

```
$ curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/filter_missing_ind.sh
$ chmod +x filter_missing_ind.sh
$ ./filter_missing_ind.sh TRS6newdp10maf.recode.vcf.gz TRS6newdp10mafMI

output:


                                          Histogram of % missing data per individual
       12 ++-------------+--------------+--------------+---------------+--------------+--------------+-------------++
          +              +              +              +       'totalmissing' using (bin($1,binwidth)):(1.0) ****** +
          |                                                                           ********                      |
          |                                                                           *      *                      |
       10 ++                                                                          *      *                     ++
          |                                                                           *      *                      |
          |                                                    *********              *      *                      |
          |                                                    *       *              *      *                      |
        8 ++                                                   *       *              *      *                     ++
          |                                                    *       *              *      *                      |
          |                             *********              *       *              *      *                      |
        6 ++                            *       *              *       *      *********      *                     ++
          |                             *       *              *       *      *       *      *                      |
          |                             *       *      *********       *      *       *      *********              |
          |                             *       *      *       *       *      *       *      *       *              |
        4 ++                            *       ********       *       ********       *      *       *      *********
          |                             *       *      *       *       *      *       *      *       *      *       *
          |                             *       *      *       *       *      *       *      *       *      *       *
          |                             *       *      *       *       *      *       *      *       *      *       *
        2 ++      ********              *       *      *       *       *      *       *      *       *      *      +*
          |       *      *              *       *      *       *       *      *       *      *       *      *       *
          |       *      ****************       *      *       *       *      *       *      *       ********       *
          +       *      *       *      *       *      *       *       *      *       *      *       *      *       *
        0 ++------***************************************************************************************************
         0.86           0.88           0.9            0.92            0.94           0.96           0.98            1
                                                       % of missing data

The 85% cutoff would be 0.972706
Would you like to set a different cutoff, yes or no
1

VCFtools - 0.1.17
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
        --gzvcf TRS6newdp10maf.recode.vcf.gz
        --remove lowDP.indv
        --recode-INFO-all
        --out TRS6newdp10mafMI
        --recode


Excluding individuals in 'exclude' list
After filtering, kept 52 out of 60 Individuals
Outputting VCF file...
After filtering, kept 32756214 out of a possible 32756214 Sites                             
Run Time = 6193.00 seconds
```

Here is a list of the worst 12 samples:

```
$ cat <(head -1 TRS6newdp10mafMI.imiss ) <(mawk '!/F_MI/' TRS6newdp10mafMI.imiss | sort -k5 -r -n ) | head -13 | column -t

output: 

INDV             N_DATA    N_GENOTYPES_FILTERED  N_MISS    F_MISS
PVD_4_combined   32756214  0                     32508284  0.992431
PVD_2_combined   32756214  0                     32460113  0.99096
PVD_1_combined   32756214  0                     32456357  0.990846
PVD_5_combined   32756214  0                     32446781  0.990553
PVD_3_combined   32756214  0                     32397576  0.989051
KIC_10_combined  32756214  0                     32007418  0.97714
MCD_1_combined   32756214  0                     31989326  0.976588
KIC_5_combined   32756214  0                     31976544  0.976198
MCD_10_combined  32756214  0                     31862169  0.972706
MCD_9_combined   32756214  0                     31815718  0.971288
MCD_3_combined   32756214  0                     31769361  0.969873
MCD_6_combined   32756214  0                     31758913  0.969554
```

Because this step filtered out individuals, I won't use it moving forward, instead I'll keep going with `TRS6newdp10maf.recode.vcf.gz`.

Moving forward with three different filtering pathways: 

- One with a popmap filter and matching max-missing filter of 0.1 (more stringent), 
- 0.5 (medium stringency), 
- and 0.8 (less stringent).

I want to test different levels of filtering to see if there are significant differences in the results. 

**Use a second custom `script pop_missing_filter.sh` to filter loci that have high missing data values in a single population.**

This step needs a file that maps individuals to populations `popmap`. dDocent automatically generates this for you based on the file names, so make sure to use a naming convention where population name is first followed by a unique identifier. 

```
$ head popmap

output:
BAR_10_combined	BAR
BAR_1_combined	BAR
BAR_2_combined	BAR
BAR_3_combined	BAR
BAR_4_combined	BAR
BAR_5_combined	BAR
BAR_6_combined	BAR
BAR_7_combined	BAR
BAR_8_combined	BAR
BAR_9_combined	BAR
```

Moving forward with three different filtering pathways: 

One with a popmap filter and matching max-missing filter of 0.1 (more stringent), 0.5 (medium stringency), and 0.8 (less stringent). I want to test different levels of filtering to see if there are significant differences in the results. 

```
$ curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/pop_missing_filter.sh
$ chmod +x pop_missing_filter.sh

# more stringent pop filter of 0.1 
$ ./pop_missing_filter.sh TRS6newdp10maf.recode.vcf.gz popmap 0.1 1 TRS6newdp10mafp9

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 131117 out of a possible 32756214 Sites
Run Time = 593.00 seconds

# medium stringent level of 0.5 
$ ./pop_missing_filter.sh TRS6newdp10maf.recode.vcf.gz popmap 0.5 1 TRS6newdp10mafp5

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 498511 out of a possible 32756214 Sites
Run Time = 685.00 seconds

# less stringent level of 0.8 
$ ./pop_missing_filter.sh TRS6newdp10maf.recode.vcf.gz popmap 0.8 1 TRS6newdp10mafp2

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 1233427 out of a possible 32756214 Sites
Run Time = 826.00 seconds
```

**Filter out any sites with less than % overall call rate and MAF of 0.001.**

I'm going to repeat the same filtering level for the overall call rate as the popmap filter step - specific to each vcf pathway. 

```
# more stringent filtering - max-missing of 0.9 
$ vcftools --gzvcf TRS6newdp10mafp9.recode.vcf.gz --recode-INFO-all --max-missing 0.9 --maf 0.001 --out TRS6newdp10mafp9g9 --recode

output: 
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 130702 out of a possible 131117 Sites
Run Time = 29.00 seconds

# medium stringent level of 0.5 
$ vcftools --gzvcf TRS6newdp10mafp5.recode.vcf.gz --recode-INFO-all --max-missing 0.5 --maf 0.001 --out TRS6newdp10mafp5g5 --recode

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 475519 out of a possible 498511 Sites
Run Time = 110.00 seconds

# less stringent filtering - max-missing of 0.2
$ vcftools --gzvcf TRS6newdp10mafp2.recode.vcf.gz --recode-INFO-all --max-missing 0.2 --maf 0.001 --out TRS6newdp10mafp2g2 --recode

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 1042212 out of a possible 1233427 Sites
Run Time = 244.00 seconds
```

**Next, split the VCF files into nDNA and mtDNA for further filtering**

```
# More stringently filtered vcf
$ zcat TRS6newdp10mafp9g9.recode.vcf.gz | mawk '!/NC_007175.2/' > TRS6newdp10mafp9g9.nDNA.vcf

$ zcat TRS6newdp10mafp9g9.recode.vcf.gz | head -10000 | mawk '/#/' > header
$ cat header <(zcat TRS6newdp10mafp9g9.recode.vcf.gz |head -20000| mawk '/NC_007175.2/') > TRS6newdp10mafp9g9.mtDNA.vcf

# Medium stringently filtered vcf
$ zcat TRS6newdp10mafp5g5.recode.vcf.gz | mawk '!/NC_007175.2/' > TRS6newdp10mafp5g5.nDNA.vcf

$ zcat TRS6newdp10mafp5g5.recode.vcf.gz | head -10000 | mawk '/#/' > header
$ cat header <(zcat TRS6newdp10mafp5g5.recode.vcf.gz |head -20000| mawk '/NC_007175.2/') > TRS6newdp10mafp5g5.mtDNA.vcf

# Less stringently filtered vcf
$ zcat TRS6newdp10mafp2g2.recode.vcf.gz | mawk '!/NC_007175.2/' > TRS6newdp10mafp2g2.nDNA.vcf

$ zcat TRS6newdp10mafp2g2.recode.vcf.gz | head -10000 | mawk '/#/' > header
$ cat header <(zcat TRS6newdp10mafp2g2.recode.vcf.gz |head -20000| mawk '/NC_007175.2/') > TRS6newdp10mafp2g2.mtDNA.vcf
```

The remaining filtering steps will be split between the mtDNA and nDNA

#### mtDNA

**Break complex mutational events (combinations of SNPs and INDELs) into separate SNP and INDEL calls, and then remove INDELs.**

```
# More stringently filtered vcf
$ vcffilter -s -f "AB < 0.001" TRS6newdp10mafp9g9.mtDNA.vcf | vcffilter -s -f "QUAL / DP > 0.25" > TRS6newdp10mafp9g9.mtDNA.F.vcf
$ vcfallelicprimitives -k TRS6newdp10mafp9g9.mtDNA.F.vcf | sed 's:\.|\.:\.\/\.:g' > TRS6newdp10mafp9g9.F.prim
$ vcftools --vcf TRS6newdp10mafp9g9.F.prim --remove-indels --recode --recode-INFO-all --out SNP.TRS6newdp10mafp9g9mtDNA

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 593 out of a possible 609 Sites
Run Time = 0.00 seconds

# Medium stringently filtered vcf
$ vcffilter -s -f "AB < 0.001" TRS6newdp10mafp5g5.mtDNA.vcf | vcffilter -s -f "QUAL / DP > 0.25" > TRS6newdp10mafp5g5.mtDNA.F.vcf
$ vcfallelicprimitives -k TRS6newdp10mafp5g5.mtDNA.F.vcf | sed 's:\.|\.:\.\/\.:g' > TRS6newdp10mafp5g5.F.prim
$ vcftools --vcf TRS6newdp10mafp5g5.F.prim --remove-indels --recode --recode-INFO-all --out SNP.TRS6newdp10mafp5g5mtDNA

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 593 out of a possible 609 Sites
Run Time = 0.00 seconds

# Less stringently filtered vcf
$ vcffilter -s -f "AB < 0.001" TRS6newdp10mafp2g2.mtDNA.vcf | vcffilter -s -f "QUAL / DP > 0.25" > TRS6newdp10mafp2g2.mtDNA.F.vcf
$ vcfallelicprimitives -k TRS6newdp10mafp2g2.mtDNA.F.vcf | sed 's:\.|\.:\.\/\.:g' > TRS6newdp10mafp2g2.F.prim
$ vcftools --vcf TRS6newdp10mafp2g2.F.prim --remove-indels --recode --recode-INFO-all --out SNP.TRS6newdp10mafp2g2mtDNA

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 593 out of a possible 609 Sites
Run Time = 0.00 seconds
```

#### nDNA

**Use a custom filter script `dDocent_ngs_filters`.**

```
$ curl -L -O https://raw.githubusercontent.com/amyzyck/EecSeq_NB_EasternOyster/master/Scripts/dDocent_ngs_filters
$ chmod +x dDocent_ngs_filters

# More stringently filtered vcf 
$ ./dDocent_ngs_filters TRS6newdp10mafp9g9.nDNA.vcf TRS6newdp10mafp9g9nDNA

output:
This script will automatically filter a FreeBayes generated VCF file using criteria related to site depth,
quality versus depth, allelic balance at heterzygous individuals, and paired read representation.
The script assumes that loci and individuals with low call rates (or depth) have already been removed. 

Contact Jon Puritz (jpuritz@gmail.com) for questions and see script comments for more details on particular filters 

Is this from a mixture of SE and PE libraries? Enter yes or no.
no
Number of sites filtered based on allele balance at heterozygous loci, locus quality, and mapping quality / Depth
 22944 of 129989 

Number of additional sites filtered based on properly paired status
 33 of 107045 

Number of sites filtered based on high depth and lower than 2*DEPTH quality score
 17663 of 107045 



                                             Histogram of mean depth per site

  4000 ++----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----++
       +     +     +     +     +     +     +     +     +   'meandepthpersite' using (bin($1,binwidth)):(1.0) ******+|
       |                                                                                                            |
  3500 ++                      **                                                                                  ++
       |                     *********                                                                              |
       |                   ******* ***                                                                              |
  3000 ++                 ** ***** *****                                                                           ++
       |                  ** ***** *******                                                                          |
       |                 *** ***** ***** ***                                                                        |
  2500 ++               **** ***** ***** ***                                                                       ++
       |                **** ***** ***** *****                                                                      |
  2000 ++               **** ***** ***** *******                                                                   ++
       |               ***** ***** ***** ***** ****                                                                 |
       |               ***** ***** ***** ***** *****                                                                |
  1500 ++            ******* ***** ***** ***** *********                                                           ++
       |             * ***** ***** ***** ***** ***** ****                                                           |
       |             * ***** ***** ***** ***** ***** ********                                                       |
  1000 ++           ** ***** ***** ***** ***** ***** ***** ****                                                    ++
       |            ** ***** ***** ***** ***** ***** ***** ****                                                     |
       |           *** ***** ***** ***** ***** ***** ***** ****                                                     |
   500 ++          *** ***** ***** ***** ***** ***** ***** ****                                                    ++
       |          **** ***** ***** ***** ***** ***** ***** ****** ***                                               |
       +     +   ***** ***** ***** ***** ***** ***** ***** ***************************************************     +|
     0 ++----********************************************************************************************************
       10    15    20    25    30    35    40    45    50    55    60    65    70    75    80    85    90    95   100
                                                        Mean Depth

The 95% cutoff would be 81
Would you like to use a different maximum mean depth cutoff than 81, yes or no
no
Maximum mean depth cutoff is 81
Number of sites filtered based on maximum mean depth
 4489 of 89356 

Total number of sites filtered
 45122 of 129989 

Remaining sites
 84867 

Filtered VCF file is called Output_prefix.FIL.recode.vcf

Filter stats stored in TRS6newdp10mafp9g9nDNA.filterstats
```

```
# Medium stringently filtered vcf 
$ ./dDocent_ngs_filters TRS6newdp10mafp5g5.nDNA.vcf TRS6newdp10mafp5g5nDNA

output:
Is this from a mixture of SE and PE libraries? Enter yes or no.
no
Number of sites filtered based on allele balance at heterozygous loci, locus quality, and mapping quality / Depth
 96234 of 474806 

Number of additional sites filtered based on properly paired status
 191 of 378572 

Number of sites filtered based on high depth and lower than 2*DEPTH quality score
 55155 of 378572 



                                              Histogram of mean depth per site

  30000 ++--+------+-----+-----+------+-----+-----+------+-----+-----+------+-----+-----+------+-----+------+-----+++
        |   +      +     +     +      +     +     +      + 'meandepthpersite' using (bin($1,binwidth)):(1.0)+****** |
        |                                                                                                           |
        |       ****                                                                                                |
  25000 ++      *  ***                                                                                             ++
        |     ***  * ***                                                                                            |
        |     * *  * * *                                                                                            |
        |     * *  * * ***                                                                                          |
  20000 ++    * *  * * * ***                                                                                       ++
        |     * *  * * * * *                                                                                        |
        |   *** *  * * * * ***                                                                                      |
  15000 ++  * * *  * * * * * *                                                                                     ++
        |   * * *  * * * * * ***                                                                                    |
        |   * * *  * * * * * * ***                                                                                  |
        |   * * *  * * * * * * * ****                                                                               |
  10000 ++  * * *  * * * * * * * *  *****                                                                          ++
        | *** * *  * * * * * * * *  * * ***                                                                         |
        | * * * *  * * * * * * * *  * * * *****                                                                     |
        | * * * *  * * * * * * * *  * * * * * *****                                                                 |
   5000 ++* * * *  * * * * * * * *  * * * * * * * ***                                                              ++
        | * * * *  * * * * * * * *  * * * * * * * * *                                                               |
        *** * * *  * * * * * * * *  * * * * * * * * ****                                                            |
        * * * * *  * * * * * * * *  * * * * * * * * *  ******************************** +      +     +      +     + |
      0 *************************************************************************************************************
            12     15    18    21     24    27    30     33    36    39     42    45    48     51    54     57    60
                                                         Mean Depth

The 95% cutoff would be 49
Would you like to use a different maximum mean depth cutoff than 49, yes or no
no
Maximum mean depth cutoff is 49
Number of sites filtered based on maximum mean depth
 16347 of 323238 

Total number of sites filtered
 167915 of 474806 

Remaining sites
 306891 

Filtered VCF file is called Output_prefix.FIL.recode.vcf

Filter stats stored in TRS6newdp10mafp5g5nDNA.filterstats
```

```
# Less stringently filtered vcf 
$ ./dDocent_ngs_filters TRS6newdp10mafp2g2.nDNA.vcf TRS6newdp10mafp2g2nDNA

output:
Is this from a mixture of SE and PE libraries? Enter yes or no.
no
Number of sites filtered based on allele balance at heterozygous loci, locus quality, and mapping quality / Depth
 208213 of 1041499 

Number of additional sites filtered based on properly paired status
 724 of 833286 

Number of sites filtered based on high depth and lower than 2*DEPTH quality score
 108451 of 833286 



                                              Histogram of mean depth per site

  80000 ****---+------+------+------+------+------+------+------+------+------+------+------+------+------+------+-++
        *  *   +      +      +      +      +      +      + 'meandepthpersite' using (bin($1,binwidth)):(1.0) ****** |
        *  *                                                                                                        |
  70000 *+ *                                                                                                       ++
        *  *****                                                                                                    |
        *  *   *                                                                                                    |
  60000 *+ *   *                                                                                                   ++
        *  *   ****                                                                                                 |
        *  *   *  *                                                                                                 |
  50000 *+ *   *  *                                                                                                ++
        *  *   *  *****                                                                                             |
  40000 *+ *   *  *   *                                                                                            ++
        *  *   *  *   ****                                                                                          |
        *  *   *  *   *  *                                                                                          |
  30000 *+ *   *  *   *  *****                                                                                     ++
        *  *   *  *   *  *   ****                                                                                   |
        *  *   *  *   *  *   *  *****                                                                               |
  20000 *+ *   *  *   *  *   *  *   ****                                                                           ++
        *  *   *  *   *  *   *  *   *  *****                                                                        |
        *  *   *  *   *  *   *  *   *  *   ****                                                                     |
  10000 *+ *   *  *   *  *   *  *   *  *   *  *                                                                    ++
        *  *   *  *   *  *   *  *   *  *   *  *****                                                                 |
        *  *   *  *   *  *   *  *   *  *   *  *   *******************************************      +      +      +  |
      0 *************************************************************************************************************
        10     12     14     16     18     20     22     24     26     28     30     32     34     36     38     40
                                                         Mean Depth

The 95% cutoff would be 33
Would you like to use a different maximum mean depth cutoff than 33, yes or no
no
Maximum mean depth cutoff is 33
Number of sites filtered based on maximum mean depth
 37656 of 724146 

Total number of sites filtered
 355009 of 1041499 

Remaining sites
 686490 

Filtered VCF file is called Output_prefix.FIL.recode.vcf

Filter stats stored in TRS6newdp10mafp2g2nDNA.filterstats
```

**Break complex mutational events (combinations of SNPs and INDELs) into separate SNP and INDEL calls, and then remove INDELs.**

```
# More stringently filtered vcf 
$ vcfallelicprimitives TRS6newdp10mafp9g9nDNA.FIL.recode.vcf --keep-info --keep-geno > TRS6newdp10mafp9g9nDNA.prim.vcf
$ vcftools --vcf TRS6newdp10mafp9g9nDNA.prim.vcf --remove-indels --recode --recode-INFO-all --out SNP.TRS6newdp10mafp9g9nDNA

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 63948 out of a possible 76865 Sites
Run Time = 13.00 seconds

# Medium stringently filtered vcf 
$ vcfallelicprimitives TRS6newdp10mafp5g5nDNA.FIL.recode.vcf --keep-info --keep-geno > TRS6newdp10mafp5g5nDNA.prim.vcf
$ vcftools --vcf TRS6newdp10mafp5g5nDNA.prim.vcf --remove-indels --recode --recode-INFO-all --out SNP.TRS6newdp10mafp5g5nDNA

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 216398 out of a possible 266479 Sites
Run Time = 44.00 seconds

# Less stringently filtered vcf 
$ vcfallelicprimitives TRS6newdp10mafp2g2nDNA.FIL.recode.vcf --keep-info --keep-geno > TRS6newdp10mafp2g2nDNA.prim.vcf
$ vcftools --vcf TRS6newdp10mafp2g2nDNA.prim.vcf --remove-indels --recode --recode-INFO-all --out SNP.TRS6newdp10mafp2g2nDNA

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 466067 out of a possible 569535 Sites
Run Time = 92.00 seconds
```

***
**Filtered by Minor Allele Frequencies**

Minor Allele Frequency greater than 5%

```
# More stringently filtered vcf 
$ vcftools --vcf SNP.TRS6newdp10mafp9g9nDNA.recode.vcf --maf 0.05 --recode --recode-INFO-all --out SNP.TRS6newdp10mafp9g9nDNAmaf05

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 23231 out of a possible 63948 Sites
Run Time = 6.00 seconds

# Medium stringently filtered vcf
$ vcftools --vcf SNP.TRS6newdp10mafp5g5nDNA.recode.vcf --maf 0.05 --recode --recode-INFO-all --out SNP.TRS6newdp10mafp5g5nDNAmaf05

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 83242 out of a possible 216398 Sites
Run Time = 19.00 seconds

# Less stringently filtered vcf
$ vcftools --vcf SNP.TRS6newdp10mafp2g2nDNA.recode.vcf --maf 0.05 --recode --recode-INFO-all --out SNP.TRS6newdp10mafp2g2nDNAmaf05

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 215841 out of a possible 466067 Sites
Run Time = 50.00 seconds
```

**Max 2 alleles**

```
# More stringently filtered vcf
$ vcftools --vcf SNP.TRS6newdp10mafp9g9nDNAmaf05.recode.vcf --recode --recode-INFO-all --out SNP.TRS6newdp10mafp9g9nDNAmaf052A --max-alleles 2

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 22874 out of a possible 23231 Sites
Run Time = 5.00 seconds

# Medium stringently filtered vcf
$ vcftools --vcf SNP.TRS6newdp10mafp5g5nDNAmaf05.recode.vcf --recode --recode-INFO-all --out SNP.TRS6newdp10mafp5g5nDNAmaf052A --max-alleles 2

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 81872 out of a possible 83242 Sites
Run Time = 17.00 seconds

# Less stringently filtered vcf
$ vcftools --vcf SNP.TRS6newdp10mafp2g2nDNAmaf05.recode.vcf --recode --recode-INFO-all --out SNP.TRS6newdp10mafp2g2nDNAmaf052A --max-alleles 2

output: 
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 211290 out of a possible 215841 Sites
Run Time = 41.00 seconds
```

There are 2-3 known large inversions in these samples. I'm going to make a vcf file with SNPs inside these inverted regions removed. I will perform all downstream analyses on this dataset and compare to the full SNP dataset. 

I have a bed file with the known inverted regions: 

```
$ cat Detected_large_inversions.bed 

output:
NC_035780.1	38560000	44960000
NC_035784.1	61560000	80150000
NC_035785.1	29630000	44000000
```

SNPs that fall within these ranges on these three chromosomes will be removed. 

```
$ vcftools --vcf SNP.TRS6newdp10mafp9g9nDNAmaf052A.recode.vcf --exclude-bed Detected_large_inversions.bed --recode-INFO-all --out SNP.TRS6newdp10mafp9g9nDNAmaf052A_noinvert --recode

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
	Read 3 BED file entries.
After filtering, kept 21070 out of a possible 22874 Sites
Run Time = 5.00 seconds
```

I moved onto Outlier Detection using vcf file `SNP.TRS6newdp10mafp9g9nDNAmaf05.recode.vcf`. Steps for outlier detection are documented [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/Analysis_Final/NB_OutlierDetection_Final.Rmd).