# Amplicon Sequence Analysis
This code is designed to implement the dada2 pipeline in R to pick ASVs and assign taxonomy, this pipeline has been adapted to 16S and ITS sequencing.

Many of these scripts can take a while to run depending on the size of your data and power of your computer. If you have a lot of samples, it will save you time to utilize a high performace computer. 
The RunR.sh script is an example of how to batch R scripts on an IBM hpc system. (But the R scripts can be run locally or on any platform that supports R)

## Step1: Setup
Before you get coding see if you can find out if your fastq file have adapter or sequencing primers present or already removed, either way no worries, we can easily trim the adaptors before ASV picking, for accurate trimming find out the exact length of the sequencing primers used.

We'll be running this using R on the hpc, I recommend using a custom conda environment to keep all package versions stable, can install packages via conda as shown:
```
module load conda
conda create -p /path/to/environments/myR r-essentials r-base
conda activate myR
conda install gridExtra
```
some packages (like dada2) can't be installed via conda, will need to begin an R session to install it:
```
R
install.packages("Bioconductor")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2")
```
We need to generate a meta data file that includes the relative file location of forward and reverse reads:
```
files<-list.files("./RawData")
R1<-files[grep("R1",files)]
R2<-files[grep("R2",files)]
AmpliconTest<-as.data.frame(cbind(R1,R2))
#I typically add any other metadata I have to this file and export it as a csv to use in down stream analysis
#Write a file for any additional customization in excel an usage in downstream analysis:  
write.table(AmpliconTest,"AmpliconTest.csv",sep=",",row.names = FALSE,)
```
# Visualize read quality and choose trimming values
Edit the RunR.sh script so that it only runs MakeMeta_Quality.R 
This will generate a quality plot like this:
![AveQuality.pdf](https://github.com/NateKorth/AmpliconSeqAnalysis/files/14921532/AveQuality.pdf)




