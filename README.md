## License
**ViCroSeq** is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
( at your option ) any later version.

**ViCroSeq** is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ViCroSeq.  Otherwise, see <http://www.gnu.org/licenses/>.  


## Authors

* **M. Cristina Rodr&iacute;guez** - crodriguez _at_ irsicaixa.es - 
* **Marc Noguera-Julian**  - mnoguera _at_ irsicaixa.es - 

# ViCroSeq

**ViCroSeq** is a tool to check cross contamination between different samples from the same HIV, HBV or HCV virus sequencing run.  
  
* Input: Paired-end or single-end fastq files, could be compressed by gzip or bz2.    

* Output:   
  crossCont/**"PROJECT"_CrosscPlot.pdf** , using chord diagram in R package, a 2D circular data track plot is generated showing the contaminant-contaminated samples attached with arrows.  

  crossCont/**CcEnd.csv** is a percentage table where columns are the contaminants and rows are the contaminated samples.  

  crossCont/**lowcow/** files with low coverage (lower than 500) will not be processed and will be moved to this directory.  

  **consensus files** ( crossCont/ clean_consdp.fasta ).  

  **cleaned fastq files** ( crossCont/ _clean.fq files ), only sequences that had mapped to Virus Reference.   

  **decontaminated fastq files** ( crossCont/finalDir/"sample"/"sample"clean_consdp.fq ), files without crosscontamination sequences.   

## Description  

**ViCroSeq-master** directory has the main script  **pCrossCont.sh** and two folders:  
   * **General** with the References files.  
   * **scripts** with secondary scripts.   

## Prerequisities

Create **RawData** directory.  
Copy fastq files to "/**RawData**/" directory.   
    
**Software that should be installed and be in $PATH:**   
* **Samtools** Tools for alignments in the SAM format.  
     Version: samtools 1.3.1 Using htslib 1.3.1   
     https://sourceforge.net/projects/samtools/files/samtools/1.2/   

* **bcftools** Is a set of utilities that manipulate variant calls in the Variant Call Format ( VCF) and its binary counterpart BCF. It contains vcfutils.pl.  
     Version: 1.2 ( using htslib 1.2.1 )  
     https://sourceforge.net/projects/samtools/files/samtools/1.2/  

* **BWA** Burrows-Wheeler Aligner is a software package for mapping low-divergent sequences against a large reference genome.  
     Version: 0.7.10_x86_64  
     https://sourceforge.net/projects/bio-bwa/files/  

* **Bowtie**  Is an ultrafast, memory-efficient short read aligner.  
     Version: 1.1.2  
     http://bowtie-bio.sourceforge.net/index.shtml  

* **bbsplit.sh**  Is a tool that bins reads by mapping to multiple references simultaneously, using **BBMap** ( aligner for DNA/RNAseq ).   
     Version: BBMap_v36.14  
     https://sourceforge.net/projects/bbmap/  

* **Rscript**  Scripting front-end. ( libraries: reshape2, RCircos; circlize )  
     Version 3.3.1   
     https://cran.r-project.org/

* **fastq-stats** from **ea-utils** package. Produces statistics for the files listed.  
     Version: ea-utils.1.1.2-537  
     https://code.google.com/archive/p/ea-utils/downloads  

## Installation  
Download "ViCroSeq-master" folder.  
Inside will be pCrossCont.sh script.  

## Running the test

**pCrossCont.sh**   arg1: "path where RawData directory is located"      arg2: "Sequenced virus: HIV, HBV or HCV"      arg3: Minimum sequence length to filter ( between 100 and 500 )    arg4: "Paired (1) or single-end files (2)"       arg5: "Project's name ( please, only one word )."


`Ex:`   
`sh ./pCrossCont.sh "/home/user/" "HIV" "100" "1" "TestHiv"`  





