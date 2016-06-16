# Cross Contamination

**CrossCont** is a tool to check cross contamination between different samples from the same HIV, HBV or HCV virus sequencing run.  
  
* Input: Paired-end fastq files ( R1,R2 ), could be compressed by gzip or bz2.    

* Output:   
  crossCont/**"PROJECT"_CrosscPlot.pdf** , using chord diagram in R package, we generate 2D circular data track plot where we can see if there is contamination.  

  crossCont/**CcEnd.csv** is a percentage table where columns are the contaminants and rows are the contaminated samples.  

  **decontaminated fastq files** ( crossCont/ _clean.fq files ), with sequences that had mapped to Virus Reference.   

  **consensus files** (crossCont/ _cons.fastq ).  

## Description  

** crossC ** directory has the main script  ** pCrossCont.sh ** and two folders:  
   * ** General ** with the References files.  
   * ** scripts ** with secondary scripts.   

## Prerequisities

Create **RawData** directory.  
Copy fastq files to "/**RawData**/" directory.   
    
**Software that should be installed and be in $PATH:**   
* ** Samtools ** Tools for alignments in the SAM format.  
     Version: samtools 1.2 Using htslib 1.2.1   
     http://www.htslib.org/download/   

* ** bcftools ** Is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF.  
     Version: 1.2 (using htslib 1.2.1).  
     http://www.htslib.org/doc/#manual-pages  

* ** vcfutils.pl ** It's part of the samtools distribution. If you download samtools, it should contain a "bcftools" folder, which contains vcfutils.pl. 

* ** BWA ** Burrows-Wheeler Aligner is a software package for mapping low-divergent sequences against a large reference genome.  
     Version: 0.7.13-r1126  
     https://sourceforge.net/projects/bio-bwa/files/  

* ** bbsplit.sh **  Is a tool that bins reads by mapping to multiple references simultaneously, using BBMap (aligner for DNA/RNAseq).   
     Version: BBSplit / BBMap v35.x.  
     http://seqanswers.com/forums/showthread.php?t=41288

* ** Rscript **  Scripting front-end. ( libraries: reshape2, RCircos; circlize )  
     Version 3.2.0 (2015-04-16)  
     https://cran.r-project.org/

* ** seqtk ** Toolkit for processing sequences in FASTA/Q formats.  
     Version: 1.1-r91.  
     https://github.com/lh3/seqtk  

* ** fastq-stats ** Produces statistics for the files listed.  
     Version: 1.01.  
     http://manpages.ubuntu.com/manpages/wily/man1/fastq-stats.1.html  

## Installation  
Download "crossC" folder.  
Inside will be pCrossCont.sh script.  

## Running the test

** pCrossCont.sh**   arg1: "path where RawData directory is located"      arg2: "path where crossC is downloaded"      arg3: "Sequenced virus HIV, HBV or HCV"      arg4: Minimum sequence length to filter (between 100 and 500)    arg5: "Project's name (please, only one word)."


`Ex:`   
`sh ./pCrossCont.sh "/home/user/" "/Work/downloaded/" "HIV" 100 "TestHiv"`  

## License
ViCroSeq is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ViCroSeq is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ViCroSeq.  If not, see <http://www.gnu.org/licenses/>.  

## Authors

* **Marc Noguera-Julian**  - mnoguera _at_ irsicaixa.es - 
* **M. Cristina Rodr&iacute;guez** - crodriguez _at_ irsicaixa.es - 


