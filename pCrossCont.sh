#!/bin/bash
#===============================================================================
#
#       FILE:  pCrossCont.sh
#
#  	DESCRIPTION:  Check cross contamination between samples.
#				1-. Filter: 
#					fastq Nº of sequences R1 is different from R2.
#					coverage is < 500 (aligning with HXB2R)
#					sequences are not HIV (bbsplit).
#				2-. Generate consensus.
#				3-. Filter low coverage positions.
#				4-. Classification with bbsplit.
#			        5-. Generate table with % of contaminant's Nº of reads per sample. 
#				6-. Plot.
#
#      AUTHOR:   M. Cristina Rodríguez(crodriguez@irsicaixa.es), 
#      COMPANY:  IrsiCaixa.
#      VERSION:  1.0
#      CREATED:  2016/05/10
#      
#===============================================================================


RUNDIR=$1  	# "/Work/data/" RawData
Virus=$2        # "HIV","HBV", "HCV"
SEQlength=$3 	# "Minimum Sequence length"
PROJECT=$4 	# "Project name"

CROSSCDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"   # path where ViCroSeq-master is downloaded

if [[ "$RUNDIR" == "" || "$Virus" == "" || "$PROJECT" == "" ]] 
   then
	echo " *---*"
 	echo " pCrossCont.sh      arg1:path where RawData directory is located      arg2: HIV or HBV or HCV            arg3: Minimum sequence length (between 100 and 500)       arg4: Project's name (one word)."
	echo ".fastq.gz (R1, R2) files must be in RawData dir".
	echo " *---*"
	exit
   fi

PIPELINEDIR=$CROSSCDIR"/scripts"
REFERENCEFILE=$CROSSCDIR"/General/HXB2R.fna"   

if [[ "$Virus" != "HIV" && "$Virus" != "HBV" && "$Virus" != "HCV" ]] 
   then
	echo " *---*"
 	echo "  arg3 must be: HIV or HBV or HCV."
	echo " *---*"
	exit
   fi

if [[ "$Virus" == "HIV" ]] 
   then
	echo "HIV "
	REFERENCEFILEALL=$CROSSCDIR"/General/Alltypes_Refe.fasta" 
	REFERENCEFILE=$CROSSCDIR"/General/HXB2R.fna"
	REFCOPIED="HXB2R.fna"
   elif [[ "$Virus" == "HBV" ]] 
   then
	echo "HBV "
	REFERENCEFILEALL=$CROSSCDIR"/General/HBV_Refe.fasta" 
	REFERENCEFILE=$CROSSCDIR"/General/HBV_Refe.fasta"
	REFCOPIED="HBV_Refe.fasta"
   elif [[ "$Virus" == "HCV" ]] 
   then
	echo "HCV"
	REFERENCEFILEALL=$CROSSCDIR"/General/HCV_Refe.fasta" 
	REFERENCEFILE==$CROSSCDIR"/General/HCV_Refe.fasta"
	REFCOPIED="HCV_Refe.fasta"
   else
        echo "** "
        echo "**  Please, check your second argument, the path is'nt correct. **"
	exit
fi

if [[ "$SEQlength" == "" ]] 
   then
	$SEQlength=100
   elif [[ "$SEQlength" -gt 500 ]] 
   then
	echo "** Minimum sequence length is $SEQlength, its value should be between 100 and 500."
	exit
   elif [[ "$SEQlength" -lt 100 ]] 
   then
	echo "** Minimum sequence length is $SEQlength, its value should be between 100 and 500."
	exit
fi

# Create directory
#.................

cd $RUNDIR
echo "**  Rundir   **"
pwd

#.......................................... functions ..............................................#
decompr (){ 		# decompress gz
	local filein=$1 
	echo $filein 
	gunzip $filein
}
decomprbz (){  		# decompress bz2
	local filein=$1 
	echo $filein 
	bzip2 -d $filein
}

linkf (){ 
	local filein=$1 
	ln -s $filein . 
}

samf (){ 
	local fileR1=$1 
	echo $fileR1  
	fileR2=`echo $fileR1 | sed s/_R1/_R2/` 
	echo $fileR2   
	samfile=${fileR1%%_R1*.fastq}_bwa.sam  
	echo " ** sam file= $samfile "  
	#bwa mem HXB2R.fna $fileR1 $fileR2 > $samfile 2>> log.txt
	bwa mem $REFCOPIED $fileR1 $fileR2 > $samfile 2>> log.txt
}

sortedf (){    #*_bwa.sam
	local fileR1=$1 
	echo $file  
	fileout=${file%%.sam}_sorted  
	echo "** sorted bam file= $fileout "  
	samtools view -Shb $file | samtools sort - $fileout 2>> log.txt   
}

depthf (){ 
	local fileb=$1  # *_bwa_sorted.bam
	echo $fileb  
	fileout=${fileb%%_bwa_sorted.bam}_bwa_depth.txt  
    	echo "** depth file= $fileout"   
	samtools depth $fileb >  $fileout    # K03455|HIVHXB2CG	2114	50  =>  id   posición   coverage
}

decontam (){
        local file1=$1  # *_R1*.fastq
	fileR1=${file1##$RUNDIR/crossCont/}
        echo " ** decontam file R1= $fileR1" 
        fileR2=`echo $fileR1 | sed s/_R1/_R2/` 
        echo " ** decontam file R2= $fileR2"  
        sample=${fileR1%%_R1*.fastq} 
        echo $sample  
        bbsplit.sh build=1 in1=$fileR1 in2=$fileR2 basename=$RUNDIR/crossCont/desContam/$sample/out_%.fq outu1=$RUNDIR/crossCont/desContam/$sample/clean1.fq outu2=$RUNDIR/crossCont/desContam/$sample/clean2.fq scafstats=$RUNDIR/crossCont/desContam/$sample/scafstats.txt 2>> log.txt
   
}

sbtypref (){
	local file1=$1  # *_R1*.fastq
	fileR1=${file1##$RUNDIR/crossCont/}
	sample=${fileR1%%_R1*.fastq}  
	fileout=$sample"_clean.fq"  
	echo " ** decontam file clean= $fileout"  
	cp $RUNDIR/crossCont/desContam/$sample/out_Alltypes_Refe.fq ./$fileout  
   #. move to dir
	cd $RUNDIR/crossCont/desContam/$sample/
# Will generate subtype reference with the majority sequence  => xxx_sbtypeRef.fasta

	intid=`head -n 2 scafstats.txt  | grep -v '#' | awk '{print $1}'` 
	fileSbref=$sample"_sbtypeRef.fasta"  
	fgrep -A 1 $intid $REFERENCEFILEALL > $RUNDIR/crossCont/SubtypeRefs/$fileSbref 
	cd $RUNDIR/crossCont/ 
	pwd 
}

sbtypind (){
	local fileRef=$1  # SubtypeRefs/*_sbtypeRef.fasta
	echo $file  
	bwa index -a is $fileRef 2>> log.txt
}


alignfq (){
	local filef=$1  # $RUNDIR/crossCont/*_clean.fq
	echo " ** clean.fq= $filef " 
	fileff=${filef##$RUNDIR/crossCont/}
	sample=${fileff%%_clean.fq}  
	fileref=$RUNDIR/crossCont/SubtypeRefs/$sample"_sbtypeRef.fasta" 
	echo " ** reference= $fileref "
	outfile=$sample"_NRclean.sam" 
	bwa mem $fileref $filef > $outfile 2>> log.txt
}

sortedsam (){
	local filesam=$1  # *_NRclean.sam
	echo $filesam  
	fileout=${filesam%%.sam}_sorted  
	echo $fileout  
	samtools view -Shb $filesam | samtools sort - $fileout 2>> log.txt  
}

consens (){
	local filebam=$1  # $RUNDIR/crossCont/*_NRclean_sorted.bam
	echo " ** _sorted.bam= $filebam"  
	filebamm=${filebam##$RUNDIR/crossCont/}
	sample=${filebamm%%_NRclean_sorted.bam}  
	fileout=$sample"_cons.fastq" 
	echo " ** output file: $fileout "   
	fileref=$RUNDIR"/crossCont/SubtypeRefs/"$sample"_sbtypeRef.fasta" 
	samtools mpileup -E -uf $fileref $filebam | bcftools call -c | vcfutils.pl vcf2fq > $fileout 2>> log.txt 
}

fqfasta (){
	local fileq=$1  # *_cons.fastq 
	echo " ** _cons.fastq= $fileq"  
	fileout=${fileq%%.fastq}.fasta  
	echo $fileout 
	seqtk seq -a $fileq > $fileout   
}

depthft (){
	local fileb=$1  # *_NRclean_sorted.bam
	echo $fileb  
	fileout=${fileb%%_NRclean_sorted.bam}_depth.txt          
	echo $fileout   
	samtools depth $fileb >  $fileout 2>> log.txt 
}

removel (){
	local filec=$1  # _cons.fasta
	echo $filec  
	filedp=${filec%%_cons.fasta}_depth.txt 
	echo " ** depth file: $filedp " 
	fileout=${filec%%_cons.fasta}clean_consdp.fasta 
	echo " ** output file: $fileout "  				
 	perl $PIPELINEDIR/filterByDepth.pl -i $filec -d $filedp > $fileout 
}

lastclassify (){     # $RUNDIR/crossCont/*_clean.fq
	local fileq=$1    
	echo " ** clean file: $fileq" 
	fileqq=${fileq##$RUNDIR/crossCont/}
	sample=${fileqq%%_clean.fq} 
	bbsplit.sh build=1 in=$fileq basename=./finalDir/$sample/out_%.fq outu=./finalDir/$sample/clean.fq 2>> log.txt 
}


#.............................    end functions   .............................................#
#..............................................................................................#

if [[ ! -d $RUNDIR/crossCont ]]
    then
    echo "** Creating crossCont Directory." 
    mkdir $RUNDIR/crossCont 
fi


# Decompress RawData files  					
#.........................
echo "** Check and decompress RawData files."
echo ""
#..............................................................................................#

cd $RUNDIR/RawData/ 
pwd 
for file in `ls`; 
do echo $file;
	if [ "${file##*.}" = "gz" ]; 
		then
	        	echo "** decompressing .gz files."
			decompr $file 
		fi
	if [ "${file##*.}" = "bz2" ]; 
		then
	        	echo "** decompressing .bz2 files."
			decomprbz $file
		fi
	if [ "${file##*.}" = "fastq" ]; 
		then
	        	echo "** There are .fastq files."
		fi
done
wait

# Check number of sequences from R1 and R2 fastq files.  					
#..............................................................................................#
echo "** Check number of sequences from R1 and R2 fastq files."
echo ""

for file in $RUNDIR/RawData/*.fastq;
do
	if [[ -e $file || ! -s $file ]]
	then
		echo $file 
		seq=`cat $file | wc -l` 
		nseq=$((seq / 4)) 
		namefile=${file#$RUNDIR/RawData/} 
		echo -e $namefile"\t"$nseq >> $RUNDIR/RawData/Nseq.txt
	else
		echo "* Err1, Checking number of sequences, /RawData/*.fastq not found or empty"
		exit
	fi
done

#for file in $RUNDIR/RawData/*.fastq; do echo $file; seq=`cat $file | wc -l`; nseq=$((seq / 4)); namefile=${file#$RUNDIR/RawData/}; echo -e $namefile"\t"$nseq >> $RUNDIR/RawData/Nseq.txt; done

fname="$RUNDIR/RawData/Nseq.txt"  
exec<$fname
while read line
do
                f1=`echo $line | awk -F" " '{print $1}'`   
                file1=`echo $f1 | awk -F"_" '{print $1}'` 
                nseq1=`echo $line | awk -F" " '{print $2}'` 
                read line
                f2=`echo $line | awk -F" " '{print $1}'`   
                file2=`echo $f2 | awk -F"_" '{print $1}'`
                nseq2=`echo $line | awk -F" " '{print $2}'`   
		if [[ "$file1" -eq "$file2" && "$nseq1" -ne "$nseq2" ]]  
                then
                        echo "** $file1 R1, R2 have different number of sequences" 
			if [[ ! -d $RUNDIR/RawData/diffR1R2NSeq ]]
    			then
   				 echo "** Creating /RawData/diffR1R2NSeq Directory." 
    				 mkdir $RUNDIR/RawData/diffR1R2NSeq 
			fi
			mv $f1 $f2 $RUNDIR/RawData/diffR1R2NSeq/      			
                fi
done

#..............................................................................................#
#..............................................................................................#

cd $RUNDIR/crossCont/ 
pwd 

# Link RawData files     					
#..............................................................................................#
echo "** Link files from RawData."
echo ""
for file in $RUNDIR/RawData/*.fastq; do echo $file; linkf $file & done 
wait 

# Align with HXB2R  					
#..............................................................................................#
echo "** Align to get coverage and filter if coverage < 500"
echo ""

cp $REFERENCEFILE $RUNDIR/crossCont/     
bwa index -a is $REFCOPIED 2> log.txt

for file in $RUNDIR/crossCont/*_R1*.fastq; 
do  
	if [[ ! -e $file || ! -s $file ]]  # file is empty
	then     	
		echo "* Err2, Decontaminating, $file not found or is empty."
		exit
	fi
done
for fileR1 in $RUNDIR/crossCont/*_R1*.fastq; do samf $fileR1 & done 
wait 

# Transform .sam into sorted.bam
echo "** Transform .sam into sorted.bam." 
echo "" 

for file in $RUNDIR/crossCont/*_bwa.sam
do
	if [[ ! -e $file || ! -s $file ]]
	then
		echo "* Err3, Transforming .sam into sorted.bam, $file not found or is empty."
		exit
	fi
done
for file in $RUNDIR/crossCont/*_bwa.sam;do sortedf $file & done 
wait

# Generate depth files.
#..............................................................................................#
echo "** Generate xxx_bwa_depth.txt files." 
echo ""

for file in $RUNDIR/crossCont/*_bwa_sorted.bam
do
	if [[ ! -e $file || ! -s $file ]]
	then
		echo "* Err4, Generate xxx_bwa_depth.txt files, $file not found or is empty."
		exit
	fi
done
for file in $RUNDIR/crossCont/*_bwa_sorted.bam; do depthf $file & done 
wait 

# Generate meanCov.txt.
#..............................................................................................#
echo "** Generate meanCov.txt." 
echo "" 

for file in $RUNDIR/crossCont/*_bwa_depth.txt
do
	if [[ ! -e $file || ! -s $file ]]
	then
		echo "* Err5, Generating meanCov.txt, $file not found or is empty."
		exit
	fi
done
for file in $RUNDIR/crossCont/*_bwa_depth.txt; do prom=`cat $file | awk '{ sum+=$3; n++} END {print sum/n}'`;echo -e "$file\t$prom" >> meanCov.txt;done

# Generate "filesWithLowCoverage.txt", sample file's names with mean coverage lower than 500.
#..............................................................................................#
echo "** Generate filesWithLowCoverage.txt, samples file's names with mean coverage lower than 500." 
echo "" 
cat meanCov.txt |  awk '{if ($2 < 500) print $1"\t"$2 } END {}'>> $RUNDIR/crossCont/filesWithLowCoverage.txt
# 3000

# Check low coverage files: lowcov.
echo "** Check low coverage files." 
echo "" 
if [[ -s filesWithLowCoverage.txt ]]      	
	then echo "There are files with low coverage."  
	if [[ ! -d $RUNDIR/crossCont/lowcov ]] 
    	then
    		echo "** Creating lowcov Directory." 
    		mkdir $RUNDIR/crossCont/lowcov 
	fi
	filename="./filesWithLowCoverage.txt"
	exec<$filename 
	while read line
   	do
		file=`echo $line | awk -F"-" '{print $1}'`   
		namef=${file%%_bwa_depth.txt}_R1*.fastq    
		mv $namef $RUNDIR/crossCont/lowcov/                                                                 
	done
	else echo "** There are no files with low coverage."  
fi  


rm *_bwa.sam 
rm *_bwa_sorted.bam 

####.............................................................................................####


# Decontaminate, 
# will generate xxx_clean.fq files, with the sequences that mapped to Alltypes_Refe.fasta.
#..............................................................................................#

if [[ ! -d $RUNDIR/crossCont/desContam ]] 
    then
    echo "** Creating desContam Directory." 
    mkdir $RUNDIR/crossCont/desContam 
fi

if [[ ! -d $RUNDIR/crossCont/SubtypeRefs ]] 
    then
    echo "** Creating SubtypeRefs Directory ." 
    mkdir $RUNDIR/crossCont/SubtypeRefs 
fi

echo "** Decontaminate =>  desContam/sample/xxxx_clean.fq" 
echo "" 
bbsplit.sh build=1 ref_Alltypes_Refe=$REFERENCEFILEALL 2>> log.txt

for fileR1 in $RUNDIR/crossCont/*_R1*.fastq; 
do  
	if [[ ! -e $fileR1 || ! -s $fileR1 ]]
	then     	
		echo "* Err6, Decontaminating, $fileR1 not found or is empty."
		exit
	fi
done
for fileR1 in $RUNDIR/crossCont/*_R1*.fastq; do decontam $fileR1 & done 
wait
for fileR1 in $RUNDIR/crossCont/*_R1*.fastq; do sbtypref $fileR1 & done
wait

####.............................................................................................####

# Generate index  
#..............................................................................................#           
echo "** Generate _sbtypeRef.fasta index"
echo "" 

for file in $RUNDIR/crossCont/SubtypeRefs/*_sbtypeRef.fasta ;
do  
	if [[ ! -e $file || ! -s $file ]]
	then     	
		echo "* Err7, Generating _sbtypeRef.fasta, $file not found or is empty."
		exit
	fi
done
for file in $RUNDIR/crossCont/SubtypeRefs/*_sbtypeRef.fasta ; do sbtypind $file & done 
wait

# Align
#..............................................................................................#
echo "** Align   => xxx_NRclean.sam." 
echo "" 

for file in $RUNDIR/crossCont/*_clean.fq ;
do  
	if [[ ! -e $file || ! -s $file ]]
	then     	
		echo "* Err8, ALigning _clean.fq, $file not found or is empty."
		exit
	fi
done
for file in $RUNDIR/crossCont/*_clean.fq; do alignfq $file & done 
wait
	
# Transform .sam into sorted.bam
echo "** Transform .sam into sorted.bam."
echo ""
#module load SAMTOOLS

for filecls in $RUNDIR/crossCont/*_NRclean.sam
do
	if [[ ! -e $filecls || ! -s $filecls ]]
	then
		echo "* Err9, Transforming sam into sorted.bam, $filecls not found or is empty."
		exit
	fi
done
for file in $RUNDIR/crossCont/*_NRclean.sam; do sortedsam $file & done
wait

# Generate consensus (from each sample) with Samtools   => xxxxx_cons.fastq
#..............................................................................................#
echo "** Generate consensus sample_cons.fastq." 
echo "" 

for fileclb in $RUNDIR/crossCont/*_NRclean_sorted.bam
do
	if [[ ! -e $fileclb || ! -s $fileclb ]]
	then
		echo "* Err10, Generating consensus sample_cons.fastq, $fileclb not found or is empty."
		exit
	fi
done
for file in $RUNDIR/crossCont/*_NRclean_sorted.bam; do consens $file & done    
wait 

# fastq => fasta
echo "** Transform fastq => fasta" 
echo "" 

for file in $RUNDIR/crossCont/*_cons.fastq ;
do  
	if [[ ! -e $file || ! -s $file ]]
	then     	
		echo "* Err11, Transforming fastq => fasta, $file not found or is empty."
		exit
	fi
done
for file in $RUNDIR/crossCont/*_cons.fastq ; do fqfasta $file & done    
wait 

# Generate depth files.
#..............................................................................................#
echo "** Generate xxx_depth.txt files." 
echo "" 

for filesb in $RUNDIR/crossCont/*_NRclean_sorted.bam
do  
	if [[ ! -e $filesb || ! -s $filesb ]]
	then
		echo "* Err12, Generating xxx_depth.txt files, $filesb not found or is empty."
		exit
	fi
done
for file in $RUNDIR/crossCont/*_NRclean_sorted.bam; do depthft $file & done   
wait 

# Remove low coverage positions.
#..............................................................................................#
echo "** Remove low coverage positions and generate xxxclean_consdp.fasta." 
echo "" 

for fcf in $RUNDIR/crossCont/*_cons.fasta
do  
	if [[ ! -e $fcf || ! -s $fcf ]]
	then
		echo "* Err13, Removig low coverage positions, $fcf not found or is empty."
		exit
	fi
done
for file in $RUNDIR/crossCont/*_cons.fasta; do removel $file & done   
wait 

# Classify with bbsplit
#..............................................................................................#

if [[ ! -d $RUNDIR/crossCont/finalDir ]]
    then
    echo "** Creating finalDir Directory" 
    mkdir $RUNDIR/crossCont/finalDir
fi

# create reference for bbsplit
for file in $RUNDIR/crossCont/*_consdp.fasta ;
do  
	if [[ ! -e $file || ! -s $file ]]
	then     	
		echo "* Err14, Creating reference for bbsplit, $file not found or is empty."
		exit
	fi
done

ref="";for file in $RUNDIR/crossCont/*_consdp.fasta;do ref=$ref","$file; echo $ref; done

# delete first "," from $ref
reference=`echo $ref |sed 's/^,//'` 
echo "ref= $reference" 

echo "** Classify with bbsplit" 
echo "" 
rm -r ref 
# create ref index
bbsplit.sh build=1 ref=$reference 2>> log.txt


for fcl in $RUNDIR/crossCont/*_clean.fq
do
	if [[ ! -e $fcl || ! -s $fcl ]]
	then
		echo "* Err15, Last classifying with bbsplit, $fcl not found or is empty."
		exit
	fi
done
for file in $RUNDIR/crossCont/*_clean.fq; do lastclassify $file & done  
wait 

#..............................................................................................#

cd $RUNDIR/crossCont/finalDir/

# Obtain the number of reads  and length from finalDir files
#..............................................................................................#
echo "** Obtain the number of reads  and length from finalDir files and write CcOutput.txt" 
echo "" 
pwd 

finalDirPath=$RUNDIR/crossCont/finalDir/
$PIPELINEDIR/results.sh $finalDirPath > $RUNDIR/crossCont/CcOutput.txt 
cd $RUNDIR/crossCont/ 

echo "** Generate CcEnd.csv." 
echo "" 
pwd 
perl $PIPELINEDIR/results.pl $RUNDIR/crossCont/CcOutput.txt $SEQlength > $RUNDIR/crossCont/CcEnd.csv 

echo "** Remove files." 
echo "" 

rm -r $RUNDIR/crossCont/SubtypeRefs;
rm -r $RUNDIR/crossCont/desContam;
rm -r $RUNDIR/crossCont/ref;
rm $RUNDIR/crossCont/*.sam;
rm $RUNDIR/crossCont/*.bam;

plotname=$RUNDIR"/crossCont/"$PROJECT"_CrosscPlot.pdf"

if [[ ! -e CcEnd.csv  || ! -s CcEnd.csv  ]]
then
	echo  "* Err16, CcEnd.csv, final file not found or is empty"	
	exit
fi

Rscript $PIPELINEDIR/cContFinal.R $RUNDIR"/crossCont/CcEnd.csv" $plotname

echo ""
echo "Plot = $plotname"
echo "End"

