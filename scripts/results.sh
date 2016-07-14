#!/bin/bash

# 2015/12/1   C.R.
# Obtain the number of reads  and length from the finalDir files.

finalDirpath=$1
cd $finalDirpath

listdir=""
for dir in `ls`
        do listdir=$listdir","$dir
done
ldir=`echo $listdir |sed 's/^,//'`

echo $ldir            
for dir in `ls`
        do echo -e "Sample\t$dir"
        cd $dir 
        for file  in *.fq
                do reads=`fastq-stats $file | sed -n 1p`
		read=`echo $reads | awk -F" " '{print $2}'`
                len=`fastq-stats $file | sed -n 3p`   # sequence length
                le=`echo $len | awk -F" " '{print $3}'`
                if  [[ "$read" == "reads" ]] # means there are no reads
                        then   le=0
                               read=0
                fi
                echo -e "$file\t$le\t$read"
        done
        cd ../
done

