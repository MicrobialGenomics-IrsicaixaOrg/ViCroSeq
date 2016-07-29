#!/usr/bin/perl

# Cristina Rodríguez  2015/12/10.
# Will transform the data from CcOutput.txt  to generate a table that will be processed in cContFinal.R
#
# Samples ...................
# Samples percentage ........
# ..
# We will not include clean.fq's data in the table, but in the calculation. 
# filter by length reads passed by parameter.
# _cons.fq =>  clean_consdp.fq

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

my $inputFile=$ARGV[0];
my $seqlength=$ARGV[1];
my $Hfilein={};
my $Htotreads={};

open(INPUT,"<$inputFile") or die("cannot open $inputFile"); 
my $headers= <INPUT>;
chomp($headers);
my @head=split(",",$headers); 
print STDOUT "Samples\t";

my $nucol=scalar(@head);
for(my $i=0; $i < scalar(@head); $i++){
	print STDOUT "$head[$i]";
	if ($i != $nucol-1 )
        {
		print STDOUT "\t";
	}
}
print STDOUT "\n";

my $sa;
# Generate hash $Hfilein->{$sample_row}->{'_totalreads'} adding the reads from each contaminated sample. 
while(<INPUT>){
	my $string=$_;
	chomp($string);
	my @line=split("\t",$string);     
	if ($line[0] eq "Sample"){        
		$sa= $line[1];           
	}
	else{
		if (($line[1] >= $seqlength) && ($line[0] ne "clean.fq"))	 
		{
			$Htotreads->{$sa}->{'_totreads'}+=$line[2];   	        
			$Hfilein->{$sa}->{$line[0]}->{'_lengthreads'}=$line[1]; 	
			$Hfilein->{$sa}->{$line[0]}->{'_numreads'}=$line[2];   
		}		
	}
}

##
#.....................
# Calculate percentage and print
foreach my $sample (sort keys %$Hfilein)
{
	print STDOUT "$sample\t";        # each sample row
	my $ncol=scalar(@head);
 	for(my $i=0; $i < scalar(@head); $i++){  	
				 
		my $hdsample = "out_".$head[$i]."clean_consdp.fq"; 
		if (defined $Hfilein->{$sample}->{$hdsample}->{'_numreads'})
		{			 
		  my $perc = ($Hfilein->{$sample}->{$hdsample}->{'_numreads'}*100)/$Htotreads->{$sample}->{'_totreads'};
		  printf STDOUT "%.2f",$perc;
		  if ($i != $ncol-1 )
                  {
			print STDOUT "\t";
		  }
		}else
		{
		  printf STDOUT "%.2f",0;
                  if ($i != $ncol-1 )
                  {
			print STDOUT "\t";
		  }
		}
	}
	print STDOUT "\n";
}

close INPUT;



