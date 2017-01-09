#!/usr/bin/perl 

# author: M. Cristina Rodríguez, 2016/10/19.
# Create consensus file from "samtools mpileup  -f xxx_sbtypeRef.fasta xxx_NRclean_sorted.bam" => result file: xxx_precons.txt.
# and filter positions with cov > 6.
# 
use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

#  perl ../scripts/getConsens.pl -i xxx_precons.txt  > xxx_cons.fasta

my %opts;

# Parameters;
# # -i input xxx_precons.txt file .

getopts('i:h', \%opts);
my $fprecons = $opts{i} or die help();

open(FILEIN,"<$fprecons") or die("cannot open $fprecons");  # $fprecons
my $x =0;
my $Np=0;
my $Ncom=0;
my $NA=0;
my $NT=0;
my $NG=0;
my $NC=0;
my $Na=0;
my $Nt=0;
my $Ng=0;
my $Nc=0;
my $Nx1=0;
my $Nx2=0;
my $seqfin;
my $finalNucleot;
my $strin;
my $header;
while(<FILEIN>){ 
	chomp;
 	my @line = split("\t", $_);  
	$header=$line[0];
	my $RefNucleotide = $line[2];
	if ($line[3] > 6){   		 # coverage should be > 6
		$strin= $line[4];
	}else{ next; }
 	
# count how many '.,AGTC' are.
	my $xp = '^].';
	my @cp = $strin =~ /\Q$xp/g;
	my $countp = $#cp;
	my $xc = '^],';
	my @cc = $strin =~ /\Q$xc/g;
	my $countc = $#cc ;

	$Np=($line[4] =~ tr/\./ /);
	$Ncom=($line[4] =~ tr/,/ /);
	$NA=($line[4] =~ tr/A/ /);
	$Na=($line[4] =~ tr/a/ /);
	$NT=($line[4] =~ tr/T/ /);
	$Nt=($line[4] =~ tr/t/ /);
	$NG=($line[4] =~ tr/G/ /);
	$Ng=($line[4] =~ tr/g/ /);
	$NC=($line[4] =~ tr/C/ /);
	$Nc=($line[4] =~ tr/c/ /);
	my $Npfinal= $Np-($countp + 1);
	my $Ncomfinal = $Ncom-($countc + 1);
	my $Nxpfinal=$countp + 1;   # '^].'
	my $Nxcfinal=$countc + 1;   # '^],'

	my $Np_cfinal=$Npfinal + $Ncomfinal;
	my $NAafinal=$NA+$Na;
	my $NTtfinal=$NT+$Nt;
	my $NGgfinal=$NG+$Ng;
	my $NCcfinal=$NC+$Nc;

	my $hNn={};
	$hNn->{"X"}=$Np_cfinal;
	$hNn->{"A"}=$NAafinal;
	$hNn->{"T"}=$NTtfinal;
	$hNn->{"G"}=$NGgfinal;
	$hNn->{"C"}=$NCcfinal;

# sort per value
	my @keys = sort { $hNn->{$b} <=> $hNn->{$a} } keys(%$hNn);
	my @vals = @{$hNn}{@keys};

# Assign nucleotide
	$finalNucleot = "";
	if ($vals[0]!= $vals[1])     # exists a maximum
	{
		if ($keys[0] eq "X")
		{
			$finalNucleot = $RefNucleotide;	
		}else{
			$finalNucleot = $keys[0];
			if ($finalNucleot eq "A")    { if ($Na > $NA) {$finalNucleot = "a"}} 
			elsif ($finalNucleot eq "T") { if ($Nt > $NT) {$finalNucleot = "t"}} 
			elsif ($finalNucleot eq "G") { if ($Ng > $NG) {$finalNucleot = "g"}}
			elsif ($finalNucleot eq "C") { if ($Nc > $NC) {$finalNucleot = "c"}} 
		}
	}else{
# Assign nucleotide IUPAC codes
		if (($vals[0] == 0) && ($vals[1] == 0)){ 
				$finalNucleot = $RefNucleotide;
		}elsif ($keys[0] eq "X") { 
			if ((($RefNucleotide eq "A") && ($keys[1]eq "G")) ||  (($RefNucleotide eq "G") && ($keys[1]eq "A"))){
				$finalNucleot = "R";
			 }
			if ((($RefNucleotide eq "A") && ($keys[1]eq "T")) ||  (($RefNucleotide eq "T") && ($keys[1]eq "A"))){
				$finalNucleot = "W";
			 }
			if ((($RefNucleotide eq "A") && ($keys[1]eq "C")) ||  (($RefNucleotide eq "C") && ($keys[1]eq "A"))){
				$finalNucleot = "M";
			 }
			if ((($RefNucleotide eq "C") && ($keys[1]eq "T")) ||  (($RefNucleotide eq "T") && ($keys[1]eq "C"))){
				$finalNucleot = "Y";
			 }
			if ((($RefNucleotide eq "G") && ($keys[1]eq "C")) ||  (($RefNucleotide eq "C") && ($keys[1]eq "G"))){
				$finalNucleot = "S";
			 }
			if ((($RefNucleotide eq "G") && ($keys[1]eq "T")) ||  (($RefNucleotide eq "T") && ($keys[1]eq "G"))){
				$finalNucleot = "K";
			 }
		}elsif ($keys[1]eq "X")
		{ 
			if ((($RefNucleotide eq "A") && ($keys[0]eq "G")) ||  (($RefNucleotide eq "G") && ($keys[0]eq "A"))){
				$finalNucleot = "R";
			 }
			if ((($RefNucleotide eq "A") && ($keys[0]eq "T")) ||  (($RefNucleotide eq "T") && ($keys[0]eq "A"))){
				$finalNucleot = "W";
			 }
			if ((($RefNucleotide eq "A") && ($keys[0]eq "C")) ||  (($RefNucleotide eq "C") && ($keys[0]eq "A"))){
				$finalNucleot = "M";
			 }
			if ((($RefNucleotide eq "C") && ($keys[0]eq "T")) ||  (($RefNucleotide eq "T") && ($keys[0]eq "C"))){
				$finalNucleot = "Y";
			 }
			if ((($RefNucleotide eq "G") && ($keys[0]eq "C")) ||  (($RefNucleotide eq "C") && ($keys[0]eq "G"))){
				$finalNucleot = "S";
			 }
			if ((($RefNucleotide eq "G") && ($keys[0]eq "T")) ||  (($RefNucleotide eq "T") && ($keys[0]eq "G"))){
				$finalNucleot = "K";
			 }
		}else
		{
			if ((($keys[0] eq "A") && ($keys[1]eq "G")) ||  (($keys[0] eq "G") && ($keys[1]eq "A"))){
				$finalNucleot = "R";
			 }
			if ((($keys[0] eq "A") && ($keys[1]eq "T")) ||  (($keys[0] eq "T") && ($keys[1]eq "A"))){
				$finalNucleot = "W";
			 }
			if ((($keys[0] eq "A") && ($keys[1]eq "C")) ||  (($keys[0] eq "C") && ($keys[1]eq "A"))){
				$finalNucleot = "M";
			 }
			if ((($keys[0] eq "C") && ($keys[1]eq "T")) ||  (($keys[0] eq "T") && ($keys[1]eq "C"))){
				$finalNucleot = "Y";
			 }
			if ((($keys[0] eq "G") && ($keys[1]eq "C")) ||  (($keys[0] eq "C") && ($keys[1]eq "G"))){
				$finalNucleot = "S";
			 }
			if ((($keys[0] eq "G") && ($keys[1]eq "T")) ||  (($keys[0] eq "T") && ($keys[1]eq "G"))){
				$finalNucleot = "K";
			 }

		}
	}
	$seqfin=$seqfin.$finalNucleot;
}
print ">$header\n";
print "$seqfin\n";
#print Dumper(@row);

close FILEIN;

#... help  function .....................
sub help{
	print STDOUT "Usage:\n";
	print STDOUT "  -i xxx_precons.txt.\n";
	print STDOUT " \n";
	exit;
}

