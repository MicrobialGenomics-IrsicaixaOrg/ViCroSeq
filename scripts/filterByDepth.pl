#!/usr/bin/perl

# Cristina Rodríguez 2016/01/19
# Write a fasta file from xxx_clean_cons.fasta and depth.txt selecting opsitions with cov > 6.

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;
use File::Basename;

# Parameters;
# # -i input _cons.fasta file
# # -d depth.txt file. 

my %opts;
getopts('i:d:h', \%opts);

my $fileIN = $opts{i} or die help();
my $fdepth = $opts{d} or die help();

open(FDEPTH,"<$fdepth") or die("cannot open $fdepth"); 

my $Hdp={};
my @dp;
my $x=0;
while(<FDEPTH>){ 
	chomp;
	@dp= split("\t",$_);     # id  posición coverage
	$Hdp->{$x}->{'_cov'}=$dp[2];
	$x++;
}

#... consensus fasta....
open(FILEIN,"<$fileIN") or die("cannot open $fileIN"); 

my $header= <FILEIN>;
my @aa;
while(<FILEIN>){ 
	chomp;
	@aa = split(//,$_); 	
}

print STDOUT "$header";	
foreach my $position (sort { $a <=> $b } keys %$Hdp)  
 	{
	if ($Hdp->{$position}->{'_cov'} > 6)
		{
	        	my $pos= $position -1;  
			print STDOUT "$aa[$pos]";	
		}
}

close FDEPTH;
close FILEIN;

#............. function help.............#

sub help{
	print STDOUT "Usage:\n";
	print STDOUT "  -i input cons.fasta .\n";
	print STDOUT "  -d depth.txt file.\n";
	print STDOUT " \n";
	exit;
}











