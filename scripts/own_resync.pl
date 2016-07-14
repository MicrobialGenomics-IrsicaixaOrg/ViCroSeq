#!/usr/bin/perl -w

###########################################################
# resync.pl
# John Garbe
# June 2012
#
# Resynchronize a pair of paired-end fastq files. Read in two unsynchronized 
# fastq files and write out two synchronized fastq files. The
# synchronized files have properly paired reads, with singelton reads removed.
# 
#
###########################################################

die "USAGE: resync.pl sample1_R1.fastq sample1_R2.fastq [sample1_R1_synced.fastq] [sample1_R2_synced.fastq]\n" unless ($#ARGV >= 1 && $#ARGV <= 3);

# open up input files
open F1, "<$ARGV[0]" or die "cannot open $ARGV[0]\n";
open F2, "<$ARGV[1]" or die "cannot open $ARGV[1]\n";

# open up output files
if ($#ARGV > 1) {
    open O1, ">$ARGV[2]" or die "cannot open $ARGV[2]\n";
} else {
    open O1, ">$ARGV[0].out" or die "cannot open $ARGV[0].out\n";
}
if ($#ARGV == 3) {
    open O2, ">$ARGV[3]" or die "cannot open $ARGV[3]\n";
} else {
    open O2, ">$ARGV[1].out" or die "cannot open $ARGV[1].out\n";
}

$readid = `head -n 1 $ARGV[0]`;
chomp $readid;
if ($readid =~ /^@\S+\/[12]$/) { # @ at the start of the line followed by non-whitespace, a /, a 1 or 2, the end of the line
    $id_type = 1;
    print STDERR "Casava 1.7 read id style\n"; # TESTING
} elsif ($readid =~ /^@\S+\W[12]\S+$/) { # @ at the start of the line followed by non-whitspace, a space, a 1 or 2, non-whitespace
    $id_type = 2;
    print STDERR "Casava 1.8 read id style\n"; # TESTING
} else {
    print STDERR "Cannot determine read id style\n";
    print STDOUT "Unknwon id style: $readid\n";
    exit 1;
}


# read in the first fastq file and store in a hash
$f1readcount = 0;
while ($f1line1 = <F1>) {
    $f1readcount++;
    if ($f1readcount % 1000000 == 0) { # print progress update
	print STDERR "$f1readcount reads processed in first file\n";
    }
    $f1line2 = <F1>;
    $f1line3 = <F1>;
    $f1line4 = <F1>;

    ($id, $junk) = split /\ /, $f1line1;
    $r1{$id} = $f1line1 . $f1line2 . $f1line3 . $f1line4;
}

# read in the second fastq file, printing out both read pairs if available, otherwise tossing the reads
$f2readcount = 0;
$f1missed = $f1readcount;
$f2missed = 0;
while ($f2line1 = <F2>) {
    $f2readcount++;
    if ($f2readcount % 1000000 == 0) { # print progress update
	print STDERR "$f2readcount reads processed in second file\n";
    }
    $f2line2 = <F2>;
    $f2line3 = <F2>;
    $f2line4 = <F2>;

    ($id, $junk) = split /\ /, $f2line1;
    if (defined($r1{$id})) {
	$f1missed--;
	print O1 "$r1{$id}";
	print O2 "$f2line1";
	print O2 "$f2line2";
	print O2 "$f2line3";
	print O2 "$f2line4";
    } else {
	$f2missed++;
    }
}

print "$f1readcount reads in $ARGV[0], $f1missed without mate\n";
print "$f2readcount reads in $ARGV[1], $f2missed without mate\n";


