#!/usr/bin/perl


# Title:        Fusionminer (wrapper)
#
# Function:     This runs the Fusionminer pipeline.
# 
# Author:       Gavin Oliver
# 
# Date:         20 Aug 2009
# 
# Version:      V1.0 
#
# Copyright (C) <2010>  <Gavin Oliver>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


use strict;
use warnings;
use Getopt::Std;

my $start = time(); #start time

my $parameters; # Fusionminer parameter file

my $bfinfile; 	# Initial input file
my $bfcommand="./blastfilter.pl "; 	# blastfilter command
my $cfcommand="./chimfind.pl "; 	# chimfind command
my $cftcommand="./chimfilt.pl "; 	# chimfilt command
my $cftscommand="./chimfiltSING.pl ";	# chimfiltsing command
my $sscommand="./splicescan.pl "; 	# splicescan command
my $ssscommand="./splicescanSING.pl ";	# splicescansing command
my $cccommand="./chimclust.pl "; 	# chimclust command
my $fccommand="./framecheck.pl ";    	# framecheck command
my $fname; # Standard naming prefix for all output files
my $outdir; # output directory
my %params; # hash for storage of parameter values
my $version = "1.0";
# Get command line options and arguments
my %opts;
&getopt('i:p:d:h',\%opts);

# If help flag given print usage
if (exists $opts{h}) {
        &usage;
        exit;
}

# Unless paramter file is specified and exists print usage and exit
if (defined $opts{p}) {
	$parameters = $opts{p};
	if (! -e $parameters) {
		print "Parameter file does not exist!";
		exit;
	}
}
else {
	&usage;
	exit;
}

# Unless bf input file specified print usage and exit
if (defined $opts{i}) {
        $bfinfile = $opts{i};
	# remove path component
	$fname =$bfinfile;
	$fname =~ s/(.*)\///;
	$fname =~ s/\..*//;
}
else {
        &usage;
        exit;
}

# check for output directory
if (defined $opts{d}) {
    $outdir=$opts{d};
    if (! -e $outdir){
        print "Problem with output directory $outdir\n";
        exit;
    }
    print "Output directory = $outdir\n";
}

# Extract parameters from file to hash
open (PARAMS, $parameters) or die $!;

while (<PARAMS>) {
	chomp;
	next if (/^$/ || /^#/);
	my @linearray = split (/\s+/, $_);
	my $flag = $linearray[0];
	my $value = $linearray[1];
	$params{$flag}=$value;	
}
close PARAMS;

# Extract parameters from hash and append to script run commands
if (defined $params{-a}) {
	$bfcommand.="-l $params{-a} ";
}
if (defined $params{-b}) {
        $bfcommand.="-s $params{-b} ";  
}	
if (defined $params{-c}) {
        $bfcommand.="-p $params{-c} ";  
}
if (defined $params{-d}) {
        $cfcommand.="-o $params{-d} ";  
}
if (defined $params{-e}) {
        $cftcommand.="-l $params{-e} ";  
	$cftscommand.="-l $params{-e} ";
	$ssscommand.="-l $params{-e} ";
}
if (defined $params{-f}) {
        $cftcommand.="-p $params{-f} ";
        $cftscommand.="-p $params{-f} ";
        $ssscommand.="-p $params{-f} ";
}
if (defined $params{-g}) {
        $cftcommand.="-o $params{-g} ";
        $cccommand.="-o $params{-g} ";
        $ssscommand.="-o $params{-g} ";
}
if (defined $params{-h}) {
        $sscommand.="-t $params{-h} ";
        $cftscommand.="-t $params{-h} ";
        $ssscommand.="-t $params{-h} ";
}
if (defined $params{-i}) {
        $sscommand.="-c $params{-i} ";
        $cftscommand.="-c $params{-i} ";
        $ssscommand.="-c $params{-i} ";
}
if (defined $params{-j}) {
	$fccommand.="-c $params{-j} ";
}

# Check input file exists 
unless (-e $bfinfile) {
	die "Input file $bfinfile does not exist!\n";
}

# Check genomic coordinates file exists
unless (-e $params{-i}) {
	die "Genomic coordinates file $params{-i} does not exist!\n";
}

# Run blastfilter
print "\nRunning blastfilter.pl\n";
system "perl $bfcommand-i $bfinfile -o $fname.filt";

# Sort blastfilter output
print "\nSorting blastfilter.pl output\n";
system "sort $fname.filt > $fname.filt.srted";

# Run chimfind on sorted blastfilter ouput
print "\nRunning chimfind.pl\n";
system "perl $cfcommand-i $fname.filt.srted";

# Run chimfilt on multi chromosome accepted candidate file
print "\nRunning chimfilt.pl on interchromosomal candidates\n";
system "perl $cftcommand-i $fname-multi.accept";

# Run chimfiltSING on single chromosome accepted candidate file
print "\nRunning chimfiltSING.pl on intrachromosomal candidates\n";
system "perl $cftscommand-i $fname-clsing.accept";

# Run splicescan on multi chromosome accepted candidate file
print "\nRunning splicescan.pl on interchromosomal candidates\n";
system "perl $sscommand-i $fname-multi.accept";

# Run splicescanSING on single chromosome accepted candidate file
print "\nRunning splicescanSING.pl on intrachromosomal candidates\n";
system "perl $ssscommand-i $fname-sing.accept";

# Run chimclust on multi chromosome accepted candidate file
print "\nRunning chimclust.pl on interchromosomal candidates\n";
system "perl $cccommand-i $fname-multi.accept";

# Run chimclust on single chromosome accepted candidate file
print "\nRunning chimclust.pl on intrachromosomal candidates\n";
system "perl $cccommand-i $fname-sing.accept";

# Run framecheck on  multi chromosome accepted and clustered candidate file
print "\nRunning framecheck.pl on clustered interchromosomal candidates\n";
system "perl $fccommand-i $fname-multi.accept.clustered -r $fname.filt.srted";

# Run framecheck on  single chromosome accepted and clustered candidate file
print "\nRunning framecheck.pl on clustered intrachromosomal candidates\n";
system "perl $fccommand-i $fname-sing.accept.clustered -r $fname.filt.srted";

# Create rejected/accepted directories and move files to them
system "mkdir Accepted";

system "mkdir Archived";

system "mv *framechecked Accepted";

#system "mv $fname-sing.accept Accepted";
#system "mv $fname-sing.accno Accepted";
#system "mv $fname-sing.accept.clustered Accepted";

#system "mv $fname-multi.accept Accepted";
#system "mv $fname-multi.accno Accepted";
#system "mv $fname-multi.accept.clustered Accepted";

system "mv $fname* Archived";

# Output run time
my ($user,$system,$cuser,$csystem) = times;
print "\n** Processing finished.  User time=$user s **\n";

# Sub to print help screen
sub usage {

    print "\n\t$0 version $version\n";
    print "\tby G.R. Oliver and A. Emerson\n";
    print "\nUsage:\n";
    print "\t $0 -i blastoutput file -p parameter file\n";
    print "\t e.g $0 -i MyBlastResults.table -p fm.parameters \n\n";
    print "Other options\n";
    print "\t -d <output-dir>\n";
    print "\t -h <print this screen>\n\n";
    print "NB: Input files should be tabular BLAST/BLAT output (-m 8 / -output=blast8) and Fusionminer format parameter file\n";

}

# Calculate and print total running time
my $end = time();
my $mins = ($end-$start)/60;
$mins = sprintf("%.2f",$mins); 
print "\nTotal time taken for analysis = $mins minutes\n\n";
