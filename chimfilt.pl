#!/usr/bin/perl

# Title:        Chimfilt (Part of FusionMiner) 
#
# Function:     Perl program to take files produced by chimfind and remove previously accepted chimeric candidates
#               when sequence alignments on different chromosomes overlap significantly or have large gaps.
#               Alignments can be rejected based on length or %id by specifying a value with the -l and -p flags
#               Multiple alignments to a single chromosome are clustered and %ID averaged in the output file
#               Output also outputs segments in 5->3 prime order for each candidate
# 
# Author:       Gavin Oliver
# 
# Date:         25 Feb 2009
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
#use warnings;
use Getopt::Std;


my $version = "25-feb-2009 1.0"; # version number
my %opts=(); # hash for storage of command line options
my $input; # input file
my $thresh; # overlap/gap threshold
my $outdir="."; # output directory
my $id; # %id filter threshold
my $length; # alignment length filter threshold

# Set command line option flags
getopts('ho:i:d:p:l:', \%opts);

# if help requested print and exit
if (exists $opts{h}) {
    &usage;
    exit;
}

# check for output directory and redefine if necessary - quit if requested dir doesn't exist
if (defined $opts{d}) {
    $outdir=$opts{d};
    if (! -e $outdir){
        print "Problem with output directory $outdir\n";
        exit;
    }
    print "Output directory = $outdir\n";
}

# Quit unless input file is defined on command line
if (exists $opts{i}) {
	$input = $opts{i};
}
	
else {&usage;
	print "\n\nInput filename must be specified on command line!\n\n";
	exit;}

# If overlap defined on command line accept value otherwise use a default of 10
if (exists $opts{o}) {
	$thresh = $opts{o};
}

else {$thresh=10;}

# If %id defined on command line accept value otherwise set to 0 to avoid filtering
if (exists $opts{p}) {
	$id = $opts{p};
}

else {$id =0;}

# If length defined on command line accept value otherwise set to 0 to avoid filtering
if (exists $opts{l}) {
        $length = $opts{l};
}

else {$length =0;}

# Set up ouput filenames based on input filename
my $rootname = $input;
$rootname=~s/(.*)\///;
$rootname=~s/\..*//;
my $accepted = "$outdir/$rootname.accept.tmp";
my $rejected = "$outdir/$rootname.reject.cfilt";
my $accs = "$outdir/$rootname.accno";

my %stats=(); # hash for storage of genome alignment positions
my %record=(); # hash to keep a record of process stats
my $accept=0; # accept or reject result
my $linebuffer; # string for storage of result lines for later printing
my %acchash=(); # hash for storage of accepted accessions
my %uniqacc=(); # hash for storage of unique accessions processed
my %rejhash=(); # hash for storage of unique accessions rejected

print "\n\t$0: version $version. Use -h for options\n\n";

open (INFILE, $input) or die $!;
open (ACC, ">$accs") or die $!;
open (ACCEPTED, ">$accepted") or die $!;
open (REJECTED, ">$rejected") or die $!;

print "% ID filter set to $id bases\n";
print "Length filter set to $length bases\n";
print "Overlap / Gap filter set to $thresh bases\n";
print "Reading alignments from $input\n";
print "Writing accepted sequences to $accepted\n";
print "Writing list of accepted accession numbers to $accs\n";
print "Writing newly rejected sequences to $rejected\n\n";

print ACCEPTED "#query_id chromosome %identity query_start query_end chrom_start chrom_end\n";
print REJECTED "#query_id chromosome %identity query_start query_end chrom_start chrom_end\n";

while (<INFILE>) {
	
	$record{nlines}++;

	# calculate overlap/gap or lack of - only does this when a blank line is reached i.e. between results
	if (/^$/) {

		# Set up keys for stat hash access
		my ($a, $b) = keys %stats;
		
		# If specified on command line, first check that %ids are sufficient otherwise reject
		if ($id) {
			my ($pid1, $pid2) = ($stats{$a}[2], $stats{$b}[2]);
			if (($pid1 < $id) || ($pid2 < $id)) {
				print REJECTED "$linebuffer\n";
				# Write rejected accessions to a hash for stats
				my @linearray = split(/\s+/, $linebuffer);
				my $acc = $linearray[0];
				# Remove asterisks so we only have unique entries
				$acc=~s/\*//;
				$rejhash{$acc}++;
				# Reset variables for processing next result
				$linebuffer="";
				%stats=();
				$accept = 0;
				next;
			}
		}

		# Extract alignment positions 		
		my ($mina, $minb, $maxa, $maxb) = ($stats{$a}[0], $stats{$b}[0], $stats{$a}[1], $stats{$b}[1]);
		
		# If specified on command line, first check that alignment lengths are sufficient otherwise reject
		if ($length) {
			
			# Extract maximum alignment lengths
			my $maxrangea = $stats{$a}[3];
			my $maxrangeb = $stats{$b}[3];

			if (($maxrangea < $length) || ($maxrangeb < $length)) {
				print REJECTED "$linebuffer\n";
				# Write rejected accessions to a hash for stats
				my @linearray = split(/\s+/, $linebuffer);
				my $acc = $linearray[0];
				# Remove asterisks so we only have unique entries
				$acc=~s/\*//;
				$rejhash{$acc}++;
				# Reset variables for processing next result
				$linebuffer="";
				%stats=();
				$accept = 0;
				next;
			}
		}

		# if sequence a's 5` end shows overlap with b
		if ($mina > $minb && $mina <= $maxb) {
			# but is not totally contained within it
			if ($maxa > $maxb) {
				my $overlap = $maxb - $mina;
				if ($overlap < $thresh) {$accept=1;}
				else {$accept=0;}
			}
			
			else {$accept=0;}
		}
		
		# if sequence a aligns 3` of sequence b
		elsif ($mina > $minb && $mina > $maxb) {
			my $gap = $mina - $maxb;
			# but the gap between them is 10 or less
			if ($gap <= $thresh) {$accept=1;}
			else {$accept=0;}
		}
	
		# # if sequence b's 3` end shows overlap with a
		elsif ($minb > $mina && $minb <= $maxa) {
                        # but is not totally contained within it
			if ($maxb > $maxa) {
                                my $overlap = $maxa - $minb;
                                if ($overlap < $thresh) {$accept=1;}
				else {$accept=0;}
			} 
	
			else {$accept=0;}
		}
		
		# if sequence b aligns 3` of sequence a
		elsif ($minb > $mina && $minb > $maxa) {
                        my $gap = $minb - $maxa;
                        # but the gap between them is 10 or less
			if ($gap <= $thresh) {$accept=1;}
			else {$accept=0;}
                }
		
		
		if ($accept) {
			# Create output string for chromosome a (multiple lines clustered)
			my $outpid;
			my $outputa = $stats{$a}[9]." ".$a." ";
			if ($stats{$a}[8]>1) {$outpid = ($stats{$a}[7])/($stats{$a}[6]);}
			else {$outpid = $stats{$a}[2];}
			$outpid = sprintf "%.2f", $outpid;
			$outputa.=$outpid." ".$mina." ".$maxa." ".$stats{$a}[4]." ".$stats{$a}[5]."\n";
			
			# Create output string for chromosome b (multiple lines clustered)
                        my $outputb = $stats{$b}[9]." ".$b." ";
                        if ($stats{$b}[8]>1) {$outpid = ($stats{$b}[7])/($stats{$b}[6]);}
                        else {$outpid = $stats{$b}[2];}
			$outpid = sprintf "%.2f", $outpid;
                        $outputb.=$outpid." ".$minb." ".$maxb." ".$stats{$b}[4]." ".$stats{$b}[5]."\n";
			
			# Ensure 5 prime segment output first
			if ($maxa > $maxb) {print ACCEPTED "$outputb"."$outputa"."\n";}
			else {print ACCEPTED "$outputa"."$outputb"."\n";}
	
			# Write newly accepted accessions to a hash for printing to acc file later
                	my @linearray = split(/\s+/, $linebuffer);
                	my $acc = $linearray[0];
                	# Remove asterisks so we only have unique entries
			$acc=~s/\*//;
                	$acchash{$acc}++;
		}

		else {
			print REJECTED "$linebuffer\n";
			# Write rejected accessions to a hash for stats
                        my @linearray = split(/\s+/, $linebuffer);
                        my $acc = $linearray[0];
                        # Remove asterisks so we only have unique entries
                        $acc=~s/\*//;
                        $rejhash{$acc}++;
		}

	
		# Reset variables for processing next result	
		$linebuffer="";
		%stats=();
		$accept = 0;
        }

	# Prior to reaching a blank line, extract result from file for processing
	else {	
		next if (/^#/); # ignore comment lines
		# Concatenate lines into a string for printing to accept or reject file
		$linebuffer.=$_;
		chomp;
		my @linearray = split (/\s+/); # split results into an array
		# Check input file format
		unless (%uniqacc) {
			unless ($#linearray ==6) {
				die "***Quitting!*** Incorrect input file format supplied to chimfilt!\n\n";
			}
		}
		# Assign their own variables for ease of use
		my ($acc, $chr, $pid, $min, $max, $cmin, $cmax)=
		($linearray[0], $linearray[1], $linearray[2], $linearray[3], $linearray[4], $linearray[5], $linearray[6]);
		# Clean accession and increment unique accession counter
		$acc=~s/\*//;
		$uniqacc{$acc}++;
		# Store max and min hit positions for each of the two chromosomes
		unless (defined $stats{$chr}[0] && $stats{$chr}[0] < $min) {$stats{$chr}[0] = $min; $stats{$chr}[4] = $cmin;}
		unless (defined $stats{$chr}[1] && $stats{$chr}[1] > $max) {$stats{$chr}[1] = $max; $stats{$chr}[5] = $cmax;}
		# Store minimum %id for each chromosomal alignment
		$stats{$chr}[2] = $pid unless (defined $stats{$chr}[2] && $stats{$chr}[2] < $pid);
		# Store maximum alignment length for each chromosomal alignment
		my $range = ($max - $min)+1;
		$stats{$chr}[3] = $range unless (defined $stats{$chr}[3] && $stats{$chr}[3] > $range);
		# Store total number of bases aligning
		$stats{$chr}[6] += $range;
		# Keep total of %id x bases to caculate overall %id where there are multiple alignments
		my $avid = $pid*$range;
		$stats{$chr}[7] += $avid;
		# Keep record of number of alignmnets for a given chromosome
		$stats{$chr}[8]++;
		# Store accession 
		$stats{$chr}[9]=$acc;
	}
}

close INFILE;
close ACCEPTED;
close REJECTED;

# Overwrite old accept file with new version
print "Replacing $input file with new version ($accepted). (Original will be saved as .old.cfilt)\n";
system "cp $input $outdir/$rootname.accept.old.cfilt";
system "mv $accepted $outdir/$rootname.accept";

# Overwrite old acc file with new version
foreach my $key (keys %acchash) {
	print ACC "$key\n";
}

close ACC;

# Print processing time info.

my ($user,$system,$cuser,$csystem) = times;
print "\n** Processing finished.  User time=$user s **\n";

# Make counts of accession numbers from their hashes for stats
$record{naccno} = keys %uniqacc;
$record{nrej} = keys %rejhash;
$record{nacc} = keys %acchash;

# Print output stats
print "\n\tchimfilt statistics\n\n";
print " No. input lines processed\t\t\t:$record{nlines}\n";
print " No. distinct accession numbers\t\t\t:$record{naccno}\n";
print " No. acc numbers rejected\t\t\t:$record{nrej}\n";
print " No. acc numbers accepted\t\t\t:$record{nacc}\n\n";

# Warn if discrepancy in stats
if ($record{nacc}+$record{nrej}!=$record{naccno}) {
   print "**WARNING: Naccepted + Nrejected <> N accession numbers\n\n";
}

# Subroutine to print help to screen
sub usage {

    print "\n\t$0 version $version\n";
    print "\tby G.R. Oliver\n";
    print "\nUsage:\n";
    print "\t $0 -i input filename (.accept file created by chimfind)\n";
    print "\t e.g $0 -i example.accept\n\n";
    print "Other options\n";
    print "\t -o <overlap/gap permitted at breakpoint (default 10)>\n";
    print "\t -p <% id filter threshold (default 0)>\n";
    print "\t -l <alignment length filter threshold (default 0)>\n";
    print "\t -d <output-dir>\n";
    print "\t -h <print this screen>\n";
    print "NB: Input should have been filtered (blastfilter), sorted (UNIX sort) and processed (chimfind)\n";

}

