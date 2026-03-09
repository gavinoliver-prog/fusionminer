#!/usr/bin/perl

# Title:        Splicescan (Part of FusionMiner)
# 
# Function:     Program to filter chimeric candidates output by chimfilt.
#               Utilises genomic exon coordinates to confirm or reject
#               chimeric candidates from chimfilt.
#
# Author:       Gavin Oliver
# 
# Date:         07 Apr 2009
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

my $version = "07-apr-2009 1.0"; #version number
my %opts=(); # hash for storage of command line options
my %stats=(); # hash for storage of genome alignment positions
my $outdir="."; # output directory
my $accept = 0; # accept or reject result
my %acchash=(); # hash for storage of accepted accessions
my %uniqacc=(); # hash for storage of unique accessions processed
my %rejhash=(); # hash for storage of unique accessions rejected
my %record=(); # hash to keep a record of process stats
my $input; # input file
my $chrexons; #file containing chromosomal exon start/stop positions
my $linebuffer;# buffer for storage of input lines 
my $threshold; # threshold distance (in bp) allowed upstream/downstream of alignment start/finish to search for splice signal 
my %exons=(); # hash to store chromosomal exon start/stop positions


# Set command line option flags
getopts('ht:i:d:c:', \%opts);

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
	print "\n\nInput filename must be specified!\n\n";
	exit;
}

# Quit unless chromosome exon position file defined on command line
if (exists $opts{c}) {
	$chrexons = $opts{c};
}

else {&usage;
	print "\n\nFile containing chromosome exon start/stop sites must be specified!\n\n";
	exit;
}

# If threshold defined on command line accept value otherwise use a default of 3
if (exists $opts{t}) {
	$threshold = $opts{t};
}

else {$threshold=3;}


# Set up ouput filenames based on input filename
my $rootname = $input;
$rootname=~s/(.*)\///;
$rootname=~s/\..*//;
my $accepted = "$outdir/$rootname.accept.tmp";
my $rejected = "$outdir/$rootname.reject.sscan";
my $accs = "$outdir/$rootname.accno";


open (CHIMERS, $input) or die $!;
open (ACC, ">$accs") or die $!;
open (ACCEPTED, ">$accepted") or die $!;
open (REJECTED, ">$rejected") or die $!;


print "Search around alignment ends set to $threshold bases\n";
print "Reading genomic coordinates from $chrexons\n";
print "Reading alignments from $input\n";
print "Writing accepted sequences to $accepted\n";
print "Writing list of accepted accession numbers to $accs\n";
print "Writing newly rejected sequences to $rejected\n\n";

# Open file of exon starts/stops and store in a hash of arrays
open (EXONS, $chrexons) or die $!;
my $flag=0; # flag to check file format
while (<EXONS>) {	
	next if (/^#/);
	chomp;
	my @linearray = split (/\s+/);
	unless ($flag) {
		unless ($#linearray == 7) {
			die "***Quitting!*** Incorrect exon coordinate file format supplied to splicescan!\n\n";
		}
	}
	$flag=1;
	my $chr = $linearray[0];
	my $orient = $linearray[1];
	my $start = $linearray[2];
	my $end = $linearray[3];
	my $band = $linearray[4];
	my $exon = $linearray[5];
	my $gene = $linearray[6];
	
	if ($orient eq "-") {
		$exons{$chr}[0][0]{$start} = "$band $gene";
		$exons{$chr}[0][1]{$end} = "$band $gene";
	}
	
	elsif ($orient eq "+") {
		$exons{$chr}[1][0]{$start} = "$band $gene";
		$exons{$chr}[1][1]{$end} = "$band $gene";
	}
}

close EXONS;
undef $flag;

print ACCEPTED "#query_id chromosome %identity query_start query_end chrom_start chrom_end band gene\n";
print REJECTED "#query_id chromosome %identity query_start query_end chrom_start chrom_end\n";


while (<CHIMERS>) {

	$record{nlines}++;

	# extract chromsomal regions and process - only does this when a blank line is reached i.e. between results
	if (/^$/) {

		# Set up keys for stat hash access
		my ($a, $b) = sort keys %stats;
		
		my $outputa;
		my $outputb;
		my $afirst;
	
		# Extract alignment positions
		my ($mina, $minb, $maxa, $maxb, $acca, $accb) = 
		($stats{$a}[0], $stats{$b}[0], $stats{$a}[1], $stats{$b}[1], $stats{$a}[2], $stats{$b}[2]);
		
		my ($aorient, $borient, $cmina, $cmaxa, $cminb, $cmaxb) = 
		($stats{$a}[3], $stats{$b}[3], $stats{$a}[4], $stats{$a}[5], $stats{$b}[4], $stats{$b}[5]);

		my ($pida, $pidb) = ($stats{$a}[6], $stats{$b}[6]);

		# Create variables for chromosomal inpection regions
		my ($lowera, $uppera, $lowerb, $upperb, $chromalow, $chromahigh, $chromblow, $chrombhigh);
		my ($lower, $upper, $pholder);

                $chromalow = $cmina;

		$chromahigh = $cmaxa;
		
		$chromblow = $cminb;

		$chrombhigh = $cmaxb;
		

		# if sequence a's 5` end overlaps with b or is 3` of it
		if ($mina > $minb) {

			$afirst=0;

			# Calculate chromosomal inspection regions based on orientation
			$lowera = $chromalow - $threshold;
			$uppera = $chromalow + $threshold;
			$lowerb = $chrombhigh - $threshold;
			$upperb = $chrombhigh + $threshold;


			# Index values for searching exon position array based on orientation
			my ($aindexa, $aindexb, $bindexa, $bindexb);
			
			# Search appropriate array dependent on orientation
			if ($aorient==0) {
				$aindexa = 0;
				$aindexb = 1;}
			else {
				$aindexa = 1;
				$aindexb = 0;}

			if ($borient==0) {
				$bindexa = 0;
				$bindexb =0;
			}
			
			else {
				$bindexa = 1;
				$bindexb = 1;
			}

			# Create range array for areas being considered for exon boundaries
			my @arange = ($lowera..$uppera);
			my @brange = ($lowerb..$upperb);

			# Accept if boundary found in our allowable range for both sequences
			my $aaccept;
			my $baccept;
			my $element;

			foreach $element (@arange) {
				if (exists $exons{$a}[$aindexa][$aindexb]{$element}){
					$outputa = $exons{$a}[$aindexa][$aindexb]{$element};
					$aaccept=1;
					last;
				}
			}

			foreach $element (@brange) {
				if (exists $exons{$b}[$bindexa][$bindexb]{$element}){
					$outputb =$exons{$b}[$bindexa][$bindexb]{$element};
					$baccept=1;
					last;
				}
			}

			if (($aaccept)&&($baccept)) {$accept = 1;}

		}

		# if sequence b's 5` end shows overlap with a or is 3` of it
		elsif ($minb > $mina) {
			
			$afirst=1;
		
			# Calculate chromosomal inspection regions based on orientation
			$lowerb = $chromblow - $threshold;
			$upperb = $chromblow + $threshold;
			$lowera = $chromahigh - $threshold;
			$uppera = $chromahigh + $threshold;


			# Index values for searching exon position array based on orientation
			my ($aindexa, $aindexb, $bindexa, $bindexb);

			# Search appropriate array dependent on orientation
			if ($borient==0) {
				$bindexa = 0;
				$bindexb = 1;}
			else {
				$bindexa = 1;
				$bindexb = 0;}

			if ($aorient==0) {
				$aindexa = 0;
				$aindexb =0;
			}

			else {
				$aindexa = 1;
				$aindexb = 1;
			}

			# Create range array for areas being considered for exon boundaries
			my @arange = ($lowera..$uppera);
			my @brange = ($lowerb..$upperb);

			# Accept if boundary found in our allowable range for both sequences
			my $aaccept;
			my $baccept;
			my $element;

			foreach $element (@arange) {
					if (exists $exons{$a}[$aindexa][$aindexb]{$element}){
					$outputa = $exons{$a}[$aindexa][$aindexb]{$element};
					$aaccept=1;
					last;
				}
			}

			foreach $element (@brange) {
				if (exists $exons{$b}[$bindexa][$bindexb]{$element}){
					$outputb =$exons{$b}[$bindexa][$bindexb]{$element};
					$baccept=1;
					last;
				}
			}

			if (($aaccept)&&($baccept)) {$accept = 1;}


		}

		# Incase sequence was originally in wrong orientation	
		# if sequence a's 5` end overlaps with b or is 3` of it
		if (!$accept && $mina > $minb) {

			$afirst=1;

			# Calculate chromosomal inspection regions based on orientation
			$lowera = $chromalow - $threshold;
			$uppera = $chromalow + $threshold;
			$lowerb = $chrombhigh - $threshold;
			$upperb = $chrombhigh + $threshold;


			# Index values for searching exon position array based on orientation
			my ($aindexa, $aindexb, $bindexa, $bindexb);
			
			# Search appropriate array dependent on orientation
			if ($aorient==0) {
				$aindexa = 1;
				$aindexb = 1;}
			else {
				$aindexa = 0;
				$aindexb = 0;}

			if ($borient==0) {
				$bindexa = 1;
				$bindexb =0;
			}
			
			else {
				$bindexa = 0;
				$bindexb = 1;
			}

			# Create range array for areas being considered for exon boundaries
			my @arange = ($lowera..$uppera);
			my @brange = ($lowerb..$upperb);

			# Accept if boundary found in our allowable range for both sequences
			my $aaccept;
			my $baccept;
			my $element;

			foreach $element (@arange) {
				if (exists $exons{$a}[$aindexa][$aindexb]{$element}){
					$outputa = $exons{$a}[$aindexa][$aindexb]{$element};
					$aaccept=1;
					last;
				}
			}

			foreach $element (@brange) {
				if (exists $exons{$b}[$bindexa][$bindexb]{$element}){
					$outputb =$exons{$b}[$bindexa][$bindexb]{$element};
					$baccept=1;
					last;
				}
			}

			if (($aaccept)&&($baccept)) {$accept = 2;}
		}
		
		# Incase sequence was originally in wrong orientation
		# if sequence b's 5` end shows overlap with a or is 3` of it
		elsif (!$accept && $minb > $mina) {
			
			$afirst=0;
		
			# Calculate chromosomal inspection regions based on orientation
			$lowerb = $chromblow - $threshold;
			$upperb = $chromblow + $threshold;
			$lowera = $chromahigh - $threshold;
			$uppera = $chromahigh + $threshold;


			# Index values for searching exon position array based on orientation
			my ($aindexa, $aindexb, $bindexa, $bindexb);

			# Search appropriate array dependent on orientation
			if ($borient==0) { 
                                $bindexa = 1;
                                $bindexb = 1;}
                        else {
                                $bindexa = 0;
                                $bindexb = 0;}

                        if ($aorient==0) {
                                $aindexa = 1;
                                $aindexb =0;
                        }

                        else {
                                $aindexa = 0;
                                $aindexb = 1;
                        }

			# Create range array for areas being considered for exon boundaries
			my @arange = ($lowera..$uppera);
			my @brange = ($lowerb..$upperb);

			# Accept if boundary found in our allowable range for both sequences
			my $aaccept;
			my $baccept;
			my $element;

			foreach $element (@arange) {
					if (exists $exons{$a}[$aindexa][$aindexb]{$element}){
					$outputa = $exons{$a}[$aindexa][$aindexb]{$element};
					$aaccept=1;
					last;
				}
			}

			foreach $element (@brange) {
				if (exists $exons{$b}[$bindexa][$bindexb]{$element}){
					$outputb =$exons{$b}[$bindexa][$bindexb]{$element};
					$baccept=1;
					last;
				}
			}

			if (($aaccept)&&($baccept)) {$accept = 2;}

		}
		

		if ($accept) {
			
			# Annotate if original sequence was in reverse orientation
			if ($accept==2) {
				$acca="(RO)-$acca";
				$accb="(RO)-$accb";
			}

			# Prepare and print output data (5 prime first)
			$outputa = "$acca $a $pida $mina $maxa $cmina $cmaxa ".$outputa."\n";
			$outputb = "$accb $b $pidb $minb $maxb $cminb $cmaxb ".$outputb."\n";


			if ($afirst) {print ACCEPTED "$outputa"."$outputb\n";}
			else  {print ACCEPTED "$outputb"."$outputa\n";}
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
		$accept=0;
		$outputa="";
		$outputb="";
		$afirst="";
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
				die "***Quitting!*** Incorrect input file format supplied to splicescan!\n\n";
			}
		}
		# Assign their own variables for ease of use
		my ($acc, $chr, $pid, $min, $max, $cmin, $cmax)=($linearray[0], $linearray[1], $linearray[2], $linearray[3], $linearray[4], $linearray[5], $linearray[6]);
		# Clean accession and increment unique accession counter
		$acc=~s/\*//;
		$uniqacc{$acc}++;
		# Store max and min hit positions for each of the two chromosomes
		$stats{$chr}[0] = $min unless (defined $stats{$chr}[0] && $stats{$chr}[0] < $min);
		$stats{$chr}[1] = $max unless (defined $stats{$chr}[1] && $stats{$chr}[1] > $max);
		
		# Store accession number 
		$stats{$chr}[2] = $acc;

		# Store chromosome positions corresponding to the max and min hit positions and determine orientation
		if ($cmin > $cmax) {
			$stats{$chr}[3] = 0;
			$stats{$chr}[4] = $cmin unless (defined $stats{$chr}[4] && $stats{$chr}[4] > $cmin);
			$stats{$chr}[5] = $cmax unless (defined $stats{$chr}[5] && $stats{$chr}[5] < $cmax);
		}
		
		else {
			$stats{$chr}[3] = 1;
			$stats{$chr}[4] = $cmin unless (defined $stats{$chr}[4] && $stats{$chr}[4] < $cmin);
			$stats{$chr}[5] = $cmax unless (defined $stats{$chr}[5] && $stats{$chr}[5] > $cmax);
		}
		# Store %id
                $stats{$chr}[6] = $pid;
		
	}

}

close CHIMERS;
close ACCEPTED;
close REJECTED;

# Overwrite old accept file with new version
print "Replacing $input file with new version ($accepted). (Original will be saved as .old.sscan)\n";
system "cp $input $outdir/$rootname.accept.old.sscan";
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
print "\n\tsplicescan statistics\n\n";
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
    print "\nUsage:\n\n";
    print "\t $0 -i input filename (.accept file created by chimfilt.pl) -c file containing exon boundary locations\n";
    print "\t e.g $0 -i example.accept -c exon.boundaries\n\n";
    print "Options:\n";
    print "\t -t <threshold bp to scan around alignment breaks for known exon boundaries (default 3)>\n";
    print "\t -c <file containing known chromosomal exon boundaries\n";
    print "\t -d <output-dir>\n";
    print "\t -h <print this screen>\n";
   # print "\t -s screen mode \n";
   # print "\t -v verbose\n";
    print "NB: Input should have been filtered (blastfilter), sorted (UNIX sort), processed (chimfind) and filtered (chimfilt)\n";

}
