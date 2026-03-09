#!/usr/bin/perl

# Title:        ChimfiltSING (Part of FusionMiner)
#
# Function:     Perl program to take single chromosome files produced by chimfind and remove previously accepted chimeric candidates
#               where candidates have insufficient % ID, minimum alignment length of fall completely inside a single gene
#               Alignments can be rejected based on length or %id by specifying a value with the -l and -p flags
#               Gene coordinate file is specified with the -c flag
#
# Author:       Gavin Oliver
# 
# Date:         04 Jun 2009
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


my $version = "04 Jun 2009 1.0"; # version number
my %opts=(); # hash for storage of command line options
my $input; # input file
my $coords; # gene coordinate file
my $outdir="."; # output directory
my $id; # %id filter threshold
my $length; # alignment length filter threshold
my $flag=0;
my %generange=(); # hash for gene coordinates
my $threshold; # error margin allowed for gene boundaries

# Set command line option flags
getopts('hi:t:d:p:l:c:', \%opts);

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

if (exists $opts{c}) {
        $coords = $opts{c};
}

else {&usage;
        print "\n\nGene coordinates filename must be specified on command line!\n\n";
        exit;}

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

# If acceptable boundary error defined on command line accept value otherwise set to default of 3
if (exists $opts{t}) {
	$threshold = $opts{t};
}

else {$threshold=3;}

# Set up ouput filenames based on input filename
my $rootname = $input;
$rootname=~s/(.*)\///;
$rootname=~s/\..*//;
my $accepted = "$outdir/$rootname.accept.tmp";
my $rejected = "$outdir/$rootname.reject.cfilt";
my $accs = "$outdir/$rootname.accno";
# unclustered file
my $rootname2=$rootname;
$rootname2=~s/-.*//;
my $unclrej = "$outdir/$rootname2-sing.reject.cfilt";
my $unclfile = "$outdir/$rootname2-sing.accept";
my $uncltemp = "$unclfile.tmp";



my @stats=(); # hash for storage of genome alignment positions
my %record=(); # hash to keep a record of process stats
my $accept=0; # accept or reject result
my $linebuffer; # string for storage of result lines for later printing
my %acchash=(); # hash for storage of accepted accessions
my %uniqacc=(); # hash for storage of unique accessions processed
my %rejhash=(); # hash for storage of unique accessions rejected

print "\n\t$0: version $version. Use -h for options\n\n";

open (INFILE, $input) or die $!;
open (COORDS, $coords) or die $!;
open (ACC, ">$accs") or die $!;
open (ACCEPTED, ">$accepted") or die $!;
open (UNCLTEMP, ">$uncltemp") or die $!;
open (REJECTED, ">$rejected") or die $!;
open (UNCLREJ, ">$unclrej") or die $!;

print "% ID filter set to $id bases\n";
print "Length filter set to $length bases\n";
print "Search around alignment ends set to $threshold bases\n";
print "Reading genomic coordinates from $coords\n";
print "Reading alignments from $input\n";
print "Writing accepted sequences (clustered) to $accepted\n";
print "Writing accepted sequences (unclustered) to $uncltemp\n";
print "Writing list of accepted accession numbers to $accs\n";
print "Writing newly rejected sequences (clustered) to $rejected\n";
print "Writing newly rejected sequences (unclustered) to $unclrej\n\n";

print ACCEPTED "#query_id chromosome %identity query_start query_end chrom_start chrom_end\n";
print UNCLTEMP "#query_id chromosome %identity query_start query_end chrom_start chrom_end\n";
print REJECTED "#query_id chromosome %identity query_start query_end chrom_start chrom_end\n";
print UNCLREJ "#query_id chromosome %identity query_start query_end chrom_start chrom_end\n";

# Read in genomic coordinates of genes and store
while (<COORDS>) {
	next if (/^#/); # ignore comment lines
        chomp;
        my @linearray = split (/\s+/); # split results into an array
        # Check input file format
	unless ($flag) {
		unless ($#linearray == 7) {
			die "***Quitting!*** Incorrect gene coordinate file format supplied!\n\n";
		}
	}
        $flag=1;

	my $chr = $linearray[0];
	#my $orient = $linearray[1];
	my $start = $linearray[2];
	my $end = $linearray[3];
	#my $band = $linearray[4];
	#my $exon = $linearray[5];
	my $gene = $linearray[6];
	my $ensgene = $linearray[7];

	$generange{$chr}{$ensgene}[0]=$start unless (defined $generange{$chr}{$ensgene}[0] && $generange{$chr}{$ensgene}[0]<$start);
	$generange{$chr}{$ensgene}[1]=$end unless (defined $generange{$chr}{$ensgene}[1] && $generange{$chr}{$ensgene}[1]>$end);
}

close COORDS;

	
while (<INFILE>) {
	
	$record{nlines}++;

	# calculate overlap/gap or lack of - only does this when a blank line is reached i.e. between results
	if (/^$/) {
	
		$accept =1;
		

		# If specified on command line, first check that %ids are sufficient otherwise reject
		if ($id) {
			my $pid1 = $stats[2];
			if ($pid1 < $id) {
				print REJECTED "$linebuffer\n";
				# Write rejected accessions to a hash for stats
				my @linearray = split(/\s+/, $linebuffer);
				my $acc = $linearray[0];
				# Remove asterisks so we only have unique entries
				$acc=~s/\*//;
				$rejhash{$acc}++;
				# Reset variables for processing next result
				$linebuffer="";
				@stats=();
				$accept = 0;
				###########
				#print REJECTED "ID rejection\n";
				next;
			}
		}

		# Extract alignment positions 		
		my ($mina, $maxa) = ($stats[0], $stats[1]);

		# If specified on command line, first check that alignment lengths are sufficient otherwise reject
		if ($length) {
			
			# Extract maximum alignment lengths
			my $minrangea = $stats[3];

			if ($minrangea < $length){
				print REJECTED "$linebuffer\n";
				# Write rejected accessions to a hash for stats
				my @linearray = split(/\s+/, $linebuffer);
				my $acc = $linearray[0];
				# Remove asterisks so we only have unique entries
				$acc=~s/\*//;
				$rejhash{$acc}++;
				# Reset variables for processing next result
				$linebuffer="";
				@stats=();
				$accept = 0;
				##########
				#print REJECTED "Length rejection\n";
				next;
			}
		}

		# Check if gene boundaries are crossed
		my ($chr, $cmin, $cmax) = ($stats[10], $stats[4], $stats[5]);
		
		foreach my $key (keys %{$generange{$chr}}) {
	
			my $genemin = $generange{$chr}{$key}[0];
			my $genemax = $generange{$chr}{$key}[1];
			$genemin-=$threshold;
			$genemax+=$threshold;

			if (($cmin >= $genemin && $cmin <= $genemax) && ($cmax >= $genemin && $cmax <= $genemax)) {
				$accept = 0;
				##########
				#print "$key\t$chr\t$genemin\t$genemax\t$cmin\t$cmax\n";
				#print REJECTED "Within a single gene\n";
			}

			# Exit foreach loop if we find a gene we're inside
			last if ($accept==0);		
		}

		# Check both start and end of sequence exist within a gene (not necessarily in exons at this stage)
		###
		if ($accept) {
			my $i=0;
			foreach my $key (keys %{$generange{$chr}}) {

				my $genemin = $generange{$chr}{$key}[0];
				my $genemax = $generange{$chr}{$key}[1];
				$genemin-=$threshold;
				$genemax+=$threshold;

				if (($cmin >= $genemin && $cmin <= $genemax) || ($cmax >= $genemin && $cmax <= $genemax)) {
					$i++;
				#	print "$chr $key\n";
				}
				if ($i==2) {last;}
			}
			#########
			if ($i!=2) {$accept=0;}#print REJECTED "Start and end\n";}
		}

			
		###
	 	
		if ($accept) {
			print ACCEPTED "$linebuffer\n";
	
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
		@stats=();
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
		unless (defined $stats[0] && $stats[0] < $min) {$stats[0] = $min; $stats[4] = $cmin;}
		unless (defined $stats[1] && $stats[1] > $max) {$stats[1] = $max; $stats[5] = $cmax;}
		# Store minimum %id for each chromosomal alignment
		$stats[2] = $pid unless (defined $stats[2] && $stats[2] < $pid);
		# Store MINimum alignment length for each chromosomal alignment
                my $range = ($max - $min)+1;
                $stats[3] = $range unless (defined $stats[3] && $stats[3] < $range);
		# Store total number of bases aligning
		$stats[6] += $range;
		# Keep total of %id x bases to caculate overall %id where there are multiple alignments
		my $avid = $pid*$range;
		$stats[7] += $avid;
		# Keep record of number of alignmnets for a given chromosome
		$stats[8]++;
		# Store accession 
		$stats[9]=$acc;
		# Store chromosome
		$stats[10]=$chr;
	}
}

close INFILE;
close ACCEPTED;
close REJECTED;

# Overwrite old accept file with new version
print "Replacing $input file with new version ($accepted). (Original will be saved as .old.cfilt)\n";
system "cp $input $outdir/$rootname.accept.old.cfilt";
system "mv $accepted $outdir/$rootname.accept";

# Create new unclustered file
my $uncflag=0;
open (UNCLFILE, "$unclfile") or die $!;
while (<UNCLFILE>) {
	next if (/^#/);
	if (/^$/ && $uncflag) {print UNCLTEMP $_; next;}
	elsif (/^$/ && !$uncflag) {print UNCLREJ $_; next;}
	my @linearray = split (/\s+/);
	my $uncid = $linearray[0];
	if (exists $acchash{$uncid}) {print UNCLTEMP $_; $uncflag=1;}
	else {print UNCLREJ $_; $uncflag=0;}
}

close UNCLTEMP;
close UNCLFILE;
# Overwrite old unclustered file with new version
print "Replacing $unclfile file with new version ($uncltemp).  (Original will be saved as .old.cfilt)\n";
system "mv $unclfile $outdir/$unclfile.old.cfilt";
system "mv $uncltemp $outdir/$unclfile";

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
    print "\t $0 -i input filename (-clsing.accept file created by chimfind)\n";
    print "\t e.g $0 -i example.accept\n\n";
    print "Other options\n";
    print "\t -c <file containing gene coordinates>\n";
    print "\t -t <error margin in bp to scan around known gene boundaries (default 3)>\n";
    print "\t -p <% id filter threshold (default 0)>\n";
    print "\t -l <alignment length filter threshold (default 0)>\n";
    print "\t -d <output-dir>\n";
    print "\t -h <print this screen>\n";
    print "NB: Input should have been filtered (blastfilter), sorted (UNIX sort) and processed (chimfind)\n";

}

