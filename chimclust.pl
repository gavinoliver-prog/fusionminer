#!/usr/bin/perl

# Title:        Chimclust (Part of FusionMiner)
#
# Function:     Program to cluster chimeric candidates output by splicescan
#               Considers individual accessions and their chromosomal alignment positions
#               Clusters together accessions whose breakpoints are within user-definable proximity
#
# Author:       Gavin Oliver
# 
# Date:         22 Apr 2009
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

my $input; # input file
my %opts=(); # hash for storage of command line arguments
my %chimers=(); # hash for storage of chimer alignment info
my $overlap; # Overlap/gap allowance for clustering by breakpoint location
my %stats=(); # hash for storage of general processing stats
my $outdir="."; # output directory
my $version = "22-apr-2009 1.0"; # version


# Set command line option flags
getopts('ho:i:d:', \%opts);

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

# If overlap defined on command line accept value otherwise use a default of 10
if (exists $opts{o}) {
        $overlap = $opts{o};
}

else {$overlap=10;}

# Set up output filename based on input filename
my $clustered = $input;
# Remove path
$clustered =~ s/(.*)\///;
$clustered = "$outdir/$clustered.clustered";

# Open chimers file and output file
open (CHIMERS, $input) or die $!;
open (CLUSTERED, ">$clustered") or die $!;

print CLUSTERED "#query_id chromosome %identity query_start query_end chrom_start chrom_end band exon\n";

print "\n\t$0: version $version. Use -h for options\n\n";

print "Permitted gap/overlap set to $overlap bases\n";
print "Reading alignments from $input\n";
print "Writing clustered results to $clustered\n"; 

while (<CHIMERS>) {
	$stats{lines}++; # Increment no. of lines processed
	next if ((/^#/) || (/^$/));

	# Split lines into fields and extract to variables
	my @linearray = split (/\s+/);
	my $line = $_;
	
	my $acc = $linearray[0];
	$acc=~s/\*//;
	my $chr = $linearray[1];
	my $qstart = $linearray[3];
	my $qend = $linearray[4];
	my $cstart = $linearray[5];
	my $cend = $linearray[6];

	$chimers{$acc}[0].=$line; # Store output string (input from chimer file)
	$chimers{$acc}[1]=0; # Already compared flag (prevent chimers being analysed more than once)
	###
	if (defined $chimers{$acc}[2]{$chr}) {$chr.=".";}

	# Store start position of alignment, unless smaller start position already stored from previous line	
	unless ((defined $chimers{$acc}[2]{$chr}[0]) && ($chimers{$acc}[2]{$chr}[0] < $qstart)) {
		$chimers{$acc}[2]{$chr}[0] = $qstart;
		$chimers{$acc}[2]{$chr}[2] = $cstart;
	}
	# Store end position of alignment, unless larger end position already stored from previous line
	unless ((defined  $chimers{$acc}[2]{$chr}[1]) &&  ($chimers{$acc}[2]{$chr}[1] > $qend)) {
		$chimers{$acc}[2]{$chr}[1] = $qend;
		$chimers{$acc}[2]{$chr}[3] = $cend;
	}
}

close CHIMERS;

# Foreach chimer
foreach my $acc (sort keys %chimers) {
	# Increment number of accessions considered
	$stats{accs}++;
	# Increment counter to denote seq has been considered once
	$chimers{$acc}[1]++ unless ($chimers{$acc}[1]); 
	# foreach chimer
	foreach my $acc2 (sort keys %chimers) {
		# unless comparing self with self or already considered
		unless ($acc eq $acc2 || $chimers{$acc2}[1]) { 
			# Extract and sort chromsomes for each chimer being compared
			my @array1 = sort keys %{$chimers{$acc}[2]};
			my @array2 = sort keys %{$chimers{$acc2}[2]};
			# If each hits same chromosomes
			if (defined $array1[0] && defined $array2[0] && defined $array1[1] && defined $array2[1] && $array1[0] eq $array2[0] && $array1[1] eq $array2[1]) {
				# Extract chromsomes and query end positions
				my $chr1 = $array1[0];
				my $chr2 = $array1[1];
				my $qend11 = $chimers{$acc}[2]{$chr1}[1];
				my $qend12 = $chimers{$acc}[2]{$chr2}[1];
				my $qend21 = $chimers{$acc2}[2]{$chr1}[1];
				my $qend22 = $chimers{$acc2}[2]{$chr2}[1];
		
				# If chromsome 'a' aligns 3` of chromsome 'b' for each sequence being considered
				if ($qend11 > $qend12 && $qend21 > $qend22) {
					
					# Extract the chromosomal start positions for 'a'
					my $cstart11 = $chimers{$acc}[2]{$chr1}[2];
					my $cstart21 = $chimers{$acc2}[2]{$chr1}[2];
					# Extract the chromosomal end positions for 'b'
					my $cend12 = $chimers{$acc}[2]{$chr2}[3];
					my $cend22 = $chimers{$acc2}[2]{$chr2}[3];
					
					# If the difference between the start/end postions for each alignment are <= overlap
					if ((abs($cstart11-$cstart21) <= $overlap) && (abs($cend12-$cend22) <= $overlap)) {
						# Increment the count of sequences clustered with for each alignment
						$chimers{$acc2}[1]++;
						$chimers{$acc}[1]++;
						# Append the output string of sequence 2 to sequence 1 and undef 2
						$chimers{$acc}[0].=$chimers{$acc2}[0];
						undef $chimers{$acc2}[0];
					}			
					
				}
				# If chromsome 'b' aligns 3` of chromsome 'a' for each sequence being considered
				elsif ($qend12 > $qend11 && $qend22 > $qend21) {
					
					# Extract the chromosomal start positions for 'b'
					my $cstart12 = $chimers{$acc}[2]{$chr2}[2];
					my $cstart22 = $chimers{$acc2}[2]{$chr2}[2];
					# Extract the chromosomal end positions for 'a'
					my $cend11 = $chimers{$acc}[2]{$chr1}[3];
                                        my $cend21 = $chimers{$acc2}[2]{$chr1}[3];

					# If the difference between the start/end postions for each alignment are <= overlap
					if ((abs($cstart12-$cstart22) <= $overlap) && (abs($cend11-$cend21) <= $overlap)) {
						# Increment the count of sequences clustered with for each alignment
						$chimers{$acc2}[1]++;
						$chimers{$acc}[1]++;
						# Append the output string of sequence 2 to sequence 1 and undef 2
						$chimers{$acc}[0].=$chimers{$acc2}[0];
						undef $chimers{$acc2}[0];
					}

				}
			}
		}
	}
} 

# Print computer times
my ($user,$system,$cuser,$csystem) = times;
print "\n** Processing finished.  User time=$user s **\n";

# Print output stats
print "\n\tchimclust statistics\n\n";
print " No. input lines processed\t\t\t:$stats{lines}\n";
print " No. distinct accession numbers\t\t\t:$stats{accs}\n\n";

# Get all cluster sizes in an array and remove the lowest and highest to use as an index
my %sizes=();
foreach my $acc (keys %chimers) {
	my $size = $chimers{$acc}[1];
	$sizes{$size}++;
}

# Print cluster sizes to screen and results to file using the index values perviously calculated
foreach my $index (reverse sort keys(%sizes)) {
	foreach my $acc (keys %chimers) {
		if ($chimers{$acc}[0] && $chimers{$acc}[1] == $index) {
			$stats{naccs}+=$index;
			if ($index >1) {
				print CLUSTERED "#$index\_accession_cluster\n$chimers{$acc}[0]\n";
			}
			else {	
				print CLUSTERED "#singlet_chimer\n$chimers{$acc}[0]\n";
			}
			$stats{$index}++;
		}
	}
	if ($index >1) {
		print " No. of chimer clusters with $index accessions:\t$stats{$index}\n";
	}
	
	elsif ($index ==1) {
		print " No. of singlet chimers:\t\t\t$stats{$index}\n" unless (!$stats{$index});
	}
}

# Warn if discrepancy in stats
if ($stats{naccs}!=$stats{accs}) {
   print "**WARNING: Naccepted + Nrejected <> N accession numbers\n\n";
}

# Subroutine to print help to screen
sub usage {

	print "\n\t$0 version $version\n";
	print "\tby G.R. Oliver\n";
	print "\nUsage:\n";
	print "\t $0 -i input filename (.accept file created by splicescan)\n";
	print "\t e.g $0 -i example.accept\n\n";
	print "Other options\n";
	print "\t -o <overlap/gap allowed at breakpoint (default 10)>\n";
	print "\t -d <output-dir>\n";
	print "\t -h <print this screen>\n";
	print "NB: Input should have been filtered (blastfilter), sorted (UNIX sort), processed (chimfind), refiltered (chimfilt) and reprocessed (chimsplice)\n";

}

	
