#!/usr/bin/perl

# Title:        SplicescanSING (Part of FusionMiner)
#
# Function:     Program to filter single chromosome chimeric candidates output by chimfilt.sing
#               Utilises genomic exon coordinates to confirm or reject
#               chimeric candidates from chimfilt.sing
#               Rejects sequences hitting 1 or more than 2 genes
#               Rejects chimers with 2 similar gene names as these are usually related
#               Also filters on % id and length
#
# Author:       Gavin Oliver
# 
# Date:         18 Aug 2009
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

my $version = "18-aug-2009 1.0"; #version number
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
my @lines; # array for storage of input lines 
my @linesclean; # array for cleaned input 
my $threshold; # threshold distance (in bp) allowed upstream/downstream of alignment start/finish to search for splice signal 
my %exons=(); # hash to store chromosomal exon start/stop positions
my %generange=(); # hash to store whole gene ranges
my $overlap; # overlap/gap threshold for chimer join point 
my $id; # % id filter
my $length; # length filter


# Set command line option flags
getopts('ht:i:d:c:o:l:p:', \%opts);

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

# If overlap defined on command line accept value otherwise use a default of 10
if (exists $opts{o}) {
        $overlap = $opts{o};
}

else {$overlap=10;}

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
my $rejected = "$outdir/$rootname.reject.sscan";
my $accs = "$outdir/$rootname.accno";


open (CHIMERS, $input) or die $!;
open (ACC, ">$accs") or die $!;
open (ACCEPTED, ">$accepted") or die $!;
open (REJECTED, ">$rejected") or die $!;

print "% ID filter set to $id bases\n";
print "Length filter set to $length bases\n";
print "Overlap / Gap filter set to $overlap bases\n";
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
	# Check input file format for coordinates file
	unless ($flag) {
		unless ($#linearray == 7) {
			die "***Quitting!*** Incorrect exon coordinate file format supplied to splicescan!\n\n";
		}
	}
	$flag=1;
	# Break lines into variables
	my $chr = $linearray[0];
	my $orient = $linearray[1];
	my $start = $linearray[2];
	my $end = $linearray[3];
	my $band = $linearray[4];
	my $exon = $linearray[5];
	my $gene = $linearray[6];
	my $ensgene = $linearray[7];
	
	# Store upper and lower limits for gene coordinates
	$generange{$chr}{$orient}{$ensgene}[0]=$start 
	unless (defined $generange{$chr}{$orient}{$ensgene}[0] && $generange{$chr}{$orient}{$ensgene}[0]<$start);

	$generange{$chr}{$orient}{$ensgene}[1]=$end 
	unless (defined $generange{$chr}{$orient}{$ensgene}[1] && $generange{$chr}{$orient}{$ensgene}[1]>$end);
	
	# Store exon coordinates and associated annotation info
	if ($orient eq "-") {

		
		$exons{$chr}[0][0]{$start}[0]= "$band";
		$exons{$chr}[0][0]{$start}[1]= "$gene";
		$exons{$chr}[0][0]{$start}[2]= "exon_$exon";

		$exons{$chr}[0][1]{$end}[0] = "$band";
		$exons{$chr}[0][1]{$end}[1] = "$gene";
		$exons{$chr}[0][1]{$end}[2] = "exon_$exon";
		
	}
	
	elsif ($orient eq "+") {

		$exons{$chr}[1][0]{$start}[0] = "$band";
		$exons{$chr}[1][0]{$start}[1] = "$gene";
		$exons{$chr}[1][0]{$start}[2] = "exon_$exon";

		$exons{$chr}[1][1]{$end}[0] = "$band";
		$exons{$chr}[1][1]{$end}[1] = "$gene";
		$exons{$chr}[1][1]{$end}[2] = "exon_$exon";

	}
}

close EXONS;
undef $flag;

print ACCEPTED "#query_id chromosome %identity query_start query_end chrom_start chrom_end band gene\n";
print REJECTED "#query_id chromosome %identity query_start query_end chrom_start chrom_end\n";

# register for how many genes matched
my %genereg=();
my $acc;

while (<CHIMERS>) {

	$record{nlines}++;

	# extract chromsomal regions and process - only does this when a blank line is reached i.e. between results
	if (/^$/) {

		# Variables for final results output
		my $outputa;
		my $outputb;
	
		# Process each line of a result
		foreach my $line (@lines) {
			my @linearray = split (/\s+/, $line);
			#my $acc = $linearray[0];
			my $chr = $linearray[1];
			my $id = $linearray[2];
			my $min = $linearray[3];
			my $max = $linearray[4];
			my $cmin = $linearray[5];
			my $cmax = $linearray[6];
		
			my $orient;
	
			if ($cmin > $cmax) {$orient="-";}
			else {$orient = "+";}
			
			
			# Check how many genes are hit by the sequence being considered
			foreach my $key (keys %{$generange{$chr}{$orient}}) {
				my $genemin = $generange{$chr}{$orient}{$key}[0];
				my $genemax = $generange{$chr}{$orient}{$key}[1];
		
				# Allow for some error
				$genemin-=$threshold;
				$genemax+=$threshold;

				# Only consider alignments over 50 bases?
				my $testrange = $max-$min;
				if ($testrange > 50) {
				if (($cmin >= $genemin && $cmin <= $genemax) && ($cmax >= $genemin && $cmax <= $genemax)) {
					# Increment number of genes every time an alignment falls within one
					$genereg{$key}++;
					$line.=" $key";
					last;
				}
				}
			}
		
			
		}

		# Test again incase sequence was in wrong orientation
		my $geneno = scalar (keys %genereg);
		if ($geneno ne 2) {
			%genereg=();
			foreach my $line (@lines) {
			my @linearray = split (/\s+/, $line);
                       	#my $acc = $linearray[0];
                       	my $chr = $linearray[1];
                        my $id = $linearray[2];
                        my $min = $linearray[3];
                        my $max = $linearray[4];
                        my $cmin = $linearray[5];
                        my $cmax = $linearray[6];

                        my $orient;

                        if ($cmin > $cmax) {$orient="+";}
                        else {$orient = "-";}


                        # Check how many genes are hit by the sequence being considered
                        foreach my $key (keys %{$generange{$chr}{$orient}}) {
                                my $genemin = $generange{$chr}{$orient}{$key}[0];
                                my $genemax = $generange{$chr}{$orient}{$key}[1];

                                # Allow for some error
                                $genemin-=$threshold;
                                $genemax+=$threshold;

                                # Only consider alignments over 50 bases?
                                my $testrange = $max-$min;
                                if ($testrange > 50) {
                                	if (($cmin >= $genemin && $cmin <= $genemax) && ($cmax >= $genemin && $cmax <= $genemax)) {
                                        # Increment number of genes every time an alignment falls within one
                                        $genereg{$key}++;
                                        $line.=" $key";
                                        last;
                                	}
                                }
                       	 }
			}
		}


			
			
		$geneno = scalar (keys %genereg);
		my $accept1;
		my $accept2;
		# Only accept if two genes are hit
		if ($geneno == 2) {$accept1=1;}
		else {$accept1=0;}
		
		# Hash and variables for next round of checks
		my %testhash=();
		my $a;
                my $b;
		my $rev;

		# If accepted as hitting two genes
		if ($accept1) {
			foreach my $line (@lines) {
				# Remove lines that didn't previously hit a gene
				my @linearray = split (/\s+/, $line);
				next unless ($linearray[7]);
				push @linesclean, $line;
			}

			# For lines hitting genes
			foreach my $line (@linesclean) {
                                my @linearray = split (/\s+/, $line);
                                my $acc = $linearray[0];
                                my $chr = $linearray[1];
                                my $id = $linearray[2];
                                my $min = $linearray[3];
                                my $max = $linearray[4];
                                my $cmin = $linearray[5];
                                my $cmax = $linearray[6];
				my $ensgene = $linearray[7];
				my $orient;	
				if ($cmin > $cmax) {$orient="-";}
                                else {$orient = "+";}

				# Cluster results to form one for each gene
				unless (defined $testhash{$ensgene} && $testhash{$ensgene}[1] > $max) {
				
					$testhash{$ensgene}[0] = $min unless 
					(defined $testhash{$ensgene}[0] && $testhash{$ensgene}[0] < $min);
					$testhash{$ensgene}[1] = $max;
					
					$testhash{$ensgene}[2] = $cmin unless 
                                        (defined $testhash{$ensgene}[0] && $testhash{$ensgene}[0] < $min);
					$testhash{$ensgene}[3] = $cmax;	
			
					$testhash{$ensgene}[4] = $orient;

					$testhash{$ensgene}[5] = $chr;
	
					$testhash{$ensgene}[6] = $acc;
					# Variables to produce average weighted %id
					my $range = ($max - $min)+1;
					my $avid = $id*$range;

					$testhash{$ensgene}[7] +=$avid;
					$testhash{$ensgene}[8] +=$range;
				}

                        }

			
			# Extract gene names and figure which one is 3` - make $a the 5` gene and $b the 3` gene
			my @genes = keys %testhash;
			my $gene1 = $genes[0];
			my $gene2 = $genes[1];
			
		
			if ($testhash{$gene1}[0] > $testhash{$gene2}[0]) {
				$a = $gene2;
				$b = $gene1;
			}
			else {
				$a = $gene1;
				$b = $gene2;
			}

			# Check length and %ids and reject if necessary
			if ($length) {
				my $alen = ($testhash{$a}[1] - $testhash{$a}[0])+1;
				my $blen = ($testhash{$b}[1] - $testhash{$b}[0])+1;
				if ($alen < $length || $blen < $length) {
					$accept1=0;
				}
			}
			# id check if length check passed only
			if ($id && $accept1) {
				my $aid = sprintf("%.2f",$testhash{$a}[7]/$testhash{$a}[8]);
				my $bid = sprintf("%.2f",$testhash{$b}[7]/$testhash{$b}[8]);
				if ($aid < $id || $bid < $id) {
					$accept1=0;
				}
			}
			
			
			# Check the level of overlap/gap at the fusion point - reject if too great
			if (abs($testhash{$b}[0]-$testhash{$a}[1])>$overlap ||!$accept1) {$accept1=0;}
			
			else {
				# Create variables for chromosomal inpection regions
                		my ($lowera, $uppera, $lowerb, $upperb, $chromalow, $chromahigh, $chromblow, $chrombhigh);
                		#my ($lower, $upper, $pholder);

				$chromalow = $testhash{$a}[2];
				$chromahigh = $testhash{$a}[3];

				$chromblow = $testhash{$b}[2];
				$chrombhigh = $testhash{$b}[3];

				my $aorient = $testhash{$a}[4];
				my $borient = $testhash{$b}[4];
			
				my $chra = $testhash{$a}[5];
				my $chrb = $testhash{$b}[5];
			
				my ($aindexa, $aindexb, $bindexa, $bindexb);

				# Preliminary acceptance if break agrees with exon boundaries
                                my $aaccept;
                                my $baccept;
                                my $element;
                                my $aband;
                                my $agene;
                                my $aexon;
                                my $bband;
                                my $bgene;
                                my $bexon;
				
				my @arange=();
				my @brange =();

				# Figure chromosomal inspection regions based on orientation and position
				if ($aorient eq "+") {
					$aindexa = 1;
					$aindexb = 1;
					$uppera = $chromahigh + $threshold;
					$lowera = $chromahigh - $threshold;
				}
				else {
					$aindexa = 0;
					$aindexb = 0;
					$uppera = $chromahigh + $threshold;
                                	$lowera = $chromahigh - $threshold;
				}

				if ($borient eq "+") {
					$bindexa = 1;
					$bindexb = 0;
					$upperb = $chromblow + $threshold;
                                	$lowerb = $chromblow - $threshold;
				}
				else {
					$bindexa = 0;
					$bindexb = 1;
					$upperb = $chromblow + $threshold;
                                	$lowerb = $chromblow - $threshold;
				}

				@arange = ($lowera..$uppera);
                        	@brange = ($lowerb..$upperb);
		
				foreach $element (@arange) {
                                        if (exists $exons{$chra}[$aindexa][$aindexb]{$element}){
                                        	$aband=$exons{$chra}[$aindexa][$aindexb]{$element}[0];
						$agene=$exons{$chra}[$aindexa][$aindexb]{$element}[1];
						$aexon=$exons{$chra}[$aindexa][$aindexb]{$element}[2];
						$outputa = "$aband $agene";
                                        	$aaccept=1;
                                        	last;
                                	}
                        	}
			
				foreach $element (@brange) {
                                	if (exists $exons{$chrb}[$bindexa][$bindexb]{$element}){
						$bband=$exons{$chrb}[$bindexa][$bindexb]{$element}[0];
						$bgene=$exons{$chrb}[$bindexa][$bindexb]{$element}[1];
						$bexon=$exons{$chrb}[$bindexa][$bindexb]{$element}[2];
                                        	$outputb="$bband $bgene";
                                        	$baccept=1;
                                        	last;
                                	}
                        	}
						
				if (($aaccept)&&($baccept)) {$accept2 = 1;}

				####
				####Test again incase sequence was in wrong orientation
				else {
		
						
					# Figure chromosomal inspection regions based on orientation and position
					if ($aorient eq "+") {
						$aindexa = 0;
						$aindexb = 1;
						$uppera = $chromahigh + $threshold;
						$lowera = $chromahigh - $threshold;
					}
					else {
						$aindexa = 1;
						$aindexb = 0;
						$uppera = $chromahigh + $threshold;
                                		$lowera = $chromahigh - $threshold;
					}

					if ($borient eq "+") {
						$bindexa = 0;
						$bindexb = 0;
						$upperb = $chromblow + $threshold;
                                		$lowerb = $chromblow - $threshold;
					}
					else {
						$bindexa = 1;
						$bindexb = 1;
						$upperb = $chromblow + $threshold;
                                		$lowerb = $chromblow - $threshold;
					}

					@arange = ($lowera..$uppera);
                        		@brange = ($lowerb..$upperb);
		
					foreach $element (@arange) {
                                        	if (exists $exons{$chra}[$aindexa][$aindexb]{$element}){
                                        		$aband=$exons{$chra}[$aindexa][$aindexb]{$element}[0];
							$agene=$exons{$chra}[$aindexa][$aindexb]{$element}[1];
							$aexon=$exons{$chra}[$aindexa][$aindexb]{$element}[2];
							$outputa = "$aband $agene";
                                        		$aaccept=1;
                                        		last;
                                		}
                        		}
			
					foreach $element (@brange) {
                                		if (exists $exons{$chrb}[$bindexa][$bindexb]{$element}){
							$bband=$exons{$chrb}[$bindexa][$bindexb]{$element}[0];
							$bgene=$exons{$chrb}[$bindexa][$bindexb]{$element}[1];
							$bexon=$exons{$chrb}[$bindexa][$bindexb]{$element}[2];
                                        		$outputb="$bband $bgene";
                                        		$baccept=1;
                                        		last;
                                		}
                        		}#######################
					if (($aaccept)&&($baccept)) {
						$accept2 = 1;
						$testhash{$a}[6]="(RO)-$testhash{$a}[6]";
                                		$testhash{$b}[6]="(RO)-$testhash{$b}[6]";
						$rev=1;
						#($a, $b) = ($b, $a);
					}
				}
				
				# Reject similar gene names
				if ($accept2) {	
					my $gtesta = substr $agene, 0, 3;
					my $gtestb = substr $bgene, 0, 3;
					if ($gtesta eq $gtestb) {
						$accept2=0;
					}
				}

					
			}

		}
		
		# Prepare and print accepted results	
		if ($accept1 && $accept2) {
			my $chra = $testhash{$a}[5];
			my $chrb = $testhash{$b}[5];
			my $acca = $testhash{$a}[6];
			my $accb = $testhash{$b}[6];
			my $mina = $testhash{$a}[0];
			my $minb = $testhash{$b}[0];
			my $maxa = $testhash{$a}[1];
			my $maxb = $testhash{$b}[1];
			my $cmina = $testhash{$a}[2];
                        my $cminb = $testhash{$b}[2];
			my $cmaxa = $testhash{$a}[3];
			my $cmaxb = $testhash{$b}[3];
			my $pida = sprintf("%.2f",$testhash{$a}[7]/$testhash{$a}[8]);
			my $pidb = sprintf("%.2f",$testhash{$b}[7]/$testhash{$b}[8]);
			
			$outputa = "$acca $chra $pida $mina $maxa $cmina $cmaxa ".$outputa."\n";
			$outputb = "$accb $chrb $pidb $minb $maxb $cminb $cmaxb ".$outputb."\n";	
			
			($outputa, $outputb) = ($outputb, $outputa) if ($rev);

			print ACCEPTED "$outputa";
			print ACCEPTED "$outputb\n";
			$acchash{$acc}++;
		}
		
		# Print rejected results
		else {
			print REJECTED join "\n", @lines;
			print REJECTED "\n\n";
			$rejhash{$acc}++;
		}
	
		# Reset variables
		@lines=();
		@linesclean=();
		undef $acc;
		%genereg=();
		

	}

	# Prior to reaching a blank line, extract result from file for processing
        else {
		next if (/^#/); # ignore comment lines
		# Place line into an array for processing
		chomp;
		push @lines, $_;
		my @linearray = split (/\s+/); # split results into an array
		# Check input file format
		unless (%uniqacc) {
			unless ($#linearray ==6) {
				die "***Quitting!*** Incorrect input file format supplied to splicescan!\n\n";
			}
		}
		# Extract accession
		$acc = $linearray[0];
		# Clean accession and increment unique accession counter
		$acc=~s/\*//;
		$uniqacc{$acc}++;
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
    print "\t -p <% id filter threshold (default 0)>\n";
    print "\t -l <alignment length filter threshold (default 0)>\n";
    print "\t -d <output-dir>\n";
    print "\t -h <print this screen>\n";
   # print "\t -s screen mode \n";
   # print "\t -v verbose\n";
    print "NB: Input should have been filtered (blastfilter), sorted (UNIX sort), processed (chimfind) and filtered (chimfilt)\n";

}
