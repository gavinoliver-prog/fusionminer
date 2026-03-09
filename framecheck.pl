#!/usr/bin/perl


# Title:        Framecheck (Part of FusionMiner)
# 
# Function:     Considers the alignment upstream and downstream of the fusion point and attempts to determine
#               if reading frame is preseved at the fusion point.       
#
# Author:       Gavin Oliver
#
# Date:         20 Oct 2010
# 
# Version:      V1.1 
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

my $infile; # input file
my %opts=(); # hash for storage of command line arguments
my $threshold; # discrepancy allowed in bp
my $outdir="."; # output directory
my $version = "20-oct-2010 1.1"; # version
my $coordfile; # genomic coordinate file
my $resultfile; # alignment results
my %lines=(); # hash for line number stats
my %uniqaccs=(); # hash for no. of unique accessions

# Set command line option flags
getopts('ht:i:d:c:r:', \%opts);

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
        $infile = $opts{i};
}

else {&usage;
        print "\n\nInput filename must be specified!\n\n";
        exit;
}

# Quit unless coord file is defined on command line
if (exists $opts{c}) {
        $coordfile = $opts{c};
}

else {&usage;
        print "\n\nGenomic coordinate filename must be specified!\n\n";
        exit;
}

# Quit unless result file is defined on command line
if (exists $opts{r}) {
        $resultfile = $opts{r};
}

else {&usage;
        print "\n\nInput results filename must be specified!\n\n";
        exit;
}

# If threshold defined on command line accept value otherwise use a default of 3
if (exists $opts{t}) {
        $threshold = $opts{t};
}

else {$threshold=3;}

# Set up output filename based on input filename
my $outfile = $infile;
# Remove path
$outfile =~ s/(.*)\///;
$outfile = "$outdir/$outfile.framechecked";

# Open coordinates file and process
open (COORDS, $coordfile) or die $!;


print "Permitted exon boundary discrepancy set to $threshold bases\n";
print "Reading genomic coordinates from $coordfile\n";
print "Reading chimers from $infile\n";
print "Writing framechecked chimers to $outfile\n";

# Hashes to store genomic coordinates and transcript ids
my %coords=();
my %ids=();

# Read in and store genomic coordinates
while (<COORDS>) {
	chomp;

	# Extract required values on each line
	my @linearray = split (/\s/, $_);
	my $tranid = $linearray[0];
	my $cdsst = $linearray[1];
	my $cdsend = $linearray[2];
	my $genename = $linearray[3];
	my $exonrank = $linearray[4];
	my $exonstart = $linearray[5];
	my $exonend = $linearray[6];
	my $utr5start = $linearray[7];
	my $utr5end = $linearray[8];
	my $utr3start = $linearray[9];
	my $utr3end = $linearray[10];
	
	# CDS start and end for individual exons
	if ($cdsst) {$coords{$tranid}[0]{$exonrank}[0] = $cdsst;}
	if ($cdsend) {$coords{$tranid}[0]{$exonrank}[1] = $cdsend;}

	# Store 5` and 3` utr starts and ends 
	if ($utr5start) {$coords{$tranid}[1] = $utr5start;}
	if ($utr5end) {$coords{$tranid}[2] = $utr5end;}
	if ($utr3start) {$coords{$tranid}[3] = $utr3start;}
	if ($utr3end) {$coords{$tranid}[4] = $utr3end;}
	
	# Start and end of individual exons
	unless ($coords{$tranid}[0]{$exonrank}[2] 
	&& $coords{$tranid}[0]{$exonrank}[2] < $exonstart) {$coords{$tranid}[0]{$exonrank}[2]=$exonstart;}
        unless ($coords{$tranid}[0]{$exonrank}[3] 
	&& $coords{$tranid}[0]{$exonrank}[3] > $exonstart) {$coords{$tranid}[0]{$exonrank}[3]=$exonend;}

	# Store all transcript ids for each gene
	$ids{$genename}{$tranid}++;
}

close COORDS;


# Open output file and print header
open (OUTFILE, ">$outfile") or die $!;
print OUTFILE "#query_id chromosome %identity query_start query_end chrom_start chrom_end band gene fusion_bases in_frame_fusion?\n";

# Variables for analysing rf
my %chimers=();
my $segment=0;
my $length=0;

# Get sequence ids for extraction of ful results from .filt.srted file and put in hash
open (INFILE,  $infile) or die $!;
while (<INFILE>) {

        next if (/^#/ || /^\s$/);
	my @linearray = split (/\s+/, $_);
	my $seqid = $linearray[0];
	if ($seqid =~/^\(RO\)\-/) {$seqid=~s/^\(RO\)\-//;}
	$chimers{$seqid}++;
}
close INFILE;

# Open result file to parse those pertaining to the chimers
open (RESULTS, $resultfile) or die $!;
while (<RESULTS>) {
	next if (/^#/);
	my @linearray = split (/\s+/, $_);
	my $seqid = $linearray[0];
	if ($chimers{$seqid}) {$chimers{$seqid}.=" $_";}
}

# Open infile to process chimers
open (INFILE, $infile) or die $!;

# Array for storage of output and coding tests
my @outtest=();
my $fusion=0; # no of fusion proteins
my $junk=0; # non coding fusions
my $unknown=0; #unable to determine

while (<INFILE>) {
	
	$lines{nlines}++;
	my %outhash=();
	next if (/^#query/);
	if (/^#/) {print OUTFILE $_;next;}
	if (/^\s$/) {
		print OUTFILE "$_";
		# Reset segment
		$segment = 0;
		##
		%outhash=();
	}

	else {	chomp;
		my $line = $_;
		# Increment to show which segment is being processed
		$segment++;
		$length=0;

		my @linearray = split (/\s+/, $_);
		my $seqid = $linearray[0];
		my $chr = $linearray[1];
		my $chrstart = $linearray[5];
		my $chrend = $linearray[6];
		my $gene = $linearray[8];
		my $upperexon=$1;
		
		if ($seqid =~/^\(RO\)\-/) {
			($chrstart, $chrend) = ($chrend, $chrstart);
			$seqid=~s/^\(RO\)\-//;
		}

		unless ($ids{$gene}) {print OUTFILE "$line Gene_not_defined\n";next;}

		$uniqaccs{$seqid}++;


		# For processing 5` segment	
		if ($segment%2) {
			my $upperindex;
			my $flag=0;
			my $agree;
			my $exception;
			my $rem;
		
	
			# Analyse each transcript corresponding to current gene
			foreach my $tranid (keys %{$ids{$gene}}) {
				#print "$tranid\n";
				$length=0;
				$agree=0;
				$exception=0;
				$rem = undef;

				# If alignment ends in 3` utr
				if (defined $coords{$tranid}[3] && $chrend>=$coords{$tranid}[3]
				&& $chrend<=$coords{$tranid}[4]) {
					$rem="Ends_in_3`_utr";
					$exception=1;
				}
				# If alignmnet ends immediately prior to 3`utr 
				elsif (defined $coords{$tranid}[3] && $chrend==($coords{$tranid}[3]-1)) {
					$rem= "Ends_at_utilised_stop_codon";
					$exception=1;
				}

				
				# Calculate joining exon number (upperindex) by determining which exon corresponds to the breakpoint
				foreach my $exonno (sort {$a<=>$b} keys %{$coords{$tranid}[0]}) {
					my $start = $coords{$tranid}[0]{$exonno}[2];
					my $end = $coords{$tranid}[0]{$exonno}[3];
					if ( defined $coords{$tranid}[0]{$exonno}[3] && abs($coords{$tranid}[0]{$exonno}[3]-$chrend)<=$threshold || 
					defined $coords{$tranid}[0]{$exonno}[2] && abs($coords{$tranid}[0]{$exonno}[2]-$chrend)<=$threshold) {
						$upperindex = $exonno;
						last;
					}
					
				}
				# Move to process next transcript if no exons in the current one correspond to the breakpoint
				next if (! defined $upperindex);

				# Figure out number of bases donated and number of exon boundaries that agree
				foreach my $exonno (1..$upperindex) {
					if (defined $coords{$tranid}[0]{$exonno}[0]) {	
						$flag=1;
						$length +=(($coords{$tranid}[0]{$exonno}[1] - $coords{$tranid}[0]{$exonno}[0])+1);
						my @resultarray = split (/\n/, $chimers{$seqid});
						foreach my $result (@resultarray){
							my @linearray = split (/\s+/, $result);
							my $chrom = $linearray[2];
							my $cstart=$linearray[9];
							my $cend=$linearray[10];
							if ($cstart > $cend) {
								($cstart, $cend) = ($cend, $cstart);
							}
							if ($chr eq $chrom) {
								if (defined $coords{$tranid}[0]{$exonno}[2]
								&& abs($coords{$tranid}[0]{$exonno}[2]-$cstart)<=$threshold){
									$agree+=0.5;
								}
								if (defined $coords{$tranid}[0]{$exonno}[3]
								&& abs($coords{$tranid}[0]{$exonno}[3]-$cend)<=$threshold) {
									$agree+=0.5;
								}
							}
								
							
						}
					}
					elsif ($coords{$tranid}) {
                                                $length+=0;
                                                $flag=1;
                                        }
				}
				$rem = $length % 3 unless ($exception);
				if ($flag) {
					# Store donated base number for transcripts with highest alignment agreement
					unless (exists $outhash{$gene}[0] && $outhash{$gene}[0]>$agree) {
						if (exists $outhash{$gene}[0] && $outhash{$gene}[0]==$agree) {
							$outhash{$gene}[0] = $agree;
							$outhash{$gene}[1] .= ",$rem";
							$outhash{$gene}[1] .= "_bp_donated" unless ($exception);
						}

						else {
							$outhash{$gene}[0] = $agree;
							$outhash{$gene}[1] = "$rem";
							$outhash{$gene}[1] .= "_bp_donated" unless ($exception);
						}
							
					}
				}
				$flag=0;
			}
		}


		# For processing 3` segment
		if ($segment%2==0) {
                        my $upperindex;
			my $exon;
			my $flag=0;
			my $agree;
			my $exception;
			my $rem;
				
			# Analyse each transcript corresponding to current gene
			foreach my $tranid (keys %{$ids{$gene}}) {
				$length=0;
				$agree=0;	
				$exception=0;
				$rem = undef;
				if (defined $coords{$tranid}[1] && $chrstart<=$coords{$tranid}[2] 
				&& $chrstart >= $coords{$tranid}[1]) {
                                        $rem= "$gene\_starts_in_5`_utr";
					$exception=1;
                                }


                                # Calculate joining exon number (upperindex) by determining which exon corresponds to the breakpoint
                                foreach my $exonno (sort {$a<=>$b} keys %{$coords{$tranid}[0]}) {
                                        my $start = $coords{$tranid}[0]{$exonno}[2];
                                        my $end = $coords{$tranid}[0]{$exonno}[3];
                                        if ( defined $coords{$tranid}[0]{$exonno}[3] && abs($coords{$tranid}[0]{$exonno}[3]-$chrstart)<=$threshold ||
                                        defined $coords{$tranid}[0]{$exonno}[2] && abs($coords{$tranid}[0]{$exonno}[2]-$chrstart)<=$threshold) {
                                                $upperindex = $exonno;
						$exon = $exonno;
                        			unless ($upperindex ==1) {$upperindex-=1;}
                                                last;
                                        }

                                }
                                # Move to process next transcript if no exons in the current one correspond to the breakpoint
                                next if (! defined $upperindex);


				# Figure out number of bases donated and number of exon boundaries that agree	
                        	foreach my $exonno (1..$upperindex) {
					if (defined $coords{$tranid}[0]{$exonno}[0]) {
						$flag=1;
						$length+= (($coords{$tranid}[0]{$exonno}[1] - $coords{$tranid}[0]{$exonno}[0])+1);
					}
					
					elsif ($coords{$tranid}) {
						$length+=0;
						$flag=1;
					}
				}

				foreach my $exonno (sort {$a<=>$b} keys %{$coords{$tranid}[0]}) {
					my @resultarray = split (/\n/, $chimers{$seqid});
					foreach my $result (@resultarray){
						my @linearray = split (/\s+/, $result);
						my $chrom = $linearray[2];
						my $cstart=$linearray[9];
						my $cend=$linearray[10];
						if ($cstart > $cend) {
                                                	($cstart, $cend) = ($cend, $cstart);
						}
						if ($chr eq $chrom) {
							if (defined $coords{$tranid}[0]{$exonno}[2]
							&& abs($coords{$tranid}[0]{$exonno}[2]-$cstart)<=$threshold) {
								$agree+=0.5;
							}
							if (defined $coords{$tranid}[0]{$exonno}[3]
							&& abs($coords{$tranid}[0]{$exonno}[3]-$cend)<=$threshold) {
								$agree+=0.5;
							}
						}


					}
				}
				unless ($exception) {
					if ($exon==1) {$rem = 0;}
					else {$rem = $length % 3;}
				}

				if ($flag) {
					# Store donated base number for transcripts with highest alignment agreement
                                	unless (exists $outhash{$gene}[0] && $outhash{$gene}[0]>$agree) {
                                                if (exists $outhash{$gene}[0] && $outhash{$gene}[0]==$agree) {
                                                        $outhash{$gene}[0] = $agree;
                                                        $outhash{$gene}[1] .= ",$rem";
							$outhash{$gene}[1] .= "_bp_expected" unless ($exception);
                                                }

                                                else {
                                                        $outhash{$gene}[0] = $agree;
                                                        $outhash{$gene}[1] = "$rem";
							$outhash{$gene}[1] .= "_bp_expected" unless ($exception);
                                                }

                                        }
				}
				$flag=0;
                        }
		}

	
		if (defined $outhash{$gene}[0]) {
			my @testarray = split (/,/, $outhash{$gene}[1]);
			my %testhash=();
			foreach my $element (@testarray) {
				$testhash{$element}++;
			}
			@testarray = keys %testhash;
			# If multiple transcripts match the alignment and don't all donate same no. of coding bases, can't determine frame prediction
			if (scalar @testarray > 1) { 
				if ($segment%2) {
					push @outtest, "$line Alignment_insufficient";
				}
				if ($segment%2==0) {
					push @outtest, "$line Alignment_insufficient";
					$outtest[0].=" Can't_determine\n";
					$outtest[1].=" Can't_determine\n";
					$unknown++;
					print OUTFILE @outtest;
					@outtest=();
				}
			}
			# If only a single transcript matches the alignment, output number of coding bases and frame prediction
			else {
				if ($segment%2) {
					push @outtest, "$line $testarray[0]";
				}
				if ($segment%2==0) {
					push @outtest, "$line $testarray[0]";
					my @line1 = split (/\s+/, $outtest[0]);
					my $comp1 = $line1[0-1];
					my @line2 = split (/\s+/, $outtest[1]);
                                        my $comp2 = $line2[0-1];
					if (($comp1=~/_bp_donated/)&&($comp2=~/_bp_expected/)) {
						$comp1=~s/_bp_donated//;
						$comp2=~s/_bp_expected//;
						if ($comp1==$comp2) {
							$outtest[0].=" Yes\n";
                                       		 	$outtest[1].=" Yes\n";
							$fusion++;
                                        		print OUTFILE @outtest;
                                        		@outtest=();
						}
						else {
							$outtest[0].=" No\n";
                                                        $outtest[1].=" No\n";
							$junk++;
                                                        print OUTFILE @outtest;
                                                        @outtest=();
						}
					}
					else {
						$outtest[0].=" Can't_determine\n";
                                        	$outtest[1].=" Can't_determine\n";
						$unknown++;
                                        	print OUTFILE @outtest;
                                        	@outtest=();
					}
				}

			}
		}
		# Can't determine result if gene is unknown
		else {
			if ($segment%2) {
				push @outtest, "$line Unknown gene";
			}
			if ($segment%2==0) {
				push @outtest, "$line Unknown gene";
				$outtest[0].=" Can't_determine\n";
                                $outtest[1].=" Can't_determine\n";
				$unknown++;
                                print OUTFILE @outtest;
                                @outtest=();
			}
		}
	}

}


# Print processing time info.
my ($user,$system,$cuser,$csystem) = times;
print "\n** Processing finished.  User time=$user s **\n";

# Print run stats
my $accs = keys %uniqaccs;
print "\n\tFrameChecker statistics\n\n";
print " No. input lines processed\t\t\t\t:$lines{nlines}\n";
print " No. distinct accession numbers\t\t\t\t:$accs\n";
print " No. of in-frame fusions predicted\t\t\t:$fusion\n";
print " No. of out-of-frame fusions predicted\t\t\t:$junk\n";
print " No. with unknown result (insufficient alignment data)\t:$unknown\n\n";

unless (($accs) == ($fusion+$unknown+$junk)) {print "Number discrepancy!\n\n";}


# Subroutine to print help to screen
sub usage {

    print "\n\t$0 version $version\n";
    print "\tby G.R. Oliver\n";
    print "\nUsage:\n\n";
    print "\t $0 -i input filename (.accept file) -c coordinate file -r results (.filt.srted file)\n";
    print "\t e.g $0 -i example.accept -c framecheck.coords -r results.filt.srted\n\n";
    print "Options:\n";
    print "\t -t <maximum discrepancy (in bp) allowed when checking exon boundaries (default 3)>\n";
    print "\t -d <output-dir>\n";
    print "\t -h <print this screen>\n";
    print "NB: Input is .accept or .accept.clustered file from FusionMiner\n";

}

