#!/usr/bin/perl

# Title:        Blastfilter (Part of FusionMiner)
#
# Function:     Perl program to filter blast/blat tabular output according to score, length or % identity
#               Any combination of the above three allowed.
#
# Authors:      A. Emerson, G.R. Oliver 
#
# Date:         17 April 2009
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



use warnings;
use strict;
use Getopt::Std;

my ($line,$lines_excluded, $lines_processed); # Current line, no of lines excluded and no of lines processed
my ($percent_cutoff,$score_cutoff,$length_cutoff); # filter cutoffs
my $print_line; # Flag to print or reject line
my @fields; # array to split columns into
my %blast_col;# hash tro specify column positions of filter criteria
my $infile; # Input filename
my $outfile; # Output filename
my $version = "17-apr-2009 1.0"; # version
my $outdir; # output directory

# Get command line options and arguments
my %opts;
&getopt('pslhi:o:d:',\%opts);

# If help flag given print usage
if (exists $opts{h}) {
	&usage;
	exit;
}

# Unless input file specified print usage and exit
if (defined $opts{i} && defined $opts{o}) {
	$infile = $opts{i};
	$outfile = $opts{o};
}
else {
	&usage;
	exit;
}

$outdir=".";
# check for output directory
if (defined $opts{d}) {
    $outdir=$opts{d};
    if (! -e $outdir){
        print "Problem with output directory $outdir\n";
        exit;
    }
    print "Output directory = $outdir\n";
}

# Set thresholds and lines rejected to zero	
$lines_excluded=0;
$percent_cutoff =0;
$score_cutoff=0;
$length_cutoff=0;
   
# Accept command line definitions for thresholds if given
$percent_cutoff = $opts{p} if defined($opts{p});
$score_cutoff = $opts{s} if defined($opts{s});
$length_cutoff = $opts{l} if defined($opts{l});

# Set column numbers for filter criteria
%blast_col = ( percent => 2,score => 11, length => 3 );

# Open input and output files
open (INFILE, $infile) or die $!;
open (OUTFILE, ">$outdir/$outfile") or die $!;

# Print general info messages to screen

print "\n\t$0: version $version. Use -h for options\n\n";

print "Reading from file $infile\n";
print "Writing to file $outfile\n";
print "Percent ID cutoff set to $percent_cutoff%\n";
print "Score cutoff set to $score_cutoff\n";
print "Length cutoff set to $length_cutoff bases\n\n";


# Filter result lines based on cutoffs
while (<INFILE>) {
	$line = $_;
	chomp($line); 
	next if ($line =~ /^#/);    # skip comment lines
	@fields=split(/\s+/,$line);
	# Check inpiut file format
	unless ($lines_processed) {
		unless ($#fields==11) {
		die "***Quitting!***  Input file must be tabular BLAST output (-m 8 for BLAST / -out=blast8 for BLAT)\n\n";}
	}
	$print_line=1;
	if ($percent_cutoff && $fields[$blast_col{percent}]<$percent_cutoff) {
		$print_line=0;
	} 
	elsif ($score_cutoff && $fields[$blast_col{score}]<$score_cutoff) {
		$print_line=0;
	} 
	elsif ($length_cutoff && $fields[$blast_col{length}]<$length_cutoff) {
		$print_line=0;
	} 
	else {	
		print OUTFILE "$line\n";
	}
	$lines_processed++;
	$lines_excluded++ if (!$print_line);
}  # end while


# Print processing time info.

my ($user,$system,$cuser,$csystem) = times;
print "\n** Processing finished.  User time=$user s **\n";

print "\n\tblastfilter statistics\n\n";
print "No. lines processed\t\t\t:$lines_processed\n";
print "No. lines filtered\t\t\t:$lines_excluded\n";
print "No. lines accepted\t\t\t:",$lines_processed - $lines_excluded,"\n\n";

sub usage {

	print "usage: $0 -i <blast/blat tabular output file> -p <per-cent> -l <length> -s <score> -o <output file>\n\n";
	print "Examples:\n\n";
	print "$0 -i blast.dat -o filtered_blast.out -p 100\n";
	print "$0 -i blast.dat -l 10 -s 100 -o out1\n";


}
