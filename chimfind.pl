#!/usr/bin/perl

#  
# Title:        Chimfind (Part of FusionMiner)
#
# Function:     In Silico transplicing program - selects initial (loose) chimer candidates
#
# Author:       Andy Emerson, Gavin Oliver
# 
# Date:         27 May 2009
#
# Version:      V1.0 
#
# Copyright (C) <2010>  <Andy Emerson, Gavin Oliver>
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

my $version="27-may-2009 1.0";


my ($blastfile,@linebuffer,$lineno,$line,$linefield,$finished);
my ($query1,$query2);
my (@id_lines,$id);
my (@linefields,@fields,@start,@sorted_lines,@accepted_lines);
####
my @linefields2;
####
my ($i,$lo,$hi,$start,$end,$start2,$end2,$reject);
my (@chrs,$chr,$chr1,$chr0,@all_chrs);
my (@indices);
my ($cstart, $cend, $chi, $cend2);
my (@sorted_chr_lines,@cat_lines,@all_cat_lines,@lines);

# positions of the blast fields when using blast -m8
my %blast = (acc=>0, chr=>1, pcid=>2, start => 6, end =>7, cstart => 8, cend =>9 );
my ($pcid,$acc);
my $def_overlap=10; # default overlap in bases
my $overlap=$def_overlap;
    my ($key,$value);     


my %chrs;


my ($outdir,$mult_accept_file,$sing_accept_file,$clsing_accept_file,$reject_file,$fname,$stem,$ext);
my ($mult_accno_file,$sing_accno_file);

my $accept;

# stats variables
my %stats = (nlines =>0,  # no. of input lines
             naccno =>0, # no. distinct acession numbers processed
             nrej1 =>0, # no. acc rejected due to 1 match on 1 chromosone
             nrej3 =>0, # no. reject on 3 or more chrs
             nacc =>0,   # no. accepted with hits on 2 chromosomes
	     nacc2=>0 	# no. accepted with hits on 1 chromosome
             );   
my ($user,$system,$cuser,$csystem);

#   get command line options
my %opts;

&getopts('ho:d:i:',\%opts);


# if help requested print and exit
if (exists $opts{h}) {
    &usage;
    exit;
}

# Quit unless input file is defined on command line
if (defined $opts{i}) {
	$blastfile = $opts{i};
}

else {&usage;
	print "\n\nInput filename must be specified!\n\n";
	exit;
}


print "\n\t$0: version $version. Use -h for options\n\n";


# if overlap has been specified change from default
if (defined $opts{o} ){
    $overlap=$opts{o};
}
print "Overlap set to $overlap bases\n";

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

# reject and accept files should use the file stem of the blastfile

# remove path component
$fname =$blastfile;
$fname =~ s/(.*)\///;

#remove extension
($stem,$ext)=split(/\./,$fname,2);

$mult_accept_file = $stem."-multi.accept";
$sing_accept_file = $stem."-sing.accept";
$clsing_accept_file = $stem."-clsing.accept";
$reject_file = $stem.".reject.cfind";
$mult_accno_file  = $stem."-multi.accno";
###Remove this file?
$sing_accno_file = $stem."-clsing.accno";

open(ACCEPT,">$outdir/$mult_accept_file") or die "Can't open $outdir/$mult_accept_file for writing";
open(ACCEPT2,">$outdir/$sing_accept_file") or die "Can't open $outdir/$sing_accept_file for writing";
open(ACCEPT3,">$outdir/$clsing_accept_file") or die "Can't open $outdir/$clsing_accept_file for writing";
open(REJECT,">$outdir/$reject_file") or die "Can't open $outdir/$reject_file for writing";
open(ACCNO,">$outdir/$mult_accno_file") or die "Can't open $mult_accno_file for writing\n";
open(ACCNO2,">$outdir/$sing_accno_file") or die "Can't open $sing_accno_file for writing\n";
print "Writing multi-chromosome accepted alignments to $outdir/$mult_accept_file\n";
print "Writing single-chromosome accepted alignments (clustered) to $outdir/$sing_accept_file\n";
print "Writing single-chromosome accepted alignments (unclustered) to $outdir/$clsing_accept_file\n";
print "Writing list of multi-chromosome accepted accession numbers to $outdir/$mult_accno_file\n";
print "Writing list of single-chromosome accepted accession numbers to $outdir/$sing_accno_file\n";
print "Writing rejected sequences to $outdir/$reject_file\n";
print ACCEPT "#query_id chromosome %identity query_start query_end chrom_start chrom_end\n";
print ACCEPT2 "#query_id chromosome %identity query_start query_end chrom_start chrom_end\n";
print ACCEPT3 "#query_id chromosome %identity query_start query_end chrom_start chrom_end\n";
#print REJECT "#query_id chromosome  %identity  query_start  query_end chrom_start chrom_end\n";
  

open (FILTERED,"$blastfile") or die "Cannot find/read blast output file=$blastfile\n";

print "Analysing  file $blastfile ..\n";
print "(In output files * indicates a concatenated match, %identities are weighted averages)\n\n";
@linebuffer=();

$lineno=0;
$finished=0;

while (!$finished && ($line=&getline(\*FILTERED)) ) {
    chomp($line);
    $stats{nlines}++;

    # discard any comment lines that might be lurking
    next if ($line =~ /^\#/);
    $lineno++;

    # reset the arrays which store the line and fields for each identifier
    @linefields=();  # this will be an array of arrays
####    
    @linefields2=(); # this array initially stores the same as @linefields but is not altered by concatenation etc
####   
    # separate current line into fields
    @fields=split(/\s+/,$line);
    # Check input file format
    unless ($stats{naccno}) {
        unless ($#fields==11) {
                die "***Quitting!***  Input file must be tabular BLAST output (-m 8 for BLAST / -out=blast8 for BLAT)\n\n";
        }
    }
    $query1 = $fields[0];   # this is the query identifier
    $stats{naccno}++;

    # store the fields array in @linefields/@linefields2
    push @linefields,[@fields];
    push @linefields2, [@fields];

    while ($line=&getline(\*FILTERED)) {
       chomp($line);
       $stats{nlines}++;
       @fields= split(/\s+/,$line);
       $query2 = $fields[0];
       
       # compare : if not the same identifier get out, 
       # remembering to put line into linebuffer because we have read 1 line too many
       if ($query2 ne $query1) {
           push (@linebuffer,$line);
           $stats{nlines}--;
           last;
       }
       # store the fields of this line
       push @linefields, [@fields];
####
       push @linefields2, [@fields];
####
    }
    $finished =1 if (!$line); # finish if end of file

    # check to see if we have found 2 or more records with the same 
    # identifier. If not carry on looping through file
 

   if (@linefields==1) {  # should this be <=1 
      # no, so put line back and start comparing again
      #push (@linebuffer,$line);
      $stats{nrej1}++;    
	####
	foreach my $entry (@linefields) {
                print REJECT "@$entry\n";
        }
	####
      next;  
    }

	
    # ok now we have a set of matches for an identifier (one or more chromosomes).
    # Now separate lines according to chromosome, collect first all the chromosomes.
   
    @chrs=();
    foreach $line (@linefields) {
       $chr =($line->[1]);    
       push @chrs,$chr;
   }


    # sort (i.e. group together) linefields according to chromosone
    @sorted_chr_lines =&field_sort(\@linefields,\@chrs,"string");
  

 # now collect lines for each chromosone:
  
    # reset %chr first
    %chrs=();
    foreach $line (@sorted_chr_lines) { 

       $chr=($line->[1]);
       push @{$chrs{$chr}},$line;

   }


    @all_cat_lines=();
    foreach $chr (keys %chrs) {
	@lines=@{$chrs{$chr}};

        # if only 1 line store the line and forget about it
        if (@lines==1) {
           push (@all_cat_lines,shift @lines);
           next;
        }
        # otherwise try to concatenate
        @cat_lines= &concat(\@lines,$overlap);
        push(@all_cat_lines,@cat_lines);
    }      

        
    # now perform an elimination step between all chrs
    @cat_lines=@all_cat_lines;
    @all_cat_lines=&eliminate(\@cat_lines);


    # now re-collect all chromosones 
    %chrs=();
    @all_chrs=();
    foreach $line (@all_cat_lines) {
       $chr =($line->[1]); 
       push @all_chrs,$chr;   
       push @{$chrs{$chr}},$line;
   }
    @chrs =keys %chrs;

# Check if we have matches on only one chromosome

    if (&one_chr(\@linefields)) {
	print ACCNO2 "$query1\n";
        $stats{nacc2}++;
####    
        #foreach my $entry (@linefields2) {
        #        $acc = $entry->[$blast{acc}];
        #        $pcid  = $entry->[$blast{pcid}];
        #        $chr = $entry->[$blast{chr}];
        #        $start = $entry->[$blast{start}];
        #        $end = $entry->[$blast{end}];
        #        $cstart = $entry->[$blast{cstart}];
        #        $cend = $entry->[$blast{cend}];
        #        print ACCEPT2 "$acc $chr $pcid $start $end $cstart $cend\n";
        #}

	#print ACCEPT2 "\n";
	my @lf2starts=();
	my @linefields2srtd=();
	
	foreach my $entry (@linefields2) {
		push @lf2starts, $entry->[$blast{start}];
	}
	
	@linefields2srtd = &field_sort(\@linefields2,\@lf2starts,"numeric");

	foreach  my $entry (@linefields2srtd) {
                $acc = $entry->[$blast{acc}];
                $pcid  = $entry->[$blast{pcid}];
                $chr = $entry->[$blast{chr}];
                $start = $entry->[$blast{start}];
                $end = $entry->[$blast{end}];
                $cstart = $entry->[$blast{cstart}];
                $cend = $entry->[$blast{cend}];
                print ACCEPT2 "$acc $chr $pcid $start $end $cstart $cend\n";
        }
	
	print ACCEPT2 "\n";


	foreach my $entry (@all_cat_lines) {
                       my $chrm = $entry->[$blast{chr}];
                        if ($chrs[0] eq $chrm) {
                                $acc = $entry->[$blast{acc}];
                                $pcid  = $entry->[$blast{pcid}];
                                $start = $entry->[$blast{start}];
                                $end = $entry->[$blast{end}];
                                $cstart = $entry->[$blast{cstart}];
                                $cend = $entry->[$blast{cend}];
                                print ACCEPT3 "$acc $chr $pcid $start $end $cstart $cend\n";
                        }
        }
        print ACCEPT3 "\n";
####
        next;
    }


   #  chimera criterion:
   #  1 chr -> reject
   #  2 chrs -> accept
   # >2 chrs -> reject
    if ($#chrs ==0) { # Matches on only 1 chromosone
        $stats{nacc2}++;
        $accept=0;
	####
	print ACCNO2 "$query1\n";

	my @lf2starts=();
        my @linefields2srtd=();

        foreach my $entry (@linefields2) {
                push @lf2starts, $entry->[$blast{start}];
        }

        @linefields2srtd = &field_sort(\@linefields2,\@lf2starts,"numeric");

	foreach my $entry (@linefields2srtd) {
			my $chrm = $entry->[$blast{chr}];
			if ($chrs[0] eq $chrm) {
				$acc = $entry->[$blast{acc}];
				$pcid  = $entry->[$blast{pcid}];
				$start = $entry->[$blast{start}];
				$end = $entry->[$blast{end}];
				$cstart = $entry->[$blast{cstart}];
				$cend = $entry->[$blast{cend}];
				print ACCEPT2 "$acc $chr $pcid $start $end $cstart $cend\n";
			}
	}

	print ACCEPT2 "\n";

	foreach my $entry (@all_cat_lines) {
                       my $chrm = $entry->[$blast{chr}];
                        if ($chrs[0] eq $chrm) {
                                $acc = $entry->[$blast{acc}];
                                $pcid  = $entry->[$blast{pcid}];
                                $start = $entry->[$blast{start}];
                                $end = $entry->[$blast{end}];
                                $cstart = $entry->[$blast{cstart}];
                                $cend = $entry->[$blast{cend}];
                                print ACCEPT3 "$acc $chr $pcid $start $end $cstart $cend\n";
                        }
        }
	print ACCEPT3 "\n";
	####
	
    }
	elsif ($#chrs ==1){ # Matches on 2 chromosones
        $stats{nacc}++;
	$accept=1;
    } 
    	else { # Matches on 3 or more chromosomes
    	$stats{nrej3}++;
    	$accept=0;
	foreach my $entry (@linefields2) {
                        print REJECT "@$entry\n";
                }
}
    # re-sort according to chr (Not really necessary?)
     @all_cat_lines=&field_sort(\@all_cat_lines,\@all_chrs,"string");  


       #  finally print out
    print ACCNO "$query1\n" if ($accept);
        foreach $line (@all_cat_lines) {
            $acc = $line->[$blast{acc}];
            $pcid  = $line->[$blast{pcid}];
            $chr = $line->[$blast{chr}];
            $start = $line->[$blast{start}];
            $end = $line->[$blast{end}];
            $cstart = $line->[$blast{cstart}];
            $cend = $line->[$blast{cend}];
	    print ACCEPT "$acc $chr $pcid $start $end $cstart $cend\n" if ($accept);
        }
        print ACCEPT "\n" if ($accept);
	#print REJECT "\n" if (!$accept);

}  #  end main while over input file

($user,$system,$cuser,$csystem) = times;
print "\n** Processing finished.  User time=$user s **\n";

#  print statistics

print "\n\tchimfind statistics\n\n";
print " No. input lines processed\t\t\t\t:$stats{nlines}\n";
print " No. distinct accession numbers\t\t\t\t:$stats{naccno}\n";
print " No. acc numbers with hits on 1 chr rejected\t\t:$stats{nrej1}\n";
print " No. acc numbers with hits on 3 or more chr rejected\t:$stats{nrej3}\n";
print " No. acc numbers with hits on 1 chromosome accepted\t:$stats{nacc2}\n";
print " No. acc numbers with hits on 2 chromosomes accepted\t:$stats{nacc}\n\n";

if ( ($stats{nacc}+$stats{nrej1}+$stats{nrej3}+$stats{nacc2})!=$stats{naccno}) {
   print "**WARNING: Naccepted + Nrejected <> N accession numbers\n\n";
} 


close(ACCEPT);
close(REJECT);
close(ACCNO);

close(FILTERED);


sub field_sort {

# sorts an array according to a set of keys

    my @arr =@{$_[0]};
    my @keys =@{$_[1]};
    my $flag=$_[2]; # set to "string" for a string comparison (othewise numeric)

    my %tmp;
    my @sorted_lines=();

    for $i (0..@keys-1) {
	$tmp{$i}=$keys[$i];
    }
    if ($flag eq "string") {
       @indices = sort {$tmp{$a} cmp $tmp{$b}} keys %tmp;
    } else {
       # numeric comparison
       @indices = sort {$tmp{$a} <=> $tmp{$b}} keys %tmp;
    }
    for $i (0..@indices-1) {
       $sorted_lines[$i] = $arr[$indices[$i]];
   }
    @sorted_lines;
}


sub getline {

    my $fh=shift;
    my $line;
    if (@linebuffer == 0 ) {
        $line =<$fh>;      
        return 0 if (!$line);
        push (@linebuffer,$line);
    }
    $line = pop (@linebuffer);
    #print "3. getline: line=$line";
    $line;
}


sub one_chr {
    use strict;
    my @linefields=@{$_[0]};
    my ($line,$chr,$chr0,$n);
    my $onechr=1;
    my @chr_fields;

    $line = shift @linefields;
    $chr0= ($line->[1]);
    
    foreach $line (@linefields) {
	$chr = ($line->[1]);
         if ($chr ne $chr0) {
	    $onechr=0;
	    return $onechr;
	}
    }
    return $onechr;

}


sub concat {

    my @lines=@{$_[0]};
    my $overlap=$_[1];

    my @accepted_lines;
    my ($i,$line,$hi,$start2,$end2,$lo, $chi, $cend2);
    my ($pcid,$spcid,$avpcid,$npcid);
    my ($nbases,$snbases);

    my @concat_lines;
    my $reject;
    my $acc;


     # eliminate matches conatined wholly within other matches   
     @accepted_lines=&eliminate(\@lines);

 
     # @accepted_lines now contains a set of matches, possibly overlapping.
     # Using the overlap criterion try to fit them together
     # NB: the %identity of the concatenated lines will be a weighted match
     # of the component %ids, e.g.
     # <%> = (nb1*%1 + nb2*%2 + nb3*%3..)/(nb1+nb2+nb3+..)
     # where nb1=no.of bases for fragment1
 
     @concat_lines=();
     $line = $accepted_lines[0];
     push @concat_lines,$line;
    $acc = $line->[$blast{acc}];  # in fact this is the same for all lines here
     $hi = $line->[$blast{end}];
    ###
    $chi = $line->[$blast{cend}];
    $lo = $line->[$blast{start}];
    $nbases = $hi-$lo+1;
    $snbases = $nbases;
    $pcid =$line->[$blast{'pcid'}]*$nbases; # weighted %
    $spcid=$pcid;

     for $i (1..$#accepted_lines) {
	 $line = $accepted_lines[$i];
	 $start2 = $line->[$blast{start}];
	 $end2 = $line->[$blast{end}];
	 ###
	 $cend2 = $line->[$blast{cend}];
         $nbases = $end2-$start2+1;
         $pcid = $line->[$blast{pcid}]*$nbases;         
 
         if ( abs($start2-$hi)<=$overlap)  {
           # concatenate the matches
           $hi=$end2;
	   ###	
	   $chi=$cend2;
           $snbases += $nbases;
           $spcid +=$pcid;
	   $concat_lines[$#concat_lines]->[$blast{end}]="$hi";
	   ###
	   $concat_lines[$#concat_lines]->[$blast{cend}]="$chi";
           $avpcid=sprintf("%.2f",$spcid/$snbases);          
           $concat_lines[$#concat_lines]->[$blast{pcid}]=$avpcid;
          
           # mark this line as concatenated
           $concat_lines[$#concat_lines]->[$blast{acc}] = $acc."*";
           next;
         } 
        # otherwise store new line, reset hi and spcid,snbases
	 push @concat_lines,$line;
         $hi=$end2;
	 ###
	 $chi=$cend2;
         $spcid=$pcid;
         $snbases=$nbases;

     }
     @concat_lines;     
}


# eliminate matches contained within other matches
#
sub eliminate {

    my @lines=@{$_[0]};

    my ($line,$lo,$hi,$i,$start2,$end2);
    my (@accepted_lines,@sorted_lines);
    my @start=();
    my $reject;

    # collect start positions and sort lines according to start
    foreach $line (@lines) {
	push @start,$line->[$blast{start}];
    }
    @sorted_lines=&field_sort(\@lines,\@start,"numeric");

   # now eliminate matches totally within other matches

     @accepted_lines=();
     $line = $sorted_lines[0];
     push @accepted_lines,$line;
     $lo = $line->[$blast{start}];
     $hi = $line->[$blast{end}];
 
     for $i (1..$#sorted_lines) {
	 $line = $sorted_lines[$i];
	 $start2 = $line->[$blast{start}];
	 $end2 = $line->[$blast{end}];
	  
            # Backwards check:
            # check if previously accepted sequence sits within
            # current sequence (NB: we know $start2 >= $lo) 

	    if ($start2 == $lo && $end2 >=$hi) {
		# remove line corresponding to lo and hi
		# make new lo,hi = start,end
                $reject = pop @accepted_lines;

            # Forwards check:
            # check if current sequence sits within previous sequence
            # (we know $start2 >$lo)
	    } elsif ($start2<$hi && $end2 <= $hi) {
		# reject current line, keep lo and hi as before
                next;
	    }
	    # otherwise accept line
            push @accepted_lines,$line;
            $lo = $start2;
            $hi = $end2; 

	}        
    @accepted_lines;
}



sub usage {

    print "\n\t$0 version $version\n";
    print "\tby A. Emerson & updated by G.R. Oliver\n";
    print "\nUsage:\n";
    print "\t $0 -i blastoutput file\n";
    print "\t e.g $0 -i example.filtered.clean.sorted \n\n";
    print "Other options\n";
    print "\t -o <gap/overlap threshold for hit joining (default 10)>\n";
    print "\t -d <output-dir>\n";
    print "NB: blast output should have been filtered (blastfilter) and sorted (UNIX sort)\n";

}


sub show_version {
    print "$0: Version $version \n\n";
}



sub filter {


    use strict;
    my $cutoff=shift;
    my $field_number = shift;
    my $infile=shift;
    my $outfile = shift;

    my $line;
    my ($query,@fields);

    open(IN,$infile) or die "Cannot find/read blast results file=$infile\n";    
    open(OUT,">$outfile") or die "Cannot open filter file=$outfile\n";    
    while ($line=<IN>) {
       chomp($line);
       @fields=split(/\s+/,$line);
       $query=$fields[0];
       if ($fields[$field_number]<$cutoff) {
#	   print "Eliminating $query with $fields[$field_number] < $cutoff \n";
       } else {
	   print OUT "$line\n";
       }
   } # end while
    close (IN);
    close (OUT);
}


#
# routine for debugging only
#
sub debug {

    my @arr = @{$_[0]};
    my $str = $_[1];
    my $id  = $_[2];

    my $a;
    print "$str\n";
    foreach $a (@arr) {
	print "$id @$a \n";
    }

}

