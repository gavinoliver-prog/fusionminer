Fusionminer Source README


1)  The only script that needs to be run is the "FusionMiner.pl" script.  Others are run automatically.

2)  Execution paramters are stored in the FusionMiner.param file.  These can be altered but it is not recommended for novice use.

3)  The program accepts genomic alignments in BLAST/BLAT tabular ouput (BLAST command line option –m=8; BLAT command line option out=blast8).  

4)  Coordinate files are based on human genome build GRCh37 (ucsc hg19).  Therefore alignments to this genome version are expected.
    Furthermore the program assumes chromsomes are named in the format chr1, chr22, chrX etc.

5)  The source version is written in Perl therefore Perl is required.  The Getopt::Std modules is needed also.

6)  A test set of results (THC.result.zip) are provided. This should be uncompressed before use.

7)  Genomic coordinates used as input to the program are provided in the "genomic.coords" file.
    The coordinates required for reading frame detrmination are provided in the "framecheck.coords" file.
    Both are accessed automatically by the program.

8)  Run the program from within the FusionMiner directory.
 
9)  Typing ./FusionMiner.pl -i THC.result -p FusionMiner.param will execute the program.

10)  Final accepted chimer result files will be found in the 'Accepted' directory following program completion.  
    One file will exist for intrachromosomal chimers and another for interchromosomal chimers.
    These will be suffixed *-sing.accept.clustered.framechecked and *-multi.accept.clustered.framechecked respectively.

11) Transient files are stored in the 'Archived' directory.  The identity of each transient file is contained in the program's run-time screen output.
    These should not be required in most instances.

12) The framecheck contained in the final result files considers the alignment upstream and downstream of the fusion point and attempts to determine
    if reading frame is preseved at the fusion point.

    As numerous different transcripts can arise from a single gene and these can be similar in sequence but different in coding potential, the software
    attempts to count how many exon boundaries the alignment agrees with for each transcript.  When a favorite is found (i.e. most exon boundary agreements)
    this transcript is assumed to be the one we represent.  If the alignment is equally accurate for a number of transcripts and each of these occurs in the
    same reading frame then this value is used.  Otherwise the software does not make a prediction.  For the situations where a favorite is found, the number 
    of upstream coding bases are counted.  If this number is divisible by 3, we know the upstream fusion partner ‘donates’ 0 bases as coding occurs in triplets.  
    If the number is not divisible by 3, a remainder of 1 or 2 bases means that the upstream partner ‘donates’ 1 or 2 bases respectively.

    The missing part of the downstream fusion partner (i.e. the sequence that is replaced by its upstream fusion partner) is similarly assessed to determine how
    many bases it naturally donates at this point.  If the number matches the number donated by the upstream fusion partner, we know that reading frame is
    preserved, and therefore whatever translation occurs downstream of the fusion will not be junk translation.

    If the upstream fusion ends at that transcript’s known, utilized stop codon, the program recognizes and outputs this.  Similarly it recognizes instances
    where the join occurs within 5` or 3` UTRs. 

    If 0 bases are donated/expected at the fusion point, no stop codon should occur until the downstream partner's known canonical stop codon.

    If 1 or 2 bases are donated/expected respectively, there is some potential for a premature stop codon to be introduced, but a higher liklihood
    of a single erroneous amino acid to be created and translation to continue.

    This added functionality may be of interest to researchers looking specifically for chimers with/without protein-coding potential.

    Further assessment of precisely what happens in terms of translation can easily be performed with a 6 frame-translator program.

 
13) Alignments rejected at the exon checking stage will be mostly junk but may contain non-canonical chimers.  
    These files are stored in the 'Archived' directory and are suffixed *-multi.reject.sscan and *-sing.reject.sscan.
    Users interested in these results may wish to inspect these files.  
    

    One recommended way of looking for genuine non-canonical fusions amongst the junk alignments is to look for fusions that occur more than once.
    Currently this functionality is only implemented for interchromosmal chimers.

    This is most easily done by running the chimclust.pl script as follows:

    ./chimclust.pl -i Archived/filename-multi.reject.sscan	(for the interchromosmal candidates)

    and inspecting the output file.  Results will be ordered by number of repeat occurrences of a potential chimer.

    Currently, these non-canonical candidates are not annotated in terms of gene name etc.  This functionality may be added in a later version.


Questions can be directed at gavin.oliver@almacgroup.com
