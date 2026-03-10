# FusionMiner

**FusionMiner** is a Perl-based pipeline for detecting and characterizing candidate fusion transcripts (chimeric RNAs) from BLAST or BLAT-created genomic alignment data, originally written way back in 2009 when NGS was only beginning to gain traction.

The software was developed at the **beginning of the next-generation sequencing era**, but was designed to operate on **Sanger sequencing–derived transcript data and alignment output**, rather than modern RNA-seq read datasets.

This repository preserves the original code as a **historical reference implementation** associated with the following publication:

**Plebani, R., Oliver, G. R., Trerotola, M., Guerra, E., Cantanelli, P., Apicella, L., Emerson, A., Albiero, A., Harkin, P. D., Kennedy, R. D., & Alberti, S.** (2012).  
*Long-range transcriptome sequencing reveals cancer cell growth regulatory chimeric mRNA.*  
**Neoplasia**, 14(11), 1087–1099.  
https://doi.org/10.1593/neo.121342

\*Authors contributed equally.


---

# Historical Context

FusionMiner reflects an early approach to fusion transcript discovery using transcript-derived sequence alignments. Rather than operating directly on sequencing reads, the pipeline processes **BLAST or BLAT alignment output** to identify candidate chimeric transcripts.

The software was developed during a period when computational pipelines for detecting fusion transcripts were still evolving. FusionMiner implemented deterministic filtering and annotation steps designed to identify biologically plausible fusion events and characterize their potential coding consequences.

The pipeline was used in the Neoplasia study referenced above, which used long-range transcriptome sequencing to identify and characterize chimeric mRNAs associated with cancer cell growth.

---

## Repository Purpose

This repository preserves the original FusionMiner source code as it existed during its use in the associated research publication. The software is provided primarily for historical reference, reproducibility, and archival purposes.

FusionMiner was designed for transcript-derived sequence alignments generated from long transcript sequencing approaches and Sanger-based workflows common at the time. As sequencing technologies evolved, fusion detection pipelines shifted toward tools designed specifically for high-throughput RNA-seq data and large-scale datasets.

This repository therefore represents an example of an early computational approach to systematic fusion transcript discovery from alignment data and may be of interest to researchers studying the evolution of fusion detection methods or reproducing analyses from the original study.

---

# Key Features

Although developed more than a decade and a half ago, FusionMiner included several analytical capabilities that addressed challenges in early fusion transcript discovery.  Miraculously I believe it still works out-of-the-box. 

## Fusion Candidate Detection

FusionMiner processes genomic alignments in tabular BLAST/BLAT format and identifies candidate chimeric alignments consistent with transcript fusion events.

## Reading Frame Preservation Analysis

The pipeline attempts to determine whether predicted fusion transcripts maintain a protein reading frame.

This is done by:

- Mapping alignment positions relative to exon boundaries
- Selecting the transcript model with the greatest agreement with exon boundaries
- Calculating codon offsets between upstream and downstream fusion partners

This allows users to prioritize fusions with potential protein-coding capability.

## Transcript Model Evaluation

Because genes may produce multiple transcripts with different exon structures, the pipeline evaluates alignment agreement across transcript models and selects the transcript with the most exon boundary matches before performing frame calculations.

## Fusion Clustering

Repeated occurrences of similar fusion events can be clustered to highlight potentially recurrent fusion candidates.

## Classification of Fusion Types

FusionMiner separates results into:

- **Intrachromosomal chimeras**
- **Interchromosomal chimeras**

These are written to separate output files for easier interpretation.

## Inspection of Non-Canonical Fusions

Alignments rejected during exon boundary filtering are retained so users can manually inspect potential **non-canonical or rare fusion events**.

A clustering script is provided to identify repeated candidate events among rejected alignments.

---

# Requirements

FusionMiner was designed for Linux environments typical of early bioinformatics pipelines.

Required software:

- **Perl**
- Perl module: `Getopt::Std`

---

# Input Data

FusionMiner accepts genomic alignments in tabular format produced by:

### BLAST

```
-m 8
```

### BLAT

```
out=blast8
```

Genome reference assumptions:

- Reference genome: **GRCh37 (UCSC hg19)**
- Chromosome naming convention:

```
chr1, chr2, chr3 ... chrX
```

---

# Example Dataset

A test dataset is included:

```
THC.result.zip (https://sourceforge.net/projects/fusionminer/files/latest/download)
```

Uncompress the file before running the pipeline.

---

# Running the Pipeline

The main executable script is:

```
FusionMiner.pl
```

Other scripts are executed automatically by the pipeline.

Run the program from within the FusionMiner directory.

Example:

```
./FusionMiner.pl -i THC.result -p FusionMiner.param
```

Parameters:

| Option | Description |
|------|-------------|
| `-i` | Input alignment result set |
| `-p` | Parameter configuration file |

---

# Configuration

Execution parameters are stored in:

```
FusionMiner.param
```

These parameters can be modified but were not intended for novice users.

---

# Coordinate Files

Two coordinate files are used internally:

```
genomic.coords
framecheck.coords
```

These provide genomic coordinate and transcript annotation information used during fusion frame determination.

---

# Output

Final accepted fusion candidates are written to the:

```
Accepted/
```

directory.

Output files include:

```
*-sing.accept.clustered.framechecked
*-multi.accept.clustered.framechecked
```

These correspond to:

- intrachromosomal fusion candidates
- interchromosomal fusion candidates

---

# Intermediate and Archived Files

Intermediate and transient files are stored in:

```
Archived/
```

Rejected alignments are also stored here and may contain non-canonical fusion candidates.

Examples:

```
*-multi.reject.sscan
*-sing.reject.sscan
```

---

# Reading Frame Determination

FusionMiner evaluates whether a fusion preserves protein coding potential by examining the upstream and downstream transcripts.

The process includes:

1. Identifying the transcript model that best matches exon boundary alignments.
2. Counting upstream coding bases relative to the fusion point.
3. Determining how many bases (0, 1, or 2) are donated by the upstream partner.
4. Determining how many bases are expected from the downstream partner.
5. Comparing the two offsets.

If the offsets match, the reading frame is predicted to be preserved.

The pipeline also identifies cases where the fusion junction occurs within:

- 5′ UTR
- 3′ UTR
- the canonical stop codon of the upstream transcript

These annotations allow researchers to identify fusion transcripts with potential biological or translational consequences.

---

# Inspecting Rejected Alignments

Rejected alignments may still contain genuine non-canonical fusion events.

Users can inspect these using the clustering script:

```
chimclust.pl
```

Example:

```
./chimclust.pl -i Archived/filename-multi.reject.sscan
```

This clusters candidate events by repeat occurrence, helping identify recurrent signals among rejected alignments.

Currently this clustering functionality is implemented for interchromosomal chimeras.

---

# Project Status

This repository preserves the original FusionMiner code primarily for **historical reference and reproducibility**.

The pipeline reflects computational approaches used in early fusion transcript discovery from Sanger sequenced transcript-derived sequence alignments prior to widespread RNA-seq–based fusion detection workflows.

---

# Original Project Location

The original project was hosted on SourceForge:

https://sourceforge.net/projects/fusionminer/

---

# Citation

If referencing FusionMiner or work derived from it, please cite:

**Plebani, R., Oliver, G. R., Trerotola, M., Guerra, E., Cantanelli, P., Apicella, L., Emerson, A., Albiero, A., Harkin, P. D., Kennedy, R. D., & Alberti, S.** (2012).  
*Long-range transcriptome sequencing reveals cancer cell growth regulatory chimeric mRNA.*  
**Neoplasia**, 14(11), 1087–1099.  
https://doi.org/10.1593/neo.121342

\*Authors contributed equally.

---

# License

This repository preserves the original source code for archival and reference purposes.  
Please refer to the original SourceForge project for historical licensing context.
