<h1 align="center"><img width="300px" src="img/sv2.png"/></h1>
<h1 align="center">pbsv2</h1>
<p align="center">SMRT Structural Variation Caller</p>

***

## Scope
*pbsv2* contains the newest tools to identify structural variations in
PacBio single-molecule sequencing data.
Starting in SMRT Link v6.0.0, those tools power the
*Structural Variation GUI-based analysis* application.

## Availability
The latest pre-release, developers-only linux binaries can be found under
[releases](https://github.com/PacificBiosciences/pbsv2/releases).
These binaries are not ISO compliant.
For research only.
Not for use in diagnostics procedures.

Official support is only provided for official and stable SMRT Link builds
provided by PacBio.

Unofficial support for binary pre-releases is provided via github issues,
not via mail to developers.

Binaries require **SSE4.1 CPU support**; CPUs after 2008 (Penryn) include it.

## Algorithm

ToDo: Explain
* core algo incl realignment
* the idea behind the `sv` file
* how joint calling is a native feature

## Workflow
`pbsv2` accepts aligned reads in the BAM format as input.
The general workflow looks as following, align each movie, convert each `.bam`
file to an intermediate sparse representation `.sv.gz` using `pbsv2 parse`,
and then jointly call all `*.sv.gz` together with `pbsv2 call`.

### Convert movie to fasta
For each movie, we are extracting the subreads from unaligned `.bam` files to
`.fasta` file. The trick employed here, only take ONE subread per ZMW, possibly
one that spans the full-length molecule:

    pbsv2 filter movie1.subreads.bam > movie1.fasta

### Alignment
For each movie, align them with your favorite aligner,
we recommend [minimap2](https://github.com/lh3/minimap2)
and [ngmlr](). Example call for `minimap`:

    minimap2 -a -O 5,56 -E 4,1 --secondary=no -z 400,50 -r2k -Y -R "@RG\\tID:sv" hg38.fasta movie1.fasta |\
    samtools view -bS > movie1.aligned.bam

Options explained:
* `-a` for SAM output
* `-O 5,56 -E 4,1` approximates the convex gap costs of `ngmlr`
* `--secondary=no` Only primary hits are of interest
* `-r2k` increase bandwidth to not break indels
* `-Y` do not hard clip
* `-R` add a read group, required for PacBio BAM files

### Parse
Each aligned `.bam` file has to be converted to a `.sv.gz` file. For this,
`pbsv2 parse` requires three arguments:

    pbsv2 parse movie1.aligned.bam HG0733 movie1.sv.gz

The first argument is the aligned `.bam` input file, followed by the bio sample name,
and third is the output `.sv.gz` file name. Multiple `.sv.gz` can have the share
a bio sample name.

### Call
After all `.sv.gz` have been created, they can be jointly called, together
with the used reference in `.fasta` format. For this,
`pbsv2 call` requires at least three arguments, different ways to call it:

    pbsv2 call hg38.fasta movie1.sv.gz mysample.vcf
    pbsv2 call hg38.fasta movie1.sv.gz movie2.sv.gz movie3.sv.gz mysample.vcf
    pbsv2 call hg38.fasta *.sv.gz mysample.vcf

### Speed ups
In order to minimize IO, stream `pbsv2 filter` into your aligner:

    minimap2 ADDITIONAL_OPTIONS hg38.fasta <(pbsv2 filter movie1.subreads.bam) | samtools view -bS > movie1.aligned.bam

## DISCLAIMER

THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
