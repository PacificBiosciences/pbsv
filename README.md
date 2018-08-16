<h1 align="center"><img width="300px" src="img/sv2.png"/></h1>
<h1 align="center">pbsv</h1>
<p align="center">PacBio structural variant (SV) calling and analysis tools</p>

***

`pbsv` is a suite of tools to call and analyze structural variants
in diploid genomes from PacBio single molecule real-time sequencing (SMRT) reads.
The tools power the *Structural Variant Calling* analysis workflow in
PacBio's SMRT Link GUI.

`pbsv` calls insertions, deletions, inversions, and translocations.
Both single-sample calling and joint (multi-sample) calling are provided. `pbsv` is most
effective for:
* insertions 20 bp to 5 kb
* deletions 20 bp to 100 kb
* inversions 200 bp to 5 kb
* translocations between different chromosomes or further than 100kb apart on a single chromosome

## Availability
Stable versions of `pbsv` with official support from PacBio are only
available through [SMRT Link releases](https://www.pacb.com/support/software-downloads/).
* SMRT Link 5.1.0 includes `pbsv` 1.0.0

The current pre-release version of the `pbsv` binary is available for Linux from
[bioconda](https://bioconda.github.io/).  All necessary dependencies are installed
automatically.  A modern (post-2008) CPU with support for [SSE4.1 instructions](https://en.wikipedia.org/wiki/SSE4#SSE4.1) is required.

    conda install pbsv

## Workflow
<p align="center"><img width="700px" src="img/pbsv-stage-workflow.png"/></p>

The general `pbsv` workflow is:
1. Align PacBio reads to a reference genome, per movie. (`.subreads.bam` to `.bam`)
2. Discover signatures of structural variation, per movie or per sample. (`.bam` to `.svsig.gz`)
3. Call structural variants and assign genotypes, all samples. (`.svsig.gz` to `.vcf`)

### 1. Align PacBio reads to a reference genome
For each movie (`.subreads.bam`), extract subreads with `pbsv fasta` and align
to a reference genome (`ref.fa`) with [`minimap2`](https://github.com/lh3/minimap2) or
[`ngmlr`](https://github.com/philres/ngmlr).

```sh
pbsv fasta movie1.subreads.bam |
    minimap2 -x map-pb -a --eqx -L -O 5,56 -E 4,1 -B 5 --secondary=no -z 400,50 -r 2k -Y -R "@RG\tID:rg1a\tSM:sample1" ref.fa - |
    samtools sort > ref.movie1.bam
samtools index ref.movie1.bam
```

It is recommended to use multithreading (`-t` for `minimap2`, `-@` for `samtools`).

The recommended `minimap2` options are:
* `-x map-pb` provides preset seeding parameters optimized for PacBio reads
* `-a` produces SAM output (**required** for `pbsv`)
* `--eqx` uses `X`/`=` extended CIGAR strings
* `-L` supports long alignments
* `-O 5,56 -E 4,1 -B 5` approximates the convex gap costs of `ngmlr`
* `--secondary=no` outputs only primary and supplementary alignments
* `-z 400,50` enables alignment of short inversions
* `-r 2k` increases alignment bandwidth to span large insertions and deletions
* `-Y` prevents hard clipping (**required** for `pbsv`)
* `-R` defines a read group with the sample name (**required** for joint calling with `pbsv`)

The sample name, stored in the `SM` tag of the read groups, associates
aligned reads with a particular sample.  It it required for downstream
joint calling.

### 2. Discover signatures of structural variation

For each aligned BAM or set of aligned BAMs for a single sample, identify signatures
of structural variation.  This reduces all aligned reads to those that are relevant
to calling structural variants.  The signatures are stored in a `.svsig.gz` file.

```sh
pbsv discover ref.movie1.bam ref.sample1.svsig.gz
pbsv discover ref.movie2.bam ref.sample2.svsig.gz
```

Sample names are transferred from the `RG` headers to the `.svsig.gz` file.

### 3. Call structural variants and assign genotypes

Call structural variants from structural variant signatures, jointly for all
samples of interest. One or more `.svsig.gz` files are accepted, including multiple
`.svsig.gz` for a single sample and/or `svsig.gz` for multiple samples.

```sh
pbsv call ref.fa ref.sample1.svsig.gz ref.sample2.svsig.gz ref.var.vcf
```

Variant calls for all samples are output in a single `.vcf` file.

### Parallel processing per chromosome
For large genomes with high sequencing coverage, it is recommended to process
chromosomes separately.  After aligning each movie:

#### Generate separate `.svsig.gz` files per chromosome

```sh
for i in {chr1,chr2,chr3,chr4,chr5,...}; do
    pbsv discover --region $i hg38.movie1.bam hg38.sample1.$i.svsig.gz
done
```

#### Call insertions, deletions, and inversion per chromosome
```sh
for i in {chr1,chr2,chr3,chr4,chr5,...}; do
    pbsv call --types INS,DEL,INV hg38.fa hg38.sample1.${i}.svsig.gz hg38.${i}.ins+del+inv.vcf
done
```

#### Call translocations for whole genome
Translocation (breakend) calling requires all chromosomes.  It is not currently supported
to call translocations separately per chromosome.

```sh
pbsv call --types BND hg38.fa hg38.sample1.*.svsig.gz hg38.*.bnd.vcf
```

## Algorithm Overview and Advanced Parameters
### Deletions
<p align="center"><img width="700px" src="img/pbsv-deletion-workflow.png"/></p>

Cluster options used during `pbsv call`:

```
SV Signature Cluster Options:
  --cluster-max-length-perc-diff   Do not cluster signatures with difference in length > P%. [25]
  --cluster-max-ref-pos-diff       Do not cluster signatures > N bp apart in reference. [200]
```

Number of flanks used for consensus generation:
```
Consensus Options:
  -x,--max-consensus-coverage      Limit to N reads for variant consensus. [20]
```

**Split deletions**: Deletions that are not fully aligned using the `D` cigar are recovered up to a
size of 100kb. Deletions greater than 100kb are currently called as translocations.

### Insertions
Insertion calling workflow is identical to the above described deletion workflow,
except for one additional criteria, the inserted sequence similarity check
during clustering:

```
SV Signature Cluster Options:
  --cluster-min-basepair-perc-id   Do not cluster signatures with basepair identity < P%. [10]
```

### Inversions
<p align="center"><img width="800px" src="img/pbsv-insertion-criteria.png"/></p>

An inversion signature is detected if a single read is split into three
alignments with different orientations / strands, either `+-+` or `-+-`.
The maximum permitted reference gap or overlap between consecutive alignments is
configured in `pbsv call`:

```
Alignment Connection Options:
 --max-inversion-gap   Do not link inverted alignments with > N bp gap or overlap with flanking alignments. [1000]
```

Clustering is performed on the inverted segment and uses the same criteria as deletion clustering.

The VCF call marks the most likely position and size of the inverted segment, as shown in this IGV screenshot:
<p align="center"><img width="800px" src="img/pbsv-inversion-igv.png"/></p>

### Translocations
Translocations are identified using breakends of individual reads.  All four breakend combinations are supported:
<p align="center"><img width="700px" src="img/pbsv-breakends.png"/></p>

### Calling and Genotyping
An variant is output if it passes all of the following criteria:
* supported by at least `--call-min-reads-all-samples [2]` reads total across samples,
* supported by at least `--call-min-reads-one-samples [2]` in a sample,
* supported by at least `--call-min-read-per-one-sample [20]` percent of reads in a sample
* assigned a non-reference genotype in at least one sample
** A sample is assigned a non-reference genotype for a variant if at least `--gt-min-reads [1]` reads
   support the variant.

### Filtering
The VCF filter column is

1) **PASS**
2) **NearReferenceGap**: variant is near (< `--filter-near-reference-gap [1000]`) from a gap (run of >= 50 Ns in the reference assembly)
3) **Decoy**: variant involves a decoy sequence, where the chromosome name contains `decoy`, `hs37d5`, or `hs38d1`

## Known Issues
`pbsv` is under active development and will continue to improve in future release.  Currently known issues and limitations are:

* Sensitivity is limited for insertions > 5kb due to difficulty aligning with current minimap2 parameters.
* Large deletions are sometimes called multiple times.
* The REF allele sequence is shifted by one base for insertions and deletions and does not match the reported genomic position.
* The REF allele representation for deletions is one base too short.  The reported POS and SVLEN are correct.
* The output VCF is not properly sorted in edge cases.

## FAQ

### To where do I report bugs and ask questions about the pre-release version of `pbsv`?
Report bugs and questions using GitHub Issues.  The pre-release version of `pbsv` is not
officially supported, but feedback from users is appreciated and will be addressed as possible.

## Change Log
 * 2.0.0: Drop RC for conda release
 * 2.0.0-RC2: First public release candidate for SMRT Link 6.0.0


## DISCLAIMER

THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.

FOR RESEARCH USE ONLY. NOT FOR USE IN DIAGNOSTICS PROCEDURES.
