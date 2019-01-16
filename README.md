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
Latest version can be installed via bioconda package `pbsv`.

Please refer to our [official pbbioconda page](https://github.com/PacificBiosciences/pbbioconda)
for information on Installation, Support, License, Copyright, and Disclaimer.

## Latest Version
Version **2.1.1**: [Full changelog here](#full-changelog)

## Workflow
<p align="center"><img width="700px" src="img/pbsv-stage-workflow.png"/></p>

The general `pbsv` workflow is:
1. Align PacBio reads to a reference genome, per movie. (`.subreads.bam`/`.ccs.fastq` to `.bam`)
2. Discover signatures of structural variation. (`.bam` to `.svsig.gz`)
3. Call structural variants and assign genotypes, all samples. (`.svsig.gz` to `.vcf`)

### 1. Align PacBio reads to a reference genome
The recommended aligner is `pbmm2` that can be installed via `conda install pbmm2`.
For each movie (`.subreads.bam` or `.ccs.fq`) align records to a
reference genome (`ref.fa`).

```sh
#For .subreads.bam input:
pbmm2 align ref.fa movie1.subreads.bam ref.movie1.bam --sort --sample 'sample1' --median-filter

#For .ccs.fq input:
pbmm2 align ref.fa movie1.ccs.fq ref.movie1.bam --sort --rg '@RG\tID:movie1\tSM:sample1' --preset CCS
```

The sample name, stored in the `SM` tag of the read groups, associates
aligned reads with a particular sample.  It is required for downstream
joint calling.

### 2. Discover signatures of structural variation

For each aligned BAM or set of aligned BAMs, identify signatures
of structural variation. This reduces all aligned reads to those that are relevant
to calling structural variants.  The signatures are stored in a `.svsig.gz` file.

```sh
pbsv discover ref.movie1.bam ref.sample1.svsig.gz
pbsv discover ref.movie2.bam ref.sample2.svsig.gz
```

It is highly recommended to provide one tandem repeat annotation `.bed` file
of your reference to `pbsv discover` via `--tandem-repeats`. This increases
sensitivity and recall.

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
for i in $(samtools view -H hg38.movie1.bam | grep '^@SQ' | cut -f2 | cut -d':' -f2); do
    pbsv discover --region $i hg38.movie1.bam hg38.sample1.$i.svsig.gz
done
```

#### Call SVs per chromosome
```sh
for i in {chr1,chr2,chr3,chr4,chr5,...}; do
    pbsv call hg38.fa hg38.sample1.${i}.svsig.gz hg38.${i}.vcf
done
```
Be aware that each translocation will get called twice, when run per chromosome.
Even though the IDs will match, the `DP` information won't necessarily be identical.

To avoid that, call insertions, deletions, and inversions independent from
translocations:

```sh
for i in {chr1,chr2,chr3,chr4,chr5,...}; do
    pbsv call --types INS,DEL,INV hg38.fa hg38.sample1.${i}.svsig.gz hg38.${i}.ins+del+inv.vcf
done
pbsv call --types BND hg38.fa hg38.sample1.*.svsig.gz hg38.bnd.vcf
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
 --max-inversion-gap   Do not link inverted alignments with > N bp gap or overlap with flanking alignments. [1K]
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
* supported by at least `--call-min-read-per-one-sample [20]` percent of reads in a sample,
* supported by at least `--call-min-reads-per-strand-all-samples [1]` reads per strand total across samples,
* assigned a non-reference genotype in at least one sample;
  a sample is assigned a non-reference genotype for a variant if at least `--gt-min-reads [1]` reads
  support the variant.

### Filtering
The VCF filter column is

1) **PASS**
2) **NearReferenceGap**: variant is near (< `--filter-near-reference-gap [1K]`) from a gap (run of >= 50 Ns in the reference assembly)
3) **Decoy**: variant involves a decoy sequence, where the chromosome name contains `decoy`, `hs37d5`, or `hs38d1`
4) **NearContigEnd**: variant is near (< `--filter-near-contig-end [1K]`) from a contig end
5) **InsufficientStrandEvidence**: variant is not supported by at least (`--call-min-reads-per-strand-all-samples [1]`) reads in forward and reverse orientation

## Performance benchmarks
Using the [publicly available HG002 15kb CCS dataset](https://bit.ly/2RW1b3I),
we are tracking `pbsv` performance with respect to the genome in a bottle annotation version 0.6.

### Step by step
#### Prepare and gather data
1) Create directory structure:
```sh
mkdir fastqs ref alns svsigs giab
```

2) Download genome in a bottle annotations:
```sh
FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/
curl -s ${FTPDIR}/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed > giab/HG002_SVs_Tier1_v0.6.bed
curl -s ${FTPDIR}/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz > giab/HG002_SVs_Tier1_v0.6.vcf.gz
```

3) Make sure that you have exact truvari version (check that all python dependencies are met, not described):
```sh
git clone https://github.com/spiralgenetics/truvari
(cd truvari; git reset --hard 600b4ed7)
```

4) Download hg19 reference with decoys and map non-ACGT characters to N:
```sh
curl -s ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz > ref/human_hs37d5.fasta.gz
gunzip ref/human_hs37d5.fasta.gz
sed -i '/^[^>]/ y/BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvwxyz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' human_hs37d5.fasta
```

5) Download hg19 tandem repeat annotations:
```sh
curl -s https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_hs37d5.trf.bed > ref/human_hs37d5.trf.bed
```

6) Download all `.fastq` files:
```sh
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180901_011437.Q20.fastq > fastqs/m54238_180901_011437.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180902_013549.Q20.fastq > fastqs/m54238_180902_013549.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180903_015530.Q20.fastq > fastqs/m54238_180903_015530.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180904_021549.Q20.fastq > fastqs/m54238_180904_021549.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180907_170406.Q20.fastq > fastqs/m54238_180907_170406.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180908_172515.Q20.fastq > fastqs/m54238_180908_172515.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180909_174539.Q20.fastq > fastqs/m54238_180909_174539.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180910_180559.Q20.fastq > fastqs/m54238_180910_180559.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180913_181445.Q20.fastq > fastqs/m54238_180913_181445.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180914_183539.Q20.fastq > fastqs/m54238_180914_183539.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180915_185611.Q20.fastq > fastqs/m54238_180915_185611.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180916_191625.Q20.fastq > fastqs/m54238_180916_191625.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180919_165326.Q20.fastq > fastqs/m54238_180919_165326.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180920_171425.Q20.fastq > fastqs/m54238_180920_171425.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180921_173448.Q20.fastq > fastqs/m54238_180921_173448.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180922_175520.Q20.fastq > fastqs/m54238_180922_175520.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180925_225123.Q20.fastq > fastqs/m54238_180925_225123.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180926_231301.Q20.fastq > fastqs/m54238_180926_231301.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180927_233325.Q20.fastq > fastqs/m54238_180927_233325.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54328_180921_232856.Q20.fastq > fastqs/m54328_180921_232856.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54328_180922_235017.Q20.fastq > fastqs/m54328_180922_235017.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54328_180924_001027.Q20.fastq > fastqs/m54328_180924_001027.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54328_180925_003051.Q20.fastq > fastqs/m54328_180925_003051.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54328_180926_222309.Q20.fastq > fastqs/m54328_180926_222309.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54328_180927_224427.Q20.fastq > fastqs/m54328_180927_224427.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54328_180928_230446.Q20.fastq > fastqs/m54328_180928_230446.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54328_180929_232524.Q20.fastq > fastqs/m54328_180929_232524.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54329_180924_222717.Q20.fastq > fastqs/m54329_180924_222717.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54329_180925_224838.Q20.fastq > fastqs/m54329_180925_224838.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54329_180926_230856.Q20.fastq > fastqs/m54329_180926_230856.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54329_180927_232921.Q20.fastq > fastqs/m54329_180927_232921.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54334_180924_221206.Q20.fastq > fastqs/m54334_180924_221206.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54334_180925_223328.Q20.fastq > fastqs/m54334_180925_223328.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54334_180926_225337.Q20.fastq > fastqs/m54334_180926_225337.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54334_180927_231334.Q20.fastq > fastqs/m54334_180927_231334.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54335_180924_221150.Q20.fastq > fastqs/m54335_180924_221150.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54335_180925_223313.Q20.fastq > fastqs/m54335_180925_223313.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54335_180926_225328.Q20.fastq > fastqs/m54335_180926_225328.Q20.fastq
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54335_180927_231344.Q20.fastq > fastqs/m54335_180927_231344.Q20.fastq
```

6) Index reference
```sh
pbmm2 index ref/human_hs37d5.fasta ref/human_hs37d5.mmi --preset CCS
```

### Process data
7) Align each movie, can be done in parallel
```sh
for i in fastqs/*.fastq; do
    FILENAME="${i#*fastqs/}"
    FILEPREFIX="${FILENAME%.*}"
    pbmm2 align ref/human_hs37d5.mmi $i "alns/hg19.${FILEPREFIX}.bam" --preset CCS \
                --sort --rg '@RG\tID:${FILEPREFIX}' --sample HG2
done
```

8) Discover SV signatures for each alignment, can be done in parallel
```sh
for i in alns/*.bam; do
    FILENAME="${i#*alns/}"
    FILEPREFIX="${FILENAME%.*}"
    pbsv discover $i "svsigs/hg19.${FILEPREFIX}.svsig.gz" --tandem-repeats ref/human_hs37d5.trf.bed
done
```

9) Call and polish SVs
```sh
pbsv call ref/human_hs37d5.fasta svsigs/*.svsig.gz hg2.pbsv.vcf --log-level INFO -z 1G\
          --call-min-read-perc-one-sample 10
bgzip hg2.pbsv.vcf
tabix hg2.pbsv.vcf.gz
```

10) Compare to ground truth
```sh
truvari/truvari.py -f ref/human_hs37d5.fasta -b  giab/HG002_SVs_Tier1_v0.6.vcf.gz\
                   --includebed giab/HG002_SVs_Tier1_v0.6.bed -o testrun --passonly\
                   --giabreport -r 1000 -p 0.01 -c hg2.pbsv.vcf.gz
```

11) Parse results
```sh
function sumsv() { cat $1 | grep ':' | tr -d ',' |sed "s/^[ \t]*//"| tr -d '"' |\
                   tr -d ' ' | tr ':' '\t' | awk '{ for (i=1; i<=NF; i++)  {
                     a[NR,i] = $i } } NF>p { p = NF } END { for(j=1; j<=p; j++)
                     { str=a[1,j]; for(i=2; i<=NR; i++){ str=str" "a[i,j]; } print str } }' |\
                   tail -n 1 | awk '{ printf "%1.4f\t%1.4f\t%1.4f\t%10.0f\n", $2,$4,$11,$1+$8 }';}
cat <(echo -e "Run\tF1\tPrecision\tRecall\tFP+FN")\
    <(for i in testrun; do printf $i"\t";sumsv $i/summary.txt;done) |\
    sed 's/testrun\///g;' | sort -k 2 -n | column -t
```

## Known Issues
`pbsv` is under active development and will continue to improve in future release.  Currently known issues and limitations are:

* Some LINE elements are not detected due to difficulties in alignment. This will improve in upcoming versions of `pbmm2`.

## FAQ

### To where do I report bugs and ask questions about the pre-release version of `pbsv`?
Please refer to our [official pbbioconda page](https://github.com/PacificBiosciences/pbbioconda)
to report bugs and ask questions.

### Where can I find an example dataset to try `pbsv`?
10-fold coverage of the Genome in a Bottle sample HG002 is [available](https://downloads.pacbcloud.com/public/dataset/HG002/Sequel-201810/).

### The binary does not work on my linux system!
If you get `Illegal instruction` upon execution of `pbsv`, then your CPU is not supported.
A modern (post-2008) CPU with support for [SSE4.1 instructions](https://en.wikipedia.org/wiki/SSE4#SSE4.1) is required.

## Full Changelog
 * **2.1.1**:
   * Improve error output if reference and svsig contigs do not match
   * Add `IMPRECISE` to VCF header
   * Separate filter with semicolon

 * 2.1.0:
   * Algorithmic improvements to increase recall and sensitivity across all SV lengths
   * Add `SAC` Stranded Allel Counts for subread input
   * Remove `pbsv fasta` and rely on `pbmm2` input
   * Add `MATEDIST`,  distance between two breakends of a translocation if the two breakends are on the same contig
   * Fix VCF POS/REF/ALT and sorting to pass GATK `ValidateVariants` and allow VCF `bgzip` and `tabix`
   * Allow multiple sample names and multiple BAM inputs for `pbsv discover`
   * Allow length `KMG` suffix
   * Flag **NearContigEnd**
   * Allow breakend calling from single chromosome `svsig.gz`
   * Allow custom `--annotations` for known sequence types
   * Loop termination fix to properly fix the issue in the 2.0.2 patch
   * Decrease memory overhead in `pbsv discover`

 * 2.0.2: Fix rare `pbsv call` abort, because of missing coverage
 * 2.0.0: Drop RC for conda release
 * 2.0.0-RC2: First public release candidate for SMRT Link 6.0.0


## DISCLAIMER

THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.

FOR RESEARCH USE ONLY. NOT FOR USE IN DIAGNOSTICS PROCEDURES.
