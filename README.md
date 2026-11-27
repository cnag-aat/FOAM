# FOAM

FOAM (Flye ONT Assembly of Mitogenomes) assembles and finishes complete
circular mitogenomes from long reads, or finishes a mitogenome candidate that
an assembler (flye, hifiasm) already produced. Each sample poses its own
challenges, so FOAM is under continuous development.

## What it does

```
 mode: reads     pass 1: seed ─▶ candidate reads ─▶ select reads ─▶ flye ─▶ candidates ─▶ annotate* (draft)
                 pass 2: own draft as bait ─▶ candidate reads ─▶ select reads ─▶ flye ─▶ candidates ─▶ annotate*
                 ─▶ polish ─▶ annotate* ─▶ final
 mode: flye      flye contig(s) = bait ─▶ candidate reads ─▶ select reads ─▶ annotate* ─▶ polish ─▶ annotate* ─▶ final
 mode: hifiasm   hifiasm contig(s) / whole assembly ─▶ candidates ─▶ annotate* ─▶ final   (no reads, no polishing)

 annotate* = split tandem copies ─▶ MITOS2 on each unique copy ─▶ repair annotation
             ─▶ keep the most complete copy ─▶ rotate/reverse-complement to the start gene
```

### Read fishing

Fishing is built to depend as little as possible on fixed identity or length
cut-offs, or on how close a relative's mitogenome is.

* **Seed (first pass)**, set by `fishing.seed`:
  * `bait`: a related mitogenome, given directly or fetched by taxid. FOAM
    walks up the NCBI taxonomy (species → genus → family → order) and takes
    up to three complete mitogenomes, preferring RefSeq.
  * `genes`: conserved mitochondrial proteins (cox1–3, cob, nad1, nad4, nad5)
    from the MITOS reference set, searched with frameshift-aware
    `diamond blastx`. No related mitogenome is needed, and gene order and
    control-region differences don't matter.
  * `auto` (default): `bait` when a bait or taxid is given, `genes`
    otherwise.
* **Candidate reads.** One streaming pass over all reads (`minimap2` against
  the bait concatenated with itself, or DIAMOND against the seed proteins)
  keeps any read with ≥ 300 bp of hit. This is deliberately permissive.
  Only read IDs and the small candidate set are written; the full read set
  is never rewritten.
* **Selection by k-mer coverage** (`fishing.method: kmer`, the default).
  Within the candidate set, mitogenome k-mers occur hundreds to thousands of
  times; the nuclear flanks of a NUMT read and sequencing errors occur only a
  few times. FOAM:
  1. finds the mitochondrial coverage peak in the k-mer abundance spectrum;
  2. scores each read by the fraction of its k-mers at mitochondrial
     coverage;
  3. keeps reads scoring at least half of a typical mito read.

  Both cut-offs are relative to the data, so they adapt to read accuracy and
  bait divergence, and NUMT-carrying reads are rejected by coverage rather
  than by an identity threshold.

  Counting uses [FastK](https://github.com/thegenemyers/FASTK) (in the
  `mapping` env): every k-mer is counted, and each read's count profile is
  read back with `Profex`, in parallel over read ranges. 300 Mb of
  candidate reads takes about 12 s on 4 threads. Without FastK on PATH (or
  with `fishing.kmer.backend: python`), a pure-Python counter samples
  1/8 of the k-mers instead; it is about 8× slower even so. `results/pass*/read_selection.json` reports
  the spectrum, the peak and the thresholds. `fishing.method: alignment`
  restores the old behaviour: keep reads ≥ 80 % covered by bait alignments.
* **Second pass, self-baited** (reads mode, `fishing.second_pass: true`).
  The pass-1 draft, which is the sample's own mitogenome, becomes the bait
  for a second fishing pass and a fresh flye assembly. Those reads also
  polish the result. With a near-identical bait the prefilter no longer
  depends on how divergent the first seed was.

### Finishing

* **Candidates.** Circular contigs (flye `assembly_info.txt`; hifiasm `…c`
  names) whose length is between 0.75× the expected length and 4 × 1.25× it,
  to allow for tandem copies. A whole hifiasm assembly can be screened
  against the bait. Use `contig_names` to pick contigs by hand.
* **Tandem copies.** Assemblers often output a "circular" contig that is
  2–3 head-to-tail copies. A `minimap2 -DP` self-alignment (k-mer fallback)
  gives the period, and an anchor k-mer cuts all copies at the same
  biological position. Identical copies are annotated once; the copy with
  the most complete annotation wins.
* **Annotation repair.** MITOS drops genes its own searches found. FOAM
  checks the result against the expected gene set and repairs it. Every
  change is listed in `*.report.json` and as a `Note=` in the GFF3:
  * rRNAs MITOS discarded for overlapping neighbours (`ignored.mitos`) are
    restored and trimmed to the neighbouring genes.
  * Missing proteins are rebuilt as the complete ORF around MITOS's own
    BLAST hits (`blast/prot`), preferring ATG/ATA starts.
  * Missing tRNAs are searched for with `cmsearch --max` and the tRNA's
    covariance model. Hits on other tRNAs/rRNAs, or deep inside proteins,
    are rejected.
  * trnL1/L2 and trnS1/S2 are relabelled by their anticodon.
  * If a tRNA covers a protein's 5′ end by ≥ 20 bp, the start moves to the
    next in-frame start codon that clears the tRNA.
  * MITOS's OH/OL predictions are dropped (unreliable outside vertebrates);
    `mitos.keep_origins: true` keeps them.
* **Orientation.** Vertebrates start at trnF; all other taxa start at the
  5′ end of cox1, on the + strand.
* **Polishing** (reads/flye modes): `flye --polish-target` on the selected
  copy, rotated by half its length so the poorly polished ends are not at the
  start gene. The final step re-orients.

## Taxon presets

`--taxon` sets the genetic code, the start gene and the expected gene set.
Any of them can be overridden (`--genetic-code`, `--anchor`, `--expected-genes`).

| preset | code | start | genes |
|---|---|---|---|
| vertebrata | 2 | trnF | 37 |
| mollusca, arthropoda, annelida, nematoda, invertebrate | 5 | cox1 | 37 |
| echinodermata, platyhelminthes | 9 | cox1 | 37 |
| tunicata | 13 | cox1 | 37 |
| cnidaria, porifera | 4 | cox1 | 13 PCG + 2 rRNA + trnM, trnW |

With `--taxid` and `--email`, `create_config.py` picks the preset from the NCBI
lineage and takes the mitochondrial genetic code from NCBI Taxonomy.

## Running it

FOAM needs Snakemake ≥ 9. Snakemake builds the tool environments itself
(`--use-conda`, from `envs/`): minimap2, seqkit, DIAMOND and FastK
(`mapping`), flye (`flye`), and MITOS 2.1.10 with Infernal (`mitos`, Linux
only). FOAM's own Python code uses only the standard library.

1. Clone this repository and create a driver environment for Snakemake:

   ```bash
   conda create -n snakemake -c conda-forge -c bioconda \
       "snakemake>=9" snakemake-executor-plugin-slurm
   ```
   Download the MITOS refseq81 databases from
   https://zenodo.org/records/2672835 into `utils/refseq_dir` and unpack them
   (`tar jxf refseq81m.tar.bz2`).
2. Write a config:

   ```bash
   # a hifiasm candidate (no reads needed)
   ./create_config.py --mode hifiasm --tolid xgXerMont --taxon mollusca \
       --contigs xgXerMont_mitogenome_candidate.fa

   # a whole hifiasm assembly, screened with a bait fetched from NCBI
   ./create_config.py --mode hifiasm --tolid xgXerMont --taxid 1234 --email you@inst.org \
       --contigs asm.p_ctg.fa --bait related_mt.fa

   # flye circular contig(s) + reads (polished)
   ./create_config.py --mode flye --tolid ibLobCana --taxon arthropoda \
       --contigs flye/assembly.fasta --flye-info flye/assembly_info.txt \
       --reads reads.ont.fastq.gz --read-type ont-hq

   # from reads only, bait fetched from NCBI
   ./create_config.py --mode reads --tolid ibLobCana --taxid 6973 --email you@inst.org \
       --reads reads.ont.fastq.gz --read-type ont-hq --expected-length 15025

   # from reads only, no bait at all (conserved-protein seed)
   ./create_config.py --mode reads --tolid ibLobCana --taxon arthropoda \
       --reads reads.ont.fastq.gz --read-type ont-hq
   ```
   `mitos.runmitos: auto` (the default) runs `runmitos` (MITOS ≥ 2.1.7) or
   `runmitos.py` (older), whichever is on PATH; `--runmitos` in
   `create_config.py` sets another command.
3. Build the conda environments once, on a node with internet access
   (compute nodes often have none). `--conda-prefix` shares one copy
   between runs:

   ```bash
   snakemake -s /path/to/FOAM/bin/foam.smk --configfile xgXerMont.foam.yaml \
       --use-conda --conda-prefix ~/foam-conda-envs --conda-create-envs-only
   ```
4. Run, locally or on SLURM with one of the example profiles:

   ```bash
   # locally / inside one job
   snakemake -s /path/to/FOAM/bin/foam.smk --configfile xgXerMont.foam.yaml \
       --use-conda --conda-prefix ~/foam-conda-envs -j 16

   # SLURM (CNAG partitions all / research / research_long)
   snakemake -s /path/to/FOAM/bin/foam.smk --configfile xgXerMont.foam.yaml \
       --conda-prefix ~/foam-conda-envs --profile /path/to/FOAM/examples/profiles/sm9-slurm-cnag
   ```

   | profile | partitions |
   |---|---|
   | `sm9-slurm-cnag` | each rule on `all` (≤ 1 h), `research` (≤ 12 h) or `research_long` with a fitting runtime |
   | `sm9-slurm-sp_gpus` | every job on `sp_gpus`, no time limit (no GPUs requested) |

   Both set memory per rule and the account `sa_assembly`; override any of
   it with `--set-resources rule:key=value`. In a profile, partition names
   are written `"'all'"`: Snakemake evaluates resource values as Python, and
   a bare `all` is Python's `all()` function.

Outputs go to `results/final/<tolid>_MT.{fa,gff3,bed,genes.tsv,report.json}`;
each step logs to `logs/<rule>.log`. Intermediate files are deleted unless
`keep_intermediate: true`.

The steps can also be run by hand, e.g. to (re)annotate a candidate:

```bash
bin/foam.py annotate -i candidate.fa -o out/sp_MT --taxon mollusca \
    --refdir utils/refseq_dir --refseq refseq81m
```

`bin/foam.py --help` lists all subcommands (`annotate`, `split-copies`,
`candidates`, `fish-filter`, `kmer-select`, `seed-proteins`, `double-bait`,
`rotate`, `fetch-bait`).

## Tests

```bash
pip install pytest pyyaml        # snakemake, minimap2 and FastK on PATH enable more tests
python -m pytest tests
```

The annotate tests use a fake `runmitos` (`tests/data/fake_runmitos.py`), so
MITOS is not needed. `foam/` is standard-library-only Python and stays
compatible with Python 3.6, so it also runs in an older MITOS 2.1.3 environment.

## Possible further improvements

1. **Assembly-graph check.** Use the flye/hifiasm graph (GFA) to confirm
   circularity and to detect real heteroplasmic/duplicated regions before
   collapsing tandem copies.
2. **Polish on a circular target.** Polish a copy extended by a few kb on each
   side (or a 2-copy dimer) and cut the polished monomer back out with the
   tandem splitter, so the linearisation point gets full read support.
3. **Second annotator.** Compare with MitoZ/MitoFinder or a tRNAscan-SE/ARWEN
   tRNA pass, and report disagreements, for divergent taxa where MITOS2's
   reference set is thin.
## Contributors

Fernando Cruz, Jèssica Gómez-Garrido and Tyler Alioto. Genome Assembly and
Annotation Team (CNAG). https://denovo.cnag.cat/about

## References

A chromosome-level reference genome for the common octopus, *Octopus vulgaris*
(Cuvier, 1797). Destanović D, Schultz DT, Styfhals R, Cruz F, Gómez-Garrido J,
Gut M, Gut I, Fiorito G, Simakov O, Alioto TS, Ponte G, Seuntjens E.
https://doi.org/10.1101/2023.05.16.540928
