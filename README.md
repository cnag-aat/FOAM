# FOAM

FOAM (Flye ONT Assembly of Mitogenomes) assembles and finishes complete
circular mitogenomes from long reads, or finishes a mitogenome candidate that
an assembler (flye, hifiasm) already produced. Each sample poses its own
challenges, so FOAM is under continuous development.

## What it does

```
 mode: reads     bait ─▶ fish reads ─▶ flye ─▶ candidates ─▶ annotate* ─▶ polish ─▶ annotate* ─▶ final
 mode: flye      flye contig(s) ─▶ candidates ─▶ fish reads (contigs = bait) ─▶ annotate* ─▶ polish ─▶ annotate* ─▶ final
 mode: hifiasm   hifiasm contig(s) / whole assembly ─▶ candidates ─▶ annotate* ─▶ final   (no reads, no polishing)

 annotate* = split tandem copies ─▶ MITOS2 on each unique copy ─▶ repair annotation
             ─▶ keep the most complete copy ─▶ rotate/reverse-complement to the start gene
```

* **Bait.** Give a related mitogenome, or a taxid. With a taxid, FOAM walks up
  the NCBI taxonomy (species → genus → family → order) and takes up to three
  complete mitogenomes, preferring RefSeq. This is the former
  `get_closest_mitogenome.py`, now `bin/foam.py fetch-bait`.
* **Read fishing.** Reads are streamed through `minimap2` against the bait
  concatenated with itself (so reads crossing the bait's linearisation point
  align end-to-end). A read is kept only if its alignments cover ≥ 80 % of
  it, so nuclear reads that merely contain a NUMT are dropped. Only the read
  IDs and the selected reads are written: the full read set is never
  rewritten.
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

1. Clone this repository. Download the MITOS refseq81 databases from
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

   # from reads only
   ./create_config.py --mode reads --tolid ibLobCana --taxid 6973 --email you@inst.org \
       --reads reads.ont.fastq.gz --read-type ont-hq --expected-length 15025
   ```
   Set `mitos.runmitos` to `runmitos.py` for MITOS 2.1.3 or `runmitos` for
   ≥ 2.1.7 (`--runmitos` in `create_config.py`).
3. Run:

   ```bash
   snakemake -s /path/to/FOAM/bin/foam.smk --configfile xgXerMont.foam.yaml --use-conda -j 16
   # SLURM: --profile /path/to/FOAM/examples/profiles/slurm
   ```

Outputs go to `results/final/<tolid>_MT.{fa,gff3,bed,genes.tsv,report.json}`;
each step logs to `logs/<rule>.log`. Intermediate files are deleted unless
`keep_intermediate: true`.

The steps can also be run by hand, e.g. to (re)annotate a candidate:

```bash
bin/foam.py annotate -i candidate.fa -o out/sp_MT --taxon mollusca \
    --refdir utils/refseq_dir --refseq refseq81m --runmitos runmitos.py
```

`bin/foam.py --help` lists all subcommands (`annotate`, `split-copies`,
`candidates`, `fish-filter`, `double-bait`, `rotate`, `fetch-bait`).

## Tests

```bash
pip install pytest pyyaml        # snakemake and minimap2 on PATH enable more tests
python -m pytest tests
```

The annotate tests use a fake `runmitos` (`tests/data/fake_runmitos.py`), so
MITOS is not needed. `foam/` is standard-library-only Python and stays
compatible with Python 3.6, the version in the MITOS 2.1.3 environment.

## Ideas for making it more robust (not implemented yet)

Read fishing still depends on cut-offs (`min_read_frac`, `min_identity`,
`min_aln`) and on how close the bait is. Possible next steps:

1. **Two-pass, self-baited fishing.** Fish loosely with the related bait,
   assemble, then re-fish with the sample's own draft mitogenome. The second
   pass is near-identical, so thresholds hardly matter and NUMTs separate
   cleanly by identity.
2. **k-mer based selection.** Count k-mers in the reads (or a subsample).
   Mito k-mers have 100–1000× nuclear coverage, so picking reads rich in
   high-abundance k-mers that also occur in the draft/bait needs no
   alignment-identity thresholds. It is also faster than minimap2 on the full
   read set, and it can be the bait-free fallback when no relative is
   sequenced.
3. **Gene-based seeds.** Instead of a whole mitogenome from another genus,
   seed with conserved protein-coding genes (cox1, cob, nad5) and fish by
   translated alignment (e.g. miniprot/diamond on reads). Order and
   control-region differences between taxa then stop mattering.
4. **Coverage-based NUMT rejection.** Mito reads have a coverage peak far
   above nuclear reads. Cluster candidate reads by depth of their best
   contig and drop the low-depth cluster instead of using identity cut-offs.
5. **Assembly-graph check.** Use the flye/hifiasm graph (GFA) to confirm
   circularity and to detect real heteroplasmic/duplicated regions before
   collapsing tandem copies.
6. **Polish on a circular target.** Polish a copy extended by a few kb on each
   side (or a 2-copy dimer) and cut the polished monomer back out with the
   tandem splitter, so the linearisation point gets full read support.
7. **Second annotator.** Compare with MitoZ/MitoFinder or a tRNAscan-SE/ARWEN
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
