import os
import subprocess
import sys

import pytest

from conftest import needs_snakemake

yaml = pytest.importorskip("yaml")


def make_config(root, tmp_path, *args):
    out = str(tmp_path / "cfg.yaml")
    subprocess.check_call([sys.executable, os.path.join(root, "create_config.py"),
                           "-o", out] + list(args))
    return out, yaml.safe_load(open(out))


@pytest.fixture
def files(tmp_path):
    d = {}
    for n in ("reads.fq.gz", "bait.fa", "contigs.fa", "assembly_info.txt"):
        p = tmp_path / n
        p.write_text(">x\nACGT\n" if n.endswith(".fa") else "")
        d[n] = str(p)
    return d


def test_presets_in_config(root, tmp_path, files):
    _, c = make_config(root, tmp_path, "--mode", "hifiasm", "--tolid", "xgTest",
                       "--taxon", "mollusca", "--contigs", files["contigs.fa"])
    assert c["genetic_code"] == 5 and c["anchor"] == ["cox1"] and c["mode"] == "hifiasm"
    _, c = make_config(root, tmp_path, "--mode", "hifiasm", "--tolid", "mTest",
                       "--taxon", "vertebrata", "--contigs", files["contigs.fa"])
    assert c["genetic_code"] == 2 and c["anchor"][0] == "trnF"
    _, c = make_config(root, tmp_path, "--mode", "flye", "--tolid", "iTest", "--taxon",
                       "arthropoda", "--contigs", files["contigs.fa"], "--reads",
                       files["reads.fq.gz"], "--read-type", "hifi", "--genetic-code", "9")
    assert c["genetic_code"] == 9
    assert c["flye"]["read_flag"] == "--pacbio-hifi"
    assert c["fishing"]["minimap2_preset"] == "map-hifi"


@pytest.mark.parametrize("args", [
    ["--mode", "reads", "--tolid", "x", "--taxon", "mollusca"],                 # no reads
    ["--mode", "reads", "--tolid", "x", "--taxon", "mollusca", "--reads", "CONTIGS",
     "--seed", "bait"],                                                         # seed bait, no bait
    ["--mode", "hifiasm", "--tolid", "x", "--taxon", "mollusca"],               # no contigs
    ["--mode", "hifiasm", "--tolid", "x", "--contigs", "/does/not/exist.fa", "--taxon", "mollusca"],
    ["--mode", "hifiasm", "--tolid", "x", "--contigs", "CONTIGS"],              # no taxon
])
def test_config_errors(root, tmp_path, files, args):
    args = [files["contigs.fa"] if a == "CONTIGS" else a for a in args]
    r = subprocess.run([sys.executable, os.path.join(root, "create_config.py"), "-o",
                        str(tmp_path / "c.yaml")] + args, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
    assert r.returncode != 0


def dryrun(root, cfg, tmp_path, jobs=False):
    """Dry-run; returns stdout, or {rule: job count} with jobs=True."""
    r = subprocess.run(["snakemake", "-n", "-s", os.path.join(root, "bin", "foam.smk"),
                        "--configfile", cfg, "-d", str(tmp_path / "run")],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    assert r.returncode == 0, r.stdout
    if not jobs:
        return r.stdout
    counts, in_table = {}, False
    for line in r.stdout.splitlines():
        f = line.split()
        if f[:2] == ["job", "count"]:
            in_table = True
            continue
        if in_table:
            if len(f) == 2 and f[1].isdigit() and f[0] != "total":
                counts[f[0]] = int(f[1])
            elif f and f[0] == "total":
                break
    return counts


@needs_snakemake
def test_dryrun_hifiasm(root, tmp_path, files):
    cfg, _ = make_config(root, tmp_path, "--mode", "hifiasm", "--tolid", "xgTest",
                         "--taxon", "mollusca", "--contigs", files["contigs.fa"])
    out = dryrun(root, cfg, tmp_path)
    assert "annotate_final" in out and "candidates" in out
    for r in ("prefilter_reads", "flye_assembly", "polish", "annotate_draft"):
        assert r not in out


@needs_snakemake
def test_dryrun_flye(root, tmp_path, files):
    cfg, _ = make_config(root, tmp_path, "--mode", "flye", "--tolid", "iTest",
                         "--taxon", "arthropoda", "--contigs", files["contigs.fa"],
                         "--flye-info", files["assembly_info.txt"], "--reads", files["reads.fq.gz"])
    out = dryrun(root, cfg, tmp_path)
    for r in ("candidates", "double_bait", "prefilter_reads", "select_reads",
              "annotate_draft", "polish", "annotate_final"):
        assert r in out
    for r in ("flye_assembly", "seed_proteins"):
        assert r not in out


@needs_snakemake
def test_dryrun_reads(root, tmp_path, files):
    cfg, _ = make_config(root, tmp_path, "--mode", "reads", "--tolid", "iTest",
                         "--taxon", "arthropoda", "--bait", files["bait.fa"],
                         "--reads", files["reads.fq.gz"])
    out = dryrun(root, cfg, tmp_path, jobs=True)
    # two fishing passes (bait, then own draft), one polish
    for r, n in (("double_bait", 2), ("prefilter_reads", 2), ("select_reads", 2),
                 ("flye_assembly", 2), ("candidates", 2), ("annotate_draft", 2),
                 ("polish", 1), ("annotate_final", 1)):
        assert out.get(r) == n, (r, out)
    assert "fetch_bait" not in out and "seed_proteins" not in out


@needs_snakemake
def test_dryrun_reads_gene_seed_single_pass(root, tmp_path, files):
    cfg, c = make_config(root, tmp_path, "--mode", "reads", "--tolid", "iTest",
                         "--taxon", "arthropoda", "--reads", files["reads.fq.gz"],
                         "--no-second-pass")
    assert c["fishing"]["seed"] == "auto"
    out = dryrun(root, cfg, tmp_path, jobs=True)
    assert out.get("seed_proteins") == 1 and "double_bait" not in out
    assert out.get("prefilter_reads") == 1 and out.get("flye_assembly") == 1


@needs_snakemake
def test_dryrun_reads_fetches_bait(root, tmp_path, files):
    cfg, _ = make_config(root, tmp_path, "--mode", "reads", "--tolid", "iTest",
                         "--taxon", "arthropoda", "--taxid", "6973", "--email", "x@y.org",
                         "--offline", "--reads", files["reads.fq.gz"])
    assert "fetch_bait" in dryrun(root, cfg, tmp_path)
