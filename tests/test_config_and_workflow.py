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


def dryrun(root, cfg, tmp_path):
    r = subprocess.run(["snakemake", "-n", "-s", os.path.join(root, "bin", "foam.smk"),
                        "--configfile", cfg, "-d", str(tmp_path / "run"), "--quiet", "rules"],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    assert r.returncode == 0, r.stdout
    return r.stdout


@needs_snakemake
def test_dryrun_hifiasm(root, tmp_path, files):
    cfg, _ = make_config(root, tmp_path, "--mode", "hifiasm", "--tolid", "xgTest",
                         "--taxon", "mollusca", "--contigs", files["contigs.fa"])
    out = dryrun(root, cfg, tmp_path)
    assert "annotate_final" in out and "candidates" in out
    for r in ("fish_reads", "flye_assembly", "polish", "annotate_prepolish"):
        assert r not in out


@needs_snakemake
def test_dryrun_flye(root, tmp_path, files):
    cfg, _ = make_config(root, tmp_path, "--mode", "flye", "--tolid", "iTest",
                         "--taxon", "arthropoda", "--contigs", files["contigs.fa"],
                         "--flye-info", files["assembly_info.txt"], "--reads", files["reads.fq.gz"])
    out = dryrun(root, cfg, tmp_path)
    for r in ("candidates", "double_bait", "fish_reads", "annotate_prepolish", "polish",
              "annotate_final"):
        assert r in out
    assert "flye_assembly" not in out


@needs_snakemake
def test_dryrun_reads(root, tmp_path, files):
    cfg, _ = make_config(root, tmp_path, "--mode", "reads", "--tolid", "iTest",
                         "--taxon", "arthropoda", "--bait", files["bait.fa"],
                         "--reads", files["reads.fq.gz"])
    out = dryrun(root, cfg, tmp_path)
    for r in ("double_bait", "fish_reads", "flye_assembly", "candidates",
              "annotate_prepolish", "polish", "annotate_final"):
        assert r in out
    assert "fetch_bait" not in out


@needs_snakemake
def test_dryrun_reads_fetches_bait(root, tmp_path, files):
    cfg, _ = make_config(root, tmp_path, "--mode", "reads", "--tolid", "iTest",
                         "--taxon", "arthropoda", "--taxid", "6973", "--email", "x@y.org",
                         "--offline", "--reads", files["reads.fq.gz"])
    assert "fetch_bait" in dryrun(root, cfg, tmp_path)
