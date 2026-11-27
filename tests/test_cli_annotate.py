import json
import os
import sys

import pytest

from conftest import rand_seq
from foam import seqio
from foam.cli import main

FAKE = os.path.join(os.path.dirname(__file__), "data", "fake_runmitos.py")


@pytest.fixture
def genome(tmp_path, monkeypatch):
    """6 kb toy mitogenome: cox1 on the minus strand, trnF and nad1 on plus."""
    cox1, trnf, nad1 = rand_seq(600, 21), rand_seq(70, 22), rand_seq(300, 23)
    spacer = [rand_seq(n, 30 + i) for i, n in enumerate([1000, 1500, 800, 1730])]
    mono = spacer[0] + seqio.revcomp(cox1) + spacer[1] + trnf + spacer[2] + nad1 + spacer[3]
    markers = tmp_path / "markers.fa"
    seqio.write_fasta([("cox1", "gene", cox1), ("trnF", "tRNA", trnf), ("nad1", "gene", nad1)],
                      str(markers))
    monkeypatch.setenv("FAKE_MITOS_MARKERS", str(markers))
    return {"mono": mono, "cox1": cox1, "trnF": trnf}


def run(tmp_path, contig_fa, *extra):
    out = str(tmp_path / "out" / "sp_MT")
    argv = ["annotate", "-i", contig_fa, "-o", out, "--refdir", str(tmp_path),
            "--runmitos", "%s %s" % (sys.executable, FAKE), "--workdir", str(tmp_path / "w"),
            "--min-period", "2000", "--expected-genes", "cox1,trnF,nad1",
            "--no-rescue"] + list(extra)
    assert main(argv) == 0
    return out


def test_tandem_contig_mollusc_starts_at_cox1(tmp_path, genome):
    fa = str(tmp_path / "c.fa")
    seqio.write_fasta([("ptg000001c", seqio.rotate(genome["mono"] * 2, 777))], fa)
    out = run(tmp_path, fa, "--taxon", "mollusca")
    seq = seqio.read_fasta(out + ".fa")[0][2]
    assert len(seq) == len(genome["mono"])
    assert seq.startswith(genome["cox1"])            # 5' end of cox1, + strand
    rep = json.load(open(out + ".report.json"))
    assert rep["orientation"]["anchor"] == "cox1"
    assert rep["orientation"]["reverse_complemented"] is True
    assert rep["qc"]["n_present"] == 3 and rep["qc"]["missing"] == []
    assert rep["contigs"][0]["n_copies"] == 2
    # identical copies are annotated once
    assert sum(1 for c in rep["candidates"] if c.get("same_as")) == 1
    gff = open(out + ".gff3").read()
    assert "transl_table=5" in gff and "\tCDS\t1\t600\t" in gff
    assert not os.path.exists(str(tmp_path / "w"))     # temp files removed


def test_vertebrate_starts_at_trnF(tmp_path, genome):
    fa = str(tmp_path / "c.fa")
    seqio.write_fasta([("contig_1", genome["mono"])], fa)
    out = run(tmp_path, fa, "--taxon", "vertebrata")
    seq = seqio.read_fasta(out + ".fa")[0][2]
    assert seq.startswith(genome["trnF"])
    assert "transl_table=2" in open(out + ".gff3").read()


def test_best_copy_selected_among_contigs(tmp_path, genome):
    broken = genome["mono"].replace(genome["trnF"], rand_seq(70, 99))
    fa = str(tmp_path / "c.fa")
    seqio.write_fasta([("bad", broken), ("good", seqio.rotate(genome["mono"], 100))], fa)
    out = run(tmp_path, fa, "--taxon", "arthropoda", "--keep-tmp")
    rep = json.load(open(out + ".report.json"))
    assert rep["selected"].startswith("good")


def test_mitos_command_auto(tmp_path, monkeypatch):
    from foam import mitos as mt
    bindir = tmp_path / "bin"
    bindir.mkdir()
    exe = bindir / "runmitos"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", str(bindir))
    assert mt.mitos_command("auto") == ["runmitos"]
    assert mt.mitos_command("runmitos.py") == ["runmitos"]      # old config, new MITOS
    assert mt.mitos_command("python x.py") == ["python", "x.py"]
    exe.rename(bindir / "runmitos.py")
    assert mt.mitos_command("runmitos") == ["runmitos.py"]
    (bindir / "runmitos.py").unlink()
    import pytest
    with pytest.raises(mt.MitosError):
        mt.mitos_command("auto")
