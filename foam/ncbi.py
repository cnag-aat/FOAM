"""Closest available mitogenome(s) for a taxon (port of get_closest_mitogenome.py).

Walks up the NCBI taxonomy from the given taxid (species -> genus -> family ->
order) and returns up to ``max_results`` complete mitochondrial genomes from
the first rank that has any, preferring RefSeq records.  Uses only the
standard library (E-utilities over HTTPS).  Also reports the taxon's
mitochondrial genetic code and a matching FOAM taxon preset.
"""

import time
import xml.etree.ElementTree as ET

try:
    from urllib.parse import urlencode
    from urllib.request import urlopen
except ImportError:  # pragma: no cover
    raise

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
RANKS = ["species", "genus", "family", "order"]

# lineage name (lower case) -> FOAM taxon preset, checked most specific first
LINEAGE_PRESETS = [
    ("vertebrata", "vertebrata"), ("tunicata", "tunicata"),
    ("echinodermata", "echinodermata"), ("platyhelminthes", "platyhelminthes"),
    ("mollusca", "mollusca"), ("arthropoda", "arthropoda"),
    ("annelida", "annelida"), ("nematoda", "nematoda"),
    ("cnidaria", "cnidaria"), ("porifera", "porifera"),
]


class Entrez(object):
    def __init__(self, email, api_key=None, tool="foam", opener=urlopen, delay=0.4):
        self.base = {"tool": tool, "email": email}
        if api_key:
            self.base["api_key"] = api_key
            delay = 0.12
        self.opener = opener
        self.delay = delay

    def get(self, util, **params):
        params.update(self.base)
        url = EUTILS + util + "?" + urlencode(params)
        for attempt in range(3):
            try:
                fh = self.opener(url)
                data = fh.read()
                fh.close()
                time.sleep(self.delay)
                return data.decode() if isinstance(data, bytes) else data
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(2 * (attempt + 1))


def parse_taxonomy(xml_text):
    root = ET.fromstring(xml_text)
    taxon = root.find("Taxon")
    if taxon is None:
        raise ValueError("no taxonomy record")
    lineage = [{"TaxId": t.findtext("TaxId"), "ScientificName": t.findtext("ScientificName"),
                "Rank": t.findtext("Rank")} for t in taxon.findall("LineageEx/Taxon")]
    me = {"TaxId": taxon.findtext("TaxId"), "ScientificName": taxon.findtext("ScientificName"),
          "Rank": taxon.findtext("Rank")}
    mgc = taxon.findtext("MitoGeneticCode/MGCId")
    return {"taxon": me, "lineage": [me] + list(reversed(lineage)),
            "mito_genetic_code": int(mgc) if mgc else None}


def preset_for_lineage(lineage):
    names = set((n["ScientificName"] or "").lower() for n in lineage)
    for lname, preset in LINEAGE_PRESETS:
        if lname in names:
            return preset
    return "invertebrate"


def parse_fasta_text(text):
    recs, name, title, buf = [], None, "", []
    for line in text.splitlines():
        if line.startswith(">"):
            if name:
                recs.append((name, title, "".join(buf)))
            parts = line[1:].split(None, 1)
            name, title, buf = parts[0], parts[1] if len(parts) > 1 else "", []
        elif line.strip():
            buf.append(line.strip())
    if name:
        recs.append((name, title, "".join(buf)))
    return recs


def taxonomy(taxid, entrez):
    return parse_taxonomy(entrez.get("efetch.fcgi", db="taxonomy", id=str(taxid), retmode="xml"))


def closest_mitogenomes(taxid, email, max_results=3, api_key=None, entrez=None):
    entrez = entrez or Entrez(email, api_key)
    tax = taxonomy(taxid, entrez)
    out = {"taxid": str(taxid), "name": tax["taxon"]["ScientificName"],
           "mito_genetic_code": tax["mito_genetic_code"],
           "taxon_preset": preset_for_lineage(tax["lineage"]), "records": []}
    for node in tax["lineage"]:
        if node["Rank"] not in RANKS:
            continue
        base = ("txid%s[Organism:exp] AND mitochondrion[Filter] AND complete genome[Title]"
                % node["TaxId"])
        ids = []
        for term in (base + " AND refseq[Filter]", base):
            root = ET.fromstring(entrez.get("esearch.fcgi", db="nuccore", term=term, retmax="100"))
            ids = [e.text for e in root.findall("IdList/Id")]
            if ids:
                break
        if not ids:
            continue
        if len(ids) > max_results:   # spread picks across the hit list
            step = len(ids) // max_results
            ids = [ids[i * step] for i in range(max_results)]
        fasta = entrez.get("efetch.fcgi", db="nuccore", id=",".join(ids),
                           rettype="fasta", retmode="text")
        for acc, title, seq in parse_fasta_text(fasta):
            out["records"].append({"accession": acc, "title": title, "sequence": seq,
                                   "length": len(seq)})
        out["rank"] = node["Rank"]
        out["rank_name"] = node["ScientificName"]
        return out
    raise RuntimeError("no complete mitogenome found for taxid %s up to order level" % taxid)
