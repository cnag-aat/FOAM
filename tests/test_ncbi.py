import io

from foam import ncbi

TAX_XML = """<?xml version="1.0"?>
<TaxaSet><Taxon>
  <TaxId>999</TaxId><ScientificName>Xerocrassa testis</ScientificName><Rank>species</Rank>
  <MitoGeneticCode><MGCId>5</MGCId><MGCName>Invertebrate Mitochondrial</MGCName></MitoGeneticCode>
  <LineageEx>
    <Taxon><TaxId>6447</TaxId><ScientificName>Mollusca</ScientificName><Rank>phylum</Rank></Taxon>
    <Taxon><TaxId>40</TaxId><ScientificName>Geomitridae</ScientificName><Rank>family</Rank></Taxon>
    <Taxon><TaxId>41</TaxId><ScientificName>Xerocrassa</ScientificName><Rank>genus</Rank></Taxon>
  </LineageEx>
</Taxon></TaxaSet>"""


def esearch(ids):
    return "<eSearchResult><Count>%d</Count><IdList>%s</IdList></eSearchResult>" % (
        len(ids), "".join("<Id>%s</Id>" % i for i in ids))


class FakeOpener(object):
    def __init__(self):
        self.urls = []

    def __call__(self, url):
        self.urls.append(url)
        if "db=taxonomy" in url:
            body = TAX_XML
        elif "esearch" in url and "txid999" in url:
            body = esearch([])                       # nothing at species level
        elif "esearch" in url and "txid41" in url and "refseq" in url:
            body = esearch([])                       # genus: no RefSeq ...
        elif "esearch" in url and "txid41" in url:
            body = esearch(["11", "12"])             # ... but GenBank records
        elif "efetch" in url and "db=nuccore" in url:
            body = ">MZ000011.1 Xerocrassa sp. mitochondrion, complete genome\nACGT\nAC\n" \
                   ">MZ000012.1 Xerocrassa sp2 mitochondrion, complete genome\nGGGG\n"
        else:
            raise AssertionError(url)
        return io.BytesIO(body.encode())


def test_closest_mitogenomes_walks_up_to_genus():
    op = FakeOpener()
    res = ncbi.closest_mitogenomes("999", "x@y.org", entrez=ncbi.Entrez("x@y.org", opener=op, delay=0))
    assert res["rank"] == "genus" and res["rank_name"] == "Xerocrassa"
    assert [r["accession"] for r in res["records"]] == ["MZ000011.1", "MZ000012.1"]
    assert res["records"][0]["sequence"] == "ACGTAC"
    assert res["mito_genetic_code"] == 5 and res["taxon_preset"] == "mollusca"
    assert all("email=x%40y.org" in u for u in op.urls)


def test_preset_for_lineage():
    lin = [{"ScientificName": "Vertebrata"}, {"ScientificName": "Chordata"}]
    assert ncbi.preset_for_lineage(lin) == "vertebrata"
    assert ncbi.preset_for_lineage([{"ScientificName": "Chordata"}]) == "invertebrate"
