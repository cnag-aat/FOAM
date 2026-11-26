"""NCBI genetic codes relevant to mitochondria (translation tables)."""

_B1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
_B2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
_B3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

# table: (amino acids, start codons marked 'M')
_TABLES = {
    1: ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "---M------**--*----M---------------M----------------------------"),
    2: ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        "--------------------------------MMMM---------------M------------"),
    4: ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "--MM---------------M------------MMMM---------------M------------"),
    5: ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        "---M----------------------------MMMM---------------M------------"),
    9: ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "-----------------------------------M---------------M------------"),
    13: ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
         "---M------------------------------MM---------------M------------"),
    14: ("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
         "-----------------------------------M----------------------------"),
}


class GeneticCode(object):
    def __init__(self, table):
        if table not in _TABLES:
            raise ValueError("genetic code %s not supported" % table)
        aas, starts = _TABLES[table]
        self.table = table
        self.codon = {}
        self.starts = set()
        self.stops = set()
        for i in range(64):
            c = _B1[i] + _B2[i] + _B3[i]
            self.codon[c] = aas[i]
            if starts[i] == "M":
                self.starts.add(c)
            if aas[i] == "*":
                self.stops.add(c)

    def translate(self, seq):
        seq = seq.upper()
        return "".join(self.codon.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))
