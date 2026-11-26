"""Select mitochondrial reads from a minimap2 PAF stream.

The bait is doubled (``seq+seq``) so that reads crossing the linearisation
point of the circular bait still align end-to-end.  A read is kept only when
its alignments (merged) cover most of the read; reads that merely *contain* a
mito-like segment (e.g. nuclear reads carrying a NUMT) are rejected.  With a
bait from a related species, alignments can be split at divergent regions,
which is why coverage is the union over all alignments of the read.
"""


def double_bait(records):
    """records: (name, desc, seq) -> doubled (name, desc, seq)."""
    return [(n, "doubled circular bait", s + s) for n, _d, s in records]


def _union_len(intervals):
    total = 0
    cur_s = cur_e = None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                total += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s
    return total


def read_passes(qlen, intervals, matches, block, min_frac=0.8, min_aln=1000,
                min_identity=0.7, min_len=1000, max_len=0):
    if qlen < min_len or (max_len and qlen > max_len):
        return False
    covered = _union_len(intervals)
    if covered < min_aln or covered < min_frac * qlen:
        return False
    if block and float(matches) / block < min_identity:
        return False
    return True


def filter_paf(lines, **kw):
    """Yield IDs of reads passing the filter.

    minimap2 writes all alignments of a read consecutively, so reads are
    processed one group at a time with O(1) memory.
    """
    cur = None
    qlen = matches = block = 0
    ivs = []
    for line in lines:
        c = line.split("\t", 12)
        if len(c) < 12:
            continue
        rid = c[0]
        if rid != cur:
            if cur is not None and read_passes(qlen, ivs, matches, block, **kw):
                yield cur
            cur, qlen, matches, block, ivs = rid, int(c[1]), 0, 0, []
        ivs.append((int(c[2]), int(c[3])))
        matches += int(c[9])
        block += int(c[10])
    if cur is not None and read_passes(qlen, ivs, matches, block, **kw):
        yield cur
