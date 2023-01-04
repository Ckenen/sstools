import pysam
from pyBioInfo.Utils import SegmentTools

def infer_library_layout(obj):
    is_paired = False
    if isinstance(obj, str):
        with pysam.AlignmentFile(obj) as f:
            for segment in f.fetch(until_eof=True):
                if segment.is_unmapped:
                    continue
                is_paired = segment.is_paired
                break
    elif isinstance(obj, pysam.AlignmentFile):
        for segment in obj.fetch(until_eof=True):
            if segment.is_unmapped:
                continue
            is_paired = segment.is_paired
            break
    else:
        raise ValueError("Unsupported input object.")
    return is_paired

def load_segments(bamfile, chrom, is_paired=False):
    if isinstance(bamfile, str):
        f = pysam.AlignmentFile(bamfile)
        need_close = True
    elif isinstance(bamfile, pysam.AlignmentFile):
        f = bamfile
        need_close = False
    else:
        raise RuntimeError()
    for segment in f.fetch(chrom):
        if segment.is_unmapped:
            continue
        if is_paired and (not segment.is_proper_pair):
            continue
        if segment.is_secondary:
            continue
        if segment.is_supplementary:
            continue
        yield segment
    if need_close:
        f.close()


def load_fragments(bamfile, chrom, is_paired):
    if is_paired:
        segments = list(sorted(load_segments(bamfile, chrom, True), key=lambda item: item.query_name))
        i = 0
        while i < len(segments) - 1:
            segment1 = segments[i]
            segment2 = segments[i + 1]
            if segment1.query_name == segment2.query_name:
                assert segment1.next_reference_start == segment2.reference_start
                assert segment2.next_reference_start == segment1.reference_start
                start = min(segment1.reference_start, segment2.reference_start)
                end = max(segment1.reference_end, segment2.reference_end)
                if segment1.is_read2:
                    segment1, segment2 = segment2, segment1
                yield [chrom, start, end, segment1, segment2]                        
                i += 2
            else:
                i += 1
    else:
        for segment in load_segments(bamfile, chrom, False):
            start = segment.reference_start
            end = segment.reference_end
            yield [chrom, start, end, segment]
    

def add_pg_to_header(header, cl):
    header = header.copy()
    if "PG" not in header:
        header["PG"] = list()
    pg_id_exists = [pg["ID"] for pg in header["PG"]]
    pg_id = "sstools"
    i = 1
    while pg_id in pg_id_exists:
        pg_id = "sstools.%d" % i
        i += 1
    header["PG"].append({'ID': pg_id, 'PN': 'sstools', 'VN': '1.0.0', 'CL': cl})
    return header


def divide_zero(v1, v2):
    # assert v1 <= v2
    if v2 > 0:
        return v1 / v2
    else:
        return 0
    
def load_bed_records(obj, chrom=None, start=None, end=None):
    raise NotImplementedError()
    # if isinstance(str, obj):
    #     pass
    # elif isinstance(pysam.TabixFile, obj):
    #     pass
    # else:
    #     raise RuntimeError()
    
# def get_mismatch_events(segment):
#     bases = SegmentTools.get_reference_sequence_base(segment)
#     rc = ";".join(["%s,%d" % (k, bases[k]) for k in sorted(bases.keys())])
#     items = []
#     for e in SegmentTools.get_events(segment):
#         if isinstance(e[3], list):
#             e[3] = "/".join(map(str, e[3]))
#         items.append(",".join(map(str, e)))
#     me = ";".join(items)
#     segment.set_tag("RC", rc)
#     segment.set_tag("ME", me)


class BaseMatrix(object):
    def __init__(self):
        self.pos = None # position
        self.ref = None # reference base
        self.allele = None # determine allele
        self.details = None # json detail
        # pos, ref, allele, allele_count, cell_count, read_count, score, conf, s2, s1

    def get_cell_count(self):
        pass

    def get_read_count(self):
        pass

    def get_score(self):
        pass


class PairedBaseMatrix(object):
    def __init__(self):
        self.pos = None # position
        self.ref = None # reference base
        self.ref1 = None # paternal allele
        self.ref2 = None # maternal allele
        self.conf = None # is high confidence
        self.allele1 = None
        self.allele2 = None
        self.details1 = None
        self.details2 = None
        # pos, ref, base_pat, base_mat, hc, base1, base2, s1, s2
