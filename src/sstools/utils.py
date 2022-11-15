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
    elif isinstance(obj. pysam.AlignmentFile):
        for segment in f.fetch(until_eof=True):
            if segment.is_unmapped:
                continue
            is_paired = segment.is_paired
            break
    else:
        raise ValueError("Unsupported input object.")
    return is_paired

def add_pg_to_header(header, cl):
    header = header.copy()
    if "PG" not in header:
        header["PG"] = list()
    header["PG"].append({'ID': 'sstools', 'PN': 'sstools', 'VN': '1.0.0', 'CL': cl})
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
    
def get_mismatch_events(segment):
    bases = SegmentTools.get_reference_sequence_base(segment)
    rc = ";".join(["%s,%d" % (k, bases[k]) for k in sorted(bases.keys())])
    items = []
    for e in SegmentTools.get_events(segment):
        if isinstance(e[3], list):
            e[3] = "/".join(map(str, e[3]))
        items.append(",".join(map(str, e)))
    me = ";".join(items)
    segment.set_tag("RC", rc)
    segment.set_tag("ME", me)
