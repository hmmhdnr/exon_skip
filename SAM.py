import re

def _CIGAR_pos(cigarstring):
    cigar = re.split("([a-zA-Z]+)", cigarstring)[:-1]
    rel_pos = []
    start_pos = 0
    x = []
    for l, stat in zip(cigar[0::2], cigar[1::2]):
        if stat == "M":
            rel_pos.append(start_pos)
        if stat in ["M", "D", "N"]:
            x.append(start_pos)
            start_pos += int(l)
    return rel_pos

def CIGAR_pos(cigarstring):
    cigar = re.split("([a-zA-Z]+)", cigarstring)[:-1]
    rel_pos = []
    start_pos = 0
    x = []
    for l, stat in zip(cigar[0::2], cigar[1::2]):
        if stat == "M":
            rel_pos.append(start_pos)
        if stat in ["M", "D", "N"]:
            start_pos += int(l)
    return rel_pos

