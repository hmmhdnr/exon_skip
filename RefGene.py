import sys

class RefGene:
    def __init__(self, line_num, line):
        self.line_id = line_num
        elem = line.split("\t")
        acc = elem[1]
        contig = elem[2]
        strand = elem[3]
        tr_start = int(elem[4]) + 1
        tr_end = int(elem[5])
        cds_start = int(elem[6]) + 1
        cds_end = int(elem[7])
        exon_count = int(elem[8])
        exon_starts = [int(x) + 1 for x in elem[9].split(",")[:-1]]
        exon_ends = [int(x) for x in elem[10].split(",")[:-1]]
        _ = elem[11]
        acc2 = elem[12]
        self.accession = acc
        self.contig = contig
        self.strand = strand
        self.tr_start = tr_start
        self.tr_end = tr_end
        self.cds_start = cds_start
        self.cds_end = cds_end
        self.exon_count = exon_count
        self.exon_starts = exon_starts
        self.exon_ends = exon_ends
        self.accession2 = acc2

    def show(self):
        return "\t".join([self.accession, self.contig, self.strand,
                         str(self.tr_start-1), str(self.tr_end),
                         str(self.cds_start-1), str(self.cds_end), str(self.exon_count),
                         ",".join([str(x) for x in self.exon_starts]),
                         ",".join([str(x) for x in self.exon_ends]),
                          self.accession2])

def getSingleRefGeneRecord(rg_file, line_num):
    with open(rg_file, "r") as f:
        lines = f.read().splitlines()
    return RefGene(line_num, lines[line_num])

def getRefGeneRecords(rg_file, start, stop):
    if start > stop:
        print("Error: invalid line number of refGene.txt was specified.")
        sys.exit()

    with open(rg_file, "r") as f:
        lines = f.read().splitlines()
        stop = min(stop + 1, len(lines))
        res = [RefGene(line_num, lines[line_num]) for line_num in range(start, stop, 1)]
    return res
