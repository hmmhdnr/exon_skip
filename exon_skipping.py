import os
import re
import pysam
import numpy as np
import pandas as pd
#from concurrent import futures
#from itertools import repeat
from tqdm import tqdm
import timeit
import argparse

from RefGene import getRefGeneRecords
from SAM import CIGAR_pos

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--refgene_path", help="path of refGene.txt", type=str, default="refGene.txt")
parser.add_argument("-s", "--refgene_start", help="line number to start searching refGene.txt", type=int, default=1)
parser.add_argument("-e", "--refgene_stop", help="line number to stop searching refGene.txt", type=int, default=1)
parser.add_argument("-n", "--sample_name", help="name of the sample analyzed", type=str, default="")
parser.add_argument("-i", "--input_bam", help="path of the input BAM file", type=str)
parser.add_argument("-o", "--output_dir", help="path of a directory where result files will be stored", type=str, default="result")
parser.add_argument("-c", "--cores", help="number of cores", type=int, default=1)
parser.add_argument("-t", "--threads", help="number of threads", type=int, default=2)
args = parser.parse_args()

def getReadsFromBAM(bam_file, chrom):
    bam = pysam.AlignmentFile(bam_file, "rb")
    if chrom == "all":
        return [x for x in bam]
    else:
        return [x for x in bam.fetch(contig=chrom)]

## -->
def check_mapped(read, exon_start, exon_end):
    if read.reference_end < exon_start:
        return 0
    elif read.reference_end <= exon_end:
        return 1
    else:
        if exon_end < read.reference_start:
            return 0
        elif exon_start <= read.reference_start:
            return 1
        else:
            cigar = CIGAR_pos(read.cigarstring)
            for c_pos in cigar:
                pos = read.reference_start + c_pos
                if exon_start <= pos and pos <= exon_end:
                    return 1
            return 0
    return 0

def _same_as_pl_check_mapped(read, exon_start, exon_end, i):
    if read.reference_start <= exon_start:
        if read.reference_end <= exon_start:
            return 0
        elif read.reference_end <= exon_end:
            return 1
        else:
            cigar_pos = CIGAR_pos(read.cigarstring)
            #pos = read.reference_start
            for c_pos in cigar_pos:
                #pos += c_pos #
                pos = read.reference_start + c_pos
                #if exon_start <= pos and pos <= exon_end:
                if exon_start <= pos and pos < exon_end: #
                    return 1
            return 0
    else:
        if read.reference_end <= exon_end:
            return 1
        elif read.reference_start <= exon_end:
            return 1
    return 0

def check_skipping(read, exon_start, exon_end):
    if read.reference_start < exon_start and exon_end < read.reference_end:
        cigar = CIGAR_pos(read.cigarstring)
        for c_pos in cigar:
            pos = read.reference_start + c_pos
            if exon_start <= pos and pos <= exon_end:
                return 0
            else:
                pass
        return 1
    else:
        return 0

def __check_skipping(read, exon_start, exon_end):
    if not(read.reference_start < exon_start and exon_end < read.reference_end):
        return 0
    cigar = CIGAR_pos(read.cigarstring)
    #pos = read.reference_start - 1
    #pos = read.reference_start #
    for c_pos in cigar:
        #pos += c_pos
        pos = read.reference_start + c_pos
        if exon_start <= pos and pos < exon_end:
        #if exon_start <= pos and pos <= exon_end: #
            return 0
    return 1

def _check_skipping(read, exon_start, exon_end):
    if not(read.reference_start < exon_start and exon_end < read.reference_end):
        return 0
    cigar_pos = CIGAR_pos(read.cigarstring)
    pos = read.reference_start
    for c_pos in cigar_pos:
        #pos += c_pos
        pos = read.reference_start + c_pos - 1
        if exon_start <= pos and pos <= exon_end:
        #if exon_start <= pos and pos < exon_end: #
            return 0
    return 1

def count_skipping(read, refgene):
    mapped = [check_mapped(read, ex_s, ex_e) for ex_s, ex_e in zip(refgene.exon_starts, refgene.exon_ends)]
    skipping = [check_skipping(read, ex_s, ex_e) for ex_s, ex_e in zip(refgene.exon_starts, refgene.exon_ends)]
    return (mapped, skipping)
## <--

def main():
    RGFile = args.refgene_path
    RGStart = args.refgene_start - 1
    RGStop = args.refgene_stop - 1

    result = pd.DataFrame(index=[], columns=["rg_acc", "rg_symbol", "rg_contig", "rg_strand", "mapped", "skipping"])
    result_file = "skipping_result_%d.tsv" % args.refgene_start

    # loading refGene
    refGene = getRefGeneRecords(RGFile, RGStart, RGStop)

    for rg in tqdm(refGene):
        print(rg.accession2)
        bam = pysam.AlignmentFile(args.input_bam, "rb")
        if len(rg.contig.split("_")) > 1:
            r = pd.Series({"rg_acc": rg.accession, "rg_symbol": rg.accession2, "rg_contig": rg.contig, "rg_strand": rg.strand, "mapped": "", "skipping": ""})
            result = pd.concat([result, r.to_frame().T])
            print(rg.accession, rg.accession2, rg.contig)
            print(rg.tr_start, rg.tr_end)
            continue

        reads_on_gene = [x for x in bam.fetch(contig=rg.contig) if type(x.reference_end) == int and rg.tr_start <= x.reference_start and x.reference_end <= rg.tr_end]
        if len(reads_on_gene) == 0:
            #
            read_count = ",".join(["0"] * rg.exon_count)
            #
            r = pd.Series({"rg_acc": rg.accession, "rg_symbol": rg.accession2, "rg_contig": rg.contig, "rg_strand": rg.strand, "mapped": read_count, "skipping": read_count})
            result = pd.concat([result, r.to_frame().T])
            continue

        #with futures.ThreadPoolExecutor(max_workers = args.threads) as executor:
        #    future_list = [executor.submit(count_skipping, read, refgene=rg) for read in reads_on_gene]
        #    results = [x.result() for x in future_list]
        #    mapped = [np.array(x[0]) for x in results]
        #    skipping = [np.array(x[1]) for x in results]
        results = [count_skipping(read, rg) for read in reads_on_gene]
        mapped = [np.array(x[0]) for x in results]
        skipping = [np.array(x[1]) for x in results]

        if rg.exon_count > 0:
            m = ",".join([str(x) for x in np.array(sum(mapped))])
            s = ",".join([str(x) for x in np.array(sum(skipping))])
        else:
            m = ""
            s = ""
        r = pd.Series({"rg_acc": rg.accession, "rg_symbol": rg.accession2, "rg_contig": rg.contig, "rg_strand": rg.strand, "mapped": m, "skipping": s})
        result = pd.concat([result, r.to_frame().T])

    output_dir = os.path.join(".", args.output_dir, args.sample_name)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    result.to_csv(os.path.join(output_dir, result_file), sep="\t")


if __name__ == "__main__":
    print("counting exon-skipping in %s" % args.sample_name)
    t0 = timeit.default_timer()
    main()
    t1 = timeit.default_timer()
    print("elapsed time = %.3f sec." % (t1 - t0))
