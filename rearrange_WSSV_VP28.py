#!/usr/bin/env python

# the WSSV genome is circular
# thus, assemblies often start at different positions
# this is not ideal for alignments
# want to rearrange all genomes to start at VP10
# this is where the first WSSV genome started (AF369029)
# will achieve this by blasting the VP10 gene against the genome
# and rearranging from that

import sys, os
import argparse
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import glob

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("take a newly assembled WSSV genome and re-arrange it so that it starts at the VP10 "
                                 "region\nnote: requires blast+ tools to be available on the command line\nnote: "
                                 "requires that biopython is installed in python\n")

parser.add_argument('-g', '--genome', type = str,
        nargs = 1, help = "fasta file with genome to be re-arranged")
parser.add_argument('-o', '--output', type = str,
        nargs = 1, help = "new name for the re-arranged genome")

if len(sys.argv) == 1:  # if no args are given
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check required arguments are provided

if args.genome is None or args.output is None:
    print("*** error: required argument is missing"
          "*** error: input genome and an output file name is required")
    parser.print_help(sys.stderr)
    sys.exit(1)

# check that blast is loaded

try:
    tmp = subprocess.call(["blastn", "-version"])
except OSError as e:
    if e.errno == os.errno.ENOENT:
        print("*** error: blastn could not be found: try 'module load blast+'\n")
        raise

# make a blast database from the genome

subprocess.call(["makeblastdb", "-dbtype", "nucl", "-in", args.genome[0]])

# create and write out a WSSV VP10 record

VP10_record = SeqRecord(Seq("tggatctttctttcactctttcggtcgtgtcggccatcctcgccatcactgctgtgattgctgtatttattgtgatttttaggtat"
                            "cacaacactgtgaccaagaccatcgaaacccacacagacaatatcgagacaaacatggatgaaaacctccgcattcctgtgactgc"
                            "tgaggttggatcaggctacttcaagatgactgatgtgtcctttgacagcgacaccttgggcaaaatcaagatccgcaatggaaagt"
                            "ctgatgcacagatgaaggaagaagatgcggatcttgtcatcactcccgtggagggccgagcactcgaagtgactgtggggcagaat"
                            "ctcacctttgagggaacattcaaggtgtggaacaacacatcaagaaagatcaacatcactggtatgcagatggtgccaaagattaa"
                            "cccatcaaaggcctttgtcggtagctccaacacctcctccttcacccccgtctctattgatgaggatgaagttggcacctttgtgt"
                            "gtggtaccacctttggcgcaccaattgcagctaccgccggtggaaatcttttcgacatgtacgtgcacgtcacctactctggcact"
                            "gagaccgagtaaataaatcgtgcttttttatatagatagggaattttaatattacaa"), id="AAK77670.1_nucl",
                            name="VP10", description="WSSV_VP10_ORF1")

SeqIO.write(VP10_record, "VP10_record.fasta", "fasta")

# now blast the VP10 gene against the WSSV genome

blast_output = subprocess.check_output(["blastn", "-query", "VP10_record.fasta", "-db", args.genome[0], "-outfmt",
                                        "6 length pident sstart send sstrand"])

blast_output = blast_output.decode("utf-8")

if len(blast_output.split("\n")) > 2:
    print("*** warning: more than 1 blast match found")
    print("*** warning: using only the top blast match to re-arrange genome")
    blast_output = blast_output.split("\n")[0]

if not blast_output:
    print("*** error: could not find a WSSV VP10 match in the genome!")
    sys.exit(1)
else:
    blast_cols = blast_output.split()
    length = blast_cols[0]
    identity = blast_cols[1]
    start = blast_cols[2]
    end = blast_cols[3]
    strand = blast_cols[4]
    print("~~~ VP10 gene found at position {}-{} in the {} strand with {}% identity over {} bps (659 bp total "
          "length)".format(start, end, strand, identity, length))

# create genome record and reverse complement if required

genome_record = SeqIO.read(args.genome[0], "fasta")

if strand == "minus":
    print("~~~ reverse complementing genome")
    genome_record = genome_record.reverse_complement(id=True, name=True, description=True)
    # also need to change the new start position
    start = len(genome_record.seq) - (int(start) - 1) # python 0-based indexing fix

if float(identity) < 80:
    print("*** warning: identity is low for the VP10 gene")
    # TODO: maybe stop re-arrangement if identity is too low

# do the actual re-arrangement in biopython

def rearrange_genome(record, new_start):
    print("~~~ re-arranging genome with new start position")
    new_start = new_start - 1 # match blast / python indexing
    start_chunk = record[new_start:]
    end_chunk = record[:new_start]
    new_record = start_chunk + end_chunk
    return(new_record)

rearranged_genome = rearrange_genome(genome_record, int(start))
SeqIO.write(rearranged_genome, args.output[0], "fasta")

# clean up database and VP10 files

files_to_remove = [args.genome[0] + ext for ext in [".nhr", ".nin", ".nsq"]]
files_to_remove.append("VP10_record.fasta")

for fl in files_to_remove:
    os.remove(fl)


print("~~~ Genome {} was successfully re-arranged to begin at position {} and written to the file {}".format(
    args.genome[0], start, args.output[0]))
