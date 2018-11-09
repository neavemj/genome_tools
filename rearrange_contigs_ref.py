#!/usr/bin/env python

# rearrange and reverse complement contigs against a given reference genome
# Matthew J. Neave

import sys, os
import argparse
import subprocess
import genomeview
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import glob

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("take viral genome(s) and calculate coordinates for the contigs along a given "
                                 "reference genome"
                                 "\nnote: requires mummer to be available on the command line\nnote: "
                                 "requires that biopython is installed in python\n")

parser.add_argument('-g', '--genomes', type = str,
                    nargs = "*", help = "fasta files with genome to be re-arranged. Multiple fasta files can be given as a space-separated list")
parser.add_argument('-r', '--reference_genome', type = str,
                    nargs = 1, help = "fasta file with reference genome used as a backbone")
parser.add_argument('-o', '--output', type = str,
                    nargs = 1, help = "name for co-ordinate file")
parser.add_argument('-a', '--arranged_contigs', type = str,
                    nargs = 1, help = "optional: output name if re-arranged contigs required")
parser.add_argument('--plot', action = 'store_true',
                    help = "plot the resulting contig tiling?")

if len(sys.argv) == 1:  # if no args are given
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check required arguments are provided

if args.genomes is None or args.output is None or args.reference_genome is None:
    print("*** error: required argument is missing"
          "*** error: input genome, reference genome and an output file name are required")
    parser.print_help(sys.stderr)
    sys.exit(1)

# check that dependancies are loaded

try:
    tmp = subprocess.call(["nucmer", "-version"])
except:
    print("*** error: nucmer could not be found: trying to 'module load mummer'\n")
    subprocess.call(["module", "load", "mummer"])

try:
    tmp = subprocess.call(["R", "--version"])
except:
    subprocess.call(["module", "load", "R"])

# function to make a dictionary of contigs from a multi-fasta genome

def make_genome_dict(gen):
    genome_dict = {}
    for contig in SeqIO.parse(gen, "fasta"):
        genome_dict[contig.id] = contig.seq
    return(genome_dict)

# align contigs to the reference using mummer and nucmer (nucleotide alignment)
# if the genomes are more divergent they can be aligned using promer (protein alignment)

out_handle = open(args.output[0], "w")
out_handle.write("\t".join(["sample", "contig_name", "start", "end", "length", "coverage", "identity",
"orientation", "wrapped"]) + "\n")

for genome in args.genomes:
    subprocess.check_output(["nucmer", "--prefix=" + genome, args.reference_genome[0],
                                             genome])

    # create contig 'tiles' along the genome
    tiling_output = subprocess.check_output(["show-tiling", "-c", "-g", "-1", genome + ".delta"])
    tiling_output = tiling_output.decode("utf-8")

    # need to re-run mummer but without the circularisation flag (-g)
    # this will make it easier to re-arrange them appropriately
    if args.arranged_contigs:
        rearrange_output = subprocess.check_output(["show-tiling", "-g", "-1", genome + ".delta"])
        rearrange_output = rearrange_output.decode("utf-8")
        print(rearrange_output)
        genome_dict = make_genome_dict(genome)
        new_genome = open(args.arranged_contigs[0], "w")
        genome_name = genome.split(".")[0]
        new_genome.write(">" + genome_name + "_contigs_arranged" + "\n")
        for contig in tiling_output.split("\n"):
            contig = contig.strip().split("\t")
            if len(contig) < 7:
                continue
            contig_name = contig[7]
            strand = contig[6]
            if strand == "+":
                new_genome.write(str(genome_dict[contig_name]))
            if strand == "-":
                rev_comp = genome_dict[contig_name].reverse_complement()
                new_genome.write(str(rev_comp))


    # b'>AF369029.2 292967 bases\n-195054\t74406\t2360\t269461\t98.19\t99.63\t+\tscaffold_0\n76767\t83037\t14875\t6271
    # \t100.00\t99.86\t+\tscaffold_1\n'

    for contig in tiling_output.split("\n"):
        contig = contig.strip().split("\t")
        if len(contig) < 7:
            continue
        contig_name = contig[7]
        genome_name = genome.split(".")[0]
        start = int(contig[0])
        end = int(contig[1])
        length = int(contig[3])
        wrapped = "False"
        if start < 0:
            # write contig from first-half of the wrap
            wrapped = "True"
            out_handle.write("\t".join([genome_name, contig_name, "0", contig[1], contig[3], contig[4],
                                        contig[5], contig[6], wrapped]) + "\n")

            # write contig from back-half of the wrap
            ref_record = SeqIO.read(args.reference_genome[0], "fasta")
            ref_len = len(ref_record.seq)
            start = ref_len -- start
            out_handle.write("\t".join([genome_name, contig_name, str(start), str(ref_len), contig[3], contig[4],
                                        contig[5], contig[6], wrapped]) + "\n")

        else:
            out_handle.write("\t".join([genome_name, contig_name, str(start), contig[1], contig[3], contig[4], contig[5],
                                    contig[6], wrapped]) + "\n")


out_handle.close()

# plot the orientation of the draft contigs against the reference genome

if args.plot:
    subprocess.call(["/datastore/nea040/PycharmProjects/genome_tools/plot_mummer.R", args.output[0]])

