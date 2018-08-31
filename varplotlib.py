import sys
import argparse
import csv
import subprocess
import os

from collections import Counter

from variant import Variant
from sample import Sample

accepted_types = ["nonsynonymous_SNV", "frameshift_deletion", "frameshift_insertion"
                  "nonframeshift_deletion", "stopgain", "splicing"]

command = 'Rscript'

def main(args):
    # Open vcf and output mutation statistics like in test_data.csv
    type = check_arguments(args)
    if type == 'table':
        read_table(args)
    elif type == 'list':
        samples, group_names = read_list(args)
    #print(samples)
    global current_path 
    current_path = os.path.dirname(os.path.realpath(__file__))
    print(current_path)
    substitutions_csv = prepare_base_substitutions(samples, group_names, args.output)
    plot_substitution_frequencies(substitutions_csv)
    gene_matrix_filename = create_gene_matrix(samples, args)
    var_n = amount_of_variants(samples, args)
    plot_gene_matrix(gene_matrix_filename, var_n, args)
    

def check_arguments(args):
    suffix = args.input.split('.')[-1]
    if suffix == 'table':
        if args.groups == None:
            print("Error! Group-file must be provided with table-input. "
            "Exiting...\n ")
            exit()
    elif suffix == 'list':
        if args.groups != None:
            print("List input doesn't need group-file.")
    else:
        print("Error! Invalid input file. Exiting..." )
        exit()
    if args.rename == None:
        args.rename = "FALSE"
    if args.zero_variants == None:
        args.zero_variants = "FALSE"
    return suffix

def read_table(args):
    pass

def read_list(args):
    # Dict where each key is a group in the vcf list
    # Value is a list of Mutation-objects
    samples = []
    group_names = []
    with open(args.input,'r') as list_of_files:
        for line in list_of_files:
            try:
                line = line.split('\t')
                file_name = line[0]
                sample_name = line[1]
                # sample_names.append(sample_name)
                variants = []
                group = line[2].strip()
                if group not in group_names:
                    group_names.append(group)
                with open(file_name, 'r') as vcf:
                    for vline in vcf:
                        if vline[0] != '#':
                            # Input split line to parse_variants
                            site_variants = parse_variants(vline.split('\t'))
                            for variant in site_variants:
                                if variant.type in accepted_types:
                                    variants.append(variant)
                new_sample = Sample(sample_name, variants, group)
                samples.append(new_sample)
            except IndexError:
                print("Error! Badly formatted input list.")
    return samples, group_names

def parse_variants(line):
    # Parses variant information from annovar-annotated vcf-file.
    variants = []
    chromosome = line[0]
    start = line[1]
    end = line[2]
    ref = line[3]
    alt = line[4]
    info = line[7]
    genes, type, af = parse_info(info.split(';'))
    for g in genes:
        new_var = Variant(chromosome, start, end, ref, alt, g, type, af, None)
        variants.append(new_var)
    return variants

def parse_info(info):
    genes = []
    type, af = None, None
    exonic = False
    for i in info:
        i = i.split('=')
        if i[0] == 'AF':
            af = i[1]
        elif i[0] == 'Gene.refGene':
            gs = i[1].split("\\x3b")
            for g in gs:
                genes.append(g)
        elif i[0] == 'Func.refGene':
            if i[1] == "exonic":
                exonic = True
            else:
                type = i[1]
        elif i[0] == "ExonicFunc.refGene":
            if exonic:
                type = i[1]
        else:
            pass
    return genes, type, af

def create_gene_matrix(samples, args):
    # Create a list of all mutated genes
    # This could be overwritten in the future
    # by the genelist-argument
    if args.genelist == None:
        genes = []
        for sample in samples:
            for variant in sample.variants:
                if variant.gene not in genes:
                    genes.append(variant.gene)
    else:
        genes = parse_genes(args.genelist)
    
    # Write csv from variant_genes_by_sample dictionary
    with open(args.output + "_gene_matrix.csv", 'w') as ofile:
        writer = csv.writer(ofile, delimiter = ',')
        writer.writerow(["Sample"] + genes + ["Group"])
        for sample in samples:
            # Init empty row for each sample
            row = []
            for gene in genes:
                found = False
                types = []
                for variant in sample.variants:
                    if gene == variant.gene:
                        print(gene,sample.id)
                        type = check_variant_type(variant)
                        if type not in types:
                            types.append(type)
                        found = True
                if found == False:
                    row.append('-')
                else:
                    row.append(';'.join(types))
            writer.writerow([sample.id] + row + [sample.group])
    return(args.output + "_gene_matrix.csv")


def check_variant_type(var):
    type = var.type
    # Same as missense
    if type == "nonsynonymous_SNV":
        return "m"
    # Frameshift
    elif type == "frameshift_deletion" or type == "frameshift_insertion":
        return "f"
    elif type == "nonframeshift_deletion":
        return "d"
    # Stopgain is same as nonsense variant
    elif type == "stopgain":
        return "n"
    elif type == "splicing":
        return "s"
    else:
        return '-'

def parse_genes(genelist):
    genes = []
    with open(genelist,'r') as g:
        for line in g:
            line = line.strip().split('\t')
            genes.append(line[0])
    return genes

def flatten(statlist):
    flatlist = []
    for s in statlist:
        if isinstance(s,list): flatlist.extend(flatten(s))
        else: flatlist.append(s)
    return flatlist

def plot_substitution_frequencies(filename):
    path_to_sub_freq = current_path + '/RScripts/substitution_frequencies.R'
    cmd = [command, path_to_sub_freq, filename]
    subprocess.call(cmd)

def plot_gene_matrix(filename, var_n, args):
    path_to_sub_freq = current_path + '/RScripts/gene_matrix.R'
    cmd = [command, path_to_sub_freq, filename, args.genelist, var_n, args.rename, args.zero_variants]
    subprocess.call(cmd)

def prepare_base_substitutions(samples, group_names, output):
    with open(output + "_substitutions.csv", 'w') as output_file:
        writer = csv.writer(output_file, delimiter = ',')
        writer.writerow(["Amount", "Group", "Mutation"])
        for group in group_names:
            variants = []
            for sample in samples:
                if sample.group == group:
                    for var in sample.variants:
                        sub = parse_base_substution(var)
                        variants.append(sub)
            dist = Counter(variants)
            for i in dist.items():
                row = [i[1], group, i[0]]
                writer.writerow(row)
    return (output + "_substitutions.csv")

def parse_base_substution(variant):
    print(variant.ref, variant.alt)
    if (variant.ref != 'A'  and variant.ref != 'C' and variant.ref != 'G' and variant.ref != 'T') \
        or (variant.alt != 'A'  and variant.alt != 'C' and variant.alt != 'G' and variant.alt != 'T'):
        # Indel
        return 'indel'
    else:
        sub = str(variant.ref + "-" + variant.alt)
        sub = correct_substitution_type(sub)
        return sub

def correct_substitution_type(sub):
    if sub == "C-A" or sub == "G-T":
        return "C-A"
    elif sub == "C-G" or sub == "G-C":
        return "C-G"
    elif sub == "C-T" or sub == "G-A":
        return "C-T"
    elif sub == "A-C" or sub == "T-G":
        return "A-C"
    elif sub == "A-G" or sub == "T-C":
        return "A-G"
    elif sub == "A-T" or sub == "T-A":
        return "A-T"
    else:
        return sub

def amount_of_variants(samples, args):
    with open(args.output + "_amount_of_variants.csv", 'w') as output_file:
        writer = csv.writer(output_file, delimiter = ',')
        writer.writerow(["Sample","Amount_of_variants"])
        for sample in samples:
            writer.writerow([sample.id,len(sample.variants)])
    return args.output + "_amount_of_variants.csv"


def print_all_variants(sample_id, samples):
    for sample in samples:
        if sample.id.strip() == sample_id.strip():
            print(sample.id)
            print(sample.group)
            for variant in sample.variants:
                print(variant.gene)

def run():
    parser=argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest = "input", 
                        help = "Input file. Can be either a list of vcf-files in "
                        "a certain format or a table (See table details in the "
                        "documentation).")
    parser.add_argument("-o", "--output", dest = "output", 
                        help = "Filename base for plot and other output files.")
    parser.add_argument("-g","--groups", dest = "groups", 
                        help = "If table is provided then group file is needed. "
                        "Tab-separated file which has sample and group.")
    parser.add_argument("-gl", "--genelist", dest = "genelist",
                        help = "Script will visualize variants on the given genes by samples. \
                        If you want to group genes by pathway etc. give the gene group in the same \
                        row as the gene and separated by tab.")
    parser.add_argument("--remove_zero_variants", dest = "zero_variants", action = 'store_const', const="TRUE")
    parser.add_argument("--rename_samples", dest = "rename", action = 'store_const', const="TRUE")
    args = parser.parse_args()
    main(args)

if __name__ == "__main__":
    run()