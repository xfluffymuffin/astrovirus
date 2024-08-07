import argparse
import pandas as pd
import os
import re
from Bio import SeqIO
from matplotlib import pyplot as plt


def cds_distr(input_file, output_dir, out_name):
    '''
    Parses GenBank file, prints sorted list of annotations found in 
    product and gene qualifiers in CDS field
    '''
    out_file = open(os.path.join(output_dir, out_name), 'w')
    
    with open(input_file) as handle:
        annot = []
        records = list(SeqIO.parse(handle, 'gb'))
        rec_annot_cds = 0
        cds_prod_count = 0
        cds_gene_count = 0
        for rec in records:
            #print(rec.id)
            rec_flag = 0
            for feature in rec.features:
                if feature.type == 'CDS':
                    if 'product' in feature.qualifiers.keys():
                        cds_prod_count += 1
                        rec_flag = 1
                        #print('product ' + feature.qualifiers['product'][0])
                        annot.append(feature.qualifiers['product'][0])
                    if 'gene' in feature.qualifiers.keys():
                        cds_gene_count += 1
                        rec_flag = 1
                        #print('gene ' + feature.qualifiers['gene'][0])
                        annot.append(feature.qualifiers['gene'][0])
            rec_annot_cds += rec_flag
        out_file.write('Number of records in file: {}\n'.format(len(records)))
        out_file.write('Records with annotated CDS: {}\n'.format(rec_annot_cds))
        out_file.write('CDS with product fields {}\n'.format(cds_prod_count))
        out_file.write('CDS with gene fields {}\n'.format(cds_gene_count))
        
        print('Number of records in file: {}'.format(len(records)))
        print('Records with annotated CDS: {}'.format(rec_annot_cds))
        print('CDS with product fields {}'.format(cds_prod_count))
        print('CDS with gene fields {}'.format(cds_gene_count))
        #print('Annotations:')
        annot = list(set(annot))
        annot.sort()
        for a in annot:
            out_file.write(a + '\n')
            #print(a)
        



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file in GenBank format",
                        required=True)
    parser.add_argument("-odir", "--output_dir", type=str,
                        help="Output directory to save the output file",
                        )
    parser.add_argument("-oname", "--output_name", type=str,
                        help="Name of output file",
                        required=True)
    args = parser.parse_args()
    
    if not args.output_dir:
            args.output_dir = ''

    cds_distr(args.input_file, args.output_dir, args.output_name)
