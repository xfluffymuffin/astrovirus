from Bio import SeqIO
import re
import pandas as pd
import argparse
import os


def ref_del_gaps(inp_aln, ref_name, save_dir):
    for record in SeqIO.parse(inp_aln, 'fasta'):
        if record.id == ref_name:
            with open(f'{save_dir}/{ref_name}_no_gaps.fas', 'w') as ref_seq:
                ref_seq.write(f'>{ref_name}_no_gaps\n{''.join(str(record.seq).split('-')).upper()}')


def parse_ref_coords(inp_ref_coords_csv):
    df_ref_coords = pd.read_csv(inp_ref_coords_csv)
    for index, row in df_ref_coords.iterrows():
        find_aln_region_by_ref_coords(args.reference, row.iloc[0], row.iloc[1], args.input_aln, args.output)


def find_aln_region_by_ref_coords(ref_name, left_border, right_border, inp_aln, save_dir):
    global df_aln_coords
    for ref_fasta in SeqIO.parse(f'{save_dir}/{ref_name}_no_gaps.fas', 'fasta'):
        subseq = str(ref_fasta.seq)[int(left_border) - 1:int(right_border)].upper()
        for ref in SeqIO.parse(inp_aln, 'fasta'):
            if ref.id == ref_name:
                aln_ref = re.findall(f"[{''.join('-')}]*".join(subseq), str(ref.seq))
                new_left_border = f"{str(ref.seq).index(aln_ref[0]) + 1}"
                new_right_border = f"{str(ref.seq).index(aln_ref[0]) + len(aln_ref[0])}"

                print(f"Reference {ref_name}: {left_border}-{right_border}, "
                      f"alignment: {new_left_border}-{new_right_border}")

                df_aln_coords = df_aln_coords._append({'new_left_border': new_left_border,
                                                       'new_right_border': new_right_border},
                                                      ignore_index=True)


def merge_csv(inp_ref_coords_csv, ref_name, save_dir):
    merged_csv = pd.concat((pd.read_csv(inp_ref_coords_csv), df_aln_coords), axis=1)
    merged_csv.to_csv(f'{save_dir}/{ref_name}_ref_aln_coords.csv', index=False)

    # Additionally saving aln coords only
    df_aln_coords.to_csv(f'{save_dir}/{ref_name}_aln_coords.csv', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-ia", "--input_aln", type=str,
                        help="Path to fasta file with alignment", required=True),
    parser.add_argument("-ic", "--input_coord", type=str,
                        help="Path to csv file with coordinates. Must contain two columns with"
                             "left and right borders of the genomic regions you'd like to cut", required=True),
    parser.add_argument("-ref", "--reference", type=str,
                        help="Full name of the reference sequence from your alignment")
    parser.add_argument("-o", "--output", type=str,
                        help="Saving directory. Folder created if doesn't exist. Script directory unless specified")

    args = parser.parse_args()

df_aln_coords = pd.DataFrame(columns=['new_left_border', 'new_right_border'])

if not args.output:
    args.output = os.path.dirname(__file__)

if not os.path.exists(args.output):
    os.makedirs(args.output)

ref_del_gaps(args.input_aln, args.reference, args.output)
parse_ref_coords(args.input_coord)
merge_csv(args.input_coord, args.reference, args.output)
