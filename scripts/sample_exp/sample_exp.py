import os
import shutil
import argparse
import Bio
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import time
from time import localtime, strftime


def create_query(input_file, cut_left, cut_right):
    try:
        os.remove('query.fa')
    except OSError:
        pass

    # Save records from the input file to another fasta file with cuts
    for record in SeqIO.parse(input_file, 'fasta'):
        with open('query.fa', 'a') as query:
            query.write(">" + str(record.id) + '\n' +
                        str(record.seq[cut_left - 1:cut_right].upper()) + '\n')


def launch_blast():
    # Perform the alignment
    res_handle = Bio.Blast.NCBIWWW.qblast("blastn", "nt",
                                          open('query.fa').read(),
                                          ungapped_alignment=True,
                                          megablast=True)

    with open('blast_results.xml', 'w') as save_file:
        blast_res = res_handle.read()
        save_file.write(blast_res)


def show_save_blast_results(cut_left, cut_right):
    if cut_right is None:
        cut_right = 'end'

    e_val_thr = 1e-20

    # Adding queries and matches to a new manually created file - a workaround
    # (No pandas because the resulting dataframe somehow had 49 matches for every query. There was
    # no error saying 'all arrays must be of the same length' which meant the dictionary I used to
    # create the dataframe was broken)

    try:
        os.remove('blast_query_match.txt')
    except OSError:
        pass

    # Queries and matches are separated by "//", matches themselves are also separated by ";"

    with open('blast_query_match.txt', 'a') as blast_txt:
        blast_txt.write(f"query_{cut_left}_to_{cut_right}//matches")

        for record in NCBIXML.parse(open("blast_results.xml")):
            if record.alignments:
                print("\n")
            print(f"query: {record.query[:150]}, positions: {cut_left} to {cut_right}")
            blast_txt.write(f"\n{record.query}//")

            for aln in record.alignments:
                for hsp in aln.hsps:
                    if hsp.expect < e_val_thr:
                        print(f"match: {aln.title[:150]}")
                        blast_txt.write(f"{aln.title};")


def move_output_files(new_dir_name):
    if new_dir_name is None:
        new_dir_name = strftime("%Y%m%d_%H%M", localtime())

    if not os.path.exists(new_dir_name):
        os.mkdir(f"{new_dir_name}")
    try:
        shutil.move("blast_query_match.txt", f"{new_dir_name}\\blast_query_match.txt")
        shutil.move("query.fa", f"{new_dir_name}\\query.fa")
        shutil.move("blast_results.xml", f"{new_dir_name}\\blast_results.xml")
    except FileNotFoundError:
        pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str,
                        help="Path to file with records/alignment (fasta)", required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Saving directory. Folder named by current time by default")
    parser.add_argument("-l", "--left", type=int, default=1,
                        help="Left cutting point. Beginning of sequence(s) by default")
    parser.add_argument("-r", "--right", type=int,
                        help="Right cutting point")
    parser.add_argument("-si", "--save_intermediate", default=False,
                        help="Save intermediate files (blast query w/gaps "
                             "(in case of alignment) and results in *.xml format), "
                             "False by default; True or '1' to activate")

    args = parser.parse_args()

create_query(args.input, args.left, args.right)
launch_blast()
show_save_blast_results(args.left, args.right)

if not args.save_intermediate:
    os.remove('query.fa')
    os.remove('blast_results.xml')

move_output_files(args.output)
