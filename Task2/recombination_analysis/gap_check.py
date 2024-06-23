from Bio import SeqIO
import pandas as pd
import argparse
import os


def check_gaps(old_path, new_path_csv, new_path_fasta, threshold, exception):
    seqs, names, res_len, clean = [], [], [], [] #
    # seqs - contains ALL records from the input file
    # names - contains names of ALL sequences
    # res_dict - contains lengths of cont. gaps
    # clean - contains names of sequences that are clean enough to not get removed, such seqs are written to a fasta file
    count = 0

    for record in SeqIO.parse(old_path, 'fasta'):
        seqs.append(record.seq)
        names.append(record.id)

    for record in range(len(seqs)):
        inter = []
        for letter in seqs[record]:
            if letter == "-":
                count += 1
            elif letter in "ATGC":
                inter.append(count)
                count = 0

        res_len.append([x for x in inter if x != 0 and x != exception])
    df = pd.DataFrame({'Names': names, "Count": res_len})
    df.to_csv(new_path_csv, index=False)

    print('\nFollowing sequences exceeded threshold: \n '
          '(sequence name --> number of continuous gaps)\n')
    for index, row in df.iterrows():
        if max(row['Count']) >= threshold:
            print(row['Names'] + " -->", max(row['Count']))
        else:
            clean.append(row['Names'])

    for record in SeqIO.parse(old_path, "fasta"):
        for x in clean:
            if record.id.startswith(x.split('_')[0]):
                record.id = x
                with open(new_path_fasta, "a") as new:
                    new.write(">" + str(record.id) + "\n" +
                              str(record.seq.upper()) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str,
                        help="Path to file with alignment (fasta)", required=True)
    parser.add_argument("-oc", "--output_csv", type=str,
                        help="Path to the output *.csv file - directory of input file by default")
    parser.add_argument("-of", "--output_fasta", type=str,
                        help="Path to the output *.fasta file containing only records not "
                             "exceeding the threshold - directory of input file by default")
    parser.add_argument("-t", "--threshold", type=int, default=50,
                        help="Sequences exceeding the threshold (max number of continuous gaps in a sequence) "
                             "will be listed on the screen")
    parser.add_argument("-e", "--exception", type=int, default=0,
                        help="If you have a certain gap length to exclude from analysis, this is your opportunity")

    args = parser.parse_args()

    args.input = os.path.realpath(args.input)
    converted_path = args.input.replace('\\', "/")
    cut_input_path = converted_path.split('/')[-1].split('.fa')[0]
    if not args.output_csv:
        args.output_csv = os.path.split(args.input)[0] + rf'\{cut_input_path}_gap_info.csv'
    if not args.output_fasta:
        args.output_fasta = os.path.split(args.input)[0] + rf'\{cut_input_path}_fewer_gaps.fas'

    check_gaps(args.input, args.output_csv, args.output_fasta, args.threshold, args.exception)
