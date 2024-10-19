from Bio import SeqIO
import argparse


def cut_seq(input_file, output_file, cut_left, cut_right):

    for record in SeqIO.parse(input_file, 'fasta'):
        with open(output_file, 'a') as inp:
            inp.write(">" + str(record.id) + '\n' +
                      str(record.seq[cut_left - 1:cut_right].upper()) + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str,
                        help="Path to file with records/alignment (fasta)", required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Saving directory.")
    parser.add_argument("-l", "--left", type=int, default=1,
                        help="Left cutting point. Beginning of sequence(s) by default")
    parser.add_argument("-r", "--right", type=int,
                        help="Right cutting point")

    args = parser.parse_args()

cut_seq(args.input, args.output, args.left, args.right)
# cut_seq(r"C:\Users\gdvov\OneDrive\Desktop\useless.fas", "test.fas", 2, 5)
