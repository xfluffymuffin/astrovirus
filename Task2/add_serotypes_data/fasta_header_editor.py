from Bio import SeqIO
import pandas as pd
import numpy as np

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)


def edit_fasta(path_old, path_new):
    # Файл с аннотацией
    ser = pd.read_csv(r"D:\PyCharm\BioPython\tree_DNA_alignment\iTOL\ORF2_DNA_tree_ann.txt",
                      sep="\t",
                      header=None,
                      skiprows=3) \
          .drop(columns=[1, 3, 4]) \
          .rename({0: 'leaf_name', 2: 'color'}, axis=1) \
          .head(104)
    # Первые строки с названиями записей содержат оригинальные названия (все остальное - сочетания как в дереве)

    # Посмотреть уникальные значения RGB цветов
    # print(ser.color.unique())

    conditions = [ser['color'] == 'rgb(191, 144, 0)',  # yellow
                  ser['color'] == 'rgb(116, 27, 71)',  # magenta
                  ser['color'] == 'rgb(0, 211, 32)',  # green
                  ser['color'] == 'rgb(204, 0, 0)',  # red
                  ser['color'] == 'rgb(255, 126, 193)',  # pink
                  ser['color'] == 'rgb(0, 0, 255)',  # blue
                  ser['color'] == 'rgb(111, 168, 220)',  # light blue
                  ser['color'] == 'rgb(180, 95, 6)']  # brown

    choices = ['serotype_8', 'serotype_4', 'serotype_6', 'serotype_5',
               'serotype_2', 'serotype_7', 'serotype_3', 'serotype_1']  # same order

    ser['serotype'] = np.select(conditions, choices)
    ser['full_name'] = ser['leaf_name'] + "_" + ser['serotype']

    mod = [str(x[0]) + '_serotype_' + str(x[-1]) for x in list(ser.full_name.str.split('_'))]
    # print(mod)

    for record in SeqIO.parse(path_old, "fasta"):
        for x in mod:
            with open(path_new, "a") as new:
                for i in range(9):
                    if x.endswith(f"_{i}") and x.split('_')[0] in record.id:
                        record.id = record.id + f"_serotype_{i}"
                        new.write(">" + str(record.id) + "\n" +
                                  str(record.seq.upper()) + "\n")


edit_fasta(r"D:\PyCharm\BioPython\final_alignments\ORF2_cleared_alignment.fas",
           r"D:\PyCharm\BioPython\add_serotypes_data\ORF2_cleared_alignment_serotypes.fas")
