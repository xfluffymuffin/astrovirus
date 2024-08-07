from Bio import SeqIO


def fix_rdp_output(path_old, path_new):
    mod = []
    # Названия записей в ORF_1A_1B_2_conc_DNA_no_amb.fas и ORF_1B_2_conc_DNA_no_amb.fas совпадают, поэтому как эталон
    # для называния заголовков выравниваний можно использовать любой из этих файлов
    with open(r"ORF_1A_1B_2_conc_RNA_no_amb.fas", 'r') as file:
        for line in file.readlines():
            if ">" in line:
                cut = line[1:].strip('\n')
                mod.append(cut)

        for record in SeqIO.parse(path_old, "fasta"):
            for x in mod:
                if record.id.startswith(x.split('_')[0]):
                    record.id = x
                    with open(path_new, "a") as new:
                        new.write(">" + str(record.id) + "\n" +
                                        str(record.seq.upper()) + "\n")

# Названия файлов были впоследствии возвращены к исходным 
fix_rdp_output(r"ORF_1A_1B_2_serotype_1.fas",
               r"ORF_1A_1B_2_serotype_1_desc_fix.fas")
