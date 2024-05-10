from Bio import SeqIO

def data_editor(path_old, path_new):
    ser = {"1": ["FJ755404", "KF211475", "AY720892", "JN887820", "FJ755403", "FJ755402", "FJ755405", "KC342249",
                 "ON571594", "NC_030922", "OQ630461", "LC694996", "LC694997", "MW485040", "MW485039", "MW485038",
                 "MW485041", "L23513", "MK059949", "Z25771", "HQ398856", "MT832896", "MT832892", "MT832895", "MT832894",
                 "MT832893", "MT832897", "MN433703", "MN433704", "LC694990", "LC694987", "LC694992", "LC694993",
                 "LC694994", "LC694986", "LC694995", "ON571596", "ON571597", "ON571597", "LC694988", "LC694991", "LC694989"],
           "2": ["KF039910", "MK059950", "MN433705", "MW485042", "KF039911", "KC285152"],
           "3": ["MT267483", "MK059951", "GU223905", "GU732187", "JF491430", "MW485043", "MZ603079", "MG571777",
                 "MK296753", "LC732124", "MN444721", "LC694985", "LC732122", "LC732123", "OP432297", "OP432300",
                 "OP432298", "OP432294", "OP432296", "OP432301"],
           "4": ["DQ070852", "AY720891", "MK059952", "KC285113", "KF039912", "MT906853", "MT906854", "MT906855",
                 "KF039913", "DQ344027", "OP432293", "OP432299", "OR371570", "OP432295", "ON571598"],
           "5": ["DQ028633", "MK059953", "JQ403108", "MT267476", "MN433706", "MW485044", "MW485045", "MF684776",
                 "ON571595", "OR222910", "MT906858", "MT906859", "MT906856", "MT906857"],
           "6": ["HM237363", "MK059954", "GQ495608", "GQ901902"], "8": ["KP862744", "AF260508", "MK059956"]}

    for record in SeqIO.parse(path_old, "fasta"):
        for key, values in ser.items():
            for id in values:
                with open(path_new, "a") as new:
                    if key == "1" and id in record.id:
                        record.id = record.id + "_serotype_1"
                        new.write(">" + str(record.id) + "\n" +
                                  str(record.seq.upper()) + "\n")
                    if key == "2" and id in record.id:
                        record.id = record.id + "_serotype_2"
                        new.write(">" + str(record.id) + "\n" +
                                  str(record.seq.upper()) + "\n")
                    if key == "3" and id in record.id:
                        record.id = record.id + "_serotype_3"
                        new.write(">" + str(record.id) + "\n" +
                                  str(record.seq.upper()) + "\n")
                    if key == "4" and id in record.id:
                        record.id = record.id + "_serotype_4"
                        new.write(">" + str(record.id) + "\n" +
                                  str(record.seq.upper()) + "\n")
                    if key == "5" and id in record.id:
                        record.id = record.id + "_serotype_5"
                        new.write(">" + str(record.id) + "\n" +
                                  str(record.seq.upper()) + "\n")
                    if key == "6" and id in record.id:
                        record.id = record.id + "_serotype_6"
                        new.write(">" + str(record.id) + "\n" +
                                  str(record.seq.upper()) + "\n")
                    if record.id == "MK059955_human_USA_NA_Oxford_NA_cell-culture":
                        record.id = record.id + "_serotype_7"
                        new.write(">" + str(record.id) + "\n" +
                                  str(record.seq.upper()) + "\n")
                    if key == "8" and id in record.id:
                        record.id = record.id + "_serotype_8"
                        new.write(">" + str(record.id) + "\n" +
                                  str(record.seq.upper()) + "\n")

data_editor(r"D:\PyCharm\BioPython\final_alignments\ORF2_cleared_alignment.fas",
            r"D:\PyCharm\BioPython\add_serotypes_data\ORF2_cleared_alignment_serotypes.fas")
