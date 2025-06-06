from Bio import SeqIO
import os


# Функция пробегает по записям fasta-файла и сначала формирует из них список, а затем достает только те, которые
# не содержат слова "sewage", т.е. убираем образцы, взятые из сточных вод
def clean_records(fasta_file):
    seqs = []
    for record, seq in enumerate(SeqIO.parse(fasta_file, "fasta")):
        seqs.append(seq)
    with open("intermediate.fasta", 'w') as file:
        for i in seqs:
            if "sewage" not in i.id:
                # Записываем в промежуточный файл нужные записи
                file.write(">" + str(i.id) + "\n" +
                           str(i.seq.upper()) + "\n")


# Переписываем fasta-файл так, чтобы в каждой строке с последовательностью содержалось 60 символов
def fix_format(fasta):
    SeqIO.write(SeqIO.parse(fasta, "fasta"),
                'Mamastrovirus_1_complete_genome_records_CLEARED.fasta',
                'fasta')
    os.remove("intermediate.fasta")


# Проверяем, есть ли в заголовках оставшихся записей дата взятия образца
def check_collection_date(fasta, output):
    coll_date, all_rec = [], []
    with open(fasta, 'r') as file:
        for line in file:
            for i in range(1900, 2025):
                if ">" in line:
                    # Все заголовки
                    all_rec.append(line.strip(">\n"))

                if ">" in line and str(i) in line:
                    # Заголовки с датами
                    coll_date.append(line.strip(">\n"))

    with open(output, 'w') as writer:
        # Если даты нет, записываем весь заголовок в отдельный файл для дальнейшей работы
        for fasta_header in list(set([x for x in all_rec if
                                      x not in coll_date])):  # Превращение в set(no_coll_date) и обратно - для отсеивания повторяющихся значений списка
            writer.write(fasta_header + "\n")



clean_records("Mamastrovirus_1_complete_genome_records.fasta")
fix_format("intermediate.fasta")
check_collection_date("Mamastrovirus_1_complete_genome_records_CLEARED.fasta",
                      "records_no_collection_date.txt")
# В последней функции:
# 1 - файл с очищенными последовательностями, 2 - выходной файл с именами последовательностей без дат взятия образца
