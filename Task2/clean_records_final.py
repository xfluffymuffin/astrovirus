from Bio import SeqIO
import os

# Функция пробегает по записям fasta-файла и сначала формирует из них список, а затем достает только те, которые
# не содержат слова "sewage", т.е. убираем образцы, взятые из сточных вод
def clean_records(fasta_file):
    seqs_cleaned = []
    for rec_num, seq in enumerate(SeqIO.parse(fasta_file, "fasta")):
        seqs_cleaned.append(seq)

    # Создание очищенного файла с форматированием в 60 символов на строку
    with open("D:/PyCharm/BioPython/intermediate.fasta", 'w') as file:
        SeqIO.write(seqs_cleaned,
                    file,
                    'fasta')
    SeqIO.write(SeqIO.parse("D:/PyCharm/BioPython/intermediate.fasta", "fasta"),
                'D:/PyCharm/BioPython/Mamastrovirus_1_complete_genome_records_CLEARED.fasta',
                'fasta')
    os.remove("D:/PyCharm/BioPython/intermediate.fasta")

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


clean_records("D:/PyCharm/BioPython/Mamastrovirus_1_complete_genome_records.fasta")
check_collection_date("Mamastrovirus_1_complete_genome_records_CLEARED.fasta", "records_no_collection_date.txt")
# 1 - очищенный файл без последовательностей сточных вод, 2 - выходной файл с последовательностями без даты выделения образца
