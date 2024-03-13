import os.path
from Bio import Entrez, SeqIO
import pandas as pd
from collections import Counter
import pickle


# Функция для скачивания и сохранения полных GenBank-записей в файл records.gb
# и частичной записи данных в словарь res
def get_gb_records(query):
    # Загружаем записи из GenBank по запросу из переменной 'query'
    print("Загрузка записей...")
    Entrez.email = "yakrit2013@yandex.ru"

    # Поиск в базе NCBI Nucleotide с учетом длины последовательности (полноты генома),
    # сохранение результатов в record.
    handle = Entrez.esearch(db="nucleotide",
                            term=query,
                            retmax=10000000)
    record = Entrez.read(handle)

    # Из record берем id результатов поиска и для каждого из них загружаем данные в формате GenBank, сохраняем в rec.
    for i in record["IdList"]:
        download_for_writing = Entrez.efetch(db="nucleotide",
                                             id=i,
                                             rettype="gb",
                                             retmode="text")
        # Создаем выходной файл с записями GenBank – records.gb
        with open("Mamastrovirus_1_complete_genome_records.gb", "a") as outfile:
            outfile.write(download_for_writing.read())

        download_for_parsing = Entrez.efetch(db="nucleotide",
                                             id=i,
                                             rettype="gb",
                                             retmode="text")

        # Загруженные данные также частично сохраняем в словарь res.
        # Если информация по какому-либо квалификатору отсутствует, пишем NaN
        rec = SeqIO.read(download_for_parsing, "genbank")
        for key in res.keys():
            try:
                if key == "GBID":
                    res[key].append(*rec.annotations["accessions"])
                else:
                    res[key].append(*rec.features[0].qualifiers.get(key))
            except TypeError:
                res[key].append("NaN")

        with open("dict_data.txt", "wb") as save_res:
            pickle.dump(res, save_res)


# Функция для подсчета и вывода встречаемости каждого квалификатора
def count_occurrences(qualifier):
    return [x[0] + ': ' + str(x[1]) for x in dict(Counter(res[qualifier]).most_common()).items()]


def file_writer():
    # Сохраняем данные из словаря res в датафрейм pandas и файл metadata.csv
    data = pd.DataFrame.from_dict(res)
    data.to_csv("metadata.csv", sep=",", index=False)

    # Выводим часть данных в файл main_output.txt
    with open("main_output.txt", "w") as outfile:
        print(f'Загружено {data.GBID.count()} записей',
              f'\nОбразцы выделены в следующих точках мира (чаще всего '
              f'{(count := Counter(res["country"]).most_common()[0])[0]}: {count[1]} раз(а)):',
              *count_occurrences('country'),
              '\nВиды-хозяева вирусов: ', *count_occurrences('host'),
              '\nДаты взятия образцов: ',
              *count_occurrences('collection_date'),
              '\nШтаммы: ', *count_occurrences('strain'),
              '\nИзвестная информация о серотипах/пациентах: ',
              *count_occurrences('isolate'),
              '\nИсточники выделения: ',
              *count_occurrences('isolation_source'),
              sep='\n', file=outfile)


# Словарь для упорядоченного хранения необходимых данных из записей
res = {"GBID": [], "country": [], "host": [], "collection_date": [], "strain": [], "isolate": [],
       "isolation_source": []}

# Если dict_data.txt существует, то скрипт запускался ранее, данные снова можно не загружать.
# Чтобы заставить программу скачивать данные и сохранять их в Mamastrovirus_1_complete_genome_records.gb заново,
# необходимо удалить файл dict_data.txt
try:
    with open("dict_data.txt", "rb") as save_res:
        res = pickle.load(save_res)
    file_writer()
except FileNotFoundError:
    # Необходимо во избежание добавления GenBank записей в конец .gb файла при удалении dict_data.txt
    # Если .gb файл существует, удаляем: get_gb_records создаст его заново
    if os.path.exists("Mamastrovirus_1_complete_genome_records.gb"):
        os.remove("Mamastrovirus_1_complete_genome_records.gb")
    else:
        pass
    get_gb_records('"Mamastrovirus 1"[porgn] AND "6000"[SLEN] : "8000"[SLEN]')
    file_writer()
