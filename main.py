from Bio import Entrez, SeqIO
import pandas as pd
from collections import Counter

# Функция для скачивания и сохранения полных GenBank-записей в файл records.gb
# и частичной записи данных в словарь res
def get_gb_records(query):
    # Загружаем записи из GenBank по запросу из переменной 'query'
    print('Загрузка записей...')
    Entrez.email = 'yakrit2013@yandex.ru'

    # Поиск в базе NCBI Nucleotide с учетом длины последовательности (полноты генома),
    # сохранение результатов в record.
    handle = Entrez.esearch(db='nucleotide',
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
        with open("records.gb", "a") as outfile:
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
                if key == 'GBID':
                    res[key].append(*rec.annotations["accessions"])
                else:
                    res[key].append(*rec.features[0].qualifiers.get(key))
            except TypeError:
                res[key].append('NaN')

# Функция для подсчета и вывода встречаемости каждого квалификатора
def count_occurences(qualifier):
    return [x[0] + ': ' + str(x[1]) for x in dict(Counter(res[qualifier])).items() if x[0] != 'NaN']

# Словарь для упорядоченного хранения необходимых данных из записей
res = {'GBID': [], 'country': [], 'host': [], 'collection_date': [], 'strain': [], 'isolate': [],
       'isolation_source': []}

get_gb_records('"Mamastrovirus 1"[porgn] AND "6000"[SLEN] : "8000"[SLEN]')

# Сохраняем данные из словаря res в датафрейм pandas и файл metadata.csv
data = pd.DataFrame.from_dict(res)
data.to_csv('metadata.csv', sep=',', index=False)

# Выводим часть данных в файл main_output.txt
with open('main_output.txt', 'w') as outfile:
    print(f'Загружено {data.GBID.count()} записей',
      f'\nОбразцы выделены в следующих точках мира (чаще всего '
      f'{(count := Counter(res["country"]).most_common()[0])[0]}: {count[1]} раз(а)):',
      *count_occurences('country'),
      '\nВиды-хозяева вирусов: ',  *count_occurences('host'),
      f'\nДаты взятия образцов (известны у '
      f'{data.collection_date.count() - data.collection_date.str.contains("N/A").sum()} '
      f'из {data.collection_date.count()})): ',
      *count_occurences('collection_date'),
      '\nШтаммы: ', *count_occurences('strain'),
      '\nИзвестная информация о серотипах/пациентах: ',
      *count_occurences('isolate'),
      '\nИсточники выделения: ',
      *count_occurences('isolation_source'),
      '\nБолее полное описание - в файле output.gb',
      sep='\n', file=outfile)