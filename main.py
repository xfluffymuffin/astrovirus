from Bio import Entrez, SeqIO
import pandas as pd
from collections import Counter


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

    # Из record берем id результатов поиска и для каждого из них загружаем данные в формате GenBank, сохраняем в rec
    for i in record["IdList"]:
        download = Entrez.efetch(db="nucleotide",
                                 id=i,
                                 rettype="gb",
                                 retmode="text")
        rec = SeqIO.read(download, "genbank")

        # Можно дополнительно выводить ID загружаемых записей
        # print(rec.annotations['accessions'])

        # Завершаем загрузку
        download.close()

        # Записываем все загруженное в файл output.txt. Если значение для какого-либо поля отсутствует, пишем N/A
        outfile = open('output.txt', 'a')
        outfile.write(f'\nGBID: {str(*rec.annotations["accessions"])}\n')
        for key in res.keys():
            try:
                if key == 'GBID':
                    res[key].append(*rec.annotations["accessions"])
                else:
                    res[key].append(*rec.features[0].qualifiers.get(key))
                    outfile.write(f'{key}: {str(*rec.features[0].qualifiers.get(key)).replace(": ", ":")}\n')
            except TypeError:
                res[key].append('N/A')
                outfile.write(f'{key}: N/A\n')
        outfile.close()


def parse():
    # Создаем датафрейм из записей в output.txt
    data = pd.read_csv('output.txt', sep=': ', header=None, engine='python')

    # Выводим важное на экран (дату и страну сбора, виды-хозяева, штаммы, серотипы, источники выделения, информация о
    # пациенте
    print(f'Загружено {len((total := data[data[0] == "collection_date"][1].to_list()))} записей'
          f'\nОбразцы выделены в следующих точках мира (чаще всего '
          f'{(countries := Counter(data[data[0] == "country"][1].to_list()).most_common()[0])[0]}, '
          f'{countries[1]} раз(а)):',
          *data[data[0] == 'country'][1].dropna().unique(),
          '\nВиды-хозяева вирусов: ', *data[data[0] == "host"][1].dropna().unique(),
          f'\nДаты взятия образцов (известны у '
          f'{len([x for x in total if len(str(x)) > 3])} '  # фильтр на длину даты, т.е. 'NaN' не выводятся
          f'из {len(total)})): ',
          *data[data[0] == "collection_date"][1].dropna().unique(),
          '\nШтаммы: ', *data[data[0] == "strain"][1].dropna().to_list(),
          '\nИзвестная информация о серотипах/пациентах: ', *data[data[0] == 'isolate'][1].dropna().unique(),
          '\nИсточники выделения: ', *data[data[0] == 'isolation_source'][1].dropna().unique(),
          '\nБолее полное описание - в файле output.txt',
          sep='\n')


# Словарь для упорядоченного хранения необходимых данных из записей
res = {'GBID': [], 'country': [], 'host': [], 'collection_date': [], 'strain': [], 'isolate': [],
       'isolation_source': []}

# Если файла output.txt не существует, то и записи отсутствуют. Тогда сначала скачиваем их. В противном случае работаем
# с готовым файлом
try:
    parse()
except FileNotFoundError:
    get_gb_records('"Mamastrovirus 1"[porgn] AND "6000"[SLEN] : "8000"[SLEN]')
    parse()
