import pandas as pd
from pandas import DataFrame
from Bio import Entrez, SeqIO
from collections import Counter

# Счетчик записей, для которых известна дата выделения
count = 0

# Словарь для упорядоченного хранения необходимых данных из записей
res = {'country': [], 'host': [], 'collection_date': [], 'strain': [], 'isolate': []}

Entrez.email = 'yakrit2013@yandex.ru'

# 1) Поиск в базе NCBI Nucleotide с учетом длины последовательности (полноты генома), сохранение результатов в record.
# Параметр retmax ограничивает максимальное количество результатов поиска
handle = Entrez.esearch(db='nucleotide',
                        term='"Mamastrovirus 1"[porgn] AND complete[All Fields] '
                             'AND genome[All Fields] AND "6000"[SLEN] : "8000"[SLEN]',
                        retmax='1000000')
record = Entrez.read(handle)

# 2) Из record берем id результатов поиска и для каждого из них загружаем данные в формате GenBank, сохраняем в rec
for i in record["IdList"]:
    download = Entrez.efetch(db="nucleotide", id=i, rettype="gb", retmode="text")
    rec = SeqIO.read(download, "genbank")

    # Завершаем загрузку
    download.close()

    # Здесь пригождаются названия ключей в словаре res, совпадающие с извлекаемыми из записей данными.
    # Например, если ключ 'country', то туда записываем исключительно информацию по странам, где выделяли образцы
    # Если в конкретной записи отсутствует необходимое значение, игнорируем ошибку (конструкция try/except)
    for key in res.keys():
        try:
            res[key].append(*rec.features[0].qualifiers.get(key))
        except TypeError:
            pass

    # 3) Создаем датафрейм с полученной из записей информацией (например, в 'Keys' будет country, а в 'Values' - China)
    df = DataFrame({'Keys': rec.features[0].qualifiers.keys(),
                    'Values': rec.features[0].qualifiers.values()})
    # print(df)

    # Считаем кол-во записей, для которых известна дата выделения
    count += len(df[df['Keys'] == 'collection_date'])

    # 4) Запись датафрейма в файл Excel. У меня df.to_csv игнорирует попытки передать ему какой-либо sep, и
    # в итоге в Excel наблюдается каша в одном столбце. С df.to_excel тоже есть проблема: if_sheet_exists='overlay'
    # не дает дополнять датафреймы на одном листе (получается информация только для одной записи,
    # т.е. все перезаписывается). Приходится добавлять каждый датафрейм на новый лист
    with pd.ExcelWriter('output.xlsx',
                        mode='a', if_sheet_exists='new') as writer:
        df.to_excel(writer, sheet_name="Main", index=False)

    # 5) В поле 'strain' почему-то может попадать информация как о серотипе, так и о пациенте
print(f'\nОбразцы выделены в следующих точках мира, чаще всего {Counter(res["country"]).most_common()[0]}:',
      *res["country"],
      '\nВиды-хозяева вирусов: ', *res['host'],
      f'\nДаты взятия образцов (известны у {count} из {len(record["IdList"])}): ', *res['collection_date'],
      '\nШтаммы: ', *res['strain'],
      '\nИнформация о серотипе или пациенте (?): ', *res['isolate'], sep='\n')
