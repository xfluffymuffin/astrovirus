## **aln/**

Содержит выравнивания различного происхождения:

* Подпапка main/: основные выравнивания, актуальные на данный момент
	* Подпапка tree_aln/: выравнивания, использовавшиеся для построения деревьев в iqtree, заканчиваются на 'no_amb'. Файлы iqtree (см. папку trees/iqtree_output/) впоследствии пропущены через iTOL для окрашивания по серотипам (см. папку trees/iTOL_output/)
		- Подпапка recomb_analysis_aln/: выравнивания, использовавшиеся для анализа событий рекомбинации - на них использован скрипт gap_check.py (папка scripts/). Включает папку serotypes_aln/: содержит выравнивания по отдельным серотипам, полученные в RDP4 (из выравниваний ORF_1A_1B_2_conc_RNA_no_amb.fas, ORF_1B_2_conc_RNA_no_amb.fas в папке aln/main/recomb_analysis_aln/)
		- Подпапка pair_dist_aln/: выравнивания для нахождения расстояния между последовательностями - результаты в папаке pairwise_distance/

	
*выравнивания в подпапках main/recomb_analysis_aln/ и main/pair_dist_aln/ частично совпадают - одинаковые файлы были использованы для решения разных задач

------------------------------------------------------------------------------------------------------------------------------------------

## **genbank/**

Содержит gb-файлы (записи с последовательностями) либо ассоциированные с ними файлы (названия записей без даты взятия образца)

------------------------------------------------------------------------------------------------------------------------------------------

## **genome_coverage/**

Содержит созданные скриптом genbank_coverage_complex.py файлы для задания по определению покрытия генома

reference.fas - референсные последовательности серотипов
reference_aln.fas - выравнивание для референсов
Final_cov.txt - данные для построения графика Final_cov_plot.png из подпапки plots/original
file.fasta - скачанные скриптом выравненные на референсы последовательности

* Подпапка plots/: содержит графики
	- Подпапка original/: исходные графики покрытия, включая графики для каждого серотипа в 			отдельности
	- Подпапка inkscape_modified/: графики (Final_plot_cov, а также breakpoint distribution plot и recombination region count matrix, полученные в RDP по выравниванию ORF_1A_1B_2_conc_RNA_no_amb.fas), переведенные в векторный формат и модифицированные с помощью Inkscape - данные с графиков сопоставлены с добавленными схемами генома

- *.svg-файл recombination_region_count_matrix.svg из папки genome_coverage/plots/inkscape_modified/ не уместился в лимит в 25 MB, поэтому на него есть ссылка в текстовом файле

------------------------------------------------------------------------------------------------------------------------------------------

## **add_info/**

Содержит вспомогательные файлы:

* gbid_serotype_corresp.csv - соответствия gb_id и серотпиов
* orf_data.csv - координаты рамок из gb-файла с записями

------------------------------------------------------------------------------------------------------------------------------------------

## **trees/**

Содержит выходные файлы после построения деревьев:

* Подпапка iqtree_output/: файлы iqtree, без визуализации; для построения деревьев использованы выравнивания из папки aln/main/tree_aln/
* Подпапка iTOL_output/: сгенерированные iTOL файлы - деревья в нескольких форматах (newick, phyloxml), аннотация деревьев
		- Подпапки pics/ для папок с деревьями: изображения деревьев в *.png и *.pdf
	
------------------------------------------------------------------------------------------------------------------------------------------

## **parwise_distance/**

Содержит полученные в MEGA v11 данные о расстояниях между белковыми и нуклеотидными последовательностями по двум последним и по всем трем рамкам (подпапки aa/, nucl/), а также соответствующие графики (подпапки aa_plot/, nucl_plot/)

Выравнивания для выполнения этого задания можно найти в папке aln/main/pair_dist_aln/.
В заголовках графиков имеется информация о выравниваниях, а также о выходных *.csv файлах.

------------------------------------------------------------------------------------------------------------------------------------------

## **RDP/**

Содержит проанализированные в RDP проверенные события рекомбинации:

Для выравнивания всей выборки
(исходные выравнивания - ORF_1A_1B_2_conc_RNA_less_amb_fewer_gaps.fas, ORF_1B_2_conc_RNA_frame_fixed_less_amb_fewer_gaps.fas из папки aln/main/recomb_analysis_aln/):

* Подпапка 1a_1b_2: по полногеномному выравниванию
* Подпапка 1b_2: по выравниванию двух последних рамок

Для выравнивания по серотипам
(исходные выравнивания - aln/main/recomb_analysis_aln/serotypes_aln/): 

* Подпапка 1a_1b_2_serotypes: по полногеномному выравниванию
* Подпапка 1b_2_serotypes: по выравниванию двух последних рамок

------------------------------------------------------------------------------------------------------------------------------------------

## **old/**

Содержит устаревшие файлы проекта, замененные обновленными (могут пригодиться больше, чем не перекочевавшие сюда из старого репозитория)

* Подпапка aln_old/: промежуточные, более не актуальные выравнивания; перемещаются сюда из подпапки aln/main/
* Подпапка rdp_old/: старые RDP-проекты и выходную информацию о них в формате *.csv
* Подпапка trees_old/: старые файлы, связанные с построением деревьев. Папка old_annot/ -  с устаревшей аннотацией

## **scripts/**

Содержит написанные для проекта скрипты вместе с входными и выходными файлами, даже если это дубликаты уже имеющихся файлов в других папках (так папку scripts можно использовать независимо от ее расположения)

### _genecds_create_csv.py_

Программа для создания csv-файла по записям из полей CDS, genes записей GenBank

Последовательно запускаются две функции:

1) launch_gene_cds_field_distribution()
Функция запускает скрипт gene_cds_field_distribution.py с помощью следующей внутренней команды:
python genecds_field_distribution.py -input Mamastrovirus_1_complete_genome_records.gb -odir alignment_of_orfs_output -oname Mamastrovirus_1_gene_product_output.txt
Выходной файл такой же, как и у запускаемого изнутри скрипта: Mamastrovirus_1_gene_product_output.txt

2) gene_cds_field_distribution_output_handle()
Функция непосредственно создает необходимый csv-файл на основании выходного файла предыдущей функции
Если строка поля CDS/gene – это слово orf1/orf1a, либо в строке имеется “unknown”, “protease”, “NS”, то после запятой добавляется “1A”.
Строка – “orf1b” или в строке “polymerase” – добавляется “1B”/
В строке “Non” – добавляется “1AB”.
В строке “non” – добавляется “1A”, “1B” и “1AB”.
Иначе добавляется “2”.
Скрипт genecds_field_distribution.py и файл с GenBank-записями Mamastrovirus_1_complete_genome_records.gb должны находиться в той же директории, что и сам genecds_create_csv.py.

					
### _parse_gb_records.py_

Программа для получения отдельных данных записей GenBank из NCBI

Использует четыре библиотеки: pandas, biopython, collections, pickle, а также модуль os.path.
- Вначале создается словарь res, который в будущем будет содержать часть данных из записей GenBank.
- Затем проверяется наличие файла dict_data.txt, содержащего данные словаря res. Файл создается после первого запуска программы и необходим, чтобы избегать повторного скачивания GenBank-записей при каждом новом запуске программы – удобно, если нужно, к примеру, поменять что-либо в формате вывода данных. Для повторного скачивания данных необходимо перед запуском скрипта удалить dict_data.txt.
- Далее при отсутствии dict_data.txt выполняется функция get_gb_records(query), где query – запрос, отправляемый на сервер NCBI (к базе данных NCBI Nucleotide). Функция выполняет три задачи: скачивание данных из GenBank, запись данных некоторых квалификаторов в словарь res и сохранение всех скачанных записей GenBank в файл records.gb.
Если dict_data.txt уже существует, то из него подгружаются данные GenBank записей, сохраненные ранее в словарь res, и функция get_gb_records игнорируется.
- Затем функцией file_writer() данные из словаря res сохраняются в датафрейм pandas, оттуда – в файл metadata.csv (ныне файл заменен аналогом metadata_serotypes_fixed.csv)
- Часть данных, сохраненная в словаре res, выводится в файл main_ output.txt.

Переменной query было передано следующее значение:
```
"Mamastrovirus 1"[porgn] AND "6000"[SLEN] : "8000"[SLEN]
```
Содержимое main_output.txt частично приведено ниже (результат скачивания данных по запросу):
```
Загружено 125 записей

Образцы выделены в следующих точках мира (чаще всего China: 19 раз(а)): 
China: 19
France: 15
Brazil: 15
...

Виды-хозяева вирусов: 
Homo sapiens: 84
NaN: 25
...

Даты взятия образцов: 
NaN: 17
2015: 10
2014: 5
...

Штаммы: 
NaN: 75
Oxford: 8
Dresden: 2
...

Известная информация о серотипах/пациентах: 
NaN: 35
Child with acute gastroenteritis: 3
V1347: 2
...

Источники выделения: 
NaN: 47
stool: 20
sewage: 20
...
```

### _clean_records.py_

Удаляет из *.fasta-файла последовательности, содержащие "sewage" в заголовке - для избавления от образцов из сточных вод (использован на примере Mamastrovirus_1_complete_genome_records.fasta из папки genbank/)
					
### _fasta_header_editor.py_

Добавляет в заголовки *.fasta-файлов информацию о серотипах (например, _serotype_1), используя файл с аннотацией iTOL, где каждому серотипу присвоен свой цвет

При запуске из папки скрипт сработает на файлах из старого репозитория, на которых работал изначально:
* ORF2_cleared_alignment.fas из папки old/aln_old/
* ORF2_RNA_tree_ann.txt из папки old/trees_old/old_annot/
					
### _metadata_editor.py_

Аналогичен fasta_header_editor.py, но корректирует metadata.csv, превращая его в metadata_serotypes_fixed.csv
					
### _fix_rdp_csv.py_

*.fasta-файлы, полученные на выходе RDP4, содержали неполые заголовки, этот скрипт их восстанавливает до полной длины по выравниванию ORF_1A_1B_2_conc_RNA_no_amb.fas из папки aln/main/
					
### _gap_check.py_

Запускается через командную строку. Использовался на выравниваниях ORF_1A_1B_2_conc_RNA_less_amb_fewer_gaps.fas и ORF_1B_2_conc_RNA_frame_fixed_less_amb_fewer_gaps.fas из папки aln/main/recomb_analysis_aln/.
Папка со скриптом также содержит выходные файлы.

Удаляет последовательность из выравнивания, если в ней длина цепочки гэпов превышает пороговое значение
Помимо обязательного пути к исходному файлу (-i) программе можно подать на вход аргументы -oc (путь к выходному csv-файлу с длинами всех гэпов всех  последовательностей в выравнивании), -of (путь к выходному fasta-файлу, очищенному от замусоренных гэпами последовательностей), -t (пороговое значение для числа подряд идущих гэпов, 50 по умолчанию), -e (позволяет исключить из рассмотрения фиксированную длину гэпов, когда мы знаем, что именно этот участок определенной длины учитывать не хотим)

### _sample_exp.py_

С помощью алгоритма BLAST выравнивает пользовательские последовательности на последовательности из GenBank. Запускается из командной строки.

Аргументы: 
*-i - путь к подаваемому на вход файлу .fasta;*
*-o - путь к папке с выходными файлами (будет создана), необязательный аргумент: по умолчанию присваивается название с текущими датой и временем;*
*-l - левая граница выбранного для анализа участка (опционально), по умолчанию - начало последовательности*
*-r - правая граница участка (опционально), по умолчанию - конец последовательности*
*-si - сохранение промежуточных файлов (query.fa, blast_results.xml), по умолчанию - False*

На вход подается выравнивание либо обычная последовательность/несколько последовательностей. На основе этого файла создается *query.fa*, где содержатся последовательности, готовые к прогону через веб-сервер BLAST. Результаты прогона сохраняются в файл *blast_results.xml* - он далее обрабатывается до состояния итогового *blast_query_match.txt*, в котором дается соответствие *query-match* (разделитель - "//", т.к. запятых в поле *match* и так много).

------------------------------------------------------------------------

## **aln/**

Contains alignments of various origins:

* Subfolder main/: main alignments currently in use
  * Subfolder tree_aln/: alignments used for tree building in IQ-TREE, ending with 'no_amb'. IQ-TREE files (see folder trees/iqtree_output/) were later processed through iTOL for serotype coloring (see folder trees/iTOL_output/)
    - Subfolder recomb_analysis_aln/: alignments used for recombination event analysis, processed with the gap_check.py script (folder scripts/). Includes the folder serotypes_aln/: contains serotype-specific alignments obtained in RDP4 (from alignments ORF_1A_1B_2_conc_RNA_no_amb.fas, ORF_1B_2_conc_RNA_no_amb.fas in folder aln/main/recomb_analysis_aln/)
    - Subfolder pair_dist_aln/: alignments for calculating sequence distances - results are in folder pairwise_distance/
	
*Alignments in subfolders main/recomb_analysis_aln/ and main/pair_dist_aln/ partially overlap - the same files were used for different tasks

------------------------------------------------------------------------------------------------------------------------------------------

## **genbank/**

Contains GenBank files (sequence records) or associated files (record names without sample collection date)

------------------------------------------------------------------------------------------------------------------------------------------

## **genome_coverage/**

Contains files created by the genbank_coverage_complex.py script for genome coverage analysis

reference.fas - reference sequences of serotypes
reference_aln.fas - reference alignments
Final_cov.txt - data for plotting Final_cov_plot.png in folder plots/original
file.fasta - sequences aligned to the references by the script

* Subfolder plots/: contains plots
  - Subfolder original/: original coverage plots, including individual serotype plots
  - Subfolder inkscape_modified/: plots (Final_plot_cov, as well as breakpoint distribution plot and recombination region count matrix, obtained in RDP from alignment ORF_1A_1B_2_conc_RNA_no_amb.fas), converted to vector format and modified with Inkscape - plot data matched with added genome schemes

- The *.svg file recombination_region_count_matrix.svg from folder genome_coverage/plots/inkscape_modified/ exceeded the 25 MB limit, so a link to it is provided in a text file

------------------------------------------------------------------------------------------------------------------------------------------

## **add_info/**

Contains auxiliary files:

* gbid_serotype_corresp.csv - correspondences of gb_id and serotypes
* orf_data.csv - ORF coordinates from the gb record

------------------------------------------------------------------------------------------------------------------------------------------

## **trees/**

Contains tree building output files:

* Subfolder iqtree_output/: IQ-TREE files, without visualization; alignments from folder aln/main/tree_aln/ were used for tree building
* Subfolder iTOL_output/: generated iTOL files - trees in various formats (newick, phyloxml), tree annotations
    - Subfolders pics/ for tree folders: tree images in *.png and *.pdf
	
------------------------------------------------------------------------------------------------------------------------------------------

## **parwise_distance/**

Contains MEGA v11 data on distances between protein and nucleotide sequences for the last two and all three ORFs (subfolders aa/, nucl/), as well as corresponding plots (subfolders aa_plot/, nucl_plot/)

Alignments for this task can be found in folder aln/main/pair_dist_aln/.
Plot headers contain information about the alignments and output *.csv files.

------------------------------------------------------------------------------------------------------------------------------------------

## **RDP/**

Contains RDP-analyzed recombination events:

For the entire dataset alignment
(original alignments - ORF_1A_1B_2_conc_RNA_less_amb_fewer_gaps.fas, ORF_1B_2_conc_RNA_frame_fixed_less_amb_fewer_gaps.fas from folder aln/main/recomb_analysis_aln/):

* Subfolder 1a_1b_2: based on whole-genome alignment
* Subfolder 1b_2: based on two last ORFs alignment

For serotype-specific alignments
(original alignments - aln/main/recomb_analysis_aln/serotypes_aln/): 

* Subfolder 1a_1b_2_serotypes: based on whole-genome alignment
* Subfolder 1b_2_serotypes: based on two last ORFs alignment

------------------------------------------------------------------------------------------------------------------------------------------

## **old/**

Contains obsolete project files, replaced by updated ones (may be more useful than files not transferred here from the old repository)

* Subfolder aln_old/: intermediate, no longer relevant alignments; moved here from subfolder aln/main/
* Subfolder rdp_old/: old RDP projects and output information in *.csv format
* Subfolder trees_old/: old tree-related files. Folder old_annot/: with outdated annotations

## **scripts/**

Contains project scripts along with input and output files, even if they duplicate files in other folders (so the scripts folder can be used independently of its location)

### _genecds_create_csv.py_

Program for creating a csv file from GenBank record CDS, genes fields

Sequentially runs two functions:

1) launch_gene_cds_field_distribution()
The function runs the script gene_cds_field_distribution.py with the following internal command:
python genecds_field_distribution.py -input Mamastrovirus_1_complete_genome_records.gb -odir alignment_of_orfs_output -oname Mamastrovirus_1_gene_product_output.txt
The output file is the same as the one generated by the script: Mamastrovirus_1_gene_product_output.txt

2) gene_cds_field_distribution_output_handle()
The function creates the required csv file based on the output file from the previous function
If the CDS/gene field string is orf1/orf1a, or if it contains "unknown", "protease", "NS", then "1A" is added after the comma.
String – "orf1b" or contains "polymerase" – "1B" is added.
In the string "Non" – "1AB" is added.
In the string "non" – "1A", "1B", and "1AB" are added.
Otherwise, "2" is added.
The script genecds_field_distribution.py and the GenBank record file Mamastrovirus_1_complete_genome_records.gb must be in the same directory as genecds_create_csv.py.

					
### _parse_gb_records.py_

Program for extracting specific data from GenBank records from NCBI

Uses four libraries: pandas, biopython, collections, pickle, and the os.path module.
- Initially, a dictionary res is created, which will contain part of the GenBank record data.
- Then the presence of the file dict_data.txt, containing the res dictionary data, is checked. The file is created after the first run of the program and is necessary to avoid re-downloading GenBank records on each new run of the program – convenient if, for example, you need to change the data output format. To re-download data, delete dict_data.txt before running the script.
- Next, if dict_data.txt is absent, the get_gb_records(query) function is executed, where query is the request sent to the NCBI server (to the NCBI Nucleotide database). The function performs three tasks: downloading data from GenBank, recording some qualifiers' data into the res dictionary, and saving all downloaded GenBank records to the records.gb file.
If dict_data.txt already exists, the GenBank record data previously saved in the res dictionary is loaded from it, and the get_gb_records function is ignored.
- Then, the file_writer() function saves the res dictionary data to a pandas dataframe, from where it is saved to the metadata.csv file (now replaced by metadata_serotypes_fixed.csv)
- Part of the data saved in the res dictionary is output to the main_output.txt file.

The query variable was assigned the following value:
```
"Mamastrovirus 1"[porgn] AND "6000"[SLEN] : "8000"[SLEN]
```
The content of main_output.txt is partially shown below (result of data download by query, translated from Russian):
```
125 records loaded

Samples were collected from the following locations (most frequently China: 19 times): 
China: 19
France: 15
Brazil: 15
...

Virus host species: 
Homo sapiens: 84
NaN: 25
...

Sample collection dates: 
NaN: 17
2015: 10
2014: 5
...

Strains: 
NaN: 75
Oxford: 8
Dresden: 2
...

Known information about serotypes/patients: 
NaN: 35
Child with acute gastroenteritis: 3
V1347: 2
...

Sample sources: 
NaN: 47
stool: 20
sewage: 20
...
```

### _clean_records.py_

Removes sequences with "sewage" in the header from a *.fasta file - used to discard sewage samples (used with Mamastrovirus_1_complete_genome_records.fasta from the genbank/ folder)
					
### _fasta_header_editor.py_

Adds serotype information (e.g., _serotype_1) to *.fasta file headers using an iTOL annotation file where each serotype is assigned a color

When run from the folder, the script works on files from the old repository it originally processed:
* ORF2_cleared_alignment.fas from folder old/aln_old/
* ORF2_RNA_tree_ann.txt from folder old/trees_old/old_annot/
					
### _metadata_editor.py_

Similar to fasta_header_editor.py but corrects metadata.csv, turning it into metadata_serotypes_fixed.csv
					
### _fix_rdp_csv.py_

*.fasta files output by RDP4 contained incomplete headers

, this script restores them to full length based on alignment ORF_1A_1B_2_conc_RNA_no_amb.fas from folder aln/main/
					
### _gap_check.py_

Run via the command line. Used on alignments ORF_1A_1B_2_conc_RNA_less_amb_fewer_gaps.fas and ORF_1B_2_conc_RNA_frame_fixed_less_amb_fewer_gaps.fas from folder aln/main/recomb_analysis_aln/.
The script folder also contains output files.

Removes a sequence from alignment if the length of its gap chain exceeds the threshold value.
Besides the mandatory input file path (-i), the program can accept arguments -oc (path to the output csv file with gap lengths of all sequences in the alignment), -of (path to the output fasta file cleaned of gap-heavy sequences), -t (threshold value for the number of consecutive gaps, 50 by default), -e (excludes a fixed gap length from consideration, when we know that this specific length should not be counted).

### _sample_exp.py_

Uses the BLAST algorithm to align user-provided sequences with sequences from GenBank. Designed to be executed from the command line.

Arguments:
*-i - path to the input .fasta file;*
*-o - path to the output directory (will be created), optional argument: defaults to a name with the current date and time;*
*-l - left boundary of the region selected for analysis (optional), defaults to the start of the sequence;*
*-r - right boundary of the region (optional), defaults to the end of the sequence;*
*-si - save intermediate files (query.fa, blast_results.xml), defaults to False.*

The input can be either an alignment or a single/multiple sequence(s). Based on this file, a *query.fa* file is generated, containing the sequences ready to be run through the BLAST web server. The results of the run are saved in *blast_results.xml*, which is further processed into the final *blast_query_match.txt* file. This file contains the *query-match* correspondence, with "//" as the delimiter (since the *match* field already contains many commas).
