# -*- coding: utf-8 -*-
#
"""
ARGEOS programm
writen by Gleb Gavrish
email: ggavrish@fbb.msu.ru

INPUT: file with a set of requests
Intermediate files: input_ArEx.txt; input_GSE.txt
OUTPUT: output_argeos.tsv; errors_argeos.txt; table_term.tsv; table_term_ArEx.tsv
Optional output files: output_argeos.txt

The program takes a file with ID GSE as input, and outputs all the necessary information to output files
The output is three files: a human-readable text file, a table -
suitable for sorting results, as well as a file with errors, where bugged IDs are written
"""
import click
import csv
import tarfile
import shutil
import requests
from bs4 import BeautifulSoup
from tqdm import tqdm
import entrezpy.esummary.esummarizer
import entrezpy.esearch.esearcher
import entrezpy.efetch.efetcher
import os
from contextlib import contextmanager, redirect_stderr, redirect_stdout
import logging
from copy import deepcopy
from pymed import PubMed
import sys
import subprocess

# ГЛОБАЛЬНЫЕ ПЕРЕМЕННЫЕ
if True:
    # NCBI требует сообщать им почту разраба и название проги
    email = "ggavrish@fbb.msu.ru"
    tool = "ARGEOS"
    # Тут определяются адреса директорий, отностительно расположения файла программы
    dirname = os.path.dirname(os.path.abspath(__file__))
    # А это директория для временных файлов
    tmp_dir = os.path.join(dirname, "argeos_tmp")
    Error_List = []  # задаем пустой лист, для повторного анализа плохих ID (Глобальная переменная!)
    reader = csv.reader(open(os.path.join(dirname, 'dict_if_final.csv'), 'r'))
    dict_if = {}
    for k, v in filter(None, reader):
        dict_if[k] = v
    VerboseG = False  # G от слова global
    Cell_size_for_tsv = 50000
    with open(os.path.join(dirname, 'bad_ids.txt'), 'r') as badids:
        lines = badids.readlines()
    Bad_IDs = [x.strip() for x in lines]


class SeriesInfo:
    # класс для удобной записи переменных датасета
    GSE = None
    Authors = None
    organism = None
    samples = None
    Type = None
    platform = None
    title = None
    sub_date = None
    Summary = None
    Overall_design = None
    BioProject = None
    BioProj_EBI = None
    SRA = None


class PubMedInfo:
    # класс для работы с данными из PubMed
    journal = None  # полное название журнала
    doi = None  # doi айдишник
    impfact = None  # импакт фактор журнала
    pbid = None  # id для поиска статьи в pubmed и вытаскивания
    title = None  # Название статьи
    full_refs = None  # Все статьи данной записи


class GsmInfo:
    Cell_type = None
    Treatment = None
    Growth = None
    Type_mol = None
    Extr_prot = None
    Characteristics = None
    All_protocols = None


class GseInfo:
    # И будет создан один класс, что бы править ими всеми
    series_info = SeriesInfo()
    PubMed_info = PubMedInfo()
    gsm_info = GsmInfo()


# Технические функции
def chunks(lst, chunk_size):
    """
    Вход: лист и переменая int - по сколько переменныых делим лист
    Выход: лист листов, в котороых фиксированное (или меньшее) кол-во переменных
    Я использую это для скаивания данных 'пакетами', а не индивидуально - часть оптимизации связи с сервером
    (больше данных за один сеанс связи)
    """
    return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]


def mydeepcopy(gse_mega):
    """
    Какие-то баги возникали из-за обычного копирования, пытался решить дипкопией.
    Не помню помогло по итогу или нет, но лучше не трогать :)
    """
    new_gse_mega = GseInfo()
    new_gse_mega.series_info.Authors = deepcopy(gse_mega.series_info.Authors)
    new_gse_mega.series_info.organism = deepcopy(gse_mega.series_info.organism)
    new_gse_mega.series_info.samples = deepcopy(gse_mega.series_info.samples)
    new_gse_mega.series_info.Type = deepcopy(gse_mega.series_info.Type)
    new_gse_mega.series_info.platform = deepcopy(gse_mega.series_info.platform)
    new_gse_mega.series_info.title = deepcopy(gse_mega.series_info.title)
    new_gse_mega.series_info.sub_date = deepcopy(gse_mega.series_info.sub_date)
    new_gse_mega.series_info.Summary = deepcopy(gse_mega.series_info.Summary)
    new_gse_mega.series_info.Overall_design = deepcopy(gse_mega.series_info.Overall_design)
    new_gse_mega.series_info.BioProject = deepcopy(gse_mega.series_info.BioProject)
    new_gse_mega.series_info.BioProj_EBI = deepcopy(gse_mega.series_info.BioProj_EBI)
    new_gse_mega.series_info.SRA = deepcopy(gse_mega.series_info.SRA)
    new_gse_mega.PubMed_info.journal = deepcopy(gse_mega.PubMed_info.journal)
    new_gse_mega.PubMed_info.doi = deepcopy(gse_mega.PubMed_info.doi)
    new_gse_mega.PubMed_info.impfact = deepcopy(gse_mega.PubMed_info.impfact)
    new_gse_mega.PubMed_info.pbid = deepcopy(gse_mega.PubMed_info.pbid)
    new_gse_mega.PubMed_info.title = deepcopy(gse_mega.PubMed_info.title)
    new_gse_mega.PubMed_info.full_refs = deepcopy(gse_mega.PubMed_info.full_refs)
    new_gse_mega.gsm_info.Cell_type = deepcopy(gse_mega.gsm_info.Cell_type)
    new_gse_mega.gsm_info.Treatment = deepcopy(gse_mega.gsm_info.Treatment)
    new_gse_mega.gsm_info.Growth = deepcopy(gse_mega.gsm_info.Growth)
    new_gse_mega.gsm_info.Type_mol = deepcopy(gse_mega.gsm_info.Type_mol)
    new_gse_mega.gsm_info.Extr_prot = deepcopy(gse_mega.gsm_info.Extr_prot)
    new_gse_mega.gsm_info.Characteristics = deepcopy(gse_mega.gsm_info.Characteristics)
    new_gse_mega.gsm_info.All_protocols = deepcopy(gse_mega.gsm_info.All_protocols)
    return new_gse_mega


def cell_splitter(enter_string):
    """
    Прото разбивает строку на несколько разделенных табом
    """
    # немного логики
    if enter_string is not None:
        n = Cell_size_for_tsv
        if n == 0:
            return enter_string
        elif n < 100:
            return enter_string
        # когда все проверили, просто разбиваем
        exit_str = ""
        if len(enter_string) > n:
            numb = len(enter_string) // n
            for i in range(numb + 1):
                exit_str = exit_str + enter_string[:n - 1] + "\t"
                enter_string = enter_string[n - 1:]
        else:
            exit_str = enter_string
        return exit_str.strip("\t")
    else:
        print("error! string is None type", file=sys.stderr)
        return " "


@contextmanager
def suppressor(a):
    a = True
    if a:
        with open(os.devnull, 'w') as fnull:
            with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
                yield (err, out)
    else:
        yield


# Основные функции
def sist_search(filename, maxterms, output_dir):
    """
    Вход: файл input_terms.txt - построчная запись запросов
    Выход: запись в два файла: input_GSE.txt - уникальные GSE ID для анализа,
    table_term.tsv - таблица запрос - кол-во найденных датасетов
    Внимание! К запросам добавлять GSE не нужно, это делается тут
    Функция считывает запросы из файла, добавляет пункт для поиска только датасетов, и возвращает список датасетов.
    Попутно она записывает в отдельный файл таблицу с количесвом найденных датасетов по каждому запросу
    """

    def id_to_gse(id):
        """
        вход: id
        выход: GSE id
        функция нужная для преобразлвания выдачи в формат GSE id, ибо ftp ищется по нему
        """
        number = int(id[1:])  # удаляем первый символ ('2') и удаляем нули в начале перобразованием в int
        gse = "GSE" + str(number)  # склеиваем на id
        return gse

    with open(os.path.join(dirname, filename), "r") as input_file:  # открываем файл с запросами
        term_list = input_file.readlines()  # читаем запросы в лист
    table_term = open(os.path.join(output_dir, "table_term.tsv"), "w")  # открываем файл для записи таблицы
    table_term.write("term" + '\t' + "number of found datasets (GEO)" + '\n')  # записываем шапку таблицы
    glob_list = []  # лист куда сливаю все GSE. Причем слив происходит при работе самой утилиты, т.е. в конце она выдает
    # общую сумму всех GSE по всем запросам. Для этого ввожу переменную i
    i = 0  # переменная для подсчета кол-ва датасетов найденный в конкретном запросе (см коммент выше)

    for term in tqdm(term_list):  # проходимся по всем запросам
        # if True:
        if "|" in term:  # если есть организм
            term = "(" + term.strip() + " AND gse[ETYP]"  # добавляем фильтрацию по датасетам (только GSE)
            term = term.replace(" |", "|")  # если перед палкой пробел, то это его удалит
            term = term.replace("|",
                                ") AND")  # Тут заменяем | у ORGN на AND (в других бд orgn в других разделах пишут)
        else:  # если нет организма
            term = "(" + term.strip() + ") AND gse[ETYP]"
        with suppressor(True):
            logging.disable(logging.CRITICAL)  # Выключаю @#%^*$ логи entrezpy (ЭТУ СТРОЧКУ Я ИСКАЛ 4 МЕСЯЦА!!!)
            e = entrezpy.esearch.esearcher.Esearcher(tool, email)
            a = e.inquire({"db": "gds", "term": term, "retmax": maxterms, "rettype": "uilist"})
            # НЕ ЗАБЫТЬ! " (кавычка) = %22
            id_list = a.get_result().uids  # считываю результат
            logging.disable(logging.NOTSET)  # возвращаем логи
        i_new = len(id_list)  # узнаем длинну текущего массива айдишников
        table_term.write(term + '\t' + str(i_new - i) + '\n')  # записываем инфу по запросу, вычитаем предыдущие
        i = i_new
        glob_list = id_list

    gse_list = []  # сервер выдает  результаты в виде просто id, но для ftp мне нужен формат GSE id, для этого лист
    for unique_id in set(glob_list):  # чтобы оптимизировать работу, проходимся только по уникальным GSE ID
        gse_id = id_to_gse(unique_id)  # переводим id в gse спец функцией
        gse_list.append(gse_id)
    num = 0
    with open(os.path.join(output_dir, "input_GSE.txt"), "w") as input_gse:  # открываем файл на запись GSE ID
        for gse in gse_list:
            input_gse.write(gse + '\n')  # записываем GSE ID в файл
            num = num + 1  # подсчитываем количестов уникальных gse_id

    table_term.write("TOTAL GEO (UNIQUE)\t" + str(num) + '\n')  # добавляем строку со суммарными находками(уникальными!)
    table_term.close()  # не забываем закрыть файл с результатами!


def gsm_analizator(gsm_list):
    """
    Вход: Лист из GSM ID
    Выход: Вся нужная инфа (тип клеток, treatment protocol) записанная в специальный класс (чтоб на выход одна
    переменная)
    Функция проходится по блокам сэмплов и вытаскивает из каждой нужную инфу.
    Дополнительно фильтрует ее, на выходе получаются только уникальные знаения.
    По сути работает как основной код, но просто смотрит много однотипных страниц, и выдает только уникальные значения.
    """
    tr_list = []  # лист для Treatment protocol
    gr_list = []  # лист для Growth protocol
    cell_types = []  # лист для типа клеток
    tp_list = []  # лист для типа экстрагируемой молекулы (?)
    exp_list = []  # лист для Extraction protocol (протокола выделения пробы)
    char_list = []  # лист для Characteristics (вся прочая инфа из образцов)
    # Упоминая переменные
    growth = "None"
    # cells = "None"
    # mol_type = "None"
    treatment = "None"
    for GSM in gsm_list:
        for chanel in GSM.find_all("channel"):
            try:
                treatment = chanel.find_all("treatment-protocol")[0].get_text().strip()
            except Exception:
                pass
            # Тут вначале обращаюсь к блоку Chanle (видимо бывают разные каналы для одного образца, нужно это учесть)
            try:
                growth = chanel.find_all("growth-protocol")[
                    0].get_text().strip()  # если нет нормального Treatment протокола,
                # то обычно инфа записанна здесь. А если уж тут нет, то скорее всего вообще другой тип эксперимента.
            except Exception:
                pass
            try:
                cells = chanel.find_all("characteristics", attrs={"tag": "cell type"})[0].get_text().strip()
            except Exception:
                cells = -1
                pass
            # а потом через characteristics с тэгом "cell type" нахожу нужное значение.
            if cells == -1:
                cells = "Source: " + chanel.find_all("source")[0].get_text().strip()
                # иногда тупо нет графы Cell type, но эксперимент на клетках.
                # Тогда обычно они записывают инфу тут, но помечаю что это из Source
            mol_type = chanel.find_all("molecule")[
                0].get_text().strip()  # малоинформативно, но может пригодиться
            extprot = chanel.find_all("extract-protocol")[0].get_text().strip()
            try:
                charact_tmplist = chanel.find_all("characteristics")  # Тут выцепляем Charact, но с ним сложнее -- их
                # много и у них часть инфы внутри тэга
                charact = ""
                for char in charact_tmplist:
                    charact = charact + char.get(
                        "tag") + ': ' + char.get_text().strip() + '; '  # объединяем tag и сам текст
            except Exception:
                charact = ""
                pass
            # далее в листы добавляются только уникальные значения
            if treatment not in tr_list:
                tr_list.append(treatment)
            if growth not in gr_list:
                gr_list.append(growth)
            if cells not in cell_types:
                cell_types.append(cells)
            if mol_type not in tp_list:
                tp_list.append(mol_type)
            if extprot not in exp_list:
                exp_list.append(extprot)
            if charact not in char_list:
                char_list.append(charact)
    out_info = GsmInfo()
    out_info.Treatment = str(tr_list).strip('[]')
    out_info.Cell_type = str(cell_types).strip('[]')
    out_info.Growth = str(gr_list).strip('[]')
    out_info.Type_mol = str(tp_list).strip('[]')
    out_info.Extr_prot = str(exp_list).strip('[]')
    out_info.Characteristics = str(char_list).strip('[]')
    # Произвожу проверку на пустые параметры, тогда прописываю None
    if len(out_info.Type_mol) == 0:
        out_info.Type_mol = "None"
    if len(out_info.Treatment) == 0:
        out_info.Treatment = "None"
    if len(out_info.Cell_type) == 0:
        out_info.Cell_type = "None"
    if len(out_info.Growth) == 0:
        out_info.Growth = "None"
    if len(out_info.Extr_prot) == 0:
        out_info.Extr_prot = "None"
    if len(out_info.Characteristics) == 0:
        out_info.Characteristics = "None"
    return out_info


def geo_xml_parser(name):
    """
    Программа производит непосредственный анализ XML файла полученного с GEO.
    ТРУБУЮТСЯ ФУНКЦИИ: pub_med_by_id, GSM_analizator
    Вход: имя файла .xml
    Выход: три переменные, каждая своего класса, записанные в переменую мегакласса:
    Series -  инфа по датасэту
    pubmed - инфа из PubMed
    GSM_info - инфа по сэмплам
    """
    with open(name, 'r', encoding='utf-8') as XML_file:
        xml = XML_file.read()
        soup = BeautifulSoup(xml, features="html.parser")
        # Варим суп из файла.
        # МНОГОБУКАФ: Изначально пытался юзать пакет lxml, причем его можно и как паресер для супа использовать.
        # Но с ним были проблемы, тупо не запускался. Так что юзаю стандартный парсер. ОДНАКО, в супе
        # я немного разочеровался: дело в том, что все парсеры раскладывают файл на дерево. И скажем нам нужен
        # конкретный лист из этого дерева,и у него есть уникальный тэг. Хотелось бы тупо по этому тэгу его и выцепить.
        # Но нет, в супе нужно пропсывать полный путь по тэгам :( А уверенности что путь один для всех файлов нет :'(
        # Кароч надеюсь все норм будет. Но если что надо менять пакет для парсенья.
        block1 = soup.find_all('series')[0]
        # Блок1 : тут инфа по датасету, то что можно в саммари найти
        block2 = soup.find_all('platform')
        # тут инфа по платформе, например ее ID
        samples_block = soup.find_all('sample')
        # а тут получаем лист из блоков по каждому образцу. Его потом анализируем функцией GSM_analizator

        # ! Блок серии (датасета)
        if True:
            series = SeriesInfo()
            # с платформой сложности, так что проходимся циклом
            i = 0
            for blochechek in block2:
                if i == 0:
                    platf = blochechek.find_all('accession')[0].get_text().strip()
                else:
                    platf = platf + "; " + blochechek.find_all('accession')[0].get_text().strip()
                i = i + 1
            org_list = []
            for blochechek in block2:
                organ = blochechek.find_all('organism')[0].get_text().strip()
                org_list.append(organ)
            first = True
            for new_org in set(org_list):
                if first:
                    organ = new_org
                    first = False
                else:
                    organ = organ + "; " + new_org
            # platform
            series.platform = platf
            # organism
            series.organism = organ
            # Далее выцепляем инфу по каждому интересуещему параметру в отдельную переменную
            # GSE ID
            series.GSE = block1.find_all('accession')[0].get_text().strip()
            # samples
            series.samples = len(samples_block)
            # Type
            types = block1.find_all('type')
            str_type = ""
            for one_type in types:
                str_type = str_type + one_type.get_text().strip() + "; "
            series.Type = str_type.strip("; ")  # стрип нужен чтоб красиво выводилось
            # title
            series.title = block1.find_all('title')[0].get_text().strip()
            # sub_date
            series.sub_date = block1.find_all('submission-date')[0].get_text().strip()
            # Summary
            series.Summary = block1.find_all('summary')[0].get_text().strip()
            series.Summary = ' '.join(series.Summary.split())
            # Overall_design
            series.Overall_design = block1.find_all('overall-design')[0].get_text().strip()
            series.Overall_design = ' '.join(series.Overall_design.split())  # удалю перенос строки внутри текста

            for links in block1.find_all('relation'):
                if links.get('type') == "BioProject":
                    # BioProject link
                    series.BioProject = links.get('target')
                elif links.get('type') == "SRA":
                    # SRA link
                    series.SRA = links.get('target')
            if series.BioProject != "None" or series.BioProject is not None:  # убедиться, что при отсутсвии get() возвращает None
                bioProj = series.BioProject.split("bioproject/")[1]
                series.BioProj_EBI = "https://www.ebi.ac.uk/ena/browser/view/" + bioProj
                series.BioProject = "https://www.ncbi.nlm.nih.gov/bioproject/" + bioProj
        # ! Блок PUBMED
        if True:
            pubmed = PubMedInfo()
            pub_med_id = block1.find_all('pubmed-id')  # получаю лист со ВСЕМИ сатьями
            # doi
            if type(pub_med_id) != "NoneType":
                lt = []
                for poob in pub_med_id:
                    lt.append(poob.get_text())
                pubmed.pbid = lt
            else:
                pubmed.pbid = "None"
            # pubmed.journal = "None"

        # ! Блок сэмплов
        # Анализируем все сэмплы и получаем с них инфу в виде листов. Всю выдачу функции записываем в GSM_info
        gsm_info = gsm_analizator(samples_block)

        # Удаляем GSM_info после записи его компонентов в переменные (мне удобнее работать с ними)
        all_info = GseInfo()
        all_info.gsm_info = gsm_info
        all_info.PubMed_info = pubmed
        all_info.series_info = series

        all_info.gsm_info.All_protocols = "[Overal design]" + str(
            all_info.series_info.Overall_design) + "; [Treatment]" + str(
            all_info.gsm_info.Treatment) + "; " + "[Growth]" + str(all_info.gsm_info.Growth) + "; [Extraction]" + str(
            all_info.gsm_info.Extr_prot) + "; [Cell type]" + str(all_info.gsm_info.Cell_type) + str(
            all_info.gsm_info.Characteristics)
    return all_info


def pubmed_parser(gse_info_list, pubmed_id_list, ERRORS):
    """"
    Вход: лист из переменных формата GSE_id, но в которых пока нет инфы из pubmed (а только pbid) + файл с ошибками
    Выход: Инфа с pubmed + импакт фактор
    # import metapub - пытался юзать этот пакет, но он битый и не ставится
    """
    pubmed_id_list = [int(i) for i in pubmed_id_list if i != "None"]  # удаляем None значения
    if VerboseG:
        print("Getting info from pubmed: " + str(pubmed_id_list))
    exit_list = []
    if len(pubmed_id_list) > 0:  # это защита от пустого листа, если все ID выдали ошибку
        with suppressor(True):
            logging.disable(logging.CRITICAL)  # отключаем логи. ЭТУ СТРОЧКУ Я ИСКАЛ 4 МЕСЯЦА!!
            e = entrezpy.esummary.esummarizer.Esummarizer(tool, email)
            analyzer = e.inquire({'db': 'pubmed', 'id': pubmed_id_list})
            res = analyzer.get_result().summaries
            logging.disable(logging.NOTSET)  # возвращаем логи
        jurs = list(dict_if.keys())
        new_jurs = []

        for pubmed in gse_info_list:
            pubmed_ret = pubmed
            if pubmed.PubMed_info.pbid is None:  # если нет статьи
                pubmed_ret.PubMed_info.impfact = "None"
                pubmed_ret.PubMed_info.journal = "None"
                pubmed_ret.PubMed_info.full_refs = "None"

            elif len(pubmed.PubMed_info.pbid) == 1:  # если отдна статья
                u_id = int(pubmed.PubMed_info.pbid[0])
                pubmed_ret.PubMed_info.pbid = u_id  # тут тупо вывожу инфу по одной статье
                ful_jur = res.get(u_id).get('fulljournalname')
                imft = dict_if.get(ful_jur.lower())
                if imft is None:
                    imft = dict_if.get(res.get(u_id).get('source').lower())  # пробуем краткое название
                if imft is None:
                    imft = dict_if.get(ful_jur.lower().replace('and', '&'))
                if imft is None:
                    imft = dict_if.get(ful_jur.lower().replace('&', 'and'))
                # все перепробовали, теперь с высокой вероятностью есть impactfactor
                if imft == "None" or imft == "Not Available":
                    imft = "None"
                elif imft is None:
                    imft = "None"
                    new_jurs.append(ful_jur)  # Если все же его нет, то записываем в логи
                pubmed_ret.PubMed_info.journal = ful_jur
                pubmed_ret.PubMed_info.impfact = imft

                doi_ex = False
                for d in res.get(u_id).get('articleids'):
                    if d.get('idtype') == 'doi':
                        u_doi = d.get('value')
                        doi_ex = True
                if doi_ex:
                    # pubmed_ret.PubMed_info.doi = "doi: " + u_doi
                    pubmed_ret.PubMed_info.doi = "https://doi.org/" + u_doi
                else:
                    pubmed_ret.PubMed_info.doi = "pubmed_id = " + str(u_id)
                    ul_doi = "pubmed_id = " + str(u_id)
                pubmed_ret.PubMed_info.title = res.get(u_id).get('title')
                pubmed_ret.PubMed_info.full_refs = "{" + ful_jur + " (" + imft + "): " + u_doi + "}"
            elif len(pubmed.PubMed_info.pbid) > 1:  # если несколько статей
                # Прохожусь по всем pubmed id, выявляю тот, у которого наибольший импакт фактор
                best = ["", -1.0]
                all_refs = ""
                for uniq_pb in pubmed.PubMed_info.pbid:
                    u_id = int(uniq_pb)
                    ful_jur = res.get(u_id).get('fulljournalname')
                    imft = dict_if.get(ful_jur.lower())
                    mb_best = [u_id, imft]
                    if imft == "None" or imft == "Not Available" or imft is None:
                        mb_best[1] = -0.5
                    else:
                        mb_best[1] = float(mb_best[1])
                    if best[1] < mb_best[1]:
                        best = mb_best

                    # тут формируем дополение в полный списко всех ссылок
                    doi_ex = False
                    for d in res.get(u_id).get('articleids'):
                        if d.get('idtype') == 'doi':
                            u_doi = d.get('value')
                            doi_ex = True
                    if imft is None:
                        imft = "None"
                        new_jurs.append(ful_jur)  # to ERRORs
                    if doi_ex:
                        out_str = "{" + ful_jur + " (" + imft + ") doi:" + u_doi + "}; "
                    else:
                        out_str = "{" + ful_jur + " (" + imft + ") pubmed_id = " + uniq_pb + "}; "
                    all_refs = all_refs + out_str
                if best[1] == -0.5:
                    best[1] = "None"
                pubmed_id = best[0]
                doi_ex = False
                for d in res.get(pubmed_id).get('articleids'):
                    if d.get('idtype') == 'doi':
                        ul_doi = d.get('value')
                        doi_ex = True
                if doi_ex:
                    # pubmed_ret.PubMed_info.doi = "doi: " + ul_doi
                    pubmed_ret.PubMed_info.doi = "https://doi.org/" + ul_doi
                else:
                    pubmed_ret.PubMed_info.doi = "pubmed_id = " + str(pubmed_id)
                pubmed_ret.PubMed_info.title = res.get(pubmed_id).get('title')
                # jur = res.get(pubmed_id).get('source')  # Тут краткое имя журнала
                # pubmed_ret.PubMed_info.journal = jur
                jur = res.get(pubmed_id).get('fulljournalname')
                if (jur not in jurs) and (jur not in new_jurs):
                    new_jurs.append(jur)
                pubmed_ret.PubMed_info.impfact = dict_if.get(jur.lower())
                pubmed_ret.PubMed_info.journal = jur
                pubmed_ret.PubMed_info.full_refs = all_refs
            exit_list.append(pubmed_ret)

        # return block
        if len(new_jurs) > 0:
            ERRORS.write("new jurs! " + str(new_jurs) + "\n")
    return exit_list


def xml_by_id(gse_list, ERRORS):
    """
    Вход: лист из GSE ID
    Выход: скаччанные .xml файлы в tmp директории
    Функция берет набор ID и через FTP скачивает архивы в tmp папку, после чего разархивирует их до .xml
    """

    def GSE_to_nnn(gse_id):
        """
        Вход: GSE_id
        ВЫХОД: название наддиректории в виде GSE...nnn (-3 последние символа)
        """
        if len(gse_id) > 5:
            out_n = gse_id[:-3] + "nnn/"
        else:
            out_n = "GSEnnn/"
        return out_n

        # компонеты для ftp подключения к GEO
        host = "ftp.ncbi.nlm.nih.gov"
        ftp = ftplib.FTP(host)
        ftp.login(user="anonymous")
        for GSE_id in tqdm(gse_list):
            if VerboseG:
                tqdm.write("Downloading " + str(GSE_id))  # Важно что использую не print(), т.к. он ломает prog.bar tqdm
            filename = "/geo/series/" + GSE_to_nnn(GSE_id) + GSE_id + "/miniml/" + GSE_id + "_family.xml.tgz"
            out = os.path.join(tmp_dir, str(GSE_id + ".xml.tgz"))
            with open(out, 'wb') as f:
                def callback(data):
                    f.write(data)

                ftp.retrbinary('RETR ' + filename, callback, blocksize=1024)  # последний параметр увелиивает скорость
            tf = tarfile.open(out)
            tf.extractall(path=tmp_dir)
            tf.close()


def arex_search(filename, output_dir):
    """
    Вход: список запросов для GEO
    Выход: лист из id ArrayExpress + таблица с кол-вом записей по каждому запросу
    """
    bhtml = "https://www.ebi.ac.uk/arrayexpress/xml/v3/experiments"  # базовая html строка
    bhtml = bhtml + "?directsub=true&"  # это позволяет отсеить данные импортированные из GEO
    with open(os.path.join(dirname, filename), "r") as input_file:  # открываем файл с запросами
        term_list = input_file.readlines()
    table_term = open(os.path.join(output_dir, "table_term.tsv"), "a")  # открываем файл для записи таблицы
    table_term.write("term" + '\t' + "number of found datasets (ArrayExpress)" + '\n')
    arexp_list = []  # объявляю лист для переменных
    for term in tqdm(term_list):
        the_term = term  # .split("|")
        tmp_list = []
        # Преобразование запроса
        # Тут преобразую запрос из формата для GEO формат для ArrayExpress
        if "|" in term:  # если в запросе есть организм, то обрабатываем его отдельно, как того требует ArrayExpress
            orgn = term.split("|")[1].strip().split("[")[0]
            term = "(" + term.split("|")[0].strip() + ")"
            norm_term = term + "&species=" + orgn
        else:
            norm_term = "(" + term.strip() + ")"
        c = []
        for i in norm_term.split(" "):  # тут для замены пробелов на плюсики (для url запроса!) сплитим по пробелам
            if len(i) > 0:  # если пользователь написал 2 пробела, то этот цикл удалит лишний
                c.append(i)
        norm_term = "+".join(c)  # наконец сливаем порезанную строку аккрутно пробелчиками
        # скачивание запроса
        req = bhtml + "keywords=" + norm_term
        r = requests.get(req)
        soup = BeautifulSoup(r.text, features="html.parser")
        # выцепление всех ArEx ID
        exp = soup.find_all("experiment")
        i = len(exp)  # кол-во находок
        for ex in exp:
            id = ex.find("accession").get_text().strip()
            tmp_list.append(id)
        table_term.write(norm_term.replace("+", " ").replace("&species=", " | ") + '\t' + str(i) + '\n')
        arexp_list = arexp_list + tmp_list
    # Фильтрация уникальных значений
    arexp_list = list(set(arexp_list))
    num = 0
    with open(os.path.join(output_dir, "input_ArEx.txt"), "w") as input_ae:  # открываем файл на запись GSE ID
        for ae in arexp_list:
            input_ae.write(ae + '\n')  # Проходимся по ID, записываем их в файл
            num = num + 1
    table_term.write("TOTAL ArrayExpress (UNIQUE)\t" + str(num) + '\n')  # дописываем кол-во уникальных
    table_term.close()  # не забываем закрыть файл


def array_express(id):
    """
    Вход: ArEx ID для анализа
    Выход: переменная для печати (в мегаформате)
    По итогу решил забить, и для каждой записи в индивидуальном порядке скачивать инфу в xml
    """
    bhtml = "https://www.ebi.ac.uk/arrayexpress/xml/v3/experiments"
    request = bhtml + '/' + id  # формируем строку запроса
    # далее получаю инфу через url запрос
    # if VerboseG:
    #     print(request, file=sys.stderr)
    r = requests.get(request)  # делаю запрос
    # with open(os.path.join(tmp_dir, 'tmp.html'), 'w', encoding='utf-8') as tmp_file:
    #     tmp_file.write(r.text)  # Кэширую страницу во временном файле, потом работаю с файлом.
    # with open(os.path.join(tmp_dir, 'tmp.html'), 'r', encoding='utf-8') as tmp_html:
    #     xml = tmp_html.read()  # читаю файл в переменную, после чего закрываю файл
    soup = BeautifulSoup(r.text, features="html.parser")  # парсим супом, формирую дерево разбора
    ae = mydeepcopy(GseInfo())  # объявляю переменную МЕГАкласса
    # тут я пока стараюсь выявить соответсвия по тегам и категориям в GEO
    ae.series_info.GSE = id
    if type(soup.find("name")) != 'NoneType':
        ae.series_info.title = soup.find("name").get_text().strip()
    else:
        ae.series_info.title = soup.find("experiment").find("name").get_text().strip()
    ae.series_info.organism = soup.find("organism").get_text().strip()
    ae.series_info.samples = soup.find("experiments").get("total-samples")
    ae.series_info.sub_date = soup.find(
        "releasedate").get_text().strip()  # Важно! формат даты аналогичен NCBI (проверил)
    if soup.find("experimenttype") is not None:
        ae.series_info.Type = soup.find("experimenttype").get_text().strip()  # Формат у них не совпадает с GEO -_-
        # А вот сейчас попытаемся перевести в формат GEO ;)
        # dict_to_geo_types = {}
        # mb_to_geo =
    # Анализ референсной статьи
    bib = soup.find_all("bibliography")
    title_list = []
    if bib is not None and len(bib) > 0:
        for bib0 in bib:
            doi_tmp = bib0.find("doi")
            if doi_tmp is not None:
                ae.PubMed_info.doi = doi_tmp.get_text().strip()
            else:
                ae.PubMed_info.doi = "None"
            title_tmp = bib0.find("title")
            if title_tmp is not None:
                ae.PubMed_info.title = title_tmp.get_text().strip()
                title_list.append(title_tmp.get_text().strip())
            else:
                ae.PubMed_info.title = "None"
            pub_tmp = bib0.find("publication")
            if pub_tmp is not None:
                ae.PubMed_info.journal = pub_tmp.get_text().strip()
                ae.PubMed_info.impfact = dict_if.get(ae.PubMed_info.journal.lower())
            else:
                ae.PubMed_info.journal = "None"
                ae.PubMed_info.impfact = "None"
    else:
        # print(soup)
        ae.PubMed_info.doi = "None"
        ae.PubMed_info.title = "None"
        ae.PubMed_info.journal = "None"
        ae.PubMed_info.impfact = "None"
    if len(title_list) > 0:
        ae.PubMed_info.title = title_list

    # ! Блок анализа протоколов
    ae.series_info.Summary = soup.find("description").get_text().strip()  # видимо по смыслу это summary (?)
    # ! Блок прохождения по всем протоколам
    protocols = soup.find_all("protocol")
    total_protocols = ""
    for prot in protocols:
        shmotocol = protocol_analyzer(str(prot.find("id").get_text().strip()))
        if shmotocol is not None:
            total_protocols = total_protocols + shmotocol
    is_gro = False
    try:
        ae.gsm_info.Growth = soup.find("experimentdesign").get_text().strip()
        is_gro = True
    except Exception:
        ae.gsm_info.Growth = "None"
    if is_gro:
        total_protocols = "[Growth] " + ae.gsm_info.Growth + "; " + total_protocols
    ae.gsm_info.All_protocols = total_protocols
    if type(soup.find_all("arraydesign")) != 'NoneType':
        if len(soup.find_all("arraydesign")) == 1:
            ae.series_info.platform = soup.find("arraydesign").find(
                "accession").get_text().strip()  # вроде только для чипов, для секов нет платформ
        elif len(soup.find_all("arraydesign")) >= 1:
            i = 0
            for arplat in soup.find_all("arraydesign"):
                if i == 0:
                    final_ar_plat = arplat.find("accession").get_text().strip()
                else:
                    final_ar_plat = final_ar_plat + "; " + arplat.find("accession").get_text().strip()
                i = i + 1
            ae.series_info.platform = final_ar_plat
        else:
            ae.series_info.platform = "None"
    # конец анализа, возвращаю готовую переменную
    return mydeepcopy(ae)


def pbid_by_title(title_list):
    """
    Вход: лист названий
    Выход: словарь название = pubmed id
    Ищет в pubmed статьи по названиям
    Вообще этим пактеом сильно удобнее вытаскивать всю инфу, но A) времени нет и Б) работает -- не трогай!
    """
    term = "(" + "[Title]) OR( ".join(title_list) + "[Title])"
    pubmed = PubMed(tool=tool, email=email)
    results = pubmed.query(term, max_results=5)
    dict_for_out = {}
    list_of_ids = []
    for article in results:
        if article.title.strip(".") in term:
            dict_for_out[article.title.strip(".")] = article.pubmed_id.split("\n")[0]
            list_of_ids.append(article.pubmed_id.split("\n")[0])
    return [dict_for_out, list_of_ids]


def protocol_analyzer(id):
    """
    Вход: id протокола из Array Express
    Выход: вся информация по протоколу, преобразованная в строку.
    """
    bhtml = "https://www.ebi.ac.uk/arrayexpress/xml/v3/protocols/"
    request = bhtml + id  # формируем строку запроса
    # далее получаю инфу через url запрос
    if VerboseG:
        print("getting protocol by id: " + id, file=sys.stderr)
        print(request)
    r = requests.get(request)  # делаю запрос
    soup = BeautifulSoup(r.text, features="html.parser")  # парсим супом, формирую дерево разбора
    is_fail = False
    text = ""
    typ = ""
    if soup.find("type") is not None:
        typ = soup.find("type").get_text()
    else:
        is_fail = True
    if soup.find("text") is not None:
        text = soup.find("text").get_text()
    else:
        is_fail = True
    if len(text) > 0 and len(typ) > 0 and not is_fail:
        total_information_string = '"[' + id + '] ' + typ + ': ' + text + '"; '
    else:
        is_fail = True
    if is_fail:
        return None
    else:
        return total_information_string


def split_to_unique(glist):
    """
    Вход и Выход: Лист МЕГАФОРМАТА
    Разбивает запись на несколько, если есть несколько организмов или типов экспериментов
    Это нужно для удобной фильтрации конечного листа записей (стандартизация формата столбца)
    """
    exit_glist = []
    prom_list = []
    # Проходимся по организмам
    for gmega in glist:
        if gmega.series_info.GSE is not None:
            if ";" in gmega.series_info.organism:
                org_list = gmega.series_info.organism.split(";")
                for org in org_list:
                    new_mega = deepcopy(gmega)  # Очень важно делать глубокую копию!
                    # Иначе у двух переменных ссылка на один и тот же объект
                    new_mega.series_info.organism = org.strip()
                    prom_list.append(new_mega)
            else:
                prom_list.append(gmega)
    # Проходимся по типам экспериментов
    for gmega2 in prom_list:
        if ";" in gmega2.series_info.Type:
            for typ in gmega2.series_info.Type.split(";"):
                new_mega = deepcopy(gmega2)
                new_mega.series_info.Type = typ.strip(" ")
                exit_glist.append(new_mega)
        else:
            exit_glist.append(gmega2)
    return exit_glist


def text_output(listochek, output_text_file):
    """
    Вход: лист из переменных мегаформата и файл для записи
    Выход: запись в файл текстового формата
    """
    for gse_info in listochek:
        if gse_info.series_info.GSE is not None:
            output_text_file.write('#\n')
            output_text_file.write('GSE\t' + str(gse_info.series_info.GSE) + '\n')  # GSE ID
            output_text_file.write('BPJ\t' + str(gse_info.series_info.BioProject) + '\n')  # BioProject
            output_text_file.write('SRA\t' + str(gse_info.series_info.SRA) + '\n')  # SRA
            output_text_file.write('ORG\t' + str(gse_info.series_info.organism) + '\n')  # Orgamism or Organisms
            output_text_file.write('TYP\t' + str(gse_info.series_info.Type) + '\n')  # тип эксперимента
            output_text_file.write('SPS\t' + str(gse_info.series_info.samples) + '\n')  # number of Samples
            output_text_file.write('PLT\t' + str(gse_info.series_info.platform) + '\n')  # Platform
            output_text_file.write('TTL\t' + str(gse_info.series_info.title) + '\n')  # Title
            output_text_file.write('SBD\t' + str(gse_info.series_info.sub_date) + '\n')  # Submison date
            output_text_file.write('DOI\t' + str(gse_info.PubMed_info.doi) + '\n')  # doi
            output_text_file.write('JUR\t' + str(gse_info.PubMed_info.journal) + '\n')  # journal name
            output_text_file.write('IMP\t' + str(gse_info.PubMed_info.impfact) + '\n')  # impact factor
            output_text_file.write('PBI\t' + str(gse_info.PubMed_info.pbid) + '\n')  # pubmed id
            output_text_file.write('PTL\t' + str(gse_info.PubMed_info.title) + '\n')  # Название статьи
            output_text_file.write('SUM\t' + str(gse_info.series_info.Summary) + '\n')  # Summary
            output_text_file.write('OVD\t' + str(gse_info.series_info.Overall_design) + '\n')  # Overall design
            output_text_file.write('CLT\t' + str(gse_info.gsm_info.Cell_type) + '\n')  #
            output_text_file.write('TRE\t' + str(gse_info.gsm_info.Treatment) + '\n')
            output_text_file.write('TYM\t' + str(gse_info.gsm_info.Type_mol) + '\n')
            output_text_file.write('GRP\t' + str(gse_info.gsm_info.Growth) + '\n')


def table_output(listochek, output_table_file):
    """
    Вход: лист из переменных мегаформата и файл для записи
    Выход: запись в файл таблицы инормации по данному GSE
    """
    # test = gse_info.gsm_info.Type_mol
    for gse_info in listochek:
        if gse_info.series_info.GSE is not None:
            if VerboseG:
                tqdm.write("writing " + gse_info.series_info.GSE, file=sys.stderr)
            all_prot = cell_splitter(deepcopy(gse_info.gsm_info.All_protocols))
            id = gse_info.series_info.GSE
            if "GSE" in id:
                link = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + id  # Если GEO
            else:
                link = "https://www.ebi.ac.uk/arrayexpress/experiments/" + id  # Если Array Express
            # посредственно запись строки
            output_table_file.write(str(gse_info.series_info.GSE) + "\t" + str(gse_info.series_info.organism) + "\t" +
                                    str(gse_info.series_info.samples) + "\t" + str(gse_info.series_info.Type) + "\t" +
                                    str(gse_info.series_info.platform) + "\t" +
                                    str(gse_info.series_info.title) + "\t" +  # str(gse_info.series_info.Authors) + "\t" +
                                    str(gse_info.series_info.sub_date) + "\t" + str(
                gse_info.series_info.Summary) + "\t" + link + "\t" + str(gse_info.PubMed_info.title) + "\t" + str(
                gse_info.PubMed_info.journal) + "\t" +
                                    str(gse_info.PubMed_info.impfact) + "\t" + str(gse_info.PubMed_info.doi) + "\t" +
                                    str(
                                        gse_info.PubMed_info.full_refs) + "\t" +
                                    str(gse_info.gsm_info.Type_mol) + "\t" + str(
                gse_info.series_info.BioProject) + "\t" +
                                    str(gse_info.series_info.BioProj_EBI) + "\t" + str(
                gse_info.series_info.SRA) + "\t" + all_prot + "\n")
        else:
            print("error with table", file=sys.stderr)


@click.command()
@click.option('--input_file', '-i', default="input_terms.txt", show_default=True, help="название файла для входа")
@click.option('--output', '-o', default="argeos_output", show_default=True,
              help="название директории с результатом")
@click.option('--verbose', '-v', is_flag=True, help="Более подробные сообщения во время работы")
@click.option('--text_out', '-t', is_flag=True, help="Добавит к выдаче текствый файл")
@click.option('--unique', '-u', is_flag=True, help="Не разбивать строки на уникальные (по типу и организмам)")
@click.option('--cell_size', '-l', default=50000, show_default=True,
              help="Ограничение ячейки для таблицы, для корректной вставики в Exel/Google sheets. "
                   "0 если нет ограничений")
@click.option('--chunk_size', '-c', default=3, show_default=True,
              help=str("Переменная, определяющая величину пакета для скачивания. Если больше - меньше сеансов связи, " +
                       "больше занимаемого места tmp директорией (и наоборот)"))
@click.option('--mode1', is_flag=True, help="Только поиск, без анализа данных")
@click.option('--mode2', is_flag=True, help="Только анализ, без поиска (входной файл input_GSE.txt)")
def main(input_file, output, text_out, chunk_size, mode1, mode2, verbose, unique, cell_size):
    """
    Программа разработанна для аннатоирования результатов поиска в базах данных GEO и ArrayExpress. На вход программа
    принимает один или нескольуо поисковых запросов, записанных на разных строках. На выходе, в output ректории
    получаются две таблицы:
    количество находок для каждого запроса и базы данных (table_term.tsv);
    таблица с подробной информацией по каждой находке;
    два файла со списком найденных ID

    Рекомендации по созданию запросов: ключевые слова писать через пробел и AND. Если нужно добавить оргнаизм, то
    добавить "|" в начале и "[ORGN]" в конце. Пример: lung AND macrophages | Rattus norvegicus[ORGN]
    """
    # -----------!!!! НАЧАЛО ОСНОВНОГО КОДА !!!!------------
    # проверка что оба мода не вызваны одновременно
    global Cell_size_for_tsv
    Cell_size_for_tsv = cell_size
    maxterms = 1000000  # формально нужно оганичение, но по факту смотрю все
    if verbose:
        global VerboseG
        VerboseG = True  # переключаем глобальную переменную флажком, сделанно для удобства написания кода
    if mode1 and mode2:
        return print("Error! Can not call mode1 and mode2 in same time!", file=sys.stderr)
    output_dir = os.path.join(dirname, output)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)  # проверяю наличие output директории, если ее нет то создаю
    # Блок инициации работы
    if (mode1 and not mode2) or (not mode1 and not mode2):
        print("Starting systematic search", file=sys.stderr)
        print("GEO search", file=sys.stderr)
        sist_search(input_file, maxterms, output_dir)  # Производим запросы, генерим файл input_GSE.txt
        print("ArrayExpress search", file=sys.stderr)
        arex_search(input_file, output_dir)
        print("End of systematic search", file=sys.stderr)
    main_true = (not mode1 and mode2) or (not mode1 and not mode2)
    if main_true:
        print("Starting main analysis", file=sys.stderr)
        # Меняю название GSE_id в зависимости от флагов
        # название для ArrayExpress пока не трогаю, пусть будет только дефолтное
        if mode2 and input_file != "input_terms.txt":
            input_GSE_name = os.path.join(dirname, input_file)
        else:
            input_GSE_name = os.path.join(output_dir, "input_GSE.txt")
        # открываем файл на чтение, считываем ID
        with open(input_GSE_name, 'r') as input_GSE:
            gse_list = input_GSE.readlines()
        # "чистим" названия, удаляя скрытый символ переноса строк
        gse_list = [x.strip() for x in gse_list]
        # тоже и для листа ArreyExpress
        with open(os.path.join(output_dir, "input_ArEx.txt"), 'r') as input_arex:
            arex_list = input_arex.readlines()
        arex_list = [x.strip() for x in arex_list]

        # Для удобства ввел переменную, чтоб оформление кусков не отличалось от text_out
        tab_out = True

        # Открываю файлы в соответсвии с флагами + сразу записываю шапки, если они нужны
        if tab_out:
            output_table = open(os.path.join(output_dir, "output_argeos.tsv"), 'w', encoding='utf-8')
            # Записываем шапку в таблицу
            output_table.write("Accession" + "\t" + "Organism" + "\t" + "Samples" + "\t" + "Type" + "\t" +
                               "Platform" + "\t" + "Title" +  # "\t" +"Authors" +
                               "\t" + "Year" + "\t" + "Summary" + "\t" + "Link" + "\t" + "Paper_title" + "\t" +
                               "Journal" + "\t" + "Impact factor 2018" + "\t" + "doi or pubmed id" + "\t" +
                               "All references" + "\t" + "Type of molecule" + "\t" + "BioProject link (NCBI)" +
                               "\t" + "BioProject link (EBI)" + "\t" + "SRA" + "\t" + "All protocols" + "\n")
        if text_out:
            output_file = open(os.path.join(output_dir, "output_argeos.txt"), 'w', encoding='utf-8')
        errors = open(os.path.join(output_dir, "errors_argeos.txt"), 'w', encoding='utf-8')

    # ---!! Начало ЦИКЛА!!---
    if main_true:
        # БЛОК GEO
        print("GEO datasets:", file=sys.stderr)
        # разбиваем наш лист на множество мелких
        gse_list_withou_bad = [y for y in gse_list if y not in Bad_IDs]
        gse_mega_list = chunks(gse_list_withou_bad, chunk_size)
        for listochek in tqdm(gse_mega_list):
            # БЛОК СКАЧИВАНИЯ ФАЙЛА
            # В начале я скачиваю пакет из нескольких листов, а потом уже спокойно прохожусь по этим файлам,
            # т.к. они уже на серваке (на компе) в tmp-шной папке
            # создание/стирание tmp директории
            try:
                shutil.rmtree(tmp_dir)  # если директория была до этого, то стираю ее
            except Exception:
                pass
            os.mkdir(tmp_dir)  # создаю заведомо  пустую директорию
            listochek = xml_by_id(listochek, errors)  # скачиваю файлы
            gse_list = []
            pubmed_id_list = []
            for GSE in listochek:
                # ! Блок АНАЛИЗА ФАЙЛА
                # парсю .xml файл, получаю всю инфу в лист
                if VerboseG:
                    print("working with " + GSE, file=sys.stderr)
                try:
                    gse_info = geo_xml_parser(os.path.join(tmp_dir, str(GSE + "_family.xml")))
                except Exception:
                    gse_info = GseInfo()
                    errors.write("error was (bad XML file) " + str(GSE) + "\n")
                gse_list.append(gse_info)
                if gse_info.PubMed_info.pbid is not None:
                    for pbid in gse_info.PubMed_info.pbid:
                        pubmed_id_list.append(pbid)
            gse_list = pubmed_parser(gse_list, pubmed_id_list,
                                     errors)  # произвожу анализ pubmed данных, обновляю лист этой инфой

            # ! Блок c разбиением на несколько строк
            if not unique:
                gse_list = split_to_unique(gse_list)
            # ! Блок записи результата
            if tab_out:
                table_output(gse_list, output_table)
            if text_out:
                text_output(gse_list, output_file)
        # БЛОК ArreyExpress
        chunk_size_lockal = 1  # Пока выключаю проход по кускам, будем индивидуально проходить
        # (иначе есть баг с записью в лист :( )
        ae_mega_list = chunks(arex_list, chunk_size_lockal)
        # Тут какая-та очень странная ошибка, вроде обращение к одному элементу.
        # Короче прохожусь по листу, а всесь лист получается заполнен копией последней переменной
        # пробовал deepcopy (даже функцию свою написал, для своего класса), но все равно не работает :(
        print("End of GEO. Starting ArrayExpress", file=sys.stderr)
        for listochek in tqdm(ae_mega_list):
            try:
                shutil.rmtree(tmp_dir)  # если директория была до этого, то стираю ее
            except Exception:
                pass
            os.mkdir(tmp_dir)  # создаю заведомо пустую tmp директорию

            pubmed_title_list = []
            ae_list = []
            # i = 0
            for aeid in listochek:
                ae_mega = array_express(aeid)
                ae_list.append(mydeepcopy(ae_mega))
                # ae_list.append(array_express(aeid))
                # i = i + 1
                if type(ae_mega.PubMed_info.title) == list:
                    for pbtit in ae_mega.PubMed_info.title:
                        pubmed_title_list.append(pbtit)
            if len(pubmed_title_list) > 0:
                tit_res = pbid_by_title(pubmed_title_list)
                tit_dict = tit_res[0]
                pb_id_list_ae = tit_res[1]
                # преобразование по словарю:
                for aesh in ae_list:
                    if type(aesh.PubMed_info.title) == list:
                        aesh.PubMed_info.pbid = []
                        all_fails = True
                        some_fails = False
                        for lt in aesh.PubMed_info.title:
                            pid = tit_dict.get(lt.strip("."))
                            if pid is not None:
                                aesh.PubMed_info.pbid.append(pid)
                                all_fails = False
                            else:
                                some_fails = True
                        if all_fails and some_fails:
                            aesh.PubMed_info.pbid = None
                ae_list = pubmed_parser(ae_list, pb_id_list_ae, errors)

            # ! Блок записи результата
            if tab_out:
                table_output(ae_list, output_table)
            if text_out:
                text_output(ae_list, output_file)

    # ---!! КОНЕЦ ЦИКЛА !!---

    # Блок терминации работы
    # Закрываю файлы на запись, удаляю tmp директорию
    if main_true:
        errors.close()
        if tab_out:
            output_table.close()
        if text_out:
            output_file.close()
        try:
            shutil.rmtree(tmp_dir)  # удаление tmp директории
        except Exception:
            pass
        print("Work finished!", file=sys.stderr)
    # Конец основного кода


if __name__ == "__main__":
    main()
