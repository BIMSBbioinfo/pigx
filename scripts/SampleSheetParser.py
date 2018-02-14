#!@PYTHON@
# -*- python -*-
# PiGx ChIPseq Pipeline.
#
# Copyright © 2018 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
#
# This file is part of the PiGx ChIPseq Pipeline.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse

import xlrd
import re
import yaml
import csv

# "PigX_ChIP sample_sheet template.excel"


def load_book(file_name):
    book = xlrd.open_workbook(file_name)
    return book


def excel_sheet_to_dict(book, sheet_name="sample_sheet"):
    sheet = book.sheet_by_name(sheet_name)
    rows = [sheet.row_values(r) for r in range(0, sheet.nrows)]
    header = rows[0]
    rows = rows[1:]
    return [dict(zip(header, row)) for row in rows]


def csv_file_to_dict(file_name):
#   # Load sample sheet
    with open(file_name, 'r') as fp:
        rows = [row for row in csv.reader(fp, delimiter=',')]
        header = rows[0]
        rows = rows[1:]
        rows = [row for row in rows if row[0]]
        return [dict(zip(header, row)) for row in rows]

#
# MAPPING_DICT = sheet_to_dict("mapping")
# PEAK_CALLING_DICT = sheet_to_dict("peak_calling")


def split_list(string_of_list):
    return string_of_list.split(";")


def split_pair(string_of_pair):
    return string_of_pair.split(",")


def split_keyValue_to_dict(string_of_keyValue):
    pair = re.split("\W+", string_of_keyValue.strip())
    if len(pair) < 2:
        pair.append('')
    return dict([pair])


def params_to_dict(sample_dict):
        tool_params = [param
                       for param
                       in sample_dict.keys()
                       if "params" in param]
        parDict = {}
        for par in tool_params:
            parStr, tool = re.split("[\W_]", par)
#            print("Found Parameters for " + tool + ".")
            toolDict = {}

            setting = sample_dict[par]
            d = {}
            for s in split_list(setting):
                d.update(split_keyValue_to_dict(s))
            toolDict[tool] = d
            parDict.update(toolDict)
            sample_dict.pop(par)
        if tool_params:
            sample_dict.update({parStr: parDict})
#            print("updated sample dictionary")


def parse_single_sample(sample_row):
    if isinstance(sample_row["ChIP"], str):
        ChIP_list = [split_pair(reads)
                     for reads
                     in split_list(sample_row["ChIP"])]
        sample_row["ChIP"] = ChIP_list
    if isinstance(sample_row["Cont"], str):
        Cont_list = [split_pair(reads)
                     for reads
                     in split_list(sample_row["Cont"])]
        sample_row["Cont"] = Cont_list
    params_to_dict(sample_row)
    return sample_row


def parse_sample_sheet(dict):
    sample_sheet = [parse_single_sample(sample_row)
                    for sample_row
                    in dict]
    return sample_sheet


def strip_ext(fastq_file):
    for ext in [".fq", ".fastq"]:
        if ext in fastq_file:
            fastq_file = fastq_file[:fastq_file.find(ext)]
    return fastq_file


def strip_pair(fastq_file):
    for pair in ["_1", "_2", "_R1", "_R2"]:
        if pair in fastq_file:
            fastq_file = fastq_file[:fastq_file.find(pair)]
    return fastq_file


def get_sample_names(samples):
    sample_names = [list(set(map(lambda x: strip_pair(strip_ext(x)), files)))
                    if len(files) == 2
                    else [strip_ext(files[0])]
                    for files in samples]
    sample_names = sum(sample_names, [])
    return sample_names


def get_lib_type(samples):
    lib_type = ["paired" if len(files) == 2
                else "single"
                for files in samples]
    return lib_type


def define_mapping(sample_dict):
    samples = [dict['ChIP'][0] for dict in sample_dict]
    samples += (dict['Cont'][0] for dict in sample_dict)
    # remove duplicates
    samples = [list(x) for x in set(tuple(x) for x in samples)]
    # remove empty strings
    samples = list(filter(None, [list(filter(None, file))
                                 for file in samples]))
    lib_type = get_lib_type(samples)
    sample_names = get_sample_names(samples)
    sample_names = ["{}.PE".format(sample)
                    if lib == "paired"
                    else sample
                    for sample, lib
                    in zip(sample_names, lib_type)]
    samples_dict = {}
    for sample, file, lib in zip(sample_names, samples, lib_type):
        samples_dict[sample] = {"fastq": file, "library": lib}
    return samples_dict


def define_peakCalling(sample_dict, mapping_dict):
    peak_dict = {}
    for dict in sample_dict:
        peak_name = dict['Peak_Name']
        peak_dict[peak_name] = {}
        chip_list = []
        cont_list = []
        for files in dict["ChIP"]:
            for key in mapping_dict.keys():
                fastq_files = mapping_dict[key]["fastq"]
                if len(dict["ChIP"]) == 1:
                    if fastq_files == files:
                        chip_list = key
                        continue
                else:
                    if fastq_files == files:
                        chip_list.append(key)
        for files in dict["Cont"]:
            for key in mapping_dict.keys():
                if len(dict["Cont"]) == 1:
                    if fastq_files == files:
                        cont_list = key
                        continue
                else:
                    if fastq_files == files:
                        cont_list.append(key)
        peak_dict[peak_name]["ChIP"] = chip_list
        peak_dict[peak_name]["Cont"] = cont_list
        peak_dict[peak_name]["params"] = dict["params"]
    return(peak_dict)


def define_idr(sample_dict):
    idr_dict = {}
    idr_dict = {}
    idr_list = [dict['idr'] for dict in sample_dict]
    idr_names = list(filter(None, set(idr_list)))
    for name in idr_names:
        if idr_list.count(name) == 2:
            peaks = [dict['Peak_Name']
                     for dict in sample_dict
                     if dict['idr'] in name]
            idr_dict.update({name: {"ChIP1": peaks[0],
                                    "ChIP2": peaks[1]}})
    return idr_dict


def define_feature(sample_dict, feature_name):
    feature_dict = {}
    feature_dict = {}
    feature_list = [dict[feature_name] for dict in sample_dict]
    feature_names = list(filter(None, set(feature_list)))
    for name in feature_names:
        peaks = [dict['Peak_Name']
                 for dict in sample_dict
                 if dict[feature_name] in name]
        feature_dict[name] = peaks
    return feature_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f',
        '--file',
        required=True,
        help="The excel or csv file to load.")
    parser.add_argument(
        '-s',
        '--sample-sheet',
        default="sample_sheet",
        help='The name of "sample_sheet" sheet.')
    args = parser.parse_args()
    # book = load_book(args.file)
    # dict = sheet_to_dict(book, args.sample_sheet)
    # print(yaml.dump(dict))
    if args.file.endswith('.xlsx'):
        with load_book(args.file) as book:
            tmp_dict = excel_sheet_to_dict(book, args.sample_sheet)
    elif args.file.endswith('.csv'):
        tmp_dict = csv_file_to_dict(args.file)
    sample_dict = parse_sample_sheet(tmp_dict)

    sample_sheet_dict = {}
    sample_sheet_dict["samples"] = define_mapping(sample_dict)
    sample_sheet_dict["peak_calling"] = define_peakCalling(sample_dict, sample_sheet_dict["samples"])
    sample_sheet_dict["idr"] = define_idr(sample_dict)
    sample_sheet_dict["feature_combination"] = define_feature(sample_dict, "feature_combination")
    sample_sheet_dict["feature_factors"] = define_feature(sample_dict, "feature_factors")

    # print(json.dumps(sample_dict, indent=4))
    print("__________")
    print(yaml.dump(sample_sheet_dict))
