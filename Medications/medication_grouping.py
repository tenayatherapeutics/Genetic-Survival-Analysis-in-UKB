#!/usr/bin/python3

import sys 
import subprocess
import os 
import re

med_groups = 'MedicationGroups.txt'
codings = 'coding4.tsv'

def create_dict(filepath):
    my_dict = {}
    with open(filepath, 'rt') as file: 
        ln = file.readline() 
        ln = file.readline() 
        while ln: 
            cols = ln.strip('\n').split('\t')
            key = cols[0].strip('\n').strip(' ')
            value = cols[1].strip('\n').strip(' ')
            if(key != ''):
                my_dict[key] = value 
            ln = file.readline()
    return my_dict

def clean_medi(medi_name): 
    medi_name = medi_name.strip('\"').lower()
    word_list = re.split(', |/| |-', medi_name)
    for word in word_list:
        if len(word) < 4:
            word_list.remove(word)
    return word_list

def clean_drug(drug_name):
    exclusion_list = ['alpha', 'flu', 'continus', 'plus', 'tablet']
    drug_name = drug_name.strip('\"').lower()
    word_list = re.split('-| |/|\\+', drug_name)
    name = word_list[0]
    for i in range(len(word_list)):
        if len(name) < 4 or name in exclusion_list:
            name = word_list[i]
    return(name)
        
med_groups_dict = create_dict(med_groups)
codings_dict = create_dict(codings)
#print(med_groups_dict.keys())
#print(list(codings_dict.keys())[0:10])

merged_list = []

for code in codings_dict: 
    for group in med_groups_dict: 
        medi_list = clean_medi(med_groups_dict[group])
        medi_word = clean_drug(codings_dict[code])
        #if any(s.startswith(medi_word) for s in medi_list):
        #    merged_list.append([code, codings_dict[code], group])
        if any(s == (medi_word) for s in medi_list):
            merged_list.append([code, codings_dict[code], group])

print(len(merged_list))
print(merged_list[0:10])

with open('merged_med_groups_stringent_anticoag.txt', 'w') as file:
    file.write('\t'.join(['Coding', 'Drug', 'Group']) + '\n')
    for i in merged_list:
        file.write('\t'.join(i) + '\n')
