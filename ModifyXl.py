# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 09:53:29 2016

@author: Wesley
"""

import pandas as pd
import numpy as np
#from pandas import ExcelWriter
from isotopomer_distribution_correction import parse_formula, Fernandez1996_correction
import metabolite_structure
import sys, traceback
#import re

def calc_fractional_labeling(molecular_formula, corrected_distribution):
    if molecular_formula:
        num_carbons = parse_formula(molecular_formula)['C']
        #num_carbons = 3
        numerator = 0
        denominator = 0
        for m in range(num_carbons + 1):
            a = corrected_distribution[m]
            numerator += (m * a)
            denominator += a
        denominator *= num_carbons
        fractional_labeling = numerator / denominator
    else:
        fractional_labeling = 0.
    return fractional_labeling
print("Starting...")
#fileName = 'ExcelExp_Short.xlsx'
#fileName = 'MCF028_AorticStenosis_Short.XLS'
#fileName = "C:\\Users\\u0062868\\Documents\\Eleonora\\MCF038_secondFocus_Short_WV.XLS"
fileName = "I:\\Wesley\\PaulineDezeeuw\\MCF3067_3034_13C_rest_Short.XLS"
#cuttoff = 2000000
cuttoff = 0


print("File Found...")
xls = pd.ExcelFile(fileName)
#xls = pd.ExcelFile(fileName)
fileName = fileName.split('.')[0]
shs = xls.sheet_names
pn = pd.Panel()
df = pd.DataFrame()
df_RT = pd.DataFrame()
for name in shs:
    if '_m' in name:
        comp = name.rsplit('_', 1)[0]
        mass = name.rsplit('_', 1)[1]
        sheet = xls.parse(name,header=4)
        index = np.where(sheet['Filename'].isnull())[0][0]
        newSheet = sheet[:index]
        newestSheet= newSheet[['Filename','Area','RT']]
        
        frame = pd.DataFrame(newestSheet)
        df[name] = newestSheet['Area'].replace('NF',float(0))
        df_RT[name] = newestSheet['RT'].replace('NF',float(0))
        
print("Processing File...")
df = df.set_index(newestSheet['Filename'])
df = df.sort_index(axis=1)
#print('ok')
compounds = list(df.columns.values)
#print(compounds)
df_perc = pd.DataFrame(index = df.index.values.tolist())
df_corr = pd.DataFrame(index = df.index.values.tolist())
df_fractLabel = pd.DataFrame(index = df.index.values.tolist())
sumAbundance = pd.DataFrame(index = df.index.values.tolist())
list_compounds = []
for comp in compounds :
    comp_name = comp.rsplit('_', 1)[0]
    list_compounds.append(comp_name)
list_uniq_compounds = list(set(list_compounds))
list_uniq_compounds = sorted(list_uniq_compounds)

for uniq in list_uniq_compounds :
#    now = df[df.columns[df.columns.to_series().str.contains(uniq)]]
    pattern = uniq + r'_m[0-9]{1,2}'
#    if '+' in pattern :
#        pattern.replace('+','\\\+')
#        print("found one")
    now=df.filter(regex=pattern)
    #print(pattern)
    df_ab = pd.DataFrame(data=now.sum(axis=1))
    try :
        df_ab.columns = ['_'.join(now.columns.values[0].split('_')[:-1])]
        sumAbundance = pd.concat([sumAbundance,df_ab],axis = 1)
    #now_max = now.max(axis=1)
        now_sum = now.sum(axis=1)
        now[now<cuttoff] = float(0)
    #now = now.div(now_max,axis=0)
        now = now.div(now_sum,axis=0)
        df_perc = pd.concat([df_perc, now], axis=1)
    except :
        print(now.columns)
        print("Unexpected error:", sys.exc_info()[0])
    
    #if '_' in uniq :
    #uniq = uniq.split('_')[0]
#        print(uniq)
    try :
        formula = metabolite_structure.metabolite[uniq]
        print(uniq+": "+formula+ "  -->  " + metabolite_structure.MEC_Synonyms.get(uniq,uniq))
        df_corrected = pd.DataFrame(index = df.index.values.tolist(),columns=now.columns.values)
        
        df_fl = pd.DataFrame(index = df.index.values.tolist(),columns=now.columns.values)
        #df_ab = pd.DataFrame(index = df.index.values.tolist(),columns=now.columns.values)
        try:
            for index, row in now.iterrows() : #index: filename, row = values & columnheads
                linelist = list(row)
                distribution = [float(v) for v in linelist]
                corrected_distribution = list(Fernandez1996_correction(formula, distribution))
                mci_or = np.max(distribution)
                distribution = [ci/mci_or for ci in distribution]
                sumCorr = np.sum(corrected_distribution)
                corrected_distribution = [ci/sumCorr for ci in corrected_distribution]
                fractional_labeling = calc_fractional_labeling(formula, distribution)
                df_fl.ix[index] = fractional_labeling
                #df_fl.columns = [df.columns.values[0].split('_')[0]]
                #df_corrected = pd.DataFrame(index = df.index.values.tolist(),columns=now.columns.values)
                #try :
                df_corrected.ix[index] = corrected_distribution
                #except :
                    #print("No C13 tracing")
        except:
            print("No C13 tracing")
            
                
        df_corr = pd.concat([df_corr, df_corrected], axis=1)
        df_fl_finish = pd.DataFrame(df_fl.ix[:,0])
        df_fl_finish.columns = ['_'.join(df_fl.columns.values[0].split('_')[:-1])]
        df_fractLabel = pd.concat([df_fractLabel, df_fl_finish], axis = 1)
    except ValueError:
        print("Could not convert data to an integer or there are two files named the same...")
        traceback.print_exc(file=sys.stdout)
    except :
        print(uniq + " not found: Is the name correct in metabolite_structure? or error during calculations...")
        print("Unexpected error:", sys.exc_info()[0])
        
        
#gly = [x.replace(' ','_') for x in metabolite_structure.Glycolysis]
#kreb = [x.replace(' ','_') for x in metabolite_structure.Krebs]
#aa = [x.replace(' ','_') for x in metabolite_structure.AminoAcids]
#ur = [x.replace(' ','_') for x in metabolite_structure.Urea]
#penpho = [x.replace(' ','_') for x in metabolite_structure.PentPhos]

df_GlyKreb = pd.concat([df_fractLabel.ix[:,metabolite_structure.MEC_Glycolysis],df_fractLabel.ix[:,metabolite_structure.MEC_Citric]],axis=1)
df_GlyKreb.columns = [metabolite_structure.MEC_Synonyms.get(x,x) for x in df_GlyKreb.columns.values]
df_Amino = df_fractLabel.ix[:,metabolite_structure.MEC_AminoAcids]
df_Amino.columns = [metabolite_structure.MEC_Synonyms.get(x,x) for x in df_Amino.columns.values]
df_PentUrea = pd.concat([df_fractLabel.ix[:,metabolite_structure.MEC_PentPhos],df_fractLabel.ix[:,metabolite_structure.MEC_Urea]],axis=1)
df_PentUrea.columns = [metabolite_structure.MEC_Synonyms.get(x,x) for x in df_PentUrea.columns.values]
df_Energy = df_fractLabel.ix[:,metabolite_structure.SpecialNucleotides]
df_Energy.columns = [metabolite_structure.MEC_Synonyms.get(x,x) for x in df_Energy.columns.values]

fileName = fileName + '_Converted.xlsx'
#fileName = fileName + '_Cut.xlsx'

def multiple_dfs(df_list, sheets, file_name, spaces):
    writer = pd.ExcelWriter(file_name,engine='xlsxwriter')   
    row = 0
    for dataframe in df_list:
        dataframe.to_excel(writer,sheet_name=sheets,startrow=row , startcol=0)   
        row = row + len(dataframe.index) + spaces + 1
    df_fractLabel.columns = [metabolite_structure.MEC_Synonyms.get(x,x) for x in df_fractLabel.columns.values]
    df_fractLabel.to_excel(writer,'FracContribution')
    df.to_excel(writer,'rawResult')
    sumAbundance.columns = [metabolite_structure.MEC_Synonyms.get(x,x) for x in sumAbundance.columns.values]
    sumAbundance.to_excel(writer,'rawAbundances')
    df_perc.to_excel(writer,'rawPercentage')
    df_corr.to_excel(writer,'Corrected&FracLabel')
    writer.save()

# list of dataframes
dfs = [df_GlyKreb,df_Amino,df_PentUrea]

# run function
multiple_dfs(dfs, 'OrderedResults', fileName, 2)
#writer.save()
print('finished')
