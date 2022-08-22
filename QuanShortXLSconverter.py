# -*- coding: utf-8 -*-

"""
Created on Fri Oct  7 09:53:29 2016

@author: Wesley
"""
import pandas as pd
import numpy as np
from isotopomer_distribution_correction import parse_formula, Fernandez1996_correction
import metabolite_structure
import re
import sys, getopt, os, traceback

np.seterr(divide='ignore', invalid='ignore')

def calc_fractional_labeling(molecular_formula, corrected_distribution, element):
    if molecular_formula:
        num_carbons = parse_formula(molecular_formula)[element]
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
    
def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print ('python QuanShortXLSconverter.py -i <inputfile> -o <outputfile,optional>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('python QuanShortXLSconverter.py -i <inputfile> -o <outputfile,optional>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    #fileName = "I:\\Wesley\\TestFolder\\testing.xlsx"
    if len(sys.argv) <= 1:
        print('You have to tell me where the excel file is, stoopy... like this: \r\n \r\n>python QuanShortXLSconverter.py -i <inputfile> -o <outputfile,optional>')
        sys.exit(1)
    fileName = inputfile
    dir_path = os.path.dirname(os.path.realpath(inputfile))
    if outputfile == '': outputfile = inputfile.split('.')[0]+ '_Converted.xlsx'
    elif '\\' not in outputfile : 
        if '.xlsx' not in outputfile:
            outputfile = dir_path+'\\'+outputfile + '.xlsx'
        else:
            outputfile = dir_path+'\\'+outputfile
    
    
    print("File Found...")
    xls = pd.ExcelFile(fileName)
    fileName = fileName.split('.')[0]
    shs = xls.sheet_names
    if 'Component' in shs:
        shs.remove('Component')
    df = pd.DataFrame()
    df_RT = pd.DataFrame()
    for name in shs:
        if '_m0' or '_m00'  in name:
            CompleteSheetFromExcel = xls.parse(name,header=4)
            index_NameOfFiles = np.where(CompleteSheetFromExcel['Filename'].isnull())[0][0]
            ReducedSheetOnlyFileNames = CompleteSheetFromExcel[:index_NameOfFiles]
            SheetFileNameAreaRT= ReducedSheetOnlyFileNames[['Filename','Area','RT']]
            df[name] = SheetFileNameAreaRT['Area'].replace('NF',float(0))
            df_RT[name] = SheetFileNameAreaRT['RT'].replace('NF',float(0))
            
    print("Processing File...")
    if any(SheetFileNameAreaRT['Filename'].duplicated(keep=False)):
        print()
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("!!!WARNING: There are duplicate rows in your excel file, RAW-files with the same name!!!")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print()
        print("Removing one of them")
        print()
        df = df.set_index(SheetFileNameAreaRT['Filename'])
        df = df.reset_index().drop_duplicates(subset='Filename', keep='first').set_index('Filename')
        df = df.drop_duplicates(keep="last")
    else:
        df = df.set_index(SheetFileNameAreaRT['Filename'])
    df = df.sort_index(axis=1)
    compoundsWithIsotopologues = list(df.columns.values)
    df_rawValuesNoChange = pd.DataFrame(index = df.index.values.tolist())
    df_raw_M_Complete = pd.DataFrame(index = df.index.values.tolist())
    df_raw_N_Complete = pd.DataFrame(index = df.index.values.tolist())
    df_percentageM_Complete = pd.DataFrame(index = df.index.values.tolist())
    df_percentageN_Complete = pd.DataFrame(index = df.index.values.tolist())
    df_corrected_Complete = pd.DataFrame(index = df.index.values.tolist())
    df_fractContribution = pd.DataFrame(index = df.index.values.tolist())
    df_fractContributionN = pd.DataFrame(index = df.index.values.tolist())
    df_sumAbundances_Complete = pd.DataFrame(index = df.index.values.tolist())
#    tracers = False
    list_compounds = []
    for compoundWIT in compoundsWithIsotopologues :
#        tracerpattern = r'[A-Z]*[a-z]*_m0?1n?0?$'
#        if(re.match(compoundWIT,tracerpattern)):
#            tracers = True
        compound_name = compoundWIT.rsplit('_', 1)[0]
        list_compounds.append(compound_name)
        
    list_uniq_compounds = list(set(list_compounds))
    list_uniq_compounds = sorted(list_uniq_compounds)
    #print(list_uniq_compounds)
#    if(tracers):
#        print("Tracers were found")
#    else:
#        print("No tracers were found")
    for uniq in list_uniq_compounds :
        try:
            formula = metabolite_structure.metabolite[uniq]
            print(uniq+": "+formula+ "  -->  Found")
            if uniq  + '_m0n0' in compoundsWithIsotopologues or uniq+'_m00n0' in compoundsWithIsotopologues:
                patternN = '^' + uniq + r'_m[0-9]{1,2}n[0-9]{1,2}$'
                UniqueCompoundTable = df.filter(regex=patternN)
                numberOfIndividualElementsInFormula = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
                if not numberOfIndividualElementsInFormula [0][1]:
                    Carbon = 1
                else:
                    Carbon = int(numberOfIndividualElementsInFormula [0][1])
                if not numberOfIndividualElementsInFormula [2][1]:
                    Nitrogen = 1
                else:
                    Nitrogen = int(numberOfIndividualElementsInFormula [2][1])
                df_M = pd.DataFrame(index = df.index.values.tolist())
                df_N = pd.DataFrame(index = df.index.values.tolist())
                
                for m in range(Carbon+1):
                    if m<10:
                        if Carbon<10:
                            temperaryPatternC = '^'+uniq+r'_m'+str(m)+'n[0-9]$'
                            df_M[uniq+'_m'+str(m)] = UniqueCompoundTable.filter(regex=temperaryPatternC).sum(axis=1)
                        else:
                            temperaryPatternC = '^'+uniq+r'_m0'+str(m)+'n[0-9]$'
                            df_M[uniq+'_m0'+str(m)] = UniqueCompoundTable.filter(regex=temperaryPatternC).sum(axis=1)
                    else:
                        temperaryPatternC = '^'+uniq+r'_m'+str(m)+'n[0-9]$'
                        df_M[uniq+'_m'+str(m)] = UniqueCompoundTable.filter(regex=temperaryPatternC).sum(axis=1)
                df_raw_M_Complete = pd.concat([df_raw_M_Complete,df_M],axis = 1)
                
                for n in range(Nitrogen+1):
                    temperaryPatternN = '^'+uniq+r'_m[0-9]{1,2}n'+str(n)+'$'
                    df_N[uniq+'_n'+str(n)] = UniqueCompoundTable.filter(regex=temperaryPatternN).sum(axis=1)
                df_raw_N_Complete = pd.concat([df_raw_N_Complete,df_N],axis = 1)
                df_Perc_N = pd.DataFrame(index = df.index.values.tolist())
                df_Perc_N = df_N.div(df_N.sum(axis=1),axis=0)
                df_percentageN_Complete = pd.concat([df_percentageN_Complete,df_Perc_N],axis = 1)
                df_FrConN = pd.DataFrame(index = df.index.values.tolist(),columns=df_Perc_N.columns.values)
                try:
                    print("Start Correction for Nitrogen")
                    for index, row in df_Perc_N.iterrows() : #index: filename, row = values & columnheads
                        linelist = list(row)
                        distribution = [float(v) for v in linelist]
                        corrected_distribution = distribution
                        mci_or = np.max(distribution)
                        distribution = [ci/mci_or for ci in distribution]
                        sumCorr = np.sum(corrected_distribution)
                        corrected_distribution = [ci/sumCorr for ci in corrected_distribution]
                        fractional_labeling = calc_fractional_labeling(formula, distribution, 'N')
                        df_FrConN.loc[index] = fractional_labeling
                except ValueError:
                    print("\tError during the correction of isotopologues (ValueError)...")
                except KeyError:
                    print("Not good --> KeyError")			
                except:
                    print("\t\tI cannnot calculate the fractional Contribution for Nitrogen!")
                    print("Unexpected error:", sys.exc_info()[0])
                
                df_fConN_finish = pd.DataFrame(df_FrConN.iloc[:,0])
                df_fConN_finish.columns = ['_'.join(df_FrConN.columns.values[0].split('_')[:-1])]
                df_fractContributionN = pd.concat([df_fractContributionN, df_fConN_finish], axis = 1)
            else:
                patternC = '^' + uniq + r'_m[0-9]{1,2}$'
                UniqueCompoundTable = df.filter(regex=patternC)
                numberOfIndividualElementsInFormula = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
                Carbon = int(numberOfIndividualElementsInFormula [0][1])
                #print("C: "+Carbon)
                df_M = UniqueCompoundTable
            
        except TypeError:
             print("TypeError:something wrong with higher Carbon count: " + str(Carbon))
        except ValueError:
            print("ValueError: Can't find a value-> ", sys.exc_info()[0])
             
        except:
             print(uniq +": the compound name does not exist in the databank, so I cannot find the formula")
             print('_' * (len(uniq)+1) + "  Unexpected error:", sys.exc_info()[0])
      
        try :
            df_rawValuesNoChange = pd.concat([UniqueCompoundTable,df_rawValuesNoChange],axis = 1)
            
            df_SumCompoundIsotopologues = pd.DataFrame(data=UniqueCompoundTable.sum(axis=1),columns=[uniq])
            df_sumAbundances_Complete = pd.concat([df_sumAbundances_Complete,df_SumCompoundIsotopologues],axis = 1)
            
            df_Perc_M = pd.DataFrame(index = df.index.values.tolist())
            df_Perc_M = df_M.div(df_M.sum(axis=1),axis=0)
            df_percentageM_Complete = pd.concat([df_percentageM_Complete,df_Perc_M],axis = 1)
        except :
            print("\tSomething wrong with the calculations, and it ain't me...\n\tError trying to add up numbers for RawAbundance for Carbons...")
#        
        try :
            df_corrected = pd.DataFrame(index = df.index.values.tolist(),columns=df_Perc_M.columns.values)
            
            df_fl = pd.DataFrame(index = df.index.values.tolist(),columns=df_Perc_M.columns.values)
            try:
                print("Start Correction for Carbon")
                for index, row in df_Perc_M.iterrows() : #index: filename, row = values & columnheads
                    linelist = list(row)
                    distribution = [float(v) for v in linelist]
                    corrected_distribution = list(Fernandez1996_correction(formula, distribution))
                    mci_or = np.max(distribution)
                    distribution = [ci/mci_or for ci in distribution]
                    sumCorr = np.sum(corrected_distribution)
                    corrected_distribution = [ci/sumCorr for ci in corrected_distribution]
                    fractional_labeling = calc_fractional_labeling(formula, distribution,'C')
                    df_fl.loc[index] = fractional_labeling
                    df_corrected.loc[index] = corrected_distribution
            except ValueError:
                print("\tError during the correction of isotopologues or no tracers found...")
            except KeyError:
                print("Not good")				
            except:
                print("\t\tI cannnot calculate the fractional Contributrion for the Carbons!")
            df_corrected_Complete = pd.concat([df_corrected_Complete, df_corrected], axis=1)
            df_fl_finish = pd.DataFrame(df_fl.iloc[:,0])
            df_fl_finish.columns = ['_'.join(df_fl.columns.values[0].split('_')[:-1])]
            df_fractContribution = pd.concat([df_fractContribution, df_fl_finish], axis = 1)
            print("End Compound")
        except ValueError:
            print("Could not convert data to an integer or there are two files named the same...")
            traceback.print_exc(file=sys.stdout)
        except :
            print(uniq + " not found: Is the name correct in metabolite_structure? or error during calculations...")
            print("Unexpected error:", sys.exc_info()[0])
            
            
    try:

        writer = pd.ExcelWriter(outputfile,engine='xlsxwriter')
        if not df_rawValuesNoChange.empty:
            df_rawValuesNoChange.reindex(sorted(df.columns), axis=1).to_excel(writer,'rawResult')
            #df_rawValuesNoChange.to_excel(writer,'rawResult')
        if not df_sumAbundances_Complete.empty:
            df_sumAbundances_Complete.to_excel(writer,'rawAbundances')
        if not df_raw_M_Complete.empty:
            df_raw_M_Complete.to_excel(writer,'RawCarbons')
        if not df_percentageM_Complete.empty:
            df_percentageM_Complete.to_excel(writer,'PercentageCarbons')
        if not df_corrected_Complete.empty:
            df_corrected_Complete.to_excel(writer,'CorrectedCarbonsOnly')
        if not df_fractContribution.empty:
            df_fractContribution.to_excel(writer,'FracContributionCarbons')
        if not df_raw_N_Complete.empty:
            df_raw_N_Complete.to_excel(writer,'RawNitrogens')
        if not df_percentageN_Complete.empty:
            df_percentageN_Complete.to_excel(writer,'PercentageNitrogens')
        if not df_fractContributionN.empty:
            df_fractContributionN.to_excel(writer,'FracContributionNitrogen')
        #writer.save()
    except ValueError as v:
        print()
        print("Value Error: " , v)
        print("First tab won't be alphabetized...  why?   good question...")
        print("Maybe check for missing compound and compound formulas in the list above.")
        print("Add that missing compound with formula to the metabolite_structure.py file.")
        df_rawValuesNoChange.to_excel(writer,'rawResult')
        
    except :
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("Unexpected error:", sys.exc_info()[0])
        
    else :
        print()
        print("No major problems detected :)")
        print(outputfile)
        
    finally :
        print()
        print('finished')
        
        writer.save()
    
if __name__ == "__main__":
    main(sys.argv[1:])
