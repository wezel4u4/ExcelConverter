'''
Created on Sep 4, 2013

@author: Dries Verdegem
'''

import re
import numpy

# The isotopic mass data is from G. Audi, A. H. Wapstra: Nucl. Phys A. 1993, 565, 1-65 and G. Audi, A. H. Wapstra: Nucl. Phys A. 1995, 595, 409-480.

# The percent natural abundance data is from the 1997 report of the IUPAC Subcommittee for Isotopic Abundance Measurements
# by K.J.R. Rosman & P.D.P. Taylor: Pure Appl. Chem. 1999, 71, 1593-1607.
natural_abundances = {}

natural_abundances['H']  = ((1,1.007825,99.9885),(2,2.014102,0.0115),(3,3.016049,0.))
natural_abundances['He'] = ((3,3.016029,0.000137),(4,4.002603,99.999863))
natural_abundances['Li'] = ((6,6.015122,7.59),(7,7.016004,92.41))
natural_abundances['Be'] = ((9,9.012182,100),)
natural_abundances['B']  = ((10,10.012937,19.9),(11,11.009305,80.1))
natural_abundances['C']  = ((12,12.000000,98.93),(13,13.003355,1.07),(14,14.003242,0.))
natural_abundances['N']  = ((14,14.003074,99.632),(15,15.000109,0.368))
natural_abundances['O']  = ((16,15.994915,99.757),(17,16.999132,0.038),(18,17.999160,0.205))
natural_abundances['F']  = ((19,18.998403,100),)
natural_abundances['Ne'] = ((20,19.992440,90.48),(21,20.993847,0.27),(22,21.991386,9.25))
natural_abundances['Na'] = ((23,22.989770,100),)
natural_abundances['Mg'] = ((24,23.985042,78.99),(25,24.985837,10.00),(26,25.982593,11.01))
natural_abundances['Al'] = ((27,26.981538,100),)
natural_abundances['Si'] = ((28,27.976927,92.230),(29,28.976495,4.683),(30,29.973770,3.087))
natural_abundances['P']  = ((31,30.973762,100),)
natural_abundances['S']  = ((32,31.972071,94.93),(33,32.971458,0.76),(34,33.967867,4.29),(36,35.967081,0.02))
natural_abundances['Cl'] = ((35,34.968853,75.78),(37,36.965903,24.22))
natural_abundances['Ar'] = ((36,35.967546,0.3365),(38,37.962732,0.0632),(40,39.962383,99.6003))
natural_abundances['K']  = ((39,38.963707,93.2581),(40,39.963999,0.0117),(41,40.961826,6.7302))
natural_abundances['Ca'] = ((40,39.962591,96.941),(42,41.958618,0.647),(43,42.958767,0.135),(44,43.955481,2.086),(46,45.953693,0.004),(48,47.952534,0.187))
natural_abundances['Sc'] = ((45,44.955910,100),)
natural_abundances['Ti'] = ((46,45.952629,8.25),(47,46.951764,7.44),(48,47.947947,73.72),(49,48.947871,5.41),(50,49.944792,5.18))
natural_abundances['V']  = ((50,49.947163,0.250),(51,50.943964,99.750))
natural_abundances['Cr'] = ((50,49.946050,4.345),(52,51.940512,83.789),(53,52.940654,9.501),(54,53.938885,2.365))
natural_abundances['Mn'] = ((55,54.938050,100),)
natural_abundances['Fe'] = ((54,53.939615,5.845),(56,55.934942,91.754),(57,56.935399,2.119),(58,57.933280,0.282))
natural_abundances['Co'] = ((59,58.933200,100),)
natural_abundances['Ni'] = ((58,57.935348,68.0769),(60,59.930791,26.2231),(61,60.931060,1.1399),(62,61.928349,3.6345),(64,63.927970,0.9256))
natural_abundances['Cu'] = ((63,62.929601,69.17),(65,64.927794,30.83))
natural_abundances['Zn'] = ((64,63.929147,48.63),(66,65.926037,27.90),(67,66.927131,4.10),(68,67.924848,18.75),(70,69.925325,0.62))
natural_abundances['Ga'] = ((69,68.925581,60.108),(71,70.924705,39.892))
natural_abundances['Ge'] = ((70,69.924250,20.84),(72,71.922076,27.54),(73,72.923459,7.73),(74,73.921178,36.28),(76,75.921403,7.61))
natural_abundances['As'] = ((75,74.921596,100),)
natural_abundances['Se'] = ((74,73.922477,0.89),(76,75.919214,9.37),(77,76.919915,7.63),(78,77.917310,23.77),(80,79.916522,49.61),(82,81.916700,8.73))
natural_abundances['Br'] = ((79,78.918338,50.69),(81,80.916291,49.31))
natural_abundances['Kr'] = ((78,77.920386,0.35),(80,79.916378,2.28),(82,81.913485,11.58),(83,82.914136,11.49),(84,83.911507,57.00),(86,85.910610,17.30))
natural_abundances['Rb'] = ((85,84.911789,72.17),(87,86.909183,27.83))
natural_abundances['Sr'] = ((84,83.913425,0.56),(86,85.909262,9.86),(87,86.908879,7.00),(88,87.905614,82.58))
natural_abundances['Y']  = ((89,88.905848,100),)
natural_abundances['Zr'] = ((90,89.904704,51.45),(91,90.905645,11.22),(92,91.905040,17.15),(94,93.906316,17.38),(96,95.908276,2.80))
natural_abundances['Nb'] = ((93,92.906378,100),)
natural_abundances['Mo'] = ((92,91.906810,14.84),(94,93.905088,9.25),(95,94.905841,15.92),(96,95.904679,16.68),(97,96.906021,9.55),(98,97.905408,24.13),(100,99.907477,9.63))
natural_abundances['Tc'] = ((98,97.907216,100),)
natural_abundances['Ru'] = ((96,95.907598,5.54),(98,97.905287,1.87),(99,98.905939,12.76),(100,99.904220,12.60),(101,100.905582,17.06),(102,101.904350,31.55),(104,103.905430,18.62))
natural_abundances['Rh'] = ((103,102.905504,100),)
natural_abundances['Pd'] = ((102,101.905608,1.02),(104,103.904035,11.14),(105,104.905084,22.33),(106,105.903483,27.33),(108,107.903894,26.46),(110,109.905152,11.72))
natural_abundances['Ag'] = ((107,106.905093,51.839),(109,108.904756,48.161))
natural_abundances['Cd'] = ((106,105.906458,1.25),(108,107.904183,0.89),(110,109.903006,12.49),(111,110.904182,12.80),(112,111.902757,24.13),(113,112.904401,12.22),(114,113.903358,28.73),(116,115.904755,7.49))
natural_abundances['In'] = ((113,112.904061,4.29),(115,114.903878,95.71))
natural_abundances['Sn'] = ((112,111.904821,0.97),(114,113.902782,0.66),(115,114.903346,0.34),(116,115.901744,14.54),(117,116.902954,7.68),(118,117.901606,24.22),(119,118.903309,8.59),(120,119.902197,32.58),(122,121.903440,4.63),(124,123.905275,5.79))
natural_abundances['Sb'] = ((121,120.903818,57.21),(123,122.904216,42.79))
natural_abundances['Te'] = ((120,119.904020,0.09),(122,121.903047,2.55),(123,122.904273,0.89),(124,123.902819,4.74),(125,124.904425,7.07),(126,125.903306,18.84),(128,127.904461,31.74),(130,129.906223,34.08))
natural_abundances['I']  = ((127,126.904468,100),)
natural_abundances['Xe'] = ((124,123.905896,0.09),(126,125.904269,0.09),(128,127.903530,1.92),(129,128.904779,26.44),(130,129.903508,4.08),(131,130.905082,21.18),(132,131.904154,26.89),(134,133.905395,10.44),(136,135.907220,8.87))
natural_abundances['Cs'] = ((133,132.905447,100),)
natural_abundances['Ba'] = ((130,129.906310,0.106),(132,131.905056,0.101),(134,133.904503,2.417),(135,134.905683,6.592),(136,135.904570,7.854),(137,136.905821,11.232),(138,137.905241,71.698))
natural_abundances['La'] = ((138,137.907107,0.090),(139,138.906348,99.910))
natural_abundances['Ce'] = ((136,135.907144,0.185),(138,137.905986,0.251),(140,139.905434,88.450),(142,141.909240,11.114))
natural_abundances['Pr'] = ((141,140.907648,100),)
natural_abundances['Nd'] = ((142,141.907719,27.2),(143,142.909810,12.2),(144,143.910083,23.8),(145,144.912569,8.3),(146,145.913112,17.2),(148,147.916889,5.7),(150,149.920887,5.6))
natural_abundances['Pm'] = ((145,144.912744,100),)
natural_abundances['Sm'] = ((144,143.911995,3.07),(147,146.914893,14.99),(148,147.914818,11.24),(149,148.917180,13.82),(150,149.917271,7.38),(152,151.919728,26.75),(154,153.922205,22.75))
natural_abundances['Eu'] = ((151,150.919846,47.81),(153,152.921226,52.19))
natural_abundances['Gd'] = ((152,151.919788,0.20),(154,153.920862,2.18),(155,154.922619,14.80),(156,155.922120,20.47),(157,156.923957,15.65),(158,157.924101,24.84),(160,159.927051,21.86))
natural_abundances['Tb'] = ((159,158.925343,100),)
natural_abundances['Dy'] = ((156,155.924278,0.06),(158,157.924405,0.10),(160,159.925194,2.34),(161,160.926930,18.91),(162,161.926795,25.51),(163,162.928728,24.90),(164,163.929171,28.18))
natural_abundances['Ho'] = ((165,164.930319,100),)
natural_abundances['Er'] = ((162,161.928775,0.14),(164,163.929197,1.61),(166,165.930290,33.61),(167,166.932045,22.93),(168,167.932368,26.78),(170,169.935460,14.93))
natural_abundances['Tm'] = ((169,168.934211,100),)
natural_abundances['Yb'] = ((168,167.933894,0.13),(170,169.934759,3.04),(171,170.936322,14.28),(172,171.936378,21.83),(173,172.938207,16.13),(174,173.938858,31.83),(176,175.942568,12.76))
natural_abundances['Lu'] = ((175,174.940768,97.41),(176,175.942682,2.59))
natural_abundances['Hf'] = ((174,173.940040,0.16),(176,175.941402,5.26),(177,176.943220,18.60),(178,177.943698,27.28),(179,178.945815,13.62),(180,179.946549,35.08))
natural_abundances['Ta'] = ((180,179.947466,0.012),(181,180.947996,99.988))
natural_abundances['W']  = ((180,179.946706,0.12),(182,181.948206,26.50),(183,182.950224,14.31),(184,183.950933,30.64),(186,185.954362,28.43))
natural_abundances['Re'] = ((185,184.952956,37.40),(187,186.955751,62.60))
natural_abundances['Os'] = ((184,183.952491,0.02),(186,185.953838,1.59),(187,186.955748,1.96),(188,187.955836,13.24),(189,188.958145,16.15),(190,189.958445,26.26),(192,191.961479,40.78))
natural_abundances['Ir'] = ((191,190.960591,37.3),(193,192.962924,62.7))
natural_abundances['Pt'] = ((190,189.959930,0.014),(192,191.961035,0.782),(194,193.962664,32.967),(195,194.964774,33.832),(196,195.964935,25.242),(198,197.967876,7.163))
natural_abundances['Au'] = ((197,196.966552,100),)
natural_abundances['Hg'] = ((196,195.965815,0.15),(198,197.966752,9.97),(199,198.968262,16.87),(200,199.968309,23.10),(201,200.970285,13.18),(202,201.970626,29.86),(204,203.973476,6.87))
natural_abundances['Tl'] = ((203,202.972329,29.524),(205,204.974412,70.476))
natural_abundances['Pb'] = ((204,203.973029,1.4),(206,205.974449,24.1),(207,206.975881,22.1),(208,207.976636,52.4))
natural_abundances['Bi'] = ((209,208.980383,100),)


def parse_formula(elemental_formula):
    elemental_formula_parsed = []
    elemental_formula_regex = re.compile(r'(\(?)([a-zA-Z]+)(\d*)(\)?)(\d*)')
    elemental_formula_parts = elemental_formula_regex.split(elemental_formula)
    
    while '' in elemental_formula_parts: elemental_formula_parts.remove('')
    
    group_ends = [item for item in range(len(elemental_formula_parts)) if elemental_formula_parts[item] == ')']
    group_counts = []
    for ge in group_ends:
        try:
            gc = int(elemental_formula_parts[ge+1])
        except:
            gc = 1
        group_counts.append(gc)
    
    # parsing out concatenated elements
    ele = None
    freq = None
    group = False
    group_done = False
    while elemental_formula_parts:
        part = elemental_formula_parts.pop(0)
        if part == '(':
            group = True
            group_count = group_counts.pop(0)
            continue
        elif part == ')':
            group = False
            group_done = True
            continue
        try:
            freq = int(part)
            if group:
                freq *= group_count
            if not group_done:
                elemental_formula_parsed.append(freq)
        except:
            ele = part
            while ele:
                i = len(ele)
                while ele not in natural_abundances:
                    ele = ele[:-1]
                    i -= 1
                    if i < 1: return
                elemental_formula_parsed.append(ele)
                if group:
                    freq = group_count
                else:
                    freq = 1
                elemental_formula_parsed.append(freq)
                ele = part = part[i:]

        group_done = False
    
    # removing redundant frequencies
    elemental_formula_parsed2 = []
    previous_part_type = type(0)
    for part in elemental_formula_parsed:
        if type(part) == previous_part_type:
            elemental_formula_parsed2.pop()
        elemental_formula_parsed2.append(part)
        previous_part_type = type(part)
    
    elemental_formula_parsed3 = {}
    for i in range(0,len(elemental_formula_parsed2),2):
        ele = elemental_formula_parsed2[i]
        freq = elemental_formula_parsed2[i+1]
        if ele not in elemental_formula_parsed3:
            elemental_formula_parsed3[ele] = freq
        else:
            elemental_formula_parsed3[ele] += freq
    
    #print elemental_formula_parsed3
    return elemental_formula_parsed3


def Fernandez1996_correction(elemental_formula,input_mass_isotopomer_distribution):
    elements = parse_formula(elemental_formula)
    CM_temp = []
    counter = 0
    elements['C'] += 1
    while elements['C'] > 0:
        elements['C'] -= 1
        mass_isotopomer_distribution = [1]
        for e,f in elements.items():
            na_s = [(info[0],info[2]/100.) for info in natural_abundances[e]]
            element_mass_distribution = []
            prev_m = na_s[0][0] - 1
            while na_s:
                m,a = na_s[0]
                if m == prev_m + 1:
                    element_mass_distribution.append(a)
                    na_s.pop(0)
                else:
                    element_mass_distribution.append(0.)
                prev_m += 1
            for _ in range(f):
                mass_isotopomer_distribution = numpy.convolve(mass_isotopomer_distribution,element_mass_distribution)
        CM_temp.append(([0.]*counter) + list(mass_isotopomer_distribution))
        counter += 1
    
    # CM needs to be square
    num_rows = len(CM_temp)
    CM = [row[:num_rows] for row in CM_temp]
    CM = numpy.array(CM)
    
    # input_mass_isotopomer_distribution needs to be of the same dimension
    if len(input_mass_isotopomer_distribution) < num_rows:
        input_mass_isotopomer_distribution.extend([0.]*(num_rows-len(input_mass_isotopomer_distribution)))
    elif len(input_mass_isotopomer_distribution) > num_rows:
        input_mass_isotopomer_distribution = input_mass_isotopomer_distribution[:num_rows]
    input_mass_isotopomer_distribution = numpy.array(input_mass_isotopomer_distribution)
    
    return numpy.dot(input_mass_isotopomer_distribution,numpy.linalg.inv(CM))
