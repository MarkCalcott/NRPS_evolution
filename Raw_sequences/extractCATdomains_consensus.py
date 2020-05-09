'''
This program identifies all gbk files in a folder. It then
searches within each protein for a LCL domain immediately upstream
to an A domain. These are then extracted and added to a fasta file.

By looking at CDS features, the C domains are always upstream to an
A domain.
'''


from Bio import SeqIO
import os
import re

def getFileNames(path):
    '''
    gets the file names for all gbk files in the specific folder
    '''
    files = []
    for file in os.listdir(os.getcwd()+path):
        if file.endswith(".gbk"):
            files.append(os.getcwd()+path+file)
    return files

def findDomains(filename):
    '''
    Takes a genbank file and finds the regions annotated as domains
    '''
    rec = SeqIO.read(filename, 'genbank')
    results = []
    for feature in rec.features:
        
        if feature.type == "CDS":
            #Each CDS refers to a protein
            if 'sec_met' in feature.qualifiers:
                
                type_sec_met_qualifiers = [feat for feat in feature.qualifiers['sec_met'] if "NRPS/PKS Domain: " in feat]
                C_region = None #reset to None
                A_region = None
                T_region = None
                
                for qualifier in type_sec_met_qualifiers:
                    Cdomain = 'Condensation_LCL'
                    Adomain = 'AMP-binding'
                    Tdomain = 'PCP'
                    if Cdomain in qualifier:
                        C_region = [int(x) for x in re.findall(Cdomain +' \(([0-9]+)-([0-9]+)\)', qualifier)[0]] #Finds a tuple for the region
                    elif Adomain in qualifier:
                        A_region = [int(x) for x in re.findall(Adomain +' \(([0-9]+)-([0-9]+)\)', qualifier)[0]] #Finds a tuple for the region
                        specificity = re.findall(r'Minowa\), (.+) \(consensus\)', qualifier)[0] #Extract consensus specificity
                    elif Tdomain in qualifier:
                        T_region = [int(x) for x in re.findall(Tdomain +' \(([0-9]+)-([0-9]+)\)', qualifier)[0]]
                        #Extract the C and A domain sequences
                        if C_region and A_region and T_region:
                            if not specificity == "nrp":
                                if A_region[0]-C_region[1] < 250 and T_region[0]-A_region[1] < 250: # C, A and T domain separated by less than 250, as it is a protein, the C is always first
                                    
                                    sequence = feature.extract(rec.seq)

                                    extract_region = [C_region[0], T_region[1]]

                                    sequenceCA = sequence[3*extract_region[0]:3*extract_region[1]]
                                    if len(sequenceCA) == 0:
                                        print ('Error', extract_region[0], extract_region[1])
                                    results.append((specificity, sequenceCA))

    return results

def wrapper(species, path):
    savefilename = species+'CATdomains_consensus.fasta'
    saveFileCA = open(savefilename, 'w')

    for fileName in getFileNames(path):
        domains = findDomains(fileName)
        for index in range(len(domains)):
            clustername = fileName.split('\\')[-1][:-4]
            saveFileCA.write('> %s_CA%d_%s\n' % (clustername, index+1, domains[index][0]))
            saveFileCA.write(str(domains[index][1]))
            saveFileCA.write('\n')
        
    saveFileCA.close()


wrapper('Bacillus_2', '\\Bacillus\\')
wrapper('Pseudomonas_2', '\\Pseudomonas\\')
wrapper('Streptomyces_2', '\\Streptomyces\\')
