import schema
import numpy as np
import matplotlib.pyplot as plt


def makeBarGraph(data, variation, name):
    '''
    Draws a graph containing the ordered data and OD600 from each well
    Saves based on sample number
    '''
    plt.clf()
    plt.figure(figsize=(15,5))

    #Add motif lines - the location of the centre of each motif
    plt.axvline(x=11, color='grey', alpha = 0.3, linestyle='dashed') #C1
    plt.axvline(x=57, color='grey', alpha = 0.3, linestyle='dashed') #C2
    plt.axvline(x=142, color='grey', alpha = 0.3, linestyle='dashed') #C3
    plt.axvline(x=179, color='grey', alpha = 0.3, linestyle='dashed') #C4
    plt.axvline(x=300, color='grey', alpha = 0.3, linestyle='dashed') #C5
    plt.axvline(x=334, color='grey', alpha = 0.3, linestyle='dashed') #C6
    plt.axvline(x=351, color='grey', alpha = 0.3, linestyle='dashed') #C7
    plt.axvline(x=501, color='grey', alpha = 0.3, linestyle='dashed') #A1
    plt.axvline(x=548, color='grey', alpha = 0.3, linestyle='dashed') #A2
    plt.axvline(x=627, color='grey', alpha = 0.3, linestyle='dashed') #A3
    plt.axvline(x=675, color='grey', alpha = 0.3, linestyle='dashed') #A4
    plt.axvline(x=778, color='grey', alpha = 0.3, linestyle='dashed') #A5
    plt.axvline(x=872, color='grey', alpha = 0.3, linestyle='dashed') #A6
    plt.axvline(x=906, color='grey', alpha = 0.3, linestyle='dashed') #A7
    plt.axvline(x=931, color='grey', alpha = 0.3, linestyle='dashed') #A8
    plt.axvline(x=1005, color='grey', alpha = 0.3, linestyle='dashed') #A9
    plt.axvline(x=1025, color='grey', alpha = 0.3, linestyle='dashed') #A10

    #Add data and error bars
    wells = range(0, len(data))
    plt.bar(wells, data, width = 1, alpha=0.5, align='center')
    plt.errorbar(wells, data, variation, linestyle='None')
    axes = plt.gca()
    axes.set_xlim([0, len(data)])

    #Add vertical lines
    plt.axvline(x=430, color='k', alpha = 0.8) #New A fwd - halfway between the C7 and A1 motifs
    plt.axvline(x=1034, color='k', alpha = 0.8) #A Rev
    
    #Add labels
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)
    axes.set_xlabel('Alignment position')
    axes.set_ylabel('Average number of clashes')
    
    plt.savefig('Profile_' + str(name))


def wrapperProfileGraph(parentFile, contactFile):
    '''
    Draws a graph of the number of clashes at each recombination point
    '''
    pdbName = contactFile.split('_')[0][-4:]
    parent_list = schema.readMultipleSequenceAlignmentFile(file(parentFile, 'r'))
    parents = [p for (k,p) in parent_list]

    pdb_contacts = schema.readContactFile(file(contactFile, 'r'))

    clash_data = [[] for x in parents[0]]
    for i in range(1, len(parents)):
        print i
        #This reshuffles the alignment to make the first and second sequences the ones analysed. It was needed as SCHEMA is limited to 9 sequences.
        newList = [parents[0], parents[i]]
        for x in range(1, len(parents)):
            if not i==x:
                newList.append(parents[x])    
        #Graphs for hotspots
        for residue in range(0, len(parents[0])):
            crossovers = [residue]  
            contacts = schema.getSCHEMAContactsWithCrossovers(pdb_contacts, newList, crossovers)
            fragments = schema.getFragments(crossovers, parents[0])
            
            clash_data[residue].append(schema.getChimeraDisruption('21', contacts, fragments, newList))
    means = [np.mean(values) for values in clash_data]
    StDev = [np.std(values) for values in clash_data]
    makeBarGraph(means, StDev, pdbName)

wrapperProfileGraph('Alignments\Pa11_aligned_SCHEMA.aln', 'Contacts\\2VSQ_contacts.txt')


