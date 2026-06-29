import os
import pandas as pd
workingDirectory = os.path.dirname(__file__)
BedFile = "WT_sorted_combined_Trinuc.bed"    
strainName = "WT"
output = "clusters/WT_clusters.bed"

with open(os.path.join(workingDirectory,"UV_yeast_Merged_SNP_DIPs_sorted_anz2_tandem.txt")) as f:
    data = f.read()
lines = [s.strip().split() for s in data.splitlines()]
genotypes = {'WT': 0, 'rad16':0, 'rad26':0, 'rad30': 0}
revComplementDict  = {
        'A' : 'T',
        'G' : 'C',
        'T' : 'A',
        'C' : 'G'
}
def GetReverseComplement(dna):
        dna2 = ""
        for base in dna: 
                dna2 += revComplementDict[base]
        return dna2[::-1]

WT_chromosomes = {}
rad16_chromosomes = {}
rad26_chromosomes = {}
rad30_chromosomes = {}
WT_indels = 0
rad16_indels = 0
rad26_indels = 0
rad30_indels = 0
genotype = []
background = []
chromosome = []
position1 = []
position2 = []
mutation_type = []
length = []
allele = []
mutation = []
frequency = []
left_reference = []
right_reference = []
left_consensus = []
right_consensus = []
trinucleotide = []
strand = []
isolate_list = {}
lines.pop(0)
for line in lines:
   
    if 'Mito' not in line[6] and 'ySR128' not in line[3]:
    
        genotype.append(line[4])
    
        mutation_type.append(line[10])
        if 'WT' in line[4]:
            genotypes['WT'] = genotypes['WT'] + 1
        if 'rad16' in line[4]:
            genotypes['rad16'] = genotypes['rad16'] + 1
        if 'rad26' in line[4]:
            genotypes['rad26'] = genotypes['rad26'] + 1
        if 'rad30' in line[4]:
            genotypes['rad30'] = genotypes['rad30'] + 1
        chromosome.append('chr' + str(line[6]))
        if 'WT' in line[4]:
            if str(line[6]) in WT_chromosomes.keys():
                WT_chromosomes[str(line[6])] = WT_chromosomes[str(line[6])] + 1

            else:
                WT_chromosomes[str(line[6])] = 1
            if line[10] == "DIP" or line[10] == "Complex DIP":
                WT_indels = WT_indels + 1
        if 'rad16' in line[4]:
            if str(line[6]) in rad16_chromosomes.keys():
                rad16_chromosomes[str(line[6])] = rad16_chromosomes[str(line[6])] + 1
            else:
                rad16_chromosomes[str(line[6])] = 1
            if line[10] == "DIP" or line[10] == "Complex DIP":
                rad16_indels = rad16_indels + 1
        if 'rad26' in line[4]:
            if str(line[6]) in rad26_chromosomes.keys():
                rad26_chromosomes[str(line[6])] = rad26_chromosomes[str(line[6])] + 1
            else:
                rad26_chromosomes[str(line[6])] = 1
            if line[10] == "DIP" or line[10] == "Complex DIP":
                rad26_indels = rad26_indels + 1
        if 'rad30' in line[4]:
            if str(line[6]) in rad30_chromosomes.keys():
                rad30_chromosomes[str(line[6])] = rad30_chromosomes[str(line[6])] + 1
            else:
                rad30_chromosomes[str(line[6])] = 1
            if line[10] == "DIP" or line[10] == "Complex DIP":
                rad30_indels = rad30_indels + 1
        position1.append(line[8])
        position2.append(line[9])
        length.append(line[11])
        allele.append(line[12])
        if len(line) >= 22: 
            if '-' not in  line[18] and line[18].isalpha() == False:
                if float(line[18]) >= 90:
                    frequency.append('Homozygous')
                else:
                    frequency.append('Heterozygous')
            elif float(line[19]) >= 90:
                frequency.append('Homozygous')
            else:
                frequency.append('Heterozygous')
            left_reference.append(line[23])
            right_reference.append(line[24])
            left_consensus.append(line[25])
            right_consensus.append(line[26])
            if line[27].isalpha() == False:
                if line[24][1] == "G" or line[24][1] == "A":
                    trinucleotide.append(GetReverseComplement(line[24]))
                    strand.append('-')
                    if '-' not in line[13]:
                        mutation.append(GetReverseComplement(line[13]))
                    else:
                        mutation.append(line[13])
                else:
                    trinucleotide.append(line[24])
                    mutation.append(line[13])
                    strand.append('+')
            else:
                if line[27][1] == "A" or line[27][1] == "G":
                    trinucleotide.append(GetReverseComplement(line[27]))
                    strand.append('-')
                    if '-' not in line[13]:
                        mutation.append(GetReverseComplement(line[13]))
                    else:
                        mutation.append(line[13])
                else:
                    trinucleotide.append(line[27])
                    mutation.append(line[13])
                    strand.append('+')

    
        else:
            frequency.append('NA')
            left_reference.append('NA')
            right_reference.append('NA')
            left_consensus.append('NA')
            right_consensus.append('NA')
            trinucleotide.append('NA')
        if (line[0], line[4]) not in isolate_list.keys():
            isolate_list[(line[0],  line[4])] = [(line[9], 'chr' + str(line[6]), trinucleotide[len(trinucleotide) - 1], mutation[len(mutation) - 1], strand[len(strand) - 1], mutation_type[len(mutation_type) - 1], line[0], line[10])]
        else:
            isolate_list[(line[0],  line[4])].append((line[9], 'chr'+ str(line[6]),  trinucleotide[len(trinucleotide) - 1], mutation[len(mutation) - 1], strand[len(strand) - 1], mutation_type[len(mutation_type) - 1], line[0], line[10]))

print(genotypes)
print(WT_indels/genotypes['WT'])
print(rad16_indels/genotypes['rad16'])
print(rad26_indels/genotypes['rad26'])
print(rad30_indels/genotypes['rad30'])




#cluster_file = open('clusters.txt', 'w+')

def FindClusters( isolate_list, strain):
    lengths = {}
    clusters = {}
    cluster_numbers = {}
    for key in isolate_list.keys():
        number = 0
        if strain in key[1]:
            for j in range (0,len (isolate_list[key]) - 1):
                flag = True
                n = 1
                while flag == True:
                    if j + n <= len(isolate_list[key]) - 1:
                        if int(isolate_list[key][j + n][0]) <= int((isolate_list[key][j + n - 1][0])) + 10:
                            if isolate_list[key][j + n][1] == isolate_list[key][j + n - 1][1]:
                                if isolate_list[key][j] not in clusters.keys() and isolate_list[key][j][7] == "DIP":
                                    if j == 0 or isolate_list[key][j - n] not in clusters.keys():
                                        clusters[isolate_list[key][j]] = [isolate_list[key][j]]
                                        number = number + 1
            
                                        if isolate_list[key][j + n][7] == "DIP":
                                            clusters[isolate_list[key][j]].append(isolate_list[key][j + n])
                                    
                                      
                                    else:
                                        n = n + 1
                                elif isolate_list[key][j + n] == "DIP":
                                    clusters[isolate_list[key][j]].append(isolate_list[key][j + n])
                                    number = number + 1
                                n = n + 1
                            else:
                                flag = False
                        else:
                            flag = False
                    else:
                        flag = False
        cluster_numbers[key] = number
    for key in clusters.keys():
        if len(clusters[key]) in lengths.keys():
            lengths[len(clusters[key])] = lengths[len(clusters[key])] + 1
        else:
            lengths[len(clusters[key])] = 1
        if len(clusters[key]) > 5:
            print(str(clusters[key][0]) + ": " + str(len(clusters[key])))
    for key in lengths.keys():
        print(str(key) + ": " + str(lengths[key]))
    for key in cluster_numbers.keys():
        print(str(key) + ": " + str(cluster_numbers[key]/ len(isolate_list[key])))
    new_clusters = {}
    for key in clusters.keys():
        for i in range(0, len(clusters[key])):
            if (clusters[key][i][6], clusters[key][i][1], clusters[key][i][0]) not in new_clusters.keys():
                new_clusters[(clusters[key][i][6], clusters[key][i][1], clusters[key][i][0])] = [(clusters[key][i][6], clusters[key][i][1], clusters[key][i][0])]



    return new_clusters


cluster = FindClusters(isolate_list, strainName)
#print(cluster, file = cluster_file)




dataBedFile = {
        "Chromosome" : [],
        "Start(0)" : [],
        "Start(1)" : [],
        "Trinucleotide" : [],
        "Insertion" : [],
        "Strand" : [],
        "Isolate" : []

}



myDfBed = pd.DataFrame(dataBedFile) #make dictionary into dataframe
clusterDF = pd.DataFrame(dataBedFile)

def parseLine(bedLine, dataframe):
       
        cleanSplitLine = bedLine.split()

        newRow = {"Chromosome" : cleanSplitLine[0], "Start(0)" : cleanSplitLine[1], "Start(1)" : 
        cleanSplitLine[2], "Trinucleotide" : cleanSplitLine[3], 
        "Insertion" : cleanSplitLine[4], "Strand" : cleanSplitLine[5], "Isolate" : cleanSplitLine[6]}
     
        newRowDf = pd.DataFrame([newRow], columns = dataframe.columns) #make the new row in dataframe format
        dataframe = pd.concat([dataframe, newRowDf], ignore_index = True) #combine the two dataframes
     
        return dataframe
       
def fileIntoDf(myDfBed):
        with open(os.path.join(workingDirectory, BedFile), "r") as file:
                myDfBed = pd.DataFrame(dataBedFile) #make dictionary into dataframe

                for line in file: #iterate through each line

                        if line.strip(): #if line not empty:
                                
                                myDfBed = parseLine(line, myDfBed)

        return myDfBed

trinucDF = fileIntoDf(myDfBed).copy()





for key in cluster.keys():
    for x in range(len(trinucDF)):
        if key[1] == trinucDF.at[x, "Chromosome"]:
            if key[2] == trinucDF.at[x, "Start(1)"]:
                if key[0] == trinucDF.at[x, "Isolate"]:
                    row = trinucDF.loc[x].to_frame().T
                    clusterDF = pd.concat([clusterDF, row])
            

clusterDF.to_csv(output, sep = '\t', header = False, index = False)