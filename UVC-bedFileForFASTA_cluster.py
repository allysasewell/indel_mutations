import os
import pandas as pd
workingDirectory = os.path.dirname(__file__)

strain= "WT"

BedFile = "Combined_Clusters/" + strain +"_CombinedAll_NONcluster.bed"    

output1 = "UVC/Combined_Clusters/nonCluster_fasta/del fasta/UVC-" + strain +"_del_NONclusters.bed"
output2 = "UVC/Combined_Clusters/nonCluster_fasta/ins fasta/UVC-" + strain +"_ins_NONclusters.bed"


#BedFile = "Combined_Clusters/" + strain +"_CombinedIndel_clusters.bed"    

#output1 = "UVC/Combined_Clusters/cluster_fasta/del fasta/UVC-" + strain +"_del_clusters.bed"
#output2 = "UVC/Combined_Clusters/cluster_fasta/ins fasta/UVC-" + strain +"_ins_clusters.bed"


#new FASTA file structure
isolate = {
        "Chromosome" : [],        
        "Start(0)" : [],
        "Start(1)" : [],
        "Isolate" : [],
        "HP Count" : [],
        "Strand" : []
}

def checkLen(chr): #chromosome lengths
        if chr == 'chrI':
                return 230218
        if chr == 'chrII':
                return 813184
        if chr == 'chrIII':
                return 316620
        if chr == 'chrIV':
                return 1531933
        if chr == 'chrV':
                return 576874
        if chr == 'chrVI':
                return 270161
        if chr == 'chrVII':
                return 1090940
        if chr == 'chrVIII':
                return 562643
        if chr == 'chrIX':
                return 439888
        if chr == 'chrX':
                return 745751
        if chr == 'chrXI':
                return 666816
        if chr == 'chrXII':
                return 1078177
        if chr == 'chrXIII':
                return 924431
        if chr == 'chrXIV':
                return 784333
        if chr == 'chrXV':
                return 1091291
        if chr == 'chrXVI':
                return 948066

def checkStartandEnd(line):  #check that is not near ends of chromosome
        if int(line[1]) - 27 < 500:
                return False
        if int(line[2]) > (checkLen(line[0]) - 500):
                return False
        else:
                return True


def parseLine(bedLine, dataframe):

        newRow = { "Chromosome" : bedLine[0], "Start(0)" : int(bedLine[1]) - 26, "Start(1)" : int(bedLine[2])  + 26, "Isolate" : bedLine[5],
                 "HP Count" : bedLine[7], "Strand" : bedLine[6], "Indel" : bedLine[4]}
     
        newRowDf = pd.DataFrame([newRow], columns = dataframe.columns) #make the new row in dataframe format
        dataframe = pd.concat([dataframe, newRowDf], ignore_index = True) #combine the two dataframes
     
        return dataframe


with open(os.path.join(workingDirectory, BedFile), "r") as file:
                dfIns = pd.DataFrame(isolate) #make dictionary into dataframe
                dfDel = pd.DataFrame(isolate)

                for line in file: #iterate through each line
                        if line.strip(): #if line not empty:
                                line = line.split()
                                if int(line[4]) == 1: #if single indel
                                        if '-' not in line[3]:
                                               if checkStartandEnd(line):
                                                       dfIns = parseLine(line, dfIns) 
                                        else:
                                                if checkStartandEnd(line):
                                                        dfDel = parseLine(line, dfDel)




dfDel.to_csv(output1, sep = '\t', header = False, index = False, float_format= "%.0f")
dfIns.to_csv(output2, sep = '\t', header = False, index = False, float_format= "%.0f")



#bedtools getfasta -fi genome.fasta -bed UVB-rad30.bed -nameOnly -s -fo UVB-rad30.fasta