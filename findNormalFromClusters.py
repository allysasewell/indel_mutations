import os
import pandas as pd

bedFilePolynucHP = "rad16_combined_HP_54nuc.bed"
bedFileClusters = "UVB-Clusters/UVB-rad16_CombinedIndels_clusters.bed"
output = "clusters/rad16_NONclusters.bed"

workingDirectory = os.path.dirname(__file__)

polynucHPBedFile = {
        "Chromosome" : [],
        "Start(0)" : [],
        "Start(1)" : [],
        "Nucleotides" : [],
        "InDel" : [], 
        "Strand" : [],
        "Isolate" : [],
        "HP Count" : [],
        "Type" : []
}
clusterBedFile = {
        "Chromosome" : [],
        "Start(0)" : [],
        "Start(1)" : [],
        "Trinucleotides" : [],
        "InDel" : [], 
        "Strand" : [],
        "Isolate" : []
}

def parseLinePolynucHP(bedLine, dataframe): #parse HP bed file
       
        cleanSplitLine = bedLine.split()

        newRow = {"Chromosome" : cleanSplitLine[0], "Start(0)" : cleanSplitLine[1], "Start(1)" : 
        cleanSplitLine[2], "Nucleotides" : cleanSplitLine[3], 
        "InDel" : cleanSplitLine[4], "Strand" : cleanSplitLine[5], "Isolate" : cleanSplitLine[6], 
        "HP Count" : cleanSplitLine[7], "Type" : cleanSplitLine[8]}
     
        newRowDf = pd.DataFrame([newRow], columns = dataframe.columns) #make the new row in dataframe format
        dataframe = pd.concat([dataframe, newRowDf], ignore_index = True) #combine the two dataframes
     
        return dataframe

def fileIntoDfPolynuc(): #HP into bed file
    
    with open (os.path.join(workingDirectory, bedFilePolynucHP), "r") as file:
        PolynucHP_DF = pd.DataFrame(polynucHPBedFile) #make dictionary into dataframe

        for line in file: #iterate through each line

            if line.strip(): #if line not empty:
                                    
                PolynucHP_DF = parseLinePolynucHP(line, PolynucHP_DF)
                
    return PolynucHP_DF


def parseLineCluster(bedLine, dataframe): #parse cluster bed file
       
        cleanSplitLine = bedLine.split()

        newRow = {"Chromosome" : cleanSplitLine[0], "Start(0)" : cleanSplitLine[1], "Start(1)" : 
        cleanSplitLine[2], "Trinucleotides" : cleanSplitLine[3], 
        "InDel" : cleanSplitLine[4], "Strand" : cleanSplitLine[5], "Isolate" : cleanSplitLine[6]}
     
        newRowDf = pd.DataFrame([newRow], columns = dataframe.columns) #make the new row in dataframe format
        dataframe = pd.concat([dataframe, newRowDf], ignore_index = True) #combine the two dataframes
     
        return dataframe

def fileIntoDfCluster(): #HP into bed file
    
    with open (os.path.join(workingDirectory, bedFileClusters), "r") as file:
        clusterDF = pd.DataFrame(clusterBedFile)

        for line in file: #iterate through each line

            if line.strip(): #if line not empty:
                                    
                clusterDF = parseLineCluster(line, clusterDF)
                
    return clusterDF

#Parse files and place into data frames
polynucHP_DF = fileIntoDfPolynuc()
cluster_DF = fileIntoDfCluster()

clusterHP_DF = pd.merge(polynucHP_DF, cluster_DF, how = 'inner', on = ["Chromosome", "Start(0)", "Start(1)", "InDel", "Strand", "Isolate"]) #find the similar records
clusterHP_DF  = clusterHP_DF.drop(columns= ["Trinucleotides"])


clusterHP_DF.to_csv(output, sep = '\t', header = False, index = False)