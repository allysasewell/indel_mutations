import os
import pandas as pd
workingDirectory = os.path.dirname(__file__)

strain = "WT"

BedFile = "UVB-"+ strain+"_Combined_all_HP.bed" 
BedFile2= "UVB-Clusters/UVB-"+ strain+"_CombinedIndel_clusters.bed"
output = "UVB/UVB-Clusters/non_cluster/UVB-"+ strain+"_NONclusters.bed"

isolate = {
        "Chromosome" : [],        
        "Start(0)" : [],
        "Start(1)" : [],
        "Nucleotides" : [],
        "InDel" : [],
        "Isolate" : [],
        "Indel Length": [],
        "Strand" : [],
        "HP Count" : [],
        "Type" : [],
        "HP Start" : [],
        "HP End" : []
}

def parseLine(bedLine, dataframe):
       
        line = bedLine.split()

        newRow = { "Chromosome" : line[0], "Start(0)" : line[1], "Start(1)" : line[2], "Nucleotides" : line[3], 
                  "InDel" : line[4],  "Isolate" : line[6], "Indel Length": line[5], "Strand" : line[7], 
                  "HP Count" : line[8], "Type" : line[9], "HP Start" : line[10], "HP End" : line[11]}
     
        newRowDf = pd.DataFrame([newRow], columns = dataframe.columns) #make the new row in dataframe format
        dataframe = pd.concat([dataframe, newRowDf], ignore_index = True) #combine the two dataframes
     
        return dataframe
 


def parseLineCLuster(bedLine, dataframe):
       
        line = bedLine.split()

        newRow = { "Chromosome" : line[0], "Start(0)" : line[1], "Start(1)" : line[2], "Nucleotides" : line[3], 
                  "InDel" : line[4],  "Isolate" : line[5], "Indel Length": line[6], "Strand" : line[7], 
                  "HP Count" : line[8], "Type" : line[9], "HP Start" : line[10], "HP End" : line[11]}
     
        newRowDf = pd.DataFrame([newRow], columns = dataframe.columns) #make the new row in dataframe format
        dataframe = pd.concat([dataframe, newRowDf], ignore_index = True) #combine the two dataframes
     
        return dataframe
 
def fileIntoDfCluster(myDfBed, BedFile):
        with open(os.path.join(workingDirectory, BedFile), "r") as file:
                myDfBed = pd.DataFrame(myDfBed) #make dictionary into dataframe

                for line in file: #iterate through each line
                        if line.strip(): #if line not empty:
                                myDfBed = parseLineCLuster(line, myDfBed)

        return myDfBed


def fileIntoDf(myDfBed, BedFile):
        with open(os.path.join(workingDirectory, BedFile), "r") as file:
                myDfBed = pd.DataFrame(myDfBed) #make dictionary into dataframe

                for line in file: #iterate through each line
                        if line.strip(): #if line not empty:
                                myDfBed = parseLine(line, myDfBed)

        return myDfBed

allDf = fileIntoDf(isolate, BedFile).copy()
clusterDf = fileIntoDfCluster(isolate, BedFile2).copy()

combDf = pd.concat([allDf, clusterDf])
combDf = combDf.drop_duplicates(keep = False)


combDf.to_csv(output, sep = '\t', header = False, index = False)