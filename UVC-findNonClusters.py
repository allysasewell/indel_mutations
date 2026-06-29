import os
import pandas as pd
workingDirectory = os.path.dirname(__file__)

strain = "rad30"

BedFile = strain +"_combined_all_HP.bed"  
BedFile2=  strain +"_CombinedIndel_clusters.bed"  
output = "UVC/Combined_Clusters/" + strain +"_CombinedAll_NONcluster.bed"



isolate = {
        "Chromosome" : [],        
        "Start(0)" : [],
        "Start(1)" : [],
        "InDel" : [],
        "Number" : [],
        "Isolate" : [],
        "Strain" : [],
        "HP Count" : []
}

def parseLineCluster(bedLine, dataframe):
       
        line = bedLine.split()

        newRow = {"Chromosome" : line[0], "Start(0)" : line[1], "Start(1)" : line[2],
                  "InDel" : line[3], "Number" : line[4],  "Isolate" : line[5], "Strain" :line[6],
                  "HP Count" : line[7]}
     
        newRowDf = pd.DataFrame([newRow], columns = dataframe.columns) #make the new row in dataframe format
        dataframe = pd.concat([dataframe, newRowDf], ignore_index = True) #combine the two dataframes
     
        return dataframe
       
def fileIntoDf(myDfBed, BedFile):
        with open(os.path.join(workingDirectory, BedFile), "r") as file:
                myDfBed = pd.DataFrame(myDfBed) #make dictionary into dataframe

                for line in file: #iterate through each line
                        if line.strip(): #if line not empty:
                                myDfBed = parseLineCluster(line, myDfBed)

        return myDfBed

def parseLineHP(bedLine, dataframe):
       
        line = bedLine.split()

        newRow = {"Chromosome" : line[0], "Start(0)" : line[1], "Start(1)" : line[2],
                  "InDel" : line[4], "Number" : line[5],  "Isolate" : line[6], "Strain" : line[7],
                  "HP Count" : line[8]}
     
        newRowDf = pd.DataFrame([newRow], columns = dataframe.columns) #make the new row in dataframe format
        dataframe = pd.concat([dataframe, newRowDf], ignore_index = True) #combine the two dataframes
     
        return dataframe
       
def fileIntoDf2(myDfBed, BedFile):
        with open(os.path.join(workingDirectory, BedFile), "r") as file:
                myDfBed = pd.DataFrame(myDfBed) #make dictionary into dataframe

                for line in file: #iterate through each line
                        if line.strip(): #if line not empty:
                                myDfBed = parseLineHP(line, myDfBed)

        return myDfBed



allDf = fileIntoDf2(isolate, BedFile).copy()
clusterDf = fileIntoDf(isolate, BedFile2).copy()

combDf = pd.concat([allDf, clusterDf])
combDf = combDf.drop_duplicates(keep = False)


combDf.to_csv(output, sep = '\t', header = False, index = False)                  


