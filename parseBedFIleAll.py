import os
import pandas as pd

workingDirectory = os.path.dirname(__file__)
allBedFiles = ["WT_DIPs_sorted_tandem_all.bed"]

dataBedFile = {
        "Chromosome" : [],
        "Start(0)" : [],
        "Start(1)" : [],
        "Trinucleotide" : [],
        "Deletion" : [],
        "Insertion" : [], 
        "Strand" : []

}

def parseLine(bedLine, dataframe):
       
        cleanSplitLine = bedLine.split()  #split by white space into strings

        newRow = {"Chromosome" : cleanSplitLine[0], "Start(0)" : cleanSplitLine[1], "Start(1)" : 
        cleanSplitLine[2], "Trinucleotide" : cleanSplitLine[3], 
         "Deletion" : cleanSplitLine[4], "Insertion" : cleanSplitLine[5], "Strand" : cleanSplitLine[6]}
     
        newRowDf = pd.DataFrame([newRow], columns = dataframe.columns) #make the new row in dataframe format
        dataframe = pd.concat([dataframe, newRowDf], ignore_index = True) #combine the two dataframes
     
        return dataframe
       
def fileIntoDf(x):
        with open(os.path.join(workingDirectory, allBedFiles[x]), "r") as file:
                myDfBed = pd.DataFrame(dataBedFile) #make dictionary into dataframe

                for line in file: #iterate through each line

                        if line.strip(): #if line not empty:
                                
                                myDfBed = parseLine(line, myDfBed)

        return myDfBed

wtAll = fileIntoDf(0).copy() 

print(wtAll.head())