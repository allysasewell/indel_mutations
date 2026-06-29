import os
import pandas as pd

workingDirectory = os.path.dirname(__file__)
allBedFiles = ["rad16_sorted_insertions_Trinuc.bed", "rad26_sorted_insertions_Trinuc.bed", 
               "rad30_sorted_insertions_Trinuc.bed", "WT_sorted_insertions_Trinuc.bed"]

dataBedFile = {
        "Chromosome" : [],
        "Start(0)" : [],
        "Start(1)" : [],
        "Trinucleotide" : [],
        "Insertion" : [],
        "Strand" : []

}

myDfBed = pd.DataFrame(dataBedFile) #make dictionary into dataframe

def parseLine(bedLine, dataframe):
       
        cleanSplitLine = bedLine.split()

        newRow = {"Chromosome" : cleanSplitLine[0], "Start(0)" : cleanSplitLine[1], "Start(1)" : 
        cleanSplitLine[2], "Trinucleotide" : cleanSplitLine[3], 
        "Insertion" : cleanSplitLine[4], "Strand" : cleanSplitLine[5]}
     
        newRowDf = pd.DataFrame([newRow], columns = dataframe.columns) #make the new row in dataframe format
        dataframe = pd.concat([dataframe, newRowDf], ignore_index = True) #combine the two dataframes
     
        return dataframe
       
def fileIntoDf(x, myDfBed):
        with open(os.path.join(workingDirectory, allBedFiles[x]), "r") as file:
                myDfBed = pd.DataFrame(dataBedFile) #make dictionary into dataframe

                for line in file: #iterate through each line

                        if line.strip(): #if line not empty:
                                
                                myDfBed = parseLine(line, myDfBed)

        return myDfBed

rad16Del = fileIntoDf(0, myDfBed).copy() #copy dataframe into variable

rad26Del = fileIntoDf(1, myDfBed).copy()

rad30Del = fileIntoDf(2, myDfBed).copy()                        #deep = True?

wtDel = fileIntoDf(3, myDfBed).copy()


print(rad16Del.head())