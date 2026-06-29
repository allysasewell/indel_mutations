import os
import pandas as pd
workingDirectory = os.path.dirname(__file__)

name = "rad16"
in_or_del = "ins"
BedFile = "UVB-BedForFASTA/UVB-" + name +"_"+ in_or_del+".bed"
output = "UVB-BedForFASTA/UVB-" + name + "_" + in_or_del+"_+-_norm.bed"

fasta = {
    "Chromosome" : [],
    "Start(0)" : [],
    "Start(1)" : [],
    "Isolate" : [],
    "HP Count" : [],
    "Strand" : []
}


def parseLine(line, dataframe):

        newRow = { "Chromosome" : line[0], "Start(0)" : line[1], "Start(1)" : line[2], "Isolate" : line[3],
                 "HP Count" : line[4], "Strand" : line[5]}
     
        newRowDf = pd.DataFrame([newRow], columns = dataframe.columns) #make the new row in dataframe format
        dataframe = pd.concat([dataframe, newRowDf], ignore_index = True) #combine the two dataframes
     
        return dataframe

def parseLine_minus(line, dataframe): #-1 for insertion?

        newRow = { "Chromosome" : line[0], "Start(0)" : int(line[1]), "Start(1)" : int(line[2]), "Isolate" : line[3],
                 "HP Count" : line[4], "Strand" : line[5]}
     
        newRowDf = pd.DataFrame([newRow], columns = dataframe.columns) #make the new row in dataframe format
        dataframe = pd.concat([dataframe, newRowDf], ignore_index = True) #combine the two dataframes
     
        return dataframe

with open(os.path.join(workingDirectory, BedFile), "r") as file:
                df = pd.DataFrame(fasta) #make dictionary into dataframe

                for line in file: #iterate through each line
                        if line.strip(): #if line not empty:
                                line = line.split()
                                if line[5] == '+':
                                    df = parseLine(line, df)

                                else:
                                    df = parseLine_minus(line, df)

df["Start(0)"] = df["Start(0)"].fillna(0).astype(int) #addition buggy with int, fixes issue
df["Start(1)"] = df["Start(1)"].fillna(0).astype(int)


df.to_csv(output, sep = '\t', header = False, index = False, float_format= "%.0f")