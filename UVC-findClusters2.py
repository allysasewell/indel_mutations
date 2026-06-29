import os
import pandas as pd
workingDirectory = os.path.dirname(__file__)
strainName = "WT"

BedFile = "Combined_clusters/"+strainName+"_Combined_all_HP.bed"    
output = "UVC/Combined_clusters/"+strainName+"_CombinedIndel_clusters.bed"

coverage = 0

#ALL SNP file
UVYeast = {
        "Isolate" : [],        
        "Strain Name" : [],
        "Chromosome" : [],
        "Start(1)" : [],
        "End(1)" : [],
        "Type" : [],
        "HP Count" : []
}
#Bed File with HP
isolate = {
        "Chromosome" : [],        
        "Start(0)" : [],
        "Start(1)" : [],
        "Indel" : [],
        "Number" : [],
        "Isolate" : [],
        "Strand" : [],
        "HP Count": []
}

wtDF = pd.DataFrame(UVYeast)
rad16DF = pd.DataFrame(UVYeast)
rad26DF = pd.DataFrame(UVYeast)
rad30DF = pd.DataFrame(UVYeast)

def createRowUVYeast(line):
     newRow = { "Isolate" : line[0], "Strain Name" : line[4], "Chromosome" : "chr" +  line[6], 
                        "Start(1)" : line[8], "End(1)" : line[9], "Type" : line[10]}
     
     return newRow

with open(os.path.join(workingDirectory,"UV_yeast_Merged_SNP_DIPs_sorted_anz2_tandem.txt")) as f:
    data = f.read()
    lines = [s.strip().split() for s in data.splitlines()]

    for line in lines:
        if 'Mito' not in line[6] and 'ySR128' not in line[3]:
            
            if line[4] == "WT":
                if strainName == "WT":
                    newRow = createRowUVYeast(line)
            
                    newRowDF = pd.DataFrame([newRow], columns = wtDF.columns)
                    wtDF = pd.concat([wtDF, newRowDF], ignore_index = True)

            elif line[4] == "rad16D/rad16D":
                    if strainName == "rad16":
                        newRow = createRowUVYeast(line)

                        newRowDF = pd.DataFrame([newRow], columns = wtDF.columns)
                        rad16DF = pd.concat([rad16DF, newRowDF], ignore_index = True)
                    
            elif line[4] == "rad26D/rad26D" :
                    if strainName == "rad26":
                        newRow = createRowUVYeast(line)

                        newRowDF = pd.DataFrame([newRow], columns = wtDF.columns)
                        rad26DF = pd.concat([rad26DF, newRowDF], ignore_index = True)
            elif line[4] == "rad30D/rad30D":
                    if strainName == "rad30" :
                        newRow = createRowUVYeast(line)

                        newRowDF = pd.DataFrame([newRow], columns = wtDF.columns)
                        rad30DF = pd.concat([rad30DF, newRowDF], ignore_index = True)

def parseLine(bedLine, dataframe):
       
        line = bedLine.split()

        newRow =  { "Chromosome" : line[0],  
                "Start(0)" : line[1], "Start(1)" : line[2], "Indel" : line[4], "Number" : line[5],  "Isolate" : line[6], 
                "HP Count" : line[8], "Strand" : line[7]}
     
        newRowDf = pd.DataFrame([newRow], columns = dataframe.columns) #make the new row in dataframe format
        dataframe = pd.concat([dataframe, newRowDf], ignore_index = True) #combine the two dataframes
     
        return dataframe
       
def fileIntoDf(myDfBed):
        with open(os.path.join(workingDirectory, BedFile), "r") as file:
                myDfBed = pd.DataFrame(myDfBed) #make dictionary into dataframe

                for line in file: #iterate through each line
                        if line.strip(): #if line not empty:
                                myDfBed = parseLine(line, myDfBed)

        return myDfBed




def checkHP(strainDF, isolateDF):

    hpRow = {"HP Count" : []}

    for y in range(len(strainDF)):
        hpRow[y] = 0
        if strainDF.at[y, "Type"] != "SNP":
            for x in range(len(isolateDF)):
                if int(strainDF.at[y, "Start(1)"])== int(isolateDF.at[x, "Start(1)"]): #if the same
                    if strainDF.at[y, "Chromosome"] == isolateDF.at[x, "Chromosome"] and strainDF.at[y, "Isolate"] == isolateDF.at[x, "Isolate"]:
                        hpRow[y] = int(isolateDF.at[x, "HP Count"]) - 1
                        break

    strainDF["HP Count"] = hpRow
    
    return strainDF

def findClusters(strainDF, isolateDF):
    clusterDF = pd.DataFrame(isolate)

    for y in range(len(strainDF)):
        for x in range(len(isolateDF)):
            if strainDF.at[y, "Chromosome"] == isolateDF.at[x, "Chromosome"]:

               if strainDF.at[y, "Isolate"] == isolateDF.at[x, "Isolate"]:
                    
                    isolateStart = int(isolateDF.at[x, "Start(1)"])
                    hpLengthIsolate = int(isolateDF.at[x, "HP Count"])
                    strainStart = int(strainDF.at[y, "Start(1)"])
                    hpLengthStrain = int(strainDF.at[y, "HP Count"])

                    if  isolateStart > strainStart : #if to the left
                        leftRange = isolateStart - (strainStart + hpLengthStrain)
                        rightRange = 11
                    elif isolateStart + hpLengthIsolate < strainStart: #if to the right of the whole HP
                        rightRange = abs((isolateStart + hpLengthIsolate) - strainStart)
                        leftRange = 11
                    elif isolateStart == strainStart: 
                        rightRange = 11      #is the HP itself
                        leftRange = 11
                    else: #it is inside the homopolymer
                         leftRange = 1
                         rightRange = 1

                    if leftRange <= 10 or rightRange <= 10: #if within desired range, it is complex
                        row = isolateDF.loc[x].to_frame().T

                        isolateDF.drop([x], axis = 0, inplace = True)
                        isolateDF.reset_index(inplace = True, drop = True)

                        clusterDF = pd.concat([clusterDF, row])

                        break

    return clusterDF



isolateDF = fileIntoDf(isolate).copy()



if(strainName == "WT"):   
    wtDF = checkHP(wtDF, isolateDF)           
    clusterDF = findClusters(wtDF, isolateDF)
elif(strainName == "rad16"):  
    rad16DF= checkHP(rad16DF, isolateDF)            
    clusterDF = findClusters(rad16DF, isolateDF)
elif(strainName == "rad26"):   
    rad26DF= checkHP(rad26DF, isolateDF)               
    clusterDF = findClusters(rad26DF, isolateDF)
elif(strainName == "rad30"):
    rad30DF= checkHP(rad30DF, isolateDF)               
    clusterDF = findClusters(rad30DF, isolateDF)

clusterDF.to_csv(output, sep = '\t', header = False, index = False)