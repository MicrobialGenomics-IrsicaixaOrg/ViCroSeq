
#Load ShortRead library for working with FastQ
library('ShortRead');

#Load functions 
source("simContFunc.R")

#Arguments
#Directory FastQ files
dirFastq<- "RawData/RawData"

#qaSummary FastQ files
fls <- dir(dirFastq, "*fastq", full=TRUE)
qaSummary <- qa(fls, type="fastq")
qaSummary[['readCounts']]

#Percentage of contamination
contPerct<-as.integer(c(5,10,15,5))

#Get samples
samples<-getSamples(qaSummary)
print(samples)

#get sampleMatrix
sampleMatrix<- getMatrix(samples,contPerct)
print(sampleMatrix)

#get ShortRead objects
sampleFastq<-getFastqObj(samples,dirFastq)
print(sampleFastq)

#Contaminate
simulateContamination(sampleMatrix,sampleFastq)
