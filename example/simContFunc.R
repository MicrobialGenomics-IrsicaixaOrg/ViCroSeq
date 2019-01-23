getSamples<- function(qaSummary,totalReads=215000,n=4){
  # Get n samples to perform the contamination
  # The samples are filtered according to a certain number of reads(totalReads)
  #
  # Returns a vector with the sample names ordered in ascendent order
  
  readCounts<-qaSummary[['readCounts']]
  smallSamples<-readCounts[readCounts$read<totalReads,]
  if (length(smallSamples) > n) {
    smallSamples<-sample(smallSamples,n)
  }
  pairEndSamples<-row.names(smallSamples)
  samples<-sort(unique(sapply(pairEndSamples, function(x){sub('_L001_R[1|2]_001.fastq.gz','',x)},USE.NAMES = FALSE)))
  return (samples)
}

getMatrix<-function(samples,contPerct) {
  # Get the matrix with relationship between the samples to be contaminated and the samples used to contaminate them
  # and the percentage of the contamination
  #
  # Returns a dataframe with 2 rows (the samples source of the contamination and the percentage of contamination)
  # and n columns (the samples that are going to be contaminated)
  
  crossSamples<-c()
  samplesCont<-c()
  for (i in seq_along(samples)){
    j= if(i!= length(samples)) i+1 else 1
    crossSamples[i]<-samples[i]
    samplesCont[i]<-samples[j]
    
  }
  sampleMatrix<-rbind(crossSamples, contPerct)
  colnames(sampleMatrix)<-samplesCont
  return (sampleMatrix)
}

getFastqObj<- function(samples,dirFastq) {
  # Load the pair-end files for each sample into ShortRead objects
  # 
  # Returns a matrix of lists that relates each sample name with the corresponding pair of ShortRead objects
  return (sapply(samples, function(sampleName){fileToFastq(sampleName,dirFastq)}, USE.NAMES = TRUE))
}

fileToFastq<- function(sampleName,dirFastq) {
  # Obtain the pair-ended files (R1 and R2) from a directory and loading them into ShortRead Objects
  #
  # Returns a list with the structure list(R[1|2]=> ShortRead object)
  r1_r2<-list.files(dirFastq,pattern=paste(sampleName,'*'),full.names = TRUE)
  names(r1_r2)<-sapply(strsplit(r1_r2,'_'),`[`, 4)
  r1_r2_fastq<-sapply(r1_r2,readFastq,USE.NAMES=TRUE)

 
  return (r1_r2_fastq)
}

contR1 <- function(samplesR1,nReads) {
    # Performs the contamination of the sample reads R1 selecting a random subset of read ids from the sample to be contaminated
    # and from the sample for contamination. Calls the function contRn to perform the contamination.
    # Convert the ids of these reads to ids in R2 to substitute the same in the R2 ShortRead object.
    #
    # Returns a list with the sample reads R1 contaminated and the ids for R2 contamination
  
    idsR1<-sapply(names(samplesR1), function(x){ sample(id(samplesR1[[x]]),nReads)},USE.NAMES=TRUE)
    idsR2<-sapply(idsR1, function(x){ BStringSet(gsub('1:N','2:N',x))},USE.NAMES=TRUE)
    #Substitute the reads in the sampleToCont that match the id of reads$sampleToCont with the reads in reads$sampleForCont
   
    samplesR1cont<-contRn(samplesR1,idsR1)
    
    return (list('samplesContR1'=samplesR1cont, 'idsR2'=idsR2))
}


contRn <- function(samplesRn,idsRn) {
  #  Performs the contamination of the sample reads Rn by substituting the reads with ids in idsRn$sampleToCont for the reads
  #   with ids in idsRn$sampleForCont
  #  
  #  Returns the contaminated sample reads Rn as a ShortRead object
  masks<-sapply(names(samplesRn),function(name) which(id(samplesRn[[name]]) %in% idsRn[[name]]), USE.NAMES = TRUE)
  samplesRn$sampleToCont[masks[,'sampleToCont']]=samplesRn$sampleForCont[masks[,'sampleForCont']]
  return (samplesRn$sampleToCont)
}


contPairEndReads<-function(samples,nReads,txt){
  # Perform the contamination of a pair-end sample
  #
  # Returns a list with the R1 and R2 ShortObject contaminated and the text output generated

  result<-contR1(samples$R1,nReads)
  samplesContR2<-contRn(samples$R2,result$idsR2)
  
  #Generate output
  txt<-c(txt,'$sampleToCont ->ids corresponding of the corresponding R2 samples to be substituted with other samples')
  txt<-c(txt,as.character(result$idsR2$sampleToCont))
  txt<-c(txt,'$sampleForCont ->ids corresponding of the corresponding R2 samples to be inserted in other smple')
  txt<-c(txt,as.character(result$idsR2$sampleForCont))
  
  return (list('R1'=result$samplesContR1,'R2'=samplesContR2,'txt'=txt))
  
}

simulateContamination <- function(sampleMatrix,sampleFastq){
 
  for (sampleName in colnames(sampleMatrix)) {
    print(sampleName)
    #obtain R1-R2 ShortRead objects for each sample
    r1r2ToCont<-sampleFastq[,sampleName]
    r1r2ForCont<-sampleFastq[,sampleMatrix['crossSamples',sampleName]]
    samples<-list()
    for (key in c('R1','R2')){
      samples[[key]]<-list('sampleToCont'=r1r2ToCont[[key]], 'sampleForCont'=r1r2ForCont[[key]])
    }
    
    #obtain the number of reads to be sunstituted
    nReads<-as.integer(as.integer(sampleMatrix['contPerct',sampleName])*length(r1r2ToCont$R1)/100)
    #nReads<-2
    
    #Generate output
    txt<-c()
    txt<-c(txt,'Sample contaminated\tSource contamination')
    txt<-c(txt,paste(sampleName,sampleMatrix['crossSamples',sampleName],sep='\t'))
    txt<-c(txt,paste('Reads',nReads,sep='\t'))
    
    #run simulation
    samplesCont<-contPairEndReads(samples,nReads,txt)
    
    txt<-samplesCont$txt
    samplesCont$txt<-NULL
    
    #write Fastq contaminated files
    lapply(names(samplesCont),function(name) {
        fn=paste(sampleName,'cont',paste0(name,'.fastq.gz'),sep='_')
        if (file.exists(fn)) 
          #Delete file if it exist
          file.remove(fn)
        writeFastq(samplesCont[[name]],fn)
      })
    
    #write output file
    fileConn<-file(paste0("output",sampleName,".txt"))
    writeLines(txt,fileConn)
    close(fileConn)
  }
  

}