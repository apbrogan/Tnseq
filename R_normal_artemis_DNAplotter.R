args = commandArgs(trailingOnly=TRUE)

Location_of_Coding_folder=args[1]
SAMPLE_name=args[2]
Total_reads=as.numeric(args[3])

library("data.table")
setwd(Location_of_Coding_folder)
SAMPLE=fread(paste("Data/",SAMPLE_name,"-trimmed-bowtieMap-TA.txt",sep=""), sep='\t')

#Figure out how many reads per strain
SAMPLEsum=sum(SAMPLE[,c(5,6)])

#Figure out which is better and firgure out multiplying factor
Multiply_factor=Total_reads/SAMPLEsum


#Normalize the datasets so they all have the same number of reads
SAMPLE[,c(5:6)]=round(SAMPLE[,c(5:6)]*Multiply_factor[1],digit=0)

TA=dim(SAMPLE)[1]
Reads=sum(SAMPLE$V5)+sum(SAMPLE$V6)
Hit=sum(SAMPLE$V5>0|SAMPLE$V6>0)
PercentTA=round(Hit/TA*100,0)
AverageCount=round(Reads/TA,0)

cat(paste(SAMPLE_name), sep="\n")
cat(paste("TA= ", TA), sep="\n")
cat(paste("Reads=", Reads), sep="\n")
cat(paste("Sites Hit = ",Hit), sep="\n")
cat(paste("Percent Tas  Hit = ",PercentTA), sep="\n")
cat(paste("Average Read Count = ", AverageCount), sep="\n")


fwrite(SAMPLE,paste("Data/",SAMPLE_name,"-trimmed-bowtieMap-TA-normalized.txt",sep=""),
       sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)

Artemis_sample=SAMPLE[,c(3,5:6)]
Artemis_sample[,3]=-1*Artemis_sample[,3]
Header=data.frame(matrix(c("colour","5:150:55","225:0:0","#BASE","Finsert","Rinsert"), 
                         nrow=2, ncol=3, byrow=TRUE))
fwrite(Header,paste("Data/",SAMPLE_name,"_artemismap.txt",sep=""),
       sep='\t', row.names=FALSE, col.names = FALSE, quote=FALSE)
fwrite(Artemis_sample,paste("Data/",SAMPLE_name,"_artemismap.txt",sep=""),
       sep='\t', row.names=FALSE, col.names = FALSE, quote=FALSE, append=TRUE)


REVERSE=as.numeric(unlist(Artemis_sample[,3]))*-1
FORWARD=as.numeric(unlist(Artemis_sample[,2]))
TOTAL=FORWARD+REVERSE
fwrite(as.data.frame(TOTAL),paste("Data/",SAMPLE_name,"_DNAplotter.txt",sep=""), sep='\t', row.names=FALSE, col.names = FALSE, quote=FALSE)





