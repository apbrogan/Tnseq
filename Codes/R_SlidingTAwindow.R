args = commandArgs(trailingOnly=TRUE)

Location_of_Coding_folder=args[1]
GENELIST_file=args[2]

if (is.na(GENELIST_file)){
  GENELIST_file='/Users/jsher/Desktop/Generic_TN/B_scratch_v2/Codes/BC_CP034551_genelist.txt'
  Location_of_Coding_folder="/Users/jsher/Desktop/Generic_TN/B_scratch_v2/"
}

library("data.table")
setwd(Location_of_Coding_folder)
SAMPLES_List=fread("Data/Samples_To_Compare.csv", sep=',')

for (i in 1:dim(SAMPLES_List)[1]){
  SAMPLE_name=SAMPLES_List[i,1]
  Control_name=SAMPLES_List[i,2]
  SAMPLE=fread(paste("Data/",SAMPLE_name,"-trimmed-bowtieMap-TA.txt",sep=""), sep='\t')
  Control=fread(paste("Data/",Control_name,"-trimmed-bowtieMap-TA.txt",sep=""), sep='\t')

  Temp_Combined=SAMPLE[,5:6]+Control[,5:6]
  Temp_Combined[,"sum"]=Temp_Combined[,1]+Temp_Combined[,2]
  Temp_ind=Temp_Combined[,3]>1
  SAMPLE_gzero=SAMPLE[as.logical(unlist(Temp_ind)),]
  Control_gzero=Control[as.logical(unlist(Temp_ind)),]

  #Make the Results file, can comment out if already exists
  SAMPLE_results=data.frame(matrix(nrow=123091,ncol=11))

  Start=1
  Stop=13
  for (i in 1:(dim(SAMPLE_gzero)[1]-13)){
    Temp=Start:Stop 
    #Find all of the "TA" sites with that gene locus in mutant
    Data_temp=SAMPLE_gzero[Temp,]
    #Find all of the "TA" sites with that gene locus in wild type/control
    Control_temp=Control_gzero[Temp,]
    
    SAMPLE_results[i,1]=Start
    SAMPLE_results[i,2]=Stop
    SAMPLE_results[i,3]=Data_temp[1,2] #TA n+13
    SAMPLE_results[i,4]=Data_temp[13,2] #TA n+13
    SAMPLE_results[i,5]=Data_temp[13,3]-Data_temp[1,2]  #Length
    SAMPLE_results[i,c(6,7)]=colSums(Data_temp[,c(5,6)]) #Mutant F and R reads
    SAMPLE_results[i,8]=sum(colSums(Data_temp[,c(5,6)])) #Mutant total
    SAMPLE_results[i,c(9,10)]=colSums(Control_temp[,c(5,6)])      #WT F and R reads
    SAMPLE_results[i,11]=sum(colSums(Control_temp[,c(5,6)]))      #Mutant total
    SAMPLE_results[i,12]=round((SAMPLE_results[i,11]+1)/(SAMPLE_results[i,8]+1),4)  # Ratio
    a=suppressWarnings(wilcox.test(c(Data_temp$V5,Data_temp$V6),c(Control_temp$V5,Control_temp$V6)))
    SAMPLE_results[i,13]=round(a$p.value/2,8)   #Mann-Whitney or Wilcoxon test
    Start=Start+1
    Stop=Stop+1
  }

  SAMPLE_results$X2[is.na(SAMPLE_results$X2)]=" "

   colnames(SAMPLE_results) =  c("TA # Start", "TA # stop","TA location start", "TA location stop", "length",
                                paste(SAMPLE_name, "F reads"),paste(SAMPLE_name, "R reads"),paste(SAMPLE_name, "Total reads"),
                                paste(Control_name,"F reads"),paste(Control_name,"R reads"),paste(Control_name,"Total reads"),
                                paste('Ratio (',Control_name,"/",SAMPLE_name,")"), "Mann Pvalue")


  fwrite(SAMPLE_results,paste("Data/",SAMPLE_name,"_", Control_name,"_slidingTAwindow.txt",sep=""),
         sep='\t', row.names=FALSE, col.names = TRUE, quote=FALSE)
}