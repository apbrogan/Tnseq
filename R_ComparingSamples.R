args = commandArgs(trailingOnly=TRUE)

Location_of_Coding_folder=args[1]
GENELIST_file=args[2]

library("data.table")
setwd(Location_of_Coding_folder)
SAMPLES_List=fread("Data/Samples_To_Compare.csv", sep=',')

for (i in 1:dim(SAMPLES_List)[1]){
  SAMPLE_name=SAMPLES_List[i,1]
  Control_name=SAMPLES_List[i,2]
  SAMPLE=fread(paste("Data/",SAMPLE_name,"-trimmed-bowtieMap-TA.txt",sep=""), sep='\t')
  Control=fread(paste("Data/",Control_name,"-trimmed-bowtieMap-TA.txt",sep=""), sep='\t')

  #Get the gene names
  Names_wdirections=fread(paste(GENELIST_file,sep=""), sep='\t',quote="\"")
  Names=Names_wdirections[,c(1,2,3,4)]
  Names=lapply(Names,as.character)
  Names$V5=as.numeric(Names$V4)-as.numeric(Names$V3)

  Directions=Names_wdirections[,c(1,2,5)]
  Directions=Directions[Directions$V5!="",]
  Directions=droplevels(Directions)

  Control=Control[-grep("IG_",Control$V4),]
  SAMPLE=SAMPLE[-grep("IG_",SAMPLE$V4),]

  #Make the Results file, can comment out if already exists
  SAMPLE_results=data.frame(matrix(nrow=length(unique(SAMPLE$V4)),ncol=15))
  for (i in 1:dim(SAMPLE_results)[1]){
    #Select the next gene locus
    Temp=unique(SAMPLE$V4)[i]
    #Find all of the "TA" sites with that gene locus in mutant
    Data_temp=SAMPLE[SAMPLE$V4==Temp,]
    #Find all of the "TA" sites with that gene locus in wild type/control
    Control_temp=Control[Control$V4==Temp,]
    
    SAMPLE_results[i,1]=as.character(Temp)     #Gene locus
    SAMPLE_results[i,2]=Names$V2[Names$V1==trimws(Temp)]  #Gene Name
    SAMPLE_results[i,3]=Names$V5[Names$V1==trimws(Temp)]  #Length
    SAMPLE_results[i,c(4,5)]=colSums(Data_temp[,c(5,6)]) #Mutant F and R reads
    SAMPLE_results[i,6]=sum(colSums(Data_temp[,c(5,6)])) #Mutant F and R reads
    SAMPLE_results[i,c(7,8)]=colSums(Control_temp[,c(5,6)])      #WT F and R reads
    SAMPLE_results[i,9]=sum(colSums(Control_temp[,c(5,6)])) #Mutant F and R reads
    SAMPLE_results[i,10]=round((SAMPLE_results[i,6]+1)/(SAMPLE_results[i,9]+1),4)  # Ratio
    a=suppressWarnings(wilcox.test(c(Data_temp$V5,Data_temp$V6),c(Control_temp$V5,Control_temp$V6)))
    SAMPLE_results[i,11]=round(a$p.value/2,8)   #Mann-Whitney or Wilcoxon test
    SAMPLE_results[i,12]=sum(Data_temp$V5>0)+sum(Data_temp$V6>0)         #Unique insertions Mutant
    SAMPLE_results[i,13]=SAMPLE_results[i,12]/SAMPLE_results[i,3]        #Unique insertions/length Mutant
    SAMPLE_results[i,14]=sum(Control_temp$V5>0)+sum(Control_temp$V6>0)   #Unique insertions WT
    SAMPLE_results[i,15]=SAMPLE_results[i,14]/SAMPLE_results[i,3]        #Unique insertions/length WT
  }

  SAMPLE_results$X2[is.na(SAMPLE_results$X2)]=" "

  colnames(SAMPLE_results) =  c('Loci', 'Name','length',
                               paste(SAMPLE_name, "F reads"),paste(SAMPLE_name, "R reads"),paste(SAMPLE_name, "Total reads"),
                               paste(Control_name,"F reads"),paste(Control_name,"R reads"),paste(Control_name,"Total reads"),
                               paste('Ratio (',SAMPLE_name,'/',Control_name,")"), "Mann Pvalue",
                               paste(SAMPLE_name,"Unique Insertions"),paste(SAMPLE_name,"Unique Insertions/Length"),
                               paste(Control_name,"Unique Insertions"),paste(SAMPLE_name,"Unique Insertions/Length"))

  #Remove those with no hits in either mutant or WT
  #SAMPLE_results=na.omit(SAMPLE_results)
  #Write the gene list
  fwrite(SAMPLE_results,paste("Data/",SAMPLE_name,"_", Control_name, "_results_genes.txt",sep=""),
         sep='\t', row.names=FALSE, col.names = TRUE, quote=FALSE)
}
