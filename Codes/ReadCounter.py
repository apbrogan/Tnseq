#Updated by JWS on 1_18_20
#Works for Any organism

#argv[1] is TA list (generated from TAfinder.py)
#argv[2] is mapped reads file (bowtie)
#argv[3] is gene list
#argv[4] is genome in fastA

import sys
#create dictionary of chrI TA sites

chrITA={}

for line in open(sys.argv[1]):  # opens files line by line
    split= line.split()         #splits string into columns--treates consecutive spaces as 1 single tab
    TAsite=int(split[0])
    chrITA[TAsite]=[0,0]

#counts reads at each TA site
for line in open(sys.argv[2]):
    split = line.split('\t')  #splits string into 'columns'
    strand = split[1]       #designate a marker for column 1 info
    if strand == '+':
        site = int(split[3]) + 15  #if on + strand, take column 3 position and add 15bp
        if site in chrITA: chrITA[site][0] += 1  #if read has been found before, tally 1 more in F reads
        site = int(split[3]) + 16  #if on + strand, take column 3 position and add 16bp
        if site in chrITA: chrITA[site][0] += 1  #if read has been found before, tally 1 more in F reads
    if strand == '-':
        site = int(split[3]) + 1 #if on - strand, add col 3 position to length of read in position 4 and add 1
        if site in chrITA: chrITA[site][1] += 1  #if read has been found before, tally 1 more in R reads
        site = int(split[3]) + 0 #if on - strand, add col 3 position to length of read in position 4
        if site in chrITA: chrITA[site][1] += 1  #if read has been found before, tally 1 more in R reads

# assign gene names

gene_dict = {}; prev_end = 0; IG = 'IG_1';

for line in open(sys.argv[3]):
    split = line.split('\t')
    gene = split[0]; start = int(split[2]); end = int(split[3]) #the 3rd and 4th column have the start and stop BP
    gene_dict[gene]=[start, end]
    IG = 'IG_' + gene
    if start > prev_end+1:
        gene_dict[IG]=[prev_end+1,start-1]
    prev_end = end

#Read the last line to get the last spot of the genome
#fileHandle = open (sys.argv[3],"r" )
#lineList = fileHandle.readlines()
#fileHandle.close()

#Read fasta to get genome size
fasta={}
fasta[1]=[]
for line in open(sys.argv[4]):
    split=line.split('\t')
    if split[0][0]!=">":
        read=split[0].rstrip()
        fasta[1].append(read)

genome="".join(fasta[1])
length_genome=len(genome)

gene_dict['IG_chrmEnd'] = [end-1,length_genome]

# store gene names into TA site dictionary

for k,v in gene_dict.iteritems():
    start = v[0]; end = v[1];
    for i in range(start,end+1):
        if i in chrITA: chrITA[i].append(k)


CI_sites = chrITA.keys()
CI_sites.sort()

for i in CI_sites:
    print 'B sub', '\t',  i, '\t',  i+1, '\t', chrITA[i][2], '\t', chrITA[i][0], '\t', chrITA[i][1]

