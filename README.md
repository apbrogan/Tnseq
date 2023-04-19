# Mariner TN-Seq pipeline by Joel Sher

Updated 1_22_19

Instructions for using pipeline.  If this is your first time start here, if not skip to step 3

1.	Programs needed:
a.	Python2 – download the latest version of Python2 from HERE (right now it is Python 2.7)
b.	R from HERE 
c.	Install the R package called “data.table”
i.	Once you have R installed, open R and type the command “install.packages(“data.table”)” 
d.	Need bowtie-1.0.0 from HERE 

2.	Need the following files for your organism
a.	Fasta file of your genome in the folder called Codes– must have extension .fa
b.	Bowtie indexes for your organism inside the bowtie-1.0.0/indexes
i.	If you do not have them can build using the following commands
ii.	Download FASTA sequence of your reference genome. Change the extension manually from .fasta to .fa (For E coli., I use NC_000913.fa)
iii.	Open Terminal (open spotlight search and type terminal)
iv.	Into the command line type: 
v.	jsher/Desktop/Generic_TN/bowtie-1.0.0/bowtie-build -f / jsher/Desktop/Generic_TN/Ecoli/Codes/NC_000913.fa MG1655Genome
c.	Gene list for your organism with the following 5 columns in a text file located in a folder called Codes
i.	1 unique identifier for each gene
ii.	Name of name if available/can copy the unique identifier if nothing else is known about the gene
iii.	Start of gene
iv.	End of the gene


3.	Make sure you have all of the required files
a.	Copy and paste the top level folder.  It should contain
i.	1 file called Universal_Shell.sh
ii.	1 folder called bowtie-1.0.0
1.	Inside should be a folder called indexes which contains the 6 files created in step 2B above
2.	The 6 files should end in 1.ebwt, 2.ebwt, 3.ebwt, 4.ebwt, rev.1.ebwt, rev.2.ewt
iii.	1 folder called Codes
1.	Must contain your fastA genome from step 2A above
2.	Must contain your genelist from step 2C above
iv.	A folder called Data containing your trimmed sequencing files (note: the files must end in -trimmed.fastq”

4.	Editing the Universal_Shell.sh
a.	Open the Universal_Shell.sh in a text editor program (I like sublime text)
b.	You must edit line 2 to the path of the top level folder described in 3A above
c.	You must edit line 3 to the location of your Rscript file (see my example for where it is commonly located)
d.	Edit lines 6-8 to reflect the name of your bowtie indexes, fastA and genelist file
i.	Double check the extensions match exactly for the fastA and genelist
e.	On line 17 add the names of all your samples after the “in”
i.	For example if you have 2 samples called WT and mutant the line of code would read for SAMPLE in WT mutant
f.	Save and close this file

5.	If you want the code to perform pairwise comparisons of samples (i.e. get ratios and p-values) then
a.	Open the file called “Samples_To_Compare.csv” located in the Data/ folder using excel.
b.	Write any pair of samples you want compared on a single row
i.	For example, if you want to compare WT to mut1 and mut2 the file would look like
Sample1	Sample2
WT	mut1
WT	mu2




c.	Save and close this file

6.	Running pipeline
a.	Open Terminal (open spotlight search and type terminal) 
b.	Type “cd” and a space.  Then drag the top level folder into the terminal.  The file path leading to the terminal should appear
i.	i.e. “cd /Users/jsher/Desktop/Generic_TN/Ecoli”
c.	type “Bash Universal_Shell.sh” and hit enter

