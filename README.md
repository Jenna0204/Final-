
Bioinformatics Assessment Task 3 (Part 1)
Collaborators- Janhavi Gowda and Hityshi Naidu


Gene Expression


The gene_expression.tsv file was downloaded from the following link: 

https://github.com/markziemann/SLE712_files/tree/master/assessment_task3/bioinfo_asst3_part1_files 

The following code was used to read the .tsv file into R
URL=“https://raw.githubusercontent.com/markziemann/SLE712_files/master/assessment_task3/bioinfo_asst3_part1_files/gene_expression.tsv”
NAME= "gene_expression.tsv"
download.file(URL,destfile=NAME)

To change the name of the rows to the names of the gene identifiers
head(read.delim("gene_expression.tsv"))  
str(read.delim("gene_expression.tsv"))
cds = read.delim("gene_expression.tsv")
head(cds)
row.names(cds)
row.names(cds) = cds$Name_Description
row.names(cds)

To show the values of the first 6 genes
head(cds,6)

To add a column with the means of the other columns 
head(rowMeans(cds[2:4])) 
cds$mean = rowMeans(cds[2:4])
head(cds,6)

To list 10 genes with the highest means 
order(-cds$mean)
sorted <- cds[order(-cds$mean),]
sorted[,c(4,ncol(sorted))]
head(sorted[,c(4,ncol(sorted))],10)

To list number of genes with a mean less than 10 
subset(cds,mean < 10)
nrow(subset(cds,mean < 10))
nrow(cds)

To plot a histogram 
hist(cds$mean,xlab="Mean Values", main="Histogram of Mean Values")



Growth Data

The file was downloaded for growth data using the below link

URL= "https://raw.githubusercontent.com/markziemann/SLE712_files/master/assessment_task3/bioinfo_asst3_part1_files/growth_data.csv"
NAME= "growth_data.csv"
download.file(URL,destfile=NAME)

This code read in the .csv file into R

read.csv("growth_data.csv")


For the following csv file downloaded the object name given as cdsA

cdsA = read.csv("growth_data.csv")
head(cdsA)

To display column names
colnames(cdsA)



Measuring circumference at the two sites northeast and southwest at the beginning and the end of the study

cdsA[cdsA$Site=="northeast",]
site1 = cdsA[cdsA$Site=="northeast",]
cdsA[cdsA$Site=="southwest",]
site2 = cdsA[cdsA$Site=="southwest",]
mean(site1[,3])
sd(site1[,3])
mean(site1[,6])
sd(site1[,6])
mean(site2[,3])
sd(site2[,3])
mean(site2[,6])
sd(site2[,6])




Boxplot for tree circumference at the start and end of the study at both sites.


boxplot(site1$Circumf_2005_cm, site1$Circumf_2020_cm,
        ylab="circumference",
        main="The northeast position starts and ends with the data",
        names = c("start","end"))

boxplot(site2$Circumf_2005_cm, site2$Circumf_2020_cm,
        ylab="circumference",
        main="The southwest position starts and ends with the data",
        names = c("start","end"))



site1gr = mean(site1[,6])-mean(site1[,4])
site1gr/mean(site1[,4])


site2gr = mean(site2[,6])-mean(site2[,4])
site2gr/mean(site2[,4])



Calculating the t-test and wilcox test

t.test(cdsA$Circumf_2020_cm - cdsA$Circumf_2010_cm~cdsA$Site)
boxplot(cdsA$Circumf_2020_cm - cdsA$Circumf_2010_cm~cdsA$Site,xlab = "location",ylab = "growth")
PVAL = t.test(cdsA$Circumf_2020_cm - cdsA$Circumf_2010_cm~cdsA$Site)$p.value
HEADER=paste("P:",PVAL)
mtext(HEADER)

wilcox.test(cdsA$Circumf_2020_cm - cdsA$Circumf_2010_cm~cdsA$Site)
boxplot(cdsA$Circumf_2020_cm - cdsA$Circumf_2010_cm~cdsA$Site,xlab = "location",ylab = "growth")
PVAL = wilcox.test(cdsA$Circumf_2020_cm - cdsA$Circumf_2010_cm~cdsA$Site)$p.value
HEADER=paste("P:",PVAL)
mtext(HEADER)







Bioinformatics Assessment Task 3 (Part 2)

#Ecoli

Downloads the ecoli file from ensemble bacteria website
library("R.utils")
URL="http://ftp.ensemblgenomes.org/pub/bacteria/release-53/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz"
download.file(URL,destfile="ecoli_cds.fa.gz")

Unzipsand lists the zipped files
gunzip("ecoli_cds.fa.gz")
list.files()

library("seqinr")
cds <- seqinr::read.fasta("ecoli_cds.fa")

Displays structure of the data
str(head(cds))

Gives length of coding sequence for ecoli
length(cds)

head(summary(cds))
str(summary(cds))
head(summary(cds)[,1])
len <- as.numeric(summary(cds)[,1])
sum(len)
lene <- sapply(X=cds,FUN=length)
sum(lene) 

Gives the mean of ecoli coding sequences
mean(len)

Gives the median of ecoli coding sequences
median(len)

Creates boxplot for ecoli
boxplot(len)
boxplot(len,ylab="sequence length (bp)")








#Acetobacter

library("R.utils")



#downloading acetobacter file 

URL="http://ftp.ensemblgenomes.org/pub/bacteria/release-53/fasta/bacteria_80_collection/acetobacter_aceti_gca_002005445/cds/Acetobacter_aceti_gca_002005445.ASM200544v1.cds.all.fa.gz"
download.file(URL,destfile="acetobacteraceti_cds1.fa.gz")
gunzip("acetobacteraceti_cds1.fa.gz")
list.files()

library("seqinr")
cds1 <- seqinr::read.fasta("acetobacteraceti_cds1.fa")
str(head(cds1))

#length of acetobaceter coding sequences 

length(cds1)

head(summary(cds1))

str(summary(cds1))

head(summary(cds1)[,1])

len1 <- as.numeric(summary(cds1)[,1])

sum(len1)

lena <- sapply(X=cds1,FUN=length)

sum(lena) 

#mean of acetobacter coding sequences 

mean(len1)

#median of acetobacter coding sequences 

median(len1)

boxplot(len1)

#boxlot for acetobacter

boxplot(len1,ylab="sequence length (bp)")
```




```{r}

#frequency of dna bases

GC(cds[[1]])

count(cds[[1]],1)

count(cds[[1]],2)

count(cds[[1]],3)

summary(cds[1:3])

sum(sapply(cds[1:3],length))

length(unlist(cds[1:3]))

dna <- unlist(cds)

GC(dna)

dna_composition <- count(dna,1)

#barplot for nucleotides

barplot(dna_composition,xlab="nucleotides",ylab="frequency", main="E coli CDS composition")

dna_composition/sum(dna_composition)

dna_proportion <- dna_composition/sum(dna_composition)

barplot(dna_proportion,xlab="nucleotides",ylab="proportion", main="E coli CDS composition")

grid()

translate(cds[[1]])

prot <- lapply(cds,translate)

# define the amino acid alphabet

aa <- unique(prot[[2]])
aa <- aa[aa != "*"]
length(aa)

count(prot[[1]],wordsize=1,alphabet=aa)

count(prot[[2]],wordsize=1,alphabet=aa)

#codon usage

uco(cds[[2]])

uco(cds[[2]],index="rscu")

uco(cds[[2]],index="rscu",as.data.frame=TRUE)

# k mer profiling

prots <- unlist(prot)

mycounts <- count(prots,wordsize=3,alphabet=aa)

str(mycounts)

head(mycounts)

myfreq <- count(prots,wordsize=3,alphabet=aa,freq=TRUE)
```



```{r}
GC(cds1[[1]])

count(cds1[[1]],1)

count(cds1[[1]],2)

count(cds1[[1]],3)

summary(cds1[1:3])

sum(sapply(cds1[1:3],length))

length(unlist(cds1[1:3]))

dna1 <- unlist(cds1)

GC(dna1)

dna1_composition <- count(dna1,1)

barplot(dna1_composition,xlab="nucleotides",ylab="frequency", main="Acetobacteraceti CDS composition")

dna1_composition/sum(dna1_composition)

dna1_proportion <- dna1_composition/sum(dna1_composition)

barplot(dna1_proportion,xlab="nucleotides",ylab="proportion", main="Acetobacteraceti CDS composition")

grid()

translate(cds1[[1]])

prot1 <- lapply(cds1,translate)

# define the amino acid alphabet
bb <- unique(prot1[[2]])



bb <- bb[bb != "*"]
length(bb)

count(prot1[[1]],wordsize=1,alphabet=bb)

count(prot1[[2]],wordsize=1,alphabet=bb)

uco(cds1[[2]])

uco(cds1[[2]],index="rscu")

uco(cds1[[2]],index="rscu",as.data.frame=TRUE)

#k mer profiling

protsAA <- unlist(prot1)

mycountsAA <- count(protsAA,wordsize=3,alphabet=bb)

str(mycountsAA)

head(mycountsAA)

myfreqAA <- count(protsAA, wordsize = 3, alphabet=bb, freq = TRUE)


