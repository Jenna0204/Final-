
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





