


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
