


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
