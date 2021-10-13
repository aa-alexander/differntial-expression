#deseq for 3 conditions 
#control vs treated 1 vs treated 2

library(DESeq2)
library(dplyr)

#reading count table
count.table <- read.table(file = "control_treat1_treat2_gene_count.tsv", header=T, row.names=1)

#creating colData
count.table %>% colnames -> sample
conditions <- c("control", "control", "control", "treated1", "treated1", "treated1", "treated2", "treated2", "treated2")
colData <- cbind(sample, conditions)
rownames(colData) <- colData[, 1]
colData <- as.data.frame(colData)

colData$conditions <- as.factor(colData$conditions)

#importing in deseq
dds <- DESeqDataSetFromMatrix(countData=count.table, colData=colData, design= ~conditions)

#releveling with control as reference
dds$conditions <- relevel(dds$conditions, ref="control")

#performing DESeq2
dds <- DESeq(dds)

#getting results for 3 probabilities
controlvstreat1.res <- results(dds, contrast=c("conditions", "control", "treated1"))
controlvstreat2.res <- results(dds, contrast=c("conditions", "control", "treated2"))
treat1vstreat2.res <- results(dds, contrast=c("conditions", "treated1", "treated2"))

# remove any rows with NA
controlvstreat1.res <- controlvstreat1.res[complete.cases(controlvstreat1.res), ]
controlvstreat2.res <- controlvstreat2.res[complete.cases(controlvstreat2.res), ]
treat1vstreat2.res <- treat1vstreat2.res[complete.cases(treat1vstreat2.res), ]

#making them as dataframe to sort them by pvalue 
controlvstreat1.res %>% as.data.frame -> dea.controlvstreat1
controlvstreat2.res %>% as.data.frame -> dea.controlvstreat2
treat1vstreat2.res %>% as.data.frame -> dea.treat1vstreat2

# sort the table: ascending of pvalue then descending of absolute valued of logFC
dea.controlvstreat1 <- dea.controlvstreat1[order(dea.controlvstreat1$pvalue, -abs(dea.controlvstreat1$log2FoldChange), decreasing = FALSE), ]
dea.controlvstreat2 <- dea.controlvstreat2[order(dea.controlvstreat2$pvalue, -abs(dea.controlvstreat2$log2FoldChange), decreasing = FALSE), ]
dea.treat1vstreat2 <- dea.treat1vstreat2[order(dea.treat1vstreat2$pvalue, -abs(dea.treat1vstreat2$log2FoldChange), decreasing = FALSE), ]

#saving the results to a file
write.table(dea.controlvstreat1,"dea_controlvstreat1.tsv", row.names = T, quote = FALSE, sep = '\t')
write.table(dea.controlvstreat2,"dea_controlvstreat2.tsv", row.names = T, quote = FALSE, sep = '\t')
write.table(dea.treat1vstreat2,"dea_treat1vstreat2.tsv", row.names = T, quote = FALSE, sep = '\t')
