```R
##加载数据##
a=read.table('gene_fpkm.xls',sep='\t',quote = "",fill = T,comment.char = "!",header = T)#提取表达矩阵
rownames(a)=a[,1]
a <- a[,-1]

##将FPKM转换为TPM##
expMatrix <- a[,1:12]
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpms <- apply(expMatrix,2,fpkmToTpm)
tpms[1:3,]
colSums(tpms)

#输出结果：
> tpms[1:3,]
                        Si2      Si3      Si1     BNC1     BNC3     BNC2      NC3      NC2      NC1     BSi2     BSi1     BSi3
ENSMUSG00000022483 36003.43 40483.79 40296.43 18316.13 18665.45 22250.52 37316.75 39536.82 35459.51 22794.09 23689.20 17668.00
ENSMUSG00000040152 25208.74 26544.18 27070.71 22036.57 21364.16 22121.49 25923.29 26411.96 25881.98 23351.89 24449.19 24353.89
ENSMUSG00000026193 23157.16 20559.33 21928.39 19370.38 18832.59 19066.85 23447.31 21852.63 24317.04 19992.68 19183.88 18835.18
> colSums(tpms)
  Si2   Si3   Si1  BNC1  BNC3  BNC2   NC3   NC2   NC1  BSi2  BSi1  BSi3 
1e+06 1e+06 1e+06 1e+06 1e+06 1e+06 1e+06 1e+06 1e+06 1e+06 1e+06 1e+06 

##差异分析##

group_list=c(rep('Si',3),rep('BNC',3),rep('NC',3),rep('BSi',3))
## 强制限定顺序
group_list <- factor(group_list,levels = c("Si","BNC","NC","BSi"),ordered = F)
#表达矩阵数据校正
exprSet <- tpms
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
##![jepg1](D:\typora\jepg1.jpeg)##
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
#判断数据是否需要转换
exprSet <- log2(exprSet+1)
#差异分析：
dat <- exprSet
design=model.matrix(~0 + group_list)
y <- voom(dat, design, plot = T)
# lmFit fits a linear model using weighted least squares for each gene:
fit <- lmFit(y, design)
head(coef(fit))
#Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models:
# Comparison between times 6 and 9 for cultivar I5
# makeContrasts实际就是定义比较分组信息
contr <- makeContrasts(group_listSi - group_listBSi, levels = colnames(coef(fit)))
# 比较每个基因
tmp <- contrasts.fit(fit, contr)
# Empirical Bayes smoothing of standard errors ， (shrinks standard errors that are much larger or smaller than those from other genes towards the average standard error) 
# (see https://www.degruyter.com/doi/10.2202/1544-6115.1027)
fit <- eBayes(tmp)
# 使用plotSA 绘制了log2残差标准差与log-CPM均值的关系。平均log2残差标准差由水平蓝线标出
plotSA(fit, main="Final model: Mean-variance trend")
# topTable 列出差异显著基因
top.table <- topTable(fit, sort.by = "P", n = Inf)
# logFC: log2 fold change of I5.9/I5.6
# AveExpr: Average expression across all samples, in log2 CPM
# t: logFC divided by its standard error
# P.Value: Raw p-value (based on t) from test that logFC differs from 0
# adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
# B: log-odds that gene is DE (arguably less useful than the other columns)
head(top.table, 20)

# p值<0.05的基因有多少个
length(which(top.table$adj.P.Val < 0.05))

#Write top.table to a file
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
library(clusterProfiler)
library(org.Mm.eg.db)
df <- bitr(unique(top.table$Gene), fromType = "ENSEMBL",
     toType = c( "SYMBOL"),
     OrgDb = org.Mm.eg.db)
head(df)
DEG=top.table
head(DEG)

DEG=merge(DEG,df,by.y='ENSEMBL',by.x='Gene')
head(DEG)
write.table(DEG, file = "Si_v_BSi.txt", row.names = F, sep = "\t", quote = F)

#
contr <- makeContrasts(group_listSi - group_listBSi, levels = colnames(coef(fit)))
```