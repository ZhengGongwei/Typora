```R
metadata <- pbmc@meta.data
metadata$cellname <- rownames(metadata)
write.table(metadata,file="metadata.txt")

#"/share2/pub/zhenggw/zhenggw/project1/AdultHemOrgans/data/RData"
metadata <- read.table("metadata.txt")

path <- "/share2/pub/zhenggw/zhenggw/project1/intest/data/RData/"
fileNames <- dir(path,"^Bac")
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')
})

for(i in 1:length(filePath)){
    result <- strsplit(fileNames[i],split=".genes")[[1]][1]
    print(result)
    name <- paste(result,"RData",sep=".")
    result <- read.table(filePath[i])
    result$cellname <- rownames(result)
    mergeinfo <- merge(result,metadata,by="cellname")
    finalinfo <- mergeinfo[order(mergeinfo$pval),]
    finalinfo$fdr <- p.adjust(finalinfo$pval,method="fdr")
    print(head(finalinfo))
    save(finalinfo,file=name)
}
```

