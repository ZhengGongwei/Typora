```R
library(reshape2)
library(ggplot2)

#breaks
bk <- c(seq(0,0.2,by=0.01))
# 做热图：
x <- read.csv("fdr_prop.csv",header=T,row.names=1)
colnames(x)
x$ID <- rownames(x)
df_m <- melt(x, id.vars=c("ID"))

p<- ggplot(data=df_m,aes(x=variable,y=ID))
p+geom_tile(alpha=0.5,aes(fill=value))+       #geom_tile绘制热图，使用value列作为颜色填充
  theme(axis.text.x=element_text(angle=45,vjust=0.5,size=14))+
theme(panel.grid = element_blank())+		#去除网格
 scale_fill_gradient(low = "white", high = "red")

```

