setwd("/projectnb/bf528/students/srija/final_project")

library(tidyverse)

expression_data <- read.csv("/projectnb/bf528/users/frazzled/project_1/project_1.csv",sep = ",",header = TRUE,row.names = 1)

threshold_value <- log2(15)

filter_expressed <- rowSums(expression_data > threshold_value) >= (0.2*ncol(expression_data))
filter_set<- expression_data[filter_expressed, ]

row_var<- apply(filter_set, 1, var)

med_var <- median(row_var)
dof<- ncol(expression_data) - 1

filter1_chi <- filter_set[((dof*row_var)/med_var)>qchisq(((1 - 0.99)/2),dof,lower.tail = FALSE),]

filter2_var <-subset(filter1_chi, apply(filter1_chi, 1, function(x) sd(x)/mean(x)) > 0.186)

write.csv(filter_set, sep=',', file="data.csv")
write.csv(filter1_chi, sep=',', file="data.1.csv")
write.csv(filter2_var, sep=',', file="data.2.csv")

cluster<-hclust(dist(t(filter2_var)))
plot(cluster,main = "HCluster Dendrogram")

clusters<-cutree(h_cluster, 2)

print(c("Number of samples in cluster1:", sum(clusters==1)))
print(c("Number of samples in cluster2:", sum(clusters==2)))

metadata<- read.csv('/project/bf528/project_1/doc/proj_metadata.csv')
metadata<- subset(metadata, select=c(geo_accession, cit.coloncancermolecularsubtype))

color <- ifelse(metadata$cit.coloncancermoleculars == "C3", "red", "blue")

jpeg("/projectnb/bf528/students/srija/final_project/heatmap.jpeg")
heatmap(as.matrix(filter2_var), ColSideColors = color, main = "Heatmap of Gene Expression across all Samples") 
legend(x = "topright", legend = c("C3","Other"), fill = c("red", "blue"), title = "Subtype")
dev.off()

cluster1<-filter2_var[, clusters==1]
cluster2<-filter2_var[, clusters==2]
plot(cluster, hang = -1, cex = 0.8, xlab = "GEO Accession")
rect.hclust(cluster, k = 2, border = "red")

ttest <- apply(as.matrix(filter2_var),MARGIN=1,function(x) t.test(x=x[clusters==1],y=x[clusters==2]))

pvals <- sapply(ttest,function(x) x$p.value)
tstats <- sapply(ttest,function(x) x$statistic)
p_adjusted <- p.adjust(pvals ,method = "fdr")
t_df <- data.frame("Probeset_ID" = c(row.names(filter2_var)),
                     tstats,pvals,p_adjusted)

write.csv(t_df,sep = ",",row.names = FALSE,file = "analysis5.4.csv")

diff_expressed <- t_df$Probeset_ID[t_df$p_adjusted<0.05]
print(c("Number of probes left with adjusted p value < 0.05:",length(diff_expressed)))


#For biologist: perform same analysis on expression matrix from 4.5 and provide as csv
t_bio <- apply(as.matrix(filter1_chi),MARGIN=1,function(x) t.test(x=x[clusters==1],y=x[clusters==2]))
p_bio <- sapply(t_bio,function(x) x$p.value)
t_stat_bio <- sapply(t_bio,function(x) x$statistic)
p_adj_bio <- p.adjust(p_bio,"fdr")
biologist_file <- data.frame("Probeset_ID" = c(row.names(filter1_chi)),
                            t_stat_bio,p_bio,p_adj_bio)

diff_expressed_bio <- biologist_file$Probeset_ID[biologist_file$p_adj_bio<0.05]
numdiff <- length(diff_expressed_bio)
write.csv(biologist_file,row.names = F,file = "analysis5.6.csv",sep=",")
