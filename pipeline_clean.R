setwd("/Users/macbook/Documents/demo_R")
library(GEOquery)
library(edgeR)
library(PCAtools)
library(biomaRt)
library(ReactomePA)
library(enrichplot)
library(umap)
library(ggplot2)
library(dbscan)
library(WGCNA)
library(genefilter)
library(flashClust)
library(RColorBrewer)
library(igraph)
#============================================================================
#
#  Download the files
#
#============================================================================
gse<-getGEOSuppFiles("GSE98212")
genes_count_all <- read_delim("~/Documents/demo_R/GSE98212/GSE98212_H_DE_genes_count.txt.gz",delim = "\t", escape_double = FALSE,trim_ws = TRUE)
##### selecting raw counts 
raw_counts<--GSE98212_H_DE_genes_count_txt[,2:31]
#============================================================================
#
#  DE analysis
#
#============================================================================

###edgeR
d=DGEList(counts=data_raw,group=factor(meta_high$Timepoint))

#####filtering no expressed genes
design.mat<-model.matrix(~0 + meta_high$Timepoint)
colnames(design.mat)=c("D0","D1")
keep<-filterByExpr(d,design.mat)
d_clean=d[keep,]

####calc norm factor 
d_clean<-calcNormFactors(d_clean)

###data exploration
BiocManager::install("PCAtools")
library(PCAtools)
norm_count<-cpm.DGEList(d_clean,log=TRUE)
p<-pca(norm_count,metadata = meta_high, removeVar = 0.1)
biplot(p,colby ='Timepoint',legendPosition = "right",lab=NULL)

#####LIMMA VOOM
y<-voom(d_clean,design.mat,plot=T)
fit_voom<-lmFit(y,design.mat) ###fit linear model
contr<-makeContrasts(D1-D0,levels=colnames(coef(fit_voom)))
tmp<-contrasts.fit(fit_voom,contr) ####contrast for each gene
tmp<-eBayes(tmp) 
table<-topTable(tmp,sort.by="P",n=5000)
####select genes FC>1.5 and FDR<0.05
tt <- abs(table$logFC) >= 0.6 & table$adj.P.Val<0.05
genes=table[tt,]
write.csv(genes,file="DEG_cl3vscl1.csv")


#####gene annotation

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
listAttributes(ensembl)
ens_name<-rownames(genes)
gene_name<-getBM(attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id"),filters="ensembl_gene_id",values=ens_name,mart = ensembl)
genes$ensembl_gene_id=rownames(genes)
df=left_join(gene_name,genes, by="ensembl_gene_id")


####sort based on logFC value
arrange(df, desc(logFC))

#============================================================================

 ####Reactome pathway gene set enrichment analysis

#============================================================================

geneList2<-df$logFC
names(geneList2)<-df$entrezgene_id
geneList2 = sort(geneList2, decreasing = TRUE)
de<-names(geneList2)
x <- enrichPathway(gene=de, pvalueCutoff = 0.05, readable=TRUE)
head(x)

####gse pathway
y <- gsePathway(geneList2, 
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE,)
head(y)

####visualization
dotplot(x, showCategory=30)
categorys <- c("DDX58/IFIH1-mediated induction of interferon-alpha/beta", "Negative regulators of DDX58/IFIH1 signaling",
               "ISG15 antiviral mechanism")

cnetplot(x, foldChange=geneList2,showCategory = categorys)

#============================================================================

####UMAP and K-mean clustering

#============================================================================
ump <- umap(norm_count, n_neighbors = 5, random_state = 123)

umap_df <- umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID=row_number())

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = ID,label=trial))+
  geom_point(size=3)+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot") + theme_minimal()+ theme(axis.title.x = element_text(size=12, color="white", face="bold",angle=0),axis.title.y.left = element_text(size=12, color="white", face="bold"),legend.title=element_text(color="white",face="bold",size=12),legend.text = element_text(color="white",face="bold",size=12))+ geom_text(hjust ="middle", nudge_x = 0.2)
ggsave("UMAP_plot.png")

####dbscan for k-means cluster definition

df_dbscan=umap_df[,c(1,2)]
cl <- hdbscan(df_dbscan, minPts = 5)
plot(df_dbscan, col=cl$cluster, pch=20)

####silhouette for assess cluster goodness

distance_matrix <- dist(umap_df[, 1:2])
silhouette <- silhouette(as.numeric(cl$cluster), dist = distance_matrix)
fviz_silhouette(silhouette)

###get membership of each clusters
umap_df$cluster = cl$cluster


####insert cluster column into metadata
####supplemental column 
meta_high$cluster=cl$cluster


####understand transcriptional differences between 
design.mat<-model.matrix(~0 + meta_high$cluster)
colnames(design.mat)<-c("cluster0","cluster2","cluster3")

y<-voom(d_clean,design.mat,plot=T)
fit_voom<-lmFit(y,design.mat) ###fit linear model
contr<-makeContrasts(cluster3-cluster2,levels=colnames(coef(fit_voom))) ###cluster3-cluster1
tmp<-contrasts.fit(fit_voom,contr) ####contrast for each gene
tmp<-eBayes(tmp) 
table<-topTable(tmp,sort.by="P",n=5000)
####select genes FC>1.5 and FDR<0.05
tt <- abs(table$logFC) >= 0.6 & table$adj.P.Val<0.05
genes=table[tt,]
genes$ensembl_gene_id=rownames(genes)
df=left_join(gene_name,genes, by="ensembl_gene_id")

write.csv(df,file="DEG_cl3vscl2.csv")

####perform Reactome 

library(ReactomePA)
geneList<-deg_cl3vs1$logFC
names(geneList)<-deg_cl3vs1$entrezgene_id
geneList = sort(geneList, decreasing = TRUE)
de<-names(geneList)
x <- enrichPathway(gene=de, pvalueCutoff = 0.05, readable=TRUE)
head(x)

####pathway visualization
dotplot(x, showCategory=200)
##identify more variable genes for 
##get norm count from dgelist
##transform matrix into expression set
nc<-log2(cpm(d_clean, normalized.lib.sizes=FALSE)+1)
exset=ExpressionSet(nc)
exset<-varFilter(exset, var.func=IQR, var.cutoff=0.7, filterByQuantile=TRUE)
mat=exprs(exset)

####ensemble gene annotation, see above
ens_name<-rownames(mat)
gene_name_wgcna<-getBM(attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id"),filters="ensembl_gene_id",values=ens_name,mart = ensembl)
mat_WGCNA=mat[gene_name_wgcna$ensembl_gene_id,]
all.equal(rownames(mat_WGCNA),gene_name_wgcna$ensembl_gene_id)

####WGCNA
allowWGCNAThreads() 

powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

sft = pickSoftThreshold(
  t(mat_WGCNA),# <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5 )

par(mfrow = c(1,2));
cex1 = 1;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type="n",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.80, col = "red",lty=3)
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

#####WGCNA 1

#build a adjacency "correlation" matrix
####construct TOM , Transforming the adjacency matrix in a topological overlap

softPower = 10
adjacency = adjacency(t(prova),power = softPower, type = "unsigned")  
TOM = TOMsimilarity(adjacency)                                  
dissTOM = 1-TOM                                                 


###Generate a clustered gene tree
library(flashClust)
par(mfrow = c(1,1))
geneTree = flashClust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

#This sets the minimum number of genes to cluster into a module
minModuleSize = 100

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

dynamicColors= labels2colors(dynamicMods)

####match color with eigengenes

genesMatchcolors <- cbind(dynamicColors,rownames(prova))

MEList= moduleEigengenes(t(prova), colors= dynamicColors,softPower = softPower)
MEs= MEList$eigengenes

MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")

#plots tree showing how the eigengenes cluster together
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

#set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0

MEDissThres = 0.6
merge = mergeCloseModules(t(prova), dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors = mergedColors
MEs = mergedMEs

####Step 3: Relate modules to external traits
#Define number of genes and samples
nGenes = nrow(prova)
nSamples = ncol(prova)
#Recalculate MEs with color labels
MEs0 = moduleEigengenes(t(prova), moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs,meta_for_corr, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)

#display the corelation values with a heatmap plot
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= colnames(moduleTraitCor),
               yLabels= names(MEs),
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.6,
               main= paste("Module-trait relationships"))

########get names of genes for each modules, from merged MES
rownames(mat_WGCNA)[moduleColors=="pink"] 

######forse meglio cosi
Group <- factor(metadataJ0$cluster, levels=c("1","2"))
design <- model.matrix(~0 +Group)
colnames(design)=c("cluster1","cluster2")
contrast_matrix <- makeContrasts(cluster1-cluster2, levels = design)


# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = design)

###apply contrast
tmp <- contrasts.fit(fit,contrast_matrix)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(tmp)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>% tibble::rownames_to_column("module")

head(stats_df)

##As a sanity check, let's use ggplot to see what module 18's eigengene looks like between treatment groups.
mod_boxplot=module_eigengenes
mod_boxplot$cluster=as.factor(metadata$time)
ggplot(mod_boxplot,aes(x = cluster,y = MEbrown,color=cluster))+geom_boxplot(width = 0.2, outlier.shape = NA) +theme_classic()
####network
nctable = as.matrix((prova[genes_net,]))
HARD.CUTOFF=0.8; VELDI=5
AdjMat <- cor(t(nctable), method="pearson") 
AdjMat[which(is.na(AdjMat))]=0
plot(density(AdjMat))
AdjMat[AdjMat <= HARD.CUTOFF]=0; AdjMat=AdjMat^VELDI
g <- graph.adjacency(AdjMat,mode="undirected",weighted=TRUE,diag=FALSE)
g <- delete_vertices(g, igraph::degree(g) == 0); length(V(g)$name)


##FIND CLUSTERS AND DELETE SMALL CLUSTERS
cluster.tmp <- cluster_louvain(g); table(cluster.tmp$membership)
cluster.rm <- which(table(cluster.tmp$membership) <= 8)
remove.tmp = V(g)$name[cluster.tmp$membership %in% cluster.rm]
g <- delete_vertices(g, remove.tmp)
length(V(g)$name)
##FIND CLUSTERS(LOUVAIN METHOD)
cluster <- cluster_louvain(g); table(cluster$membership)

# mst.communities cluster
library(RColorBrewer)
cluster.colors.tmp=brewer.pal(n = length(table(clustering$membership)), name = "Paired")
l <-layout.fruchterman.reingold(g)

###plot network
plot(g, mark.col=cluster.colors.tmp,
     vertex.color=cluster.colors.tmp[cluster$membership],
     layout=l,vertex.size=3, vertex.frame.color="gray10",
     vertex.label.color="transparent",
     edge.color="gray20",
     asp=FALSE,vertex.label.cex=1.0)
legend("topleft", xjust= 1,legend=paste0("cluster", 1:length(table(clustering$membership)),"  ", table(cluster$membership)), fill=cluster.colors.tmp, cex=0.5, title="Network Module   Nr. genes",bty = "n")

###genes in clusters
V(g2)$name[cluster$membership==6]


####PATHWAY ENRICHMENT
ORA.quick.results <- data.frame(
  module=factor(),term_id=character(), source=character(), term_name=character(),
  p_value=numeric(),intersection_size=integer(),term_size=integer(),query_size=integer(),
  stringsAsFactors=FALSE)
#GO
for(i in 1:length(table(clustering$membership))){ print(paste0("Analyzing cluster ",i))
  GenesConnected=V(g)$name[clustering$membership==i]
  rez <- gost(query=GenesConnected,
              organism = "hsapiens", ordered_query = F,
              significant = T,exclude_iea = F,
              measure_underrepresentation = F, evcodes = T,
              sources=c("GO:MF","GO:CC","GO:BP","REAC","KEGG","MIRNA"),
              user_threshold = 0.05,correction_method = "gSCS")
  rez <- rez$result
  rez <- rez[!(rez$intersection_size <= 2),]
  rez <- rez[(rez$term_size < 1100),]
  rez <- rez[(rez$term_size > 2),]
  NoResults=F
  if(is.null(dim(rez)[1])){NoResults=T}else if(dim(rez)[1]==0){NoResults=T}
  if(!NoResults){
    rez <- rez[, c("term_id", "source", "term_name", "p_value", "intersection_size", "term_size", "query_size","intersection")]
    rez <- rez[order(rez$p_value), ]
    res=cbind(rep(paste0("cluster",i),nrow(rez)),rez)
    colnames(res)[1]=c("module")
    ORA.quick.results=rbind(ORA.quick.results,res)
    print(res[1,]) 
  }else{print("No results")}
}; sum(-log10(ORA.quick.results$p_value)); 
View(ORA.quick.results)





