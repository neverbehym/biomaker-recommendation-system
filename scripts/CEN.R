library(igraph)
library(WGCNA)
library(ggplot2)
library(topGO)
library(reshape2)
library(viridis)
library(centiserve)
library(expm)

############# read a pool of candidate biomarker panels #############

panel_log <- read.table('../results/CandidateBiomarkerList.csv', sep=",", header = TRUE,stringsAsFactors = FALSE)

panel_list <- strsplit(panel_log$genes,';')

biomarker_union<-unique(unlist(panel_list))
length(biomarker_union)

############# read gene expression data  #############
m_express <- as.matrix(read.table("../data/NicolasGeneExpression_processed.csv", sep=",", header = TRUE, stringsAsFactors = FALSE,row.names = 1))

############# Construct of Co-Expression Network #############

### Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
#  For each power the scale free topology fit index is calculated and returned along with other information on connectivity
#  [powerEstimate, fitIndices] the lowest power for which the scale free topology fit R^2 exceeds RsquaredCut is returned
sft = pickSoftThreshold(m_express, RsquaredCut = 0.85, powerVector = powers, verbose = 5)
sft$fitIndices
# Plot the results:
sizeGrWindow(12, 5)
par(mfrow = c(1,3));
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = 'Scale independence',cex.axis = 1.2,cex.lab = 1.5);
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1.2,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"),cex.axis = 1.2,cex.lab = 1.5)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=1.2,col="red")
# The function plots a log-log plot of a histogram of the given connectivities, 
# and fits a linear model plus optionally a truncated exponential model. The R^2 of the fit can be considered an index of the scale freedom of the network topology.
scaleFreePlot(k, truncated=FALSE,cex.axis = 1.5,cex.lab = 1.4)

### calculate the adjacencies
adjacency_CEN = adjacency(m_express, power = sft$powerEstimate) 
dimnames(adjacency_CEN)<-list(gene_list,gene_list)
write.table(adjacency_CEN,'../results/adj_mat.csv',col.names = TRUE,row.names = TRUE,quote=FALSE,sep=',')

############# topological properties #############

### Topological Overlap Map (TOM)
TOM_CEN = TOMsimilarity(adjacency_CEN);
dissTOM_CEN = 1-TOM_CEN
dimnames(TOM_CEN)<-list(gene_list,gene_list)
dimnames(dissTOM_CEN)<-list(gene_list,gene_list)

### connectivity centrality, i.e. node degree
k=softConnectivity(m_express,power=sft$powerEstimate) # node connectivity
names(k)<-gene_list
quantile(k, prob=seq(0,1,0.1))

### weights
weights <- as.vector(adjacency_CEN[lower.tri(adjacency_CEN, diag = FALSE)])
length(weights) 
quantile(weights, prob=seq(0,1,0.1))

### max weights of each node
node_maxweights<-rep(0,num_genes)
for(i in 1:num_genes){
  x<-adjacency_CEN[i,]
  x<-x[-i]
  node_maxweights[i]<-max(x)
}
names(node_maxweights)<-gene_list
quantile(node_maxweights, prob=seq(0,1,0.1))

# Plot a histogram of k and a scale free topology plot
par(mfrow=c(1,3))
hist(k,cex.axis = 1.5,cex.lab = 1.4,cex.main=1.4,xlab='Node connectivity (K)',ylab='Count')
hist.data = hist(weights, plot=F)
hist.data$counts = log(hist.data$counts, 10)
plot(hist.data,cex.axis = 1.5,cex.lab = 1.4,cex.main=1.4,xlab='Edge weight',ylab='log10(Count)')
hist(node_maxweights,cex.axis = 1.5,cex.lab = 1.4,cex.main=1.4,xlab='Node maximum weight',ylab='Count')


############# module detection #############
# scan parameters to maximally reduce the unassigned genes (grey), include the biomarker genes and increase the consensus scores 
regnet_regulator<-read.table("../data/GeneRegulatoryNetwork.csv",quote="", sep="\t",header = TRUE, stringsAsFactors = FALSE)
colnames(regnet_regulator)<-c('target','source','sign')

regulon_lists<-list()
for (regulator in unique(regnet_regulator$source))
  regulon_lists[[regulator]]<-regnet_regulator$target[regnet_regulator$source==regulator]

detectCutHeight_list=vector() 
mergeCutHeight_list=vector() 
minModuleSize_list=vector() 

ConsensusScore_list = vector()

detectCutHeight_scan <- c(0.98,0.985,0.99,0.995)
mergeCutHeight_scan <- c(0.15,0.2, 0.25,0.3)
minModuleSize_scan <- c(3,4,5,8,10,15,20)
maxConsensusScore<-0
for (I in 1:length(detectCutHeight_scan)){
  for (J in 1:length(mergeCutHeight_scan)){
    for (K in 1:length(minModuleSize_scan)){

      net = blockwiseModules(m_express, power = sft$powerEstimate, corType = "pearson",TOM_CENType = "unsigned",  maxBlockSize = 3000,
                             deepSplit = 2,detectCutHeight = detectCutHeight_scan[I], minModuleSize = minModuleSize_scan[K],
                             reassignThreshold = 1e-6, mergeCutHeight = mergeCutHeight_scan[J],
                             pamStage = TRUE, pamRespectsDendro = FALSE,
                             numericLabels = TRUE,  saveTOM_CENs = FALSE,verbose = 0)
      
       # module membership
      moduleLabels = net$colors
      moduleColors = labels2colors(net$colors)
      names(moduleColors)<-names(moduleLabels)
      num_modules<-length(table(moduleLabels))-1
      
      #### consensus with regulon list
      module_lists<-list()
      for (module_temp in unique(moduleColors))
        module_lists[[module_temp]]<-names(moduleColors)[moduleColors==module_temp]
      
      maxDiceScore<-rep(0,length(module_lists))
      names(maxDiceScore)<-names(module_lists)
      for (i in 1:length(module_lists)){
        for (j in 1:length(regulon_lists)){
          Consensus_Num<-length(intersect(module_lists[[i]],regulon_lists[[j]]))
          DiceScore<-2*length(intersect(module_lists[[i]],regulon_lists[[j]]))/(length(module_lists[[i]])+length(regulon_lists[[j]]))
          
          maxDiceScore[i]<-max(maxDiceScore[i],DiceScore)
        }
      }
      ConsensusScore <- mean(maxDiceScore)
      
      detectCutHeight_list = c(detectCutHeight_list,detectCutHeight_scan[I])
      mergeCutHeight_list = c(mergeCutHeight_list,mergeCutHeight_scan[J])
      minModuleSize_list = c(minModuleSize_list,minModuleSize_scan[K])
      ConsensusScore_list = c(ConsensusScore_list,ConsensusScore)
      if (ConsensusScore>maxConsensusScore){
        maxConsensusScore<-ConsensusScore
        detectCutHeight_chosen<-detectCutHeight_scan[I]
        mergeCutHeight_chosen<-mergeCutHeight_scan[J]
        minModuleSize_chosen<-minModuleSize_scan[K]
      }
    }

  }
}
print(paste0('Best consensus ccore=',maxConsensusScore,' achieved when detectCutHeight=',detectCutHeight_chosen,' mergeCutHeight_=',mergeCutHeight_chosen,' minModuleSize=',minModuleSize_chosen))
df_modIdTune<-data.frame(detectCutHeight=detectCutHeight_list,mergeCutHeight=mergeCutHeight_list,minModuleSize=minModuleSize_list,ConsensusScore=ConsensusScore_list)
write.table(df_modIdTune,'../results/modIdTune_summary.csv',col.names = TRUE,row.names = TRUE,quote=FALSE,sep=',')


net = blockwiseModules(m_express, power = sft$powerEstimate, corType = "pearson",TOM_CENType = "unsigned",  maxBlockSize = 3000,
                                              deepSplit = 2,detectCutHeight = 0.98, minModuleSize = 4,
                                              reassignThreshold = 1e-6, mergeCutHeight = 0.3,
                                              pamStage = TRUE, pamRespectsDendro = FALSE,
                                              numericLabels = TRUE,  saveTOM_CENs = FALSE,verbose = 0)
table(net$colors) # branch number assigned to each gene, 0 means nonassigned

### module membership
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
names(moduleColors)<-names(moduleLabels)
num_modules<-length(table(moduleLabels))-1
write.table(moduleColors, "../results/moduleColors.csv", append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = F)

### print biomarker genes that are not assigned to any modules
greyModule_genes<- biomarker_union[moduleColors[biomarker_union]=='grey'] 
for (i in 1:length(panel_list)){
 panel_temp<-panel_list[[i]]
 if (sum(greyModule_genes%in% panel_temp))
  print(paste0(c('In Panel',i, 'biomakrer genes ',panel_temp[panel_temp %in%greyModule_genes], 'are not assigned to any modules'), collapse=" "))
 
}

### consensus with regulon list
module_lists<-list()
for (module_temp in unique(moduleColors))
  module_lists[[module_temp]]<-names(moduleColors)[moduleColors==module_temp]                       
   
DiceScores<-matrix(0,nrow=length(module_lists),ncol=length(regulon_lists),dimnames=list(names(module_lists),names(regulon_lists)))
PercentsInModule<-matrix(0,nrow=length(module_lists),ncol=length(regulon_lists),dimnames=list(names(module_lists),names(regulon_lists)))
CountsInModule<-matrix(0,nrow=length(module_lists),ncol=length(regulon_lists),dimnames=list(names(module_lists),names(regulon_lists)))
PercentsInRegulon<-matrix(0,nrow=length(module_lists),ncol=length(regulon_lists),dimnames=list(names(module_lists),names(regulon_lists)))

ConsensusScore <- mean(matrixStats::rowMaxs(DiceScores))

# open a graphics window
sizeGrWindow(12, 6)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], colors=moduleColors[net$blockGenes[[1]]],"Module colors", dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)

############# module profiles #############

###  module size and density (average weights (correlation strength) within module#### 
table(moduleColors)
module_size<-table(moduleColors)

module_density<-rep(0,num_modules+1)
names(module_density)<-labels2colors(0:num_modules)
for (i in 0:num_modules){
  module_genes<-names(moduleLabels)[moduleLabels==i]
  module_adj<-adjacency_CEN[module_genes,module_genes]
  module_density[i+1]<-mean(module_adj[lower.tri(module_adj, diag = FALSE)])
}
module_density<-sort(module_density)
names(tail(module_density))

module_consensus<- matrixStats::rowMaxs(DiceScores)
names(module_consensus)<-names(module_lists)
module_consensus<-sort(module_consensus)

### barplot of module size and module density
par(mar=c(5.1,4.1,1,1))  #(bottom, left, top, right)
par(mfrow = c(3,1))
labels<-rep("",length(module_size))
labels[seq(2,length(x),2)]= module_size[seq(2,length(x),2)]
barp<-barplot(module_size,xaxt = 'n', xlab = '',ylab='module size',ylim=c(0,max(module_size)+200))
x <-substr(as.character(module_size),start=1,stop=4)
text(x = barp, y = module_size, label = labels, pos = 3, cex = 0.8, col = "red")
axis(1, at=barp, labels=names(module_size), tick=FALSE, las=2, line=-0.5, cex.axis=0.8)
barp<-barplot(module_density,xaxt = 'n', xlab = '',ylab='module density',ylim=c(0,max(module_density)+0.1))
x <-substr(as.character(module_density),start=1,stop=4)
labels<-rep("",length(x))
labels[seq(2,length(x),2)]<-x[seq(2,length(x),2)]
text(x = barp, y = module_density, label = labels, pos = 3, cex = 0.8, col = "red")
axis(1, at=barp, labels=names(module_density), tick=FALSE, las=2, line=-0.5, cex.axis=0.8)
barp<-barplot(module_consensus,xaxt = 'n', xlab = '',ylab='module consensus score',ylim=c(0,max(module_consensus)+0.1))
x <-substr(as.character(module_consensus),start=1,stop=4)
labels<-rep("",length(x))
labels[seq(2,length(x),2)]<-x[seq(2,length(x),2)]
text(x = barp, y = module_consensus, label = labels, pos = 3, cex = 0.8, col = "red")
axis(1, at=barp, labels=names(module_consensus), tick=FALSE, las=2, line=-0.5, cex.axis=0.8)


### module eigengenes of the found modules (given by colors)
MEs = net$MEs 
#geneTree = net$dendrograms[[1]];

### K_ME and K_IM ####  
### eigengene-based connectivity, i.e, module fuzzy membership
# return a data frame [input genes, module eigengenes]
K_ME<-signedKME(m_express, MEs)
colnames(K_ME)<- labels2colors(as.numeric(stringr::str_remove(colnames(K_ME),'kME')))
colnames(MEs)<-labels2colors(as.numeric(stringr::str_remove(colnames(MEs),'ME')))
#### intramodular connectivity
# return a 4-columns data frame of total, intramodular, extra-modular connectivity, and the difference of  intra- and extra-modular connectivities
K_IM<-intramodularConnectivity(adjacency_CEN, moduleColors, scaleByMax = FALSE)

df_modulegens<-data.frame(gene=names(moduleColors),ModuleColor=moduleColors,K_ME=rep(0,length(moduleColors)),K_IM=rep(0,length(moduleColors)))
for(i in 1:length(moduleColors)){                                                                                     
  df_modulegens$K_ME[i]<-K_ME[df_modulegens$gene[i],df_modulegens$ModuleColor[i]]
  df_modulegens$K_IM[i]<-K_IM[df_modulegens$gene[i],'kWithin']
}

write.table(df_modulegens, paste0("../results/modulegenes.csv"), append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = F)


df_moduleProfile <- data.frame(module_size)
names(df_moduleProfile)<-c('moduleColor','moduleSize')

CentralGenes = vector()
TopRegulators = vector()
PercentRegulators = vector()
for (module_name in df_moduleProfile$moduleColor){
  df_thismodule <- df_modulegens[df_modulegens$ModuleColor == module_name,]
  centralgenes = df_thismodule$gene[order(df_thismodule$K_IM,decreasing = TRUE)[1:3]] 
  CentralGenes =c(CentralGenes,paste0(centralgenes,collapse=';') )
  
  thismodule_dices <- DiceScores[module_name,]
  tops <- sort(thismodule_dices[thismodule_dices>0.15], decreasing = TRUE)
  if (length(tops)==0)
    tops <- sort(thismodule_dices[thismodule_dices>0.05], decreasing = TRUE)
  if (length(tops)==0)
    tops <- sort(thismodule_dices[thismodule_dices>0], decreasing = TRUE)
  if (length(tops)==0){
    TopRegulators =c(TopRegulators,'' )
    PercentRegulators =c(PercentRegulators,'' )
  }
  else{
    topregulators = names(tops)
    percentregulators = unname(tops)
    TopRegulators =c(TopRegulators,paste0(topregulators,collapse=';') )
    PercentRegulators =c(PercentRegulators,paste0(round(percentregulators*100,2),"%",collapse=';') )
    }
  
}

df_moduleProfile$CentralGenes <- CentralGenes
df_moduleProfile$TopRegulators <- TopRegulators
df_moduleProfile$PercentRegulators <- PercentRegulators

write.table(df_moduleProfile, "../results/module_profiles.csv", append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = F)

############# igraph of the CEN #############
adjacency_CEN_invert<-mean(weights)/adjacency_CEN
diag(adjacency_CEN_invert)<-0
write.table(adjacency_CEN_invert,'../results/adj_mat_invert.csv',col.names = TRUE,row.names = TRUE,quote=FALSE,sep=',')

g_coexpgnet<- graph_from_adjacency_matrix(adjacency_CEN_invert, weighted=TRUE,mode='undirected',diag = FALSE)
weights_list<-E(g_coexpgnet)$weight
quantile(weights_list)
length(weights_list[weights_list<0.01])

############# obtain the network metrics #############
### degree
dg<-degree(g_coexpgnet)

### distances of all pairwise nodes
distance_table<-distances(g_coexpgnet)

### betweenness centrality
betweenness_cent<-estimate_betweenness(g_coexpgnet, cutoff=max(distance_table),directed = FALSE)

### latora closeness centrality
latoraclose_cent<-closeness.latora(g_coexpgnet)

df_vp_CEN<-data.frame(betweenness=betweenness_cent,closeness=latoraclose_cent,degree=dg)

############# compute the evaluation metrics for biomarker panels #############

## metrics for each panel for size 10
num_biomarkers=lengths(panel_list)
num_modules=vector()

avg_dgs=vector() # average node degrees of biomarkers
avg_bet=vector()  # average betweenness centrality of biomarkers
avg_clo=vector() # average closeness centrality of biomarkers

avg_dts=vector() # average pairwise distance between biomarkers

directmodules=list()
 for (i in 1:length(panel_list)){
   
    panel<-panel_list[[i]]
    
    module_table <- table(moduleColors[panel])
    directmodules[[paste0('panel',i)]]= module_table
    num_modules<-c(num_modules,length(module_table))
    
    avg_dt<- sum(distance_table[panel,panel])/length(panel)/(length(panel)-1)
    avg_dts<-c(avg_dts,avg_dt)
    
    avg_dgs<-c(avg_dgs,mean(df_vp_CEN[panel,'degree']))
    avg_bet<-c(avg_bet,mean(df_vp_CEN[panel,'betweeness']))
    avg_clo<-c(avg_clo,mean(df_vp_CEN[panel,'closeness']))
    
    
  }



df_metrics_biomarkers$pairwise_distance_CEN=avg_dts
df_metrics_biomarkers$degree_centrality_CEN=avg_dgs
df_metrics_biomarkers$n_modules=num_modules
df_metrics_biomarkers$betweenness_centrality_CEN=avg_bet
df_metrics_biomarkers$closeness_centrality_CEN=avg_clo

df_metrics_CEN <- data.frame(CEN_pairwise_distance=avg_dts,CEN_degree_centrality=avg_dgs,CEN_n_modules=num_modules,CEN_betweenness_centrality=avg_bet,CEN_closeness_centrality=avg_clo)
write.table(df_metrics_CEN,file='../results/GEN_metrics.csv', sep=',',row.names=TRUE, col.names=TRUE,quote=FALSE)


