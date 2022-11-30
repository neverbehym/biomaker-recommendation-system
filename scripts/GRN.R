library(igraph)
library(ggplot2)
library(viridis)
library(centiserve)
library(expm)

args = commandArgs(trailingOnly=TRUE)
Panellog_FILE<-args[1]
GeneRegulatoryNetwork_FILE<-args[2]

############# read a pool of candidate biomarker panels #############

panel_log <- read.table('../results/CandidateBiomarkerList.csv', sep=",", header = TRUE,stringsAsFactors = FALSE)

panel_list <- strsplit(panel_log$genes,';')

biomarker_union<-unique(unlist(panel_list))
length(biomarker_union)

############# read Gene Regulatory Network #############

regnet_regulator<-read.table("../data/GeneRegulatoryNetwork.csv",quote="", sep="\t",header = TRUE, stringsAsFactors = FALSE)
colnames(regnet_regulator)<-c('target','source','sign')
node_list<-unique(c(regnet_regulator$source,regnet_regulator$target))
edge_table<-data.frame(source_id=match(regnet_regulator$source,node_list)-1,target_id=match(regnet_regulator$target,node_list)-1)
edge_table<-cbind(regnet_regulator,edge_table)

regulon_lists<-list()
for (regulator in unique(regnet_regulator$source))
  regulon_lists[[regulator]]<-regnet_regulator$target[regnet_regulator$source==regulator]
lengths(regulon_lists)


############# create the graph of regulatory network #############

edge_list<-as.vector(t(as.matrix(regnet_regulator[,1:2])))
length(edge_list)
g_regnet <- make_graph(edge_list,directed = FALSE)
vertex_list<-V(g_regnet)$name
length(vertex_list)
set_edge_attr(g_regnet,name='sign',index=E(g_regnet),value=regnet_regulator[,3])

############# obtain the network metrics #############

### degree
dg<-degree(g_regnet)

### distances of all pairwise nodes
distance_table<-distances(g_regnet)

### connected components
comps<-components(g_regnet)

### betweenness centrality
betweenness_cent<-estimate_betweenness(g_regnet, cutoff=max(distance_table),directed = FALSE)

### latora closeness centrality
latoraclose_cent<-closeness.latora(g_regnet)

###  eigenvector centrality
eigen_cent<-eigen_centrality(g_regnet, directed = FALSE, scale = TRUE)
eigen_cent<-eigen_cent$vector

### topological coefficients
topocoe_cent<-topocoefficient(g_regnet)
topocoe_cent[is.nan(topocoe_cent)]<-0


df_vp_GRN<-data.frame(component_id=as.character(comps$membership),component_type=c('Other components','Largest component')[(comps$membership==1)+1],
                            betweenness=betweenness_cent,closeness=latoraclose_cent,eigenvector=eigen_cent,topological=topocoe_cent,degree=dg)


############# compute the evaluation metrics for biomarker panels #############

num_biomarkers <- lengths(panel_list)

num_regulators=vector()
num_biomarkersGRNIncluded=vector()
num_components=vector()
avg_dgs=vector() # average node degrees of biomarkers
avg_indgs=vector() # average node in degrees of biomarkers
avg_outdgs=vector() # average node out degrees of biomarkers
avg_bet=vector()  # average betweenness centrality of biomarkers
avg_clo=vector() # average closeness centrality of biomarkers
avg_dts=vector() # average pairwise distance between biomarkers

connectedregulators_list=list()
for (i in 1:length(panel_list)){
  panel<-panel_list[[i]]
  
  panel_slim<-panel[panel %in% c(regnet_regulator$source,regnet_regulator$target)]
  num_biomarkersGRNIncluded<-c(num_biomarkersGRNIncluded,length(panel_slim))
  
  connected_regulators<-regnet_regulator$source[regnet_regulator$target%in%panel_slim]
  connectedregulators_list=append(connectedregulators_list,list(connected_regulators))
  
  num_regulators<-c(num_regulators,length(connected_regulators))
  num_components<-c(num_components,length(unique(comps$membership[panel_slim])))
  
  avg_dt<- sum(distance_table[panel_slim,panel_slim])/length(panel_slim)/(length(panel_slim)-1)
  avg_dts<-c(avg_dts,avg_dt)

  avg_dgs<-c(avg_dgs,mean(df_vp_GRN[panel_slim,'degree']))

  avg_bet<-c(avg_bet,mean(df_vp_GRN[panel_slim,'betweenness']))
  avg_clo<-c(avg_clo,mean(df_vp_GRN[panel_slim,'closeness']))
  
}

avg_dts[is.na(avg_dts)] = max(avg_dts,na.rm = TRUE)+1

df_metrics_GRN <- data.frame(GRN_n_BiomarkerIncluded=num_biomarkersGRNIncluded,GRN_n_regulators=num_regulators,GRN_n_components=num_components,GRN_Degree=avg_dgs,GRN_pairwise_distance=avg_dts,GRN_betweenness_centrality=avg_bet,closeness_centrality=avg_clo)
write.table(df_metrics_GRN,file='../results/GRN_metrics.csv', sep=',',row.names=TRUE, col.names=TRUE,quote=FALSE)
