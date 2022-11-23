
############# read a pool of candidate biomarker panels #############

panel_log <- read.table('../data/CandidateBiomarkerList.csv', sep=",", header = TRUE,stringsAsFactors = FALSE)

panel_list <- strsplit(panel_log$genes,';')

biomarker_union<-unique(unlist(panel_list))
length(biomarker_union)

### read gene expression data and multi-class (cellular states) labels ###
m_express <- as.matrix(read.table("../data/NicolasGeneExpression_processed.csv", sep=",", header = TRUE, stringsAsFactors = FALSE,row.names = 1))
multiclass_labels <- read.table('../data/MultiClassLabel.csv',sep=',',stringsAsFactors = FALSE,row.names = 1)
num_classes <- length(unique(multiclass_labels$classID))

############# 1) performance score #############

panel_log$performance

############# 2) frequency score #############

s_freqPerbio<- rep(0,length(biomarker_union))
names(s_freqPerbio)<-biomarker_union
for (i in 1:length(panel_list))
  s_freqPerbio[panel_list[[i]]]<- s_freqPerbio[panel_list[[i]]]+1
s_freqPerbio<-s_freqPerbio/length(panel_list)

s_freqPerPanel<-rep(0,length(panel_list))
names(s_freqPerPanel)<-1:length(panel_list)
for (i in 1:length(panel_list))
  s_freqPerPanel[i]<-mean(s_freqPerbio[panel_list[[i]]])

############ 3) effect size score #######

cohen_d<- matrix(0,nrow=length(biomarker_union),ncol = 10)
for (i in 1:num_classes){
  exp_groupA<-m_express[biomarker_union,which(class==i)]
  min_d<-rep(1000,length(biomarker_union))
  for (j in c(1:(i-1),(i+1):10)){
    if((i==1 & j<2) | (i==10 & j>9))
      next
    exp_groupB<-m_express[biomarker_union,which(class==j)]
    n1<-dim(exp_groupA)[2]
    n2<-dim(exp_groupB)[2]
    pooled_sd<-sqrt(((matrixStats::rowSds(exp_groupA))^2*(n1-1)+(matrixStats::rowSds(exp_groupB))^2*(n2-1))/(n1+n2-2))
    temp_d<-(rowMeans(exp_groupA)-rowMeans(exp_groupB))/pooled_sd

    min_d<-ifelse(abs(temp_d)<abs(min_d),temp_d,min_d)
  }
  cohen_d[,i]<-min_d
}
rownames(cohen_d)<-biomarker_union
colnames(cohen_d)<-paste0('class',1:10)

s_dscorePerPanel<-rep(0,length(panel_list))
names(s_dscorePerPanel)<-1:length(panel_list)
for (i in 1:length(panel_list)){
  panel<-panel_list[[i]]
  s_dscorePerPanel[i]<-mean(matrixStats::colMaxs(abs(cohen_d[panel,])))
}

############# write the evaluation metrics in csv file #############
df_metrics_BIM<-data.frame(perfermance_f1score=panel_log$performance,frequency_score=s_freqPerPanel[as.character(1:length(panel_list))],effectSize_score=s_dscorePerPanel[as.character(1:length(panel_list))])
df_metrics_BIM <- data.frame(size=panel_log$size, perfermance_f1score=panel_log$performance,frequency_score=s_freqPerPanel,effectSize_score=s_dscorePerPanel)
