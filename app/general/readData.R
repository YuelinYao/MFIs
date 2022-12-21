path_countmatrix="./data/Count_matrix.csv"
path_metadata="./data/Meta_Data.csv"


#read_data<-function(path_countmatrix,path_metadata){
  
## Count matrix
Count_matrix <- as.matrix(fread(path_countmatrix),rownames=1)
Meta_data<-read.csv(path_metadata,row.names = 1,colClasses = "character")

## Binary matrix
count<-Count_matrix
count[count>0]<-1
count[count==0]<-0

# Creat CreateSeuratObject
Count_matrix<-t(Count_matrix) 
colnames(Meta_data)<-c("Cell_State","Cell_Types")

srt<-CreateSeuratObject(counts =Count_matrix,meta.data =Meta_data)
## preprocessing srt object
srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 10000)
srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)
srt <- ScaleData(srt)

#srt <- SCTransform(srt, return.only.var.genes = FALSE)

#data<-list(count,Meta_data,srt)
#return(data)

#}
#srt = SCTransform(srt, return.only.var.genes = FALSE)


#save(Meta_data,count,srt,file="./data/Example.Rdata")



