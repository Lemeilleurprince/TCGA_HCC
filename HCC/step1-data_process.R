## 
### ---------------
###
### Create: Noe Lee
### Email: 2021585592@qq.com
###
### ---------------
options(stringsAsFactors = F)
path = './HCC/mainfest'
floders = list.files(path)#列出总目录下含有文件夹的名称
BRCA_counts = data.frame()#构建dataframe
fd1 = floders[1]#第一个文件；
file_name  =list.files(paste(path,'/',fd1,sep =''))#列出文件夹fd1中的全部文件名
file_list = substr(file_name[1],1,28)#对第一个文件名进行文字截图；
mydata = read.table(gzfile(paste(path,'/',fd1,'/',file_name[1],sep = '')))#读取解压包gz内部的数据；
names(mydata) <- c('ECSG_ID',file_list)#读取第一个gz文件之后，把文件名当作列名；
BRCA_counts = mydata#赋值为BRCA——counts矩阵；
for (fd in floders[2:179]){#循环对179个文件进行处理；
  files_name = list.files(paste(path,'/',fd,sep = ''))
  print(files_name[1])
  file_list =substr(files_name[1],1,28)
  mydata = read.table(gzfile(paste(path,'/',fd,'/',files_name[1],sep = '')))
  names(mydata) <- c('ECSG_ID',file_list)
  BRCA_counts <- merge(BRCA_counts,mydata,by ='ECSG_ID')#进行逐个以ensg编号进行合并；
  
}
write.csv(BRCA_counts,'/home/wenlong/TCGA/CA/BRCA_counts.csv')#把数据写入csv文件中

setwd("./HCC")
library(rjson)
exprSet<-read.csv('HCC_counts.csv')
pdata<-fromJSON(file='metadata.cart.2020-03-20.json')
library(biomaRt)
library(curl)
rownames(exprSet) <- exprSet[,2]
exprSet <- exprSet[c(-1,-2)]
print(rownames(exprSet))
char =substr(rownames(exprSet),1,15)
rownames(exprSet) <- substr(rownames(exprSet),1,15)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))  #构建mart
my_ensembl_gene_id<-char #需要转换的exsembl的编码；
my_ensembl_gene_id
mms_symbols<- getBM(attributes=c('ensembl_gene_id','hgnc_symbol',"description"),filters = 'ensembl_gene_id',values = my_ensembl_gene_id,mart = mart)#基因注释之后
hgnc_symbol <- mms_symbols[,2]
result_diff<-cbind(hgnc_symbol,exprSet[1:56543,])
write.csv(result_diff,'result.csv')
exprSet<-read.csv('result.csv')
exprSet<-t(exprSet)
exprSet[1:4,1:4]
colnames(exprSet)<-exprSet[2,]
exprSet<-t(exprSet[c(-1,-2),])
head(exprSet)
save(exprSet,file='exprSet.Rdata')

load('exprSet.Rdata')

##再想想如何过滤掉重复基因?
id=colnames(exprSet)
code =substr(id,14,15)
code=as.numeric(code)
group=ifelse(code %in% 1:9,"Tumor","Normal")
group_list=cbind(id,group)
exprSet<-apply(exprSet,1,as.numeric)
rownames(exprSet)=group_list[,1]
exprSet=t(exprSet)
##筛除表达量低的基因
exprSet=exprSet[apply(exprSet,1, function(x) sum(x>1) > 5),]
#查看不重复基因数
table(!duplicated(rownames(exprSet)))
{
  MAX = by( exprSet,rownames(exprSet) , 
            function(x) rownames(x)[ which.max( rowMeans(x) ) ] )
  MAX = as.character(MAX)
  exprSet = exprSet[ rownames(exprSet) %in% MAX , ]
}
exprSet[1:4,1:4]
exprSet<-log(exprSet+0.001)
save(exprSet,group_list,file='step1-output.Rdata')

## 聚类
library( "ggfortify" )
load('step1-output.Rdata')
exprSet=t(exprSet)
{
  nodePar <- list( lab.cex = 0.3, pch = c( NA, 19 ), cex = 0.3, col = "red" )
  hc = hclust( dist( exprSet ) )
  png('hclust.png')
  plot( as.dendrogram( hc ), nodePar = nodePar, horiz = TRUE )
  dev.off()
}

##主成分分析
dat=exprSet#画PCA图时要求是行名时样本名，列名时探针名
dat=as.data.frame(dat)#将matrix转换为data.frame
dat=cbind(dat,group_list[,2]) #cbind横向追加，即将分组信息追加到最后一列
library("FactoMineR")#画主成分分析图需要加载这两个包
library("factoextra") 

dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)#现在dat最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$group_list, # color by groups
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
ggsave('PCA.png')  

##limma差异基因分析
library( "limma" )
{
  design <- model.matrix( ~0 + factor( group_list[,2]) )
  colnames( design ) = levels( factor( group_list[,2] ) )
  rownames( design ) = group_list[,1]
}
  design
  contrast.matrix <- makeContrasts( "Tumor-Normal", levels = design )
  contrast.matrix
{
  fit <- lmFit( t(exprSet), design )
  fit2 <- contrasts.fit( fit, contrast.matrix ) 
  fit2 <- eBayes( fit2 )
  nrDEG = topTable( fit2, coef = 1, n = Inf )
  write.table( nrDEG, file = "nrDEG")
}
head(nrDEG)

##绘制热图
library( "pheatmap" )
exprSet=t(exprSet)
{
  #差异结果需要先根据p值挑选
  tmp = nrDEG[nrDEG$P.Value < 0.05,]
  
  nrDEG_D = tmp[ order( tmp$logFC ), ]
  nrDEG_U = tmp[ order( -tmp$logFC ), ]
  choose_gene = c( rownames( nrDEG_D )[1:100], rownames( nrDEG_U )[1:100] )
  choose_matrix = exprSet[ choose_gene, ]
  choose_matrix = t( scale( t( choose_matrix ) ) )
  
  choose_matrix[choose_matrix > 2] = 2
  choose_matrix[choose_matrix < -2] = -2
  
  annotation_col = data.frame(Participants = factor( group_list[,2] ) )
  rownames( annotation_col ) = colnames( exprSet )
  pheatmap( fontsize = 2, choose_matrix, annotation_col = annotation_col, show_rownames = T, annotation_legend = T, filename = "heatmap.png")
}

##绘制火山图

library( "ggplot2" )
logFC_cutoff <- with( nrDEG, mean( abs( logFC ) ) + 2 * sd( abs( logFC ) ) )
logFC_cutoff
{
  nrDEG$change = as.factor( ifelse( nrDEG$P.Value < 0.05 & abs(nrDEG$logFC) > logFC_cutoff,
                                    
                                    ifelse( nrDEG$logFC > logFC_cutoff , 'UP', 'DOWN' ), 'STABLE' ) )
  
  save( nrDEG, file = "nrDEG_array.Rdata" )
  
  Title <- paste0( 'Cutoff for logFC is ', round( logFC_cutoff, 3 ),
                       
                       '\n The number of up gene is ', nrow(nrDEG[ nrDEG$change =='UP', ] ),
                       
                       '\nThe number of down gene is ', nrow(nrDEG[ nrDEG$change =='DOWN', ] ) )
  
  volcano = ggplot(data = nrDEG, aes( x = logFC, y = -log10(P.Value), color = change)) +
    geom_point( alpha = 0.4, size = 1.75) +
    theme_set( theme_set( theme_bw( base_size = 15 ) ) ) +
    xlab( "log2 fold change" ) + ylab( "-log10 p-value" ) +
    ggtitle( Title ) + theme( plot.title = element_text( size = 15, hjust = 0.5)) +
    scale_colour_manual( values = c('blue','black','red') )
  ggsave( volcano, filename = 'volcano.png' )
  dev.off()
}


##注释
library( "clusterProfiler" )
library( "org.Hs.eg.db" )

df <- bitr( rownames( nrDEG ), fromType = "SYMBOL", toType = c( "ENTREZID" ), OrgDb = org.Hs.eg.db )
head( df )

{
  nrDEG$SYMBOL = rownames( nrDEG )
  nrDEG = merge( nrDEG, df, by='SYMBOL' )
}

head( nrDEG )

{
  gene_up = nrDEG[ nrDEG$change == 'UP', 'ENTREZID' ] 
  gene_down = nrDEG[ nrDEG$change == 'DOWN', 'ENTREZID' ]
  gene_diff = c( gene_up, gene_down )
  gene_all = as.character(nrDEG[ ,'ENTREZID'] )
}

{
  geneList = nrDEG$logFC
  names( geneList ) = nrDEG$ENTREZID
  geneList = sort( geneList, decreasing = T )
}

library( "ggplot2" )

# kegg  enrich 
{
  {
    ## KEGG pathway analysis
    kk.up <- enrichKEGG(   gene          =  gene_up    ,
                           
                           organism      =  'hsa'      ,
                           
                           universe      =  gene_all   ,
                           
                           pvalueCutoff  =  0.99       ,
                           
                           qvalueCutoff  =  0.99        )
    
    kk.down <- enrichKEGG( gene          =  gene_down  ,
                           
                           organism      =  'hsa'      ,
                           
                           universe      =  gene_all   ,
                           
                           pvalueCutoff  =  0.99       ,
                           
                           qvalueCutoff  =  0.99        )
    
  }
  
  head( kk.up )[ ,1:6 ]
  head( kk.down )[ ,1:6 ]
  
  kegg_down_dt <- as.data.frame( kk.down )
  kegg_up_dt <- as.data.frame( kk.up )
  
  down_kegg <- kegg_down_dt[ kegg_down_dt$pvalue < 0.05, ]
  down_kegg$group = -1
  up_kegg <- kegg_up_dt[ kegg_up_dt$pvalue < 0.05, ]
  up_kegg$group = 1

  dat = rbind( up_kegg, down_kegg )
  dat$pvalue = -log10( dat$pvalue )
  dat$pvalue = dat$pvalue * dat$group
  dat = dat[ order( dat$pvalue, decreasing = F ), ]
  g_kegg <- ggplot( dat, 
                    aes(x = reorder( Description, order( pvalue, decreasing=F ) ), y = pvalue, fill = group)) + 
    geom_bar( stat = "identity" ) + 
    scale_fill_gradient( low = "blue", high = "red", guide = FALSE ) + 
    scale_x_discrete( name = "Pathway names" ) +
    scale_y_continuous( name = "log10P-value" ) +
    coord_flip() + theme_bw() + theme( plot.title = element_text( hjust = 0.5 ) ) +
    ggtitle( "Pathway Enrichment" ) 
  print( g_kegg )
  ggsave( g_kegg, filename = 'kegg_up_down.png' )
}

{
  ###  GSEA 
  kk_gse <- gseKEGG(geneList     = geneList,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 30,
                    pvalueCutoff = 0.95,
                    verbose      = FALSE)
  head(kk_gse)[,1:6]
  gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
  down_kegg<-kk_gse[kk_gse$pvalue<0.01 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
  up_kegg<-kk_gse[kk_gse$pvalue<0.01 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
  
  dat = rbind( up_kegg, down_kegg )
  dat$pvalue = -log10( dat$pvalue )
  dat$pvalue = dat$pvalue * dat$group
  dat = dat[ order( dat$pvalue, decreasing = F ), ]
  
  g_kegg <- ggplot( dat, 
                    aes(x = reorder( Description, order( pvalue, decreasing=F ) ), y = pvalue, fill = group)) + 
    geom_bar( stat = "identity" ) + 
    scale_fill_gradient( low = "blue", high = "red", guide = FALSE ) + 
    scale_x_discrete( name = "Pathway names" ) +
    scale_y_continuous( name = "log10P-value" ) +
    coord_flip() + theme_bw() + theme( plot.title = element_text( hjust = 0.5 ) ) +
    ggtitle( "Pathway Enrichment" ) 
  print( g_kegg )
  ggsave(g_kegg,filename = 'kegg_up_down_gsea.png')
}

# GO annotation
g_list = list( gene_up = gene_up, gene_down = gene_down, gene_diff = gene_diff)
go_enrich_results <- lapply( g_list, function( gene ) {
  lapply( c( 'BP', 'MF', 'CC' ) , function( ont ) {
    cat( paste( 'Now process', ont ) )
    ego <- enrichGO( gene          =  gene,
                     universe      =  gene_all,
                     OrgDb         =  org.Hs.eg.db,
                     ont           =  ont ,
                     pAdjustMethod =  "BH",
                     pvalueCutoff  =  0.99,
                     qvalueCutoff  =  0.99,
                     readable      =  TRUE)
    print( head( ego ) )
    return( ego )
  })
  
})

save( go_enrich_results, file = 'go_enrich_results.Rdata' )
n1 = c( 'gene_up', 'gene_down', 'gene_diff' )
n2 = c( 'BP', 'MF', 'CC' ) 
for ( i in 1:3 ){
  for ( j in 1:3 ){
    fn = paste0( 'dotplot_', n1[i], '_', n2[j], '.png' )
    cat( paste0( fn, ' ' ) )
    png( fn, res = 150, width = 1080 )
    print( dotplot( go_enrich_results[[i]][[j]] ) )
    dev.off()
  }
}


