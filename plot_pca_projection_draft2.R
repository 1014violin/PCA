#PCA ANALYSIS
#Alexzandra Morris, alexzandradmorris@gmail.com 
#Nick Bayley wrote some of the existing functions
#last edited July 1, 2021
#this pca analysis can work for just a single dataset or a single dataset + a projected dataset
#this pca script has specific functionality for differential expression results
library("EnhancedVolcano")
library(dplyr)
library(ggplot2)
#https://github.com/kevinblighe/EnhancedVolcano/blob/master/R/EnhancedVolcano.R

#----------set up before running----------
data_dir<-"/Users/alexzandramorris/Desktop/Work/pca_project"
output_dir <- "/Users/alexzandramorris/Desktop/Work/pca_project"
plots_pdf_filename <- "plot_output.pdf"
setwd(data_dir)
orig <- read.table("Nathanson DE PT-XG 20210423.txt", sep = "\t", row.names = 1, header = T)
projected <- read.table("Nathanson DE PTOPCNPC_PT-XG 20210523.txt", sep = "\t", row.names = 1, header = T)
num_filename <- c("_prcomp_scores_orig.txt", "_prcomp_scores_projected.txt", "_prcomp_loadings.txt", "eigenvalues.txt" )
projected_dataset <- F
run_match.rows <- T
run_match.cols <- T
use_signed.log.p <- T #TRY FOR NOW 
colNames = c("log2FoldChange", "signed.log.p")
#colNames = c("log2FoldChange", "pvalue")
MAX.NUM.PCS <- 5
#---------------------------------------FUNCTIONS-------------------------------------
#make sure not to call match.rows() and match.cols() on pca with no projected data 
if(projected_dataset==F) {
  run_match.rows <- F
  run.match.cols <- F
}


match.columns <- function(origg, projectedd){
  x<- select(origg, !matches(colnames(projectedd))) 
  if (dim(x)[2]==0) {
    paste("colnames match in both datasets")
  } else{
    print(dim(x)[2])
    stop("ERROR: COLNAMES DO NOT MATCH IN SETS ")
    #quit(save = "default", status = 0, runLast = TRUE)
  }
    
}

match.rows <- function(data_list){
	names_list <- lapply(data_list, rownames)
	shared_genes <- Reduce(intersect, names_list) #Reduce() reduces a value to a single vector 
	out_list <- lapply(data_list, function(x){
		x[shared_genes,]
		})
	return(out_list)
}


#padj is p adjusted 
add.signed.log.p <- function(de_table, col_name = "padj", sign_name = "log2FoldChange"){ 
	de_table$signed.log.p <- -log10(de_table[,col_name]) * sign(de_table[,sign_name])
	return(de_table)
}

run.pca <- function(data, scale = F){
	pc_res <- prcomp(data, scale = scale) #prcomp() key function behind principal component analysis
	return(pc_res)
}

get.pc.score <- function(pca, full_results = F){ #only used for original data
	if(full_results){
		return(pca$x)
	} else {
		return(pca$x[,1])
	}
}

pca.projection <- function(pca, new_data, scale = F, self_center = F, full_results = F){
	if(self_center){
		projections <-scale(new_data, center = self_center, scale = scale) %*% pca$rotation
	} else {
		projections <-scale(new_data, center = pca$center, scale = scale) %*% pca$rotation
	}
	if(full_results){
		return(projections)
	} else {
		return(projections[,1])
	}
}

#filenames and all_scores are required, The other parameters are optional. 
save_num_results <- function(filenamess, all_scoress, all_projected_scoress, loadingss, evalss){
  setwd(output_dir)
  if ((dim(all_scoress)[2])>MAX.NUM.PCS) all_scoress <- select(all_scoress, 1:5)

  if(projected_dataset!=F){
    if ((dim(all_projected_scoress)[2])>MAX.NUM.PCS) all_projected_scoress <- select(all_projected_scoress, 1:5)
    write.table(all_projected_scoress,filenamess[2],sep='\t',row.names=FALSE,quote=FALSE)
  }
  
  write.table(all_scoress,filenamess[1],sep='\t',row.names=FALSE,quote=FALSE)
  

  write.table(loadingss,filenamess[3],sep='\t',row.names=FALSE,quote=FALSE)
  write.table(evalss, filenamess[4],sep='\t',row.names=FALSE,quote=FALSE)
}

make_plots <- function(de_pcaa, all_scoress, all_projected_scoress, matched_datasetss){
  pdf(file = paste0(output_dir, "/", plots_pdf_filename),   
      width = 8, 
      height = 8) 
  
  #original data: pc1 v. pc2 
  plt1 <- plot(all_scoress[,1], all_scoress[,2],
               main="Original Data: PC1 v. PC2",
               xlab = paste0("PC1  ", summary(de_pcaa)$importance[2,1]*100, "%"),
               ylab=paste0("PC2  ", summary(de_pcaa)$importance[2,2]*100, "%"),
               type="p",
               col="black")
  print(plt1)
  
  #projected data: pc1 v. pc2
  if(projected_dataset==T){
    plt2 <- plot(all_projected_scoress[,1], all_projected_scoress[,2],
                 main="Projected data: PC1 v. PC2",
                 xlab = "PC1",
                 ylab="PC2",
                 type="p",
                 col="black")
    print(plt2)
  }
  
  #pc1 scores v. pc1 scores projected
  if(projected_dataset==T){
    plt3 <- plot(all_scoress[,1], all_projected_scoress[,1],
                 main = "PC1 Scores v. PC1 Projected Scores",
                 xlab = "PC1 Scores",
                 ylab = "PC1 Projected Scores",
                 type = "p",
                 col = "black")
    print(plt3)
  }
  
  gene_up <- rownames(filter(matched_datasetss$ORIG, matched_datasetss$ORIG$log2FoldChange>=0.0))
  gene_down <- rownames(filter(matched_datasetss$ORIG, matched_datasetss$ORIG$log2FoldChange<0.0))
  
  plt4 <- EnhancedVolcano(matched_datasetss$ORIG,
                          lab = rownames(matched_datasetss$ORIG),
                          title = 'Volcano Plot: Original Data',
                          x = 'log2FoldChange',
                          y = 'padj',
                          type = 'p',
                          encircle = gene_down,
                          encircleCol = 'black',
                          encircleSize = 1.5,
                          encircleFill = 'green')
  print(plt4)
  
  if(projected_dataset==T){
    plt5 <- EnhancedVolcano(matched_datasetss$PROJECTED,
                            lab = rownames(matched_datasetss$PROJECTED),
                            title = 'Volcano Plot: Projected Data',
                            x = 'log2FoldChange',
                            y = 'padj',
                            type = 'p')
    print(plt5)
  
     if(use_signed.log.p==T){
       logp<-matched_datasetss$ORIG$signed.log.p
       logp2<-matched_datasetss$PROJECTED$signed.log.p
       plt6 <- ggplot(matched_datasetss$ORIG, aes(x=PC.Score, y=logp)) + geom_point() +
         labs(title="Original data: PC1 vs. signed.log.p",
              x="PC1", y = "signed.log.p")
       
       print(plt6)
       
       plt7<- ggplot(matched_datasetss$PROJECTED, aes(x=PC.Score, y=logp2)) + geom_point() +
         labs(title="Projected data: PC1 vs. signed.log.p",
              x="PC1", y = "signed.log.p")
       
       print(plt7)
       
     }else{
       
    
    p<-matched_datasetss$ORIG$pvalue
    p2<-matched_datasetss$PROJECTED$pvalue
    
    plt6 <- ggplot(matched_datasetss$ORIG, aes(x=PC.Score, y=p)) + geom_point() +
      labs(title="Original data: PC1 vs. pvalue",
           x="PC1", y = "pvalue")

    print(plt6)

    plt7<- ggplot(matched_datasetss$PROJECTED, aes(x=PC.Score, y=p2)) + geom_point() +
      labs(title="Projected data: PC1 vs. pvalue",
           x="PC1", y = "pvalue")

    print(plt7)
     }
   

    plt9<- ggplot(matched_datasetss$PROJECTED, aes(x=PC.Score, y=log2FoldChange)) + geom_point() +
      labs(title="Projected data: PC1 vs. log2FoldChange",
           x="PC1", y = "log2FC")

    print(plt9)

    
  }
  plt8 <- ggplot(matched_datasetss$ORIG, aes(x=PC.Score, y=log2FoldChange)) + geom_point() +
    labs(title="Original data: PC1 vs. log2FoldChange",
         x="PC1", y = "log2FC")
  
  print(plt8)
  
  while (!is.null(dev.list())) {
  dev.off()
  }
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
library(base)
if(projected_dataset==T & run_match.rows==T){
  match.columns(orig, projected)
  matched_datasets <- match.rows(list(ORIG = orig, PROJECTED = projected)) #add parameter match = T
  if(use_signed.log.p==T){
    matched_datasets <- lapply(matched_datasets, add.signed.log.p) 
  }
  de_pca <- run.pca(matched_datasets$ORIG[,colNames], scale = T) 
}
 
if(isTRUE(projected_dataset) && isFALSE(run_match.rows)){
  print("hellooooo")
  match.columns(orig, projected)
  matched_datasets <- list(ORIG = orig, PROJECTED = projected)
  if(use_signed.log.p==T){
    matched_datasets <- lapply(matched_datasets, add.signed.log.p) 
  }
  de_pca <- run.pca(matched_datasets$ORIG[,colNames], scale = T) 
}
 
if(projected_dataset==F & run_match.rows==F){
  matched_datasets <- list(ORIG = orig)
  if(use_signed.log.p==T){
    matched_datasets <- lapply(matched_datasets, add.signed.log.p) 
  }
  de_pca <- run.pca(matched_datasets$ORIG[,colNames], scale = T) 
}

all_scores<-get.pc.score(de_pca, full_results = T) #for all pcs

matched_datasets$ORIG$PC.Score <- all_scores[,1] #just pc1 for this dataset

if(projected_dataset==T){
  all_projected_scores <- pca.projection(de_pca, 
      matched_datasets$PROJECTED[,colNames], scale = de_pca$scale, full_results = T)
  matched_datasets$PROJECTED$PC.Score <- all_projected_scores[,1] #just pc1
}

loadings=de_pca$rotation 
evals <- de_pca$sdev

#SAVE NUMERICAL RESULTS 
if(projected_dataset==T){
  save_num_results(num_filename, all_scores,all_projected_scores,loadings, evals)

}else{
  save_num_results(num_filename, all_scores,NA,loadings, evals)
}

#SAVE PLOTS
if (projected_dataset==T){
make_plots(de_pca, all_scores, all_projected_scores, matched_datasets)
}else{
  make_plots(de_pca, all_scores, NA , matched_datasets)
}

while (!is.null(dev.list())) {
  dev.off()
}

