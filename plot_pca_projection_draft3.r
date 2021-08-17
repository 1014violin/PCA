#PCA ANALYSIS
#Alexzandra Morris, alexzandradmorris@gmail.com 
#Nick Bayley wrote some of the existing functions
#last edited July 20, 2021
#this pca analysis can work for just a single dataset or a single dataset + a projected dataset
#this pca script has specific functionality for differential expression results
#________________________________

library("EnhancedVolcano")
library(dplyr)
library(ggplot2)
library(base)

#----------set up before running----------
working_dir<-"/Users/alexzandramorris/Desktop/Work/pca_project"
setwd(working_dir)
filename_prefix = "GBM_proj"
orig <- read.table("Nathanson DE PT-XG 20210423.txt", sep = "\t", row.names = 1, header = T)
projected <- read.table("Nathanson DE PTOPCNPC_PT-XG 20210523.txt", sep = "\t", row.names = 1, header = T)
projected_dataset <- T
run_match.rows <- T
run_match.cols <- T

colNames = c("log2FoldChange", "signed.log.p")#used for run_pca(). 

MAX.NUM.PCS <- 100
PLOTS_PDF <- paste0(filename_prefix, "_plot_output.pdf")
FILENAMES <- c(paste0(filename_prefix,"_prcomp_scores_orig.txt"), paste0(filename_prefix,"_prcomp_scores_projected.txt"), paste0(filename_prefix,"_prcomp_loadings.txt"), paste0(filename_prefix, "_eigenvalues.txt" ))



#---------------------------------------FUNCTIONS-------------------------------------
#make sure not to call match.rows() and match.cols() on pca with no projected data 
if(projected_dataset==F) {
  run_match.rows <- F
  run.match.cols <- F
}

#Returns void 
# @param {data.frame} origg The original data set. 
# @param {data.frame} projectedd The projected data set. 
# @return {void}
match.columns <- function(origg, projectedd){
  x<- select(origg, !matches(colnames(projectedd))) 
  if (dim(x)[2]==0) {
    paste("colnames match in both datasets")
  } else{
    print(dim(x)[2])
    #stop("ERROR: COLNAMES DO NOT MATCH IN SETS ")
    quit(save = "default", status = 0, runLast = TRUE)
  }
    
}

#Returns out_list
# @param {list} data_list, A list of the two lists: orig, projected.
# @return {list} out_list, A list of two lists, now with the same number of rows.
match.rows <- function(data_list){
	names_list <- lapply(data_list, rownames)
	shared_genes <- Reduce(intersect, names_list) 
	out_list <- lapply(data_list, function(x){
		x[shared_genes,]
		})
	return(out_list)
}

#Returns de_table
# @param {list} de_table, List of the two lists: Orig and Projected
# @param {String} col_name, Column name that contains the pvalue or p adjusted value we are using 
# @param {String} sign_name,   "    "     "     "      " log2FoldChange values 
# @return {list} de_table, list of the two lists
add.signed.log.p <- function(de_table, col_name = "padj", sign_name = "log2FoldChange"){ 
	de_table$signed.log.p <- -log10(de_table[,col_name]) * sign(de_table[,sign_name])
	return(de_table)
}

#Returns pc_res
# @param {data.frame list} data, datasets$ORIG[,colNames]
# @param {Boolean} scale
# @return {prcomp list} pc_res, list of the two lists
run.pca <- function(data, scale = F){
	pc_res <- prcomp(data, scale = scale) 
	return(pc_res)
}
#@function get.pc.score() 
# @param {prcomp}{list} de_pca, List of 5 components: sdev, rotation, center, scale, x, loadings
# @param {Boolean} full_results, Option to return all the pc scores or pc1 
# @return {double}{"matrix" "array"} pca$x or pca$x[,1]
get.pc.score <- function(pca, full_results = F){ #only used for original data
	if(full_results){
		return(pca$x)
	} else {
		return(pca$x[,1])
	}
}

#@function pca.projection() 
# @param {prcomp}{list} pca
# @param {data.frame}{list} new_data, new_data is datasets$PROJECTED[,colNames]
# @param {String} scale
# @param {list} self_center
# @param {list} full_results
# @return {"matrix" "array"}{list} projections or projections[,1]
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
#@function save_num_results() 
# @param {"matrix" "array"} all_scoress
# @param {"matrix" "array"} all_projected_scoress
# @param {"matrix" "array" } loadingss
# @param {numeric} evalss
# @return {void} 
#all_scoress is required, The other parameters are optional. 
save_num_results <- function(all_scoress, all_projected_scoress, loadingss, evalss){
  setwd(working_dir)
  if ((dim(all_scoress)[2])>MAX.NUM.PCS) all_scoress <- select(all_scoress, 1:5)

  if(projected_dataset==T){
    if ((dim(all_projected_scoress)[2])>MAX.NUM.PCS) all_projected_scoress <- select(all_projected_scoress, 1:5)
    write.table(all_projected_scoress,FILENAMES[2],sep='\t',row.names=FALSE,quote=FALSE)
  }
  
  write.table(all_scoress,FILENAMES[1],sep='\t',row.names=FALSE,quote=FALSE)
  write.table(loadingss,FILENAMES[3],sep='\t',row.names=FALSE,quote=FALSE)
  write.table(evalss, FILENAMES[4],sep='\t',row.names=FALSE,quote=FALSE)
}

#@function make_plots() 
# @param {prcomp list} de_pcaa
# @param {"matrix" "array"} all_scoress
# @param {"matrix" "array" } all_projected_scoress
# @param {list} datasetss List of two lists Original set and Projected set.
# @return {void}
make_plots <- function(de_pcaa, all_scoress, all_projected_scoress, datasetss){
  pdf(file = paste0(working_dir, "/", PLOTS_PDF), onefile = TRUE,   
      width = 8, 
      height = 8) 
  
  #NO PROJECTION PLOTS 
  orig1 <- plot(all_scoress[,1], all_scoress[,2],
               main="Original Data: PC1 v. PC2",
               xlab = paste0("PC1  ", summary(de_pcaa)$importance[2,1]*100, "%"),
               ylab=paste0("PC2  ", summary(de_pcaa)$importance[2,2]*100, "%"),
               type="p",
               col="black")
 
  gene_up <- rownames(filter(datasetss$ORIG, datasetss$ORIG$log2FoldChange>=0.05))
  gene_down <- rownames(filter(datasetss$ORIG, datasetss$ORIG$log2FoldChange<0.05))
  
  orig2 <- EnhancedVolcano(datasetss$ORIG,
                          lab = rownames(datasetss$ORIG),
                          title = 'Volcano Plot: Original Data',
                          x = 'log2FoldChange',
                          y = 'padj',
                          type = 'p',
                          encircle = gene_down,
                          encircleCol = 'black',
                          encircleSize = .5,
                          encircleFill = 'pink')
  
  orig3 <- ggplot(datasetss$ORIG, aes(x=PC.Score, y=log2FoldChange)) + geom_point() +
    labs(title="Original data: PC1 vs. log2FoldChange",
         x="PC1", y = "log2FC")
  
  print(orig1)
  print(orig2)
  print(orig3)

  #################################################################################

  #PROJECTED DATA PLOTS 
  if(projected_dataset==T & run_match.rows==T){
    proj1 <- plot(all_projected_scoress[,1], all_projected_scoress[,2],
                 main="Projected data: PC1 v. PC2",
                 xlab = "PC1",
                 ylab="PC2",
                 type="p",
                 col="black")
    
    #pc1 scores v. pc1 scores projected
    proj2 <- plot(datasetss$ORIG$PC.Score, datasetss$PROJECTED.PC.Score,
                 main = "PC1 Scores v. PC1 Projected Scores",
                 xlab = "PC1 Scores",
                 ylab = "PC1 Projected Scores",
                 type = "p",
                 col = "black")

    proj3 <- EnhancedVolcano(datasetss$PROJECTED,
                            lab = rownames(datasetss$PROJECTED),
                            title = 'Volcano Plot: Projected Data',
                            x = 'log2FoldChange',
                            y = 'padj',
                            type = 'p')
    
    logp<-datasetss$ORIG$signed.log.p
    logp2<-datasetss$PROJECTED$signed.log.p
       
    proj4 <- ggplot(datasetss$ORIG, aes(x=PC.Score, y=logp)) + geom_point() +
         labs(title="Original data: PC1 vs. signed.log.p",
              x="PC1", y = "signed.log.p")

    proj5<- ggplot(datasetss$PROJECTED, aes(x=PC.Score, y=logp2)) + geom_point() +
         labs(title="Projected data: PC1 vs. signed.log.p",
              x="PC1", y = "signed.log.p")
       
    proj6<- ggplot(datasetss$PROJECTED, aes(x=PC.Score, y=log2FoldChange)) + geom_point() +
      labs(title="Projected data: PC1 vs. log2FoldChange",
           x="PC1", y = "log2FC")

    print(proj1)
    print(proj2)
    print(proj3)
    print(proj4)
    print(proj5)
    print(proj6)
    
  }
  
  
  while (!is.null(dev.list())) {
  dev.off()
  }
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


##############################
if(run_match.cols==T) match.columns(orig, projected)

if(projected_dataset==T){
  datasets <- lapply(list(ORIG = orig, PROJECTED = projected), add.signed.log.p) 
  
}else{
  datasets <- lapply(list(ORIG = orig), add.signed.log.p) 
}

de_pca <- run.pca(datasets$ORIG[,colNames], scale = T)
all_scores <- get.pc.score(de_pca, full_results = T) 
datasets$ORIG$PC.Score <- all_scores[,1]
loadings=de_pca$rotation 
evals <- de_pca$sdev

#SAVE NUMERICAL RESULTS + SAVE PLOTS
if(projected_dataset==T){
  all_projected_scores <- pca.projection(de_pca,datasets$PROJECTED[,colNames], scale = de_pca$scale, full_results = T)
  datasets$PROJECTED$PC.Score <- all_projected_scores[,1] #just pc1
  if(run_match.rows==T) datasets <- match.rows(datasets)
  
  save_num_results(all_scores,all_projected_scores,loadings, evals)
  make_plots(de_pca, all_scores, all_projected_scores, datasets)
  
}else{
  save_num_results(all_scores,NA,loadings, evals)
  make_plots(de_pca, all_scores, NA , datasets)
  
}

#ignore warnings#

while (!is.null(dev.list())) {
  print("dev.off()")
  dev.off()
}

