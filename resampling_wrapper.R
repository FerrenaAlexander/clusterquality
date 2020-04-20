library(Seurat)
library(tidyverse)


setwd('/Users/ferrenaa/Documents/tam/clusterrobustness/kpmodelling_primary')

filelocation = '/Users/ferrenaa/Documents/tam/kpmodelling/sobj_tSNE_GH1306_noceiling.rds'

seq(0.6,1.2,0.1)

rm(list = ls())
filelocation = '/Users/ferrenaa/Documents/tam/kpmodelling/sobj_tSNE_GH1306_noceiling.rds'


for(res in seq(0.6,1.2,0.1) ){
  
  sobj <- readRDS(file = filelocation)
  
  #mkdir for downsampled object storage and outfiles
  dir.create(paste0('ds/', res))
  dir.create(paste0('outs/', res))
  
  #do run clustering with input resolution
  sobj <- FindNeighbors(object = sobj, dims = 1:30)
  sobj <- FindClusters(object = sobj, resolution = res)
  sobj <- RunTSNE(object = sobj, dims = 1:30)
  
  
  #get cells, and count total
  cells <- colnames(sobj)
  samplesize <- length(cells) * 0.1
  
  #set "pool" of cells to pull from and iteratively downsample
  pool <- cells
  
  
  pool <- cells
  for(i in c(1:10)){
    
    
    #randomly pick cells w/o replacement from pool object
    cs <- sample(pool, size = samplesize, replace = F)
    
    ###keep cells NOT selected above --> downstream analysis###
    ds <- cells[cells %in% cs == F]
    
    #remove the selected cells from pool
    #allowing them to be included in ds in next iteration
    pool <- setdiff(pool, cs)
    
    #downsample
    sobjds <- subset(sobj, cells = ds)
    
    
    
    
    #normalization
    sobjds <- NormalizeData(object = sobjds, normalization.method = "LogNormalize", scale.factor = 1e4)
    
    #find DEG
    sobjds <- FindVariableFeatures(object = sobjds, selection.method = 'mean.var.plot',
                                   mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
    
    
    #scale data... warning, takes a long time
    sobjds <- ScaleData(object = sobjds, features = rownames(x = sobjds),
                        vars.to.regress = c("nCount_RNA", "percent.mito"))
    
    #pca
    sobjds <- RunPCA(object = sobjds, pc.genes = sobjds@var.genes)
    
    #graph, cluster, and tSNE
    sobjds <- FindNeighbors(object = sobjds, dims = 1:30)
    sobjds <- FindClusters(object = sobjds, resolution = res)
    sobjds <- RunTSNE(object = sobjds, dims = 1:30)
    
    saveRDS(sobjds, file = paste0("ds/", res, "/dssobj_", as.character(i), ".rds") )
    rm(sobjds)
  }
  
  
  #now find jaccard score
  
  #dir pointer
  dsdir <- paste0('ds/', res)
  
  
  jaccardscorelist <- list()
  jaccardmeans <- data.frame()
  for(ds in list.files(dsdir)){
    
    
    sobjds <- readRDS(paste(dsdir, ds, sep = '/'))
    
    
    
    jaccardscore <- data.frame()
    for(i in as.numeric(levels(sobj@active.ident))){
      
      refclust <- sobj@active.ident[sobj@active.ident==i]
      
      tmp2 <- data.frame()
      for(j in as.numeric(levels(sobj@active.ident))){
        
        dsclust <- sobjds@active.ident[sobjds@active.ident==j]
        
        #jaccard index between reference and downsampled clusters
        tmpscore <- length(intersect( names(refclust), names(dsclust) )) / length(union ( names(refclust), names(dsclust) ))
        tmp2 <- rbind(tmp2, cbind(j, tmpscore) )
        tmp2 <- tmp2[which.max(tmp2$tmpscore),]
        
      }
      
      jaccardscore <- rbind(jaccardscore, tmp2)
    }
    
    colnames(jaccardscore) <- c("ref_prefers", "refmaxjac")
    
    
    #method for counting num clusters that stay the same between ref and downsample
    # after normalizing by jaccard score ; 
    #not intuitive, seems to double count jaccard index, soon to be deprecated
    jaccardscore$same <- seq(0,length(jaccardscore[,1])-1) == jaccardscore[,1]
    jaccardscore$same[jaccardscore$same == T] <- 1
    jaccardscore$same[jaccardscore$same == F] <- -1
    jaccardscore$same <- jaccardscore$same / length(jaccardscore$ref_prefers)
    jaccardscore$same <- jaccardscore$same * jaccardscore$refmaxjac
    
    jaccardscorelist <- c(jaccardscorelist, list(jaccardscore))
    
    jaccardmeans <- rbind(jaccardmeans, mean(jaccardscore$refmaxjac) + sum(jaccardscore$same))
    
    rm(sobjds)
    
  }
  
  colnames(jaccardmeans) <- "jaccardmeans"
  jaccardscorelist <- c(jaccardscorelist, list(jaccardmeans, mean(jaccardmeans[,1])) )
  
  outrds <- paste0('outs/' ,res, '/', res, '.rds')
  outtxt <- paste0('outs/' ,res, '/', res, '.txt')

  saveRDS(jaccardscorelist ,outrds)
  
  sink(outtxt)
  print(jaccardscorelist)
  sink()
  
  
}





#read out
scores <- data.frame()
for(i in list.files('outs/')){
  tmp <- readRDS(paste0('outs/', i, '/', i, '.rds'))
  
  sd <- sd(tmp[[length(tmp)-1]][,1])
  se <- sd / sqrt(tmp[[length(tmp)-1]][,1])
  ci <- se * qt(.95/2 + 0.5, length(tmp[[length(tmp)-1]][,1]) - 1)
  
  scores <- rbind(scores, as.numeric( cbind(i, tmp[[length(tmp)]], ci)) )
  
  
  
}

colnames(scores) <- c('res', 'score', 'ci')

ggplot(scores, aes(x = res, y = score))+
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymin=score-ci, ymax=score+ci), colour="black", width=.03)+
  scale_x_continuous(breaks = scores[,1])+
  xlab('Resolution')+
  theme_minimal()





#read out; ignore weird scoring method from before
scores <- data.frame()
for(i in list.files('outs/')){
  tmp <- readRDS(paste0('outs/', i, '/', i, '.rds'))
  
  tmp[[length(tmp)-1]] <- data.frame()
  for(j in 1:length(list.files(paste0('ds/', i)))){
    
    tmp[[length(tmp)-1]] <- rbind(tmp[[length(tmp)-1]], 
                                  mean(tmp[[j]][,2]))
    
    
  }
  
  tmp[[length(tmp)]] <- mean(tmp[[length(tmp)-1]][,1]) + 1/j
  
  sd <- sd(tmp[[length(tmp)-1]][,1])
  se <- sd / sqrt(length(tmp[[length(tmp)-1]][,1]))
  ci <- se * qt(.95/2 + 0.5, length(tmp[[length(tmp)-1]][,1]) - 1)
  
  scores <- rbind(scores, as.numeric( cbind(i, tmp[[length(tmp)]], ci)) )
  
  
  
}

colnames(scores) <- c('res', 'score', 'ci')

ggplot(scores, aes(x = res, y = score, color = '5-fold'))+
  geom_line()+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=score-ci, ymax=score+ci), width=.03)+
  scale_x_continuous(breaks = scores[,1])+
  xlab('Resolution')+
  ylab('Jaccard Index Score')+
  theme_minimal()+
  geom_line(data=scores10, aes(color="10-fold"))+
  geom_point(data=scores10, aes(color="10-fold"), size=2)+
  geom_errorbar(data=scores10, aes(ymin=score-ci, ymax=score+ci, color='10-fold'), 
                width=.05, alpha=0.75, size=1)



pdf('outtsne.pdf')

for(dsdir in list.files('ds') ){
  dsdir <- paste0('ds/', dsdir, '/')
  

#dsdir <- "ds/0.1/"
plotlist <- list()
lablist <- c()
for(ds in list.files(dsdir)){
  sobjds <- readRDS(paste0(dsdir, ds))
  lablist <- c(lablist, paste(str_split(dsdir, pattern = '/')[[1]][2],
               str_split(ds, pattern = '\\.')[[1]][1], sep = '_' )
               )
  plotlist <- c(plotlist, list(DimPlot(sobjds, label = T, label.size = 5)))
}

print ( cowplot::plot_grid(plotlist = plotlist, labels = lablist) ) 


}
dev.off()


