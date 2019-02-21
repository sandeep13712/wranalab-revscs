setwd("F:/Computational/SingleCellAnalysis/SetTheory/ForCluPaper/C07_CEN/")
data = read.csv("C07_raw_data.csv", header = TRUE,  row.names = 1, sep=",")


library(igraph)

#Step 2: filetering the geneList
{
  {
    listOfImportantGenes = read.csv("Differential_Markers_Crypts_Cluster_Analysis_OnlyStemCells_Final.csv", header = TRUE, sep=",")
    
    geneClassesStrength = rep(0, length(listOfImportantGenes[1,]))
    
    AllMarkers = NULL
    for (i in 1:length(listOfImportantGenes[1,]))
    {
      print(i)
      geneClassesStrength[i] = length(unique(listOfImportantGenes[,i]))#-1
      AllMarkers = union(AllMarkers, listOfImportantGenes[1:geneClassesStrength[i],i])
    }
    
    # Delete all genes from AllMarkers lists which are not present in data
    indexToDelete = NULL
    for (i in 1:length(AllMarkers))
    {
      
      indexNo = which(tolower(rownames(data)) == tolower(AllMarkers[i]))
      
      if(length(indexNo) <= 0)
      {
        indexToDelete = c(indexToDelete, i)
      }
    }
    
    if(length(indexToDelete) > 0)
    {
      AllMarkers = AllMarkers[-c(indexToDelete)]
    }
    
    cellCount = length(data[1,])
    totalGenes = length(data[,1])
  }
}

#Step 3: Correlation matrix
{
  filteredData = t(data[AllMarkers,])
  
  corMatrix = cor(filteredData)
  
  #for heatmap.2 () function
  library(gplots)
  library(RColorBrewer)
  
  Cairo(file="CorMatrix.png", type="png", units="in", width=70, height=70, pointsize=12*300/72, dpi=300)
  heatmap.2(corMatrix, col=brewer.pal(11,"RdYlBu"), density.info = "none", trace = "none", cexRow = 2.5, cexCol = 2.5, offsetRow = 0, offsetCol = 0, margin=c(10, 10), dendrogram = "row", lwd = 10)
  dev.off()
  
  corMatrixGraphMatrix = corMatrix
  diag(corMatrixGraphMatrix) = 0
}

#Step 4: Set the colors and classes
{
  classOfGenes = rep(0, length(AllMarkers))
  geneColorCode = rep(0, length(AllMarkers))
  textColorCode = rep(0, length(AllMarkers))
  
  colorSet = c("red","grey","darkkhaki","darkorchid1","green") 
  
  for (i in 1:length(AllMarkers))
  {
    classOfGenes[i] = which(listOfImportantGenes == AllMarkers[i], arr.ind = TRUE)[2]
    geneColorCode[i] = colorSet[classOfGenes[i]]
  }
  
  jpeg("colorCode.jpg")
  plot(rep(5,length(colorSet)), seq(1,length(colorSet)), pch=21, bg = colorSet, col = colorSet, cex = 3, axes = FALSE, xlab = "", ylab = "")
  text(rep(5.5,length(colorSet)), seq(1,length(colorSet)), labels = colnames(listOfImportantGenes))
  dev.off()
}

corGraph = graph_from_adjacency_matrix(corMatrixGraphMatrix, mode="upper", weighted = TRUE, diag = FALSE)

normalized = 1+(10*(E(corGraph)$weight - min(E(corGraph)$weight))/(max(E(corGraph)$weight) - min(E(corGraph)$weight)))

E(corGraph)$color = brewer.pal(11,"RdYlBu")[normalized]



Cairo(file="grahColorCode.png", type="png", units="in", width=40, height=40, pointsize=12*300/72, dpi=300)
plot(rep(1,11), pch = 21, col = brewer.pal(11,"RdYlBu"), bg = brewer.pal(11,"RdYlBu"), xlim = c(0,12), ylim = c(0,2), cex = 3)
dev.off()

Cairo(file="CorelationGraph.png", type="png", units="in", width=40, height=40, pointsize=12*300/72, dpi=300)
tkplot(corGraph, vertex.color = geneColorCode, vertex.label=AllMarkers, vertex.label.color= "black", vertex.size= 20, edge.width = normalized/2)
dev.off()


Cairo(file="MSTCorelationGraph.png", type="png", units="in", width=40, height=40, pointsize=12*300/72, dpi=300)
tkplot(mst(corGraph), vertex.color = geneColorCode, vertex.label=AllMarkers, vertex.label.color= "black", vertex.size= 20, edge.width = normalized/2)
dev.off()