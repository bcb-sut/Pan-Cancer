sig <- c('1.1', '1.2', '1.3', '1.4','1.5.1', '1.5.2','1.5.3', '1.5.4', '1.5.5', '1.5.6','1.5.7', '1.5.8', '1.6', '1.7', '1.8', '2.1', '2.2')
for (i in sig[7]){
  signature_data <- read.table(paste0("mutational signatures analysis/output/resulst_for_clust" , i , "/Evaluation.txt"),  sep = " ")
  signature_data$V3 <- signature_data$V3 /1000
  V1 <- signature_data$V1
  V2 <- signature_data$V2
  V3 <- signature_data$V3
  par(mar=c(5, 4, 4, 6) + 0.1)
  
  ## Plot first set of data and draw its axis
  plot(V1, V2,  axes=FALSE, ylim=c(0,1), xlab="", ylab="", type="b",col="black")
  axis(2, ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
  mtext("RECONSTRUCTION ERROR",side=2,line=2.5)
  box()
  
  ## Allow a second plot on the same graph
  par(new=TRUE)
  
  ## Plot the second plot and put axis scale on right
  plot(V1, V3, xlab="", ylab="", ylim=c(0.5,1.5), 
       axes=FALSE, type="b", col="red")
  ## a little farther out (line=4) to make room for labels
  mtext("SIGNATURE REPRODUCIBILITY",side=4,col="red",line=4) 
  axis(4, ylim=c(1,4), col="red",col.axis="red",las=1)
  ## Draw the time axis
  axis(1,1:15, 1:15)
  mtext("N",side=1,col="black",line=2.5)  
}

