RPqa.boot = function(dLV, min.length=2, n.samples=99){
  #
  # No checking is performed
  #
  calcIndex = function(dati, mlen, indices){
  distB= dati[indices]
  index1 = mean(distB[distB >= mlen]) # mean value
  index2 = sum(distB[distB >= mlen])/(sum(distB)) # proportion 
  return(c(index1, index2))
}
  bootVector = rep(1:dim(dLV)[2], dLV[1,]) # build sample of L diagonals
  risT = boot(bootVector, calcIndex, R=n.samples, mlen=min.length)
  ris = rbind(risT$t0, risT$t)
  bootVector = rep(1:dim(dLV)[2], dLV[2,]) # build sample of vertical lines
  risT = boot(bootVector, calcIndex, R=n.samples, mlen=min.length)
  ris = cbind(ris, rbind(risT$t0, risT$t))
  rownames(ris) = c("Orig", paste("s", 1:n.samples, sep = ""))
  colnames(ris) = c("Lmean", "DET", "Vmean", "LAM")

  cat("ORDINARY NONPARAMETRIC BOOTSTRAP \n")
  cat("\n")
  cat("Bootstrap statistics:\n")
  cat("        original     bias       std.error\n")
  cat("Lmean ", ris[1, 1]," ", mean(ris[2:(n.samples+1), 1])-ris[1, 1]," ", sd(ris[2:(n.samples+1), 1]), "\n")
  cat("DET   ", ris[1, 2]," ", mean(ris[2:(n.samples+1), 2])-ris[1, 2]," ", sd(ris[2:(n.samples+1), 2]), "\n")
  cat("Vmean ", ris[1, 3]," ", mean(ris[2:(n.samples+1), 3])-ris[1, 3]," ", sd(ris[2:(n.samples+1), 3]), "\n")
  cat("LAM   ", ris[1, 4]," ", mean(ris[2:(n.samples+1), 4])-ris[1, 4]," ", sd(ris[2:(n.samples+1), 4]), "\n")
  
return(ris)}
