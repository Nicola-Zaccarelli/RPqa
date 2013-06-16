RPqamw = function(dist.mat, eps=1, mlen= 2, mw= 0, step=0, theiler= 0, plotG=FALSE){
    #
    # version: 10-09-2012
    #
    # Moving window version of the analysis
    # Theiler window is now implemented
    # Insert some checking
    stopifnot( is.null(dist.mat) | (class(dist.mat)=="dist"), 
               is.null(eps) | (length(eps) == 1 & eps > 0),
               is.null(mlen) | (length(mlen) == 1 & mlen > 1),
               is.null(mw) | (length(mw) == 1 & mw >=0 ),
               is.null(step) | (length(step) == 1 & step >= 0),
               is.null(theiler) | (length(theiler) == 1 & theiler >= 0),
               plotG %in% c("FALSE", "F", "T", "TRUE"),
               theiler = floor(theiler), mw= floor(mw), eps= floor(eps),
               mlen = floor(mlen), step = floor(step), mw > theiler
    )
    # Calculate length
    n = (1+(1+8*length(dist.mat))^0.5)/2
    # Calculate series
    D = as.matrix(dist.mat)
    # Calculate recurrance matrix
    D = matrix(as.integer(D <= eps), n, n)
    # Calculate number of moving windows
    Nmw= floor(n/mw)
    if (Nmw == 0) stop("Moving window too big!")
    RepMW = floor((n -mw)/step)
    if (RepMW == 0) stop("Step too big with this moving window!")
    #
    # Build the container of the results: a matrix
    #
    ris = matrix(0, 14, RepMW)
    rownames(ris) = c("recur", "diag", "vert", "prec", "det", "lam", "ratio", "meanL", "TT", "lmax", "vmax", "DIV", "H", "Hn")
    #
    # Calculate the values
    #
    for (w in 1:RepMW) {
    # Build amatrix dLV for holding distribution results
    # We are looking at half of the matrix excluding LOI
    Dmw = D[(n-mw-w):(n-w+1), (n-mw-w):(n-w+1) ]
    dLV = matrix(0, 2, mw)
    # let's calculate dL
    for (i in (2+theiler):mw) {
      l =0
      for (j in 1:(mw-i+1))
        if (Dmw[(i+j -1), j] == 1) l=l+1 else { 
          dLV[1, 1] = dLV[1, 1] + 1
          if (l>0) { dLV[1, (l+1)] = dLV[1,(l+1)] + 1
                     l=0 } 
        }
      if (l > 0) dLV[1, (l+1)] = dLV[1,(l+1)] + 1
    }
    # let's calculate dV going through columns
    for (j in 1:(mw-theiler-1)) {
      l =0
      for (i in (j+1+theiler):mw)
        if (Dmw[i, j] == 1) l=l+1 else { 
          dLV[2, 1] = dLV[2, 1] + 1
          if (l>0) { dLV[2, (l+1)] = dLV[2,(l+1)] + 1
                     l=0 } 
        }
      if (l > 0) dLV[2, (l+1)] = dLV[2,(l+1)] + 1
    }
    #  dLV = dLV * 2 # Because we have used only half of the matrix
    #  dLV[1, n+1] = dLV[1, n+1]/2
    dLV = dLV[, -1]
    dLL = c(1:(mw-1)) #vector for calculation on length
    # Compute RQA based on Webber!!!
    ris[1, w] = sum(dLV[1, ] * dLL)
    ris[2, w] = sum(dLV[1, mlen:(mw-1)])
    ris[3, w] = sum(dLV[2, mlen:(mw-1)])
    ris[4, w] = ris[1, w]/((mw^2 -mw)/2)
    ris[5, w] = sum(dLV[1, mlen:(mw-1)]*dLL[mlen:(mw-1)])/ris[1, w]
    ris[6, w] = sum(dLV[2, mlen:(mw-1)]*dLL[mlen:(mw-1)])/ris[1, w]
    ris[7, w] = ris[5, w]/ris[4, w]
    ris[8, w] = sum(dLV[1, mlen:(mw-1)]*dLL[mlen:(mw-1)])/sum(dLV[1, mlen:(mw-1)])
    ris[9, w] = sum(dLV[2, mlen:(mw-1)]*dLL[mlen:(mw-1)])/sum(dLV[2, mlen:(mw-1)])
    ris[10, w] = max((dLV[1, mlen:(mw-1)] > 0)*dLL[mlen:(mw-1)])
    ris[11, w] = max((dLV[2, mlen:(mw-1)] > 0)*dLL[mlen:(mw-1)])
    ris[12, w] = 1/ris[11, w]
    probL = dLV[1, mlen:(mw-1)]/sum(dLV[1, mlen:(mw-1)])
    probL[probL==0]= 1 #to get rid of zeros in entropy
    ris[13, w] = -1*sum(probL*log(probL))
    ris[14, w] = ris[13, w]/log(length(probL))
    }
    if (plotG) {
      def.par =  par(no.readonly = TRUE)
      layout(matrix(c(1, 2, 3, 4), 2, 2), c(1, 1), c(1, 1), respect=FALSE)
      par(mar=c(3,4,1,1))
      plot(c(1:RepMW), ris[1, ], type="l", xlab="", ylab="RR")
      abline(h=mean(ris[1,]), lty=3)
      par(mar=c(4,4,0,1))
      plot(c(1:RepMW), ris[5, ], type="l", xlab="", ylab="DET")
      abline(h=mean(ris[5,]), lty=3)
      mtext("moving window", side=1, line=2)
      par(mar=c(3,4,1,1))
      plot(c(1:RepMW), ris[6, ], type="l", xlab="", ylab="LAM")
      abline(h=mean(ris[6,]), lty=3)
      par(mar=c(4,4,0,1))
      plot(c(1:RepMW), ris[9, ], type="l", xlab="", ylab="TT")
      abline(h=mean(ris[9,]), lty=3)
      mtext("moving window", side=1, line=2)
      mtext(paste("Moving window dimension of", mw, "with step of", step), side=1, line=3, cex=0.8, font=2)
      par(def.par)
    }
    return(ris)
  }
