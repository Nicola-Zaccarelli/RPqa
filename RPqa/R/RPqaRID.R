RPqaRID = function(dist.mat, eps=1, mlen= 2, theiler= 0, plotRP=TRUE){
    #
    # version: 11-09-2012
    #
    #
    # dist.amt: distance metrix
    # eps     : radiuos of reserach for neibourghing points
    # mlen    : minimum length of lines
    # theiler : window for correletaed points
    #           by default the LOI is excluded
    #
    # here the LOI goes from N[1, 1] to N[n ,n]
    #
    # Insert some checking
    stopifnot( is.null(eps) | (length(eps) == 1 & eps > 0),
               is.null(mlen) | (length(mlen) == 1 & mlen > 1),
               is.null(theiler) | (length(theiler) == 1 & theiler >= 0),
               mlen == floor(mlen), theiler == floor(theiler),
               plotRP %in% c("TRUE", "T", "FALSE", "F")
    )
    # Calculate length
    n = (1+(1+8*length(dist.mat))^0.5)/2
    # Calculate series
    D = as.matrix(dist.mat)
    #
    # Calculate recurrance matrix
    D = matrix(as.integer(D <= eps), n, n)
    #
    # Calculate RQAs
    # Build amatrix dLV for holding distribution results
    # We are looking at half of the matrix excluding LOI
    dLV = matrix(0, 2, n)
    row.names(dLV) = c("dL", "dV")
    # let's calculate dL
    for (i in (2+theiler):n) {
        l =0
        for (j in 1:(n-i+1))
            if (D[(i+j -1), j] == 1) l=l+1 else { 
                dLV[1, 1] = dLV[1, 1] + 1
                if (l>0) { dLV[1, (l+1)] = dLV[1,(l+1)] + 1
                           l=0 } 
            }
        if (l > 0) dLV[1, (l+1)] = dLV[1,(l+1)] + 1
    }
    # let's calculate dV going through columns
    # the Thieler correction is in
    # the -1 is to exclude the LOI
    for (j in 1:(n-theiler-1)) {
        l =0
        for (i in (j+1+theiler):n)
            if (D[i, j] == 1) l=l+1 else { 
                dLV[2, 1] = dLV[2, 1] + 1
                if (l>0) { dLV[2, (l+1)] = dLV[2,(l+1)] + 1
                           l=0 } 
            }
        if (l > 0) dLV[2, (l+1)] = dLV[2,(l+1)] + 1
    }
    dLV = dLV[, -1]
    dLL = c(1:(n-1)) #vector for calculation on length
    # Compute RQA based on Webber!!!
    recur = sum(dLV[1, ] * dLL)
    cat("#recur: ", recur, "\n")
    diag = sum(dLV[1, mlen:(n-1)])
    cat("#diag: ", diag, "\n")
    vert = sum(dLV[2, mlen:(n-1)])
    cat("#vert: ", vert, "\n")
    prec = recur/((n^2 -n)/2)
    cat("%rec (RR): ",prec, "\n")
    det = sum(dLV[1, mlen:(n-1)]*dLL[mlen:(n-1)])/recur
    cat("%DET: ", det, "\n")
    lam = sum(dLV[2, mlen:(n-1)]*dLL[mlen:(n-1)])/recur
    cat("%lam: ", lam, "\n")
    ratio = det/prec
    cat("ratio: ", ratio, "\n")
    meanL = sum(dLV[1, mlen:(n-1)]*dLL[mlen:(n-1)])/sum(dLV[1, mlen:(n-1)])
    cat("Average L: ", meanL, "\n")
    TT = sum(dLV[2, mlen:(n-1)]*dLL[mlen:(n-1)])/sum(dLV[2, mlen:(n-1)])
    cat("Trapping Time (TT): ", TT, "\n")
    lmax = max((dLV[1, mlen:(n-1)] > 0)*dLL[mlen:(n-1)])
    cat("l max: ", lmax, "\n")
    vmax = max((dLV[2, mlen:(n-1)] > 0)*dLL[mlen:(n-1)])
    cat("v max: ", vmax, "\n")
    DIV = 1/lmax
    cat("Divergence (DIV): ", DIV, "\n")
    probL = dLV[1, mlen:(n-1)]/sum(dLV[1, mlen:(n-1)])
    probL[probL==0]= 1 #to get rid of zeros in entropy
    H = -1*sum(probL*log(probL, base = 2))
    Hn = H / log(length(probL), base = 2)
    cat("Entropy (base 2):", H, " and its normalized form: ", Hn,"\n")
    # TREND =
    # it is missing
    #
    # Graph
    #
    if (plotRP) {
        image(c(1:n), c(1:n), D, col=c("white", "black"), xlab=expression(x[i]), ylab=expression(x[i]), main="Recurrence Plot", asp=1, cex.lab=1.5)
        lines(c(0, n), c(0, n), col="red", lwd=1)
    }
    # Repack the dLV matrix
    colnames(dLV)=dLL
    delCol= colSums(dLV)*dLL
    delCol = -1*(delCol == 0)*dLL
    dLV[, delCol] # get rid of void colums
    rownames(dLV) = c("Diagonal", "Vertical")
    #
    # Build a list of results
    #
    return(list(dLV=dLV, recur=recur, diag=diag, vert=vert, prec=prec, det=det, lam=lam, ratio=ratio,  meanL=meanL, TT=TT, lmax=lmax, vmax=vmax, DIV=DIV, H=H, Hn=Hn))
}
