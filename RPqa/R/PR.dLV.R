RP.dLV = function(dist.mat, eps=1, theiler= 0){
    #
    # version: 11-09-2012
    #
    # dist.amt: distance metrix
    # eps     : radiuos of reserach for neibourghing points
    # mlen    : minimum length of lines
    # theiler : window for correletaed points
    #
    # ASSUMES LOI from N[1, 1] to N[n,n]
    #
    # Insert some checking
    stopifnot( is.null(dist.mat) | (class(dist.mat)=="dist"),
               is.null(eps) | (length(eps) == 1 & eps > 0),
               is.null(theiler) | (length(theiler) == 1 & theiler >= 0),
               theiler == floor(theiler)
    )
    # Calculate length
    n = (1+(1+8*length(dist.mat))^0.5)/2
    # Calculate series
    D = as.matrix(dist.mat)
    #
    # Calculate recurrance matrix
    D = matrix(as.integer(D <= eps), n, n)
    #
    # Calculate length distribution
    # 
    dLV = matrix(0, 2, n)
    row.names(dLV) = c("dL", "dV")
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
    colnames(dLV)=dLL
    delCol= colSums(dLV)*dLL
    delCol = -1*(delCol == 0)*dLL
    dLV[, delCol] # get rid of void colums
    rownames(dLV) = c("Diagonal", "Vertical")
    return(dLV)
}
