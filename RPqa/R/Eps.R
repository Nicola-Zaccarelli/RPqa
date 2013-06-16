Eps = function(dist.mat, type="perc.pd", param = 0.1){
  #
  # version: 10-09-2012
  #
  # dist.mat: distance metrix
  # type    : type of algorithm for threshold value, possible options
  #   * max.pd, maximum phase diameter 
  #   * mean.pd, mean phase diameter
  #   * median.pd, media phase diameter
  #   * perc.pd, percentiles of phase diameter param in val
  #   * perc.rr, percentiles of rr param in (0.01 - 0.1)
  #   * fix.nn, fixed number of neighbour param = nn
  #   * noise, 5*base noise  NOT THERE
  #   * range, explore the effect of changing eps within a range on RR or DET
  #
  # Some checking on the arguments
  # To be reviewed
  #
  stopifnot( class(dist.mat)=="dist",
             type %in% c("max.pd", "mean.pd", "median.pd", "perc.pd", "perc.rr", "fix.nn", "noise", "range"),
             sum(param >0) == length(param))
  # Calculate length
  n = attr(dist.mat, "Size")
  #
  if (type=="perc.pd") stopifnot(length(param) == 1,param > 0,  param <= 0.15)
  if (type=="perc.rr") stopifnot(length(param) == 1,param > 0,  param <= 1)
  if (type=="fix.nn")  stopifnot(length(param) == 1,param > 0,  param <= n/10, param == floor(param))
  if (type=="range")   stopifnot(length(param) == 3, sum(param > 0) == 3,  param[1] < param[2])
  #
  # Evaluation part
  #
  max.pd = max(dist.mat)
  mean.pd=mean(dist.mat)
  # median.pd=median(dist.mat)
  cat("Summary statistics on distances: \n")
  cat("Values for the x series: \n")
  cat("min: ", min(dist.mat), "max", max.pd, "\n")
  Dsd= sd(as.vector(dist.mat))
  cat("mean", mean.pd, "\n")
  cat("sd: ", Dsd, "\n")
  ris = list(max.pd=max.pd, mean.pd=mean.pd, type=type)
  if (type=="median.pd") {
    median.pd=media(dist.mat)
    ris = list(max.pd=max.pd, mean.pd=mean.pd, median.pd=median.pd, type=type)
    cat("Median phase state diameter: ", median.pd, "\n")
  }
  if (type=="perc.pd") {
    prob = c(param, seq(0.01, 0.1, 0.01))
    perc = max.pd * prob
    ris = cbind(prob, perc)
  }
  if (type=="perc.rr") {
    ecdfD = ecdf(as.vector(dist.mat))
    prob = c(param, seq(0.05, 0.25, 0.05))
    perc = quantile(as.vector(dist.mat), prob=prob)
    ris = cbind(prob, perc)
  }
  if (type=="fix.nn") {
    fix.nn = NULL
    for (i in 1:(n-1)){ 
      vect.dist = NULL
      for (j in 1:(n-1)){
        vect.dist = c(vect.dist, dist.mat[n*(i-1) - i*(i-1)/2 + j-i])
      }
      vect.dist = sort(vect.dist, na.last=TRUE, method="quick")
      fix.nn = c(fix.nn, vect.dist[(param-1)])
    }
    ris = fix.nn
  }
  # if(length(eps)>1) {
  #   for (i in 1:(n-1)) 
  #      for (j in 1:(n-1)) if(D[n*(i-1) - i*(i-1)/2 + j-i] <= eps(i)) D=1 else D=0
  #   D = as.matrix(D)
  # }
  if (type=="noise") {
    noise = 1
  }
  if (type=="range") {
    range = seq(param[1], param[2], param[3])
    ris = matrix(0, length(range), 7)
    colnames(ris) =c("eps", "recur", "prec", "det", "lam", "TT", "DIV")
    for (i in 1:length(range)) {
    risTMP = RPqaRID(dist.mat, eps=range[i], plotRP=FALSE)
    ris[i, 1] = range[i]
    ris[i, 2] = risTMP$recur
    ris[i, 3] = risTMP$prec
    ris[i, 4] = risTMP$det
    ris[i, 5] = risTMP$lam
    ris[i, 6] = risTMP$TT
    ris[i, 7] = risTMP$DIV
    }
    ris = list(ris=ris, type=type)
  }
  return(ris)
}
