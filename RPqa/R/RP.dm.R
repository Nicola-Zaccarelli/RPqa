RP.dm = function(x, y=NULL, name="data", type= "RP", m=2, d=1, start.time= 1, end.time=c(length(x), length(y)),
           dist.type = "euclidean", normalize= TRUE, plotD=FALSE)
    {
    #
    # version: 10-09-2012
    #
    #---------------------------------------------------------------------------------
    # This should be for computing RP and CRP
    # when using CRP:
    # - length of time series can be different
    # - the biggest embedding dimension is use
    # - deley can be different
    # - for start and end time recycling is on if length is = 1
    #
    # DIFFERENT PLOTs are produced to account for distances and distance distributions
    #
    # For JRP runs for each time series and then use the specific function to visualise
    #
    # PS: the package peoxy is required for CRP
    #
#######################################################################################
    #
    # Some checking on the arguments
    # To be reviewed
    #
    stopifnot( is.null(x) | (length(x) > 10 & class(x) %in% c("numeric", "ts")),
               !is.null(name),
               dist.type %in% c("euclidean", "maximum", "minimum", "normalized"),
               plotD %in% c("FALSE", "F", "T", "TRUE"),
               normalize %in% c("FALSE", "F", "T", "TRUE"), 
               type %in% c("RP", "CRP"))
    
    if(is.null(y)) {
      if(type !="RP") stop("Wrong type of analysis, use type=RP")
      end.time = end.time[-2] # delete last zero
      stopifnot(length(m) == 1, m >= 2, m == floor(m),
                length(d) == 1, d >= 1, d == floor(d),
                length(start.time) == 1, start.time >= 1, start.time == floor(start.time),
                length(end.time) == 1, end.time >= 1, end.time == floor(end.time),
                end.time > start.time)
      if((length(x) <= m * d) | (end.time - start.time <= m * d)) stop("The x series is TOO SHORT with such m and d!")
    } else {
      if(type !="CRP") stop("Wrong type of analysis, use type= CRP.")
      if((length(y) > 10 & class(y) %in% c("numeric", "ts"))== FALSE) stop("Wrong length or data type for y.")
      if(((length(m) %in% c(1, 2)) & class(m)=="numeric")==FALSE) stop("Embedding dimension uncorrect.")
      if(length(m)==2) {
        cat("For CRP the maximum embedding dimension is used. \n")
        m=max(m)
        if((m >= 2 & m == floor(m))==FALSE) stop("Embedding dimension uncorrect.")
      } else {if((m >= 2 & m == floor(m))==FALSE) stop("Embedding dimension uncorrect.")}
      if((length(d) %in% c(1, 2) & class(d)=="numeric")==FALSE) stop("Delay dimension uncorrect.")
      if((d >= 1 && d == floor(d))==FALSE) stop("Delay dimension uncorrect.")
      if(length(d) == 1) d = c(d,d) # to help futher analysis
      if(((length(start.time) %in% c(1, 2)) && (start.time >= 1) && (start.time == floor(start.time)))==FALSE) stop("Wrong start.time value.")
      if(((length(end.time) %in% c(1, 2)) && (end.time >= 1) && (end.time == floor(end.time)))==FALSE) stop("Wrong end.time value.")
      if(prod(start.time < end.time) != 1) stop("Check start or end time as they are not compattible.")
      if(length(start.time)==1) start.time=c(start.time, start.time)
      if(length(end.time)==1) end.time=c(end.time, end.time)
      if((length(x) <= m * d) || (end.time - start.time <= m * d)) stop("The x series is TOO SHORT with such m and d!")
      if((length(y) <= m * d) || (end.time - start.time <= m * d)) stop("The y series is TOO SHORT with such m and d!")
   }
############################################################################
# Calculate series
#
   x = x[start.time[1]:end.time[1]]
   if (type =="CRP") y = y[start.time[2]:end.time[2]]
############################################################################
# Normalize data if requested
    if (normalize) {
        x = (x - mean(x))/sd(x)
        if (type =="CRP") y = (y - mean(y))/sd(y)
    }
#
# Calculate embendding series using a specific function
#
   n = length(x) - (m - 1) * d[1]
   embSeriesX = matrix(0, n, m)
   if (type =="CRP") {
     ny = length(y) - (m - 1) * d[2]
     embSeriesY = matrix(0, ny, m)
     for (i in 1:m) embSeriesY[, i] = y[((i - 1) * d[2] + 1):(ny + (i - 1) * d[2])]
   }
   for (i in 1:m) embSeriesX[, i] = x[((i - 1) * d[1] + 1):(n + (i - 1) * d[1])]
#
# Calculate distance matrix
  if (type =="CRP") {
    cat("The package proxy is required. Loading... \n")
    require("proxy", quietly=TRUE, warn.conflicts=FALSE)
    if (dist.type=="minimum") dist.type="manhattan"
    if (dist.type=="maximum") dist.type="supremum"
    if (dist.type=="normalized") {
        for (i in 1:dim(embSeriesX)[1]) embSeriesX[i, ] = embSeriesX[i, ] / sqrt(sum((embSeriesX[i, ])^2)) ## Check
        for (i in 1:dim(embSeriesY)[1]) embSeriesY[i, ] = embSeriesY[i, ] / sqrt(sum((embSeriesY[i, ])^2)) ## Check
        dist.type="euclidean"}
    D = dist(embSeriesX, embSeriesY, method=dist.type)
  } else {
      if (dist.type=="minimum") dist.type="manhattan"
      if (dist.type=="normalized") {
          # Not an elegant solution but it works
          for (i in 1:dim(embSeriesX)[1]) embSeriesX[i, ] = embSeriesX[i, ] / sqrt(sum((embSeriesX[i, ])^2)) ## Check
          dist.type="euclidean"}
      D = dist(embSeriesX , method=dist.type)}
# Provide some basic info on distances
  cat("Settings used: \n")
  cat(" Type of recurrence plot required:", type, "\n")
  if (type=="CRP") cat("Values for the x vs. y series: \n") 
     else  cat("Values for the x series: \n")
  cat("Embedding: ",m, "Delay:", d, "\n")
  cat("Start time: ", start.time, " and end time: ", end.time, "\n")
  cat("Distance type:", dist.type, "\n")
  cat("Summary statistics on distances: \n")
  cat("Values for the x series: \n")
  cat("min: ", min(D), "max", max(D), "\n")
  Dsd= sd(as.vector(D))
  Dmean = mean(D)
  Dmed  = median(D)
  cat("median: ", Dmed, "mean", Dmean, "\n")
  cat("sd: ", Dsd, "\n")
  ecdfD = ecdf(as.vector(D))
  cat("Proportion of distances below 1 sd value:", ecdfD(Dsd) , "\n")
  quaD10 = quantile(as.vector(D), prob=0.1)
  cat("Percentile of 0.1 of length ditribution:", quaD10 , "\n")
  if (plotD & type=="RP") {
      def.par =  par(no.readonly = TRUE)
      layout(matrix(c(1, 3, 2, 3), 2, 2), c(3, 1), c(12, 1), respect=FALSE)
  #
  # Distance image in grey colors
  #
      par(mar=c(4,5,3,0))
      image(c(start.time:(start.time + n)), c(start.time:(start.time + n)), as.matrix(D), col=grey(1:11 /11),
              xlab=expression(x[i]), ylab=expression(x[i]), main="Recurrence Distance Plot (RP)", asp=1, cex.lab=1.5)
      lines(c(start.time:(start.time + n)), c(start.time:(start.time + n)), col="green")
  #
  # Distribution of distances d with statistics
  #
      histD=hist(D, plot=FALSE)
      midpoints = histD$mids
      n=length(midpoints)
      midpoints= c(histD$breaks[1], midpoints, histD$breaks[n+1])
      histDens= c(0, histD$counts/sum(histD$counts), 0)
      par(mar=c(4,3,3,1))
      plot(histDens, midpoints, type="b", main="Dist. distrib.", xlab="pdf", ylab="") #, ylab="distance")
      abline(h=Dmed, lty=3, col="red") # median distance value
      text(max(histDens)/2, 1.1*Dmed, label="median(d)", col="red", cex=1)
      abline(h=Dsd, lty=4, col="blue")    # sd distance value
      text(0.8*max(histDens), 1.1*Dsd, label="sd(d)", col="blue", cex=1)
      abline(h=quaD10, lty=4, col="darkolivegreen")    # 0.1 percentile distance value
      text(0.3*max(histDens), 0.7*quaD10, label="10%", col="darkolivegreen", cex=1)
  #
  # Legend of distance image
  #
      par(mar=c(0,0,0,0))
      plot(c(0, 12), c(0, 2), type="n", axes=F, xlab="", ylab="")
      xl=c(0:10)
      rect(xl, 1, (xl+1), 2, col=grey(1:11 /11), border="black")
      text(0.5, 0.5, label="smaller dist. (0)")
      text(11, 0.5, label=paste("bigger dist. (", round(max(D), 3), ")"))
      par(def.par)
    }
if (plotD & type=="CRP") {
  def.par =  par(no.readonly = TRUE)
  layout(matrix(c(1, 3, 2, 3), 2, 2), c(3, 1), c(12, 1), respect=FALSE)
  #
  # Distance image in grey colors
  #
  par(mar=c(4,5,3,0))
  image(c(start.time[1]:(start.time[1] + n)), c(start.time[2]:(start.time[2] + ny)), as.matrix(D), col=grey(1:11 /11),
        xlab=expression(x[i]), ylab=expression(x[i]), main="Recurrence Cross Distance Plot (CRP)", asp=1, cex.lab=1.5)
  abline(b=1, col="green")
  #
  # Distribution of distances d with statistics
  #
  histD=hist(D, plot=FALSE)
  midpoints = histD$mids
  n=length(midpoints)
  midpoints= c(histD$breaks[1], midpoints, histD$breaks[n+1])
  histDens= c(0, histD$counts/sum(histD$counts), 0)
  par(mar=c(4,3,3,1))
  plot(histDens, midpoints, type="b", main="Cross Dist. distrib.", xlab="pdf", ylab="") #, ylab="distance")
  abline(h=Dmed, lty=3, col="red") # median distance value
  text(max(histDens)/2, 1.1*Dmed, label="median(d)", col="red", cex=1)
  abline(h=Dsd, lty=4, col="blue")    # sd distance value
  text(0.8*max(histDens), 1.1*Dsd, label="sd(d)", col="blue", cex=1)
  abline(h=quaD10, lty=4, col="darkolivegreen")    # 0.1 percentile distance value
  text(0.3*max(histDens), 0.7*quaD10, label="10%", col="darkolivegreen", cex=1)
  #
  # Legend of distance image
  #
  par(mar=c(0,0,0,0))
  plot(c(0, 12), c(0, 2), type="n", axes=F, xlab="", ylab="")
  xl=c(0:10)
  rect(xl, 1, (xl+1), 2, col=grey(1:11 /11), border="black")
  text(0.5, 0.5, label="smaller dist. (0)")
  text(11, 0.5, label=paste("bigger dist. (", round(max(D), 3), ")"))
  par(def.par)
}
  if(type=="CRP") detach(package:proxy, unload=TRUE) # to have back the original dist function
  return(D)
}
