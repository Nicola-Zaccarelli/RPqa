dVL <-
function(RPMatrix, thieler=0){
  #
  # We calculate the distribution of diagonal and vetical lines
  # in a RP.
  # FOR NOW THE MATRIX HAS TO BE SQUARE!
  # Input: RPmatrix: a matrix of a recurrence plot or a binary square matrix
  #
  # IMPORTANT: the LOI must go from [N,1] to [1, N]
  #
  ##########################################################################
  #
  # Perform checking
  #
  stopifnot(
    length(thieler) == 1, thieler >= 0, thieler == floor(thieler),
    class(RPMatrix) %in% c("matrix"), dim(RMatrix)[1]==dim(RMatrix)[2]
    )
  if(is.numeric(RMatrix[1,1])==FALSE) stop("The matrix has to be numeric!")
  cat("IMPORTANT: the LOI must go from [N,1] to [1, N] \n")
  num = dim(RPMatrix)[1]
  dLV= matrix(0, 2, num)
  distrLV <- .C("dVL", D=as.integer(RPMatrix), N=as.integer(num), T=as.integer(theiler),
                dLV=as.integer(dLV), PACKAGE="RPqa")$dLV
  
  distrLV = as.matrix(distrLV, 2, num)
  colnames(distrLV)=0:num
  delCol= colSums(distrLV)*(0:num)
  distrLV[, -delCol] # get rid of void colums
  distrLV[, -1]      # get rid of zero colums
  rownames(distrLV) = c("Diagonal", "Vertical")
  return(distrLV)
}
