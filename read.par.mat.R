read.par.mat <- function(filename, eegfile)
{
  params <- readMat(filename)$params
  names(params) <- sprintf("%s", attr(params, 'dimnames')[[1]])
  attr(params , "dimnames") <- NULL
  attr(params , "dim") <- NULL
  
  W <- as.numeric(params$W)
  th <- as.numeric(params$th)
  ufeats <- params$feats
    
  
  ch <- to.channels(input(1))
  p0 <- pipe.decimate(ch, 1, 20 , coef_10000_to_500)
  p1 <- pipe.references(p0, c(42,15))
  p2 <- pipe.spatial(p1, M) # online centering with std devision
  
  output(1)$connect(p2)
  
}