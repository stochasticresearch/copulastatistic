# Conditional CoS

conditional_cosdv3d <- function(x, y, z, conditional_flag, num_bins=10) {
  if(conditional_flag==1) {
    cc = x
    x1 = y
    x2 = z
  }
  else if(conditional_flag==2) {
    cc = y
    x1 = x
    x2 = z
  }
  else if(conditional_flag==3) {
    cc = z
    x1 = x
    x2 = y
  }
  
  hist_obj = hist(cc,breaks=seq(min(cc),max(cc),(max(cc)-min(cc))/num_bins))
  br = hist_obj$breaks
  dep_total = 0
  total_segments_processed = 0
  for (ii in 1:(length(br)-1)) {
    rngLo = br[ii]
    rngHi = br[ii+1]
    
    # find the indices of all the elements that match this range
    idxs = which(cc>=rngLo & cc<rngHi)
    x_subset = x1[idxs]
    y_subset = x2[idxs]
    
    if(length(x_subset)>2 & length(y_subset)>2) {
      dep = cosdv(x_subset,y_subset)
      dep_total = dep_total + dep
      total_segments_processed = total_segments_processed + 1
    }
    
  }
  dep_total = dep_total/total_segments_processed
  return(dep_total)
}