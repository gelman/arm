# function to append array  

array.append<-function(array1, array2, d = 3){
  if(any(dim(array1)[-d]!= dim(array2)[-d])){
    stop(message="array dimention must be same for all the dimention except for the one that you are trying to append")
  } else{
    newdim <-  dim(array1)
    newdim[d] <- ifelse(is.na( dim(array1)[d]),1,dim(array1)[d])+ 
      ifelse(is.na( dim(array2)[d]),1,dim(array2)[d])
    newarray <- c(array1,array2)
    dim(newarray) = newdim
  }
  return(newarray)
}
