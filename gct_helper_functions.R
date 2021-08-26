require(cmapR)

#gct_helper_functions.R
#Author: Pierre M Jean-Beltran
#Set of auxiliary functions for GCT object creation and manipulation
#Please feel free to contribute through github

create_gct <- function(mat,rdesc = NULL,cdesc = NULL){
  #Function to create a GCT object. This function will check for the order of the
  #ids in both rdesc and cdesc
  
  if(!is.matrix(mat)){
    error("mat must be a matrix")
  }
  
  if(!is.null(rdesc)){
    if(is.null(rdesc$id)){
      error("rdesc object requires an 'id' column")
    }
    
    if(nrow(mat)!=nrow(rdesc)){
      error("Number of rdesc rows is not equal to number of mat rows")
    }
    
    if(sum(rownames(mat) %in% rdesc$id) != nrow(mat)){
      error("Ids in rdesc do not match rownames of mat")
    }
    
    mat <- mat[rdesc$id,]
  }
  
  if(!is.null(cdesc)){
    if(is.null(cdesc$id)){
      error("cdesc object requires an 'id' column")
    }
    
    if(ncol(mat)!=nrow(cdesc)){
      error("Number of cdesc rows is not equal to number of mat columns")
    }
    
    if(sum(colnames(mat) %in% cdesc$id) != ncol(mat)){
      error("Ids in rdesc do not match rownames of mat")
    }
    
    mat <- mat[,cdesc$id]
  }
  
  return(new('GCT',mat=mat,cdesc=cdesc,rdesc=rdesc))
  
  
}
