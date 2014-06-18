#'
#' @title Adds some misclassification error to binary data
#' @description Introduces some misclassification in binary data. 
#' The number of values misclassified is determined by the misclassification rates.
#' @param binary.vector a vector of binary values.
#' @param error.1.0 1 to 0 misclassification rate.
#' @param error.0.1 0 to 1 misclassification rate.
#' @return a binary vector
#' @keywords internal
#' @author Gaye A.
#'
misclassify <- function (binary.vector, error.1.0=0.05, error.0.1=0.05){
  
  A <- binary.vector
  A.initial <- A
  misclassification.rate.1.to.0 <- error.1.0
  misclassification.rate.0.to.1 <- error.0.1
  
  # DISPLAY AN ALERT AND STOP THE PROCESS IF THE VECTOR IS NOT BINARY
  if(length(table(A)) != 2)
  {
    cat("\n\n ALERT!\n")
    cat(" The vector you provided is not a binary vector\n")
    stop(" Process ended\n\n", call.=FALSE)
  }
  
  # IF VALUES ARE NOT 0S AND 1S CONVERT TO 0 AND 1 VALUES
  if(names(table(A.initial)[1]) != "0" || names(table(A.initial)[2]) != "1")
  { 
    val1 <- as.numeric(names(table(A))[1])
    val2 <- as.numeric(names(table(A))[2])
    x1 <- which(A == names(table(A))[1])
    y2 <- which(A == names(table(A))[2])
    A[x1] <- 0
    A[y2] <- 1
  }
  
  # GET THE NUMBER OF "1" VALUES IN THE VECTOR
  num.A.1 <- length(which(A == 1))
  
  # GET THE NUMBER OF "1" POSITIONS TO MISCLASSIFY (DEPENDS ON MISCLASSIFICATION RATE)
  num.A.1.missed <- round(misclassification.rate.1.to.0 * num.A.1)
  
  # GET THE INDICES OF THE "1" VALUES 
  A.1.indx <- which(A == 1)
  
  # GET RANDOMLY THE INDICES OF THE "1" VALUES TO MISCLASSIFY
  A.1.to.change <- sample(A.1.indx, num.A.1.missed, replace=FALSE)
  
  # GET THE NUMBER OF "0" VALUES IN THE VECTOR
  num.A.0 <- length(which(A == 0))
  
  # GET THE NUMBER OF "0" POSITIONS TO MISCLASSIFY (DEPENDS ON MISCLASSIFICATION RATE)
  num.A.0.missed <- round(misclassification.rate.0.to.1 * num.A.0)
  
  # GET THE INDICES OF THE "0" VALUES 
  A.0.indx <- which(A == 0)
  
  # GET RANDOMLY THE INDICES OF THE "0" VALUES TO MISCLASSIFY
  A.0.to.change <- sample(A.0.indx, num.A.0.missed, replace=FALSE)
  
  # MISCLASSIFY (1 to 0 and 0 to 1)
  A[A.1.to.change] <- 0
  A[A.0.to.change] <- 1
  
  # IF INITIAL VALUES WERE NOT 0S AND 1S CONVERT BACK
  if(names(table(A.initial)[1]) != "0" || names(table(A.initial)[2]) != "1")
  { 
    x2 <- which(A == 0)
    y2 <- which(A == 1)
    A[x2] <- val1
    A[y2] <- val2
  }
  A.new <- A
  
  # RETURN THE MISCLASSIFIED VECTOR
  misclassified.vector <- A.new
}
