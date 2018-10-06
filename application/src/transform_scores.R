#' Transform scores to range from -1 to 1.
#'
#' @param scores a p by p matrix of scores
#' @return a p by p matrix of transformed scores
transform_scores <- function(scores) {
  #S7: transform s so that the scores are in [-1, 1].
  ranges <- range(scores[lower.tri(scores)]) #Ignore diagonal values.
  scores <- 2 * (scores - ranges[1])/(ranges[2] - ranges[1]) - 1
  diag(scores) <- 0
  
  return(scores)
}

scale_to_zero_one <- function(scores) {
  scores <- (transform_scores(scores) + 1) / 2
  diag(scores) <- 0
  return(scores)
}

scores_to_zvalues <- function(scores) {
  if(is.matrix(scores)) {
    scores <- scores[lower.tri(scores)]
  }
  
  scores <- (rank(scores) - 0.5) / length(scores)
  scores <- qnorm(scores)
  return(scores)
}

#' Transform scores to have mean 0 and variance 1.
#'
#' @param scores a vector or symmetric matrix of scores
#' @param ignore_zeroes if true, entries with zeros will be ignored when 
#' centering and scaling. 
#' @param robust if true, robust estimates for mean and variance will be used.
#' @return a vector or symmetric matrix of centered and scaled scores.
standardize_scores <- function(scores, ignore_zeroes = TRUE, robust = FALSE) {
  if(is.matrix(scores)) {
    vals <- scores[lower.tri(scores)]
  } else {
    vals <- scores
  }
  
  if(ignore_zeroes) {
    index <- which(vals != 0)
  } else {
    index <- 1:length(vals)
  }
  if(length(index) > 1) {
    if(robust) {
      mu_est <- median(vals[index])
      sd_est <- 1.4826 * median(abs(vals[index] - median(vals[index]))) 
    } else {
      mu_est <- mean(vals[index])
      sd_est <- sd(vals[index])
    }
    vals[index] <- (vals[index] - mu_est) / sd_est
  }
  
  if(is.matrix(scores)) {
    std_scores <- matrix(0, nrow(scores), ncol(scores))
    std_scores[lower.tri(std_scores)] <- vals # Put values in lower triangle.
    std_scores <- std_scores + t(std_scores) # Symmetrize
    diag(std_scores) <- 0
  } else {
    std_scores <- vals
  }
  
  return(std_scores)
}

#' Transform scores to be in [0, 1]
unsign_scores <- function(scores) {
  return(scale_to_zero_one(abs(scores)))
}


#' Normalize scores
#'
#' @param scores a matrix of cPLS association scores
#' @param kurt_threshold If the kurtosis of scores is above this 
#' value, then the transformation f(x) = sign(x)log(|x| + 1) is applied
#' to the standardized scores.
#' @return a matrix of normalized scores
#' @export
#' @note The scores are first standardized to mean zero and unit variance.
#' The transformation is applied if the scores have very heavy tails. This can
#' increase power when using the fdr procedure. However, if the kurtosis threshold
#' is too low, power will be lost if the transformation is applied. The
#' default value was found in simulation testing to help the severe cases.
normalize_cpls_scores <- function(scores, kurt_threshold = 30) {
  kurtosis <- function(scores) {
    scores <- scores[lower.tri(scores)]
    n <- length(scores)
    mu <- mean(scores)
    n * sum((scores - mu)^4) / (sum((scores - mu)^2)^2)
  }
  
  scores <- log(scores + 1)
  while(kurtosis(scores) > kurt_threshold) {
    scores <- sign(scores) * log(abs(scores) + 1)
    scores <- log(scores + 1)
  }
  return(scores)
}
