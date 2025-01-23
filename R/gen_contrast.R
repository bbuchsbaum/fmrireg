#' @keywords internal
#' @noRd
#' @examples
#' des <- data.frame(time=factor(rep(1:4, 1)), condition=factor(rep(c("face", "scene"), each=2)))
#' con <- generate_interaction_contrast(des, c("time", "condition"))
generate_interaction_contrast <- function(des, factors) {
  assert_that(all(factors %in% names(des)), msg=paste("Contrast factors:", factors, "--> not found in design table."))
  
  # Create lists of contrast vectors/matrices for each factor
  Clist <- lapply(names(des), function(ev) rep(1, nlevels(des[[ev]])))
  Dlist <- lapply(names(des), function(ev) t(-diff(diag(nlevels(des[[ev]])))))
  
  nfac <- length(Clist)
  mats <- vector("list", nfac)
  ind <- match(factors, names(des))
  
  # Assign difference contrasts to specified factors
  mats[ind] <- Dlist[ind]
  mats[-ind] <- Clist[-ind]  # Assign ones to other factors
  
  # Compute the Kronecker product
  cmat <- Reduce(kronecker, mats)
  
  # Correct dimension check
  if (nrow(cmat) != prod(sapply(des, nlevels))) {
    stop("Contrasts: design table and contrast matrix have different number of rows.")
  }
  
  return(cmat)
}

#' @keywords internal
#' @noRd
#' @examples
#' des = data.frame(time=factor(c(1,1,2,2,3,3)), condition=factor(rep(c("face", "scene"), each=3)))
#' con <- generate_main_effect_contrast(des, "time")
generate_main_effect_contrast <- function(des, factor) {
  assert_that(factor %in% names(des), msg=paste("Contrast factor:", factor, " --> not found in design table."))
  
  # Create lists of contrast vectors/matrices for each factor
  Clist <- lapply(names(des), function(ev) rep(1, nlevels(des[[ev]])))
  Dlist <- lapply(names(des), function(ev) t(-diff(diag(nlevels(des[[ev]])))))
  
  nfac <- length(Clist)
  i <- match(factor, names(des))
  
  mats <- vector("list", nfac)
  mats[[i]] <- Dlist[[i]]    # Assign difference contrast to the main effect factor
  mats[-i] <- Clist[-i]      # Assign ones to other factors
  
  # Compute the Kronecker product
  ret <- Reduce(kronecker, mats)
  
  # Correct dimension check
  if (nrow(ret) != prod(sapply(des, nlevels))) {
    stop("Contrasts: design table and contrast matrix have different number of rows.")
  }
  
  return(ret)
}
