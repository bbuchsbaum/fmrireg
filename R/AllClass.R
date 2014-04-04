
#' HRF
#' 
#' A class that encapuslates a hemodynamic response function
#' @name HRF-class
#' @slot hrf the underlying function that maps from \empahsis{time} to a response value
#' @slot name the name of the hemodynamic response function
#' @slot nbasis the number of basis functions 
#' @export
setClass("HRF", representation(hrf="function", name="character", nbasis="integer"))
         

#' Regressor
#' 
#' A class representing a regressor variable that is formed by convolving a hemodynamic response function with a set of events.
#' @name Regressor-class
#' @slot onsets the event onsets since the start of the scan in seconds
#' @slot hrf the hemodynamic response function
#' @slot amplitude a vector of values used to scale the response function for each event (default to 1)
#' @slot duration a vector of values indicating the duration of each event
#' @slot span the period of time after the event onset for which the hrf will be evaluated
#' @export
setClass("Regressor",
         representation(onsets="numeric", hrf="HRF", amplitude="vector", duration="vector", span="numeric"),
         contains="function")


#' ConstantVector
#' 
#' A class in which a single value is returned forat all vector indices from 1 to N
#' @name ConstantVector-class
#' @slot length the length of the vector
#' @export
setClass("ConstantVector",
         representation(length="integer"),
         contains="vector")