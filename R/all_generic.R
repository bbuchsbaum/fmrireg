parentTerms <- function(x) UseMethod("parentTerms")

cells <- function(x, ...) UseMethod("cells")

conditions <- function(x, ...) UseMethod("conditions")

convolve <- function(x, hrf, samplingFrame,...) UseMethod("convolve")

isContinuous <- function(x) UseMethod("isContinuous")

levels <- function(x) UseMethod("levels")

columns <- function(x) UseMethod("columns")

designMatrix <- function(x, ...) UseMethod("designMatrix")

elements <- function(x, ...) UseMethod("elements")

evaluate <-  function(x, samplingGrid, ...) UseMethod("evaluate")

nbasis <-  function(x) UseMethod("nbasis")

onsets <-  function(x) UseMethod("onsets")

durations <-  function(x) UseMethod("durations")

amplitudes <-  function(x) UseMethod("amplitudes")

samples <-  function(x, ...) UseMethod("samples")