% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmri_dataset.R
\name{data_chunks.fmri_mem_dataset}
\alias{data_chunks.fmri_mem_dataset}
\title{Create Data Chunks for fmri_mem_dataset Objects}
\usage{
\method{data_chunks}{fmri_mem_dataset}(x, nchunks = 1, runwise = FALSE, ...)
}
\arguments{
\item{x}{An object of class 'fmri_mem_dataset'.}

\item{nchunks}{The number of data chunks to create. Default is 1.}

\item{runwise}{If TRUE, the data chunks are created run-wise. Default is FALSE.}

\item{...}{Additional arguments.}
}
\value{
A list of data chunks, with each chunk containing the data, voxel indices, row indices, and chunk number.
}
\description{
This function creates data chunks for fmri_mem_dataset objects. It allows for the retrieval of run-wise or sequence-wise data chunks, as well as arbitrary chunks.
}
\examples{
# Create an fmri_mem_dataset
# ... (see example for fmri_mem_dataset)

# Create an iterator with 5 chunks
iter <- data_chunks(dset, nchunks=5)
`\%do\%` <- foreach::`\%do\%`
y <- foreach::foreach(chunk = iter) \%do\% { colMeans(chunk$data) }
length(y) == 5

# Create an iterator with 100 chunks
iter <- data_chunks(dset, nchunks=100)
y <- foreach::foreach(chunk = iter) \%do\% { colMeans(chunk$data) }
length(y) == 100

# Create a "runwise" iterator
iter <- data_chunks(dset, runwise=TRUE)
y <- foreach::foreach(chunk = iter) \%do\% { colMeans(chunk$data) }
length(y) == 2
}
