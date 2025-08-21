context("runwise_lm rank deficiency")

options(mc.cores=1)

library(testthat)
library(fmrireg)

set.seed(123)

sframe <- fmrihrf::sampling_frame(blocklens = c(10,10), TR = 1)

# simple event table with one event per run
etab <- data.frame(onset = c(1,1),
                   condition = factor(c("A","A")),
                   run = c(1,2))

Y <- matrix(rnorm(sum(fmrihrf::blocklens(sframe)) * 2),
            sum(fmrihrf::blocklens(sframe)), 2)

dset <- matrix_dataset(Y, TR = 1,
                       run_length = fmrihrf::blocklens(sframe),
                       event_table = etab)

espec <- event_model(onset ~ hrf(condition), data = etab,
                     block = ~run, sampling_frame = sframe)

# nuisance columns duplicate the runwise intercept
nlist <- list(matrix(1,10,1), matrix(1,10,1))
bspec <- baseline_model(basis = "poly", degree = 1,
                        sframe = sframe, intercept = "runwise",
                        nuisance_list = nlist)

fmod <- fmri_model(espec, bspec, dset)

expect_warning(
  fit <- fmri_lm(onset ~ hrf(condition), block = ~run,
                 dataset = dset, baseline_model = bspec,
                 strategy = "runwise", use_fast_path = TRUE),
  regexp = "rank deficient"
)

# cov.unscaled should keep column names
expect_equal(colnames(fit$result$cov.unscaled), colnames(design_matrix(fmod)))
