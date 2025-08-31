# fmri_ttest demonstration script
# Shows AFNI 3dttest++ style workflows in fmrireg

library(fmrireg)

# Simulate some example data
set.seed(42)
n_subjects <- 20
n_voxels <- 100

# Example 1: One-sample t-test
# Testing if group mean differs from zero
cat("=== One-sample t-test ===\n")
Y_onesample <- matrix(rnorm(n_subjects * n_voxels, mean = 0.5), 
                      nrow = n_subjects, ncol = n_voxels)

gd_onesample <- list(
  blocks = list(
    list(
      Y = Y_onesample,
      V = NULL,
      covars = data.frame(subject = 1:n_subjects),
      feature = NULL
    )
  )
)
class(gd_onesample) <- "group_data"

fit1 <- fmri_ttest(gd_onesample, formula = ~ 1, engine = "classic")
cat("Mean t-value:", mean(fit1$t[1,]), "\n")
cat("Proportion significant (p < 0.05):", mean(fit1$p[1,] < 0.05), "\n\n")

# Example 2: Two-sample t-test
# Testing group differences
cat("=== Two-sample t-test ===\n")
group <- factor(rep(c("Control", "Treatment"), each = 10))
Y_twosample <- matrix(rnorm(n_subjects * n_voxels), nrow = n_subjects, ncol = n_voxels)
Y_twosample[group == "Treatment", ] <- Y_twosample[group == "Treatment", ] + 0.8

gd_twosample <- list(
  blocks = list(
    list(
      Y = Y_twosample,
      V = NULL,
      covars = data.frame(subject = 1:n_subjects, group = group),
      feature = NULL
    )
  )
)
class(gd_twosample) <- "group_data"

fit2 <- fmri_ttest(gd_twosample, formula = ~ 1 + group, engine = "classic")
cat("Group effect t-value (mean):", mean(fit2$t["groupTreatment",]), "\n")
cat("Proportion significant:", mean(fit2$p["groupTreatment",] < 0.05), "\n\n")

# Example 3: Paired t-test (via differences)
# Simulating paired A/B conditions
cat("=== Paired t-test ===\n")
Y_A <- matrix(rnorm(n_subjects * n_voxels, mean = 0.2), nrow = n_subjects, ncol = n_voxels)
Y_B <- matrix(rnorm(n_subjects * n_voxels, mean = -0.1), nrow = n_subjects, ncol = n_voxels)

blkA <- list(
  Y = Y_A,
  V = matrix(runif(n_subjects * n_voxels, 0.1, 0.3), nrow = n_subjects, ncol = n_voxels),
  meta = list(subjects = paste0("sub", 1:n_subjects), contrast = "TaskA"),
  covars = data.frame(subject = 1:n_subjects)
)

blkB <- list(
  Y = Y_B,
  V = matrix(runif(n_subjects * n_voxels, 0.1, 0.3), nrow = n_subjects, ncol = n_voxels),
  meta = list(subjects = paste0("sub", 1:n_subjects), contrast = "TaskB"),
  covars = data.frame(subject = 1:n_subjects)
)

# Compute paired difference
diff_blk <- paired_diff_block(blkA, blkB, rho = 0.3)
cat("Mean difference (A - B):", mean(diff_blk$Y), "\n")

# Test if difference differs from zero
gd_paired <- list(
  blocks = list(diff_blk)
)
class(gd_paired) <- "group_data"

fit3 <- fmri_ttest(gd_paired, formula = ~ 1, engine = "meta")
cat("Paired difference z-value (mean):", mean(fit3$z[1,]), "\n\n")

# Example 4: ANCOVA with continuous covariate
cat("=== ANCOVA with age covariate ===\n")
age <- rnorm(n_subjects, mean = 40, sd = 10)
Y_ancova <- matrix(rnorm(n_subjects * n_voxels), nrow = n_subjects, ncol = n_voxels)
# Add age effect
Y_ancova <- Y_ancova + outer(age - mean(age), rep(0.05, n_voxels))

gd_ancova <- list(
  blocks = list(
    list(
      Y = Y_ancova,
      V = NULL,
      covars = data.frame(subject = 1:n_subjects, group = group, age = age),
      feature = NULL
    )
  )
)
class(gd_ancova) <- "group_data"

fit4 <- fmri_ttest(gd_ancova, formula = ~ 1 + group + age, engine = "classic")
cat("Age effect t-value (mean):", mean(fit4$t["age",]), "\n")
cat("Group effect (adjusted for age):", mean(fit4$t["groupTreatment",]), "\n\n")

# Example 5: Welch t-test for unequal variances
cat("=== Welch t-test ===\n")
group_welch <- factor(rep(c("Low_Var", "High_Var"), each = 10))
Y_welch <- matrix(0, nrow = n_subjects, ncol = n_voxels)
Y_welch[group_welch == "Low_Var", ] <- matrix(rnorm(10 * n_voxels, mean = 0, sd = 1), nrow = 10)
Y_welch[group_welch == "High_Var", ] <- matrix(rnorm(10 * n_voxels, mean = 0.5, sd = 3), nrow = 10)

gd_welch <- list(
  blocks = list(
    list(
      Y = Y_welch,
      V = NULL,
      covars = data.frame(subject = 1:n_subjects, group = group_welch),
      feature = NULL
    )
  )
)
class(gd_welch) <- "group_data"

fit5_classic <- fmri_ttest(gd_welch, formula = ~ 1 + group, engine = "classic")
fit5_welch <- fmri_ttest(gd_welch, formula = ~ 1 + group, engine = "welch")

# Get the group coefficient row (second row after intercept)
cat("Classic t-test df:", unique(fit5_classic$df[2,]), "\n")
cat("Welch t-test df (mean):", mean(fit5_welch$df[2,]), "\n")
cat("Note: Welch df < classic df due to variance heterogeneity\n\n")

# Example 6: Multiple comparisons with spatial FDR
cat("=== Multiple comparisons correction ===\n")
fit6_uncorr <- fmri_ttest(gd_twosample, formula = ~ 1 + group, engine = "classic", mc = NULL)
fit6_bh <- fmri_ttest(gd_twosample, formula = ~ 1 + group, engine = "classic", mc = "bh")

cat("Uncorrected p < 0.05:", sum(fit6_uncorr$p[2,] < 0.05), "voxels\n")
cat("BH-corrected q < 0.05:", sum(fit6_bh$q[2,] < 0.05), "voxels\n")
cat("FDR control reduces false positives\n")