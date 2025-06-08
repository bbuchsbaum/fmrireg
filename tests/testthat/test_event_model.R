library(testthat)
testthat::local_edition(3)
#library(fmrireg)
library(rlang)

# Helper function to create consistent test data
create_test_data <- function(seed = 123, n_events = 10, n_blocks = 2) {
  set.seed(seed)
  block_len <- ceiling(n_events / n_blocks)
  blocklens <- rep(block_len, n_blocks)
  total_scans <- sum(blocklens)
  TR <- 2
  sf <- fmrihrf::sampling_frame(blocklens = blocklens, TR = TR)
  
  events <- data.frame(
    x = rnorm(n_events),
    y = rnorm(n_events),
    condition = factor(rep(letters[1:2], length.out=n_events)),
    modulator = seq(0, 1, length.out=n_events),
    # Generate realistic onsets within blocks
    onset = unlist(lapply(blocklens, function(bl) sort(sample(0:(bl*TR - TR), size=n_events/n_blocks, replace=FALSE)))),
    block = rep(1:n_blocks, each=n_events/n_blocks)
  )
  list(events=events, sf=sf)
}

# === Test Suite for the *new* event_model pipeline ===

test_that("event_model (formula) creates correct object structure and attributes", {
  test_data <- create_test_data()
  events <- test_data$events
  sf <- test_data$sf
  
  # Define formula with various terms
  form <- onset ~ hrf(condition)
  
  model <- event_model(
    formula_or_list = form,
    data = events,
    block = ~block, 
    sampling_frame = sf
  )
  
  # 1. Basic Structure Checks
  expect_s3_class(model, "event_model")
  expect_s3_class(model, "list")
  expect_named(model, c("terms", "design_matrix", "blockids", "sampling_frame", "contrasts", "model_spec"))
  expect_type(model$terms, "list")
  expect_s3_class(model$design_matrix, "tbl_df")
  expect_s3_class(model$sampling_frame, "sampling_frame")
  expect_true(is.numeric(model$blockids) || is.integer(model$blockids))
  expect_length(model$blockids, nrow(events)) # Check blockids length against original data
  
  # Check number of terms generated
  expect_equal(length(model$terms), 1)
  expect_s3_class(model$terms[[1]], "event_term") # Check class of a term
  
  # Check design matrix dimensions (rows should match total scans in sampling frame)
  expect_equal(nrow(model$design_matrix), sum(fmrihrf::blocklens(sf)))
  
  # 2. Design Matrix Attribute Checks
  dm <- model$design_matrix
  term_spans <- attr(dm, "term_spans")
  col_indices <- attr(dm, "col_indices")
  
  expect_true(!is.null(term_spans))
  expect_true(!is.null(col_indices))
  
  expect_true(is.vector(term_spans), "term_spans should be a vector")
  expect_equal(length(term_spans), length(model$terms))
  if (length(model$terms) > 0) {
    expect_equal(term_spans, cumsum(sapply(model$terms, function(term) ncol(design_matrix(term)))))
  }
  expect_type(col_indices, "list")
  expect_equal(length(col_indices), length(model$terms))
  expect_equal(names(col_indices), names(model$terms))
  expect_type(col_indices[[1]], "integer") # Indices should be integers
  # Check total number of columns matches sum of indices lengths
  expect_equal(sum(sapply(col_indices, length)), ncol(dm))
  # Check indices are within bounds
  expect_true(all(unlist(col_indices) >= 1 & unlist(col_indices) <= ncol(dm)))
})

test_that("event_model (list) creates correct object structure and attributes", {
  test_data <- create_test_data()
  events <- test_data$events
  sf <- test_data$sf
  
  # Define list of hrfspec objects
  spec_list <- list(
    hrf(condition), 
    hrf(modulator),
    hrf(x, y),
    hrf(Poly(x, 2))
  )
  
  model <- event_model(
    formula_or_list = spec_list,
    data = events,
    block = events$block, # Use vector block specification
    sampling_frame = sf
  )
  
  # Repeat the same structure and attribute checks as the formula test
  # 1. Basic Structure Checks
  expect_s3_class(model, "event_model")
  expect_named(model, c("terms", "design_matrix", "blockids", "sampling_frame", "contrasts", "model_spec"))
  expect_type(model$terms, "list")
  expect_equal(length(model$terms), length(spec_list))
  expect_s3_class(model$design_matrix, "tbl_df")
  expect_equal(nrow(model$design_matrix), sum(fmrihrf::blocklens(sf)))
  
  # 2. Design Matrix Attribute Checks
  dm <- model$design_matrix
  term_spans <- attr(dm, "term_spans")
  col_indices <- attr(dm, "col_indices")
  expect_true(!is.null(term_spans))
  expect_true(!is.null(col_indices))
  expect_true(inherits(term_spans, "matrix"))
  expect_equal(dim(term_spans), c(length(model$terms), 2))
  expect_type(col_indices, "list")
  expect_equal(length(col_indices), length(model$terms))
  expect_equal(names(col_indices), names(model$terms))
  expect_equal(sum(sapply(col_indices, length)), ncol(dm))
})

test_that("contrast_weights.event_model works with new attributes", {
  test_data <- create_test_data()
  events <- test_data$events
  sf <- test_data$sf
  
  # Define a contrast within hrfspec
  con1 <- pair_contrast(~ condition == "a", ~ condition == "b", name = "A_vs_B")
  
  # Model with one term containing a contrast
  model <- event_model(
    onset ~ hrf(condition, contrasts = list(con1)), 
    data = events, 
    block = ~block, 
    sampling_frame = sf
  )
  
  # Calculate contrast weights using the S3 method for event_model
  cw_list <- contrast_weights.event_model(model)
  
  # Checks
  expect_type(cw_list, "list")
  expect_equal(length(cw_list), 1) # Should find one contrast
  expect_named(cw_list, "condition.A_vs_B")
  
  cw <- cw_list[[1]]
  expect_s3_class(cw, "cell_contrast")
  #expect_named(cw, c("weights", "name", "term", "level1", "level2", "where", "offset_weights", "condnames"))
  
  # Check dimensions of offset_weights (should match full design matrix)
  dm <- design_matrix(model)
  expect_equal(nrow(cw$offset_weights), ncol(dm))
  expect_equal(ncol(cw$offset_weights), 1)
  expect_equal(rownames(cw$offset_weights), colnames(dm))
  
  # Check that weights are non-zero only in the columns for the relevant term
  term_indices <- attr(dm, "col_indices")[[names(model$terms)]] # Get indices for the 'condition' term
  expect_true(all(cw$offset_weights[-term_indices, ] == 0))
  expect_true(any(cw$offset_weights[term_indices, ] != 0))
  
  # Check the actual weight values (simple 1 vs -1 for pair_contrast)
  # The local weights should be [1, -1]
  expect_equal(as.vector(cw$offset_weights[term_indices, ]), c(1, -1))
})

# --- Test for B-Spline Basis --- 
test_that("event_model handles hrf_bspline basis correctly", {
  # Sampling frame with 2 blocks, 100 scans each
  sf_bspline <- fmrihrf::sampling_frame(blocklens = c(100, 100), TR = 2)
  
  # Events data spanning both blocks
  n_events_bspline <- 20
  events_bspline <- data.frame(
    onset = unlist(lapply(c(100, 100), function(bl) sort(sample(0:(bl*2 - 2), size=n_events_bspline/2, replace=FALSE)))),
    condition = factor(rep(letters[3:4], length.out=n_events_bspline)), # Use c, d
    block = rep(1:2, each=n_events_bspline/2)
  )
  
  # Define formula with bspline basis
  form_bspline <- onset ~ hrf(condition, basis="bspline", nbasis=5)
  
  # Create the event model
  model_bspline <- event_model(
    formula_or_list = form_bspline,
    data = events_bspline,
    block = ~block, 
    sampling_frame = sf_bspline
  )
  
  # Basic structural checks
  expect_s3_class(model_bspline, "event_model")
  expect_length(model_bspline$terms, 1)
  term_tag <- names(model_bspline$terms)[[1]] # Get the assigned term tag
  expect_equal(term_tag, "condition") # Check default tag
  
  # Check HRF specification within the term
  hrfspec <- attr(model_bspline$terms[[1]], "hrfspec")
  expect_s3_class(hrfspec$hrf, "HRF")
  expect_equal(fmrihrf::nbasis(hrfspec$hrf), 5, info = "Checking nbasis of the HRF object")
  expect_match(attr(hrfspec$hrf, "name"), "bspline", info = "Checking name of the HRF object")
  
  # Check design matrix dimensions
  dm_bspline <- design_matrix(model_bspline)
  expect_equal(nrow(dm_bspline), 200, info = "Checking total rows in design matrix")
  expect_equal(ncol(dm_bspline), 10, info = "Checking number of columns in design matrix")
  
  # Check attributes are present
  expect_true(!is.null(attr(dm_bspline, "term_spans")))
  expect_true(!is.null(attr(dm_bspline, "col_indices")))
  
  # Check term span corresponds to ncol
  term_spans <- attr(dm_bspline, "term_spans")
  expect_true(is.numeric(term_spans))
  expect_equal(length(term_spans), 1)
  expect_equal(as.vector(term_spans[1]), 10)
  
  # Check bspline basis column naming according to new scheme
  expected_base_conds <- c("condition.c", "condition.d")
  expected_names <- sort(c(
    paste0(term_tag, "_", expected_base_conds[1], "_b0", 1:5),
    paste0(term_tag, "_", expected_base_conds[2], "_b0", 1:5)
  ))
  expect_equal(sort(colnames(dm_bspline)), expected_names)
  
  # Test column_contrast for specific basis function pattern (should still work with new names)
  term_bspline <- model_bspline$terms[[1]]
  cc3 <- column_contrast(pattern_A = "_b03$", name = "Basis3")
  cw3 <- contrast_weights(cc3, term_bspline)
  all_rel <- conditions(term_bspline, drop.empty = FALSE, expand_basis = TRUE)
  expect_equal(rownames(cw3$weights), all_rel)
  mask <- grepl("_b03$", all_rel)
  expect_equal(sum(mask), length(unique(events_bspline$condition)), ignore_attr = TRUE)
  expect_equal(cw3$weights[mask, 1], rep(1 / sum(mask), sum(mask)), ignore_attr = TRUE)
  expect_equal(cw3$weights[!mask, 1], rep(0, length(all_rel) - sum(mask)), ignore_attr = TRUE)
  
  # Test pair_contrast broadcasting across bspline basis
  con1 <- pair_contrast(~ condition == "c", ~ condition == "d", name = "C_vs_D")
  model_contrast_bspline2 <- event_model(
    onset ~ hrf(condition, basis = "bspline", nbasis = 5, contrasts = list(con1)),
    data = events_bspline,
    block = ~block,
    sampling_frame = sf_bspline
  )
  cw_list2 <- contrast_weights.event_model(model_contrast_bspline2)
  expect_named(cw_list2, "condition.C_vs_D") # Name comes from contrast + term
  cw2 <- cw_list2[[1]]
  rown <- rownames(cw2$offset_weights)
  # Grep for the new naming patterns
  pos_idx <- grep(paste0("^",term_tag,"_condition\\.c_b"), rown)
  neg_idx <- grep(paste0("^",term_tag,"_condition\\.d_b"), rown)
  expect_true(length(pos_idx) == 5) # Should match 5 basis functions
  expect_true(length(neg_idx) == 5) # Should match 5 basis functions
  expect_true(all(cw2$offset_weights[pos_idx, 1] > 0))
  expect_true(all(cw2$offset_weights[neg_idx, 1] < 0))
})

# --- Tests for Term Naming and Clash Resolution ---
test_that("event_model handles term naming and clashes correctly", {
  test_data <- create_test_data()
  events <- test_data$events
  sf <- test_data$sf

  # Scenario 1: Simple clash, no explicit names
  # Default name for hrf(condition) is "condition"
  model1 <- event_model(
    onset ~ hrf(condition) + hrf(condition),
    data = events, block = ~block, sampling_frame = sf
  )
  # make_term_tag generates "condition", then "condition#1"
  expect_named(model1$terms, c("condition", "condition#1"))

  # Scenario 2: Clash with one explicit name identical to default
  # hrf(condition) -> default "condition"
  # hrf(condition, name="condition") -> explicit "condition"
  model2 <- event_model(
    onset ~ hrf(condition) + hrf(condition, name="condition"),
    data = events, block = ~block, sampling_frame = sf
  )
  # make_term_tag handles the clash: "condition", "condition#1"
  expect_named(model2$terms, c("condition", "condition#1"))

  # Scenario 3: No clash, one explicit name
  # hrf(condition, name="my_cond") -> explicit "my_cond"
  # hrf(modulator) -> default "modulator"
  model3 <- event_model(
    onset ~ hrf(condition, name="my_cond") + hrf(modulator),
    data = events, block = ~block, sampling_frame = sf
  )
  # No clash, names are kept as generated
  expect_named(model3$terms, c("my_cond", "modulator"))

  # Scenario 4: Clash between two explicit names
  # hrf(condition, name="term") -> explicit "term"
  # hrf(modulator, name="term") -> explicit "term"
  model4 <- event_model(
    onset ~ hrf(condition, name="term") + hrf(modulator, name="term"),
    data = events, block = ~block, sampling_frame = sf
  )
  # make_term_tag resolves the clash: "term", "term#1"
  expect_named(model4$terms, c("term", "term#1"))
  
  # Scenario 5: Multiple auto clashes mixed with unique names
  # hrf(condition) -> default "condition"
  # hrf(modulator) -> default "modulator"
  # hrf(condition) -> default "condition" (clash 1)
  # hrf(x, y) -> default "x_y"
  model5 <- event_model(
    onset ~ hrf(condition) + hrf(modulator) + hrf(condition) + hrf(x, y),
    data = events, block = ~block, sampling_frame = sf
  )
  # Expected sequence: "condition", "modulator", "condition#1", "x_y"
  expect_named(model5$terms, c("condition", "modulator", "condition#1", "x_y"))
  
  # Scenario 6: Using a list input with clashes
  # hrf(condition) -> default "condition"
  # hrf(modulator) -> default "modulator"
  # hrf(condition) -> default "condition" (clash 1)
  spec_list6 <- list(
    hrf(condition), 
    hrf(modulator),
    hrf(condition)
  )
  model6 <- event_model(spec_list6, data = events, block = ~block, sampling_frame = sf)
  # Expected sequence: "condition", "modulator", "condition#1"
  expect_named(model6$terms, c("condition", "modulator", "condition#1"))
  
  # Scenario 7: Clash where explicit name matches an auto-generated name 
  # hrf(condition) -> default "condition"
  # hrf(modulator, name="condition") -> explicit "condition" (clash 1)
  model7 <- event_model(
    onset ~ hrf(condition) + hrf(modulator, name="condition"),
    data = events, block = ~block, sampling_frame = sf
  )
  # Expected sequence: "condition", "condition#1"
  expect_named(model7$terms, c("condition", "condition#1"))
  
  # Scenario 8: Ensure subsetting doesn't affect default naming
  # hrf(condition, subset=cond_flag) -> default "condition"
  # hrf(condition, subset=!cond_flag) -> default "condition" (clash 1)
  events_subset <- events
  events_subset$cond_flag <- rep(c(TRUE, FALSE), length.out=nrow(events_subset))
  model8 <- event_model(
    onset ~ hrf(condition, subset=cond_flag) + hrf(condition, subset=!cond_flag),
    data = events_subset, block = ~block, sampling_frame = sf
  )
  # Expected sequence: "condition", "condition#1"
  expect_named(model8$terms, c("condition", "condition#1"))
  
  # Scenario 9: Explicit name that clashes with an index-disambiguated name
  # hrf(condition) -> default "condition"
  # hrf(condition) -> default "condition" -> clash 1 -> "condition#1"
  # hrf(modulator, name="condition") -> explicit "condition" -> clash 2 -> "condition#2"
  model9 <- event_model(
      onset ~ hrf(condition) + hrf(condition) + hrf(modulator, name="condition"),
      data = events, block = ~block, sampling_frame = sf
  )
  # Expected sequence: "condition", "condition#1", "condition#2" 
  expect_named(model9$terms, c("condition", "condition#1", "condition#2"))
  
  # Scenario 10: Name containing dot, sanitized
  # hrf(condition, name="term.one") -> explicit "term.one" -> sanitized "term_one"
  model10 <- event_model(
    onset ~ hrf(condition, name="term.one"),
    data = events, block = ~block, sampling_frame = sf
  )
  expect_named(model10$terms, c("term_one"))
  
  # Scenario 11: Clash after sanitization
  # hrf(condition, name="term.one") -> sanitized "term_one"
  # hrf(modulator, name="term_one") -> sanitized "term_one" -> clash -> "term_one#1"
  model11 <- event_model(
    onset ~ hrf(condition, name="term.one") + hrf(modulator, name="term_one"),
    data = events, block = ~block, sampling_frame = sf
  )
  expect_named(model11$terms, c("term_one", "term_one#1"))
})

# --- Golden Heading Snapshot Test ---
test_that("Golden heading test: diverse terms produce expected colnames", {
  td <- create_test_data(seed=555, n_events=10, n_blocks=1)
  events <- td$events
  sf     <- td$sf
  
  # Add more variables for diverse terms
  events$stim_type <- factor(rep(c("Face", "Scene"), length.out=10))
  events$group <- factor(rep(c("G1", "G2"), each=5))
  events$rt <- rnorm(10, mean=1.5, sd=0.2)
  
  # Create a model with various term types
  # Use explicit names to ensure stable term tags for snapshot
  model_golden <- event_model(
    onset ~ hrf(condition, name="Cond") + 
            hrf(condition, modulator, name="CondXMod") +
            hrf(stim_type, rt, basis="spmg3", name="StimRT_spmg3") + 
            hrf(Poly(x, 3), name="PolyX") +
            hrf(ScaleWithin(y, group), name="ScaleYbyG") +
            hrf(modulator, subset=condition=="a", name="Mod_CondA") + # Potential clash if auto-named
            hrf(modulator, subset=condition=="b", name="Mod_CondB"), 
    data = events,
    block = ~block,
    sampling_frame = sf
  )
  
  dm_golden <- design_matrix(model_golden)
  
  # Snapshot the column names
  expect_snapshot(colnames(dm_golden))
})

# --- Clash Test ---
test_that("Clash test: duplicate default term tags get unique suffixes", {
  td <- create_test_data(seed=1)
  events <- td$events
  sf <- td$sf
  
  # Create model where two terms would default to the same tag ("x")
  model_clash <- event_model(
    onset ~ hrf(x, subset = condition == "a") + 
            hrf(x, subset = condition == "b"), 
    data = events,
    block = ~block,
    sampling_frame = sf
  )
  
  dm_clash <- design_matrix(model_clash)
  cols <- colnames(dm_clash)
  
  # Check term tags assigned
  expect_named(model_clash$terms, c("x", "x#1"))
  
  # Check column name prefixes
  expect_true(all(startsWith(cols[1], "x_")))
  expect_true(all(startsWith(cols[2], "x.1_")))
  
  # Check specific column names (more precise)
  expect_equal(cols, c("x_x", "x.1_x"))
})

# --- Basis Suffix Test ---
test_that("Basis suffix test: _b## added correctly for multi-basis HRFs", {
  td <- create_test_data(seed=2)
  events <- td$events
  sf <- td$sf
  
  # HRF with 3 basis functions
  hrf_multi <- fmrihrf::HRF_SPMG3
  expect_equal(fmrihrf::nbasis(hrf_multi), 3)
  
  # HRF with 1 basis function
  hrf_single <- fmrihrf::HRF_SPMG1
  expect_equal(fmrihrf::nbasis(hrf_single), 1)
  
  # Model with multi-basis HRF
  model_multi <- event_model(
    onset ~ hrf(condition, modulator, basis=hrf_multi, name="Multi"), 
    data = events,
    block = ~block,
    sampling_frame = sf
  )
  dm_multi <- design_matrix(model_multi)
  cols_multi <- colnames(dm_multi)
  
  # Expected number of columns = nlevels(cond) * nbasis = 2 * 3 = 6
  expect_equal(ncol(dm_multi), 6)
  # All columns should end with _b01, _b02, or _b03
  expect_true(all(grepl("_b0[1-3]$", cols_multi)))
  # Check specific examples
  expect_true("Multi_condition.a_modulator_b01" %in% cols_multi)
  expect_true("Multi_condition.b_modulator_b03" %in% cols_multi)
  
  hrf_single <- fmrihrf::HRF_SPMG1
  # Model with single-basis HRF
  model_single <- event_model(
    onset ~ hrf(condition, modulator, basis=hrf_single, name="Single"), 
    data = events,
    block = ~block,
    sampling_frame = sf
  )
  dm_single <- design_matrix(model_single)
  cols_single <- colnames(dm_single)
  
  # Expected number of columns = nlevels(cond) * nbasis = 2 * 1 = 2
  expect_equal(ncol(dm_single), 2)
  # No columns should end with _b## suffix
  expect_false(any(grepl("_b[0-9]+$", cols_single)))
  # Check specific examples
  expect_equal(sort(cols_single), 
               sort(c("Single_condition.a_modulator", "Single_condition.b_modulator")))
})

# --- Interaction Test ---
test_that("Interaction test: columns named correctly for factor*factor and factor*poly", {
  td <- create_test_data(seed=3)
  events <- td$events
  sf <- td$sf
  events$facA <- factor(rep(c("A1", "A2"), each=5))
  events$facB <- factor(rep(c("B1", "B2", "B1", "B2", "B1"), 2))
  events$contX <- rnorm(10)
  
  # Factor x Factor
  model_ff <- event_model(
    onset ~ hrf(facA, facB, name="FF"), 
    data = events,
    block = ~block,
    sampling_frame = sf
  )
  dm_ff <- design_matrix(model_ff)
  cols_ff <- colnames(dm_ff)
  # Expected: FF_facA.A1_facB.B1, FF_facA.A2_facB.B1, FF_facA.A1_facB.B2, FF_facA.A2_facB.B2
  expect_equal(ncol(dm_ff), 4)
  expect_true(all(startsWith(cols_ff, "FF_facA.")))
  expect_setequal(cols_ff, c("FF_facA.A1_facB.B1", "FF_facA.A2_facB.B1", 
                           "FF_facA.A1_facB.B2", "FF_facA.A2_facB.B2"))
                           
  # Factor x Poly
  model_fp <- event_model(
    onset ~ hrf(facA, Poly(contX, 2), name="FP"), 
    data = events,
    block = ~block,
    sampling_frame = sf
  )
  dm_fp <- design_matrix(model_fp)
  cols_fp <- colnames(dm_fp)
  # Expected: FP_facA.A1_01, FP_facA.A2_01, 
  #           FP_facA.A1_02, FP_facA.A2_02
  expect_equal(ncol(dm_fp), 4)                 
  expect_true(all(startsWith(cols_fp, "FP_facA.")))
  expect_setequal(cols_fp, c("FP_facA.A1_01", "FP_facA.A2_01",
                           "FP_facA.A1_02", "FP_facA.A2_02"))

  # Factor x Poly with multi-basis HRF
  hrf_multi <- fmrihrf::HRF_SPMG3
  model_fp_basis <- event_model(
    onset ~ hrf(facA, Poly(contX, 2), basis=hrf_multi, name="FPB"), 
    data = events,
    block = ~block,
    sampling_frame = sf
  )
  dm_fpb <- design_matrix(model_fp_basis)
  cols_fpb <- colnames(dm_fpb)
  # Expected: 4 conditions * 3 basis = 12 columns
  expect_equal(ncol(dm_fpb), 12)
  expect_true(all(startsWith(cols_fpb, "FPB_facA.")))
  expect_true(all(grepl("_b0[1-3]$", cols_fpb))) # Check suffix
  expect_true("FPB_facA.A1_01_b01" %in% cols_fpb) # Corrected
  expect_true("FPB_facA.A2_02_b03" %in% cols_fpb) # Corrected
})

# --- Placeholder for future tests to be added ---

# (Original tests below are mostly for the deprecated create_event_model)
# It might be useful to adapt some scenarios (e.g., subsetting, basis functions) 
# to the new event_model interface later.

# test_that("create_event_model handles various variable specifications", { ... old code ... })
# test_that("create_event_model handles HRF specifications correctly", { ... old code ... })
# test_that("create_event_model handles subset specifications", { ... old code ... })
# test_that("create_event_model handles block-wise onsets correctly", { ... old code ... })
# test_that("event_model handles ~1 block formula correctly", { ... old code ... })

# Factorial design test: formula vs list interfaces
test_that("formula and list interfaces produce identical design matrices for 3x3 factorial design", {
  # Create events with two factors A (a,b,c) and B (x,y,z)
  n_events <- 9
  events <- data.frame(
    A = factor(rep(letters[1:3], each = 3)),
    B = factor(rep(letters[24:26], times = 3)),
    onset = seq_len(n_events),
    block = rep(1, n_events)
  )
  # Single block sampling frame
  sf <- fmrihrf::sampling_frame(blocklens = 9, TR = 1)

  # Formula-based model
  fm <- event_model(
    onset ~ hrf(A) + hrf(B),
    data = events,
    block = ~block,
    sampling_frame = sf
  )
  
  # List-based model with same hrfspecs
  specs <- list(hrf(A), hrf(B))
  lm <- event_model(
    formula_or_list = specs,
    data = events,
    block = ~block,
    sampling_frame = sf
  )
  
  # Terms and design matrices should match
  expect_named(fm$terms, names(lm$terms))
  expect_equal(colnames(fm$design_matrix), colnames(lm$design_matrix))
  expect_equal(fm$design_matrix, lm$design_matrix, ignore_attr = TRUE)
})

# 1. Factor × Continuous Interaction with Custom Durations
test_that("event_model handles factor × continuous interaction with custom durations", {
  td <- create_test_data(seed=42, n_events=12, n_blocks=3)
  events <- td$events
  sf     <- td$sf

  # add a random durations vector
  events$durations <- runif(nrow(events), 0.1, 2)

  # Interaction term: condition × modulator, custom durations
  model_ic <- event_model(
    onset ~ hrf(condition, modulator),
    data = events,
    block = ~block,
    sampling_frame = sf,
    durations = events$durations
  )
  dm_ic <- design_matrix(model_ic)

  # Expect two columns: one per level of 'condition' × modulator
  expect_equal(ncol(dm_ic), length(levels(events$condition)))
  

  # Confirm the durations argument didn't break row counts
  expect_equal(nrow(dm_ic), sum(fmrihrf::blocklens(sf)))
})

# 2. Fourier Basis on Continuous Modulator
test_that("event_model supports Fourier basis on continuous modulator", {
  td <- create_test_data(seed=7, n_events=8, n_blocks=2)
  events <- td$events; sf <- td$sf

  # Fourier basis of order 4
  model_fourier <- event_model(
    onset ~ hrf(modulator, basis = "fourier", nbasis = 4),
    data = events,
    block = ~block,
    sampling_frame = sf
  )
  dm_f <- design_matrix(model_fourier)
  term_tag <- names(model_fourier$terms)[[1]] # This is the name of the list element
  # The actual term_tag attribute on the event_term object is what matters for naming
  # but for a simple hrf(modulator), both should be "modulator".
  expect_equal(term_tag, "modulator")

  # Should have exactly 4 basis columns
  expect_equal(ncol(dm_f), 4)
  
  # Corrected expectation based on naming_refactoring.md:
  # term_tag = "modulator"
  # condition_tag for simple continuous variable "modulator" is "modulator"
  # HRF basis suffix = "_b##"
  # Final: modulator_modulator_b01, modulator_modulator_b02, ...
  expected_names <- paste0(term_tag, "_modulator", "_b", zeropad(1:4, 4))
  expect_setequal(colnames(dm_f), expected_names)
})

# 3. Mixed Bases via List Interface
test_that("event_model list interface supports mixed bases (Poly + BSpline + Std)", {
  td <- create_test_data(seed=101, n_events=15, n_blocks=3)
  events <- td$events; sf <- td$sf

  specs <- list(
    hrf(condition, name="COND"),
    hrf(modulator, basis = "bspline", nbasis = 4, name="MOD"),
    hrf(Poly(x, 2)), 
    hrf(Standardized(y)) 
  )
  model_mix <- event_model(
    formula_or_list = specs,
    data = events,
    block = ~block,
    sampling_frame = sf
  )
  dm_mix <- design_matrix(model_mix)
 
  # Check term tags based on new rules
  expected_term_tags <- c("COND", "MOD", "Poly_x", "Standardized_y")
  expect_named(model_mix$terms, expected_term_tags)
 
  # Check column names for each term
  actual_colnames <- colnames(dm_mix)

  # COND: COND_condition.a, COND_condition.b
  expect_true("COND_condition.a" %in% actual_colnames)
  expect_true("COND_condition.b" %in% actual_colnames)
  # MOD: MOD_modulator_b01 ... MOD_modulator_b04 (HRF bspline basis applied to modulator event)
  expect_true("MOD_modulator_b01" %in% actual_colnames)
  expect_true("MOD_modulator_b04" %in% actual_colnames)
  # Poly(x,2): Poly_x_01, Poly_x_02
  expect_true("Poly_x_01" %in% actual_colnames)
  expect_true("Poly_x_02" %in% actual_colnames)
  # Standardized(y): Standardized_y_y (Now returns argname 'y')
  expect_true("Standardized_y_y" %in% actual_colnames)
  
  # --- Check Total Columns --- (This calculation relies on the corrected columns/nbasis interaction)
  expected_cols_term1 <- length(levels(events$condition)) * 1 # nbasis=1 default for COND
  expected_cols_term2 <- 1 * 4 # modulator is 1 column (event), HRF nbasis=4 for MOD
  expected_cols_term3 <- 2 * 1 # Poly(x, 2) is 2 columns, nbasis=1 default
  expected_cols_term4 <- 1 * 1 # Standardized(y) is 1 column, nbasis=1 default
  expect_equal(ncol(dm_mix),
               expected_cols_term1 + expected_cols_term2 + expected_cols_term3 + expected_cols_term4)
})

# 4. Duration Recycling and Zero-Duration Events
test_that("event_model durations vector of length 1 is recycled", {
  td <- create_test_data(seed=99, n_events=6, n_blocks=2)
  events <- td$events; sf <- td$sf

  # Single scalar duration
  model_dur <- event_model(
    onset ~ hrf(condition),
    data = events,
    block = ~block,
    sampling_frame = sf,
    durations = 1.5
  )
  dm_dur <- design_matrix(model_dur)

  expect_equal(nrow(dm_dur), sum(fmrihrf::blocklens(sf)))
  # We know duration doesn't affect the column count
  expect_equal(ncol(dm_dur), length(levels(events$condition)))
})

# 5. Amplitude Modulation: Factor × Identity Basis
test_that("event_model supports amplitude modulation via identity basis", {
  td <- create_test_data(seed=202, n_events=10, n_blocks=2)
  events <- td$events; sf <- td$sf

  # identity basis of modulator: treat modulator as its own basis
  # For hrf(Ident(modulator)) with no id, term_tag is NULL.
  # The term name in model_amp$terms might be auto-generated if make_term_tag returns NULL.
  # We care about the final column name and the term_tag attribute.
  model_amp <- event_model(
    onset ~ hrf(Ident(modulator)), # Wrap in Ident
    data = events,
    block = ~block,
    sampling_frame = sf
  )
  dm_amp <- design_matrix(model_amp)
  
  # Check term_tag attribute on the event_term object itself
  term_obj_name <- names(model_amp$terms)[[1]] 
  term_tag_attr <- attr(model_amp$terms[[term_obj_name]], "term_tag")
  expect_null(term_tag_attr) # term_tag attribute should be NULL for Ident() only with no id
  # Single continuous variable Identity(modulator) -> 1 column
  expect_equal(ncol(dm_amp), 1)
  # Check final name: condition_tag (modulator) directly
  expect_equal(colnames(dm_amp), "modulator") # Corrected
})

test_that("event_model correctly names columns for Ident() and combined basis functions", {
  td <- create_test_data(n_events = 10, n_blocks = 1) # Simpler data for these tests
  events <- td$events
  sf <- td$sf
  events$RT1 <- rnorm(10)
  events$RT2 <- rnorm(10)
  events$Age <- rnorm(10)

  # --- Ident() Naming Tests ---
  # 1. Simple Ident
  model_ident_simple <- event_model(onset ~ hrf(Ident(RT1, RT2)), data = events, block = ~block, sampling_frame = sf)
  cols_ident_simple <- colnames(design_matrix(model_ident_simple))
  expect_setequal(cols_ident_simple, c("RT1", "RT2"))

  # 2. Ident with HRF basis
  model_ident_hrf <- event_model(onset ~ hrf(Ident(RT1, RT2), basis="spmg2"), data = events, block = ~block, sampling_frame = sf)
  cols_ident_hrf <- colnames(design_matrix(model_ident_hrf))
  expect_setequal(cols_ident_hrf, c("RT1_b01", "RT1_b02", "RT2_b01", "RT2_b02"))

  # 3. Ident with id
  model_ident_id <- event_model(onset ~ hrf(Ident(RT1, RT2), id="myvars"), data = events, block = ~block, sampling_frame = sf)
  cols_ident_id <- colnames(design_matrix(model_ident_id))
  expect_setequal(cols_ident_id, c("myvars_RT1", "myvars_RT2"))

  # 4. Ident with id and HRF basis
  model_ident_id_hrf <- event_model(onset ~ hrf(Ident(RT1, RT2), id="myvars", basis="spmg2"), data = events, block = ~block, sampling_frame = sf)
  cols_ident_id_hrf <- colnames(design_matrix(model_ident_id_hrf))
  expect_setequal(cols_ident_id_hrf, c("myvars_RT1_b01", "myvars_RT1_b02", "myvars_RT2_b01", "myvars_RT2_b02"))

  # --- Ident() Clash Test ---
  # Expect warning because 'RT1' will be generated by both hrf() terms initially
  expect_warning(
    model_ident_clash <- event_model(onset ~ hrf(Ident(RT1, Age)) + hrf(Ident(RT1, RT2)), 
                                   data = events, block = ~block, sampling_frame = sf),
    regexp = "Duplicate column names detected and automatically resolved"
  )
  cols_ident_clash <- colnames(design_matrix(model_ident_clash))
  # Order might matter here if make.names is predictable, but setequal is safer for now.
  # We expect RT1, Age, RT1.1 (or similar), RT2
  expect_true("RT1" %in% cols_ident_clash)
  expect_true("Age" %in% cols_ident_clash)
  expect_true(any(grepl("^RT1\\.\\d+$", cols_ident_clash))) # Checks for RT1.1, RT1.2 etc.
  expect_true("RT2" %in% cols_ident_clash)
  expect_equal(length(cols_ident_clash), 4) # Ensure no extra columns
  
  # --- Combined Basis Naming Tests (Poly + HRF) ---
  # Test Poly() without explicit HRF basis (should use default single HRF basis, no _b## suffix)
  model_poly_plain <- event_model(onset ~ hrf(Poly(Age, 2)), 
                                  data = events, block = ~block, sampling_frame = sf)
  cols_poly_plain <- colnames(design_matrix(model_poly_plain))
  # TermTag from Poly(Age,2) is Poly_Age. ConditionTags from columns.Poly are 01, 02.
  # Final: Poly_Age_01, Poly_Age_02
  expect_setequal(cols_poly_plain, c("Poly_Age_01", "Poly_Age_02"))
  
  # Default name for Poly(Age,2) -> term_tag = "Poly_Age"
  model_poly_hrf <- event_model(onset ~ hrf(Poly(Age, 2), basis="spmg2"), 
                                data = events, block = ~block, sampling_frame = sf)
  cols_poly_hrf <- colnames(design_matrix(model_poly_hrf))
  # Expected: Poly_Age_01_b01, Poly_Age_01_b02, Poly_Age_02_b01, Poly_Age_02_b02
  expected_poly_hrf_cols <- c(
    "Poly_Age_01_b01", "Poly_Age_01_b02", 
    "Poly_Age_02_b01", "Poly_Age_02_b02"
  )
  expect_setequal(cols_poly_hrf, expected_poly_hrf_cols)
  
  # Poly + HRF with explicit ID
  model_poly_id_hrf <- event_model(onset ~ hrf(Poly(Age, 2), basis="spmg2", id="age_poly"), 
                                   data = events, block = ~block, sampling_frame = sf)
  cols_poly_id_hrf <- colnames(design_matrix(model_poly_id_hrf))
  # Expected: age_poly_01_b01, age_poly_01_b02, age_poly_02_b01, age_poly_02_b02
  expected_poly_id_hrf_cols <- c(
    "age_poly_01_b01", "age_poly_01_b02", 
    "age_poly_02_b01", "age_poly_02_b02"
  )
  expect_setequal(cols_poly_id_hrf, expected_poly_id_hrf_cols)

  # --- Combined Basis Naming Tests (BSpline + HRF) ---
  # Test BSpline() without explicit HRF basis
  model_bs_plain <- event_model(onset ~ hrf(BSpline(Age, 3)), 
                                data = events, block = ~block, sampling_frame = sf)
  cols_bs_plain <- colnames(design_matrix(model_bs_plain))
  # Expected: BSpline_Age_01, BSpline_Age_02, BSpline_Age_03
  expect_setequal(cols_bs_plain, c("BSpline_Age_01", "BSpline_Age_02", "BSpline_Age_03"))
  
  # Default name for BSpline(Age,3) -> term_tag = "BSpline_Age"
  model_bs_hrf <- event_model(onset ~ hrf(BSpline(Age, 3), basis="spmg2"), 
                              data = events, block = ~block, sampling_frame = sf)
  cols_bs_hrf <- colnames(design_matrix(model_bs_hrf))
  # Expected: BSpline_Age_01_b01 ... BSpline_Age_03_b02
  expected_bs_hrf_cols <- c(
    "BSpline_Age_01_b01", "BSpline_Age_01_b02", 
    "BSpline_Age_02_b01", "BSpline_Age_02_b02",
    "BSpline_Age_03_b01", "BSpline_Age_03_b02"
  )
  expect_setequal(cols_bs_hrf, expected_bs_hrf_cols)
  
  # BSpline + HRF with explicit ID
  model_bs_id_hrf <- event_model(onset ~ hrf(BSpline(Age, 3), basis="spmg2", id="age_bs"), 
                                 data = events, block = ~block, sampling_frame = sf)
  cols_bs_id_hrf <- colnames(design_matrix(model_bs_id_hrf))
  # Expected: age_bs_01_b01 ... age_bs_03_b02
  expected_bs_id_hrf_cols <- c(
    "age_bs_01_b01", "age_bs_01_b02", 
    "age_bs_02_b01", "age_bs_02_b02",
    "age_bs_03_b01", "age_bs_03_b02"
  )
  expect_setequal(cols_bs_id_hrf, expected_bs_id_hrf_cols)

})


# New tests for prefix, subset, and row count errors

test_that("prefix works with multi-variable term", {
  td <- create_test_data(seed = 111)
  events <- td$events
  sf <- td$sf

  model_pre <- event_model(
    onset ~ hrf(condition, modulator, prefix = "pre"),
    data = events,
    block = ~block,
    sampling_frame = sf
  )
  cols <- colnames(design_matrix(model_pre))
  expect_setequal(cols, c("pre_condition.a_modulator", "pre_condition.b_modulator"))
})

test_that("subset expression filters events", {
  td <- create_test_data(seed = 222)
  events <- td$events
  sf <- td$sf
  events$flag <- rep(c(TRUE, FALSE), length.out = nrow(events))

  m_sub <- event_model(onset ~ hrf(condition, subset = flag),
                       data = events, block = ~block, sampling_frame = sf)
  m_manual <- event_model(onset ~ hrf(condition),
                          data = subset(events, flag),
                          block = ~block, sampling_frame = sf)

  expect_equal(design_matrix(m_sub), design_matrix(m_manual), ignore_attr = TRUE)
})

test_that("error when onset length mismatched with data", {
  td <- create_test_data(seed = 333)
  events <- td$events
  sf <- td$sf
  bad_onsets <- events$onset[-nrow(events)]
  expect_error(
    event_model(bad_onsets ~ hrf(condition), data = events,
                block = ~block, sampling_frame = sf),
    "Length of extracted onset variable"
  )
})

