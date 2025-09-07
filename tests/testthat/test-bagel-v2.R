# Tests for BAGEL v2.0 Core Functions

test_that("validate_segments works correctly", {
  # Valid segments
  valid_segments <- data.frame(
    Sample = c("S1", "S1"),
    Chromosome = c(1, 2),
    Start = c(1000, 2000),
    End = c(1500, 2500),
    Log2Ratios = c(0.1, -0.2)
  )
  
  expect_silent(validate_segments(valid_segments))
  
  # Missing required columns
  invalid_segments <- data.frame(
    Sample = c("S1", "S1"),
    Chromosome = c(1, 2)
  )
  
  expect_error(validate_segments(invalid_segments), "Missing required columns")
  
  # Invalid coordinates
  invalid_coords <- data.frame(
    Sample = c("S1"),
    Chromosome = c(1),
    Start = c(2000),
    End = c(1000),  # End before Start
    Log2Ratios = c(0.1)
  )
  
  expect_error(validate_segments(invalid_coords), "Invalid coordinates")
})

test_that("get_breakpoint_data function works", {
  skip_if_not(exists("consensus_arm_definitions"), "Package data not loaded")
  
  # Test consensus data access
  consensus <- get_breakpoint_data("consensus", "arm_definitions")
  expect_s3_class(consensus, "data.frame")
  expect_true("arm" %in% names(consensus))
  expect_true("chr_num" %in% names(consensus))
})

test_that("validate_and_clean_segments removes invalid data", {
  segments <- data.frame(
    Sample = c("S1", "S1", "S1", "S1"),
    Chromosome = c(1, 23, 1, 1),  # Chr 23 should be removed
    Start = c(1000, 1000, 2000, 3000),
    End = c(1500, 1500, 1500, 3500),  # Invalid: Start > End
    Log2Ratios = c(0.1, 0.2, NA, 0.4)  # NA should be removed
  )
  
  cleaned <- validate_and_clean_segments(segments)
  
  expect_equal(nrow(cleaned), 1)  # Only one valid segment should remain
  expect_true(all(cleaned$Chromosome %in% 1:22))
  expect_true(all(!is.na(cleaned$Log2Ratios)))
  expect_true(all(cleaned$Start < cleaned$End))
})

test_that("fit_background_distribution works with different inputs", {
  segments <- data.frame(
    Sample = rep("S1", 100),
    Chromosome = rep(1, 100),
    Start = 1:100 * 1000,
    End = 1:100 * 1000 + 500,
    Log2Ratios = rnorm(100, mean = 0, sd = 0.2)
  )
  
  bg_dist <- fit_background_distribution(segments)
  
  expect_type(bg_dist, "list")
  expect_true("mu" %in% names(bg_dist))
  expect_true("sigma" %in% names(bg_dist))
  expect_true("type" %in% names(bg_dist))
  expect_equal(bg_dist$type, "gaussian")
})

test_that("convert_legacy_breakpoints works", {
  # Create mock legacy format
  legacy_breakpoints <- list(
    tel_bound = data.frame(
      arm = c("1p", "1q"),
      smallest_start = c(1, 125000000),
      largest_end = c(124999999, 249250621),
      direction = c("del", "amp")
    )
  )
  
  converted <- convert_legacy_breakpoints(legacy_breakpoints)
  
  expect_s3_class(converted, "data.frame")
  expect_true("arm" %in% names(converted))
  expect_true("chr_num" %in% names(converted))
  expect_true("arm_type" %in% names(converted))
  expect_equal(nrow(converted), 2)
})

test_that("bagel_log function works", {
  # Test basic logging
  expect_output(bagel_log("Test message", "INFO"), "INFO: Test message")
  expect_output(bagel_log("Warning message", "WARN"), "WARN: Warning message")
})

test_that("process_breakpoint_input handles different formats", {
  skip_if_not(exists("consensus_arm_definitions"), "Package data not loaded")
  
  # Test character input
  result1 <- process_breakpoint_input("consensus")
  expect_s3_class(result1, "data.frame")
  
  # Test data frame input
  df_input <- data.frame(
    arm = c("1p", "1q"),
    chr_num = c(1, 1),
    arm_type = c("p", "q"),
    arm_start = c(1, 125000000),
    arm_end = c(124999999, 249250621),
    direction = c("del", "amp")
  )
  
  result2 <- process_breakpoint_input(df_input)
  expect_s3_class(result2, "data.frame")
  expect_equal(nrow(result2), 2)
})

test_that("error handling works for invalid inputs", {
  # Test invalid cancer type
  expect_error(get_breakpoint_data("invalid_cancer_type"), "not found")
  
  # Test invalid breakpoint format
  expect_error(process_breakpoint_input(list(invalid = "data")), "Invalid breakpoints format")
})

# Integration test (only run if package data is available)
test_that("basic workflow integration works", {
  skip_if_not(exists("consensus_arm_definitions"), "Package data not loaded")
  
  # Create minimal test data
  test_segments <- data.frame(
    Sample = rep(c("S1", "S2"), each = 10),
    Chromosome = rep(1:2, 10),
    Start = rep(c(1000, 125000000) * 1:10, 2),
    End = rep(c(1000, 125000000) * 1:10 + 1000, 2),
    Log2Ratios = rnorm(20, 0, 0.3)
  )
  
  # Test the processing pipeline
  expect_no_error({
    results <- process_segments_pipeline(
      segments = test_segments,
      arm_definitions = consensus_arm_definitions,
      amp_threshold = 0.3,
      del_threshold = -0.3
    )
  })
  
  expect_type(results, "list")
  expect_true("segments" %in% names(results))
  expect_true("arm_summaries" %in% names(results))
})