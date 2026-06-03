test_that("explicit control map resolves channel order from machine channels", {
  sample_names <- c(
    "1-F2-empty_vector.fcs",
    "1-F3-FL5-A.fcs",
    "1-F4-FL8-A.fcs",
    "1-F5-FL11-A.fcs",
    "1-F6-FL1-A.fcs"
  )
  machine_channels <- c("FL5-A", "FL11-A", "FL1-A", "FL8-A")
  controls_map <- list(
    unstained = 1,
    channels = list(
      "FL5-A" = 2,
      "FL11-A" = 4,
      "FL1-A" = 5,
      "FL8-A" = 3
    )
  )

  resolved <- resolve_control_indices(sample_names, machine_channels, controls_map)

  expect_equal(resolved$controls_index, c(1, 2, 4, 5, 3))
  expect_equal(resolved$mapping$channel, c(NA, machine_channels))
  expect_equal(resolved$mapping$sample_name, sample_names[c(1, 2, 4, 5, 3)])
})

test_that("explicit control map supports exact names and sample patterns", {
  sample_names <- c("Unstained.fcs", "FITC_single.fcs", "APC_single.fcs")
  controls_map <- list(
    unstained = "Unstained.fcs",
    channels = list(
      "FL1-A" = "FITC",
      "FL3-A" = "APC"
    )
  )

  resolved <- resolve_control_indices(sample_names, c("FL1-A", "FL3-A"),
                                      controls_map)

  expect_equal(resolved$controls_index, c(1, 2, 3))
})

test_that("explicit control map fails before compensation when ambiguous", {
  sample_names <- c("Unstained.fcs", "FITC_A.fcs", "FITC_B.fcs")
  controls_map <- list(
    unstained = "Unstained",
    channels = list("FL1-A" = "FITC")
  )

  expect_error(
    resolve_control_indices(sample_names, "FL1-A", controls_map),
    "matched more than one detected control"
  )
})

test_that("explicit control map requires every selected compensation channel", {
  sample_names <- c("Unstained.fcs", "FITC_A.fcs")
  controls_map <- list(
    unstained = 1,
    channels = list("FL1-A" = 2)
  )

  expect_error(
    resolve_control_indices(sample_names, c("FL1-A", "FL3-A"), controls_map),
    "missing compensation controls for: FL3-A"
  )
})

test_that("compensation channels follow the current channel_pattern", {
  channels <- c("FSC-A", "SSC-A", "FL5-A", "FL5-H", "FL11-A", "FITC-A")

  expect_equal(compensation_channel_names(channels, "-A"),
               c("FL5-A", "FL11-A", "FITC-A"))
  expect_equal(compensation_channel_names(channels, "-H"), "FL5-H")
})

test_that("explicit control map fails when no compensation channels match", {
  expect_error(
    resolve_control_indices(
      sample_names = c("Unstained.fcs", "FITC_A.fcs"),
      channels = character(),
      controls_map = list(
        unstained = 1,
        channels = list("FL1-A" = 2)
      )
    ),
    "No compensation channels matched"
  )
})

test_that("design csv identifies controls by well and resolves machine order", {
  sample_names <- c(
    "01-Controls-F2.fcs",
    "01-Controls-F3.fcs",
    "01-Controls-F4.fcs",
    "01-Controls-F5.fcs",
    "01-Controls-F6.fcs",
    "01-Sample-G1.fcs"
  )
  design <- data.table::data.table(
    ID = c("1-F2", "1-F3", "1-F4", "1-F5", "1-F6", "1-G1"),
    Channel = c("unstained", "APC", "Pacific Blue", "ECD", "FITC", "FITC")
  )
  control_design <- list(
    design = design,
    sample_id_col = "ID",
    channel_col = "Channel"
  )
  channel_metadata <- data.table::data.table(
    channel = c("FL5-A", "FL11-A", "FL1-A", "FL8-A"),
    desc = c("APC-A", "mRuby3 ECD-A", "mCitrine FITC-A",
             "Pacific Blue-A")
  )

  control_index <- grepl("Controls", sample_names)
  identified <- identify_control_samples_from_design(sample_names,
                                                     control_design,
                                                     control_index)
  expect_equal(which(identified$control_index), 1:5)

  resolved <- resolve_control_indices_from_design(
    sample_names[identified$control_index],
    channels = c("FL5-A", "FL11-A", "FL1-A", "FL8-A"),
    matched_design = identified$design,
    sample_id_col = "ID",
    channel_metadata = channel_metadata
  )

  expect_equal(resolved$controls_index, c(1, 2, 4, 5, 3))
  expect_equal(resolved$mapping$sample_name,
               sample_names[c(1, 2, 4, 5, 3)])
  expect_equal(resolved$mapping$control_value,
               c("unstained", "APC", "ECD", "FITC", "Pacific Blue"))
})

test_that("design csv can use plain well IDs", {
  sample_names <- c("Sample_A1.fcs", "Sample_A2.fcs", "Sample_A3.fcs")
  design <- data.table::data.table(
    well = c("A1", "A2", "A3"),
    control_channel = c("unstained", "FL1-A", "FL5-A")
  )
  control_design <- list(
    design = design,
    sample_id_col = "well",
    channel_col = "control_channel"
  )

  identified <- identify_control_samples_from_design(sample_names, control_design)
  resolved <- resolve_control_indices_from_design(
    sample_names,
    channels = c("FL5-A", "FL1-A"),
    matched_design = identified$design,
    sample_id_col = "well"
  )

  expect_equal(which(identified$control_index), 1:3)
  expect_equal(resolved$controls_index, c(1, 3, 2))
})

test_that("design csv control rows are selected by controls_pattern matches", {
  sample_names <- c("Controls_A1.fcs", "Controls_A2.fcs", "Sample_A3.fcs")
  design <- data.table::data.table(
    well = c("A1", "A2", "A3"),
    Channel = c("unstained", "FITC", "FITC")
  )
  control_design <- list(
    design = design,
    sample_id_col = "well",
    channel_col = "Channel"
  )

  identified <- identify_control_samples_from_design(
    sample_names,
    control_design,
    control_index = grepl("Controls", sample_names)
  )

  expect_equal(which(identified$control_index), 1:2)
  expect_equal(identified$design$well, c("A1", "A2"))
})

test_that("control design can be filtered by gated population", {
  design <- data.table::data.table(
    well = c("A1", "A2", "B1", "B2", "C1", "C2"),
    Channel = c("unstained", "FITC", "unstained", "FITC",
                "unstained", "FITC"),
    gated_population = c("HEK", "HEK", "Tcells", "Tcells", "", "")
  )
  control_design <- list(
    design = design,
    sample_id_col = "well",
    channel_col = "Channel",
    population_col = "gated_population"
  )

  hek_design <- filter_control_design_for_population(control_design, "HEK")
  tcells_design <- filter_control_design_for_population(control_design, "Tcells")
  global_design <- filter_control_design_for_population(control_design, "Bcells")

  expect_equal(hek_design$design$well, c("A1", "A2"))
  expect_equal(tcells_design$design$well, c("B1", "B2"))
  expect_equal(global_design$design$well, c("C1", "C2"))
})

test_that("control population column is optional for legacy designs", {
  design <- data.table::data.table(
    well = c("A1", "A2"),
    Channel = c("unstained", "FITC")
  )
  control_design <- list(
    design = design,
    sample_id_col = "well",
    channel_col = "Channel",
    population_col = "gated_population"
  )

  filtered <- filter_control_design_for_population(control_design, "HEK")

  expect_equal(filtered$design$well, design$well)
})

test_that("FCS channel metadata aliases fluorophores to detector channels", {
  metadata <- data.table::data.table(
    channel = c("FL1-A", "FL11-A"),
    desc = c("mCitrine FITC-A", "mRuby3 ECD-A")
  )
  aliases <- channel_alias_table(metadata)

  expect_equal(resolve_control_value_to_channel("FITC",
                                                c("FL1-A", "FL11-A"),
                                                aliases),
               "FL1-A")
  expect_equal(resolve_control_value_to_channel("ECD",
                                                c("FL1-A", "FL11-A"),
                                                aliases),
               "FL11-A")
  expect_equal(resolve_control_value_to_channel("FL1-A",
                                                c("FL1-A", "FL11-A"),
                                                aliases),
               "FL1-A")
})
