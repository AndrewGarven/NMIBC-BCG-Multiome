# ============================================================
# 00_validate_metadata.R
#
# Validate sample-, patient-, and legacy alias metadata for the
# NMIBC BCG single-cell multiome reproducibility repository.
#
# This script uses base R only so that metadata validation
# does not depend on any external packages.
#
# Run from the repository root:
#
#   Rscript scripts/00_validate_metadata.R
#
# ============================================================


# ============================================================
# 1. FILE PATHS
# ============================================================

sample_metadata_path <- file.path(
  "metadata",
  "sample_metadata.csv"
)

patient_metadata_path <- file.path(
  "metadata",
  "patient_metadata.csv"
)

sample_aliases_path <- file.path(
  "metadata",
  "sample_aliases.csv"
)


# ============================================================
# 2. CHECK FILES EXIST
# ============================================================

required_files <- c(
  sample_metadata_path,
  patient_metadata_path,
  sample_aliases_path
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  stop(
    "Required metadata file(s) not found:\n",
    paste(
      paste0("  - ", missing_files),
      collapse = "\n"
    )
  )
}


# ============================================================
# 3. READ METADATA
# ============================================================

sample_meta <- read.csv(
  sample_metadata_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

patient_meta <- read.csv(
  patient_metadata_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_aliases <- read.csv(
  sample_aliases_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 4. HELPER FUNCTIONS
# ============================================================

check_required_columns <- function(
  data,
  required_columns,
  data_name
) {

  missing_columns <- setdiff(
    required_columns,
    colnames(data)
  )

  if (length(missing_columns) > 0) {

    stop(
      data_name,
      " is missing required column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
}


check_no_missing_values <- function(
  data,
  columns,
  data_name
) {

  for (column_name in columns) {

    missing_index <- is.na(data[[column_name]]) |
      trimws(
        as.character(
          data[[column_name]]
        )
      ) == ""

    if (any(missing_index)) {

      stop(
        data_name,
        " contains missing values in column '",
        column_name,
        "'."
      )
    }
  }
}


# ============================================================
# 5. SAMPLE METADATA VALIDATION
# ============================================================

required_sample_columns <- c(
  "sample_number",
  "sample_name",
  "sample_id",
  "patient_id",
  "timepoint",
  "bcg_week",
  "sample_context",
  "treatment",
  "outcome_raw",
  "analysis_outcome",
  "nuclei_batch_id"
)

check_required_columns(
  sample_meta,
  required_sample_columns,
  "sample_metadata.csv"
)

check_no_missing_values(
  sample_meta,
  required_sample_columns,
  "sample_metadata.csv"
)


# ------------------------------------------------------------
# Basic dimensions
# ------------------------------------------------------------

stopifnot(
  nrow(sample_meta) == 32
)

stopifnot(
  length(
    unique(sample_meta$sample_number)
  ) == 32
)

stopifnot(
  length(
    unique(sample_meta$sample_name)
  ) == 32
)

stopifnot(
  length(
    unique(sample_meta$sample_id)
  ) == 32
)

stopifnot(
  length(
    unique(sample_meta$patient_id)
  ) == 16
)


# ------------------------------------------------------------
# Sample numbering
# ------------------------------------------------------------

sample_numbers <- sort(
  as.integer(sample_meta$sample_number)
)

stopifnot(
  identical(
    sample_numbers,
    1:32
  )
)


# ------------------------------------------------------------
# Canonical sample names
# ------------------------------------------------------------

expected_sample_names <- sprintf(
  "multiome_%02d",
  1:32
)

stopifnot(
  identical(
    sort(sample_meta$sample_name),
    sort(expected_sample_names)
  )
)


# ------------------------------------------------------------
# Allowed categorical values
# ------------------------------------------------------------

stopifnot(
  all(
    sample_meta$timepoint %in%
      c("pre", "post")
  )
)

stopifnot(
  all(
    sample_meta$bcg_week %in%
      c("W1", "W6")
  )
)

stopifnot(
  all(
    sample_meta$sample_context %in%
      c(
        "pre-BCG",
        "pre-sixth_BCG"
      )
  )
)

stopifnot(
  all(
    sample_meta$treatment %in%
      c(
        "untreated",
        "BCG-treated"
      )
  )
)

stopifnot(
  all(
    sample_meta$outcome_raw %in%
      c(
        "remission",
        "recurrence"
      )
  )
)

stopifnot(
  all(
    sample_meta$analysis_outcome %in%
      c(
        "recurrence_free",
        "early_recurrence"
      )
  )
)


# ------------------------------------------------------------
# Timepoint consistency
# ------------------------------------------------------------

stopifnot(
  all(
    sample_meta$bcg_week[
      sample_meta$timepoint == "pre"
    ] == "W1"
  )
)

stopifnot(
  all(
    sample_meta$bcg_week[
      sample_meta$timepoint == "post"
    ] == "W6"
  )
)

stopifnot(
  all(
    sample_meta$sample_context[
      sample_meta$timepoint == "pre"
    ] == "pre-BCG"
  )
)

stopifnot(
  all(
    sample_meta$sample_context[
      sample_meta$timepoint == "post"
    ] == "pre-sixth_BCG"
  )
)

stopifnot(
  all(
    sample_meta$treatment[
      sample_meta$timepoint == "pre"
    ] == "untreated"
  )
)

stopifnot(
  all(
    sample_meta$treatment[
      sample_meta$timepoint == "post"
    ] == "BCG-treated"
  )
)


# ------------------------------------------------------------
# Outcome mapping consistency
# ------------------------------------------------------------

stopifnot(
  all(
    sample_meta$analysis_outcome[
      sample_meta$outcome_raw == "remission"
    ] == "recurrence_free"
  )
)

stopifnot(
  all(
    sample_meta$analysis_outcome[
      sample_meta$outcome_raw == "recurrence"
    ] == "early_recurrence"
  )
)


# ------------------------------------------------------------
# Each patient must have exactly two samples
# ------------------------------------------------------------

samples_per_patient <- table(
  sample_meta$patient_id
)

stopifnot(
  all(samples_per_patient == 2)
)


# ------------------------------------------------------------
# Every patient must have exactly one pre and one post sample
# ------------------------------------------------------------

patient_timepoints <- table(
  sample_meta$patient_id,
  sample_meta$timepoint
)

stopifnot(
  all(patient_timepoints[, "pre"] == 1)
)

stopifnot(
  all(patient_timepoints[, "post"] == 1)
)


# ------------------------------------------------------------
# Outcome must remain constant within each patient
# ------------------------------------------------------------

outcomes_per_patient <- tapply(
  sample_meta$analysis_outcome,
  sample_meta$patient_id,
  function(x) {
    length(unique(x))
  }
)

stopifnot(
  all(outcomes_per_patient == 1)
)


# ------------------------------------------------------------
# Raw outcome must remain constant within each patient
# ------------------------------------------------------------

raw_outcomes_per_patient <- tapply(
  sample_meta$outcome_raw,
  sample_meta$patient_id,
  function(x) {
    length(unique(x))
  }
)

stopifnot(
  all(raw_outcomes_per_patient == 1)
)


# ============================================================
# 6. PATIENT METADATA VALIDATION
# ============================================================

required_patient_columns <- c(
  "patient_id",
  "sex",
  "age_at_consent",
  "tumor_stage",
  "tumor_grade",
  "aua_risk",
  "outcome_raw",
  "analysis_outcome"
)

check_required_columns(
  patient_meta,
  required_patient_columns,
  "patient_metadata.csv"
)

check_no_missing_values(
  patient_meta,
  required_patient_columns,
  "patient_metadata.csv"
)


# ------------------------------------------------------------
# Basic patient dimensions
# ------------------------------------------------------------

stopifnot(
  nrow(patient_meta) == 16
)

stopifnot(
  length(
    unique(patient_meta$patient_id)
  ) == 16
)


# ------------------------------------------------------------
# Same patients must occur in sample and patient metadata
# ------------------------------------------------------------

sample_patients <- sort(
  unique(sample_meta$patient_id)
)

patient_patients <- sort(
  unique(patient_meta$patient_id)
)

stopifnot(
  identical(
    sample_patients,
    patient_patients
  )
)


# ------------------------------------------------------------
# Allowed clinical values
# ------------------------------------------------------------

stopifnot(
  all(
    patient_meta$sex %in%
      c("Male", "Female")
  )
)

stopifnot(
  all(
    patient_meta$tumor_stage %in%
      c("Ta", "T1")
  )
)

stopifnot(
  all(
    patient_meta$aua_risk %in%
      c(2, 3)
  )
)

stopifnot(
  all(
    patient_meta$outcome_raw %in%
      c(
        "remission",
        "recurrence"
      )
  )
)

stopifnot(
  all(
    patient_meta$analysis_outcome %in%
      c(
        "recurrence_free",
        "early_recurrence"
      )
  )
)


# ------------------------------------------------------------
# Patient-level outcome distribution
# ------------------------------------------------------------

patient_outcome_counts <- table(
  patient_meta$analysis_outcome
)

stopifnot(
  patient_outcome_counts[
    "recurrence_free"
  ] == 8
)

stopifnot(
  patient_outcome_counts[
    "early_recurrence"
  ] == 8
)


# ------------------------------------------------------------
# Clinical cohort consistency
# ------------------------------------------------------------

stopifnot(
  sum(
    patient_meta$sex == "Female"
  ) == 4
)

stopifnot(
  sum(
    patient_meta$sex == "Male"
  ) == 12
)

stopifnot(
  sum(
    patient_meta$tumor_stage == "T1"
  ) == 3
)

stopifnot(
  sum(
    patient_meta$tumor_stage == "Ta"
  ) == 13
)

stopifnot(
  sum(
    patient_meta$aua_risk == 2
  ) == 7
)

stopifnot(
  sum(
    patient_meta$aua_risk == 3
  ) == 9
)


# ------------------------------------------------------------
# Cross-check patient outcomes between sample and patient tables
# ------------------------------------------------------------

sample_patient_outcomes <- unique(
  sample_meta[
    ,
    c(
      "patient_id",
      "outcome_raw",
      "analysis_outcome"
    )
  ]
)

sample_patient_outcomes <- sample_patient_outcomes[
  order(
    sample_patient_outcomes$patient_id
  ),
]

patient_outcomes <- patient_meta[
  ,
  c(
    "patient_id",
    "outcome_raw",
    "analysis_outcome"
  )
]

patient_outcomes <- patient_outcomes[
  order(
    patient_outcomes$patient_id
  ),
]

rownames(sample_patient_outcomes) <- NULL
rownames(patient_outcomes) <- NULL

stopifnot(
  identical(
    sample_patient_outcomes,
    patient_outcomes
  )
)


# ============================================================
# 7. SAMPLE ALIAS VALIDATION
# ============================================================

required_alias_columns <- c(
  "object_sample_id",
  "sample_number",
  "sample_name",
  "sample_id"
)

check_required_columns(
  sample_aliases,
  required_alias_columns,
  "sample_aliases.csv"
)

check_no_missing_values(
  sample_aliases,
  required_alias_columns,
  "sample_aliases.csv"
)


# ------------------------------------------------------------
# Basic alias dimensions
# ------------------------------------------------------------

stopifnot(
  nrow(sample_aliases) == 32
)

stopifnot(
  length(
    unique(sample_aliases$object_sample_id)
  ) == 32
)

stopifnot(
  length(
    unique(sample_aliases$sample_number)
  ) == 32
)

stopifnot(
  length(
    unique(sample_aliases$sample_name)
  ) == 32
)

stopifnot(
  length(
    unique(sample_aliases$sample_id)
  ) == 32
)


# ------------------------------------------------------------
# Alias sample numbers must cover exactly 1:32
# ------------------------------------------------------------

alias_sample_numbers <- sort(
  as.integer(
    sample_aliases$sample_number
  )
)

stopifnot(
  identical(
    alias_sample_numbers,
    1:32
  )
)


# ------------------------------------------------------------
# Alias canonical sample names must match sample metadata
# ------------------------------------------------------------

stopifnot(
  identical(
    sort(sample_aliases$sample_name),
    sort(sample_meta$sample_name)
  )
)


# ------------------------------------------------------------
# Alias canonical sample IDs must match sample metadata
# ------------------------------------------------------------

stopifnot(
  identical(
    sort(sample_aliases$sample_id),
    sort(sample_meta$sample_id)
  )
)


# ------------------------------------------------------------
# Alias table must map exactly onto canonical sample metadata
# ------------------------------------------------------------

alias_canonical <- sample_aliases[
  ,
  c(
    "sample_number",
    "sample_name",
    "sample_id"
  )
]

sample_canonical <- sample_meta[
  ,
  c(
    "sample_number",
    "sample_name",
    "sample_id"
  )
]

alias_canonical$sample_number <- as.integer(
  alias_canonical$sample_number
)

sample_canonical$sample_number <- as.integer(
  sample_canonical$sample_number
)

alias_canonical <- alias_canonical[
  order(
    alias_canonical$sample_number
  ),
]

sample_canonical <- sample_canonical[
  order(
    sample_canonical$sample_number
  ),
]

rownames(alias_canonical) <- NULL
rownames(sample_canonical) <- NULL

stopifnot(
  identical(
    alias_canonical,
    sample_canonical
  )
)


# ------------------------------------------------------------
# Confirm expected legacy aliases for samples 9-16
# ------------------------------------------------------------

expected_legacy_w6_aliases <- data.frame(
  sample_number = 9:16,
  object_sample_id = c(
    "Sample_A_P76W6",
    "Sample_B_P80W6",
    "Sample_C_P54W6",
    "Sample_D_P43W6",
    "Sample_E_P87W6",
    "Sample_F_P29W6",
    "Sample_G_P46W6",
    "Sample_H_P42W6"
  ),
  stringsAsFactors = FALSE
)

observed_legacy_w6_aliases <- sample_aliases[
  sample_aliases$sample_number %in% 9:16,
  c(
    "sample_number",
    "object_sample_id"
  )
]

observed_legacy_w6_aliases$sample_number <-
  as.integer(
    observed_legacy_w6_aliases$sample_number
  )

observed_legacy_w6_aliases <-
  observed_legacy_w6_aliases[
    order(
      observed_legacy_w6_aliases$sample_number
    ),
  ]

rownames(
  observed_legacy_w6_aliases
) <- NULL

rownames(
  expected_legacy_w6_aliases
) <- NULL

stopifnot(
  identical(
    observed_legacy_w6_aliases,
    expected_legacy_w6_aliases
  )
)


# ------------------------------------------------------------
# Confirm remaining object IDs follow multiome_sample_X naming
# ------------------------------------------------------------

nonlegacy_aliases <- sample_aliases[
  !sample_aliases$sample_number %in% 9:16,
]

expected_object_ids <- paste0(
  "multiome_sample_",
  nonlegacy_aliases$sample_number
)

stopifnot(
  all(
    nonlegacy_aliases$object_sample_id ==
      expected_object_ids
  )
)


# ============================================================
# 8. THREE-FILE CROSS-CHECK
# ============================================================

# Join aliases to canonical sample metadata using sample_number.
# This ensures that every legacy object identifier resolves to
# exactly one canonical sample and patient.

alias_full <- merge(
  sample_aliases,
  sample_meta[
    ,
    c(
      "sample_number",
      "patient_id",
      "timepoint",
      "analysis_outcome"
    )
  ],
  by = "sample_number",
  all.x = TRUE,
  all.y = FALSE
)

stopifnot(
  nrow(alias_full) == 32
)

stopifnot(
  !any(
    is.na(alias_full$patient_id)
  )
)

stopifnot(
  !any(
    is.na(alias_full$timepoint)
  )
)

stopifnot(
  !any(
    is.na(alias_full$analysis_outcome)
  )
)


# ============================================================
# 9. REPORT
# ============================================================

cat("\n")
cat("============================================\n")
cat("Metadata validation successful\n")
cat("============================================\n\n")


cat(
  "Samples:",
  nrow(sample_meta),
  "\n"
)

cat(
  "Patients:",
  nrow(patient_meta),
  "\n"
)

cat(
  "Legacy sample aliases:",
  nrow(sample_aliases),
  "\n\n"
)


# ------------------------------------------------------------
# Sample summaries
# ------------------------------------------------------------

cat("Samples by timepoint:\n")
print(
  table(
    sample_meta$timepoint
  )
)

cat("\nPatients by outcome:\n")
print(
  patient_outcome_counts
)

cat("\nPatient sex:\n")
print(
  table(
    patient_meta$sex
  )
)

cat("\nTumor stage:\n")
print(
  table(
    patient_meta$tumor_stage
  )
)

cat("\nAUA risk:\n")
print(
  table(
    patient_meta$aua_risk
  )
)


# ------------------------------------------------------------
# Paired patient mapping
# ------------------------------------------------------------

cat("\nPaired sample mapping:\n")

paired <- reshape(
  sample_meta[
    ,
    c(
      "patient_id",
      "sample_id",
      "timepoint",
      "analysis_outcome"
    )
  ],
  idvar = c(
    "patient_id",
    "analysis_outcome"
  ),
  timevar = "timepoint",
  direction = "wide"
)

paired <- paired[
  order(
    paired$analysis_outcome,
    paired$patient_id
  ),
]

rownames(paired) <- NULL

print(
  paired,
  row.names = FALSE
)


# ------------------------------------------------------------
# Legacy alias mapping
# ------------------------------------------------------------

cat("\nLegacy object alias mapping:\n")

alias_report <- merge(
  sample_aliases,
  sample_meta[
    ,
    c(
      "sample_number",
      "patient_id",
      "timepoint"
    )
  ],
  by = "sample_number",
  all.x = TRUE
)

alias_report <- alias_report[
  order(
    alias_report$sample_number
  ),
]

alias_report <- alias_report[
  ,
  c(
    "sample_number",
    "object_sample_id",
    "sample_name",
    "sample_id",
    "patient_id",
    "timepoint"
  )
]

rownames(alias_report) <- NULL

print(
  alias_report,
  row.names = FALSE
)


# ------------------------------------------------------------
# Final confirmation
# ------------------------------------------------------------

cat("\n")
cat("All metadata checks passed.\n")
cat("\n")