# ============================================================
# 00_validate_metadata.R
#
# Validate sample- and patient-level metadata for the
# NMIBC BCG single-cell multiome reproducibility repository.
#
# This script uses base R only so that metadata validation
# does not depend on any external packages.
# ============================================================


# ------------------------------------------------------------
# File paths
# ------------------------------------------------------------

sample_metadata_path <- file.path(
  "metadata",
  "sample_metadata.csv"
)

patient_metadata_path <- file.path(
  "metadata",
  "patient_metadata.csv"
)


# ------------------------------------------------------------
# Check files exist
# ------------------------------------------------------------

if (!file.exists(sample_metadata_path)) {
  stop(
    "Sample metadata file not found: ",
    sample_metadata_path
  )
}

if (!file.exists(patient_metadata_path)) {
  stop(
    "Patient metadata file not found: ",
    patient_metadata_path
  )
}


# ------------------------------------------------------------
# Read metadata
# ------------------------------------------------------------

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


# ============================================================
# SAMPLE METADATA VALIDATION
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

missing_sample_columns <- setdiff(
  required_sample_columns,
  colnames(sample_meta)
)

if (length(missing_sample_columns) > 0) {
  stop(
    "Missing required sample metadata columns: ",
    paste(missing_sample_columns, collapse = ", ")
  )
}


# ------------------------------------------------------------
# Basic sample counts
# ------------------------------------------------------------

stopifnot(nrow(sample_meta) == 32)

stopifnot(
  length(unique(sample_meta$sample_number)) == 32
)

stopifnot(
  length(unique(sample_meta$sample_name)) == 32
)

stopifnot(
  length(unique(sample_meta$sample_id)) == 32
)

stopifnot(
  length(unique(sample_meta$patient_id)) == 16
)


# ------------------------------------------------------------
# Expected values
# ------------------------------------------------------------

stopifnot(
  all(sample_meta$timepoint %in% c("pre", "post"))
)

stopifnot(
  all(sample_meta$bcg_week %in% c("W1", "W6"))
)

stopifnot(
  all(
    sample_meta$sample_context %in%
      c("pre-BCG", "pre-sixth_BCG")
  )
)

stopifnot(
  all(
    sample_meta$treatment %in%
      c("untreated", "BCG-treated")
  )
)

stopifnot(
  all(
    sample_meta$outcome_raw %in%
      c("remission", "recurrence")
  )
)

stopifnot(
  all(
    sample_meta$analysis_outcome %in%
      c("recurrence_free", "early_recurrence")
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
# Each patient should have exactly two samples
# ------------------------------------------------------------

samples_per_patient <- table(
  sample_meta$patient_id
)

stopifnot(
  all(samples_per_patient == 2)
)


# ------------------------------------------------------------
# Every patient should have exactly one pre and one post sample
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
# Outcome should be identical across paired samples
# ------------------------------------------------------------

outcomes_per_patient <- tapply(
  sample_meta$analysis_outcome,
  sample_meta$patient_id,
  function(x) length(unique(x))
)

stopifnot(
  all(outcomes_per_patient == 1)
)


# ============================================================
# PATIENT METADATA VALIDATION
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

missing_patient_columns <- setdiff(
  required_patient_columns,
  colnames(patient_meta)
)

if (length(missing_patient_columns) > 0) {
  stop(
    "Missing required patient metadata columns: ",
    paste(missing_patient_columns, collapse = ", ")
  )
}


# ------------------------------------------------------------
# Patient-level counts
# ------------------------------------------------------------

stopifnot(nrow(patient_meta) == 16)

stopifnot(
  length(unique(patient_meta$patient_id)) == 16
)


# ------------------------------------------------------------
# Ensure the same patients exist in both metadata files
# ------------------------------------------------------------

sample_patients <- sort(
  unique(sample_meta$patient_id)
)

patient_patients <- sort(
  unique(patient_meta$patient_id)
)

stopifnot(
  identical(sample_patients, patient_patients)
)


# ------------------------------------------------------------
# Patient outcome distribution
# ------------------------------------------------------------

patient_outcome_counts <- table(
  patient_meta$analysis_outcome
)

stopifnot(
  patient_outcome_counts["recurrence_free"] == 8
)

stopifnot(
  patient_outcome_counts["early_recurrence"] == 8
)


# ------------------------------------------------------------
# Cross-check outcome between sample and patient tables
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
  order(sample_patient_outcomes$patient_id),
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
  order(patient_outcomes$patient_id),
]

rownames(sample_patient_outcomes) <- NULL
rownames(patient_outcomes) <- NULL

stopifnot(
  identical(
    sample_patient_outcomes,
    patient_outcomes
  )
)


# ------------------------------------------------------------
# Clinical consistency checks
# ------------------------------------------------------------

stopifnot(
  all(patient_meta$sex %in% c("Male", "Female"))
)

stopifnot(
  all(patient_meta$tumor_stage %in% c("Ta", "T1"))
)

stopifnot(
  all(patient_meta$aua_risk %in% c(2, 3))
)

stopifnot(
  sum(patient_meta$tumor_stage == "T1") == 3
)

stopifnot(
  sum(patient_meta$tumor_stage == "Ta") == 13
)

stopifnot(
  sum(patient_meta$sex == "Female") == 4
)

stopifnot(
  sum(patient_meta$sex == "Male") == 12
)


# ============================================================
# REPORT
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
  "\n\n"
)


cat("Samples by timepoint:\n")
print(
  table(sample_meta$timepoint)
)

cat("\nPatients by outcome:\n")
print(
  patient_outcome_counts
)


cat("\nPatient sex:\n")
print(
  table(patient_meta$sex)
)


cat("\nTumor stage:\n")
print(
  table(patient_meta$tumor_stage)
)


cat("\nAUA risk:\n")
print(
  table(patient_meta$aua_risk)
)


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

print(
  paired,
  row.names = FALSE
)

cat("\n")
cat("All metadata checks passed.\n")
cat("\n")