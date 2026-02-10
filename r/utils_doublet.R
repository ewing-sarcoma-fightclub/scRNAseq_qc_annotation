coerce_logical_vec <- function(x) {
  if (is.null(x)) return(logical())
  if (is.logical(x)) return(x)
  if (is.numeric(x)) {
    out <- x != 0
    out[is.na(x)] <- NA
    return(out)
  }
  val <- tolower(as.character(x))
  out <- rep(NA, length(val))
  out[val %in% c("true", "t", "1", "yes", "y", "doublet")] <- TRUE
  out[val %in% c("false", "f", "0", "no", "n", "singlet")] <- FALSE
  out
}

resolve_doublet_calls <- function(meta) {
  df_col <- NULL
  if ("doubletfinder_call" %in% colnames(meta)) {
    df_col <- "doubletfinder_call"
  } else {
    df_cols <- grep("^DF\\.classifications", colnames(meta), value = TRUE)
    if (length(df_cols) >= 1) df_col <- df_cols[[length(df_cols)]]
  }
  df_call <- if (!is.null(df_col)) coerce_logical_vec(meta[[df_col]]) else rep(NA, nrow(meta))

  scrub_col <- NULL
  for (c in c("scrublet_call", "predicted_doublets", "predicted_doublet", "scrublet_doublet")) {
    if (c %in% colnames(meta)) {
      scrub_col <- c
      break
    }
  }
  scrub_call <- if (!is.null(scrub_col)) coerce_logical_vec(meta[[scrub_col]]) else rep(NA, nrow(meta))

  list(df_call = df_call, scrub_call = scrub_call, df_col = df_col, scrub_col = scrub_col)
}

combine_doublet_calls <- function(df_call, scrub_call, mode = "doubletfinder") {
  mode <- tolower(mode)
  if (!(mode %in% c("doubletfinder", "scrublet", "union", "intersection"))) {
    mode <- "doubletfinder"
  }
  has_df <- !all(is.na(df_call))
  has_scrub <- !all(is.na(scrub_call))

  df_call2 <- if (has_df) df_call else rep(NA, length(df_call))
  scrub_call2 <- if (has_scrub) scrub_call else rep(NA, length(scrub_call))

  df_call2[is.na(df_call2)] <- FALSE
  scrub_call2[is.na(scrub_call2)] <- FALSE

  union_call <- df_call2 | scrub_call2
  intersection_call <- df_call2 & scrub_call2

  if (mode == "scrublet") {
    doublet_call <- if (has_scrub) scrub_call2 else df_call2
  } else if (mode == "union") {
    doublet_call <- if (has_df || has_scrub) union_call else rep(FALSE, length(df_call2))
  } else if (mode == "intersection") {
    doublet_call <- if (has_df && has_scrub) intersection_call else if (has_df) df_call2 else scrub_call2
  } else {
    doublet_call <- if (has_df) df_call2 else scrub_call2
  }

  list(
    doublet_call = doublet_call,
    union_call = union_call,
    intersection_call = intersection_call,
    has_df = has_df,
    has_scrub = has_scrub
  )
}
