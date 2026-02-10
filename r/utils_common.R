parse_bool <- function(x, default = TRUE) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) return(default)
  val <- tolower(as.character(x))
  if (val %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (val %in% c("false", "f", "0", "no", "n")) return(FALSE)
  return(default)
}

get_env_num <- function(name, default, min = NULL, max = NULL, invalid_value = default) {
  raw <- Sys.getenv(name, unset = NA_character_)
  if (is.na(raw)) return(default)
  if (!nzchar(raw)) return(invalid_value)
  val <- suppressWarnings(as.numeric(raw))
  if (!is.finite(val)) return(invalid_value)
  if (!is.null(min) && val < min) return(invalid_value)
  if (!is.null(max) && val > max) return(invalid_value)
  val
}

get_env_int <- function(name, default, min = NULL, max = NULL, invalid_value = default) {
  raw <- Sys.getenv(name, unset = NA_character_)
  if (is.na(raw)) return(default)
  if (!nzchar(raw)) return(invalid_value)
  val <- suppressWarnings(as.integer(raw))
  if (!is.finite(val)) return(invalid_value)
  if (!is.null(min) && val < min) return(invalid_value)
  if (!is.null(max) && val > max) return(invalid_value)
  val
}

get_env_num_fallback <- function(name, fallback_name, fallback_default) {
  raw <- Sys.getenv(name, unset = NA_character_)
  if (!is.na(raw)) {
    return(suppressWarnings(as.numeric(raw)))
  }
  raw_fb <- Sys.getenv(fallback_name, unset = NA_character_)
  if (is.na(raw_fb)) return(fallback_default)
  suppressWarnings(as.numeric(raw_fb))
}
