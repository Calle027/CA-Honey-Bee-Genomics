###############################
## 00_rstudio_helpers.R
## Helpers for locating scripts when run from RStudio or Rscript.
###############################

resolve_script_path <- function(fallback_name = NULL) {
  script_path <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
    if (length(file_arg) > 0) {
      return(normalizePath(file_arg[1], mustWork = TRUE))
    }

    frame_files <- vapply(
      sys.frames(),
      function(frame) if (!is.null(frame$ofile)) frame$ofile else NA_character_,
      character(1)
    )
    frame_files <- frame_files[!is.na(frame_files)]
    if (length(frame_files) > 0) {
      return(normalizePath(tail(frame_files, 1), mustWork = TRUE))
    }

    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      context_path <- rstudioapi::getActiveDocumentContext()$path
      if (!is.null(context_path) && nzchar(context_path)) {
        return(normalizePath(context_path, mustWork = TRUE))
      }
    }

    NA_character_
  }, error = function(e) NA_character_)

  if (!is.na(script_path)) {
    return(script_path)
  }

  if (!is.null(fallback_name)) {
    fallback_path <- file.path(getwd(), "honeybee_scripts", fallback_name)
    if (file.exists(fallback_path)) {
      return(normalizePath(fallback_path, mustWork = TRUE))
    }
  }

  NA_character_
}

resolve_project_root <- function(script_path = NA_character_) {
  if (!is.na(script_path)) {
    script_dir <- dirname(script_path)
    if (basename(script_dir) == "honeybee_scripts") {
      return(normalizePath(file.path(script_dir, ".."), mustWork = TRUE))
    }
    return(normalizePath(script_dir, mustWork = TRUE))
  }

  normalizePath(getwd(), mustWork = TRUE)
}
