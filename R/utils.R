`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    y
  } else {
    x
  }
}

check_scalar_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1 || is.na(x) || !nzchar(trimws(x))) {
    stop(name, " must be a non-empty character string", call. = FALSE)
  }
  trimws(x)
}

normalize_optional_text <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return("")
  }
  value_chr <- trimws(as.character(value)[1])
  if (toupper(value_chr) %in% c("", "-", "NA", "NULL", "NONE")) {
    return("")
  }
  value_chr
}

parse_optional_numeric <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return(NA_real_)
  }
  value_chr <- trimws(as.character(value)[1])
  if (toupper(value_chr) %in% c("", "-", "NA", "NULL", "NONE")) {
    return(NA_real_)
  }
  out <- suppressWarnings(as.numeric(value_chr))
  if (length(out) == 0) {
    NA_real_
  } else {
    out[1]
  }
}

parse_gene_list <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return(character())
  }
  if (length(value) == 1) {
    value <- unlist(strsplit(as.character(value), "[;,]", perl = TRUE), use.names = FALSE)
  }
  value <- trimws(as.character(value))
  unique(value[nzchar(value) & !is.na(value)])
}

validate_text_side <- function(value, name) {
  value <- tolower(check_scalar_string(value, name))
  if (!value %in% c("left", "right", "auto")) {
    stop(name, " must be 'left', 'right', or 'auto'", call. = FALSE)
  }
  value
}

as_plain_data_frame <- function(data) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame", call. = FALSE)
  }
  as.data.frame(data, stringsAsFactors = FALSE)
}
