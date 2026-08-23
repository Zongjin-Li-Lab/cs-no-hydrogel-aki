project_root <- function() {
  root <- Sys.getenv("CS_NO_PROJECT_ROOT", unset = getwd())
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)

  if (!file.exists(file.path(root, "README.md")) ||
      !dir.exists(file.path(root, "scripts"))) {
    stop(
      "Run this script from the repository root or set CS_NO_PROJECT_ROOT.\n",
      "Current candidate root: ", root
    )
  }

  root
}

make_output_dirs <- function(root) {
  dirs <- list(
    base = root,
    raw_data = file.path(root, "data", "raw"),
    processed = file.path(root, "data", "processed"),
    results = file.path(root, "results", "generated"),
    qc = file.path(root, "results", "generated", "01_QC"),
    clustering = file.path(root, "results", "generated", "02_Clustering"),
    de = file.path(root, "results", "generated", "03_DifferentialExpression"),
    trajectory = file.path(root, "results", "generated", "04_Trajectory"),
    cellchat = file.path(root, "results", "generated", "05_CellChat"),
    advanced = file.path(root, "results", "generated", "06_Advanced"),
    nichenet = file.path(root, "results", "generated", "07_NicheNet"),
    figures = file.path(root, "figures"),
    code = file.path(root, "scripts")
  )

  output_dirs <- dirs[setdiff(names(dirs), c("base", "raw_data", "code"))]
  invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
}

assert_file <- function(path, hint = NULL) {
  if (!file.exists(path)) {
    message_text <- paste0("Required file not found: ", path)
    if (!is.null(hint)) message_text <- paste(message_text, hint, sep = "\n")
    stop(message_text)
  }
  invisible(path)
}

set.seed(20260614)
