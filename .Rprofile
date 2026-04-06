quarto_pandoc_dir <- switch(
  Sys.info()[["machine"]],
  arm64 = "/Applications/quarto/bin/tools/aarch64",
  aarch64 = "/Applications/quarto/bin/tools/aarch64",
  x86_64 = "/Applications/quarto/bin/tools/x86_64",
  NULL
)

if (!nzchar(Sys.getenv("RSTUDIO_PANDOC")) &&
    !is.null(quarto_pandoc_dir) &&
    file.exists(file.path(quarto_pandoc_dir, "pandoc"))) {
  Sys.setenv(RSTUDIO_PANDOC = quarto_pandoc_dir)
}

source("renv/activate.R")
