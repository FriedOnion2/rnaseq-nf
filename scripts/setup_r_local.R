# One-time setup: install BiocManager + DESeq2 stack into a USER-level library
# (no admin needed, avoids writing to C:/Program Files/R/.../library).
# This lets you run the differential-expression step directly on this Windows machine
# while Docker/WSL2 is being set up.

# User library location (stable, outside Program Files)
user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "library")
dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))
cat("User library:", user_lib, "\n")

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = user_lib, repos = "https://cloud.r-project.org")
}
BiocManager::install(
  c("DESeq2", "SummarizedExperiment", "airway"),
  lib = user_lib, ask = FALSE, update = FALSE, Ncpus = 4
)
install.packages(c("ggplot2", "pheatmap", "RColorBrewer", "optparse"),
                 lib = user_lib, repos = "https://cloud.r-project.org")
cat("ALL DONE\n")