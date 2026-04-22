get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  match <- grep(file_arg, cmd_args)
  if (length(match) > 0) {
    return(normalizePath(sub(file_arg, "", cmd_args[match[1]])))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
  return(NA_character_)
}

script_path <- get_script_path()
project_root <- if (!is.na(script_path)) {
  normalizePath(file.path(dirname(script_path), ".."))
} else {
  normalizePath(getwd())
}

# Load model function
cat(" Loading best_model() from MS_MixGGM_all.R...\n")
source("/hopper/home/syooseph/projects/MPLN/CombinedModelSelection/MS_MixGGM_all.R")

# Read the CLR matrix
clr_fp <- file.path(project_root, "data", "processed", "lineage_clr_matrix_soil_filtered.csv")
cat(" Reading matrix from", clr_fp, "...\n")
X <- as.matrix(read.csv(clr_fp, row.names = 1, check.names = FALSE))

# Set parameters
Kval <- 20
d <- nrow(X)
method <- 1
seed <- 1001
niter <- 150
start <- 25

# Create output directory
outdir <- file.path(project_root, "models")
if (!dir.exists(outdir)) {
  dir.create(outdir)
  cat(" Created directory:", outdir, "\n")
}

# Run the model
cat("  Running best_model()...\n")
resultX <- best_model(X = X, M = Kval, beta = 1e-6, v = d, V = d * diag(d),
                      niter = niter, seed = seed, start = start, method = method)

# Save output
outfile <- file.path(outdir, "best_model_soil_k20.RData")
cat(" Saving result to", outfile, "...\n")
save(resultX, file = outfile, ascii = FALSE)

cat(" Done. Model complete.\n")
