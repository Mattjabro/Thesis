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

args <- commandArgs(trailingOnly = TRUE)
model_fp <- if (length(args) >= 1) args[1] else {
  file.path(project_root, "models", "best_model_soil_k20.RData")
}
sra_fp <- if (length(args) >= 2) args[2] else {
  file.path(project_root, "data", "raw", "sraRunTable.csv")
}
clr_fp <- if (length(args) >= 3) args[3] else {
  file.path(project_root, "data", "processed", "lineage_clr_matrix_soil_filtered.csv")
}
out_dir <- if (length(args) >= 4) args[4] else {
  file.path(project_root, "outputs", "cluster_metadata")
}
metadata_out_fp <- if (length(args) >= 5) args[5] else {
  file.path(project_root, "data", "intermediate", "soil_metadata.tsv")
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

loaded <- load(model_fp)
model_obj <- if ("resultX" %in% loaded) {
  resultX
} else {
  get(loaded[1])
}

X <- as.matrix(read.csv(clr_fp, row.names = 1, check.names = FALSE))
sample_ids <- colnames(X)

extract_clusters <- function(obj, sample_ids) {
  n_samples <- length(sample_ids)
  summaries <- list()
  best <- list(score = -Inf, assign = NULL, path = NA_character_, note = NA_character_)

  record_summary <- function(path, kind, nrow, ncol, len, n_unique, score, note) {
    summaries[[length(summaries) + 1]] <<- data.frame(
      path = path,
      kind = kind,
      nrow = nrow,
      ncol = ncol,
      length = len,
      n_unique = n_unique,
      score = score,
      note = note,
      stringsAsFactors = FALSE
    )
  }

  normalize_assign <- function(vec) {
    if (is.factor(vec)) {
      return(as.integer(vec))
    }
    if (is.character(vec)) {
      return(as.integer(factor(vec)))
    }
    if (is.numeric(vec)) {
      ok <- !is.na(vec)
      if (all(abs(vec[ok] - round(vec[ok])) < 1e-6)) {
        return(as.integer(round(vec)))
      }
    }
    return(as.integer(factor(vec)))
  }

  consider_vector <- function(vec, path) {
    if (!is.vector(vec) || is.list(vec) || length(vec) != n_samples) {
      return(NULL)
    }

    n_unique <- length(unique(vec))
    score <- 0
    note <- ""
    if (is.factor(vec) || is.integer(vec)) {
      score <- score + 3
    }
    if (is.numeric(vec)) {
      ok <- !is.na(vec)
      if (all(abs(vec[ok] - round(vec[ok])) < 1e-6)) {
        score <- score + 2
      }
    }
    if (!is.null(names(vec)) && all(names(vec) %in% sample_ids)) {
      score <- score + 3
      note <- "names_match_samples"
    }
    if (n_unique >= 2 && n_unique <= min(100, n_samples / 2)) {
      score <- score + 2
    }
    if (n_unique == n_samples || n_unique == 1) {
      score <- score - 5
    }

    record_summary(path, "vector", NA_integer_, NA_integer_, length(vec), n_unique, score, note)

    if (score > best$score) {
      best$score <<- score
      best$assign <<- normalize_assign(vec)
      best$path <<- path
      best$note <<- note
    }
  }

  consider_matrix <- function(mat, path) {
    if (!(is.matrix(mat) || is.data.frame(mat))) {
      return(NULL)
    }
    mat <- as.matrix(mat)
    nrow_mat <- nrow(mat)
    ncol_mat <- ncol(mat)

    if (nrow_mat == n_samples && ncol_mat == 1) {
      consider_vector(mat[, 1], paste0(path, "[,1]"))
      return(NULL)
    }
    if (ncol_mat == n_samples && nrow_mat == 1) {
      consider_vector(mat[1, ], paste0(path, "[1,]"))
      return(NULL)
    }
    if (!(nrow_mat == n_samples || ncol_mat == n_samples)) {
      return(NULL)
    }
    if (!is.numeric(mat)) {
      return(NULL)
    }

    score <- 0
    note <- ""
    assign <- NULL
    tol <- 1e-2
    if (!is.null(rownames(mat)) && all(rownames(mat) %in% sample_ids)) {
      score <- score + 3
    }
    if (!is.null(colnames(mat)) && all(colnames(mat) %in% sample_ids)) {
      score <- score + 3
    }
    if (nrow_mat == n_samples) {
      row_ok <- mean(abs(rowSums(mat) - 1) < tol, na.rm = TRUE)
      if (is.finite(row_ok) && row_ok > 0.9) {
        score <- score + 5
        note <- "rowsum~1"
        assign <- max.col(mat)
      }
    }
    if (is.null(assign) && ncol_mat == n_samples) {
      col_ok <- mean(abs(colSums(mat) - 1) < tol, na.rm = TRUE)
      if (is.finite(col_ok) && col_ok > 0.9) {
        score <- score + 5
        note <- "colsum~1"
        assign <- max.col(t(mat))
      }
    }

    record_summary(path, "matrix", nrow_mat, ncol_mat, NA_integer_, NA_integer_, score, note)

    if (!is.null(assign) && score > best$score) {
      best$score <<- score
      best$assign <<- as.integer(assign)
      best$path <<- path
      best$note <<- note
    }
  }

  consider_group_list <- function(lst, path) {
    if (!is.list(lst) || length(lst) < 2) {
      return(NULL)
    }

    if (all(sapply(lst, function(x) is.atomic(x) && is.numeric(x)))) {
      vals <- unlist(lst, use.names = FALSE)
      if (length(vals) == 0) {
        return(NULL)
      }
      if (min(vals) >= 1 && max(vals) <= n_samples) {
        if (length(unique(vals)) == n_samples) {
          assign <- rep(NA_integer_, n_samples)
          for (k in seq_along(lst)) {
            assign[as.integer(lst[[k]])] <- k
          }
          if (!any(is.na(assign))) {
            record_summary(path, "index_list", NA_integer_, NA_integer_, NA_integer_, NA_integer_, 6, "index_list")
            if (6 > best$score) {
              best$score <<- 6
              best$assign <<- assign
              best$path <<- path
              best$note <<- "index_list"
            }
          }
        }
      }
    }

    if (all(sapply(lst, is.character))) {
      vals <- unlist(lst, use.names = FALSE)
      if (length(vals) == 0) {
        return(NULL)
      }
      if (all(vals %in% sample_ids) && length(unique(vals)) == n_samples) {
        assign <- rep(NA_integer_, n_samples)
        names(assign) <- sample_ids
        for (k in seq_along(lst)) {
          assign[match(lst[[k]], sample_ids)] <- k
        }
        if (!any(is.na(assign))) {
          record_summary(path, "id_list", NA_integer_, NA_integer_, NA_integer_, NA_integer_, 6, "id_list")
          if (6 > best$score) {
            best$score <<- 6
            best$assign <<- assign
            best$path <<- path
            best$note <<- "id_list"
          }
        }
      }
    }
  }

  if (is.list(obj)) {
    candidates <- c(
      "clust", "cluster", "clusters", "classification", "class", "z", "Z",
      "best_clust", "best_cluster", "best_class", "membership", "posterior", "assign"
    )
    for (nm in candidates) {
      val <- obj[[nm]]
      if (!is.null(val)) {
        if (is.matrix(val) || is.data.frame(val)) {
          mat <- as.matrix(val)
          if (nrow(mat) == n_samples) {
            return(list(assign = max.col(mat), candidates = NULL))
          }
          if (ncol(mat) == n_samples) {
            return(list(assign = max.col(t(mat)), candidates = NULL))
          }
        }
        if (is.vector(val) && length(val) == n_samples && !is.list(val)) {
          return(list(assign = normalize_assign(val), candidates = NULL))
        }
      }
    }
  }

  max_depth <- 4
  walk <- function(x, path, depth) {
    if (depth > max_depth) {
      return(NULL)
    }
    if (is.list(x)) {
      consider_group_list(x, path)
      nms <- names(x)
      for (i in seq_along(x)) {
        nm <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else paste0("[[", i, "]]")
        walk(x[[i]], paste0(path, "$", nm), depth + 1)
      }
    } else if (is.matrix(x) || is.data.frame(x)) {
      consider_matrix(x, path)
    } else {
      consider_vector(x, path)
    }
  }

  walk(obj, "model_obj", 0)

  candidates_df <- if (length(summaries) > 0) {
    do.call(rbind, summaries)
  } else {
    data.frame(
      path = character(0),
      kind = character(0),
      nrow = integer(0),
      ncol = integer(0),
      length = integer(0),
      n_unique = integer(0),
      score = numeric(0),
      note = character(0),
      stringsAsFactors = FALSE
    )
  }

  return(list(assign = best$assign, candidates = candidates_df))
}

cluster_res <- extract_clusters(model_obj, sample_ids)
cluster_assign <- cluster_res$assign
if (is.null(cluster_assign)) {
  candidates_fp <- file.path(out_dir, "model_object_candidates.tsv")
  write.table(
    cluster_res$candidates,
    candidates_fp,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  capture.output(
    str(model_obj, max.level = 2),
    file = file.path(out_dir, "model_object_structure.txt")
  )
  stop(
    "Could not find cluster assignments in model object. ",
    "Wrote candidates to ", candidates_fp, "."
  )
}
if (length(cluster_assign) != length(sample_ids)) {
  stop(
    "Cluster assignment length (", length(cluster_assign),
    ") does not match number of samples (", length(sample_ids), ")."
  )
}

cluster_df <- data.frame(
  X.SampleID = sample_ids,
  cluster = as.integer(cluster_assign),
  stringsAsFactors = FALSE
)

write.table(
  cluster_df,
  file.path(out_dir, "cluster_assignments.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

sra <- read.csv(sra_fp, stringsAsFactors = FALSE, check.names = FALSE)

core_ids <- sub("(_[FR]_filt\\.fastq\\.gz|_filt\\.fastq\\.gz)$", "", sample_ids)
map_df <- data.frame(
  X.SampleID = sample_ids,
  Run = core_ids,
  stringsAsFactors = FALSE
)

meta_keep <- c(
  "Run",
  "environment_(biome)",
  "environment_(feature)",
  "environment_(material)",
  "geo_loc_name_country",
  "geographic_location_(latitude)",
  "geographic_location_(longitude)",
  "geographic_location_(depth)",
  "geographic_location_(elevation)",
  "Collection_Date"
)
meta_keep <- meta_keep[meta_keep %in% names(sra)]
meta <- sra[, meta_keep, drop = FALSE]

names(meta) <- c(
  "Run",
  "env_biome",
  "env_feature",
  "env_material",
  "country",
  "latitude_deg",
  "longitude_deg",
  "depth_m",
  "elevation_m",
  "collection_date"
)[seq_along(names(meta))]

to_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}
for (nm in c("latitude_deg", "longitude_deg", "depth_m", "elevation_m")) {
  if (nm %in% names(meta)) {
    meta[[nm]] <- to_numeric(meta[[nm]])
  }
}

merged_meta <- merge(map_df, meta, by = "Run", all.x = TRUE)
merged_meta <- merged_meta[, c("X.SampleID", setdiff(names(merged_meta), c("X.SampleID", "Run"))), drop = FALSE]

write.table(
  merged_meta,
  metadata_out_fp,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

merged <- merge(cluster_df, merged_meta, by = "X.SampleID", all.x = TRUE)
write.table(
  merged,
  file.path(out_dir, "cluster_metadata_join.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cluster_sizes <- as.data.frame(table(cluster_df$cluster), stringsAsFactors = FALSE)
names(cluster_sizes) <- c("cluster", "n_samples")
write.table(
  cluster_sizes,
  file.path(out_dir, "cluster_sizes.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat_vars <- c("env_biome", "env_feature", "env_material", "country")
num_vars <- c("latitude_deg", "longitude_deg", "depth_m", "elevation_m")

tests <- data.frame(
  variable = character(0),
  type = character(0),
  method = character(0),
  p_value = numeric(0),
  stringsAsFactors = FALSE
)

for (var in cat_vars) {
  if (!var %in% names(merged)) {
    next
  }
  tab <- table(merged$cluster, merged[[var]])
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    next
  }
  pval <- tryCatch(chisq.test(tab)$p.value, error = function(e) NA_real_)
  tests <- rbind(
    tests,
    data.frame(
      variable = var,
      type = "categorical",
      method = "chisq",
      p_value = pval,
      stringsAsFactors = FALSE
    )
  )
}

for (var in num_vars) {
  if (!var %in% names(merged)) {
    next
  }
  vals <- merged[[var]]
  ok <- !is.na(vals) & !is.na(merged$cluster)
  if (length(unique(merged$cluster[ok])) < 2) {
    next
  }
  pval <- tryCatch(kruskal.test(vals[ok] ~ as.factor(merged$cluster[ok]))$p.value,
    error = function(e) NA_real_
  )
  tests <- rbind(
    tests,
    data.frame(
      variable = var,
      type = "numeric",
      method = "kruskal",
      p_value = pval,
      stringsAsFactors = FALSE
    )
  )
}

if (nrow(tests) > 0) {
  tests$padj <- p.adjust(tests$p_value, method = "BH")
}

write.table(
  tests,
  file.path(out_dir, "cluster_metadata_tests.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Wrote outputs to", out_dir, "\n")
