# ---- Helpers ----
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(ggplot2)
  library(forcats); library(stringr); library(scales); library(tools)
})

# case-insensitive exact match
.pick_col <- function(nms, candidates) {
  nlow <- tolower(nms); cand_low <- tolower(candidates)
  hit <- match(cand_low, nlow)
  if (any(!is.na(hit))) nms[hit[which(!is.na(hit))[1]]] else NA_character_
}

# Read BP sheet from an Excel GO file (keeps original column names)
read_go_bp <- function(file, sheet_name = "BP") {
  # tolerate ".xlsxa" typo
  if (!file.exists(file) && grepl("\\.xlsxa$", file, ignore.case = TRUE)) {
    cand <- sub("\\.xlsxa$", ".xlsx", file, ignore.case = TRUE)
    if (file.exists(cand)) file <- cand
  }
  stopifnot(file.exists(file))
  
  sheets <- readxl::excel_sheets(file)
  sheets_trim <- trimws(sheets)
  # resolve BP sheet: exact name first, else any containing "bp", else first
  sheet_resolved <- if (sheet_name %in% sheets_trim) {
    sheets[match(sheet_name, sheets_trim)]
  } else {
    cand <- sheets[grepl("(^|\\b)bp(\\b|$)", tolower(sheets_trim))]
    if (length(cand)) cand[[1]] else sheets[[1]]
  }
  
  df <- readxl::read_excel(file, sheet = sheet_resolved)
  nms <- names(df)
  
  # Parse GeneRatio to numeric [0..1]
  gr_col <- .pick_col(nms, c("GeneRatio","Gene_Ratio","Gene Ratio"))
  if (!is.na(gr_col)) {
    if (is.character(df[[gr_col]])) {
      parts <- str_split(df[[gr_col]], "/", n = 2, simplify = TRUE)
      df$GeneRatio_num <- if (ncol(parts) == 2) suppressWarnings(as.numeric(parts[,1]) / as.numeric(parts[,2]))
      else suppressWarnings(as.numeric(df[[gr_col]]))
    } else {
      df$GeneRatio_num <- suppressWarnings(as.numeric(df[[gr_col]]))
    }
  } else df$GeneRatio_num <- NA_real_
  
  # Count column or estimate from BgRatio
  count_col <- .pick_col(nms, "Count")
  if (is.na(count_col)) {
    bgr <- .pick_col(nms, c("BgRatio","Bg_Ratio","Bg Ratio"))
    if (!is.na(bgr) && is.character(df[[bgr]]) && any(!is.na(df$GeneRatio_num))) {
      bgp <- str_split(df[[bgr]], "/", n = 2, simplify = TRUE)
      denom <- suppressWarnings(as.numeric(bgp[,2]))
      df$Count <- round(df$GeneRatio_num * denom)
      count_col <- "Count"
    }
  }
  
  df$.file  <- file
  df$.sheet <- sheet_resolved
  df
}

# Plot from a prepared BP tibble; uses zScore on y-axis
plot_go_bp_df <- function(
    bp_df,
    label            = NULL,   # used in file name and subtitle
    out_dir          = "BP_plots_top5_blue",
    top_n            = 5,
    min_count        = 0,
    max_count        = Inf,
    truncate_width   = 70,
    annotate         = TRUE,
    bar_color        = "#FFC107",
    width_in         = 12,
    height_in        = 6,
    dpi_val          = 300,
    title = "GO: Biological Process"
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  nms <- names(bp_df)
  desc_col <- .pick_col(nms, c("Description"))
  z_col    <- .pick_col(nms, c("zScore"))
  
  if (is.na(desc_col)) stop("Description/Term column not found.")
  if (is.na(z_col))    stop("zScore column not found.")
  
  count_col <- .pick_col(nms, "Count")
  
  df <- bp_df %>%
    mutate(
      Term_raw = as.character(.data[[desc_col]]),
      zscore   = suppressWarnings(as.numeric(.data[[z_col]])),
      CountVal = if (!is.na(count_col)) suppressWarnings(as.numeric(.data[[count_col]])) else NA_real_,
      GeneRatio_num = if ("GeneRatio_num" %in% names(.)) GeneRatio_num else NA_real_
    ) %>%
    filter(!is.na(Term_raw), !is.na(zscore)) %>%
    filter(is.na(CountVal) | (CountVal >= min_count & CountVal <= max_count)) %>%
    arrange(desc(abs(zscore))) %>%
    slice_head(n = top_n) %>%
    mutate(Term = str_trunc(Term_raw, width = truncate_width))
  
  if (!nrow(df)) stop("No rows left after filtering.")
  
  # Build plot: y = zScore (can be +/-), order by |z|
  p <- ggplot(df, aes(x = fct_reorder(Term, abs(zscore)), y = zscore)) +
    geom_col(width = 0.7, fill = bar_color) +
    coord_flip(clip = "off") +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey60") +
    labs(
      title = title,
      subtitle = if (!is.null(label)) paste0(label, " — ", unique(df$.sheet), " (Top ", nrow(df), " by |z|)") else NULL,
      x = NULL, y = "z-score"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title.position = "plot",
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(hjust = 1, lineheight = 0.95),
      plot.margin = margin(10, 30, 10, 10)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.10)))
  
  if (isTRUE(annotate) && (!all(is.na(df$CountVal)) || !all(is.na(df$GeneRatio_num)))) {
    lab <- if (!all(is.na(df$CountVal)) && !all(is.na(df$GeneRatio_num))) {
      paste0("n=", df$CountVal, "  r=", percent(df$GeneRatio_num, accuracy = 0.1))
    } else if (!all(is.na(df$CountVal))) paste0("n=", df$CountVal)
    else paste0("r=", percent(df$GeneRatio_num, accuracy = 0.1))
    p <- p + geom_text(aes(label = lab), hjust = -0.1, size = 3.4)
  }
  
  # File name
  base_name <- if (!is.null(label)) file_path_sans_ext(basename(label)) else "GO_BP"
  png_file  <- file.path(out_dir, paste0(base_name, "_", unique(df$.sheet), "_terms_top", nrow(df), "_zscore.png"))
  ggsave(png_file, p, width = width_in, height = height_in, dpi = dpi_val, bg = "white")
  
  invisible(list(plot = p, file = png_file, data = df))
}

# small inline helper (does not touch your existing functions)
add_gene_ratio <- function(df) {
  if (is.null(df) || !NROW(df)) return(df)
  if (!"GeneRatio_num" %in% names(df) && "GeneRatio" %in% names(df)) {
    parts <- strsplit(as.character(df$GeneRatio), "/", fixed = TRUE)
    num <- suppressWarnings(as.numeric(vapply(parts, `[`, "", 1)))
    den <- suppressWarnings(as.numeric(vapply(parts, `[`, "", 2)))
    df$GeneRatio_num <- ifelse(is.finite(den) & den > 0, num / den, NA_real_)
  }
  df
}

# Minimal wrappers for MF and CC
read_go_mf <- function(path) read_go_sheet(path, "MF")
read_go_cc <- function(path) read_go_sheet(path, "CC")


suppressMessages(library(readxl))

# --- Event-type -> color (match your pie) ---
evt_cols <- c(
  EX   = "#F8766D",  # coral
  INT  = "#7CAE00",  # green
  ALTA = "#00BFC4",  # cyan
  ALTD = "#C77CFF"   # purple
)

# pull EX/INT/ALTA/ALTD from any filename like:
# "GO_EX__FUS.xlsx", "genes_by_combo_INT__TAF15_GO.xlsx", etc.
get_event_from_name <- function(path) {
  b <- basename(path)
  m <- regmatches(b, regexpr("(EX|INT|ALTA|ALTD)", b, perl = TRUE))
  if (length(m) && nzchar(m)) m else NA_character_
}


input_folder <- "output/GO2"  # <- make sure this matches where your file is
folder <- list.files(input_folder, pattern = "\\.xlsx$", full.names = TRUE)

out_dir_bp <- file.path(input_folder, "BP_all_plots_top5_z")
out_dir_mf <- file.path(input_folder, "MF_all_plots_top5_z")
out_dir_cc <- file.path(input_folder, "CC_all_plots_top5_z")
dir.create(out_dir_bp, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_mf, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_cc, recursive = TRUE, showWarnings = FALSE)

top_n     <- 5
bar_color <- "#43A047"

for (file in folder) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " — ", basename(file), "\n", sep = "")
  f <- file
  evt <- get_event_from_name(f)
  bar_color <- if (!is.na(evt) && evt %in% names(evt_cols)) evt_cols[[evt]] else "#43A047"
  
  # # --- BP ---
  # bp <- tryCatch(read_go_bp(f), error = function(e) NULL)
  # if (is.null(bp)) bp <- tryCatch(read_xlsx(f, sheet = "BP"), error = function(e) NULL)
  # if (!is.null(bp) && "Description" %in% names(bp)) {
  #   plot_go_bp_df(bp, label = f, out_dir = out_dir_bp, top_n = top_n, bar_color = bar_color, title = "GO: Biological Process")
  # } else {
  #   cat("  • no BP enrichment results (Description column missing or read failed)\n")
  # }
  # 
  # # --- MF ---
  # mf <- tryCatch(read_go_mf(f), error = function(e) NULL)
  # if (is.null(mf)) mf <- tryCatch(read_xlsx(f, sheet = "MF"), error = function(e) NULL)
  # if (!is.null(mf) && "Description" %in% names(mf)) {
  #   plot_go_bp_df(mf, label = f, out_dir = out_dir_mf, top_n = top_n, bar_color = bar_color, title = "GO: Molecular function")
  # } else {
  #   cat("  • no MF enrichment results (Description column missing or read failed)\n")
  # }
  # 
  # # --- CC ---
  # cc <- tryCatch(read_go_cc(f), error = function(e) NULL)
  # if (is.null(cc)) cc <- tryCatch(read_xlsx(f, sheet = "CC"), error = function(e) NULL)
  # if (!is.null(cc) && "Description" %in% names(cc)) {
  #   plot_go_bp_df(cc, label = f, out_dir = out_dir_cc, top_n = top_n, bar_color = bar_color, title = "GO: Cellular component")
  # } else {
  #   cat("  • no CC enrichment results (Description column missing or read failed)\n")
  # }
  # --- BP ---
  bp <- tryCatch(read_go_bp(f), error = function(e) NULL)
  if (is.null(bp)) bp <- tryCatch(read_xlsx(f, sheet = "BP"), error = function(e) NULL)
  bp <- add_gene_ratio(bp)
  if (!is.null(bp) && "Description" %in% names(bp)) {
    plot_go_bp_df(bp, label = f, out_dir = out_dir_bp, top_n = top_n,
                  bar_color = bar_color, title = "GO: Biological Process")
  } else cat("  • no BP enrichment results (Description missing or read failed)\n")
  
  # --- MF ---
  mf <- tryCatch(read_go_mf(f), error = function(e) NULL)
  if (is.null(mf)) mf <- tryCatch(read_xlsx(f, sheet = "MF"), error = function(e) NULL)
  mf <- add_gene_ratio(mf)
  if (!is.null(mf) && "Description" %in% names(mf)) {
    plot_go_bp_df(mf, label = f, out_dir = out_dir_mf, top_n = top_n,
                  bar_color = bar_color, title = "GO: Molecular function")
  } else cat("  • no MF enrichment results (Description missing or read failed)\n")
  
  # --- CC ---
  cc <- tryCatch(read_go_cc(f), error = function(e) NULL)
  if (is.null(cc)) cc <- tryCatch(read_xlsx(f, sheet = "CC"), error = function(e) NULL)
  cc <- add_gene_ratio(cc)
  if (!is.null(cc) && "Description" %in% names(cc)) {
    plot_go_bp_df(cc, label = f, out_dir = out_dir_cc, top_n = top_n,
                  bar_color = bar_color, title = "GO: Cellular component")
  } else cat("  • no CC enrichment results (Description missing or read failed)\n")
  
}

