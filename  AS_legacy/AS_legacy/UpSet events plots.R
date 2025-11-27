# 
library(UpSetR)
library(grid)
for (event_name in names(partitioned_events)) {
  out_dir  <- "output"                        
  event_type <- event_name  # choose from names(partitioned_events)
  out_name <- paste0("UpSetR_", event_type, ".png")
  out_path <- file.path(out_dir, out_name)
  # --- build expression counts from your partitioned list ---
  lst <- partitioned_events[[event_type]]                 # your exclusive combos
  expr_counts <- lengths(lst)
  
  # convert names like "FUS_EWSR1_TAF15" -> "EWSR1KO&FUSKO&TAF15KO"
  to_upset_names <- function(v) {
    nm <- names(v)
    nm2 <- vapply(strsplit(nm, "_"),
                  function(parts) paste0(paste0(parts, "KO"), collapse = "&"),
                  "")
    names(v) <- nm2
    v
  }
  expr_counts <- to_upset_names(expr_counts)
  
  bm_upset <- UpSetR::fromExpression(expr_counts)
  
  # --- appearance ---
  sets     <- c("EWSR1KO","FUSKO","TAF15KO")
  set_cols <- c(EWSR1KO = "#d53031", FUSKO = "#009E73", TAF15KO = "#0072B2")
  
  qrys <- list(
    list(query = intersects, params = list("EWSR1KO"), color = set_cols["EWSR1KO"], active = TRUE),
    list(query = intersects, params = list("FUSKO"),   color = set_cols["FUSKO"],   active = TRUE),
    list(query = intersects, params = list("TAF15KO"), color = set_cols["TAF15KO"], active = TRUE)
  )
  
  # --- draw ---
  png(out_path, width = 12, height = 8, units = "in", res = 300)
  grid::grid.newpage()
  UpSetR::upset(
    bm_upset,
    sets               = sets,
    keep.order         = TRUE,
    order.by           = "degree",
    decreasing         = c(TRUE, TRUE),
    sets.bar.color     = unname(set_cols[sets]),
    mainbar.y.label    = "Intersection size",
    sets.x.label       = NULL,
    matrix.color       = "black",
    main.bar.color     = "black",
    point.size         = 5,
    line.size          = 1,
    mb.ratio           = c(0.70, 0.30),
    text.scale         = 4,
    queries            = qrys
  )
  dev.off()
}

