# used '18m_female' spleen for it 300 cells, 9 049 features
# used the seurat5 SCT fuction version (seurat5_integrated)
# only bootstrapping was measured
# the old version is not stephans, a rewritten one by me but also not parallelized at the cluster level

library(tools)
library(patchwork)
library(ggplot2)


base_dir <- "data/comparison/processing_time_tracking/logs"


parse_cores_label <- function(fname) {
  m <- stringr::str_match(fname, "core\\s*([0-9]+)")[,2]
  n <- suppressWarnings(as.integer(m))
  if (!is.na(n)) return(paste0(n, "_cores"))

  sub("\\.[^.]+$", "", fname)
}


parse_repeats <- function(fname, lines = character(0)) {
  r <- stringr::str_match(fname, "reps[_-]?([0-9]+)")[,2]
  r <- suppressWarnings(as.integer(r))
  if (!is.na(r)) return(r)
  # fallback from content if present
  if (length(lines)) {
    hit <- lines[str_detect(tolower(lines), "(repeats|n_cycles|ncycles|bootstrap)")]
    if (length(hit)) {
      nums <- str_extract_all(hit[1], "\\d+")[[1]]
      nums <- suppressWarnings(as.integer(nums))
      nums <- nums[!is.na(nums)]
      if (length(nums)) return(tail(nums, 1))
    }
  }
  NA_integer_
}

# read one file -> one row
parse_file <- function(path) {
  fname <- basename(path)
  lines <- tryCatch(readr::read_lines(path), error = function(e) character(0))

  # duration (seconds) -> hours
  dur_line <- lines[str_detect(lines, "Duration")]
  dur_val  <- suppressWarnings(as.numeric(str_extract(dur_line, "[0-9.]+")))
  dur_h    <- ifelse(length(dur_val) && !is.na(dur_val[1]), dur_val[1] / 3600, NA_real_)

  # RAM diff (assumed GB already)
  ram_line <- lines[str_detect(lines, "RAM diff")]
  ram_val  <- suppressWarnings(as.numeric(str_extract(ram_line, "[0-9.]+")))
  ram_gb   <- ifelse(length(ram_val) && !is.na(ram_val[1]), ram_val[1], NA_real_)

  tibble(
    cores_label      = parse_cores_label(fname),
    run              = sub("\\.[^.]+$", "", fname),
    repeats          = parse_repeats(fname, lines),
    `Duration in h`  = dur_h,
    `RAM used in GB` = ram_gb)
}


files <- list.files(base_dir, full.names = TRUE, recursive = FALSE)
files <- files[file.info(files)$isdir %in% c(FALSE, NA)]
if (length(files) == 0L) stop("No files found in: ", base_dir)

log_summary_df <- purrr::map_dfr(files, parse_file) %>%
  mutate(repeats = as.integer(repeats)) %>%
  filter(!is.na(repeats)) %>%
  arrange(cores_label, run, repeats)

# manually adding old rows
row1 <- tibble(
  cores_label      = "1_cores_previous",
  run              = "pipeline_log_core1_reps1_previous",
  repeats          = 1L,
  `Duration in h`  = 0.002766,
  `RAM used in GB` = 0)

row2 <- tibble(
  cores_label      = "1_cores_previous",
  run              = "pipeline_log_core1_reps100_previous",
  repeats          = 100L,
  `Duration in h`  = 0.007205,
  `RAM used in GB` = 0.08)

row3 <- tibble(
  cores_label      = "1_cores_previous",
  run              = "pipeline_log_core1_reps500_previous",
  repeats          = 500L,
  `Duration in h`  = 0.0076237,
  `RAM used in GB` = 0)

log_summary_df <- bind_rows(log_summary_df, row1)
log_summary_df <- bind_rows(log_summary_df, row2)
log_summary_df <- bind_rows(log_summary_df, row3)
log_summary_df$cores_label <- factor(log_summary_df$cores_label)

# plot
core_colors <- c(
  "1_cores"          = "#3B518B",
  "1_cores_previous" = "#30678D",
  "20_cores"        = "#268E89",
  "35_cores"        = "#42A574",
  "40_cores"        = "#9BD93C"
)

p_time <- ggplot( log_summary_df, aes(x = repeats, y = `Duration in h`, group =  cores_label, colour = cores_label)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  labs(x = "Repeats", y = "Duration (h)", colour = "Cores") +
   scale_x_log10(labels = scales::label_log()) +
     scale_color_manual(values = core_colors) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(size = 20),
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.position = "right")


ggsave("data/comparison/processing_time_tracking/tracking_time.png",
       p_time, width = 18, height = 12, units = "cm", dpi = 300)



 #RAM

p_ram <- ggplot(
  log_summary_df,
  aes(x = repeats, y = `RAM used in GB`,
      group = cores_label,
      colour = cores_label)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_text(aes(label = round(`RAM used in GB`, 2)),
            vjust = -0.7, size = 3, show.legend = FALSE) +
  labs(x = "Repeats", y = "RAM used (GB)", colour = "Cores") +
  scale_x_continuous(breaks = x_breaks) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "right")


ggsave("data/comparison/processing_time_tracking/tracking_ram.png",
       p_ram,  width = 18, height = 12, units = "cm", dpi = 300)







#  version for getting the total run directly from the processed folder


file_vec <- c("reproduced_vs3", "added_seurat5_vs3", "added_seurat5_integrated", "added_bpcells_fast", 'reproduced_fast')

log_summary_list <- list()
for (file in file_vec){
  log_file <- paste0("data/processed/", file, "/logs/pipeline_log.txt")
  lines <- read_lines(log_file)

  duration_line <- lines[str_detect(lines, "Duration")]
  duration_val  <- str_extract(duration_line, "[0-9\\.]+") %>% as.numeric()
  duration_h    <- duration_val / 3600

  ram_line <- lines[str_detect(lines, "RAM diff")]
  ram_val  <- str_extract(ram_line, "[0-9\\.]+") %>% as.numeric()

  log_summary_list[[file]] <- tibble(
    file              = file,
    `Duration in h`   = duration_h,
    `RAM used in GB`  = ram_val)}

log_summary_df <- bind_rows(log_summary_list)

log_summary_df <- log_summary_df %>%
  mutate(file = factor(file, levels = file_vec)) %>%
   mutate(version_label = dplyr::recode(file,
      "reproduced_vs3"         = "Reproduced",
      "added_seurat5_vs3"      = "Updated to Seurat 5",
      "added_seurat5_integrated" = "Seurat 5 & \n new bootstrapping",
      "added_bpcells_fast"     = "BPCells & \n new bootstrapping"))

plot <- ggplot(log_summary_df, aes(x = version_label, y = `Duration in h`, group = 1)) +
  geom_line(linewidth = 1, colour = "#4B9C79") +
  geom_point(size = 3, colour = "#4B9C79") +
  geom_text(aes(label = round(`Duration in h`, 2)), vjust = -0.7, size = 7, color = "#4B9C79") +
  labs(
    x = "Version",
    y = "Duration (h)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    panel.grid.minor = element_blank()  )

ggsave("data/comparison/processing_time_tracking/time_tracking.png", plot, width = 15, height = 10, dpi = 300)





#  enter the results manually


log_summary_df <- data.frame(
  version_label  = c("Original version", "Optimized bootstrapping", "Updated to Seurat5", "Memory optimized" ),
  `Duration in h` = c(26.55, 16.88, 3.15, 3.35 ),
  `Memory in GB`  = c(0, 15.39, 9.94, 0),
  check.names = FALSE)

# order
log_summary_df$version_label <- factor(log_summary_df$version_label,
                                        levels = log_summary_df$version_label)


# time
p_duration <- ggplot(log_summary_df, aes(x = version_label, y = `Duration in h`, group = 1)) +
  geom_line(linewidth = 1, colour = "#5FB4E5") +
  geom_point(size = 3, colour = "#5FB4E5") +
  geom_text(aes(label = round(`Duration in h`, 2)), vjust = -0.7, size = 7, color = "#5FB4E5") +
  scale_y_continuous(n.breaks = 10, limits = c(0, NA)) +
  labs(y = "Duration (h)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
      panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.3),
    panel.grid.minor   = element_blank())

# memory
p_memory <- ggplot(log_summary_df, aes(x = version_label, y = `Memory in GB`, group = 1)) +
  geom_line(linewidth = 1, colour = "#3CBFAE") +
  geom_point(size = 3, colour = "#3CBFAE") +
  geom_text(aes(label = round(`Memory in GB`, 2)), vjust = -0.7, size = 7, color = "#3CBFAE") +
 scale_y_continuous(n.breaks = 10, limits = c(0, NA)) +
  labs(y = "Memory (GB)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.3),
    panel.grid.minor   = element_blank())


plot <- p_duration + p_memory

ggsave("data/comparison/processing_time_tracking/time_tracking_combined.png", plot, width = 22, height = 12, dpi = 300)



# cells upscaling

# reproduced 100 000 there not all objects ran: load and count cells


#  enter the results manually
log_summary_df <- data.frame(
  version_label  = rep(c("10 000 cells", "50 000 cells", "100 000 cells"), 2),
  version        = rep(c("Original", "Improved"), each = 3),
  `Duration in h` = c(0, 0, 0,  # Original
                      0, 0, 8.0), # Improved 
  `Memory in GB`  = c(0, 0, 0,     # Original
                      0, 0, 16.7),    # Improved 
  check.names = FALSE)

# Order factors
log_summary_df$version_label <- factor(
  log_summary_df$version_label,
  levels = c("1 000 cells", "50 000 cells", "100 000 cells")
)


# time
p_duration <- ggplot(log_summary_df, aes(x = version_label, y = `Duration in h`, color = version, group = version)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_text(aes(label = round(`Duration in h`, 2)), vjust = -0.7, size = 5) +
  scale_y_continuous(n.breaks = 10, limits = c(0, NA)) +
  labs(y = "Duration (h)", color = "Version") +
  scale_color_manual(values = c("Original" = "#5FB4E5", "Improved" = "#FF7F50")) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    legend.position = "top"
  )

p_memory <- ggplot(log_summary_df, aes(x = version_label, y = `Memory in GB`, color = version, group = version)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_text(aes(label = round(`Memory in GB`, 2)), vjust = -0.7, size = 5) +
  scale_y_continuous(n.breaks = 10, limits = c(0, NA)) +
  labs(y = "Memory (GB)", color = "Version") +
  scale_color_manual(values = c("Original" = "#3CBFAE", "Improved" = "#FF7F50")) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    legend.position = "top"
  )



plot <- p_duration + p_memory

ggsave("data/comparison/processing_time_tracking/time_tracking_upscaling.png", plot, width = 22, height = 12, dpi = 300)



