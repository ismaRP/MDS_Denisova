
library(tidyverse)   
library(viridis)  
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)  
library(tidyr)
library(ggpubr)     
library(rstatix)     
library(cowplot)
library(forcats)
library(here)

# ================================================================
# Read in data
# ================================================================

q2e_vals_per_sample = read_csv(here('data/q2e_per_sample.csv'))
q2e_vals_per_sample_filtered = read_csv(here('data/q2e_per_sample_filtered.csv'))

mds_data = read_csv(here('data/mds_data.csv')) %>%
  mutate(layer = factor(layer, levels=sort(unique(layer))))
q2e_data = read_csv(here('data/q2e_data.csv')) %>%
  mutate(layer = factor(layer, levels=sort(unique(layer))))
q2e_mds_data = read_csv(here('data/q2e_mds_data.csv'))

# ================================================================
# Heat-map of row counts after filtering (before counts in parentheses)
# ================================================================
df_before <- q2e_vals_per_sample
df_after  <- q2e_vals_per_sample_filtered


needed <- c("zooMS_ID", "pep_number")
stopifnot(all(needed %in% names(df_before)),
          all(needed %in% names(df_after)))

counts_before <- df_before %>%
  count(zooMS_ID, pep_number, name = "count_before")

counts_after <- df_after %>%
  count(zooMS_ID, pep_number, name = "count_after")


df <- counts_before %>%
  full_join(counts_after, by = c("zooMS_ID", "pep_number")) %>%
  mutate(
    count_before = as.integer(replace_na(count_before, 0)),
    count_after  = as.integer(replace_na(count_after,  0))
  ) %>%
  arrange(zooMS_ID, pep_number)


df <- df %>%
  mutate(
    retention  = count_after / if_else(count_before == 0, NA_real_, count_before),
    label_text = paste0(count_after, " (", count_before, ")")
  )


tail_species <- c("Canidae", "Crocuta/Panthera", "Ursidae")


other_species <- df %>%
  distinct(zooMS_ID) %>%
  pull(zooMS_ID) %>%
  setdiff(tail_species)
other_species <- other_species[order(tolower(other_species))]

zoo_levels <- c(other_species, tail_species)


df <- df %>%
  mutate(
    zooMS_ID   = factor(zooMS_ID, levels = zoo_levels),
    pep_number = factor(pep_number, levels = sort(unique(pep_number)))
  )


df <- df %>%
  mutate(text_col = if_else(replace_na(retention, 0) > 0.6, "black", "white"))


ggplot(df, aes(x = pep_number, y = zooMS_ID, fill = retention)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = label_text, colour = text_col),
            size = 5, fontface = "bold") +
  scale_fill_viridis(
    name   = "Filtering\n(after / before)",
    option = "C",
    limits = c(0, 1),
    na.value = "#f0f0f0"
  ) +
  scale_colour_identity() +
  # Reverse y-axis to match zoo_levels top-to-bottom
  scale_y_discrete(limits = rev(zoo_levels)) +
  labs(
    x = "pep_number",
    y = "zooMS_ID",
    title = "Row counts after filtering (before counts in parentheses)"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 13),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 13),
    panel.grid  = element_blank(),
    plot.title  = element_text(size = 10, face = "bold", hjust = 0.5)
  )


# ================================================================
# wls q2e all species all layer facet grid
# ================================================================
q2e_wide <- q2e_mds_data %>%
  select(sample, pep_number, q2e) %>%
  pivot_wider(names_from = pep_number, values_from = q2e)


pep_cols <- setdiff(names(q2e_wide), "sample")
pep_cols <- pep_cols[order(parse_number(pep_cols), na.last = TRUE)]
q2e_wide <- q2e_wide %>% select(sample, all_of(pep_cols))


data <- mds_q2e_merged <- mds_data %>% left_join(q2e_wide, by = "sample")


# Select relevant columns for the plot
peptide_columns <- c('pep1', 'pep2', 'pep3', 'pep4', 'pep5', 'pep6', 'pep7')

# Reshape the data for easier plotting


long_data <- data %>%
  select(layer, sample, batch, plate, Q2E, zooMS_ID, all_of(peptide_columns)) %>%
  pivot_longer(cols = all_of(peptide_columns), names_to = "peptide", values_to = "value") %>%
  drop_na(value)

# Reorder zooMS_ID to move "Crocuta/Panthera" and "Ursidae" to the end
long_data$zooMS_ID <- factor(long_data$zooMS_ID, 
                             levels = c(setdiff(unique(long_data$zooMS_ID), 
                                                c('Canidae', "Crocuta/Panthera", "Ursidae")), 
                                        'Canidae', "Crocuta/Panthera", "Ursidae"))

# Set the layer order for a consistent x-axis
long_data$layer <- factor(long_data$layer, levels = sort(unique(data$layer)))

# Set the peptide order as requested
long_data$peptide <- factor(long_data$peptide, levels = c('pep1', 'pep2', 'pep3', 'pep4', 'pep5', 'pep6', 'pep7'))

# Function to modify facet labels
facet_labeller = function(variable, value) {
  if (variable=='zooMS_ID') {
    return(str_replace_all(value, '/', '\n'))
  } else {
    return(value)
  }
}

# Generate the plot with both violin and scatter overlay
long_data %>% dplyr::filter(!peptide %in% c('pep6', 'pep7')) %>%
  ggplot(aes(x = layer, y = value, fill = peptide)) +
    geom_violin(scale = "width", alpha = 0.4, drop=T) +
    geom_jitter(aes(color = peptide), width = 0.2, alpha = 0.3) +
    geom_smooth(aes(x=as.numeric(layer), y=value), alpha=0.5) +
    stat_summary(fun = median, geom = "point", size = 2, colour = "black", alpha=0.5) +
    facet_grid(zooMS_ID ~ peptide, scales = "fixed", labeller=facet_labeller) +  # 只保留一个 facet
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(
      title = "Violin + Jitter of q2e by Layer and Peptide",
      x = "Layer",
      y = "WLS q2e"
    ) +
    scale_y_continuous(limits = c(0, 1))


long_data %>% dplyr::filter(peptide %in% c('pep6', 'pep7')) %>%
  ggplot(aes(x = layer, y = value, fill = peptide)) +
  geom_violin(scale = "width", alpha = 0.4, drop=T) +
  geom_jitter(aes(color = peptide), width = 0.2, alpha = 0.3) +
  geom_smooth(aes(x=as.numeric(layer), y=value), alpha=0.5) +
  stat_summary(fun = median, geom = "point", size = 2, colour = "black", alpha=0.5) +
  facet_grid(zooMS_ID ~ peptide, scales = "fixed", labeller=facet_labeller) +  # 只保留一个 facet
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "Violin + Jitter of q2e by Layer and Peptide",
    x = "Layer",
    y = "WLS q2e"
  ) +
  scale_y_continuous(limits = c(0, 1))


# ================================================================
# wls q2e only by species violin
# ================================================================

# Create a plot ignoring the layer and breaking up by zooMS_ID
ggplot(long_data, aes(x = peptide, y = value, fill = peptide)) +
  geom_violin(scale = "width", alpha = 0.3) + # Violin plot for distribution
  geom_jitter(aes(color = peptide), width = 0.15, alpha = 0.2) + # Add jitter for points
  stat_summary(fun = median, geom = "point", size = 2, colour = "black") +
  facet_wrap(~ zooMS_ID, scales = "free") + # Break up by zooMS_ID
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(title = "Peptide Values by zooMS_ID",
       x = "Peptide",
       y = "WLS q2e") +
  scale_y_continuous(limits = c(0, 1))


# ================================================================
# equidae
# ================================================================
# Split the data: one excluding Equidae, one including only Equidae
long_data <- data %>%
  select(layer, sample, batch, plate, Q2E, zooMS_ID, all_of(peptide_columns)) %>%
  pivot_longer(cols = all_of(peptide_columns), names_to = "peptide", values_to = "value") %>%
  drop_na(value)


data_excl_equidae <- long_data %>% filter(zooMS_ID != "Equidae")
data_equidae_only <- long_data %>% filter(zooMS_ID == "Equidae")




# ---- PLOT 1: Only Equidae by batch ----

pep4_filtered <- data_equidae_only %>%
  filter(peptide == "pep4" & !is.na(value))

# Add a category column
pep4_filtered <- pep4_filtered %>%
  mutate(category = ifelse(value > 0.5, "Above 0.5", "Below 0.5"))

# Arrange data by layer and value
pep4_filtered <- pep4_filtered %>%
  arrange(layer, value)

# Print table
print(pep4_filtered)

# Optionally save the table as a CSV file
write.csv(pep4_filtered, "pep4_above_below_0.5.csv", row.names = FALSE)


pep4_filtered$batch <- as.factor(pep4_filtered$batch)

# Generate the scatter plot
pep4_equidae_plot = ggplot(pep4_filtered, aes(x = layer, y = value, color = batch)) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 5) +  # Jitter to avoid overlap
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "pep4 of Equidae",
       x = "Layer",
       y = "WLS q2e",
       color = "Batch") +
  scale_color_brewer(palette = "Set1")  # Use a color-friendly palette




# ---- PLOT 2: Only Equidae grouping mapping back to all peptides ----

color_coding <- pep4_filtered %>%
  select(sample, value) %>%
  mutate(category = ifelse(value > 0.5, "Above 0.5", "Below 0.5"))

# Step 2: Merge color categories into the full dataset
data_colored <- data_equidae_only %>%
  left_join(color_coding, by = "sample") %>%
  mutate(category = replace_na(category, "Other"))  # Fill missing categories

# Debugging: Check if value column exists after merge
print(colnames(data_colored))

# Step 3: Ensure 'category' is a factor
data_colored$category <- as.factor(data_colored$category)

# Step 4: Plot all peptides with consistent coloring from pep5
wls_q2e_vs_layer_all_pep = ggplot(data_colored, aes(x = layer, y = value.x, color = category)) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 3) +  # Scatter plot points
  facet_wrap(~ peptide) +  # Facet for each peptide
  scale_color_manual(values = c("Above 0.5" = "red", "Below 0.5" = "blue", "Other" = "gray")) +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Equidae",
       x = "Layer",
       y = "WLS q2e",
       color = "pep4 Category")+
  scale_y_continuous(limits = c(0, 1))


wls_q2e_vs_layer_combined <- pep4_equidae_plot/wls_q2e_vs_layer_all_pep + plot_annotation(tag_levels = 'A')
wls_q2e_vs_layer_combined


# ------------------------------------------------
# MDS 
# ------------------------------------------------

x_labels = levels(mds_data$layer)
x_lab_alt = levels(mds_data$layer)
x_lab_alt[seq(2, 10, 2)] = paste0('\n', x_labels[seq(2, 10, 2)])


tail_species <- c("Canidae", "Crocuta/Panthera", "Ursidae")

mds_data$zooMS_ID <- factor(
  mds_data$zooMS_ID,
  levels = c(
    sort(setdiff(unique(mds_data$zooMS_ID), tail_species)),  # 其余物种 A→Z
    tail_species                                             # 固定尾部 3 个
  )
)

#MDS vs GA q2e
ggplot(mds_data) +
  geom_point(aes(x = GA_q2e, y = MDS.Model, color = zooMS_ID), alpha = 0.4) +
  geom_smooth(method = "lm", aes(x = GA_q2e, y = MDS.Model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") + 
  xlab('GA q2e') + 
  ylab('MALDI Deamidation Score') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8))


#wls q2e vs GA q2e
ggplot(q2e_data) +
  geom_point(aes(x=GA_q2e, y=q2e, color=zooMS_ID)) +
  geom_smooth(aes(x=GA_q2e, y=q2e, color=zooMS_ID)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  ylab('WLS q2e') + xlab('GA q2e') +
  ylim(0,1) + xlim(0,1) +
  facet_wrap(~pep_number) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size=8))



#MDS vs layer
ggplot(mds_data) +
  geom_jitter(aes(x=as.numeric(layer), y=MDS.Model, fill=layer), shape=21, alpha=0.3) +
  geom_violin(aes(x=as.numeric(layer), y=MDS.Model, fill=layer, color=layer), scale='width', alpha=0.1) +
  geom_smooth(aes(x=as.numeric(layer), y=MDS.Model)) +
  scale_x_discrete(limits=levels(mds_data$layer), labels=x_labels, name='layer') +
  ylab('MALDI Deamidation Score') +
  theme_bw() +
  theme(axis.title.x.bottom = element_text(size=16),
        axis.text.x.bottom = element_text(size=16))


#GA q2e vs layer 
ggplot(mds_data) +
  geom_jitter(aes(x=layer, y=GA_q2e, fill=layer), shape=21, alpha=0.3) +
  geom_violin(aes(x=layer, y=GA_q2e, fill=layer, color=layer), scale='width', alpha=0.1) +
  geom_smooth(aes(x=as.numeric(layer), y=GA_q2e)) +
  scale_x_discrete(limits=levels(mds_data$layer), labels=x_lab_alt, name='layer') +
  scale_y_continuous(minor_breaks = NULL) +
  ylab('GA q2e') +
  theme_bw() +
  theme(axis.title.x.bottom = element_text(size=16),
        axis.text.x.bottom = element_text(size=12),
        strip.text.x = element_text(size=12))







# ------------------------------------------------
# layer 11.2-15, WLS
# ------------------------------------------------


df <- mds_q2e_merged

layers   <- c(11.2, 11.3, 11.4, 12, 13, 14, 15)   
peptides <- c("pep1", "pep2")

id_by <- "Bison/Yak"
id_cg <- "Cervidae/Gazella/Saiga"
id_eq <- "Equidae"  

plot_df_all <- df |>
  dplyr::filter(layer %in% layers,
                zooMS_ID %in% c(id_by, id_cg, id_eq)) |>
  tidyr::pivot_longer(all_of(peptides), names_to = "peptide", values_to = "value") |>
  tidyr::drop_na(value) |>
  dplyr::mutate(
    layer    = factor(layer, levels = layers, ordered = TRUE),
    peptide  = factor(peptide, levels = peptides),
    zooMS_ID = factor(zooMS_ID, levels = c(id_by, id_cg, id_eq))
  )


y_limits_pep1 <- c(-0.05, 1.05)  
y_limits_pep2 <- c(-0.05, 0.85)  


make_violin <- function(dat, pep, y_limits) {
  
  dat_pep_plot <- dat |> dplyr::filter(peptide == pep)                    
  dat_pep_stat <- dat_pep_plot |> dplyr::mutate(layer = forcats::fct_drop(layer))  
  do_stats <- nlevels(dat_pep_stat$layer) >= 2
  
  kw_p <- if (do_stats) rstatix::kruskal_test(dat_pep_stat, value ~ layer)$p else NA_real_
  
  pw_sig <- tibble::tibble()
  if (do_stats) {
    lvl <- levels(dat_pep_stat$layer)
    adj_pairs <- tibble::tibble(group1 = head(lvl, -1), group2 = tail(lvl, -1))
    
    pw <- rstatix::pairwise_wilcox_test(dat_pep_stat, value ~ layer,
                                        p.adjust.method = "bonferroni") |>
      dplyr::semi_join(adj_pairs, by = c("group1","group2")) |>
      dplyr::mutate(
        p.signif = dplyr::case_when(
          p.adj < 0.001 ~ "***",
          p.adj < 0.01  ~ "**",
          p.adj < 0.05  ~ "*",
          TRUE          ~ "ns"
        )
      )
    
    pw_sig <- pw |> dplyr::filter(p.adj < 0.05)
    if (nrow(pw_sig) > 0) {
      y_base <- y_limits[2] - 0.05   
      step   <- 0.02                
      pw_sig <- pw_sig |>
        dplyr::mutate(y.position = y_base - step * (dplyr::row_number() - 1))
    }
  }
  
  p <- ggplot2::ggplot(dat_pep_plot, ggplot2::aes(layer, value, fill = layer)) +
    ggplot2::geom_violin(alpha = .25, width = .9, trim = FALSE, colour = NA) +
    ggplot2::geom_boxplot(width = 0.20, outlier.shape = NA,
                          colour = "#555555", fill = "white", alpha = 0.9) +
    ggplot2::geom_text(
      data = dat_pep_plot |> dplyr::add_count(layer) |> dplyr::distinct(layer, n),
      ggplot2::aes(layer, y_limits[1] + 0.04 * diff(y_limits), label = paste0("n=", n)),
      inherit.aes = FALSE
    ) +
    ggplot2::coord_cartesian(ylim = y_limits, expand = FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE) +             
    ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    scale_fill_viridis_d(
      option = "cividis", begin = 0.1, end = 0.95,
      limits = levels(plot_df_all$layer), drop = FALSE
    ) +
    ggplot2::labs(
      x = "Layer", y = "WLS q2e",
      title    = paste0(pep, " — ", unique(dat_pep_plot$zooMS_ID)),
      subtitle = if (is.na(kw_p)) "Kruskal-Wallis: NA )"
      else paste0("Kruskal-Wallis p = ", formatC(kw_p, format = "e", digits = 2))
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  if (nrow(pw_sig) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(
      pw_sig, label = "p.signif",
      inherit.aes = FALSE, tip.length = 0.01, size = 4
    )
  }
  p
}


dat_BY <- plot_df_all |> dplyr::filter(zooMS_ID == id_by)
dat_CG <- plot_df_all |> dplyr::filter(zooMS_ID == id_cg)
dat_EQ <- plot_df_all |> dplyr::filter(zooMS_ID == id_eq)


p11 <- make_violin(dat_BY, "pep1", y_limits_pep1)  
p12 <- make_violin(dat_BY, "pep2", y_limits_pep2) 
p21 <- make_violin(dat_EQ, "pep1", y_limits_pep1)  
p22 <- make_violin(dat_EQ, "pep2", y_limits_pep2)  
p31 <- make_violin(dat_CG, "pep1", y_limits_pep1)  
p32 <- make_violin(dat_CG, "pep2", y_limits_pep2)  


p12 <- p12 + theme(axis.title.y = element_blank())
p22 <- p22 + theme(axis.title.y = element_blank())
p32 <- p32 + theme(axis.title.y = element_blank())


leg <- cowplot::get_legend(
  ggplot2::ggplot(plot_df_all, ggplot2::aes(layer, value)) +
    ggplot2::geom_violin(ggplot2::aes(fill = layer), alpha = .6, colour = NA) +
    scale_fill_viridis_d(
      option = "cividis", begin = 0.1, end = 0.95,
      limits = levels(plot_df_all$layer), drop = FALSE
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position   = "right",
      legend.box.margin = ggplot2::margin(0, 6, 0, 6, unit = "mm"),
      plot.margin       = ggplot2::margin(0, 0, 0, 0, unit = "mm")
    )
)


p11_nl <- p11 + theme(legend.position = "none")
p12_nl <- p12 + theme(legend.position = "none")
p21_nl <- p21 + theme(legend.position = "none")
p22_nl <- p22 + theme(legend.position = "none")
p31_nl <- p31 + theme(legend.position = "none")
p32_nl <- p32 + theme(legend.position = "none")


left_grid <- cowplot::plot_grid(
  p11_nl, p12_nl,
  p21_nl, p22_nl,
  p31_nl, p32_nl,
  ncol = 2, labels = c("A","B","C","D","E","F"), align = "hv"
)

final_fig <- cowplot::plot_grid(
  left_grid, leg, ncol = 2, rel_widths = c(1, 0.28)
)
final_fig






















