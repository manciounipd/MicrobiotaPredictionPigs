library(ggplot2)
library(dplyr)
library(tidyr)
library(ggnewscale)

COLOR_MAPPING <- c(
  # Inner ring colors
  "m[y]"   = "#a6cee3",   # pastel blue
  "e[y]"   = "#fb9a99",   # pastel pink/red
  "g[y]"   = "#b2df8a",   # pastel green
  "sys[y]" = "#fdbf6f",   # pastel orange
  
  # Outer ring colors (matched to inner)
  "g[m]"   = "#b2df8a",   # pastel green
  "res[m]" = "#fb9a99",   # pastel pink/red
  "sys[m]" = "#fdbf6f"    # pastel orange
)

ALL_EFFECTS <- names(COLOR_MAPPING)

# 2. MAKE_PLOT function
MAKE_PLOT <- function(donut_df, MODEL){
  
  inner_labels <- donut_df %>% 
    dplyr::filter(level == 1) %>% 
    dplyr::distinct(label2) %>% 
    dplyr::mutate(display = label2)
  
  outer_labels <- donut_df %>% 
    dplyr::filter(level == 2) %>% 
    dplyr::distinct(label2) %>% 
    dplyr::mutate(display = label2)
  
  ggplot() +
    # Inner ring
  
  geom_rect(
    data = donut_df %>% filter(level == 1, val > 0),
    aes(xmin = 1, xmax = 2, ymin = ymin, ymax = ymax, fill = label2),
    color = "white", linewidth = 0.3
  ) +
  scale_fill_manual(
    name = "Inner",
    values = COLOR_MAPPING,
    labels = function(x) parse(text = x)   # This converts m[y], g[m] to proper legend
  ) +
    new_scale_fill() +
    
    # Outer ring
    geom_rect(
      data = donut_df %>% filter(level == 2, val > 0),
      aes(xmin = 2, xmax = 2.5, ymin = ymin, ymax = ymax, fill = label2),
      color = "white", linewidth = 0.3
    ) +
    scale_fill_manual(
        name = "Outer",
    values = COLOR_MAPPING,
    labels = function(x) parse(text = x)   # This converts m[y], g[m] to proper legend
  ) +
    new_scale_fill() +
    
    # Segment labels
    geom_text(
      data = donut_df %>% filter(val > 0.04),
      aes(x = ifelse(level == 1, 1.5, 2.25), 
          y = label_pos, 
          label = paste0(round(val*100, 0), "%")),
      color = "black", size = 5.3
    ) +
    
    # Center label
    annotate("text", x = 0, y = 0, label = MODEL, size = 6, fontface = "italic") +
    
    coord_polar(theta = "y") +
    scale_x_continuous(limits = c(0, 2.5)) +
    theme_void() +
    theme(
          plot.margin = unit(c(0, 0, 0, 0), "null"),
    panel.margin = unit(c(0, 0, 0, 0), "null"),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 8)
    )
}

# 3. GET_DB_X_PLOT function
GET_DB_X_PLOT <- function(df) {
  
  inner <- df %>%
    dplyr::group_by(subeffect) %>%
    dplyr::summarise(val = sum(val), .groups = "drop") %>%
    dplyr::mutate(level = 1, label = subeffect)
  
  outer <- df %>%
    dplyr::filter(subeffect == "M") %>%
    dplyr::mutate(level = 2, label = effect)
  
  # Map inner labels
  inner <- inner %>%
    dplyr::mutate(label2 = dplyr::case_when(
      label == "F"   ~ "sys[y]",
      label == "RES" ~ "e[y]",
      label == "G"   ~ "g[y]",
      label == "M"   ~ "m[y]",
      TRUE ~ label
    ))
  
  # Map outer labels
  outer <- outer %>%
    dplyr::mutate(label2 = dplyr::case_when(
      effect == "G"   ~ "g[m]",
      effect == "F"   ~ "sys[m]",
      effect == "RES" ~ "res[m]",
      TRUE ~ effect
    ))
  
  donut_df <- dplyr::bind_rows(inner, outer) %>%
    # Ensure all effects are present
    tidyr::complete(label2 = ALL_EFFECTS, fill = list(val = 0, level = 1)) %>%
    # Fix factor order
    dplyr::mutate(label2 = factor(label2, levels = ALL_EFFECTS)) %>%
    # Compute positions
    dplyr::group_by(level) %>%
    dplyr::arrange(level, label2) %>%
    dplyr::mutate(
      ymax = cumsum(val),
      ymin = ymax - val,
      label_pos = (ymax + ymin) / 2
    ) %>%
    dplyr::ungroup()
  
  return(donut_df)
}
