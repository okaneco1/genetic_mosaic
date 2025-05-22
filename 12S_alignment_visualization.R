# 12S Great Lakes Fish Alignment Visual

# libraries
library(Biostrings)
library(msa)
library(tidyverse)


# Load in fasta
nt_sequences <- readDNAStringSet("./12S_no_human.fa")

# Convert to character vectors
seq_df <- data.frame(
  full_name = names(nt_sequences),  
  sequence = as.character(nt_sequences)
)

# species names for axis display
seq_df <- seq_df %>%
  mutate(species_label = sub("^[^_]+_", "", full_name))  


# split sequences into rows, label with nucleotide numbers and ungroup
seq_df_long <- seq_df %>%
  mutate(bases = strsplit(sequence, "")) %>%
  unnest_longer(bases) %>%
  group_by(full_name) %>%
  mutate(position = row_number()) %>%
  ungroup()

# merge back species names
seq_df_long <- left_join(seq_df_long, seq_df %>% select(full_name), by = "full_name")


# color palette for nucleotides
base_colors <- c(
  A = "#26547c",
  C = "#ef476f",
  G = "#ffd166",
  T = "#06d6a0",
  `-` = "white"  # gaps can be background color
)

# add color info
seq_df_long <- seq_df_long %>%
  mutate(color = base_colors[bases])

# plot
ggplot(seq_df_long, aes(x = position, y = fct_rev(full_name))) +
  geom_point(aes(fill = bases), shape = 21, size = 3.2, color="white", stroke=0.4) +
  scale_fill_manual(values = base_colors) +
  # need to perform label switch as some species are repeated
  scale_y_discrete(labels = setNames(seq_df$species_label, seq_df$full_name)) +
  theme_minimal() +
  coord_fixed(ratio = 1) +  
  labs(x = "Alignment Position", y = "Species") +
  theme(axis.text.y = element_text(size = 4),
        panel.grid = element_blank())

# there is some trial and error involved with getting the right dimensions for
# the plot you save, and how the points are spaced. The fixed ratio keeps an equal
# distance, but even spacing requires some messing around with image dimensions

# save the plot
ggsave("./12S_alignment_visualization.png", dpi=600, height = 24, width = 18)

# minor adjustments to the axes, additional titles, and layout of the design
# were then made in Adobe Illustrator with this file.


#-------------
# small bit of code to export names for proper x axis label spacing in Illustrator
clean_names <- rev(names(nt_sequences)) %>%
  sub("^[^_]+_", "", .) %>%     # Remove identifier before underscore
  gsub("_", " ", .)             # Replace underscore with space
writeLines(clean_names, "./species_labels.txt")

