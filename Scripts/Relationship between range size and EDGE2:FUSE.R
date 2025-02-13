
#Range size of each species 
load(here::here("Data","GAP analyses","grids_commonspecies_corrected_021024.Rdata"))

# Extract latitude from row names
lat <- as.numeric(sapply(strsplit(rownames(grids.fd_new), "_"), "[", 1))

# Calculate cell area for each latitude (in km²)
cell_areas <- 111.32 * 0.5 * # longitude dimension
  111.32 * 0.5 * # latitude dimension
  cos(lat * pi/180) # correction for latitude

# Calculate range size in km²
range_sizes_km2 <- colSums(grids.fd_new > 0 * cell_areas)

#Create a dataframe
range_size_df <- data.frame(
  Species = names(range_sizes_km2),
  Range_size = range_sizes_km2
)     

saveRDS(range_size_df, file = here::here("Data", "species_range_sizes.rds"))

# Open the EDGE2 and FUSE data and create dataframes
file_name <- here::here("Data", "My dataframes", "shark_conservation_metrics_no_freshwater.rds")
results <- readRDS(file = file_name)

EDGE2_df=cbind(results$EDGE2$info$Species,
         results$EDGE2$transformed)
colnames(EDGE2_df)="Species"

FUSE_df=cbind(results$FUSE$info$Species,
        results$FUSE$transformed)
colnames(FUSE_df)="Species"

# Create a complete dataframe
conservation_df <- merge(range_size_df, EDGE2_df[,c(1,3)], by = "Species", all = TRUE)
conservation_df <- merge(conservation_df, FUSE_df[,c(1,3)], by = "Species", all = TRUE)

# Rename columns for clarity
colnames(conservation_df)[2:4] <- c("Range_size", "EDGE2", "FUSE")

# 1. Correlation tests
# Spearman correlation (doesn't assume normality)
cor_edge <- cor.test(log10(conservation_df$Range_size), conservation_df$EDGE2, 
                     method = "spearman", use = "complete.obs")
cor_fuse <- cor.test(log10(conservation_df$Range_size), conservation_df$FUSE, 
                     method = "spearman", use = "complete.obs")

# Print correlation results
print("Correlation between Range Size and EDGE2:")
print(cor_edge)
print("Correlation between Range Size and FUSE:")
print(cor_fuse)

# 2. Visualizations using ggplot2
library(ggplot2)
library(gridExtra)

# Create scatter plot for EDGE2
p1 <- ggplot(conservation_df, aes(x = Range_size, y = EDGE2)) +
  geom_point(alpha = 0.6) +
  #geom_smooth(method = "lm", color = "blue") +
  scale_x_log10() +
  theme_bw() +
  labs(x = "Range Size (km²)", 
       y = "EDGE2 Score",
       title = "Range Size vs EDGE2") +
  theme(plot.title = element_text(hjust = 0.5))
p1

# Create scatter plot for FUSE
p2 <- ggplot(conservation_df, aes(x = Range_size, y = FUSE)) +
  geom_point(alpha = 0.6) +
 # geom_smooth(method = "lm", color = "red") +
  scale_x_log10() +
  theme_bw() +
  labs(x = "Range Size (km²)", 
       y = "FUSE Score",
       title = "Range Size vs FUSE") +
  theme(plot.title = element_text(hjust = 0.5))
p2
# Arrange plots side by side
grid.arrange(p1, p2, ncol = 2)
