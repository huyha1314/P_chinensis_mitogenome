# 1. Install required packages (run this once if you don't have them)
# install.packages(c("ggplot2", "sf", "rnaturalearth", "rnaturalearthdata", "ggrepel"))

# 2. Load the libraries
library(ggplot2)
library(sf)
library(rnaturalearth)
library(ggrepel)

# 3. Get the base map data (High resolution)
world_map <- ne_countries(scale = "medium", returnclass = "sf", continent = "Asia")

# 4. Create a dataframe with your EXACT Parashorea chinensis locations
sample_data <- data.frame(
  Sample = c("VN1 (Thai Nguyen)", "CH1 (Yunnan)", "CH2 (Hainan)"),
  Location = c("Thai Nguyen, Vietnam", "Yunnan, China", "Hainan, China"),
  Lon = c(106.010565, 101.255000, 108.783333), # Exact Longitudes
  Lat = c(21.843256, 21.922000, 18.700000),    # Exact Latitudes
  Type = c("Mainland", "Mainland", "Island") 
)

# Convert your dataframe into a spatial 'sf' object so ggplot understands it
sample_points <- st_as_sf(sample_data, coords = c("Lon", "Lat"), crs = 4326)

# 5. Plot the map!
my_map <- ggplot(data = world_map) +
  # Draw the country borders (light gray fill, white borders)
  geom_sf(fill = "#e0e0e0", color = "white", size = 0.3) +
  
  # Add your sample points, colored by whether they are Mainland or Island
  geom_sf(data = sample_points, aes(color = Type), size = 4, shape = 19) +
  
  # Add labels for the points, automatically repelled so they don't overlap
  geom_text_repel(data = sample_data, aes(x = Lon, y = Lat, label = Sample), 
                  fontface = "bold", size = 4.5, nudge_y = -0.5) +
  
  # Set the specific zoom box to perfectly frame the 3 samples
  coord_sf(xlim = c(98, 112), ylim = c(17, 24), expand = FALSE) +
  
  # Use theme_void() to remove all axes, grids, and numbers for a clean diagram look
  theme_void() +
  theme(
    panel.background = element_rect(fill = "#e0f3f8", color = NA), # Light blue ocean
    legend.position = "bottom",
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank() # Hide the word "Type" above the legend
  ) +
  # Set distinct colors: Green for Mainland stability, Orange for Island variation
  scale_color_manual(values = c("Island" = "#d95f02", "Mainland" = "#1b9e77"))

# View the map in R
print(my_map)

# 6. Export as a Vector PDF
ggsave("results/Fig1_Parashorea_Geographic_Map_Final.pdf", plot = my_map, width = 6, height = 5, units = "in", device = "pdf", bg = "white")

# 7. Export as a high-resolution PNG
ggsave("results/Fig1_Parashorea_Geographic_Map_Final.png", plot = my_map, width = 6, height = 5, units = "in", dpi = 600, bg = "white")

# 8. Export as a TIFF
ggsave("results/Fig1_Parashorea_Geographic_Map_Final.tiff", plot = my_map, width = 174, height = 145, units = "mm", dpi = 600, compression = "lzw", bg = "white")