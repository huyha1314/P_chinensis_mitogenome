library(ggplot2)
p <- ggplot() + 
  annotate("text", x=1, y=1, label="Hello World Shorea", size=10) +
  theme_minimal()
dir.create("scratch", showWarnings = FALSE)
ggsave("scratch/test_default.png", p)
ggsave("scratch/test_cairo.png", p, type="cairo")
ggsave("scratch/test_png_device.png", p, device=png)
ggsave("scratch/test_ragg.png", p, device=ragg::agg_png)
cat("Done testing fonts!\n")
