#textwidth = 17.5 # cm

# sizes (linewidth, point size)
lsz = 1.2
psz = 1.0
erlsz = 5*lsz/6
erlwth = lsz/12

# jitter spacings
jtwdt = 0.15
jthgt = 0

# alpha
alp =  0.3

# Condition colors
cond_cols <- scale_color_manual(values = c("orange", "blue"), labels = c("External", "Action"))

# define a 'theme'
mytheme	 <- theme_classic() + theme(
  # Set font
  text = element_text(family = "Microsoft Sans Serif"),
  # Set axis text size
  axis.text = element_text(size = rel(0.9)),
  axis.title = element_text(size = rel(1.1)),
  # Set title elements
  plot.title = element_text(size = rel(1.1), hjust = 0.5),
  plot.tag = element_text(size = rel(1.5)),
  strip.text = element_text(size = rel(1.1)),
  strip.background = element_blank(),
  # Set legend text size
  legend.title = element_text(size = rel(0.9)),
  legend.text = element_text(size = rel(0.8)),
  # Set legend spacing and position
  legend.box.spacing = unit(0, "pt"),
  legend.position = "bottom",
  legend.key = element_blank(),
  # Remove grid and background
  panel.grid = element_blank(),
  panel.background = element_blank(),
  axis.line = element_blank())
