# Knit the Rmd to produce the SVG figures (no pandoc needed for figures).
setwd("C:/Users/HagayLadany/Desktop/webpage/rmd/spatial")
knitr::knit("spatial.Rmd", output = "spatial.md", quiet = TRUE)
cat("DONE\n")
list.files("figs")
