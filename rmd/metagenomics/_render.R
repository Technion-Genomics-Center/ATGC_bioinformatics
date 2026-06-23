# Knit the Rmd to produce the SVG figures (no pandoc needed for figures).
setwd("C:/Users/HagayLadany/Desktop/webpage/rmd/metagenomics")
knitr::knit("metagenomics.Rmd", output = "metagenomics.md", quiet = TRUE)
cat("DONE\n")
list.files("figs")
