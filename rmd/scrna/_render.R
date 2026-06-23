# Knit the Rmd to produce the SVG figures (no pandoc needed for figures).
setwd("C:/Users/HagayLadany/Desktop/webpage/rmd/scrna")
knitr::knit("scrna.Rmd", output = "scrna.md", quiet = TRUE)
cat("DONE\n")
list.files("figs")
