# load ggplot2 or tidyverse
library(tidyverse)

# set the working directory
setwd("~/Dropbox/Teaching/Syracuse_University/Applied_Genomics/R_data/")

# Read in the annual CPI data (get failiar with the columns and what's in the data)
cpi <- read.csv(file = "CPI_data.txt", header = T, sep = "\t")

# Read in a different formt of this data
cpi_by_year <- read.csv("CPI_year_data.csv", header = T)

# Read in the CPI and HDI data from 2017
econ <- read.csv(file = "Rgraphics/dataSets/EconomistData.csv", sep = ",", header = T)

# Plot CPI against HDI 
ggplot(data = econ, aes(x = CPI, y = HDI)) +
  geom_point()

# ... several ways to write the above, e.g.:
ggplot(data = econ) +
  geom_point(aes(x = CPI, y = HDI))

# ..or:
ggplot() +
  geom_point(data = econ, aes(x = CPI, y = HDI))

###-----------------------------------------###

# Plot CPI against HDI and map "Region" to colour 
ggplot(data = econ, aes(x = CPI, y = HDI)) +
  geom_point(aes(colour = Region))

# Add a second layer, such as a fitted line:
ggplot(data = econ, aes(x = CPI, y = HDI)) +
  geom_point(aes(colour = Region)) +
  geom_smooth(method = "lm", formula = y ~ log(x))

# Add a third layer, such as some country names. First, make a vector of countries:
myCountries <- c("United States", "Singapore", "Brazil", "Russia")

ggplot(data = econ, aes(x = CPI, y = HDI)) +
  geom_point(aes(colour = Region)) +
  geom_smooth(method = "lm", formula = y ~ log(x)) +
  geom_text(data = subset(econ, Country %in% myCountries), aes(label = Country))

###----------------------------------------###

# select the HDI column and merge it with the annual data:
hdi <- select(econ, Country, HDI)
cpi_with_hdi <- merge(cpi, hdi, by.x = "Country", by.y = "Country", all.y = TRUE)

###----------------------------------------###
ggplot(cpi_with_hdi) +
  geom_point(aes(CPI_Score_2012, HDI), colour = "blue", alpha = 0.5) +
  geom_point(aes(CPI_score_2019, HDI), colour = "red", alpha = 0.5) +
  geom_smooth(aes(CPI_Score_2012, HDI), method = "lm", formula = y ~ log(x), colour = "blue") +
  geom_smooth(aes(CPI_score_2019, HDI), method = "lm", formula = y ~ log(x), colour = "red")

###----------------------------------------###

# Now let's try other geoms:

# bargraph:
ggplot(filter(cpi_with_hdi, Region == "AME"), aes(reorder(Country, CPI_score_2019), CPI_score_2019))+
  geom_bar(stat = "identity", fill = "green") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_errorbar(mapping = aes(ymin = CPI_score_2019 - Standard_error_2019, ymax = CPI_score_2019 + Standard_error_2019), width = 0.3, colour = "purple") +
  labs(x = "Country", y = "Corruption Index")

# boxplot

ggplot(data = cpi_by_year) +
  geom_boxplot(aes(as.character(Year), CPI_score)) +
  geom_jitter(aes(as.character(Year), CPI_score, colour = Region), width = 0.1)

###-----------------------------------------###

# Faceting is very useful for splitting data into panels
cpi_by_year_with_hdi <- merge(cpi_by_year, hdi, by.x = "Country", by.y = "Country", all.x = T)

# facet_wrap
ggplot(cpi_by_year_with_hdi, aes(CPI_score, HDI)) +
  geom_point() +
  facet_wrap(~Year)

# facet_grid
ggplot(cpi_by_year_with_hdi, aes(CPI_score, HDI)) +
  geom_point() +
  facet_grid(Region~Year) 

###-----------------------------------------###

# Position adjustments with barplots

ggplot(data = diamonds) + 
  geom_bar(mapping = aes(x = cut, fill = color), position = "dodge")

###-----------------------------------------###

# Coordinate systems can be flipped

ggplot(data = diamonds) + 
  geom_bar(mapping = aes(x = cut, fill = color), position = "dodge") +
  coord_flip()