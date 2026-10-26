
install.packages("tidyverse")
library(tidyverse)


# CREATING OBJECTS IN R --------------------------------------------------------

# arithmetic in R
3 + 5
12 / 7
2 * 3
5 - 3
2^3
2**3
sqrt(25)

# assignment operators
weight_kg <- c(55, 60)

55 -> weight_kg
  
weight_kg = 55

# arithmetic using objects
2.2 * weight_kg

# change an objects value
weight_kg <- 57.5

weight_lb <- 2.2 * weight_kg

# DATA TYPES --------------------------------------------------------------
# 6 basic data types

# 1. character
a <- "abc"
typeof(a)

# 2. numeric (real or decimal)
b <- 2
c <- 13.5
typeof(b)

# integer
d <- 2L
typeof(d)

# logical
e <- TRUE
f <- FALSE

# complex
g <- 1+4i

# Note: all these are atomic vectors (holds data of a single data type)


# DATA STRUCTURES ---------------------------------------------------------


######### vectors #################


weight_g <- c(50, 60, 65, 82) # a vector containing numbers
weight_g

animals <- c("mouse", "rat", "dog") # a vector containing characters
animals

# how many elements are in the vector
length(weight_g)
length(animals)

# get the structure of the objects
str(weight_g)
str(animals)

# adding elements to a vector using the c() function
weight_g <- c(weight_g, 90) # add to the end of the vector
weight_g <- c(30, weight_g) # add to the beginning of the vector
weight_g

# Subsetting vectors
animals <- c("mouse", "rat", "dog", "cat")
animals[2] # get the element in the second position

# Note: R indices start at 1

animals[c(3, 2)]

# Conditional subsetting
weight_g <- c(21, 34, 39, 54, 55)

weight_g[weight_g > 30 & weight_g < 50]

weight_g[weight_g <= 30 | weight_g == 55]

animals <- c("mouse", "rat", "dog", "cat")
animals[animals == "cat" | animals == "rat"] # returns both rat and cat

# test whether any elements are in the vector
animals %in% c("rat", "cat", "dog", "duck", "goat")

# subset based on  whether elements are present
animals[animals %in% c("rat", "cat", "dog", "duck", "goat")]

######## Lists ######################

# a list can have different data types
list_example <- list(1, "a", TRUE, 1+4i)
list_example

str(list_example)

list_example[[2]]


####################################
# Missing values
heights <- c(2, 4, 4, NA, 6)
mean(heights)
max(heights)
mean(heights, na.rm = TRUE)
max(heights, na.rm = TRUE)

## Extract those elements which are not missing values.
heights <- heights[!is.na(heights)]

na.omit(heights)

heights[complete.cases(heights)]


########### Factors ################
# for categorical data
data = c(1,2,2,3,1,2,3,3,1,2,3,3,1)
fdata = factor(data)
fdata

sex <- factor(c("male", "female", "female", "male"))

levels(sex)
nlevels(sex)


########## Dataframes##############
download.file(url = "https://ndownloader.figshare.com/files/2292169",
              destfile = "data/portal_data_joined.csv")

# loading data into R
surveys <- read.csv("data/portal_data_joined.csv")
View(surveys)

# View the first 6 records
head(surveys)

# structure of the dataframe
str(surveys)

nrow(surveys)
ncol(surveys)

summary(surveys)

# Indexing and subsetting data frames
# first element in the first column of the data frame (as a vector)
surveys[1, 1]

# first element in the 6th column (as a vector)
surveys[1, 6] 

# first column of the data frame (as a vector)
surveys[, 1]

# first column of the data frame (as a data.frame)
surveys[1]  

# first three elements in the 7th column (as a vector)
surveys[1:3, 7] 

# the 3rd row of the data frame (as a data.frame)
surveys[3, ] 

# equivalent to head_surveys <- head(surveys)
head_surveys <- surveys[1:6, ] 

weight <- surveys$weight

# FUNCTIONS ---------------------------------------------------------------

# calling a function
weight_kg <- sqrt(10)

# arguments of functions
round(3.14159, digits = 2)
round(digits = 2, x = 3.14159) # another way of specifying the arguments
args(round)

# getting help on a function
?round


# creating functions in R

fahrenheit_to_celsius <- function(temp_F) {
  
  temp_C <- (temp_F - 32) * 5 / 9
  
  return(temp_C)
}


# freezing point of water
fahrenheit_to_celsius(32)

# boiling point of water
fahrenheit_to_celsius(212)



###################### VISUALIZATION ##########################################

download.file("https://raw.githubusercontent.com/datacarpentry/R-genomics/gh-pages/data/Ecoli_metadata.csv", 
              "data/Ecoli_metadata.csv")

metadata <- read.csv("data/Ecoli_metadata.csv")

# get the genome size vector
genome_size <- metadata$genome_size

## Scatter plot 
plot(genome_size)

# change the shape of the data point
plot(genome_size, pch=8)

?plot

# Add title to your plot
plot(genome_size, pch=8, main="Scatter plot of genome sizes")

## Histogram
hist(genome_size)


## Boxplot
boxplot(genome_size ~ cit, metadata)

boxplot(genome_size ~ cit, metadata,  col=c("pink","purple", "darkgrey"),
        main="Average expression differences between celltypes", ylab="Expression")

### GGPLOT ---------------------------------------------------------------------
library(ggplot2)

# To initialize the basic graph structure: 
ggplot()

# scatterplot
ggplot(metadata) + 
  geom_point() # note what happens here


ggplot(metadata) +
  geom_point(aes(x = sample, y= genome_size))

ggplot(metadata) +
  geom_point(aes(x = sample, 
                 y = genome_size, 
                 color = generation, 
                 shape = cit), 
             size = rel(3.0)) +
  theme(axis.text.x = element_text(angle=45, hjust=1))


# Histogram
ggplot(metadata) +
  geom_bar(aes(x = genome_size), stat = "bin", binwidth=0.05)


# Boxplot
ggplot(metadata) +
  geom_boxplot(aes(x = cit, y = genome_size, fill = cit)) +
  ggtitle('Boxplot of genome size by citrate mutant type') +
  xlab('citrate mutant') +
  ylab('genome size') +
  theme(panel.grid.major = element_line(size = .5, color = "grey"),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))
ggsave("figures/genome-size-by-citrate-mutant.png")


#-------------------------------------------------------------------------------

# Another dataset
download.file("https://raw.githubusercontent.com/ILRI-Genomics-Platform/klebgen-africacdc-ilri/refs/heads/main/assets/genotype_vs_phenotype.tsv",
              "data/genotype_vs_phenotype-test.tsv")


genotype_vs_phenotype <- read.csv("data/genotype_vs_phenotype.tsv",
                                  sep = "\t")

# Kleborate
genotype_vs_phenotype |> 
  filter(tool == "Kleborate") |> 
  select(tool, drug, result) |> 
  # Count and compute proportions per drug
  count(drug, result) |> 
  group_by(drug) |> 
  mutate(pct = n / sum(n),
         total = sum(n)) |> 
  ungroup() |> 
  mutate(result = factor(result,
                         levels = c("agree",
                                    "major error (false resistance)",
                                    "very major error (missed resistance)"))) |> 
  ggplot(aes(x = drug, y = pct, fill = result)) +
  geom_col(width = 0.7, colour = "white") +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5),
            colour = "white", size = 4, fontface = "bold") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("agree" = "#2E8B57",
                               "major error (false resistance)" = "#E69F00",
                               "very major error (missed resistance)" = "#C0392B"),
                    name = "Result") +
  labs(title = "Kleborate genotype–phenotype concordance by drug",
       subtitle = "Numbers inside bars show isolate counts",
       x = NULL, y = "Percentage of isolates") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical",
        axis.text.x = element_text(face = "bold"))



# AMRFinderPlus
genotype_vs_phenotype |> 
  filter(tool == "AMRFinderPlus") |> 
  select(tool, drug, result) |> 
  # Count and compute proportions per drug
  count(drug, result) |> 
  group_by(drug) |> 
  mutate(pct = n / sum(n),
         total = sum(n)) |> 
  ungroup() |> 
  mutate(result = factor(result,
                         levels = c("agree",
                                    "major error (false resistance)",
                                    "very major error (missed resistance)"))) |> 
  ggplot(aes(x = drug, y = pct, fill = result)) +
  geom_col(width = 0.7, colour = "white") +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5),
            colour = "white", size = 4, fontface = "bold") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("agree" = "#2E8B57",
                               "major error (false resistance)" = "#E69F00",
                               "very major error (missed resistance)" = "#C0392B"),
                    name = "Result") +
  labs(title = "AMRFinderPlus genotype–phenotype concordance by drug",
       subtitle = "Numbers inside bars show isolate counts",
       x = NULL, y = "Percentage of isolates") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical",
        axis.text.x = element_text(face = "bold"))



# Both tools side by side
genotype_vs_phenotype |> 
  select(tool, drug, result) |> 
  count(tool, drug, result) |> 
  group_by(tool, drug) |> 
  mutate(pct = n / sum(n)) |> 
  ungroup() |> 
  mutate(result = factor(result,
                         levels = c("agree",
                                    "major error (false resistance)",
                                    "very major error (missed resistance)"))) |> 
ggplot(aes(x = drug, y = pct, fill = result)) +
  geom_col(width = 0.7, colour = "white") +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5),
            colour = "white", size = 4, fontface = "bold") +
  facet_wrap(~ tool) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("agree" = "#2E8B57",
                               "major error (false resistance)" = "#E69F00",
                               "very major error (missed resistance)" = "#C0392B"),
                    name = "Result") +
  labs(title = "Genotype–phenotype concordance by drug: AMRFinderPlus vs Kleborate",
       subtitle = "Numbers inside bars show isolate counts",
       x = NULL, y = "Percentage of isolates") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold", size = 14),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))

#-------------------------------------------------------------------------------

# A different way to visualize both Kleborate and AMRFInderPlus
genotype_vs_phenotype |> 
  select(tool, drug, result) |>
  count(tool, drug, result) |> 
  group_by(tool, drug) |> 
  mutate(pct = n / sum(n)) |> 
  ungroup() |> 
  mutate(result = factor(result,
                         levels = c("agree",
                                    "major error (false resistance)",
                                    "very major error (missed resistance)"))) |> 

ggplot(aes(x = tool, y = pct, fill = result)) +
  geom_col(width = 0.85, colour = "white") +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5),
            colour = "white", size = 3.8, fontface = "bold") +
  facet_grid(~ drug, switch = "x") +          # drug labels below, bars grouped
  scale_y_continuous(labels = scales::percent_format(),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("agree" = "#2E8B57",
                               "major error (false resistance)" = "#E69F00",
                               "very major error (missed resistance)" = "#C0392B"),
                    name = "Result") +
  labs(title = "Genotype–phenotype concordance: AMRFinderPlus vs Kleborate",
       subtitle = "Numbers inside bars show isolate counts",
       x = NULL, y = "Percentage of isolates") +
  theme_minimal(base_size = 13) +
  theme(strip.placement = "outside",
        strip.text.x = element_text(face = "bold", size = 12),
        panel.spacing.x = unit(0.8, "lines"),   # gap between drug groups
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "bottom")
ggsave("figures/concordance.png", width = 10, height = 5, dpi = 300)



