###############################################################

# A global investigation reveals limited evidence of impacts of bat exploitation of anthropogenic structures

## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

# Analysis performed with R studio (v. 4.2.2)
#setwd

# clean the workspace -----------------------------------------------------
rm(list = ls())

# Loading R package -------------------------------------------------------
libraries <- c(
  "metafor", "ggplot2", "tidyverse", "ggpubr",
  "RColorBrewer", "scatterpie", "dplyr", "gridtext",
  "maps", "ggsci", "networkD3", "irr", "lpSolve", "wrapr",
  "viridis", "stringr"
)

# Install and load missing libraries
for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    install.packages(lib, dependencies = TRUE)
  }
  library(lib, character.only = TRUE)
}



## Data preparation: Load datasets and get overall counts ----------------------------------------------------

# Loading the Database
meta <- read.csv(file = "Meta_24.csv")
data <- read.csv(file = "descriptiveES.csv")
pop_bat <- read.csv(file = "population.csv") #obtained online, see manuscript
pop_human <- read.csv("https://raw.githubusercontent.com/autistic96/project-2/main/world_population.csv")
#pop_human <- read.csv("world_population") if the online link ever fails to work, use the data file
world <- map_data("world") %>% subset(region != "Antarctica")
###############################################################

# Replace 'Habitat loss' with 'Roost loss' in the Impact and re-assigned columns
data <- data %>%
  mutate(
    Impact = str_replace_all(Impact, "Habitat loss", "Roost loss"),
    Impact.re.assigned = str_replace_all(Impact.re.assigned, "Habitat loss", "Roost loss")
  )


## Data preparation: Meta-analysis and plot (Fig 5): ----------------------------------------------------

# Meta-analysis
meta$Pearson.s_r <- as.numeric(meta$Pearson.s_r)
meta$N <- as.numeric(meta$N)
meta <- meta %>% select(ID, N, Genus, Response_Group, r = Pearson.s_r)


# Derive Fischer's Z
meta <- metafor::escalc(measure = "COR", ri = r, ni = N, data = meta)
meta <- na.omit(meta) # a couple did not calculate vi for models (Environmental Temp, Occupancy)

# meta-analysis for effects of Tag on "health" "behavior"
Bmeta <- rma.mv(yi, vi, random =  ~ 1 | ID, data = meta[meta$Response_Group == "Behavior",] )
Etmeta <- rma.mv(yi, vi, random =  ~ 1 | ID, data = meta[meta$Response_Group == "Environmental temperature",] ) 
Hmeta <- rma.mv(yi, vi, random =  ~ 1 | ID, data = meta[meta$Response_Group == "Health",] ) #might be unstable
Ometa <- rma.mv(yi, vi, random =  ~ 1 | ID, data = meta[meta$Response_Group == "Occupancy",] )

# Evaluation of publication bias via Rosenthal’s method.
failsafe_B <- fsn(yi, vi, type = "Rosenthal", data = meta[meta$Response_Group == "Behavior",])
failsafe_Et <- fsn(yi, vi, type = "Rosenthal", data = meta[meta$Response_Group == "Environmental temperature",]) 
failsafe_H <- fsn(yi, vi, type = "Rosenthal", data = meta[meta$Response_Group == "Health",]) 
failsafe_O <- fsn(yi, vi, type = "Rosenthal", data = meta[meta$Response_Group == "Occupancy",]) 

# Extracting estimates for plot # estimates, papers
plot <- data.frame(label = c( "Occupancy (7,3)", "Health (9,3)", "Environmental temperature (7,3)", "Behavior (6,2)"),
                   ES    = c( ((exp(Ometa$b)-1))/((exp(Ometa$b)+1)), ((exp(Hmeta$b)-1))/((exp(Hmeta$b)+1)),
                              ((exp(Etmeta$b)-1))/((exp(Etmeta$b)+1)), ((exp(Bmeta$b)-1))/((exp(Bmeta$b)+1))),
                   
                   L     = c( ((exp(Ometa$ci.lb)-1)/(exp(Ometa$ci.lb) +1)), ((exp(Hmeta$ci.lb)-1)/(exp(Hmeta$ci.lb)+1)), 
                              ((exp(Etmeta$ci.lb)-1)/(exp(Etmeta$ci.lb) +1)), ((exp(Bmeta$ci.lb)-1)/(exp(Bmeta$ci.lb)+1))), 
                   
                   U     = c( ((exp(Ometa$ci.ub)-1)/(exp(Ometa$ci.ub) +1)), ((exp(Hmeta$ci.ub)-1)/(exp(Hmeta$ci.ub)+1)), 
                              ((exp(Etmeta$ci.ub)-1)/(exp(Etmeta$ci.ub) +1)), ((exp(Bmeta$ci.ub)-1)/(exp(Bmeta$ci.ub)+1))) )


# Create a new variable in meta with the specified labels
meta$Specified_Label <- plot$label[match(meta$Response_Group, c("Occupancy", "Health", "Environmental temperature", "Behavior"))]


# Extracting estimates for BEHAVIOR plot
plotBehavior <- data.frame(label = c("C Use of torpor", "B Use of torpor", "A Use of torpor", "Use of torpor", "Frequency of deep torpor", "Torpid body temperature", "Behavior"),
                            
                            ES    = c( -0.022, -0.002, -0.285, -0.653, -0.161, -0.586, ((exp(Bmeta$b)-1))/((exp(Bmeta$b)+1)) ),
                            L     = c( NA, NA, NA, NA, NA, NA,  ((exp(Bmeta$ci.lb)-1)/(exp(Bmeta$ci.lb)+1)) ),
                            U     = c( NA, NA, NA, NA, NA, NA, ((exp(Bmeta$ci.ub)-1)/(exp(Bmeta$ci.ub)+1)) )
                            
)


plotBehavior$label <- factor(as.character(plotBehavior$label), levels = unique(plotBehavior$label))


B <- ggplot(plotBehavior, aes(y=label, x=ES, xmin=L, xmax=U)) + 
  geom_point(color = '#330066', shape=16, size=2.0) + 
  geom_errorbarh(height=0, size= 0.5, color = '#330066') +
  geom_vline(xintercept=0, color='black', linetype='dashed') +  
  scale_x_continuous(limits=c(-1,1), name=("")) +
  ylab('') +
  labs(subtitle="Behavior") +
  theme(text = element_text(size = 12),
        axis.text.y = element_text(face = c(rep('plain', 10), 'bold')),
        plot.margin = unit(c(0.5,0.4,0,2.0), 'lines'), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())




# Extracting estimates for ENVIRONMENTAL TEMP plot
plotEnvironment <- data.frame(label = c("A Air temperature", "Air temperature", "A Surface temperature", "Surface temperature", "Microclimate availability", "Day temperature", "Night temperature",  "Environmental temperature"),
                           
                           ES    = c( 0.726, 0.714, 0.818, 0.643, 0.818, 0.278, 0.341, ((exp(Etmeta$b)-1))/((exp(Etmeta$b)+1)) ),
                           L     = c( NA, NA, NA, NA, NA, NA, NA,  ((exp(Etmeta$ci.lb)-1)/(exp(Etmeta$ci.lb)+1)) ),
                           U     = c( NA, NA, NA, NA, NA, NA, NA, ((exp(Etmeta$ci.ub)-1)/(exp(Etmeta$ci.ub)+1)) )
                           
)


plotEnvironment$label <- factor(as.character(plotEnvironment$label), levels = unique(plotEnvironment$label))


Et <- ggplot(plotEnvironment, aes(y=label, x=ES, xmin=L, xmax=U)) + 
  geom_point(color = '#330066', shape=16, size=2.0) + 
  geom_errorbarh(height=0, size= 0.5, color = '#330066') +
  geom_vline(xintercept=0, color='black', linetype='dashed') + 
  scale_x_continuous(limits=c(-1,1), name=("")) +
  ylab('') +
  labs(subtitle="Environmental Temperature") +
  theme(text = element_text(size = 12),
        axis.text.y = element_text(face = c(rep('plain', 10), 'bold')),
        plot.margin = unit(c(0.5,0.2,0,1.5), 'lines'), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())


# Extracting estimates for OCCUPANCY plot
plotHealth <- data.frame(label = c("Juvenile parasitic load", "Female parasitic load", "A Female immune response", "Female immune response", "Juvenile body condition","Female body condition","A Spring female mass", "Spring female mass", "Antibody level", "Health"),
                            
                            ES    = c( 0.532, 0.211, -0.665, -0.333, -0.556, -0.450, 0.983, 0.999, 0.337, ((exp(Hmeta$b)-1))/((exp(Hmeta$b)+1)) ),
                            L     = c( NA, NA, NA, NA, NA, NA, NA, NA, NA, ((exp(Hmeta$ci.lb)-1)/(exp(Hmeta$ci.lb)+1)) ),
                            U     = c( NA, NA, NA, NA, NA, NA, NA, NA, NA, ((exp(Hmeta$ci.ub)-1)/(exp(Hmeta$ci.ub)+1)) )
                            
)


plotHealth$label <- factor(as.character(plotHealth$label), levels = unique(plotHealth$label))


H <- ggplot(plotHealth, aes(y=label, x=ES, xmin=L, xmax=U)) + 
  geom_point(color = '#330066', shape=16, size=2.0) + 
  geom_errorbarh(height=0, size= 0.5, color = '#330066') +
  geom_vline(xintercept=0, color='black', linetype='dashed') +
  scale_x_continuous(limits=c(-1,1), name=("")) +
  ylab('') +
  labs(subtitle="Health") +
  theme(text = element_text(size = 12),
        axis.text.y = element_text(face = c(rep('plain', 10), 'bold')),
        plot.margin = unit(c(0.5,0.4,0,1.2), 'lines'), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())




# Extracting estimates for OCCUPANCY plot
plotOccupancy <- data.frame(label = c("D Roost selection", "C Roost selection", "B Roost selection", "Roost selection", "a Number of reproductive females", "Number of reproductive females", "Consecutive days in roost",  "Occupancy"),
                         
                         ES    = c( -0.999, 0.983, 0.998, 0.997, 0.460, 0.758, 0.522, ((exp(Ometa$b)-1))/((exp(Ometa$b)+1)) ),
                         L     = c( NA, NA, NA, NA, NA, NA, NA,  ((exp(Ometa$ci.lb)-1)/(exp(Ometa$ci.lb)+1)) ),
                         U     = c( NA, NA, NA, NA, NA, NA, NA, ((exp(Ometa$ci.ub)-1)/(exp(Ometa$ci.ub)+1)) )
                         
)


plotOccupancy$label <- factor(as.character(plotOccupancy$label), levels = unique(plotOccupancy$label))


O <- ggplot(plotOccupancy, aes(y=label, x=ES, xmin=L, xmax=U)) + 
  geom_point(color = '#330066', shape=16, size=2.0) + 
  geom_errorbarh(height=0, size= 0.5, color = '#330066') +
  geom_vline(xintercept=0, color='black', linetype='dashed') + 
  scale_x_continuous(limits=c(-1,1), name=("")) +
  ylab('') +
  labs(subtitle="Occupancy") +
  theme(text = element_text(size = 12),
        axis.text.y = element_text(face = c(rep('plain', 10), 'bold')),
        plot.margin = unit(c(0.5,0.2,0,-0.5), 'lines'), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())


# View and save
# Used Inkscape to made visual edits.

library(gridExtra)
Meta <- grid.arrange(B, Et, H, O, ncol = 2)
#ggsave("Fig5_R1.svg", Meta, dpi = 600, width = 225, height = 100, units = "mm")



## qnorm (calculating Z value) for Wilcoxon Meta analyses
# these data were entered into the meta-analysis csv file

#WoS_667 for both
#qnorm(0.01/2) #-2.575829
###############################################################


## Data preparation (Fig 1a): Bar plot of different structures used by bats----------------------------------------------------

# Loading the Database
barbat <- data %>% select(ID, Higher_Geography, System_renamed)

# removing missing data
barbat <- na.omit(barbat)

# separating higher geography if multiple listed within a cell and replicating line; 
# separating systems if multiple listed within a cell and replicating line;
# group by ID, when higher geography and system are distinct
barbat <- barbat %>%
  separate_rows(Higher_Geography, sep = " ;", convert = TRUE) %>%
  mutate(Higher_Geography = trimws(Higher_Geography)) %>%
  separate_rows(System_renamed, sep = " ;", convert = TRUE) %>%
  mutate(System_renamed = trimws(System_renamed))

# group by ID, when higher geography and system are distinct
# this means that some papers will have multiple lines IF there are multiple regions and multiple systems mentioned
barbat <- barbat %>%
  group_by(ID) %>%
  distinct(Higher_Geography, System_renamed, .keep_all = TRUE)

# Get the table of counts for each factor in System_renamed
system_table <- table(barbat$System_renamed)

# Calculate percentages
percentages <- prop.table(system_table) * 100

# Create a data frame for better visualization
result_table <- data.frame(System = names(system_table), Count = as.vector(system_table), Percentage = as.vector(percentages))

# Print the result
print(result_table)

# data for the plot
summary_data <- barbat %>% group_by(Higher_Geography, System_renamed) %>% count()

# set color pallet
mixa <- c("#FF7F11", "#f0c571", "#8F4300", "#FFBE85", "#745D58", "#DFA616", "#3D1D00", "#FF9D47")

# plot
x <- ggplot(summary_data, aes(x = Higher_Geography, y = n, fill = System_renamed)) +
  geom_col(position = "dodge", width = 0.9) +  
  labs(x = "", y = "Quantity") +
  scale_fill_manual(values = mixa, name = "System") + 
  theme_minimal(base_size = 24) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"), 
        axis.ticks = element_line(color = "black"),  
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24))


#ggsave("Fig1a.png", x, dpi = 600, width = 350, height = 150, bg = "white", units = "mm")
###############################################################


## Data preparation (Fig 1b) Trend of Roost Type by Biogeographic Realm----------------------------------------------------

roostbio <- data %>% dplyr::select(ID, Year_publication, Higher_Geography, System_renamed)

# separating higher geography if multiple listed within a cell and replicating line; 
# separating systems if multiple listed within a cell and replicating line;
# group by ID, when higher geography and system are distinct
roostbio <- roostbio %>%
  separate_rows(Higher_Geography, sep = " ;", convert = TRUE) %>%
  mutate(Higher_Geography = trimws(Higher_Geography)) %>%
  separate_rows(System_renamed, sep = " ;", convert = TRUE) %>%
  mutate(System_renamed = trimws(System_renamed))

# group by ID, when higher geography and system are distinct
# this means that some papers will have multiple lines IF there are multiple regions and multiple systems mentioned
roostbio <- roostbio %>%
  group_by(ID) %>%
  distinct(Higher_Geography, System_renamed, .keep_all = TRUE)

# calculate the counts of System_renamed by Year_publication and Higher_Geography
roostbio_counts <- roostbio %>%
  group_by(Year_publication, Higher_Geography, System_renamed) %>%
  summarise(count = n(), .groups = 'drop')


# Determine the range of years in your dataset
min_year <- min(roostbio_counts$Year_publication)
max_year <- max(roostbio_counts$Year_publication)

# Create the dot plot with consistent x-axis limits and breaks
sup <- ggplot(roostbio_counts, aes(x = Year_publication, y = System_renamed, size = count, color = System_renamed)) +
  geom_point() +
  facet_wrap(~ Higher_Geography, scales = "free_y") +
  scale_color_manual(values = mixa) +  
  scale_x_continuous(
    limits = c(min_year, max_year),  
    breaks = seq(min_year, max_year, by = 5)  
  ) +
  scale_size_continuous(
    name = "Count",     
    breaks = c(1, 5, 10), 
    range = c(2, 10)    
  ) +
  labs(
    x = "Year",
    y = "System"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # Remove legends
  )

#ggsave("Fig1b.png", sup, dpi = 600, width = 350, height = 150, bg = "white", units = "mm")
###############################################################


## Data preparation (Fig 2a): World plot with pie charts for total impacts to bats (studied and documented/mentioned) ----------------------------------------------------

# Load data
batworld <- data %>% select(ID, Higher_Geography, Risk, Regress, Impact.re.assigned)

# Subset the data to remove Human 
batworld <- subset(batworld, Risk != "Human")

# Remove rows with NA values
batworld <- batworld[complete.cases(batworld), ]

# separating Impact.re.assigned & Higher_Geography if multiple listed within a cell and replicating line
batworld <- batworld %>%
  separate_rows(Impact.re.assigned, sep = " ;", convert = TRUE) %>%
  mutate(Impact.re.assigned = trimws(Impact.re.assigned)) %>%
  separate_rows(Higher_Geography, sep = " ;", convert = TRUE) %>%
  mutate(Higher_Geography = trimws(Higher_Geography))

# change No to Document/Mention
batworld <- batworld %>%
  mutate(Regress = if_else(Regress == "No", "Document/Mention", Regress))

# Group by Unique ID retaining unique Regress, and Impact; e.g., type of colony might create duplicates so this is necessary
batworld <- batworld %>%
  group_by(ID) %>%
  distinct(Higher_Geography, Regress, Impact.re.assigned, .keep_all = TRUE)

# Convert column to character type
batworld$Regress <- as.character(batworld$Regress)

# attributes levels
levels(batworld$Higher_Geography)

# remove columns not needed for the rest
batworld<- batworld[,-c(1,3,5)]

# Summary for the pie figure
summary_df <- batworld %>%
  group_by(Higher_Geography, Regress) %>%
  summarise(Frequency = n()) %>%
  ungroup()

# add coordinates for the realms
coordinates_df <- data.frame(
  Higher_Geography = c("Afrotropical", "Australasian", "Indomalayan", "Nearctic", "Neotropical", "Palearctic"),
  long = c(21, 133.4, 79, -104, -57, 45),
  lat = c(8, -23.4, 27, 44, -16.8, 56)
)

# Merge the coordinates with the summary data
merged_df <- merge(summary_df, coordinates_df, by = "Higher_Geography")

# Reshape the df
reshaped_df <- merged_df %>%
  pivot_wider(names_from = Regress, values_from = c("Frequency"), names_sep = "_") 

# Calculate the total number of impacts per realm
df_sum <- reshaped_df %>%
  group_by(Higher_Geography) %>%
  summarise(total_studies = sum(`Study`+ `Document/Mention`, na.rm = TRUE)) 

# Add the total number of impacts to the original data frame
df <- left_join(reshaped_df, df_sum, by = "Higher_Geography")

# Calculate proportional area for the pie chart
df$proportional_area <- sqrt(df$total_studies)

# Create world map with colors
worldbatmap <- ggplot(world, aes(long, lat)) +
  geom_map(map = world, aes(map_id = region), 
           color = "#ADADAD", fill = "#ADADAD", size = 0.3) +
  coord_quickmap() 

# Append the pie charts
worldbatmap <- worldbatmap + 
  geom_scatterpie(data = df, cols = c("Study", "Document/Mention"), color = NA, alpha= .9,
                  aes(x = long, y= lat, group = Higher_Geography, r = proportional_area)) 

# plot
my_colors <- c("#FB6376","#5D2A42")
worldbatmap <- worldbatmap + scale_fill_manual("",labels = c("Study", "Document/Mention"), values = my_colors)+
  labs(title = NULL, subtitle = NULL, x = NULL, y = "") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),  
        axis.ticks.y = element_blank()) +
  ggtitle("A")

# Specify min and max values for the legend based on total_studies
min_total_studies <- min(df$total_studies)
max_total_studies <- max(df$total_studies)

# Calculate proportional area for the legend based on the range of total_studies
legend_proportional_area <- sqrt(c(min_total_studies, max_total_studies))

# Add the scatterpie legend with fixed min and max values and proportional area size
worldbatmap <- worldbatmap + 
  geom_scatterpie_legend(legend_proportional_area, x = -150, y = -45, n = 2, 
                         labeller = function(x) c(min_total_studies, max_total_studies))


#ggsave("WorldbatmapPlot.svg", worldbatmap, dpi = 600, width = 350, height = 200, units = "mm")


# This makes the heat map for bat count. I used Inkscape to overlay pies.
# plot
abundance <- ggplot(pop_bat) +
  geom_map(
    data = world, map = world, aes(map_id = region),
    fill = "white", color = "#7f7f7f", size = 0.25
  ) +
  geom_map(map = world, aes(map_id = region, fill = Count), size = 0.25) +
  scale_fill_gradient(low = "#CCCCCC", high = "#474747", name = "Quantity") +
  expand_limits(x = world$long, y = world$lat) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank()    
  ) +
  ggtitle("A")

#ggsave("abundance1.svg", abundance, dpi = 600, width = 350, height = 150, units = "mm")


# Other summary info
# Create a table of counts for the combinations of Higher_Geography and Regress
table_counts <- table(batworld$Higher_Geography, batworld$Regress)

# Print the count table
print(table_counts)

# Calculate the percentage for each factor in Higher_Geography
percentage_table <- prop.table(table_counts, margin = 1) * 100

# Print the percentage table
print(percentage_table)
###############################################################

## Data preparation (Fig 2b): Trend Figure for just impacts associated with bats; not humans ----------------------------------------------------

# selecting data for the regression of proportions
trend <- data %>% dplyr::select(ID, Year_publication, Risk, Regress, Impact.re.assigned)

# Subset the data to remove Human and remove NAs
trend <- subset(trend, Risk != "Human" & !is.na(Risk))

# Separate Impact.re.assigned if multiple values listed within a cell and replicate lines
# Also, remove white space around words
trend <- trend %>%
  tidyr::separate_rows(Impact.re.assigned, sep = " ;", convert = TRUE) %>%
  mutate(Impact.re.assigned = trimws(Impact.re.assigned))

# Retain unique IDs based on unique Regress (study, docu/mention) and Impact
trend <- trend %>% 
  group_by(ID) %>%
  distinct(Regress, Impact.re.assigned, .keep_all = TRUE)

# Create a table with counts
trend_table <- data.frame(table(trend$Year_publication, trend$Regress))

# create column names
colnames(trend_table) <- c("Year", "Impact", "N")

# Using dplyr to create trend_data
trend_data <- trend_table %>%
  group_by(Year, Impact) %>%
  summarize(N = sum(N)) %>%
  spread(Impact, N, fill = 0) %>%
  ungroup()

# Convert yr to numeric
trend_data$Year <- as.numeric(as.character(trend_data$Year))

# Reshape data frame to long format for ggplot2
trend_long <- tidyr::gather(trend_data, key = "status", value = "value", -Year)

# Convert Year to factor with desired order
trend_long$Year <- factor(trend_long$Year, levels = unique(trend_long$Year))

# Set custom colors for fill
my_colors <- c("#FB6376","#5D2A42")

# Create stacked bar plot with custom colors and renamed labels, including missing years
t <- ggplot(trend_long, aes(x = Year, y = value, fill = factor(status, levels = c("Study", "No")))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors, labels = c("Study", "Document/Mention")) + 
  labs(x = NULL, y = "Quantity", fill = "") +  
  theme_minimal(base_size = 12) +
  theme(plot.margin = unit(c(1, 1, 0, 1), "lines")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  
        panel.grid = element_blank(),  
        panel.border = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.ticks = element_line(color = "black"), 
        axis.ticks.length = unit(0.2, "cm"),  
        axis.ticks.x = element_line(),  
        axis.ticks.y = element_line(),  
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +  
  scale_x_discrete(breaks = unique(trend_long$Year)[c(TRUE, rep(FALSE, 9))])  


#ggsave("TrendPlot.png", t, dpi = 600, width = 250, height = 200, units = "mm")
###############################################################


## Data preparation: Regression on proportion of impacts studied for bats between 1997-2022 ----------------------------------------------------

# select years for the regression model
regression <- trend[trend$Year_publication > 1996 & trend$Year_publication < 2023, ]

# data for glm
regress2 <- data.frame(table(regression$Year_publication, regression$Regress)) ; colnames(regress2) <- c("yr","Study","N")

# ensure year is numeric
regress2$yr <- as.numeric(as.character(regress2$yr))

# create data frame for glm data
glm_data <- data.frame(yr = unique(regress2$yr),
                       study = regress2[regress2$Study=="Study",]$N, 
                       mention = regress2[regress2$Study=="No",]$N)


# Fit the model
m1 <- glm(cbind(study, mention) ~ yr, data = glm_data, family = "binomial")

# check model performance
performance::check_overdispersion(m1) 
rsq::rsq(m1)

# calculate the proportion
glm_data$sum <- glm_data$study + glm_data$mention
glm_data$propStudy <- glm_data$study/glm_data$sum

# make numeric
glm_data$propStudy <- as.numeric(glm_data$propStudy)

# obtain model parameters (log odds [coefficient], SE, p)
(pM1 <- parameters::model_parameters(m1))
###############################################################


## Data preparation (Fig 2c): Bar plots of Studied consequences + Mentioned/Documented by Bat family----------------------------------------------------

# Load data
family <- data %>% select(ID, Family, Risk, Regress, Impact.re.assigned)

# change No to Document/Mention
family <- family %>%
  mutate(Regress = if_else(Regress == "No", "Document/Mention", Regress))

# separating Family and Impact.re.assigned if multiple listed within a cell and replicating line
family <-tidyr::separate_rows(family, Family, sep = " ;", convert = T)
family <-tidyr::separate_rows(family, Impact.re.assigned, sep = " ;", convert = T)

# remove white space around words
family$Impact.re.assigned <- trimws(family$Impact.re.assigned)
family$Family <- trimws(family$Family)

# remove all NAs to retain just Study and Document/Mention
family <- na.omit(family)

# Subset the data for Bat
family_bat <- subset(family, Risk == "Bat")

# Group by Unique ID retaining unique Family, Regress, Impact.re.assigned
family_bat <- family_bat %>%
  group_by(ID) %>%
  distinct(Family, Regress, Impact.re.assigned, .keep_all = TRUE)

# summarize data - used for percentages for where studies occurred
sum(table(family_bat$Family)) #total
table(family_bat$Family)/sum(table(family_bat$Family))*100

# Summarize data - used for percentages for Family
family_summary <- table(family_bat$Family)
family_total <- sum(family_summary)
family_percentage <- family_summary / family_total * 100

# Reorder the levels of the "Regress" variable
family_bat$Regress <- factor(family_bat$Regress, levels = c("Study", "Document/Mention"))

# Set color values
bat_colors <- c("#FB6376", "#5D2A42")

# Plot for Bat
family <- ggplot(family_bat, aes(x = Family, fill = Regress)) +
  geom_bar(position = "stack", stat = "count", width = 0.8) +
  scale_fill_manual(values = bat_colors) +
  theme_minimal(base_size = 24) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.ticks = element_line(color = "black"),  
        axis.title = element_blank(),
        plot.margin = margin(l = 60),
        legend.text = element_blank()) +
  guides(fill = guide_legend(title = NULL)) +
  guides(fill = FALSE)

print(family)


#ggsave("familyplot.svg", family, dpi = 600, width = 350, height = 150, units = "mm")
###############################################################


## Data preparation (Fig 2d & Fig 4b): Bar plots of Studied impacts + Mentioned/Documented for Bats + Humans----------------------------------------------------

# Load data
bat_data <- data %>% select(ID, Risk, Regress, Impact.re.assigned)

# Subset the data for Bat with associated Risk
bat_data <- subset(bat_data, Risk == "Bat")

# Count unique IDs for studies on bat risks
unique_id_count_bat <- bat_data %>%
  distinct(ID) %>%
  count()

# Print the count
print(unique_id_count_bat)

# change No to Document/Mention
bat_data <- bat_data %>%
  mutate(Regress = if_else(Regress == "No", "Document/Mention", Regress))

# Separate Impact.re.assigned if multiple listed within a cell and replicate lines
bat_data <- bat_data %>%
  separate_rows(Impact.re.assigned, sep = " ;", convert = TRUE) %>%
  mutate(
    Impact.re.assigned = trimws(Impact.re.assigned)
  )

#Group by Unique ID retaining unique Risk, Regress, and Impact
bat_data <- bat_data %>%
  group_by(ID) %>%
  distinct(Risk, Regress, Impact.re.assigned, .keep_all = TRUE)


# Summarize data - used for percentages for Impact.re.assigned
impact_summary <- table(bat_data$Impact.re.assigned)
impact_total <- sum(impact_summary)
impact_percentage <- impact_summary / impact_total * 100

# Summarize data - used for percentages for Regress
regress_summary <- table(bat_data$Regress)
regress_total <- sum(regress_summary)
regress_percentage <- regress_summary / regress_total * 100


# Reorder the levels of the "Regress" variable
bat_data$Regress <- factor(bat_data$Regress, levels = c("Study", "Document/Mention"))

# Set common breaks for both Y-axes, stopping at 40
common_y_breaks <- seq(0, 40, by = 10)

# Set color values
bat_colors <- c("#FB6376", "#5D2A42")

# Plot for Bat
bat_plot <- ggplot(bat_data, aes(x = Impact.re.assigned, fill = Regress)) +
  geom_bar(position = "stack", stat = "count", width = 0.8, show.legend = FALSE) +
  labs(x = "",
       y = "Count") +
  scale_fill_manual(values = bat_colors) +
  scale_y_continuous(breaks = common_y_breaks) +  
  theme_minimal(base_size = 24) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.ticks = element_line(color = "black"),  
        plot.margin = margin(l = 80),
        legend.text = element_blank()) +
  guides(fill = guide_legend(title = NULL))

print(bat_plot)
#ggsave("batplot.svg", bat_plot, dpi = 600, width = 300, height = 150, units = "mm")


# Now, repeat the same process for Human
human_data <- data %>% select(ID, Risk, Impact.re.assigned, Impact_studied_human)

# Subset the data for Human
human_data <- subset(human_data, Risk == "Human")

# Count unique IDs
unique_id_count <- human_data %>%
  distinct(ID) %>%
  count()

# Print the count
print(unique_id_count)

# Separate Impact.re.assigned if multiple listed within a cell and replicate lines
human_data <- human_data %>%
  separate_rows(Impact.re.assigned, sep = " ;", convert = TRUE) %>%
  mutate(
    Impact.re.assigned = trimws(Impact.re.assigned)
  )

# Group by Unique ID retaining unique Risk, Impact_studied_human, and Impact
human_data <- human_data %>%
  group_by(ID) %>%
  distinct(Risk, Impact_studied_human, Impact.re.assigned, .keep_all = TRUE)

# Summarize data - used for percentages for Impact.re.assigned
impact_summary <- table(human_data$Impact.re.assigned)
impact_total <- sum(impact_summary)
impact_percentage <- impact_summary / impact_total * 100

# Define your own unique colors
my_colors <- c("#5CBEFF", "#2C7CDD", "#14427B", "#B4E1FF")

# Set common breaks for both Y-axes, stopping at 40
common_y_breaks <- seq(0, 40, by = 10)

# Create the plot
human_plot <- ggplot(human_data, aes(x = Impact.re.assigned, fill = as.factor(Impact.re.assigned))) +
  geom_bar(position = "stack", stat = "count", width = 0.4) +
  labs(title = "",
       x = "",
       y = "Quantity") +
  scale_y_continuous(breaks = common_y_breaks) + 
  theme_minimal(base_size = 24) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.ticks = element_line(color = "black"),  
        plot.margin = margin(b = 61)) +
  scale_fill_manual(values = my_colors)  

# Display the plot
human_plot

#ggsave("humanplot.svg", human_plot, dpi = 600, width = 200, height = 150, units = "mm")


###############################################################


## Data preparation: Fisher's test of Regress/System due to small samples----------------------------------------------------
#testing if there is a relationship between Study/No and the System type

# Load data
fisher_system <- data %>% select(ID, Risk, Regress, System_renamed)

# Subset the data for Bat with associated Risk
fisher_system <- subset(fisher_system, Risk == "Bat")

# Separate Impact.re.assigned if multiple listed within a cell and replicate lines
fisher_system <- fisher_system %>%
  separate_rows(System_renamed, sep = " ;", convert = TRUE) %>%
  mutate(System_renamed = trimws(System_renamed)    
  )

# group by ID and retain unique regress and system
fisher_system <- fisher_system %>% 
  group_by(ID) %>%
  distinct(Regress, System_renamed, .keep_all = TRUE)

# Create the contingency table
contingency_table <- table(fisher_system$System_renamed, fisher_system$Regress)

# Perform Fisher's exact test
fisher_test <- fisher.test(contingency_table)

# Print the results
print(fisher_test)
###############################################################


## Data preparation (Fig 3): Heat plot of different structures used by bats and risks----------------------------------------------------
# these numbers won't match the trend as we are looking at specific structures and risks: that is ok that they don't match

# Loading the Database 
siterisk <- data %>% select(ID, Higher_Geography, System_renamed, Regress, Impact.re.assigned, Risk)

# Filter rows where the value in the 'Risk' column is "Bat" and remove NA
siterisk_bat <- siterisk %>%
  filter(Risk == "Bat") %>%
  na.omit()

# Separate and trim Impact.re.assigned
siterisk_bata <- siterisk_bat %>%
  separate_rows(Impact.re.assigned, sep = " ;", convert = TRUE) %>%
  mutate(Impact.re.assigned = trimws(Impact.re.assigned))

# Separate and trim System_renamed
siterisk_bata <- siterisk_bata %>%
  separate_rows(System_renamed, sep = " ;", convert = TRUE) %>%
  mutate(System_renamed = trimws(System_renamed))

# Separate and trim Higher_Geography
siterisk_bata <- siterisk_bata %>%
  separate_rows(Higher_Geography, sep = " ;", convert = TRUE) %>%
  mutate(Higher_Geography = trimws(Higher_Geography))

# Rename "No" to "Document/Mention" in the Regress column
siterisk_bata <- siterisk_bata %>%
  mutate(Regress = case_when(
    Regress == "No" ~ "Document/Mention",
    TRUE ~ Regress  # Keep the original value if it's not "No"
  ))

# Aggregating the data to get counts of studied and mentioned impacts by Higher_Geography, System_renamed, Impact.re.assigned, and Impact_studied
heatmap_data <- siterisk_bata %>%
  group_by(Regress, Impact.re.assigned, System_renamed, Higher_Geography) %>%
  summarize(Count = n()) %>%
  ungroup()

# Calculate counts and percentages for summary data
result <- heatmap_data %>%
  group_by(System_renamed, Higher_Geography, Impact.re.assigned, Regress) %>%
  summarise(
    Count = sum(Count)
  ) %>%
  group_by(System_renamed, Higher_Geography, Impact.re.assigned) %>%
  mutate(
    Total_Count = sum(Count[Regress %in% c("Study", "Document/Mention")]),
    Percent = ifelse(Regress %in% c("Study", "Document/Mention"), Count / Total_Count * 100, NA)
  )  

# Display the result
print(result, n = 160)

# Calculate count of "Study" and "Document/Mention" within each group
study_mention_counts <- heatmap_data %>%
  group_by(System_renamed, Higher_Geography, Impact.re.assigned) %>%
  summarise(
    Total_Count = sum(Count),
    Study_Count = sum(Count[Regress == "Study"]),
    Document_Mention_Count = sum(Count[Regress == "Document/Mention"])
  ) %>%
  mutate(
    Study_Percentage = (Study_Count / Total_Count) * 100,
    Document_Mention_Percentage = (Document_Mention_Count / Total_Count) * 100
  )

# Merge counts back into the original dataframe
heatmap_data <- merge(heatmap_data, study_mention_counts, by = c("System_renamed", "Higher_Geography", "Impact.re.assigned"), all.x = TRUE)

# Create the gradient plot
### VIRIDIS COLOR SCALE AND DIFFERENT TEXT COLOR
heat <- ggplot(heatmap_data, aes(x = System_renamed, y = Impact.re.assigned, fill = Study_Percentage)) +
  geom_tile() +
  scale_fill_viridis(na.value = "white", 
                     guide = guide_colorbar(reverse = FALSE, title.position = "top", title.vjust = 1, label.size = 15),
                     name = "Percent (%) studied impacts") +
  labs(x = "System", y = "Impact") +
  theme_linedraw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  facet_wrap(~ Higher_Geography, nrow = 2)


heat <- heat +
  geom_text(aes(label = Total_Count, 
                color = ifelse(Study_Percentage > 50, "black", "white")), size = 3) +
  scale_color_identity()


# Print the heatmap
print(heat)

#ggsave("Fig3_R1.png", heat, dpi = 600, width = 350, height = 150, units = "mm")
###############################################################


## Data Preparation (Fig 4a): For human population map and pie charts of impacts ----------------------------------------------------
#the goal is to plot where studies on bats and humans are 
#https://rpubs.com/kelly_eng03/1094370
humanworld <- data %>% select(ID, Higher_Geography, Risk, Impact_studied_human, Impact.re.assigned)

#separating Impact.re.assigned & Higher_Geography if multiple listed within a cell and replicating line
humanworld <- humanworld %>%
  separate_rows(Impact.re.assigned, sep = " ;", convert = TRUE) %>%
  mutate(Impact.re.assigned = trimws(Impact.re.assigned)) %>%
  separate_rows(Higher_Geography, sep = " ;", convert = TRUE) %>%
  mutate(Higher_Geography = trimws(Higher_Geography))

# Subset the data to remove Human but retain the NA studies
humanworld <- subset(humanworld, Risk != "Bat")

# Remove rows with NA values
humanworld <- humanworld[complete.cases(humanworld), ]

# Group by Unique ID retaining unique Risk, Regress, and Impact
humanworld <- humanworld %>%
  group_by(ID) %>%
  distinct(Higher_Geography, Risk, Impact_studied_human, Impact.re.assigned, .keep_all = TRUE)

# Convert column to character type
humanworld$Impact_studied_human <- as.character(humanworld$Impact_studied_human)

# Count and percentage of each category within Higher_Geography for summary data
#Impacts to humans from bats occupying anthropogenic structures
higher_geography_summary <- humanworld %>% 
  group_by(Higher_Geography) %>% 
  summarise(count = n()) %>% 
  mutate(percentage = (count / sum(count)) * 100)

# view
print(higher_geography_summary)

# create levels
levels(humanworld$Higher_Geography)

#remove columns not needed for the rest
humanworldpie <- humanworld[,-c(1,3,4)]

# summarize data 
sum(table(humanworldpie$Higher_Geography)) #total
table(humanworldpie$Higher_Geography)/sum(table(humanworldpie$Higher_Geography))

# Summary for the pie figure
summary_df <- humanworldpie %>%
  group_by(Higher_Geography, Impact.re.assigned) %>%
  summarise(Frequency = n()) %>%
  ungroup()

# add coordinates for the realms
coordinates_df <- data.frame(
  Higher_Geography = c("Afrotropical", "Australasian", "Indomalayan", "Nearctic", "Neotropical", "Palearctic"),
  long = c(21, 133.4, 79, -104, -57, 45),
  lat = c(8, -23.4, 27, 44, -16.8, 56)
)

# Merge the coordinates with the summary data
merged_df <- merge(summary_df, coordinates_df, by = "Higher_Geography")

# Reshape the df
reshaped_df <- merged_df %>%
  pivot_wider(names_from = Impact.re.assigned, values_from = "Frequency", names_sep = "_") 

# Calculate the total number of studies per region
df_sum <- reshaped_df %>%
  group_by(Higher_Geography) %>%
  summarise(total_studies = sum(`Nuisance`, `Parasitism`, `Pathogen`, `Other`, na.rm = TRUE))

# Add the total number of studies to the original data frame
df <- left_join(reshaped_df, df_sum, by = "Higher_Geography")

# Calculate proportional area for the pie chart
df$proportional_area <- sqrt(df$total_studies)

# Create world map with colors
worldhumanmap <- ggplot(world, aes(long, lat)) +
  geom_map(map = world, aes(map_id = region), 
           color = "#ADADAD", fill = "#ADADAD", size = 0.3) +
  coord_quickmap() 

# Recode NA values to 0 in relevant columns
df$Nuisance[is.na(df$Nuisance)] <- 0
df$Parasitism[is.na(df$Parasitism)] <- 0
df$Pathogen[is.na(df$Pathogen)] <- 0
df$Other[is.na(df$Other)] <- 0

# Specify min and max values for the legend based on total_studies
min_total_studies <- min(df$total_studies)
max_total_studies <- max(df$total_studies)

# Calculate proportional area for the legend based on the range of total_studies
legend_proportional_area <- sqrt(c(min_total_studies, max_total_studies))

# Define size factor
size_factor <- 2.0

# Define your own unique colors
my_colors <- c("#5CBEFF", "#2C7CDD", "#14427B", "#B4E1FF")

# Append the pie charts
worldhumanmap <- worldhumanmap + 
  geom_scatterpie(data = df, cols = c("Nuisance", "Other", "Parasitism", "Pathogen"), color = NA, alpha = 0.9,
                  aes(group = Higher_Geography, x = long, y = lat, r = proportional_area * size_factor, fill = ..group..)) +
  scale_fill_manual(values = my_colors) +
  labs(title = NULL, subtitle = NULL, x = NULL, y = "") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),  
        axis.ticks.y = element_blank()) +
  ggtitle("A")

# Display the map with pie charts and legend
worldhumanmap <- worldhumanmap + 
  geom_scatterpie_legend(legend_proportional_area, x = -150, y = -45, n = 2, 
                         labeller = function(x) c(min_total_studies, max_total_studies))

# view
worldhumanmap

#legend size changed with Inkscape to reflect the appropriate size.
#pies overlaid on human_pop map using Inkscape.

#ggsave("WorldhumanpiePlot.svg", worldhumanmap, dpi = 600, width = 350, height = 200, units = "mm")


# rename the columns for the human population map
pop_human <- pop_human %>% rename("region" = Country.Territory, "Count" = X2022.Population)

# select specific columns
pop_human <- pop_human %>% select(region, Count)

# Define replacements
replacements <- c("United States" = "USA",
                  "Antigua and Barbuda" = "Antigua",
                  "British Virgin Islands" = "Virgin Islands",
                  "DR Congo" = "Democratic Republic of the Congo",
                  "Eswatini" = "Swaziland",
                  "Vatican City" = "Vatican",
                  "Republic of the Congo" = "Republic of Congo",
                  "Saint Kitts and Nevis" = "Saint Kitts",
                  "Saint Vincent and the Grenadines" = "Saint Vincent",
                  "United Kingdom" = "UK",
                  "Trinidad and Tobago" = "Trinidad")

# Apply replacements using case_when
pop_human <- pop_human %>%
  mutate(region = case_when(
    region %in% names(replacements) ~ replacements[region],
    TRUE ~ region
  ))

# Create new rows
new_rows <- list(list(region = "Nevis", Count = 47657),
                 list(region = "Tobago", Count = 1531044))

# Add new rows
pop_human <- bind_rows(pop_human, new_rows)

#hong kong and macau - added pop value to china as they were listed as individual countries
#in the other dataset
pop_human <- pop_human %>%
  mutate(Count = ifelse(region == "China", 1434071370, Count))

# plot
humanabundance <- ggplot(pop_human) +
  geom_map(
    data = world, map = world, aes(map_id = region),
    fill = "white", color = "#7f7f7f", size = 0.25
  ) +
  geom_map(map = world, aes(map_id = region, fill = Count), size = 0.25) +
  scale_fill_gradient(
    low = "#CCCCCC", high = "#474747", name = "Quantity",
    labels = c("500M", "1B"),  # Labels for 5e+08 and 1e+09
    breaks = c(5e+08, 1e+09)  # Breaks for 5e+08 and 1e+09
  ) +
  expand_limits(x = world$long, y = world$lat) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank()
  )

humanabundance

#ggsave("humanabundancex.svg", humanabundance, dpi = 600, width = 350, height = 150, units = "mm")
###############################################################


## Data preparation: Other summary data----------------------------------------------------

## for obtaining the count of unique IDs in the full dataset
unique_id_count_data <- data %>%
  distinct(ID) %>%
  count()

# Print the count
print(unique_id_count_data)


## For determining the % of unique papers from each biogeographical realms
geog <- data %>% select(ID, Higher_Geography, Risk, Regress)

geog <- geog %>%
  separate_rows(Higher_Geography, sep = " ;", convert = TRUE) %>%
  mutate(Higher_Geography = trimws(Higher_Geography))

geog <- geog %>%
  group_by(ID) %>%
  distinct(Higher_Geography, .keep_all = TRUE)

# Count the occurrences of each factor in Higher_Geography
count_table <- table(geog$Higher_Geography)

# Calculate percentages
percentage_table <- prop.table(count_table) * 100

# Create a data frame for better visualization
result_table <- data.frame(
  Higher_Geography = names(count_table),
  Count = as.vector(count_table),
  Percentage = as.vector(percentage_table)
)

# Print the result
print(result_table)
###############################################################
