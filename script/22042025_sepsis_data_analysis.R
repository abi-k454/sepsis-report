#______________________----
#PACKAGES----
library(dplyr)
library(tidyverse)
library(janitor)
library(plotly)
#_____________________----
#IMPORT DATA----
#importing data of known sepsis patient data with healthy controls
data<-read_csv("data/sepsis_data.csv")
#_____________________----
#TIDY DATA----
#looking at format of data
glimpse(data)

#making sure raw data column names is in the right format
data <-
  data %>%
  clean_names()

#check for na values
is.na(data)%>%
  sum()

# Transpose the wide data: genes become columns, samples are rows
data_transposed <- data %>%
  column_to_rownames("gene_id") %>%
  t() %>%
  as.data.frame()

glimpse(data_transposed)

# Add sample names as a column
data_transposed <- data_transposed %>%
  rownames_to_column("sample") %>%
  mutate(group = case_when(
    str_detect(sample, "ctr") ~ "control",
    str_detect(sample, "shock") ~ "shock",
    str_detect(sample, "sepsis") ~ "sepsis",
    TRUE ~ "unknown"
  ))

# Remove zero-variance genes (columns)
expr_data <- data_transposed %>%
  select(-sample, -group) %>%
  select(where(~ sd(.) > 0))  # keep only columns with non-zero standard deviation

glimpse(expr_data)
# Run PCA
pca <- prcomp(expr_data, scale. = TRUE)

# Add metadata back
pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  mutate(sample = data_transposed$sample,
         group = data_transposed$group)
#_____________----
#PLOT----

ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(title = "PCA: Group Segregation of Samples",
       x = "Principal Component 1",
       y = "Principal Component 2",
       color = "Group") +
  theme_minimal()

summary(pca)
#____________________----
#PCA does not convey full data so increasing to 3D
# Build a dataframe with the top 3 PCs + sample info
pca_df <- as.data.frame(pca$x[, 1:3]) %>%
  mutate(sample = data_transposed$sample,
         group = data_transposed$group)

# 3D scatter plot
plot_ly(pca_df, x = ~PC1, y = ~PC2, z = ~PC3,
        color = ~group, colors = c("control" = "blue", "sepsis" = "orange", "shock" = "red"),
        text = ~sample,
        type = "scatter3d", mode = "markers") %>%
  layout(title = "3D PCA Plot of Gene Expression")

summary(pca)
