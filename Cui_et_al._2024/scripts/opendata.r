library(Seurat)
library(dplyr)

# mock data to check if rds in general is working
set.seed(123)  
data <- data.frame(RandomNumbers = runif(100, min = 0, max = 100))
saveRDS(data, "random_numbers.rds")


d <- readRDS("data/RDS_objects_new/ref_data_Macrophage.RDS")
d 
unique(d@meta.data$sample)

head(d@meta.data)



controls <- d@meta.data |>
    filter(sample == 'PBS') |>
   # group_by(rep) |>
    count()
unique(controls$rep)

controls


data <- readRDS("data/ref_data_Monocyte.RDS")
read.delim("data/ref_data_Macrophage.RDS", row.names = NULL)
source("data/ref_data_Macrophage.RDS")


simon <- readRDS("data/RDS_objects_new/simon_2024_matrix.RDS")
head(simon)

mac <- readRDS("data/RDS_objects/ref_data_Macrophage.RDS")
mono <-readRDS("data/RDS_objects/ref_data_Monocyte.RDS")
cdc1 <- readRDS("data/RDS_objects/ref_data_cDC1.RDS")
cdc2 <- readRDS("data/RDS_objects/ref_data_cDC2.RDS")
pdc <- readRDS("data/RDS_objects/ref_data_pDC.RDS")

myeloid_1 <- merge(mac, y = mono)
myeloid_2 <- merge(myeloid_1, y = cdc1)
myeloid_3 <- merge(myeloid_2, y = cdc2)
myeloid_object <- merge(myeloid_3, y = pdc)

saveRDS(myeloid_object, "data/RDS_objects/myeloid_object.RDS")

