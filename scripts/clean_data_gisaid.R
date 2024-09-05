#!/usr/bin/env Rscript

# Carga las bibliotecas necesarias
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "readxl", "lubridate")

# Obtén una lista de todos los archivos .xls en la carpeta "data"
file_list <- list.files(path = "data/raw/gisaid/", pattern = "\\.xls$", full.names = TRUE)

# Lee y combina todos los archivos en un solo dataframe
df_combined <- file_list %>%
  lapply(read_excel) %>%
  bind_rows()

# Selecciona las columnas relevantes
df_combined <- df_combined[,c(5,12,1,26,17,18,53)]

# Renombra las columnas
df_combined <- df_combined %>%
  rename(
    strain = Isolate_Name,
    date = Collection_Date,
    isolate_id = Isolate_Id,
    host = Host,
    domestic_status = Domestic_Status
  )

# Convierte la columna de fecha al formato de fecha de R
df_combined$date <- as.Date(df_combined$date, format = "%Y-%m-%d")

# Añade una columna con el valor "avian_flu"
df_combined$virus <- "avian_flu"

# Separa la columna "Location" en múltiples columnas
df_combined <- df_combined %>%
  separate(Location, into = c("region", "country", "division", "location"), sep = " / ")

# Reordena las columnas
df_combined <- df_combined[,c(1,2,11,3,4,5,6,7,8,9,10)]

# Elimina filas duplicadas en base a "strain"
df_combined <- df_combined %>%
  distinct(strain, .keep_all = TRUE)

# Exporta el dataframe a un archivo TSV
write_tsv(df_combined, "data/clean_metadata.tsv")

# Fin del script
