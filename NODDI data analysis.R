library(irr)
library(stringr)
library(blandr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggpubr)

setwd("C:/WSL2_dir/NODDISAH_11/NODDI_processing")

#load data
mni_NDI_data = read.csv('MNI_NDI_raw_intensities.csv')
mni_FWF_data = read.csv('MNI_FWF_raw_intensities.csv')
mni_ODI_data = read.csv('MNI_ODI_raw_intensities.csv')

caudate = ggdensity(mni_NDI_data$Caudate,
          main = "Density plot caudate intensities",
          xlab = "Intensity"
          )

cerebellum = ggdensity(mni_NDI_data$Cerebellum,
                    main = "Density plot Cerebellum intensities",
                    xlab = "Intensity"
                    )

cerebellum = ggdensity(mni_NDI_data$Cerebellum,
                       main = "Density plot Cerebellum intensities",
                       xlab = "Intensity"
                       )


plot_grid(Predict_QT_ABC2, Predict_QT_DM, 
          labels = c("A","B"),
          ncol = 2, nrow = 1)
