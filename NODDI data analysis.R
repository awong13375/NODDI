library(irr)
library(stringr)
library(blandr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggpubr)
library(ggeasy)

setwd("C:/WSL2_dir/NODDISAH_11/NODDI_processing")

#load data
mni_NDI_data = read.csv('MNI_NDI_raw_intensities.csv')
mni_FWF_data = read.csv('MNI_FWF_raw_intensities.csv')
mni_ODI_data = read.csv('MNI_ODI_raw_intensities.csv')

df2 <- mutate_all(mni_NDI_data, function(x) as.numeric(x))


mni_NDI_data[is.nan(mni_NDI_data)] <- NA
mni_NDI_data[mni_NDI_data == 0] <- NA


caudate = ggdensity(mni_NDI_data$Caudate,
          main = "caudate",
          xlab = "Intensity"
          ) + ggeasy::easy_center_title()
 
cerebellum = ggdensity(mni_NDI_data$Cerebellum,
                    main = "Cerebellum",
                    xlab = "Intensity"
                    ) + ggeasy::easy_center_title()

Frontal_Lobe = ggdensity(mni_NDI_data$Frontal_Lobe,
                       main = "Frontal_Lobe",
                       xlab = "Intensit y"
                       ) + ggeasy::easy_center_title()

Insula = ggdensity(mni_NDI_data$Insula,
                         main = "Insula",
                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Occipital_Lobe = ggdensity(mni_NDI_data$Occipital_Lobe,
                         main = "Occipital_Lobe",
                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Parietal_Lobe = ggdensity(mni_NDI_data$Parietal_Lobe,
                         main = "Parietal_Lobe",
                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Lobe = ggdensity(mni_NDI_data$Temporal_Lobe,
                         main = "Temporal_Lobe",
                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Thalamus = ggdensity(mni_NDI_data$Thalamus,
                         main = "Thalamus",
                         xlab = "Intensity"
) + ggeasy::easy_center_title()


plot_grid(caudate, cerebellum, Frontal_Lobe, Insula, Occipital_Lobe, Parietal_Lobe, Temporal_Lobe, Thalamus, 
          labels = c("FWF"),
          ncol = 4, nrow = 2)

