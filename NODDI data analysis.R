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

HO_cort_NDI_data = read.csv('HO_cort_NDI_raw_intensities.csv')
HO_cort_FWF_data = read.csv('HO_cort_FWF_raw_intensities.csv')
HO_cort_ODI_data = read.csv('HO_cort_ODI_raw_intensities.csv')

HO_subcort_NDI_data = read.csv('HO_subcort_NDI_raw_intensities.csv')
HO_subcort_FWF_data = read.csv('HO_subcort_FWF_raw_intensities.csv')
HO_subcort_ODI_data = read.csv('HO_subcort_ODI_raw_intensities.csv')

JHU_label_NDI_data = read.csv('JHU_label_NDI_raw_intensities.csv')
JHU_label_FWF_data = read.csv('JHU_label_FWF_raw_intensities.csv')
JHU_label_ODI_data = read.csv('JHU_label_ODI_raw_intensities.csv')

JHU_tract_NDI_data = read.csv('JHU_tract_NDI_raw_intensities.csv')
JHU_tract_FWF_data = read.csv('JHU_tract_FWF_raw_intensities.csv')
JHU_tract_ODI_data = read.csv('JHU_tract_ODI_raw_intensities.csv')

mni_NDI_data <- mutate_all(mni_NDI_data, function(x) as.numeric(as.character(x)))
mni_FWF_data <- mutate_all(mni_FWF_data, function(x) as.numeric(as.character(x)))
mni_ODI_data <- mutate_all(mni_ODI_data, function(x) as.numeric(as.character(x)))

HO_cort_NDI_data <- mutate_all(HO_cort_NDI_data, function(x) as.numeric(as.character(x)))
HO_cort_FWF_data <- mutate_all(HO_cort_FWF_data, function(x) as.numeric(as.character(x)))
HO_cort_ODI_data <- mutate_all(HO_cort_ODI_data, function(x) as.numeric(as.character(x)))

HO_subcort_NDI_data <- mutate_all(HO_subcort_NDI_data, function(x) as.numeric(as.character(x)))
HO_subcort_FWF_data <- mutate_all(HO_subcort_FWF_data, function(x) as.numeric(as.character(x)))
HO_subcort_ODI_data <- mutate_all(HO_subcort_ODI_data, function(x) as.numeric(as.character(x)))

JHU_label_NDI_data <- mutate_all(JHU_label_NDI_data, function(x) as.numeric(as.character(x)))
JHU_label_FWF_data <- mutate_all(JHU_label_FWF_data, function(x) as.numeric(as.character(x)))
JHU_label_ODI_data <- mutate_all(JHU_label_ODI_data, function(x) as.numeric(as.character(x)))

JHU_tract_NDI_data <- mutate_all(JHU_tract_NDI_data, function(x) as.numeric(as.character(x)))
JHU_tract_FWF_data <- mutate_all(JHU_tract_FWF_data, function(x) as.numeric(as.character(x)))
JHU_tract_ODI_data <- mutate_all(JHU_tract_ODI_data, function(x) as.numeric(as.character(x)))


mni_NDI_data[mni_NDI_data == 0] <- NA
mni_NDI_data <- mni_NDI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

mni_FWF_data[mni_FWF_data == 0] <- NA
mni_FWF_data <- mni_FWF_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

mni_ODI_data[mni_ODI_data == 0] <- NA
mni_ODI_data <- mni_ODI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))


HO_cort_NDI_data[HO_cort_NDI_data == 0] <- NA
HO_cort_NDI_data <- HO_cort_NDI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

HO_cort_FWF_data[HO_cort_FWF_data == 0] <- NA
HO_cort_FWF_data <- HO_cort_FWF_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

HO_cort_ODI_data[HO_cort_ODI_data == 0] <- NA
HO_cort_ODI_data <- HO_cort_ODI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))


HO_subcort_NDI_data[HO_subcort_NDI_data == 0] <- NA
HO_subcort_NDI_data <- HO_subcort_NDI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

HO_subcort_FWF_data[HO_subcort_FWF_data == 0] <- NA
HO_subcort_FWF_data <- HO_subcort_FWF_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

HO_subcort_ODI_data[HO_subcort_ODI_data == 0] <- NA
HO_subcort_ODI_data <- HO_subcort_ODI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))


JHU_label_NDI_data[JHU_label_NDI_data == 0] <- NA
JHU_label_NDI_data <- JHU_label_NDI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

JHU_label_FWF_data[JHU_label_FWF_data == 0] <- NA
JHU_label_FWF_data <- JHU_label_FWF_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

JHU_label_ODI_data[JHU_label_ODI_data == 0] <- NA
JHU_label_ODI_data <- JHU_label_ODI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))


JHU_tract_NDI_data[JHU_tract_NDI_data == 0] <- NA
JHU_tract_NDI_data <- JHU_tract_NDI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

JHU_tract_FWF_data[JHU_tract_FWF_data == 0] <- NA
JHU_tract_FWF_data <- JHU_tract_FWF_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

JHU_tract_ODI_data[JHU_tract_ODI_data == 0] <- NA
JHU_tract_ODI_data <- JHU_tract_ODI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))


# MNI

col_list = c("Caudate","Cerebellum","Frontal_Lobe","Insula","Occipital_Lobe","Parietal_Lobe","Temporal_Lobe","Thalamus")
for (n in 1:length(col_list)){
    col = col_list[n]
    print (paste("Max:", col, "=", max(mni_NDI_data[col], na.rm=T)))
    print (paste("Min:", col, "=", min(mni_NDI_data[col], na.rm=T)))
    
  
}
for (n in 1:length(col_list)){
  col = col_list[n]
  print (paste("Max:", col, "=", max(mni_FWF_data[col], na.rm=T)))
  print (paste("Min:", col, "=", min(mni_FWF_data[col], na.rm=T)))
  
  
}
for (n in 1:length(col_list)){
  col = col_list[n]
  print (paste("Max:", col, "=", max(mni_ODI_data[col], na.rm=T)))
  print (paste("Min:", col, "=", min(mni_ODI_data[col], na.rm=T)))
  
  
}


#create plots
#NDI
Caudate = gghistogram(mni_NDI_data$Caudate,
                    main = "Caudate",
                    xlab = "Intensity"
                   ) + ggeasy::easy_center_title()
 
Cerebellum = gghistogram(mni_NDI_data$Cerebellum,
                       main = "Cerebellum",
                       xlab = "Intensity"
                      ) + ggeasy::easy_center_title()

Frontal_Lobe = gghistogram(mni_NDI_data$Frontal_Lobe,
                         main = "Frontal_Lobe",
                         xlab = "Intensity"
                        ) + ggeasy::easy_center_title()

Insula = gghistogram(mni_NDI_data$Insula,
                   main = "Insula",
                   xlab = "Intensity"
                  ) + ggeasy::easy_center_title()

Occipital_Lobe = gghistogram(mni_NDI_data$Occipital_Lobe,
                           main = "Occipital_Lobe",
                           xlab = "Intensity"
                          ) + ggeasy::easy_center_title()

Parietal_Lobe = gghistogram(mni_NDI_data$Parietal_Lobe,
                          main = "Parietal_Lobe",
                          xlab = "Intensity"
                         ) + ggeasy::easy_center_title()

Temporal_Lobe = gghistogram(mni_NDI_data$Temporal_Lobe,
                          main = "Temporal_Lobe",
                          xlab = "Intensity"
                         ) + ggeasy::easy_center_title()

Thalamus = gghistogram(mni_NDI_data$Thalamus,
                     main = "Thalamus",
                     xlab = "Intensity"
                    ) + ggeasy::easy_center_title()


plot_grid(Caudate, Cerebellum, Frontal_Lobe, Insula, Occipital_Lobe, Parietal_Lobe, Temporal_Lobe, Thalamus, 
          labels = c("MNI_NDI"),
          ncol = 4, nrow = 2)

#FWF
Caudate = gghistogram(mni_FWF_data$Caudate,
                    main = "caudate",
                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Cerebellum = gghistogram(mni_FWF_data$Cerebellum,
                       main = "Cerebellum",
                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Lobe = gghistogram(mni_FWF_data$Frontal_Lobe,
                         main = "Frontal_Lobe",
                         xlab = "Intensit y"
) + ggeasy::easy_center_title()

Insula = gghistogram(mni_FWF_data$Insula,
                   main = "Insula",
                   xlab = "Intensity",
                   xlim = c(0,0.5)
) + ggeasy::easy_center_title()

Occipital_Lobe = gghistogram(mni_FWF_data$Occipital_Lobe,
                           main = "Occipital_Lobe",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Parietal_Lobe = gghistogram(mni_FWF_data$Parietal_Lobe,
                          main = "Parietal_Lobe",
                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Lobe = gghistogram(mni_FWF_data$Temporal_Lobe,
                          main = "Temporal_Lobe",
                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Thalamus = gghistogram(mni_FWF_data$Thalamus,
                     main = "Thalamus",
                     xlab = "Intensity"
) + ggeasy::easy_center_title()


plot_grid(Caudate, Cerebellum, Frontal_Lobe, Insula, Occipital_Lobe, Parietal_Lobe, Temporal_Lobe, Thalamus, 
          labels = c("MNI_FWF"),
          ncol = 4, nrow = 2)

#ODI
Caudate = gghistogram(mni_ODI_data$Caudate,
                    main = "caudate",
                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Cerebellum = gghistogram(mni_ODI_data$Cerebellum,
                       main = "Cerebellum",
                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Lobe = gghistogram(mni_ODI_data$Frontal_Lobe,
                         main = "Frontal_Lobe",
                         xlab = "Intensit y"
) + ggeasy::easy_center_title()

Insula = gghistogram(mni_ODI_data$Insula,
                   main = "Insula",
                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Occipital_Lobe = gghistogram(mni_ODI_data$Occipital_Lobe,
                           main = "Occipital_Lobe",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Parietal_Lobe = gghistogram(mni_ODI_data$Parietal_Lobe,
                          main = "Parietal_Lobe",
                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Lobe = gghistogram(mni_ODI_data$Temporal_Lobe,
                          main = "Temporal_Lobe",
                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Thalamus = gghistogram(mni_ODI_data$Thalamus,
                     main = "Thalamus",
                     xlab = "Intensity",
                     xlim = c(0,0.8)
) + ggeasy::easy_center_title()


plot_grid(Caudate, Cerebellum, Frontal_Lobe, Insula, Occipital_Lobe, Parietal_Lobe, Temporal_Lobe, Thalamus, 
          labels = c("MNI_ODI"),
          ncol = 4, nrow = 2)


col_list = c("Frontal_Pole","Parahippocampal_Gyrus_posterior_division","Parahippocampal_Gyrus_anterior_division","Frontal_Orbital_Cortex","Cuneal_Cortex","Precuneous_Cortex",
             "Cingulate_Gyrus_posterior_division","Cingulate_Gyrus_anterior_division","Paracingulate_Gyrus","Subcallosal_Cortex","Juxtapositional_Lobule_Cortex_formerly_Supplementa","Frontal_Medial_Cortex",
             "Intracalcarine_Cortex","Lateral_Occipital_Cortex_inferior_division","Lateral_Occipital_Cortex_superior_division","Angular_Gyrus","Inferior_Temporal_Gyrus_temporooccipital_part",
             "Inferior_Temporal_Gyrus_posterior_division","Inferior_Temporal_Gyrus_anterior_division","Middle_Temporal_Gyrus_temporooccipital_part","Middle_Temporal_Gyrus_posterior_division",
             "Middle_Temporal_Gyrus_anterior_division","Superior_Temporal_Gyrus_posterior_division","Superior_Temporal_Gyrus_anterior_division","Temporal_Pole","Precentral_Gyrus","Inferior_Frontal_Gyrus_pars_opercularis",
             "Inferior_Frontal_Gyrus_pars_triangularis","Middle_Frontal_Gyrus","Superior_Frontal_Gyrus","Insular_Cortex","Occipital_Pole","Supracalcarine_Cortex","Planum_Temporale","Heschls_Gyrus_includes_H1_and_H2",
             "Planum_Polare","Parietal_Operculum_Cortex","Central_Opercular_Cortex","Frontal_Operculum_Cortex","Occipital_Fusiform_Gyrus","Temporal_Occipital_Fusiform_Cortex","Temporal_Fusiform_Cortex_posterior_division",
             "Temporal_Fusiform_Cortex_anterior_division","Lingual_Gyrus")
#NDI
Caudate = gghistogram(mni_NDI_data$Caudate,
                      main = "Caudate",
                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Cerebellum = gghistogram(mni_NDI_data$Cerebellum,
                         main = "Cerebellum",
                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Lobe = gghistogram(mni_NDI_data$Frontal_Lobe,
                           main = "Frontal_Lobe",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Insula = gghistogram(mni_NDI_data$Insula,
                     main = "Insula",
                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Occipital_Lobe = gghistogram(mni_NDI_data$Occipital_Lobe,
                             main = "Occipital_Lobe",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Parietal_Lobe = gghistogram(mni_NDI_data$Parietal_Lobe,
                            main = "Parietal_Lobe",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Lobe = gghistogram(mni_NDI_data$Temporal_Lobe,
                            main = "Temporal_Lobe",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Thalamus = gghistogram(mni_NDI_data$Thalamus,
                       main = "Thalamus",
                       xlab = "Intensity"
) + ggeasy::easy_center_title()


plot_grid(Caudate, Cerebellum, Frontal_Lobe, Insula, Occipital_Lobe, Parietal_Lobe, Temporal_Lobe, Thalamus, 
          labels = c("MNI_NDI"),
          ncol = 4, nrow = 2)


