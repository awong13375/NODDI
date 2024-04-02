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

#HO Cort 

col_list = c("Frontal_Pole","Parahippocampal_Gyrus_posterior_division","Parahippocampal_Gyrus_anterior_division","Frontal_Orbital_Cortex","Cuneal_Cortex","Precuneous_Cortex",
             "Cingulate_Gyrus_posterior_division","Cingulate_Gyrus_anterior_division","Paracingulate_Gyrus","Subcallosal_Cortex","Juxtapositional_Lobule_Cortex_formerly_Supplementa","Frontal_Medial_Cortex",
             "Intracalcarine_Cortex","Lateral_Occipital_Cortex_inferior_division","Lateral_Occipital_Cortex_superior_division","Angular_Gyrus","Inferior_Temporal_Gyrus_temporooccipital_part",
             "Inferior_Temporal_Gyrus_posterior_division","Inferior_Temporal_Gyrus_anterior_division","Middle_Temporal_Gyrus_temporooccipital_part","Middle_Temporal_Gyrus_posterior_division",
             "Middle_Temporal_Gyrus_anterior_division","Superior_Temporal_Gyrus_posterior_division","Superior_Temporal_Gyrus_anterior_division","Temporal_Pole","Precentral_Gyrus","Inferior_Frontal_Gyrus_pars_opercularis",
             "Inferior_Frontal_Gyrus_pars_triangularis","Middle_Frontal_Gyrus","Superior_Frontal_Gyrus","Insular_Cortex","Occipital_Pole","Supracalcarine_Cortex","Planum_Temporale","Heschls_Gyrus_includes_H1_and_H2",
             "Planum_Polare","Parietal_Operculum_Cortex","Central_Opercular_Cortex","Frontal_Operculum_Cortex","Occipital_Fusiform_Gyrus","Temporal_Occipital_Fusiform_Cortex","Temporal_Fusiform_Cortex_posterior_division",
             "Temporal_Fusiform_Cortex_anterior_division","Lingual_Gyrus")
#NDI
Frontal_Pole = gghistogram(HO_cort_NDI_data$Frontal_Pole,
                      main = "Frontal_Pole",
                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Parahippocampal_Gyrus_posterior_division = gghistogram(HO_cort_NDI_data$Parahippocampal_Gyrus_posterior_division,
                         main = "Parahippocampal_Gyrus_posterior_division",
                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Parahippocampal_Gyrus_anterior_division = gghistogram(HO_cort_NDI_data$Parahippocampal_Gyrus_anterior_division,
                           main = "Parahippocampal_Gyrus_anterior_division",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Orbital_Cortex = gghistogram(HO_cort_NDI_data$Frontal_Orbital_Cortex,
                     main = "Frontal_Orbital_Cortex",
                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cuneal_Cortex = gghistogram(HO_cort_NDI_data$Cuneal_Cortex,
                             main = "Cuneal_Cortex",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Precuneous_Cortex = gghistogram(HO_cort_NDI_data$Precuneous_Cortex,
                            main = "Precuneous_Cortex",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulate_Gyrus_posterior_division = gghistogram(HO_cort_NDI_data$Cingulate_Gyrus_posterior_division,
                            main = "Cingulate_Gyrus_posterior_division",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulate_Gyrus_anterior_division = gghistogram(HO_cort_NDI_data$Cingulate_Gyrus_anterior_division,
                       main = "Cingulate_Gyrus_anterior_division",
                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Paracingulate_Gyrus = gghistogram(HO_cort_NDI_data$Paracingulate_Gyrus,
                                                main = "Paracingulate_Gyrus",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Subcallosal_Cortex = gghistogram(HO_cort_NDI_data$Subcallosal_Cortex,
                                                main = "Subcallosal_Cortex",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Juxtapositional_Lobule_Cortex_formerly_Supplementa = gghistogram(HO_cort_NDI_data$Juxtapositional_Lobule_Cortex_formerly_Supplementa,
                                                main = "Juxtapositional_Lobule_Cortex_formerly_Supplementa",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Medial_Cortex = gghistogram(HO_cort_NDI_data$Frontal_Medial_Cortex,
                                                main = "Frontal_Medial_Cortex",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Intracalcarine_Cortex = gghistogram(HO_cort_NDI_data$Intracalcarine_Cortex,
                                                main = "Intracalcarine_Cortex",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Lateral_Occipital_Cortex_inferior_division = gghistogram(HO_cort_NDI_data$Lateral_Occipital_Cortex_inferior_division,
                                                main = "Lateral_Occipital_Cortex_inferior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Lateral_Occipital_Cortex_superior_division = gghistogram(HO_cort_NDI_data$Lateral_Occipital_Cortex_superior_division,
                                                main = "Lateral_Occipital_Cortex_superior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Angular_Gyrus = gghistogram(HO_cort_NDI_data$Angular_Gyrus,
                                                main = "Angular_Gyrus",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Supramarginal_Gyrus_posterior_division = gghistogram(HO_cort_NDI_data$Supramarginal_Gyrus_posterior_division,
                                                main = "Supramarginal_Gyrus_posterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Supramarginal_Gyrus_anterior_division = gghistogram(HO_cort_NDI_data$Supramarginal_Gyrus_anterior_division,
                                                main = "Supramarginal_Gyrus_anterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Postcentral_Gyrus = gghistogram(HO_cort_NDI_data$Postcentral_Gyrus,
                                                main = "Postcentral_Gyrus",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_temporooccipital_part = gghistogram(HO_cort_NDI_data$Inferior_Temporal_Gyrus_temporooccipital_part,
                                                main = "Inferior_Temporal_Gyrus_temporooccipital_part",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_NDI_data$Inferior_Temporal_Gyrus_posterior_division,
                                                main = "Inferior_Temporal_Gyrus_posterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_NDI_data$Inferior_Temporal_Gyrus_anterior_division,
                                                main = "Inferior_Temporal_Gyrus_anterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_temporooccipital_part = gghistogram(HO_cort_NDI_data$Middle_Temporal_Gyrus_temporooccipital_part,
                                                main = "Middle_Temporal_Gyrus_temporooccipital_part",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_NDI_data$Middle_Temporal_Gyrus_posterior_division,
                                                main = "Cingulate_Gyrus_anterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_NDI_data$Middle_Temporal_Gyrus_anterior_division,
                                                main = "Middle_Temporal_Gyrus_anterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_NDI_data$Superior_Temporal_Gyrus_posterior_division,
                                                main = "Superior_Temporal_Gyrus_posterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_NDI_data$Superior_Temporal_Gyrus_anterior_division,
                                                main = "Superior_Temporal_Gyrus_anterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Pole = gghistogram(HO_cort_NDI_data$Temporal_Pole,
                                                main = "Temporal_Pole",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Precentral_Gyrus = gghistogram(HO_cort_NDI_data$Precentral_Gyrus,
                                                main = "Precentral_Gyrus",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Frontal_Gyrus_pars_opercularis = gghistogram(HO_cort_NDI_data$Inferior_Frontal_Gyrus_pars_opercularis,
                                                main = "Inferior_Frontal_Gyrus_pars_opercularis",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Frontal_Gyrus_pars_triangularis = gghistogram(HO_cort_NDI_data$Inferior_Frontal_Gyrus_pars_triangularis,
                                                main = "Inferior_Frontal_Gyrus_pars_triangularis",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Frontal_Gyrus = gghistogram(HO_cort_NDI_data$Middle_Frontal_Gyrus,
                                                main = "Middle_Frontal_Gyrus",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Frontal_Gyrus = gghistogram(HO_cort_NDI_data$Superior_Frontal_Gyrus,
                                                main = "Superior_Frontal_Gyrus",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Insular_Cortex = gghistogram(HO_cort_NDI_data$Insular_Cortex,
                                     main = "Insular_Cortex",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Occipital_Pole = gghistogram(HO_cort_NDI_data$Occipital_Pole,
                                     main = "Occipital_Pole",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Supracalcarine_Cortex = gghistogram(HO_cort_NDI_data$Supracalcarine_Cortex,
                                     main = "Supracalcarine_Cortex",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Planum_Temporale = gghistogram(HO_cort_NDI_data$Planum_Temporale,
                                     main = "Planum_Temporale",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Heschls_Gyrus_includes_H1_and_H2 = gghistogram(HO_cort_NDI_data$Heschls_Gyrus_includes_H1_and_H2,
                                     main = "Heschls_Gyrus_includes_H1_and_H2",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Planum_Polare = gghistogram(HO_cort_NDI_data$Planum_Polare,
                                     main = "Planum_Polare",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Parietal_Operculum_Cortex = gghistogram(HO_cort_NDI_data$Parietal_Operculum_Cortex,
                                     main = "Parietal_Operculum_Cortex",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Central_Opercular_Cortex = gghistogram(HO_cort_NDI_data$Central_Opercular_Cortex,
                                     main = "Central_Opercular_Cortex",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Operculum_Cortex = gghistogram(HO_cort_NDI_data$Frontal_Operculum_Cortex,
                                     main = "Frontal_Operculum_Cortex",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Occipital_Fusiform_Gyrus = gghistogram(HO_cort_NDI_data$Occipital_Fusiform_Gyrus,
                                     main = "Occipital_Fusiform_Gyrus",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Occipital_Fusiform_Cortex = gghistogram(HO_cort_NDI_data$Temporal_Occipital_Fusiform_Cortex,
                                     main = "Temporal_Occipital_Fusiform_Cortex",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Fusiform_Cortex_posterior_division = gghistogram(HO_cort_NDI_data$Temporal_Fusiform_Cortex_posterior_division,
                                     main = "Temporal_Fusiform_Cortex_posterior_division",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Fusiform_Cortex_anterior_division = gghistogram(HO_cort_NDI_data$Temporal_Fusiform_Cortex_anterior_division,
                                     main = "Temporal_Fusiform_Cortex_anterior_division",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Lingual_Gyrus = gghistogram(HO_cort_NDI_data$Lingual_Gyrus,
                                     main = "Lingual_Gyrus",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Frontal_Pole,Parahippocampal_Gyrus_posterior_division,Parahippocampal_Gyrus_anterior_division,Frontal_Orbital_Cortex,Cuneal_Cortex,Precuneous_Cortex,
          Cingulate_Gyrus_posterior_division,Cingulate_Gyrus_anterior_division,Paracingulate_Gyrus,Subcallosal_Cortex,Juxtapositional_Lobule_Cortex_formerly_Supplementa,Frontal_Medial_Cortex,
          Intracalcarine_Cortex,Lateral_Occipital_Cortex_inferior_division,Lateral_Occipital_Cortex_superior_division,Angular_Gyrus,Inferior_Temporal_Gyrus_temporooccipital_part,
          Inferior_Temporal_Gyrus_posterior_division,Inferior_Temporal_Gyrus_anterior_division,Middle_Temporal_Gyrus_temporooccipital_part,Middle_Temporal_Gyrus_posterior_division,
          Middle_Temporal_Gyrus_anterior_division,Superior_Temporal_Gyrus_posterior_division,Superior_Temporal_Gyrus_anterior_division,Temporal_Pole,Precentral_Gyrus,Inferior_Frontal_Gyrus_pars_opercularis,
          Inferior_Frontal_Gyrus_pars_triangularis,Middle_Frontal_Gyrus,Superior_Frontal_Gyrus,Insular_Cortex,Occipital_Pole,Supracalcarine_Cortex,Planum_Temporale,Heschls_Gyrus_includes_H1_and_H2,
          Planum_Polare,Parietal_Operculum_Cortex,Central_Opercular_Cortex,Frontal_Operculum_Cortex,Occipital_Fusiform_Gyrus,Temporal_Occipital_Fusiform_Cortex,Temporal_Fusiform_Cortex_posterior_division,
          Temporal_Fusiform_Cortex_anterior_division,Lingual_Gyrus, 
          labels = c("HO_cort_NDI"),
          ncol = 4, nrow = 12)


#FWF
Frontal_Pole = gghistogram(HO_cort_FWF_data$Frontal_Pole,
                           main = "Frontal_Pole",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Parahippocampal_Gyrus_posterior_division = gghistogram(HO_cort_FWF_data$Parahippocampal_Gyrus_posterior_division,
                                                       main = "Parahippocampal_Gyrus_posterior_division",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Parahippocampal_Gyrus_anterior_division = gghistogram(HO_cort_FWF_data$Parahippocampal_Gyrus_anterior_division,
                                                      main = "Parahippocampal_Gyrus_anterior_division",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Orbital_Cortex = gghistogram(HO_cort_FWF_data$Frontal_Orbital_Cortex,
                                     main = "Frontal_Orbital_Cortex",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cuneal_Cortex = gghistogram(HO_cort_FWF_data$Cuneal_Cortex,
                            main = "Cuneal_Cortex",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Precuneous_Cortex = gghistogram(HO_cort_FWF_data$Precuneous_Cortex,
                                main = "Precuneous_Cortex",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulate_Gyrus_posterior_division = gghistogram(HO_cort_FWF_data$Cingulate_Gyrus_posterior_division,
                                                 main = "Cingulate_Gyrus_posterior_division",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulate_Gyrus_anterior_division = gghistogram(HO_cort_FWF_data$Cingulate_Gyrus_anterior_division,
                                                main = "Cingulate_Gyrus_anterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Paracingulate_Gyrus = gghistogram(HO_cort_FWF_data$Paracingulate_Gyrus,
                                  main = "Paracingulate_Gyrus",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Subcallosal_Cortex = gghistogram(HO_cort_FWF_data$Subcallosal_Cortex,
                                 main = "Subcallosal_Cortex",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Juxtapositional_Lobule_Cortex_formerly_Supplementa = gghistogram(HO_cort_FWF_data$Juxtapositional_Lobule_Cortex_formerly_Supplementa,
                                                                 main = "Juxtapositional_Lobule_Cortex_formerly_Supplementa",
                                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Medial_Cortex = gghistogram(HO_cort_FWF_data$Frontal_Medial_Cortex,
                                    main = "Frontal_Medial_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Intracalcarine_Cortex = gghistogram(HO_cort_FWF_data$Intracalcarine_Cortex,
                                    main = "Intracalcarine_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Lateral_Occipital_Cortex_inferior_division = gghistogram(HO_cort_FWF_data$Lateral_Occipital_Cortex_inferior_division,
                                                         main = "Lateral_Occipital_Cortex_inferior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Lateral_Occipital_Cortex_superior_division = gghistogram(HO_cort_FWF_data$Lateral_Occipital_Cortex_superior_division,
                                                         main = "Lateral_Occipital_Cortex_superior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Angular_Gyrus = gghistogram(HO_cort_FWF_data$Angular_Gyrus,
                            main = "Angular_Gyrus",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Supramarginal_Gyrus_posterior_division = gghistogram(HO_cort_FWF_data$Supramarginal_Gyrus_posterior_division,
                                                     main = "Supramarginal_Gyrus_posterior_division",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Supramarginal_Gyrus_anterior_division = gghistogram(HO_cort_FWF_data$Supramarginal_Gyrus_anterior_division,
                                                    main = "Supramarginal_Gyrus_anterior_division",
                                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Postcentral_Gyrus = gghistogram(HO_cort_FWF_data$Postcentral_Gyrus,
                                main = "Postcentral_Gyrus",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_temporooccipital_part = gghistogram(HO_cort_FWF_data$Inferior_Temporal_Gyrus_temporooccipital_part,
                                                            main = "Inferior_Temporal_Gyrus_temporooccipital_part",
                                                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_FWF_data$Inferior_Temporal_Gyrus_posterior_division,
                                                         main = "Inferior_Temporal_Gyrus_posterior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_FWF_data$Inferior_Temporal_Gyrus_anterior_division,
                                                        main = "Inferior_Temporal_Gyrus_anterior_division",
                                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_temporooccipital_part = gghistogram(HO_cort_FWF_data$Middle_Temporal_Gyrus_temporooccipital_part,
                                                          main = "Middle_Temporal_Gyrus_temporooccipital_part",
                                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_FWF_data$Middle_Temporal_Gyrus_posterior_division,
                                                       main = "Cingulate_Gyrus_anterior_division",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_FWF_data$Middle_Temporal_Gyrus_anterior_division,
                                                      main = "Middle_Temporal_Gyrus_anterior_division",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_FWF_data$Superior_Temporal_Gyrus_posterior_division,
                                                         main = "Superior_Temporal_Gyrus_posterior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_FWF_data$Superior_Temporal_Gyrus_anterior_division,
                                                        main = "Superior_Temporal_Gyrus_anterior_division",
                                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Pole = gghistogram(HO_cort_FWF_data$Temporal_Pole,
                            main = "Temporal_Pole",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Precentral_Gyrus = gghistogram(HO_cort_FWF_data$Precentral_Gyrus,
                               main = "Precentral_Gyrus",
                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Frontal_Gyrus_pars_opercularis = gghistogram(HO_cort_FWF_data$Inferior_Frontal_Gyrus_pars_opercularis,
                                                      main = "Inferior_Frontal_Gyrus_pars_opercularis",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Frontal_Gyrus_pars_triangularis = gghistogram(HO_cort_FWF_data$Inferior_Frontal_Gyrus_pars_triangularis,
                                                       main = "Inferior_Frontal_Gyrus_pars_triangularis",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Frontal_Gyrus = gghistogram(HO_cort_FWF_data$Middle_Frontal_Gyrus,
                                   main = "Middle_Frontal_Gyrus",
                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Frontal_Gyrus = gghistogram(HO_cort_FWF_data$Superior_Frontal_Gyrus,
                                     main = "Superior_Frontal_Gyrus",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Insular_Cortex = gghistogram(HO_cort_FWF_data$Insular_Cortex,
                             main = "Insular_Cortex",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Occipital_Pole = gghistogram(HO_cort_FWF_data$Occipital_Pole,
                             main = "Occipital_Pole",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Supracalcarine_Cortex = gghistogram(HO_cort_FWF_data$Supracalcarine_Cortex,
                                    main = "Supracalcarine_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Planum_Temporale = gghistogram(HO_cort_FWF_data$Planum_Temporale,
                               main = "Planum_Temporale",
                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Heschls_Gyrus_includes_H1_and_H2 = gghistogram(HO_cort_FWF_data$Heschls_Gyrus_includes_H1_and_H2,
                                               main = "Heschls_Gyrus_includes_H1_and_H2",
                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Planum_Polare = gghistogram(HO_cort_FWF_data$Planum_Polare,
                            main = "Planum_Polare",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Parietal_Operculum_Cortex = gghistogram(HO_cort_FWF_data$Parietal_Operculum_Cortex,
                                        main = "Parietal_Operculum_Cortex",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Central_Opercular_Cortex = gghistogram(HO_cort_FWF_data$Central_Opercular_Cortex,
                                       main = "Central_Opercular_Cortex",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Operculum_Cortex = gghistogram(HO_cort_FWF_data$Frontal_Operculum_Cortex,
                                       main = "Frontal_Operculum_Cortex",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Occipital_Fusiform_Gyrus = gghistogram(HO_cort_FWF_data$Occipital_Fusiform_Gyrus,
                                       main = "Occipital_Fusiform_Gyrus",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Occipital_Fusiform_Cortex = gghistogram(HO_cort_FWF_data$Temporal_Occipital_Fusiform_Cortex,
                                                 main = "Temporal_Occipital_Fusiform_Cortex",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Fusiform_Cortex_posterior_division = gghistogram(HO_cort_FWF_data$Temporal_Fusiform_Cortex_posterior_division,
                                                          main = "Temporal_Fusiform_Cortex_posterior_division",
                                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Fusiform_Cortex_anterior_division = gghistogram(HO_cort_FWF_data$Temporal_Fusiform_Cortex_anterior_division,
                                                         main = "Temporal_Fusiform_Cortex_anterior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Lingual_Gyrus = gghistogram(HO_cort_FWF_data$Lingual_Gyrus,
                            main = "Lingual_Gyrus",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Frontal_Pole,Parahippocampal_Gyrus_posterior_division,Parahippocampal_Gyrus_anterior_division,Frontal_Orbital_Cortex,Cuneal_Cortex,Precuneous_Cortex,
          Cingulate_Gyrus_posterior_division,Cingulate_Gyrus_anterior_division,Paracingulate_Gyrus,Subcallosal_Cortex,Juxtapositional_Lobule_Cortex_formerly_Supplementa,Frontal_Medial_Cortex,
          Intracalcarine_Cortex,Lateral_Occipital_Cortex_inferior_division,Lateral_Occipital_Cortex_superior_division,Angular_Gyrus,Inferior_Temporal_Gyrus_temporooccipital_part,
          Inferior_Temporal_Gyrus_posterior_division,Inferior_Temporal_Gyrus_anterior_division,Middle_Temporal_Gyrus_temporooccipital_part,Middle_Temporal_Gyrus_posterior_division,
          Middle_Temporal_Gyrus_anterior_division,Superior_Temporal_Gyrus_posterior_division,Superior_Temporal_Gyrus_anterior_division,Temporal_Pole,Precentral_Gyrus,Inferior_Frontal_Gyrus_pars_opercularis,
          Inferior_Frontal_Gyrus_pars_triangularis,Middle_Frontal_Gyrus,Superior_Frontal_Gyrus,Insular_Cortex,Occipital_Pole,Supracalcarine_Cortex,Planum_Temporale,Heschls_Gyrus_includes_H1_and_H2,
          Planum_Polare,Parietal_Operculum_Cortex,Central_Opercular_Cortex,Frontal_Operculum_Cortex,Occipital_Fusiform_Gyrus,Temporal_Occipital_Fusiform_Cortex,Temporal_Fusiform_Cortex_posterior_division,
          Temporal_Fusiform_Cortex_anterior_division,Lingual_Gyrus, 
          labels = c("HO_cort_FWF"),
          ncol = 4, nrow = 12)



#ODI

Frontal_Pole = gghistogram(HO_cort_ODI_data$Frontal_Pole,
                           main = "Frontal_Pole",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Parahippocampal_Gyrus_posterior_division = gghistogram(HO_cort_ODI_data$Parahippocampal_Gyrus_posterior_division,
                                                       main = "Parahippocampal_Gyrus_posterior_division",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Parahippocampal_Gyrus_anterior_division = gghistogram(HO_cort_ODI_data$Parahippocampal_Gyrus_anterior_division,
                                                      main = "Parahippocampal_Gyrus_anterior_division",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Orbital_Cortex = gghistogram(HO_cort_ODI_data$Frontal_Orbital_Cortex,
                                     main = "Frontal_Orbital_Cortex",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cuneal_Cortex = gghistogram(HO_cort_ODI_data$Cuneal_Cortex,
                            main = "Cuneal_Cortex",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Precuneous_Cortex = gghistogram(HO_cort_ODI_data$Precuneous_Cortex,
                                main = "Precuneous_Cortex",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulate_Gyrus_posterior_division = gghistogram(HO_cort_ODI_data$Cingulate_Gyrus_posterior_division,
                                                 main = "Cingulate_Gyrus_posterior_division",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulate_Gyrus_anterior_division = gghistogram(HO_cort_ODI_data$Cingulate_Gyrus_anterior_division,
                                                main = "Cingulate_Gyrus_anterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Paracingulate_Gyrus = gghistogram(HO_cort_ODI_data$Paracingulate_Gyrus,
                                  main = "Paracingulate_Gyrus",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Subcallosal_Cortex = gghistogram(HO_cort_ODI_data$Subcallosal_Cortex,
                                 main = "Subcallosal_Cortex",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Juxtapositional_Lobule_Cortex_formerly_Supplementa = gghistogram(HO_cort_ODI_data$Juxtapositional_Lobule_Cortex_formerly_Supplementa,
                                                                 main = "Juxtapositional_Lobule_Cortex_formerly_Supplementa",
                                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Medial_Cortex = gghistogram(HO_cort_ODI_data$Frontal_Medial_Cortex,
                                    main = "Frontal_Medial_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Intracalcarine_Cortex = gghistogram(HO_cort_ODI_data$Intracalcarine_Cortex,
                                    main = "Intracalcarine_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Lateral_Occipital_Cortex_inferior_division = gghistogram(HO_cort_ODI_data$Lateral_Occipital_Cortex_inferior_division,
                                                         main = "Lateral_Occipital_Cortex_inferior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Lateral_Occipital_Cortex_superior_division = gghistogram(HO_cort_ODI_data$Lateral_Occipital_Cortex_superior_division,
                                                         main = "Lateral_Occipital_Cortex_superior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Angular_Gyrus = gghistogram(HO_cort_ODI_data$Angular_Gyrus,
                            main = "Angular_Gyrus",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Supramarginal_Gyrus_posterior_division = gghistogram(HO_cort_ODI_data$Supramarginal_Gyrus_posterior_division,
                                                     main = "Supramarginal_Gyrus_posterior_division",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Supramarginal_Gyrus_anterior_division = gghistogram(HO_cort_ODI_data$Supramarginal_Gyrus_anterior_division,
                                                    main = "Supramarginal_Gyrus_anterior_division",
                                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Postcentral_Gyrus = gghistogram(HO_cort_ODI_data$Postcentral_Gyrus,
                                main = "Postcentral_Gyrus",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_temporooccipital_part = gghistogram(HO_cort_ODI_data$Inferior_Temporal_Gyrus_temporooccipital_part,
                                                            main = "Inferior_Temporal_Gyrus_temporooccipital_part",
                                                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_ODI_data$Inferior_Temporal_Gyrus_posterior_division,
                                                         main = "Inferior_Temporal_Gyrus_posterior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_ODI_data$Inferior_Temporal_Gyrus_anterior_division,
                                                        main = "Inferior_Temporal_Gyrus_anterior_division",
                                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_temporooccipital_part = gghistogram(HO_cort_ODI_data$Middle_Temporal_Gyrus_temporooccipital_part,
                                                          main = "Middle_Temporal_Gyrus_temporooccipital_part",
                                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_ODI_data$Middle_Temporal_Gyrus_posterior_division,
                                                       main = "Cingulate_Gyrus_anterior_division",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_ODI_data$Middle_Temporal_Gyrus_anterior_division,
                                                      main = "Middle_Temporal_Gyrus_anterior_division",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_ODI_data$Superior_Temporal_Gyrus_posterior_division,
                                                         main = "Superior_Temporal_Gyrus_posterior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_ODI_data$Superior_Temporal_Gyrus_anterior_division,
                                                        main = "Superior_Temporal_Gyrus_anterior_division",
                                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Pole = gghistogram(HO_cort_ODI_data$Temporal_Pole,
                            main = "Temporal_Pole",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Precentral_Gyrus = gghistogram(HO_cort_ODI_data$Precentral_Gyrus,
                               main = "Precentral_Gyrus",
                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Frontal_Gyrus_pars_opercularis = gghistogram(HO_cort_ODI_data$Inferior_Frontal_Gyrus_pars_opercularis,
                                                      main = "Inferior_Frontal_Gyrus_pars_opercularis",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Frontal_Gyrus_pars_triangularis = gghistogram(HO_cort_ODI_data$Inferior_Frontal_Gyrus_pars_triangularis,
                                                       main = "Inferior_Frontal_Gyrus_pars_triangularis",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Frontal_Gyrus = gghistogram(HO_cort_ODI_data$Middle_Frontal_Gyrus,
                                   main = "Middle_Frontal_Gyrus",
                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Frontal_Gyrus = gghistogram(HO_cort_ODI_data$Superior_Frontal_Gyrus,
                                     main = "Superior_Frontal_Gyrus",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Insular_Cortex = gghistogram(HO_cort_ODI_data$Insular_Cortex,
                             main = "Insular_Cortex",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Occipital_Pole = gghistogram(HO_cort_ODI_data$Occipital_Pole,
                             main = "Occipital_Pole",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Supracalcarine_Cortex = gghistogram(HO_cort_ODI_data$Supracalcarine_Cortex,
                                    main = "Supracalcarine_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Planum_Temporale = gghistogram(HO_cort_ODI_data$Planum_Temporale,
                               main = "Planum_Temporale",
                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Heschls_Gyrus_includes_H1_and_H2 = gghistogram(HO_cort_ODI_data$Heschls_Gyrus_includes_H1_and_H2,
                                               main = "Heschls_Gyrus_includes_H1_and_H2",
                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Planum_Polare = gghistogram(HO_cort_ODI_data$Planum_Polare,
                            main = "Planum_Polare",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Parietal_Operculum_Cortex = gghistogram(HO_cort_ODI_data$Parietal_Operculum_Cortex,
                                        main = "Parietal_Operculum_Cortex",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Central_Opercular_Cortex = gghistogram(HO_cort_ODI_data$Central_Opercular_Cortex,
                                       main = "Central_Opercular_Cortex",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Operculum_Cortex = gghistogram(HO_cort_ODI_data$Frontal_Operculum_Cortex,
                                       main = "Frontal_Operculum_Cortex",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Occipital_Fusiform_Gyrus = gghistogram(HO_cort_ODI_data$Occipital_Fusiform_Gyrus,
                                       main = "Occipital_Fusiform_Gyrus",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Occipital_Fusiform_Cortex = gghistogram(HO_cort_ODI_data$Temporal_Occipital_Fusiform_Cortex,
                                                 main = "Temporal_Occipital_Fusiform_Cortex",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Fusiform_Cortex_posterior_division = gghistogram(HO_cort_ODI_data$Temporal_Fusiform_Cortex_posterior_division,
                                                          main = "Temporal_Fusiform_Cortex_posterior_division",
                                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Fusiform_Cortex_anterior_division = gghistogram(HO_cort_ODI_data$Temporal_Fusiform_Cortex_anterior_division,
                                                         main = "Temporal_Fusiform_Cortex_anterior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Lingual_Gyrus = gghistogram(HO_cort_ODI_data$Lingual_Gyrus,
                            main = "Lingual_Gyrus",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Frontal_Pole,Parahippocampal_Gyrus_posterior_division,Parahippocampal_Gyrus_anterior_division,Frontal_Orbital_Cortex,Cuneal_Cortex,Precuneous_Cortex,
          Cingulate_Gyrus_posterior_division,Cingulate_Gyrus_anterior_division,Paracingulate_Gyrus,Subcallosal_Cortex,Juxtapositional_Lobule_Cortex_formerly_Supplementa,Frontal_Medial_Cortex,
          Intracalcarine_Cortex,Lateral_Occipital_Cortex_inferior_division,Lateral_Occipital_Cortex_superior_division,Angular_Gyrus,Inferior_Temporal_Gyrus_temporooccipital_part,
          Inferior_Temporal_Gyrus_posterior_division,Inferior_Temporal_Gyrus_anterior_division,Middle_Temporal_Gyrus_temporooccipital_part,Middle_Temporal_Gyrus_posterior_division,
          Middle_Temporal_Gyrus_anterior_division,Superior_Temporal_Gyrus_posterior_division,Superior_Temporal_Gyrus_anterior_division,Temporal_Pole,Precentral_Gyrus,Inferior_Frontal_Gyrus_pars_opercularis,
          Inferior_Frontal_Gyrus_pars_triangularis,Middle_Frontal_Gyrus,Superior_Frontal_Gyrus,Insular_Cortex,Occipital_Pole,Supracalcarine_Cortex,Planum_Temporale,Heschls_Gyrus_includes_H1_and_H2,
          Planum_Polare,Parietal_Operculum_Cortex,Central_Opercular_Cortex,Frontal_Operculum_Cortex,Occipital_Fusiform_Gyrus,Temporal_Occipital_Fusiform_Cortex,Temporal_Fusiform_Cortex_posterior_division,
          Temporal_Fusiform_Cortex_anterior_division,Lingual_Gyrus, 
          labels = c("HO_cort_ODI"),
          ncol = 4, nrow = 12)



# HO subcort

col_list = c("Left_Cerebral_White_Matter","Left_Cerebral_Cortex","Left_Lateral_Ventricle","Left_Thalamus","Left_Caudate","Left_Putamen","Left_Pallidum","Brain.Stem","Left_Hippocampus","Left_Amygdala",
             "Left_Accumbens","Right_Cerebral_White_Matter","Right_Cerebral_Cortex","Right_Lateral_Ventricle","Right_Thalamus","Right_Caudate","Right_Putamen","Right_Pallidum","Right_Hippocampus",
             "Right_Amygdala","Right_Accumbens")
#NDI
Left_Cerebral_White_Matter = gghistogram(HO_subcort_NDI_data$Left_Cerebral_White_Matter,
                           main = "Left_Cerebral_White_Matter",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Cerebral_Cortex = gghistogram(HO_subcort_NDI_data$Left_Cerebral_Cortex,
                                                       main = "Left_Cerebral_Cortex",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Lateral_Ventricle = gghistogram(HO_subcort_NDI_data$Left_Lateral_Ventricle,
                                                      main = "Left_Lateral_Ventricle",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Thalamus = gghistogram(HO_subcort_NDI_data$Left_Thalamus,
                                     main = "Left_Thalamus",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Caudate = gghistogram(HO_subcort_NDI_data$Left_Caudate,
                            main = "Left_Caudate",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Putamen = gghistogram(HO_subcort_NDI_data$Left_Putamen,
                                main = "Left_Putamen",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Pallidum = gghistogram(HO_subcort_NDI_data$Left_Pallidum,
                                                 main = "Left_Pallidum",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Brain.Stem = gghistogram(HO_subcort_NDI_data$Brain.Stem,
                                                main = "Brain.Stem",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Hippocampus = gghistogram(HO_subcort_NDI_data$Left_Hippocampus,
                                  main = "Left_Hippocampus",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Amygdala = gghistogram(HO_subcort_NDI_data$Left_Amygdala,
                                 main = "Left_Amygdala",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Accumbens = gghistogram(HO_subcort_NDI_data$Left_Accumbens,
                                                                 main = "Left_Accumbens",
                                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Cerebral_White_Matter = gghistogram(HO_subcort_NDI_data$Right_Cerebral_White_Matter,
                                    main = "Right_Cerebral_White_Matter",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Cerebral_Cortex = gghistogram(HO_subcort_NDI_data$Right_Cerebral_Cortex,
                                    main = "Right_Cerebral_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Lateral_Ventricle = gghistogram(HO_subcort_NDI_data$Right_Lateral_Ventricle,
                                                         main = "Right_Lateral_Ventricle",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Thalamus = gghistogram(HO_subcort_NDI_data$Right_Thalamus,
                                                         main = "Right_Thalamus",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Caudate = gghistogram(HO_subcort_NDI_data$Right_Caudate,
                            main = "Right_Caudate",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Putamen = gghistogram(HO_subcort_NDI_data$Right_Putamen,
                                                     main = "Right_Putamen",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Pallidum = gghistogram(HO_subcort_NDI_data$Right_Pallidum,
                                                    main = "Right_Pallidum",
                                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Hippocampus = gghistogram(HO_subcort_NDI_data$Right_Hippocampus,
                                main = "Right_Hippocampus",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Amygdala = gghistogram(HO_subcort_NDI_data$Right_Amygdala,
                                                            main = "Right_Amygdala",
                                                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Accumbens = gghistogram(HO_subcort_NDI_data$Right_Accumbens,
                                                         main = "Right_Accumbens",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Left_Cerebral_White_Matter,Left_Cerebral_Cortex,Left_Lateral_Ventricle,Left_Thalamus,Left_Caudate,Left_Putamen,Left_Pallidum,Brain.Stem,Left_Hippocampus,Left_Amygdala,
          Left_Accumbens,Right_Cerebral_White_Matter,Right_Cerebral_Cortex,Right_Lateral_Ventricle,Right_Thalamus,Right_Caudate,Right_Putamen,Right_Pallidum,Right_Hippocampus,
          Right_Amygdala,Right_Accumbens, 
          labels = c("HO_subcort_NDI"),
          ncol = 3, nrow = 7)

#FWF
Left_Cerebral_White_Matter = gghistogram(HO_subcort_FWF_data$Left_Cerebral_White_Matter,
                                         main = "Left_Cerebral_White_Matter",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Cerebral_Cortex = gghistogram(HO_subcort_FWF_data$Left_Cerebral_Cortex,
                                   main = "Left_Cerebral_Cortex",
                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Lateral_Ventricle = gghistogram(HO_subcort_FWF_data$Left_Lateral_Ventricle,
                                     main = "Left_Lateral_Ventricle",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Thalamus = gghistogram(HO_subcort_FWF_data$Left_Thalamus,
                            main = "Left_Thalamus",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Caudate = gghistogram(HO_subcort_FWF_data$Left_Caudate,
                           main = "Left_Caudate",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Putamen = gghistogram(HO_subcort_FWF_data$Left_Putamen,
                           main = "Left_Putamen",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Pallidum = gghistogram(HO_subcort_FWF_data$Left_Pallidum,
                            main = "Left_Pallidum",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Brain.Stem = gghistogram(HO_subcort_FWF_data$Brain.Stem,
                         main = "Brain.Stem",
                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Hippocampus = gghistogram(HO_subcort_FWF_data$Left_Hippocampus,
                               main = "Left_Hippocampus",
                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Amygdala = gghistogram(HO_subcort_FWF_data$Left_Amygdala,
                            main = "Left_Amygdala",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Accumbens = gghistogram(HO_subcort_FWF_data$Left_Accumbens,
                             main = "Left_Accumbens",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Cerebral_White_Matter = gghistogram(HO_subcort_FWF_data$Right_Cerebral_White_Matter,
                                          main = "Right_Cerebral_White_Matter",
                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Cerebral_Cortex = gghistogram(HO_subcort_FWF_data$Right_Cerebral_Cortex,
                                    main = "Right_Cerebral_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Lateral_Ventricle = gghistogram(HO_subcort_FWF_data$Right_Lateral_Ventricle,
                                      main = "Right_Lateral_Ventricle",
                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Thalamus = gghistogram(HO_subcort_FWF_data$Right_Thalamus,
                             main = "Right_Thalamus",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Caudate = gghistogram(HO_subcort_FWF_data$Right_Caudate,
                            main = "Right_Caudate",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Putamen = gghistogram(HO_subcort_FWF_data$Right_Putamen,
                            main = "Right_Putamen",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Pallidum = gghistogram(HO_subcort_FWF_data$Right_Pallidum,
                             main = "Right_Pallidum",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Hippocampus = gghistogram(HO_subcort_FWF_data$Right_Hippocampus,
                                main = "Right_Hippocampus",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Amygdala = gghistogram(HO_subcort_FWF_data$Right_Amygdala,
                             main = "Right_Amygdala",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Accumbens = gghistogram(HO_subcort_FWF_data$Right_Accumbens,
                              main = "Right_Accumbens",
                              xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Left_Cerebral_White_Matter,Left_Cerebral_Cortex,Left_Lateral_Ventricle,Left_Thalamus,Left_Caudate,Left_Putamen,Left_Pallidum,Brain.Stem,Left_Hippocampus,Left_Amygdala,
          Left_Accumbens,Right_Cerebral_White_Matter,Right_Cerebral_Cortex,Right_Lateral_Ventricle,Right_Thalamus,Right_Caudate,Right_Putamen,Right_Pallidum,Right_Hippocampus,
          Right_Amygdala,Right_Accumbens, 
          labels = c("HO_subcort_FWF"),
          ncol = 3, nrow = 7)


#ODI
Left_Cerebral_White_Matter = gghistogram(HO_subcort_ODI_data$Left_Cerebral_White_Matter,
                                         main = "Left_Cerebral_White_Matter",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Cerebral_Cortex = gghistogram(HO_subcort_ODI_data$Left_Cerebral_Cortex,
                                   main = "Left_Cerebral_Cortex",
                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Lateral_Ventricle = gghistogram(HO_subcort_ODI_data$Left_Lateral_Ventricle,
                                     main = "Left_Lateral_Ventricle",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Thalamus = gghistogram(HO_subcort_ODI_data$Left_Thalamus,
                            main = "Left_Thalamus",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Caudate = gghistogram(HO_subcort_ODI_data$Left_Caudate,
                           main = "Left_Caudate",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Putamen = gghistogram(HO_subcort_ODI_data$Left_Putamen,
                           main = "Left_Putamen",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Pallidum = gghistogram(HO_subcort_ODI_data$Left_Pallidum,
                            main = "Left_Pallidum",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Brain.Stem = gghistogram(HO_subcort_ODI_data$Brain.Stem,
                         main = "Brain.Stem",
                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Hippocampus = gghistogram(HO_subcort_ODI_data$Left_Hippocampus,
                               main = "Left_Hippocampus",
                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Amygdala = gghistogram(HO_subcort_ODI_data$Left_Amygdala,
                            main = "Left_Amygdala",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Accumbens = gghistogram(HO_subcort_ODI_data$Left_Accumbens,
                             main = "Left_Accumbens",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Cerebral_White_Matter = gghistogram(HO_subcort_ODI_data$Right_Cerebral_White_Matter,
                                          main = "Right_Cerebral_White_Matter",
                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Cerebral_Cortex = gghistogram(HO_subcort_ODI_data$Right_Cerebral_Cortex,
                                    main = "Right_Cerebral_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Lateral_Ventricle = gghistogram(HO_subcort_ODI_data$Right_Lateral_Ventricle,
                                      main = "Right_Lateral_Ventricle",
                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Thalamus = gghistogram(HO_subcort_ODI_data$Right_Thalamus,
                             main = "Right_Thalamus",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Caudate = gghistogram(HO_subcort_ODI_data$Right_Caudate,
                            main = "Right_Caudate",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Putamen = gghistogram(HO_subcort_ODI_data$Right_Putamen,
                            main = "Right_Putamen",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Pallidum = gghistogram(HO_subcort_ODI_data$Right_Pallidum,
                             main = "Right_Pallidum",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Hippocampus = gghistogram(HO_subcort_ODI_data$Right_Hippocampus,
                                main = "Right_Hippocampus",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Amygdala = gghistogram(HO_subcort_ODI_data$Right_Amygdala,
                             main = "Right_Amygdala",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Accumbens = gghistogram(HO_subcort_ODI_data$Right_Accumbens,
                              main = "Right_Accumbens",
                              xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Left_Cerebral_White_Matter,Left_Cerebral_Cortex,Left_Lateral_Ventricle,Left_Thalamus,Left_Caudate,Left_Putamen,Left_Pallidum,Brain.Stem,Left_Hippocampus,Left_Amygdala,
          Left_Accumbens,Right_Cerebral_White_Matter,Right_Cerebral_Cortex,Right_Lateral_Ventricle,Right_Thalamus,Right_Caudate,Right_Putamen,Right_Pallidum,Right_Hippocampus,
          Right_Amygdala,Right_Accumbens, 
          labels = c("HO_subcort_ODI"),
          ncol = 3, nrow = 7)

#JHU labels

col_list = c("Anterior_corona_radiata_L","Anterior_corona_radiata_R","Superior_corona_radiata_R","Superior_corona_radiata_L","Retrolenticular_part_of_internal_capsule_L","Retrolenticular_part_of_internal_capsule_R",
             "Posterior_limb_of_internal_capsule_L","Posterior_limb_of_internal_capsule_R","Anterior_limb_of_internal_capsule_L","Anterior_limb_of_internal_capsule_R","Cerebral_peduncle_L","Cerebral_peduncle_R",
             "Superior_cerebellar_peduncle_L","Superior_cerebellar_peduncle_R","Inferior_cerebellar_peduncle_L","Inferior_cerebellar_peduncle_R","Medial_lemniscus_L","Medial_lemniscus_R","Corticospinal_tract_L",
             "Corticospinal_tract_R","Fornix_column_and_body_of_fornix","Splenium_of_corpus_callosum","Body_of_corpus_callosum","Genu_of_corpus_callosum","Pontine_crossing_tract_a_part_of_MCP","Middle_cerebellar_peduncle",
             "Tapetum_L","Tapetum_R","Uncinate_fasciculus_L","Uncinate_fasciculus_R","Inferior_fronto.occipital_fasciculus_L","Inferior_fronto.occipital_fasciculus_R","Superior_fronto.occipital_fasciculus_L",
             "Superior_fronto.occipital_fasciculus_R","Superior_longitudinal_fasciculus_L","Superior_longitudinal_fasciculus_R","Fornix_cres_L","Fornix_cres_R","Cingulum_hippocampus_L","Cingulum_hippocampus_R",
             "Cingulum_cingulate_gyrus_L","Cingulum_cingulate_gyrus_R","External_capsule_L","External_capsule_R","Sagittal_stratum_R","Sagittal_stratum_L","Posterior_thalamic_radiation_include_optic_radiation_L",
             "Posterior_thalamic_radiation_include_optic_radiation_R","Posterior_corona_radiata_L","Posterior_corona_radiata_R")

#NDI

Anterior_corona_radiata_L = gghistogram(JHU_label_NDI_data$Anterior_corona_radiata_L,
                           main = "Anterior_corona_radiata_L",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_corona_radiata_R = gghistogram(JHU_label_NDI_data$Anterior_corona_radiata_R,
                                                       main = "Anterior_corona_radiata_R",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_corona_radiata_R = gghistogram(JHU_label_NDI_data$Superior_corona_radiata_R,
                                                      main = "Superior_corona_radiata_R",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_corona_radiata_L = gghistogram(JHU_label_NDI_data$Superior_corona_radiata_L,
                                     main = "Superior_corona_radiata_L",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Retrolenticular_part_of_internal_capsule_L = gghistogram(JHU_label_NDI_data$Retrolenticular_part_of_internal_capsule_L,
                            main = "Retrolenticular_part_of_internal_capsule_L",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Retrolenticular_part_of_internal_capsule_R = gghistogram(JHU_label_NDI_data$Retrolenticular_part_of_internal_capsule_R,
                                main = "Retrolenticular_part_of_internal_capsule_R",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_limb_of_internal_capsule_L = gghistogram(JHU_label_NDI_data$Posterior_limb_of_internal_capsule_L,
                                                 main = "Posterior_limb_of_internal_capsule_L",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_limb_of_internal_capsule_R = gghistogram(JHU_label_NDI_data$Posterior_limb_of_internal_capsule_R,
                                                main = "Posterior_limb_of_internal_capsule_R",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_limb_of_internal_capsule_L = gghistogram(JHU_label_NDI_data$Anterior_limb_of_internal_capsule_L,
                                  main = "Anterior_limb_of_internal_capsule_L",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_limb_of_internal_capsule_R = gghistogram(JHU_label_NDI_data$Anterior_limb_of_internal_capsule_R,
                                 main = "Anterior_limb_of_internal_capsule_R",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Cerebral_peduncle_L = gghistogram(JHU_label_NDI_data$Cerebral_peduncle_L,
                                                                 main = "Cerebral_peduncle_L",
                                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Cerebral_peduncle_R = gghistogram(JHU_label_NDI_data$Cerebral_peduncle_R,
                                    main = "Cerebral_peduncle_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_cerebellar_peduncle_L = gghistogram(JHU_label_NDI_data$Superior_cerebellar_peduncle_L,
                                    main = "Superior_cerebellar_peduncle_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_cerebellar_peduncle_R = gghistogram(JHU_label_NDI_data$Superior_cerebellar_peduncle_R,
                                                         main = "Superior_cerebellar_peduncle_R",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_cerebellar_peduncle_L = gghistogram(JHU_label_NDI_data$Inferior_cerebellar_peduncle_L,
                                                         main = "Inferior_cerebellar_peduncle_L",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_cerebellar_peduncle_R = gghistogram(JHU_label_NDI_data$Inferior_cerebellar_peduncle_R,
                            main = "Inferior_cerebellar_peduncle_R",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Medial_lemniscus_L = gghistogram(JHU_label_NDI_data$Medial_lemniscus_L,
                                                     main = "Medial_lemniscus_L",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Medial_lemniscus_R = gghistogram(JHU_label_NDI_data$Medial_lemniscus_R,
                                                    main = "Medial_lemniscus_R",
                                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_L = gghistogram(JHU_label_NDI_data$Corticospinal_tract_L,
                                main = "Corticospinal_tract_L",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_R = gghistogram(JHU_label_NDI_data$Corticospinal_tract_R,
                                                            main = "Corticospinal_tract_R",
                                                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_column_and_body_of_fornix = gghistogram(JHU_label_NDI_data$Fornix_column_and_body_of_fornix,
                                                         main = "Fornix_column_and_body_of_fornix",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Splenium_of_corpus_callosum = gghistogram(JHU_label_NDI_data$Splenium_of_corpus_callosum,
                                                        main = "Splenium_of_corpus_callosum",
                                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Body_of_corpus_callosum = gghistogram(JHU_label_NDI_data$Body_of_corpus_callosum,
                                                          main = "Body_of_corpus_callosum",
                                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Genu_of_corpus_callosum = gghistogram(JHU_label_NDI_data$Genu_of_corpus_callosum,
                                                       main = "Genu_of_corpus_callosum",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Pontine_crossing_tract_a_part_of_MCP = gghistogram(JHU_label_NDI_data$Pontine_crossing_tract_a_part_of_MCP,
                                                      main = "Pontine_crossing_tract_a_part_of_MCP",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_cerebellar_peduncle = gghistogram(JHU_label_NDI_data$Middle_cerebellar_peduncle,
                                                         main = "Middle_cerebellar_peduncle",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Tapetum_L = gghistogram(JHU_label_NDI_data$Tapetum_L,
                                                        main = "Tapetum_L",
                                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Tapetum_R = gghistogram(JHU_label_NDI_data$Tapetum_R,
                            main = "Tapetum_R",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_L = gghistogram(JHU_label_NDI_data$Uncinate_fasciculus_L,
                               main = "Uncinate_fasciculus_L",
                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_R = gghistogram(JHU_label_NDI_data$Uncinate_fasciculus_R,
                                                      main = "Uncinate_fasciculus_R",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_L = gghistogram(JHU_label_NDI_data$Inferior_fronto.occipital_fasciculus_L,
                                                       main = "Inferior_fronto.occipital_fasciculus_L",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_R = gghistogram(JHU_label_NDI_data$Inferior_fronto.occipital_fasciculus_R,
                                   main = "Inferior_fronto.occipital_fasciculus_R",
                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_fronto.occipital_fasciculus_L = gghistogram(JHU_label_NDI_data$Superior_fronto.occipital_fasciculus_L,
                                     main = "Superior_fronto.occipital_fasciculus_L",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_fronto.occipital_fasciculus_R = gghistogram(JHU_label_NDI_data$Superior_fronto.occipital_fasciculus_R,
                             main = "Superior_fronto.occipital_fasciculus_R",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_L = gghistogram(JHU_label_NDI_data$Superior_longitudinal_fasciculus_L,
                             main = "Superior_longitudinal_fasciculus_L",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_R = gghistogram(JHU_label_NDI_data$Superior_longitudinal_fasciculus_R,
                                    main = "Superior_longitudinal_fasciculus_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_cres_L = gghistogram(JHU_label_NDI_data$Fornix_cres_L,
                               main = "Fornix_cres_L",
                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_cres_R = gghistogram(JHU_label_NDI_data$Fornix_cres_R,
                                               main = "Fornix_cres_R",
                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_L = gghistogram(JHU_label_NDI_data$Cingulum_hippocampus_L,
                            main = "Cingulum_hippocampus_L",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_R = gghistogram(JHU_label_NDI_data$Cingulum_hippocampus_R,
                                        main = "Cingulum_hippocampus_R",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_L = gghistogram(JHU_label_NDI_data$Cingulum_cingulate_gyrus_L,
                                       main = "Cingulum_cingulate_gyrus_L",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_R = gghistogram(JHU_label_NDI_data$Cingulum_cingulate_gyrus_R,
                                       main = "Cingulum_cingulate_gyrus_R",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

External_capsule_L = gghistogram(JHU_label_NDI_data$External_capsule_L,
                                       main = "External_capsule_L",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

External_capsule_R = gghistogram(JHU_label_NDI_data$External_capsule_R,
                                                 main = "External_capsule_R",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Sagittal_stratum_R = gghistogram(JHU_label_NDI_data$Sagittal_stratum_R,
                                                          main = "Sagittal_stratum_R",
                                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Sagittal_stratum_L = gghistogram(JHU_label_NDI_data$Sagittal_stratum_L,
                                                         main = "Sagittal_stratum_L",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_thalamic_radiation_include_optic_radiation_L = gghistogram(JHU_label_NDI_data$Posterior_thalamic_radiation_include_optic_radiation_L,
                            main = "Posterior_thalamic_radiation_include_optic_radiation_L",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_thalamic_radiation_include_optic_radiation_R = gghistogram(JHU_label_NDI_data$Posterior_thalamic_radiation_include_optic_radiation_R,
                            main = "Posterior_thalamic_radiation_include_optic_radiation_R",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()


Posterior_corona_radiata_L = gghistogram(JHU_label_NDI_data$Posterior_corona_radiata_L,
                            main = "Posterior_corona_radiata_L",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_corona_radiata_R = gghistogram(JHU_label_NDI_data$Posterior_corona_radiata_R,
                                         main = "Posterior_corona_radiata_R",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()


plot_grid(Frontal_Pole,Parahippocampal_Gyrus_posterior_division,Parahippocampal_Gyrus_anterior_division,Frontal_Orbital_Cortex,Cuneal_Cortex,Precuneous_Cortex,
          Cingulate_Gyrus_posterior_division,Cingulate_Gyrus_anterior_division,Paracingulate_Gyrus,Subcallosal_Cortex,Juxtapositional_Lobule_Cortex_formerly_Supplementa,Frontal_Medial_Cortex,
          Intracalcarine_Cortex,Lateral_Occipital_Cortex_inferior_division,Lateral_Occipital_Cortex_superior_division,Angular_Gyrus,Inferior_Temporal_Gyrus_temporooccipital_part,
          Inferior_Temporal_Gyrus_posterior_division,Inferior_Temporal_Gyrus_anterior_division,Middle_Temporal_Gyrus_temporooccipital_part,Middle_Temporal_Gyrus_posterior_division,
          Middle_Temporal_Gyrus_anterior_division,Superior_Temporal_Gyrus_posterior_division,Superior_Temporal_Gyrus_anterior_division,Temporal_Pole,Precentral_Gyrus,Inferior_Frontal_Gyrus_pars_opercularis,
          Inferior_Frontal_Gyrus_pars_triangularis,Middle_Frontal_Gyrus,Superior_Frontal_Gyrus,Insular_Cortex,Occipital_Pole,Supracalcarine_Cortex,Planum_Temporale,Heschls_Gyrus_includes_H1_and_H2,
          Planum_Polare,Parietal_Operculum_Cortex,Central_Opercular_Cortex,Frontal_Operculum_Cortex,Occipital_Fusiform_Gyrus,Temporal_Occipital_Fusiform_Cortex,Temporal_Fusiform_Cortex_posterior_division,
          Temporal_Fusiform_Cortex_anterior_division,Lingual_Gyrus, 
          labels = c("JHU_label_NDI"),
          ncol = 5, nrow = 8)

#FWF

Anterior_corona_radiata_L = gghistogram(JHU_label_FWF_data$Anterior_corona_radiata_L,
                                        main = "Anterior_corona_radiata_L",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_corona_radiata_R = gghistogram(JHU_label_FWF_data$Anterior_corona_radiata_R,
                                        main = "Anterior_corona_radiata_R",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_corona_radiata_R = gghistogram(JHU_label_FWF_data$Superior_corona_radiata_R,
                                        main = "Superior_corona_radiata_R",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_corona_radiata_L = gghistogram(JHU_label_FWF_data$Superior_corona_radiata_L,
                                        main = "Superior_corona_radiata_L",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Retrolenticular_part_of_internal_capsule_L = gghistogram(JHU_label_FWF_data$Retrolenticular_part_of_internal_capsule_L,
                                                         main = "Retrolenticular_part_of_internal_capsule_L",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Retrolenticular_part_of_internal_capsule_R = gghistogram(JHU_label_FWF_data$Retrolenticular_part_of_internal_capsule_R,
                                                         main = "Retrolenticular_part_of_internal_capsule_R",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_limb_of_internal_capsule_L = gghistogram(JHU_label_FWF_data$Posterior_limb_of_internal_capsule_L,
                                                   main = "Posterior_limb_of_internal_capsule_L",
                                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_limb_of_internal_capsule_R = gghistogram(JHU_label_FWF_data$Posterior_limb_of_internal_capsule_R,
                                                   main = "Posterior_limb_of_internal_capsule_R",
                                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_limb_of_internal_capsule_L = gghistogram(JHU_label_FWF_data$Anterior_limb_of_internal_capsule_L,
                                                  main = "Anterior_limb_of_internal_capsule_L",
                                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_limb_of_internal_capsule_R = gghistogram(JHU_label_FWF_data$Anterior_limb_of_internal_capsule_R,
                                                  main = "Anterior_limb_of_internal_capsule_R",
                                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Cerebral_peduncle_L = gghistogram(JHU_label_FWF_data$Cerebral_peduncle_L,
                                  main = "Cerebral_peduncle_L",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Cerebral_peduncle_R = gghistogram(JHU_label_FWF_data$Cerebral_peduncle_R,
                                  main = "Cerebral_peduncle_R",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_cerebellar_peduncle_L = gghistogram(JHU_label_FWF_data$Superior_cerebellar_peduncle_L,
                                             main = "Superior_cerebellar_peduncle_L",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_cerebellar_peduncle_R = gghistogram(JHU_label_FWF_data$Superior_cerebellar_peduncle_R,
                                             main = "Superior_cerebellar_peduncle_R",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_cerebellar_peduncle_L = gghistogram(JHU_label_FWF_data$Inferior_cerebellar_peduncle_L,
                                             main = "Inferior_cerebellar_peduncle_L",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_cerebellar_peduncle_R = gghistogram(JHU_label_FWF_data$Inferior_cerebellar_peduncle_R,
                                             main = "Inferior_cerebellar_peduncle_R",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Medial_lemniscus_L = gghistogram(JHU_label_FWF_data$Medial_lemniscus_L,
                                 main = "Medial_lemniscus_L",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Medial_lemniscus_R = gghistogram(JHU_label_FWF_data$Medial_lemniscus_R,
                                 main = "Medial_lemniscus_R",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_L = gghistogram(JHU_label_FWF_data$Corticospinal_tract_L,
                                    main = "Corticospinal_tract_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_R = gghistogram(JHU_label_FWF_data$Corticospinal_tract_R,
                                    main = "Corticospinal_tract_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_column_and_body_of_fornix = gghistogram(JHU_label_FWF_data$Fornix_column_and_body_of_fornix,
                                               main = "Fornix_column_and_body_of_fornix",
                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Splenium_of_corpus_callosum = gghistogram(JHU_label_FWF_data$Splenium_of_corpus_callosum,
                                          main = "Splenium_of_corpus_callosum",
                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Body_of_corpus_callosum = gghistogram(JHU_label_FWF_data$Body_of_corpus_callosum,
                                      main = "Body_of_corpus_callosum",
                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Genu_of_corpus_callosum = gghistogram(JHU_label_FWF_data$Genu_of_corpus_callosum,
                                      main = "Genu_of_corpus_callosum",
                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Pontine_crossing_tract_a_part_of_MCP = gghistogram(JHU_label_FWF_data$Pontine_crossing_tract_a_part_of_MCP,
                                                   main = "Pontine_crossing_tract_a_part_of_MCP",
                                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_cerebellar_peduncle = gghistogram(JHU_label_FWF_data$Middle_cerebellar_peduncle,
                                         main = "Middle_cerebellar_peduncle",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Tapetum_L = gghistogram(JHU_label_FWF_data$Tapetum_L,
                        main = "Tapetum_L",
                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Tapetum_R = gghistogram(JHU_label_FWF_data$Tapetum_R,
                        main = "Tapetum_R",
                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_L = gghistogram(JHU_label_FWF_data$Uncinate_fasciculus_L,
                                    main = "Uncinate_fasciculus_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_R = gghistogram(JHU_label_FWF_data$Uncinate_fasciculus_R,
                                    main = "Uncinate_fasciculus_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_L = gghistogram(JHU_label_FWF_data$Inferior_fronto.occipital_fasciculus_L,
                                                     main = "Inferior_fronto.occipital_fasciculus_L",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_R = gghistogram(JHU_label_FWF_data$Inferior_fronto.occipital_fasciculus_R,
                                                     main = "Inferior_fronto.occipital_fasciculus_R",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_fronto.occipital_fasciculus_L = gghistogram(JHU_label_FWF_data$Superior_fronto.occipital_fasciculus_L,
                                                     main = "Superior_fronto.occipital_fasciculus_L",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_fronto.occipital_fasciculus_R = gghistogram(JHU_label_FWF_data$Superior_fronto.occipital_fasciculus_R,
                                                     main = "Superior_fronto.occipital_fasciculus_R",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_L = gghistogram(JHU_label_FWF_data$Superior_longitudinal_fasciculus_L,
                                                 main = "Superior_longitudinal_fasciculus_L",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_R = gghistogram(JHU_label_FWF_data$Superior_longitudinal_fasciculus_R,
                                                 main = "Superior_longitudinal_fasciculus_R",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_cres_L = gghistogram(JHU_label_FWF_data$Fornix_cres_L,
                            main = "Fornix_cres_L",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_cres_R = gghistogram(JHU_label_FWF_data$Fornix_cres_R,
                            main = "Fornix_cres_R",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_L = gghistogram(JHU_label_FWF_data$Cingulum_hippocampus_L,
                                     main = "Cingulum_hippocampus_L",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_R = gghistogram(JHU_label_FWF_data$Cingulum_hippocampus_R,
                                     main = "Cingulum_hippocampus_R",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_L = gghistogram(JHU_label_FWF_data$Cingulum_cingulate_gyrus_L,
                                         main = "Cingulum_cingulate_gyrus_L",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_R = gghistogram(JHU_label_FWF_data$Cingulum_cingulate_gyrus_R,
                                         main = "Cingulum_cingulate_gyrus_R",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

External_capsule_L = gghistogram(JHU_label_FWF_data$External_capsule_L,
                                 main = "External_capsule_L",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

External_capsule_R = gghistogram(JHU_label_FWF_data$External_capsule_R,
                                 main = "External_capsule_R",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Sagittal_stratum_R = gghistogram(JHU_label_FWF_data$Sagittal_stratum_R,
                                 main = "Sagittal_stratum_R",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Sagittal_stratum_L = gghistogram(JHU_label_FWF_data$Sagittal_stratum_L,
                                 main = "Sagittal_stratum_L",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_thalamic_radiation_include_optic_radiation_L = gghistogram(JHU_label_FWF_data$Posterior_thalamic_radiation_include_optic_radiation_L,
                                                                     main = "Posterior_thalamic_radiation_include_optic_radiation_L",
                                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_thalamic_radiation_include_optic_radiation_R = gghistogram(JHU_label_FWF_data$Posterior_thalamic_radiation_include_optic_radiation_R,
                                                                     main = "Posterior_thalamic_radiation_include_optic_radiation_R",
                                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()


Posterior_corona_radiata_L = gghistogram(JHU_label_FWF_data$Posterior_corona_radiata_L,
                                         main = "Posterior_corona_radiata_L",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_corona_radiata_R = gghistogram(JHU_label_FWF_data$Posterior_corona_radiata_R,
                                         main = "Posterior_corona_radiata_R",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()


plot_grid(Frontal_Pole,Parahippocampal_Gyrus_posterior_division,Parahippocampal_Gyrus_anterior_division,Frontal_Orbital_Cortex,Cuneal_Cortex,Precuneous_Cortex,
          Cingulate_Gyrus_posterior_division,Cingulate_Gyrus_anterior_division,Paracingulate_Gyrus,Subcallosal_Cortex,Juxtapositional_Lobule_Cortex_formerly_Supplementa,Frontal_Medial_Cortex,
          Intracalcarine_Cortex,Lateral_Occipital_Cortex_inferior_division,Lateral_Occipital_Cortex_superior_division,Angular_Gyrus,Inferior_Temporal_Gyrus_temporooccipital_part,
          Inferior_Temporal_Gyrus_posterior_division,Inferior_Temporal_Gyrus_anterior_division,Middle_Temporal_Gyrus_temporooccipital_part,Middle_Temporal_Gyrus_posterior_division,
          Middle_Temporal_Gyrus_anterior_division,Superior_Temporal_Gyrus_posterior_division,Superior_Temporal_Gyrus_anterior_division,Temporal_Pole,Precentral_Gyrus,Inferior_Frontal_Gyrus_pars_opercularis,
          Inferior_Frontal_Gyrus_pars_triangularis,Middle_Frontal_Gyrus,Superior_Frontal_Gyrus,Insular_Cortex,Occipital_Pole,Supracalcarine_Cortex,Planum_Temporale,Heschls_Gyrus_includes_H1_and_H2,
          Planum_Polare,Parietal_Operculum_Cortex,Central_Opercular_Cortex,Frontal_Operculum_Cortex,Occipital_Fusiform_Gyrus,Temporal_Occipital_Fusiform_Cortex,Temporal_Fusiform_Cortex_posterior_division,
          Temporal_Fusiform_Cortex_anterior_division,Lingual_Gyrus, 
          labels = c("JHU_label_FWF"),
          ncol = 5, nrow = 8)


#ODI

Anterior_corona_radiata_L = gghistogram(JHU_label_ODI_data$Anterior_corona_radiata_L,
                                        main = "Anterior_corona_radiata_L",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_corona_radiata_R = gghistogram(JHU_label_ODI_data$Anterior_corona_radiata_R,
                                        main = "Anterior_corona_radiata_R",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_corona_radiata_R = gghistogram(JHU_label_ODI_data$Superior_corona_radiata_R,
                                        main = "Superior_corona_radiata_R",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_corona_radiata_L = gghistogram(JHU_label_ODI_data$Superior_corona_radiata_L,
                                        main = "Superior_corona_radiata_L",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Retrolenticular_part_of_internal_capsule_L = gghistogram(JHU_label_ODI_data$Retrolenticular_part_of_internal_capsule_L,
                                                         main = "Retrolenticular_part_of_internal_capsule_L",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Retrolenticular_part_of_internal_capsule_R = gghistogram(JHU_label_ODI_data$Retrolenticular_part_of_internal_capsule_R,
                                                         main = "Retrolenticular_part_of_internal_capsule_R",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_limb_of_internal_capsule_L = gghistogram(JHU_label_ODI_data$Posterior_limb_of_internal_capsule_L,
                                                   main = "Posterior_limb_of_internal_capsule_L",
                                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_limb_of_internal_capsule_R = gghistogram(JHU_label_ODI_data$Posterior_limb_of_internal_capsule_R,
                                                   main = "Posterior_limb_of_internal_capsule_R",
                                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_limb_of_internal_capsule_L = gghistogram(JHU_label_ODI_data$Anterior_limb_of_internal_capsule_L,
                                                  main = "Anterior_limb_of_internal_capsule_L",
                                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_limb_of_internal_capsule_R = gghistogram(JHU_label_ODI_data$Anterior_limb_of_internal_capsule_R,
                                                  main = "Anterior_limb_of_internal_capsule_R",
                                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Cerebral_peduncle_L = gghistogram(JHU_label_ODI_data$Cerebral_peduncle_L,
                                  main = "Cerebral_peduncle_L",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Cerebral_peduncle_R = gghistogram(JHU_label_ODI_data$Cerebral_peduncle_R,
                                  main = "Cerebral_peduncle_R",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_cerebellar_peduncle_L = gghistogram(JHU_label_ODI_data$Superior_cerebellar_peduncle_L,
                                             main = "Superior_cerebellar_peduncle_L",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_cerebellar_peduncle_R = gghistogram(JHU_label_ODI_data$Superior_cerebellar_peduncle_R,
                                             main = "Superior_cerebellar_peduncle_R",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_cerebellar_peduncle_L = gghistogram(JHU_label_ODI_data$Inferior_cerebellar_peduncle_L,
                                             main = "Inferior_cerebellar_peduncle_L",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_cerebellar_peduncle_R = gghistogram(JHU_label_ODI_data$Inferior_cerebellar_peduncle_R,
                                             main = "Inferior_cerebellar_peduncle_R",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Medial_lemniscus_L = gghistogram(JHU_label_ODI_data$Medial_lemniscus_L,
                                 main = "Medial_lemniscus_L",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Medial_lemniscus_R = gghistogram(JHU_label_ODI_data$Medial_lemniscus_R,
                                 main = "Medial_lemniscus_R",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_L = gghistogram(JHU_label_ODI_data$Corticospinal_tract_L,
                                    main = "Corticospinal_tract_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_R = gghistogram(JHU_label_ODI_data$Corticospinal_tract_R,
                                    main = "Corticospinal_tract_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_column_and_body_of_fornix = gghistogram(JHU_label_ODI_data$Fornix_column_and_body_of_fornix,
                                               main = "Fornix_column_and_body_of_fornix",
                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Splenium_of_corpus_callosum = gghistogram(JHU_label_ODI_data$Splenium_of_corpus_callosum,
                                          main = "Splenium_of_corpus_callosum",
                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Body_of_corpus_callosum = gghistogram(JHU_label_ODI_data$Body_of_corpus_callosum,
                                      main = "Body_of_corpus_callosum",
                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Genu_of_corpus_callosum = gghistogram(JHU_label_ODI_data$Genu_of_corpus_callosum,
                                      main = "Genu_of_corpus_callosum",
                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Pontine_crossing_tract_a_part_of_MCP = gghistogram(JHU_label_ODI_data$Pontine_crossing_tract_a_part_of_MCP,
                                                   main = "Pontine_crossing_tract_a_part_of_MCP",
                                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_cerebellar_peduncle = gghistogram(JHU_label_ODI_data$Middle_cerebellar_peduncle,
                                         main = "Middle_cerebellar_peduncle",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Tapetum_L = gghistogram(JHU_label_ODI_data$Tapetum_L,
                        main = "Tapetum_L",
                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Tapetum_R = gghistogram(JHU_label_ODI_data$Tapetum_R,
                        main = "Tapetum_R",
                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_L = gghistogram(JHU_label_ODI_data$Uncinate_fasciculus_L,
                                    main = "Uncinate_fasciculus_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_R = gghistogram(JHU_label_ODI_data$Uncinate_fasciculus_R,
                                    main = "Uncinate_fasciculus_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_L = gghistogram(JHU_label_ODI_data$Inferior_fronto.occipital_fasciculus_L,
                                                     main = "Inferior_fronto.occipital_fasciculus_L",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_R = gghistogram(JHU_label_ODI_data$Inferior_fronto.occipital_fasciculus_R,
                                                     main = "Inferior_fronto.occipital_fasciculus_R",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_fronto.occipital_fasciculus_L = gghistogram(JHU_label_ODI_data$Superior_fronto.occipital_fasciculus_L,
                                                     main = "Superior_fronto.occipital_fasciculus_L",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_fronto.occipital_fasciculus_R = gghistogram(JHU_label_ODI_data$Superior_fronto.occipital_fasciculus_R,
                                                     main = "Superior_fronto.occipital_fasciculus_R",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_L = gghistogram(JHU_label_ODI_data$Superior_longitudinal_fasciculus_L,
                                                 main = "Superior_longitudinal_fasciculus_L",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_R = gghistogram(JHU_label_ODI_data$Superior_longitudinal_fasciculus_R,
                                                 main = "Superior_longitudinal_fasciculus_R",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_cres_L = gghistogram(JHU_label_ODI_data$Fornix_cres_L,
                            main = "Fornix_cres_L",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_cres_R = gghistogram(JHU_label_ODI_data$Fornix_cres_R,
                            main = "Fornix_cres_R",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_L = gghistogram(JHU_label_ODI_data$Cingulum_hippocampus_L,
                                     main = "Cingulum_hippocampus_L",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_R = gghistogram(JHU_label_ODI_data$Cingulum_hippocampus_R,
                                     main = "Cingulum_hippocampus_R",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_L = gghistogram(JHU_label_ODI_data$Cingulum_cingulate_gyrus_L,
                                         main = "Cingulum_cingulate_gyrus_L",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_R = gghistogram(JHU_label_ODI_data$Cingulum_cingulate_gyrus_R,
                                         main = "Cingulum_cingulate_gyrus_R",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

External_capsule_L = gghistogram(JHU_label_ODI_data$External_capsule_L,
                                 main = "External_capsule_L",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

External_capsule_R = gghistogram(JHU_label_ODI_data$External_capsule_R,
                                 main = "External_capsule_R",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Sagittal_stratum_R = gghistogram(JHU_label_ODI_data$Sagittal_stratum_R,
                                 main = "Sagittal_stratum_R",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Sagittal_stratum_L = gghistogram(JHU_label_ODI_data$Sagittal_stratum_L,
                                 main = "Sagittal_stratum_L",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_thalamic_radiation_include_optic_radiation_L = gghistogram(JHU_label_ODI_data$Posterior_thalamic_radiation_include_optic_radiation_L,
                                                                     main = "Posterior_thalamic_radiation_include_optic_radiation_L",
                                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_thalamic_radiation_include_optic_radiation_R = gghistogram(JHU_label_ODI_data$Posterior_thalamic_radiation_include_optic_radiation_R,
                                                                     main = "Posterior_thalamic_radiation_include_optic_radiation_R",
                                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()


Posterior_corona_radiata_L = gghistogram(JHU_label_ODI_data$Posterior_corona_radiata_L,
                                         main = "Posterior_corona_radiata_L",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_corona_radiata_R = gghistogram(JHU_label_ODI_data$Posterior_corona_radiata_R,
                                         main = "Posterior_corona_radiata_R",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()


plot_grid(Frontal_Pole,Parahippocampal_Gyrus_posterior_division,Parahippocampal_Gyrus_anterior_division,Frontal_Orbital_Cortex,Cuneal_Cortex,Precuneous_Cortex,
          Cingulate_Gyrus_posterior_division,Cingulate_Gyrus_anterior_division,Paracingulate_Gyrus,Subcallosal_Cortex,Juxtapositional_Lobule_Cortex_formerly_Supplementa,Frontal_Medial_Cortex,
          Intracalcarine_Cortex,Lateral_Occipital_Cortex_inferior_division,Lateral_Occipital_Cortex_superior_division,Angular_Gyrus,Inferior_Temporal_Gyrus_temporooccipital_part,
          Inferior_Temporal_Gyrus_posterior_division,Inferior_Temporal_Gyrus_anterior_division,Middle_Temporal_Gyrus_temporooccipital_part,Middle_Temporal_Gyrus_posterior_division,
          Middle_Temporal_Gyrus_anterior_division,Superior_Temporal_Gyrus_posterior_division,Superior_Temporal_Gyrus_anterior_division,Temporal_Pole,Precentral_Gyrus,Inferior_Frontal_Gyrus_pars_opercularis,
          Inferior_Frontal_Gyrus_pars_triangularis,Middle_Frontal_Gyrus,Superior_Frontal_Gyrus,Insular_Cortex,Occipital_Pole,Supracalcarine_Cortex,Planum_Temporale,Heschls_Gyrus_includes_H1_and_H2,
          Planum_Polare,Parietal_Operculum_Cortex,Central_Opercular_Cortex,Frontal_Operculum_Cortex,Occipital_Fusiform_Gyrus,Temporal_Occipital_Fusiform_Cortex,Temporal_Fusiform_Cortex_posterior_division,
          Temporal_Fusiform_Cortex_anterior_division,Lingual_Gyrus, 
          labels = c("JHU_label_ODI"),
          ncol = 5, nrow = 8)


#JHU tract

col_list = c("Superior_longitudinal_fasciculus_temporal_part_R","Superior_longitudinal_fasciculus_temporal_part_L","Uncinate_fasciculus_R","Uncinate_fasciculus_L","Superior_longitudinal_fasciculus_R",
             "Superior_longitudinal_fasciculus_L","Inferior_longitudinal_fasciculus_R","Inferior_longitudinal_fasciculus_L","Inferior_fronto.occipital_fasciculus_R","Inferior_fronto.occipital_fasciculus_L",
             "Forceps_minor","Forceps_major","Cingulum_hippocampus_R","Cingulum_hippocampus_L","Cingulum_cingulate_gyrus_R","Cingulum_cingulate_gyrus_L","Corticospinal_tract_R","Corticospinal_tract_L",
             "Anterior_thalamic_radiation_R","Anterior_thalamic_radiation_L")

#NDI

Superior_longitudinal_fasciculus_temporal_part_R = gghistogram(JHU_tract_NDI_data$Superior_longitudinal_fasciculus_temporal_part_R,
                                        main = "Superior_longitudinal_fasciculus_temporal_part_R",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_temporal_part_L = gghistogram(JHU_tract_NDI_data$Superior_longitudinal_fasciculus_temporal_part_L,
                                        main = "Superior_longitudinal_fasciculus_temporal_part_L",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_R = gghistogram(JHU_tract_NDI_data$Uncinate_fasciculus_R,
                                        main = "Uncinate_fasciculus_R",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_L = gghistogram(JHU_tract_NDI_data$Uncinate_fasciculus_L,
                                        main = "Uncinate_fasciculus_L",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_R = gghistogram(JHU_tract_NDI_data$Superior_longitudinal_fasciculus_R,
                                                         main = "Superior_longitudinal_fasciculus_R",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_L = gghistogram(JHU_tract_NDI_data$Superior_longitudinal_fasciculus_L,
                                                         main = "Superior_longitudinal_fasciculus_L",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_longitudinal_fasciculus_R = gghistogram(JHU_tract_NDI_data$Inferior_longitudinal_fasciculus_R,
                                                   main = "Inferior_longitudinal_fasciculus_R",
                                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_longitudinal_fasciculus_L = gghistogram(JHU_tract_NDI_data$Inferior_longitudinal_fasciculus_L,
                                                   main = "Inferior_longitudinal_fasciculus_L",
                                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_R = gghistogram(JHU_tract_NDI_data$Inferior_fronto.occipital_fasciculus_R,
                                                  main = "Inferior_fronto.occipital_fasciculus_R",
                                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_L = gghistogram(JHU_tract_NDI_data$Inferior_fronto.occipital_fasciculus_L,
                                                  main = "Inferior_fronto.occipital_fasciculus_L",
                                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Forceps_minor = gghistogram(JHU_tract_NDI_data$Forceps_minor,
                                  main = "Forceps_minor",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Forceps_major = gghistogram(JHU_tract_NDI_data$Forceps_major,
                                  main = "Forceps_major",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_R = gghistogram(JHU_tract_NDI_data$Cingulum_hippocampus_R,
                                             main = "Cingulum_hippocampus_R",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_L = gghistogram(JHU_tract_NDI_data$Cingulum_hippocampus_L,
                                             main = "Cingulum_hippocampus_L",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_R = gghistogram(JHU_tract_NDI_data$Cingulum_cingulate_gyrus_R,
                                             main = "Cingulum_cingulate_gyrus_R",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_L = gghistogram(JHU_tract_NDI_data$Cingulum_cingulate_gyrus_L,
                                             main = "Cingulum_cingulate_gyrus_L",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_R = gghistogram(JHU_tract_NDI_data$Corticospinal_tract_R,
                                 main = "Corticospinal_tract_R",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_L = gghistogram(JHU_tract_NDI_data$Corticospinal_tract_L,
                                 main = "Corticospinal_tract_L",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_thalamic_radiation_R = gghistogram(JHU_tract_NDI_data$Anterior_thalamic_radiation_R,
                                    main = "Anterior_thalamic_radiation_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_thalamic_radiation_L = gghistogram(JHU_tract_NDI_data$Anterior_thalamic_radiation_L,
                                    main = "Anterior_thalamic_radiation_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Superior_longitudinal_fasciculus_temporal_part_R,Superior_longitudinal_fasciculus_temporal_part_L,Uncinate_fasciculus_R,Uncinate_fasciculus_L,Superior_longitudinal_fasciculus_R,
          Superior_longitudinal_fasciculus_L,Inferior_longitudinal_fasciculus_R,Inferior_longitudinal_fasciculus_L,Inferior_fronto.occipital_fasciculus_R,Inferior_fronto.occipital_fasciculus_L,
          Forceps_minor,Forceps_major,Cingulum_hippocampus_R,Cingulum_hippocampus_L,Cingulum_cingulate_gyrus_R,Cingulum_cingulate_gyrus_L,Corticospinal_tract_R,Corticospinal_tract_L,
          Anterior_thalamic_radiation_R,Anterior_thalamic_radiation_L, 
          labels = c("JHU_tract_NDI"),
          ncol = 5, nrow = 4)



# FWF

Superior_longitudinal_fasciculus_temporal_part_R = gghistogram(JHU_tract_FWF_data$Superior_longitudinal_fasciculus_temporal_part_R,
                                                               main = "Superior_longitudinal_fasciculus_temporal_part_R",
                                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_temporal_part_L = gghistogram(JHU_tract_FWF_data$Superior_longitudinal_fasciculus_temporal_part_L,
                                                               main = "Superior_longitudinal_fasciculus_temporal_part_L",
                                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_R = gghistogram(JHU_tract_FWF_data$Uncinate_fasciculus_R,
                                    main = "Uncinate_fasciculus_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_L = gghistogram(JHU_tract_FWF_data$Uncinate_fasciculus_L,
                                    main = "Uncinate_fasciculus_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_R = gghistogram(JHU_tract_FWF_data$Superior_longitudinal_fasciculus_R,
                                                 main = "Superior_longitudinal_fasciculus_R",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_L = gghistogram(JHU_tract_FWF_data$Superior_longitudinal_fasciculus_L,
                                                 main = "Superior_longitudinal_fasciculus_L",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_longitudinal_fasciculus_R = gghistogram(JHU_tract_FWF_data$Inferior_longitudinal_fasciculus_R,
                                                 main = "Inferior_longitudinal_fasciculus_R",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_longitudinal_fasciculus_L = gghistogram(JHU_tract_FWF_data$Inferior_longitudinal_fasciculus_L,
                                                 main = "Inferior_longitudinal_fasciculus_L",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_R = gghistogram(JHU_tract_FWF_data$Inferior_fronto.occipital_fasciculus_R,
                                                     main = "Inferior_fronto.occipital_fasciculus_R",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_L = gghistogram(JHU_tract_FWF_data$Inferior_fronto.occipital_fasciculus_L,
                                                     main = "Inferior_fronto.occipital_fasciculus_L",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Forceps_minor = gghistogram(JHU_tract_FWF_data$Forceps_minor,
                            main = "Forceps_minor",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Forceps_major = gghistogram(JHU_tract_FWF_data$Forceps_major,
                            main = "Forceps_major",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_R = gghistogram(JHU_tract_FWF_data$Cingulum_hippocampus_R,
                                     main = "Cingulum_hippocampus_R",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_L = gghistogram(JHU_tract_FWF_data$Cingulum_hippocampus_L,
                                     main = "Cingulum_hippocampus_L",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_R = gghistogram(JHU_tract_FWF_data$Cingulum_cingulate_gyrus_R,
                                         main = "Cingulum_cingulate_gyrus_R",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_L = gghistogram(JHU_tract_FWF_data$Cingulum_cingulate_gyrus_L,
                                         main = "Cingulum_cingulate_gyrus_L",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_R = gghistogram(JHU_tract_FWF_data$Corticospinal_tract_R,
                                    main = "Corticospinal_tract_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_L = gghistogram(JHU_tract_FWF_data$Corticospinal_tract_L,
                                    main = "Corticospinal_tract_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_thalamic_radiation_R = gghistogram(JHU_tract_FWF_data$Anterior_thalamic_radiation_R,
                                            main = "Anterior_thalamic_radiation_R",
                                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_thalamic_radiation_L = gghistogram(JHU_tract_FWF_data$Anterior_thalamic_radiation_L,
                                            main = "Anterior_thalamic_radiation_L",
                                            xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Superior_longitudinal_fasciculus_temporal_part_R,Superior_longitudinal_fasciculus_temporal_part_L,Uncinate_fasciculus_R,Uncinate_fasciculus_L,Superior_longitudinal_fasciculus_R,
          Superior_longitudinal_fasciculus_L,Inferior_longitudinal_fasciculus_R,Inferior_longitudinal_fasciculus_L,Inferior_fronto.occipital_fasciculus_R,Inferior_fronto.occipital_fasciculus_L,
          Forceps_minor,Forceps_major,Cingulum_hippocampus_R,Cingulum_hippocampus_L,Cingulum_cingulate_gyrus_R,Cingulum_cingulate_gyrus_L,Corticospinal_tract_R,Corticospinal_tract_L,
          Anterior_thalamic_radiation_R,Anterior_thalamic_radiation_L, 
          labels = c("JHU_tract_FWF"),
          ncol = 5, nrow = 4)

# ODI

Superior_longitudinal_fasciculus_temporal_part_R = gghistogram(JHU_tract_ODI_data$Superior_longitudinal_fasciculus_temporal_part_R,
                                                               main = "Superior_longitudinal_fasciculus_temporal_part_R",
                                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_temporal_part_L = gghistogram(JHU_tract_ODI_data$Superior_longitudinal_fasciculus_temporal_part_L,
                                                               main = "Superior_longitudinal_fasciculus_temporal_part_L",
                                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_R = gghistogram(JHU_tract_ODI_data$Uncinate_fasciculus_R,
                                    main = "Uncinate_fasciculus_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_L = gghistogram(JHU_tract_ODI_data$Uncinate_fasciculus_L,
                                    main = "Uncinate_fasciculus_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_R = gghistogram(JHU_tract_ODI_data$Superior_longitudinal_fasciculus_R,
                                                 main = "Superior_longitudinal_fasciculus_R",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_L = gghistogram(JHU_tract_ODI_data$Superior_longitudinal_fasciculus_L,
                                                 main = "Superior_longitudinal_fasciculus_L",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_longitudinal_fasciculus_R = gghistogram(JHU_tract_ODI_data$Inferior_longitudinal_fasciculus_R,
                                                 main = "Inferior_longitudinal_fasciculus_R",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_longitudinal_fasciculus_L = gghistogram(JHU_tract_ODI_data$Inferior_longitudinal_fasciculus_L,
                                                 main = "Inferior_longitudinal_fasciculus_L",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_R = gghistogram(JHU_tract_ODI_data$Inferior_fronto.occipital_fasciculus_R,
                                                     main = "Inferior_fronto.occipital_fasciculus_R",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_L = gghistogram(JHU_tract_ODI_data$Inferior_fronto.occipital_fasciculus_L,
                                                     main = "Inferior_fronto.occipital_fasciculus_L",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Forceps_minor = gghistogram(JHU_tract_ODI_data$Forceps_minor,
                            main = "Forceps_minor",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Forceps_major = gghistogram(JHU_tract_ODI_data$Forceps_major,
                            main = "Forceps_major",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_R = gghistogram(JHU_tract_ODI_data$Cingulum_hippocampus_R,
                                     main = "Cingulum_hippocampus_R",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_L = gghistogram(JHU_tract_ODI_data$Cingulum_hippocampus_L,
                                     main = "Cingulum_hippocampus_L",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_R = gghistogram(JHU_tract_ODI_data$Cingulum_cingulate_gyrus_R,
                                         main = "Cingulum_cingulate_gyrus_R",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_L = gghistogram(JHU_tract_ODI_data$Cingulum_cingulate_gyrus_L,
                                         main = "Cingulum_cingulate_gyrus_L",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_R = gghistogram(JHU_tract_ODI_data$Corticospinal_tract_R,
                                    main = "Corticospinal_tract_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_L = gghistogram(JHU_tract_ODI_data$Corticospinal_tract_L,
                                    main = "Corticospinal_tract_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_thalamic_radiation_R = gghistogram(JHU_tract_ODI_data$Anterior_thalamic_radiation_R,
                                            main = "Anterior_thalamic_radiation_R",
                                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_thalamic_radiation_L = gghistogram(JHU_tract_ODI_data$Anterior_thalamic_radiation_L,
                                            main = "Anterior_thalamic_radiation_L",
                                            xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Superior_longitudinal_fasciculus_temporal_part_R,Superior_longitudinal_fasciculus_temporal_part_L,Uncinate_fasciculus_R,Uncinate_fasciculus_L,Superior_longitudinal_fasciculus_R,
          Superior_longitudinal_fasciculus_L,Inferior_longitudinal_fasciculus_R,Inferior_longitudinal_fasciculus_L,Inferior_fronto.occipital_fasciculus_R,Inferior_fronto.occipital_fasciculus_L,
          Forceps_minor,Forceps_major,Cingulum_hippocampus_R,Cingulum_hippocampus_L,Cingulum_cingulate_gyrus_R,Cingulum_cingulate_gyrus_L,Corticospinal_tract_R,Corticospinal_tract_L,
          Anterior_thalamic_radiation_R,Anterior_thalamic_radiation_L, 
          labels = c("JHU_tract_ODI"),
          ncol = 5, nrow = 4)




