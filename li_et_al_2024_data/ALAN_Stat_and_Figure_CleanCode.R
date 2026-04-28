### Behavioural and transgenerational effects of artificial light at night (ALAN) of varying spectral compositions in zebrafish (Danio rerio)
### Weiwei Li，Dongxu Zhang, Qingqing Zou, Aneesh P. H. Bose, Alex Jordan, Erin S. McCallum, Jianghui Bao, Ming Duan




################################## PREPARES ################################----

### Uploading packages-----
library(tidyverse)
library(lme4)
library(glmmTMB)
library(performance)
library(DHARMa)
library(emmeans)
library(mgcv)
library(visreg) 
library(forestplot)
library(patchwork)
library(ggpubr)
library(lmerTest)
library(multcomp)
library(mgcViz)
library(gridExtra)

### Creating a work space----
# set working directory as your first step
setwd("~/GitHub_Projects/ALAN_MA_2026/li_et_al_2024_data") # replace with your working directory
### Prepare datasets-------

# adult individual datasets----
ALAN_individual <- read.csv("ALAN_individual_data.csv", stringsAsFactors = T)
# adult group data set----
ALAN_group <- read.csv("ALAN_group_data.csv", stringsAsFactors = T)
# offspring dataset
ALAN_offspring <- read.csv("ALAN_offspring_data.csv", stringsAsFactors = T)

# set factors and turn ID names to characters

ALAN_individual$lightWavelength <- factor(ALAN_individual$lightWavelength, 
                                          levels=c("dark", "white", "365nm", "420nm", 
                                                   "470nm", "500nm", "520nm", "590nm", 
                                                   "620nm", "635nm", "660nm")) # reordering the factor levels so that the 'control' (dark condition) is first

ALAN_group$lightWavelength <- factor(ALAN_group$lightWavelength, 
                                     levels=c("dark", "white", "365nm", "420nm", 
                                              "470nm", "500nm", "520nm", "590nm", 
                                              "620nm", "635nm", "660nm")) # reordering the factor levels so that the 'control' (dark condition) is first

ALAN_offspring$lightWavelength <- factor(ALAN_offspring$lightWavelength, 
                                         levels=c("dark", "white", "365nm", 
                                                  "420nm", "470nm", "500nm", 
                                                  "520nm", "590nm", "620nm", 
                                                  "635nm", "660nm")) # reordering the factor levels so that the 'control' (dark condition) is first

ALAN_individual <- ALAN_individual %>%
  mutate(individualAdultID = as.character(individualAdultID))

ALAN_group <- ALAN_group %>%
  mutate(dumyID_EachDay_EachLight = as.character(dumyID_EachDay_EachLight))

ALAN_offspring <- ALAN_offspring %>%
  mutate(lavaeID = as.character(lavaeID),
         motherID = as.character(motherID))

# Making wavelength an ordered factor
ALAN_individual <- transform(ALAN_individual, lightWavelength = ordered(lightWavelength))
ALAN_group <- transform(ALAN_group, lightWavelength = ordered(lightWavelength))
ALAN_offspring <- transform(ALAN_offspring, lightWavelength = ordered(lightWavelength))


# separate ALAN_individual to different data
individual_thigmotaxis <- subset(ALAN_individual, behaviorType == "Time_In_Thigmotaxis_Zone")
individual_transition <- subset(ALAN_individual, behaviorType == "Time_In_Transition_Zone")
individual_center <- subset(ALAN_individual, behaviorType == "Time_In_Center_Zone")
individual_distance <- subset(ALAN_individual, behaviorType == "Distance")

# separate ALAN_group to different data
group_thigmotaxis <- subset(ALAN_group, behaviorType == "Time_In_Thigmotaxis_Zone")
group_transition <- subset(ALAN_group, behaviorType == "Time_In_Transition_Zone")
group_center <- subset(ALAN_group, behaviorType == "Time_In_Center_Zone")
group_distance <- subset(ALAN_group, behaviorType == "Distance")
group_inter_individual_distance <- subset(ALAN_group, behaviorType == "Inter_Individual_Distance")

############################ STATISTICAL ANALYSES##########################---- 
### - Part 1: adult individual (GAM)-----------

## -Thigmotaxis model- ----
individual_thigmotaxis.model1 <- gam(Value ~ s(day, by = lightWavelength) + lightWavelength, data = individual_thigmotaxis, method = "REML")
summary(individual_thigmotaxis.model1) # R-sq.(adj) =  0.946

## -transition model- ----
individual_transition.model1 <- gam(Value ~ s(day, by = lightWavelength) + lightWavelength, data = individual_transition, method = "REML")
summary(individual_transition.model1)  # R-sq is 0.914

## -center model- ----
individual_center.model1 <- gam(Value ~ s(day, by = lightWavelength) + lightWavelength, data = individual_center, method = "REML")
summary(individual_center.model1)  # R-sq is 0.764

## -Distance model- ----
individual_distance.model1 <- gam(Value ~ s(day, by = lightWavelength) + lightWavelength, data = individual_distance, method = "REML")
summary(individual_distance.model1)  # R-sq is 0.932

### - Part 2: adult group (GAMM, because trialID has a potential random effect)-----------

## -Thigmotaxis model- ----
group_thigmotaxis.model1 <- gamm(value ~ s(day, by = lightWavelength) + lightWavelength,
                                 random = list(trailID = ~1),
                                 data = group_thigmotaxis, method = "REML")
summary(group_thigmotaxis.model1$gam)  ## good, R-sq.(adj) =  0.964  

## -transition model- ----
group_transition.model1 <- gamm(value ~ s(day, by = lightWavelength) + lightWavelength,
                                random = list(trailID = ~1),
                                data = group_transition, method = "REML")
summary(group_transition.model1$gam)  ## good, R-sq.(adj) =  0.946 

## -center model- ----
group_center.model1 <- gamm(value ~ s(day, by = lightWavelength) + lightWavelength,
                            random = list(trailID = ~1),
                            data = group_center, method = "REML")
summary(group_center.model1$gam)  ## okay, R-sq.(adj) =  0.792  

## -distance model- ----
group_distance.model1 <- gamm(value ~ s(day, by = lightWavelength) + lightWavelength,
                              random = list(trailID = ~1),
                              data = group_distance, method = "REML")
summary(group_distance.model1$gam)  ## okay, R-sq.(adj) =  0.804  

## -inter_individual_distance model- -----
group_inter_individual_distance.model1 <- gamm(value ~ s(day, by = lightWavelength) + lightWavelength,
                                               random = list(trailID = ~1),
                                               data = group_inter_individual_distance, method = "REML")
summary(group_inter_individual_distance.model1$gam)  ## okay, R-sq.(adj) =  0.53

### - Part 3: offspring (LMM)-----------

## -totalMovementDistance model- ----
offspring_totalMovementDistance.model1 <- lmer(totalMovementDistance ~ lightWavelength + (1|motherID), data = ALAN_offspring)
summary(offspring_totalMovementDistance.model1)   # fix effect: all light treatment has big effect; 
emmeans_totalMovementDistance_results <- emmeans(offspring_totalMovementDistance.model1, ~ lightWavelength)
print(summary(emmeans_totalMovementDistance_results))
pairwise_totalMovementDistance_comparisons <- pairs(emmeans_totalMovementDistance_results)
print(summary(pairwise_totalMovementDistance_comparisons))

## -activity model- ----
offspring_activity.model1 <- lmer(activity ~ lightWavelength + (1|motherID), data = ALAN_offspring)
summary(offspring_activity.model1)
emmeans_activity_results <- emmeans(offspring_activity.model1, ~ lightWavelength)
print(summary(emmeans_activity_results))
pairwise_activity_comparisons <- pairs(emmeans_activity_results)
print(summary(pairwise_activity_comparisons))


##################### FIGURES ##################################################

## prepare the bisic of each figure
# set colors
wavelength_colors <- c(
  "dark" = "#505050",     
  "white" = "#C0C0C0",    
  "365nm" = "#7F00FF",    
  "420nm" = "#4B0082",   
  "470nm" = "#0000FF",    
  "500nm" = "#00FFFF",   
  "520nm" = "#00FF00",   
  "590nm" = "#FFFF00",    
  "620nm" = "#FF7F00",    
  "635nm" = "#FF4500",    
  "660nm" = "#FF0000"     
)

# GAM: Use this function to generate predicted data and confidence intervals by simply replacing the dataset and model name
generate_prediction_data <- function(data, model, wavelengths) {
  pred_data_wl <- with(data[data$lightWavelength == wavelengths, ],
                       expand.grid(day = seq(min(day), max(day), length.out = 100),
                                   lightWavelength = wavelengths))
  gam_predictions_wl <- predict(model, newdata = pred_data_wl, type = "response", se.fit = TRUE)
  data.frame(day = pred_data_wl$day,
             lightWavelength = pred_data_wl$lightWavelength,
             fit = gam_predictions_wl$fit,
             se.fit = gam_predictions_wl$se.fit)
}

# GAMM: Use this function to generate predicted data and confidence intervals by simply replacing the dataset and model name
generate_prediction_data_gamm <- function(data, gamm_object, wavelengths) {

  pred_data_wl <- with(data[data$lightWavelength == wavelengths, ],
                       expand.grid(day = seq(min(day), max(day), length.out = 100),
                                   lightWavelength = wavelengths))

  gam_predictions_wl <- predict(gamm_object$gam, newdata = pred_data_wl, type = "response", se.fit = TRUE)
  
  data.frame(day = pred_data_wl$day,
             lightWavelength = pred_data_wl$lightWavelength,
             fit = gam_predictions_wl$fit,
             se.fit = gam_predictions_wl$se.fit)
}

## - Part 1: adult individual - 

## -Thigmotaxis Fig. 1.1- ----

# -1.1a- Time and trend

# Note the replacement of your_data and your_model with your actual data and model names
unique_wavelengths <- unique(individual_thigmotaxis$lightWavelength)    # replace data here
pred_list <- lapply(unique_wavelengths, generate_prediction_data, data = individual_thigmotaxis, model = individual_thigmotaxis.model1)  # replace data and model here
pred_data_all <- do.call(rbind, pred_list)
pred_data_all$CI_low <- pred_data_all$fit - 1.96 * pred_data_all$se.fit
pred_data_all$CI_high <- pred_data_all$fit + 1.96 * pred_data_all$se.fit
pred_data_all$lightWavelength <- factor(pred_data_all$lightWavelength, levels = unique_wavelengths)

# Plotting graphs, including raw data points and predicted confidence intervals
f1.1a <- ggplot() +
  geom_point(data = individual_thigmotaxis, aes(x = day, y = Value, group = lightWavelength, fill = lightWavelength, color = lightWavelength), size = 1, alpha = 0.1) +   #  replace data here
  geom_ribbon(data = pred_data_all, aes(x = day, ymin = CI_low, ymax = CI_high, fill = lightWavelength), alpha = 0.2, show.legend = FALSE) +
  geom_line(data = pred_data_all, aes(x = day, y = fit, color = lightWavelength)) +
  scale_color_manual(values = wavelength_colors) +
  scale_fill_manual(values = wavelength_colors) +
  labs(title = "Trend of Value Over Time by Light Wavelength",
       x = "Day", y = "Time in thigmotaxis zone",  # replace new name for y
       color = "Light Wavelength", fill = "Light Wavelength") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.5, colour = "black")
  ) +
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE), fill = guide_legend(nrow = 1, byrow = TRUE))

# b figure: light wavelength effect size

thigmotaxis_effectSizeLight <- data.frame(
  lightWavelength = c("dark", "white", "365nm", "420nm", "470nm", "500nm", "520nm", "590nm", "620nm", "635nm", "660nm"),
  Estimate = c(0, 55.670, 129.529, 104.974, 154.405, 67.176, 57.535, 69.063, 75.609, 68.392, 54.376), # compare to dark
  Std.Error = c(1.735, 2.454, 2.454, 2.454, 2.454, 2.454, 2.454, 2.454, 2.454, 2.454, 2.454)
)
# calculate CI
thigmotaxis_effectSizeLight$LowerCI <- thigmotaxis_effectSizeLight$Estimate - 1.96 * thigmotaxis_effectSizeLight$Std.Error
thigmotaxis_effectSizeLight$UpperCI <- thigmotaxis_effectSizeLight$Estimate + 1.96 * thigmotaxis_effectSizeLight$Std.Error

f1.1b <- ggplot(thigmotaxis_effectSizeLight, aes(x = lightWavelength, y = Estimate, ymin = LowerCI, ymax = UpperCI, color = lightWavelength)) +
  geom_point(size = 0.5) +
  geom_linerange(aes(ymin = LowerCI, ymax = UpperCI), linewidth = 0.25) +
  geom_point(data = subset(thigmotaxis_effectSizeLight, lightWavelength == "470nm"), size = 1) +  # 加粗470nm的点
  geom_linerange(data = subset(thigmotaxis_effectSizeLight, lightWavelength == "470nm"), ymin = subset(thigmotaxis_effectSizeLight, lightWavelength == "470nm")$LowerCI, ymax = subset(thigmotaxis_effectSizeLight, lightWavelength == "470nm")$UpperCI, linewidth = 0.5) +  # 加粗470nm的线
  geom_text(data = subset(thigmotaxis_effectSizeLight, lightWavelength == "470nm"), aes(label = "470nm"), vjust = 2, size = 3.5, nudge_y = -2) +  # 添加标签并微调位置
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # 添加水平虚线
  scale_color_manual(values = wavelength_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "light wavelength effect size", x = "light wavelength", y = "Partial for light wavelength") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# a and b to figure 1
f1.1 <- f1.1a + inset_element(f1.1b, 0.05, 0.6, 0.4, 0.95)
legend1.1 <- get_legend(f1.1a + theme(legend.position="bottom",
                                      legend.spacing.x = unit(0.1, "cm"),
                                      legend.key.size = unit(0.5, "cm"),
                                      legend.text = element_text(size = 10)))
Fig1.1 <- cowplot::plot_grid(f1.1, legend1.1, ncol=1, rel_heights = c(1, 0.2))
Fig1.1

## -transition Fig. 1.2- ----

# Note the replacement of your_data and your_model with your actual data and model names
unique_wavelengths <- unique(individual_transition$lightWavelength)    # replace data here
pred_list <- lapply(unique_wavelengths, generate_prediction_data, data = individual_transition, model = individual_transition.model1)  # replace data and model here
pred_data_all <- do.call(rbind, pred_list)
pred_data_all$CI_low <- pred_data_all$fit - 1.96 * pred_data_all$se.fit
pred_data_all$CI_high <- pred_data_all$fit + 1.96 * pred_data_all$se.fit
pred_data_all$lightWavelength <- factor(pred_data_all$lightWavelength, levels = unique_wavelengths)

# Plotting graphs, including raw data points and predicted confidence intervals
f1.2a <- ggplot() +
  geom_point(data = individual_transition, aes(x = day, y = Value, group = lightWavelength, fill = lightWavelength, color = lightWavelength), size = 1, alpha = 0.1) +   #  replace data here
  geom_ribbon(data = pred_data_all, aes(x = day, ymin = CI_low, ymax = CI_high, fill = lightWavelength), alpha = 0.2, show.legend = FALSE) +
  geom_line(data = pred_data_all, aes(x = day, y = fit, color = lightWavelength)) +
  scale_color_manual(values = wavelength_colors) +
  scale_fill_manual(values = wavelength_colors) +
  labs(title = "Trend of Value Over Time by Light Wavelength",
       x = "Day", y = "Time in transition zone",  # replace new name for y
       color = "Light Wavelength", fill = "Light Wavelength") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.5, colour = "black")
  ) +
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE), fill = guide_legend(nrow = 1, byrow = TRUE))

# b figure: light wavelength effect size

transition_effectSizeLight <- data.frame(
  lightWavelength = c("dark", "white", "365nm", "420nm", "470nm", "500nm", "520nm", "590nm", "620nm", "635nm", "660nm"),
  Estimate = c(0, -35.461, -50.237, -42.115, -68.952, -20.494, -20.490, -33.195, -41.080, -42.013, -39.330), 
  Std.Error = rep(1.705, 11) 
)

# calculate CI
transition_effectSizeLight$LowerCI <- transition_effectSizeLight$Estimate - 1.96 * transition_effectSizeLight$Std.Error
transition_effectSizeLight$UpperCI <- transition_effectSizeLight$Estimate + 1.96 * transition_effectSizeLight$Std.Error

f1.2b <- ggplot(transition_effectSizeLight, aes(x = lightWavelength, y = Estimate, ymin = LowerCI, ymax = UpperCI, color = lightWavelength)) +
  geom_point(size = 0.5) +
  geom_linerange(aes(ymin = LowerCI, ymax = UpperCI), linewidth = 0.25) +
  geom_point(data = subset(transition_effectSizeLight, lightWavelength == "470nm"), size = 1) +  # make 470nm thicker
  # geom_linerange(data = subset(transition_effectSizeLight, lightWavelength == "470nm"), ymin = subset(transition_effectSizeLight, lightWavelength == "470nm")$LowerCI, ymax = subset(transition_effectSizeLight, lightWavelength == "470nm")$UpperCI, linewidth = 0.5) +  # 加粗470nm的线
  # geom_text(data = subset(transition_effectSizeLight, lightWavelength == "470nm"), aes(label = "470nm"), vjust = 2, size = 3.5, nudge_y = -2) +  # add text
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # add  dash line
  scale_color_manual(values = wavelength_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "light wavelength effect size", x = "light wavelength", y = "Partial for light wavelength") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# a and b to Fig 1.2
f1.2 <- f1.2a + inset_element(f1.2b, 0.05, 0.60, 0.4, 0.95)
legend1.2 <- get_legend(f1.2a + theme(legend.position="bottom",
                                      legend.spacing.x = unit(0.1, "cm"),
                                      legend.key.size = unit(0.5, "cm"),
                                      legend.text = element_text(size = 10)))
Fig1.2 <- cowplot::plot_grid(f1.2, legend1.2, ncol=1, rel_heights = c(1, 0.2))
Fig1.2

## -Center Fig. 1.3- ----

# Note the replacement of your_data and your_model with your actual data and model names
unique_wavelengths <- unique(individual_center$lightWavelength)    # replace data here
pred_list <- lapply(unique_wavelengths, generate_prediction_data, data = individual_center, model = individual_center.model1)  # replace data and model here
pred_data_all <- do.call(rbind, pred_list)
pred_data_all$CI_low <- pred_data_all$fit - 1.96 * pred_data_all$se.fit
pred_data_all$CI_high <- pred_data_all$fit + 1.96 * pred_data_all$se.fit
pred_data_all$lightWavelength <- factor(pred_data_all$lightWavelength, levels = unique_wavelengths)

# Plotting graphs, including raw data points and predicted confidence intervals
f1.3a <- ggplot() +
  geom_point(data = individual_center, aes(x = day, y = Value, group = lightWavelength, fill = lightWavelength, color = lightWavelength), size = 1, alpha = 0.1) +   #  replace data here
  geom_ribbon(data = pred_data_all, aes(x = day, ymin = CI_low, ymax = CI_high, fill = lightWavelength), alpha = 0.2, show.legend = FALSE) +
  geom_line(data = pred_data_all, aes(x = day, y = fit, color = lightWavelength)) +
  scale_color_manual(values = wavelength_colors) +
  scale_fill_manual(values = wavelength_colors) +
  labs(title = "Trend of Value Over Time by Light Wavelength",
       x = "Day", y = "Time in center zone",  # replace new name for y
       color = "Light Wavelength", fill = "Light Wavelength") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.5, colour = "black")
  ) +
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE), fill = guide_legend(nrow = 1, byrow = TRUE))

# b figure: light wavelength effect size

center_effectSizeLight <- data.frame(
  lightWavelength = c("dark", "white", "365nm", "420nm", "470nm", "500nm", "520nm", "590nm", "620nm", "635nm", "660nm"),
  Estimate = c(0, -20.208, -79.292, -62.859, -85.453, -46.682, -37.045, -35.868, -34.529, -26.379, -15.046), 
  Std.Error = rep(2.847, 11) 
)

# calculate CI
center_effectSizeLight$LowerCI <- center_effectSizeLight$Estimate - 1.96 * center_effectSizeLight$Std.Error
center_effectSizeLight$UpperCI <- center_effectSizeLight$Estimate + 1.96 * center_effectSizeLight$Std.Error

f1.3b <- ggplot(center_effectSizeLight, aes(x = lightWavelength, y = Estimate, ymin = LowerCI, ymax = UpperCI, color = lightWavelength)) +
  geom_point(size = 0.5) +
  geom_linerange(aes(ymin = LowerCI, ymax = UpperCI), linewidth = 0.25) +
  geom_point(data = subset(center_effectSizeLight, lightWavelength == "470nm"), size = 1) +  
  # geom_linerange(data = subset(center_effectSizeLight, lightWavelength == "470nm"), ymin = subset(center_effectSizeLight, lightWavelength == "470nm")$LowerCI, ymax = subset(center_effectSizeLight, lightWavelength == "470nm")$UpperCI, linewidth = 0.5) +  # 加粗470nm的线
  # geom_text(data = subset(center_effectSizeLight, lightWavelength == "470nm"), aes(label = "470nm"), vjust = 2, size = 3.5, nudge_y = -2) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  scale_color_manual(values = wavelength_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "light wavelength effect size", x = "light wavelength", y = "Partial for light wavelength") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

f1.3b

# a and b to figure 1
f1.3 <- f1.3a + inset_element(f1.3b, 0.05, 0.6, 0.4, 0.95)
legend1.3 <- get_legend(f1.3a + theme(legend.position="bottom",
                                      legend.spacing.x = unit(0.1, "cm"),
                                      legend.key.size = unit(0.5, "cm"),
                                      legend.text = element_text(size = 10)))
Fig1.3 <- cowplot::plot_grid(f1.3, legend1.3, ncol=1, rel_heights = c(1, 0.2))
Fig1.3

## -Distance Fig. 1.4- ----

# Note the replacement of your_data and your_model with your actual data and model names
unique_wavelengths <- unique(individual_distance$lightWavelength)    # replace data here
pred_list <- lapply(unique_wavelengths, generate_prediction_data, data = individual_distance, model = individual_distance.model1)  # replace data and model here
pred_data_all <- do.call(rbind, pred_list)
pred_data_all$CI_low <- pred_data_all$fit - 1.96 * pred_data_all$se.fit
pred_data_all$CI_high <- pred_data_all$fit + 1.96 * pred_data_all$se.fit
pred_data_all$lightWavelength <- factor(pred_data_all$lightWavelength, levels = unique_wavelengths)

# Plotting graphs, including raw data points and predicted confidence intervals
f1.4a <- ggplot() +
  geom_point(data = individual_distance, aes(x = day, y = Value, group = lightWavelength, fill = lightWavelength, color = lightWavelength), size = 1, alpha = 0.1) +   #  replace data here
  geom_ribbon(data = pred_data_all, aes(x = day, ymin = CI_low, ymax = CI_high, fill = lightWavelength), alpha = 0.2, show.legend = FALSE) +
  geom_line(data = pred_data_all, aes(x = day, y = fit, color = lightWavelength)) +
  scale_color_manual(values = wavelength_colors) +
  scale_fill_manual(values = wavelength_colors) +
  labs(title = "Trend of Value Over Time by Light Wavelength",
       x = "Day", y = "Total movement distance",  # replace new name for y
       color = "Light Wavelength", fill = "Light Wavelength") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.5, colour = "black")
  ) +
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE), fill = guide_legend(nrow = 1, byrow = TRUE))

# -1.4b-: light wavelength effect size
distance_effectSizeLight <- data.frame(
  lightWavelength = c("dark", "white", "365nm", "420nm", "470nm", "500nm", "520nm", "590nm", "620nm", "635nm", "660nm"),
  Estimate = c(0, -1781.93, -1671.53, -1847.53, -1706.56, 132.02, 86.93, -784.32, -357.85, 58.78, 63.65),
  Std.Error = rep(37.10, 11)
)

# calculate CI
distance_effectSizeLight$LowerCI <- distance_effectSizeLight$Estimate - 1.96 * distance_effectSizeLight$Std.Error
distance_effectSizeLight$UpperCI <- distance_effectSizeLight$Estimate + 1.96 * distance_effectSizeLight$Std.Error

f1.4b <- ggplot(distance_effectSizeLight, aes(x = lightWavelength, y = Estimate, ymin = LowerCI, ymax = UpperCI, color = lightWavelength)) +
  geom_point(size = 0.5) +
  geom_linerange(aes(ymin = LowerCI, ymax = UpperCI), linewidth = 0.25) +
  geom_point(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), size = 1) +  
  geom_linerange(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), ymin = subset(distance_effectSizeLight, lightWavelength == "470nm")$LowerCI, ymax = subset(distance_effectSizeLight, lightWavelength == "470nm")$UpperCI, linewidth = 0.5) +  # 加粗470nm的线
  geom_text(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), aes(label = "470nm"), vjust = 2, size = 3.5, nudge_y = -2) +  # 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # 
  scale_color_manual(values = wavelength_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "light wavelength effect size", x = "light wavelength", y = "Partial for light wavelength") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

f1.4b

# a and b to figure 1
f1.4 <- f1.4a + inset_element(f1.4b, 0.05, 0.05, 0.4, 0.45)
legend1.4 <- get_legend(f1.4a + theme(legend.position="bottom",
                                      legend.spacing.x = unit(0.1, "cm"),
                                      legend.key.size = unit(0.5, "cm"),
                                      legend.text = element_text(size = 10)))
Fig1.4 <- cowplot::plot_grid(f1.4, legend1.4, ncol=1, rel_heights = c(1, 0.2))
Fig1.4

## save all figures
ggsave("Fig1.1.png", plot = Fig1.1, device = "png", width = 12, height = 9)
ggsave("Fig1.2.png", plot = Fig1.2, device = "png", width = 12, height = 9)
ggsave("Fig1.3.png", plot = Fig1.3, device = "png", width = 12, height = 9)
ggsave("Fig1.4.png", plot = Fig1.4, device = "png", width = 12, height = 9)

## - Part 2: adult group - 

## -Thigmotaxis Fig. 2.1- ----

# -2.1a- Time and trend

# Note the replacement of your_data and your_model with your actual data and model names
unique_wavelengths <- unique(group_thigmotaxis$lightWavelength)    # replace data here
pred_list <- lapply(unique_wavelengths, generate_prediction_data_gamm, data = group_thigmotaxis, gamm_object = group_thigmotaxis.model1)  # replace fuction and data and model here
pred_data_all <- do.call(rbind, pred_list)
pred_data_all$CI_low <- pred_data_all$fit - 1.96 * pred_data_all$se.fit
pred_data_all$CI_high <- pred_data_all$fit + 1.96 * pred_data_all$se.fit
pred_data_all$lightWavelength <- factor(pred_data_all$lightWavelength, levels = unique_wavelengths)

# Plotting graphs, including raw data points and predicted confidence intervals
f2.1a <- ggplot() +
  geom_point(data = group_thigmotaxis, aes(x = day, y = value, group = lightWavelength, fill = lightWavelength, color = lightWavelength), size = 1, alpha = 0.1) +   #  replace data here
  geom_ribbon(data = pred_data_all, aes(x = day, ymin = CI_low, ymax = CI_high, fill = lightWavelength), alpha = 0.2, show.legend = FALSE) +
  geom_line(data = pred_data_all, aes(x = day, y = fit, color = lightWavelength)) +
  scale_color_manual(values = wavelength_colors) +
  scale_fill_manual(values = wavelength_colors) +
  labs(title = "Trend of Value Over Time by Light Wavelength",
       x = "Day", y = "Time in thigmotaxis zone",  # replace new name for y
       color = "Light Wavelength", fill = "Light Wavelength") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.5, colour = "black")
  ) +
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE), fill = guide_legend(nrow = 1, byrow = TRUE))


# - f2.1b -
# subtract effect size and CI
fixed_effects <- summary(group_thigmotaxis.model1$lme)$tTable   # replace the model here

group_thigmotaxis_effectSizeLight <- data.frame(
  lightWavelength = c("dark", "white", "365nm", "420nm", "470nm", "500nm", "520nm", "590nm", "620nm", "635nm", "660nm"),
  Estimate = c(0, 60.747, 99.453533, 115.596253, 146.726933, 48.188733, 50.080267, 59.313133, 62.6292, 58.595, 72.0168),
  Std.Error = c(2.002732, 2.832291, 2.832291, 2.832291, 2.832291, 2.832291, 2.832291, 2.832291, 2.832291, 2.832291, 2.832291)
)

# calculate CI
group_thigmotaxis_effectSizeLight$LowerCI <- group_thigmotaxis_effectSizeLight$Estimate - 1.96 * group_thigmotaxis_effectSizeLight$Std.Error
group_thigmotaxis_effectSizeLight$UpperCI <- group_thigmotaxis_effectSizeLight$Estimate + 1.96 * group_thigmotaxis_effectSizeLight$Std.Error

f2.1b <- ggplot(group_thigmotaxis_effectSizeLight, aes(x = lightWavelength, y = Estimate, ymin = LowerCI, ymax = UpperCI, color = lightWavelength)) +
  geom_point(size = 0.5) +
  geom_linerange(aes(ymin = LowerCI, ymax = UpperCI), linewidth = 0.25) +
  # geom_point(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), size = 1) +  #
  # geom_linerange(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), ymin = subset(distance_effectSizeLight, lightWavelength == "470nm")$LowerCI, ymax = subset(distance_effectSizeLight, lightWavelength == "470nm")$UpperCI, linewidth = 0.5) +  # 加粗470nm的线
  # geom_text(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), aes(label = "470nm"), vjust = 2, size = 3.5, nudge_y = -2) +  # 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # 
  scale_color_manual(values = wavelength_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "light wavelength effect size", x = "light wavelength", y = "Partial for light wavelength") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

f2.1b

# a and b to figure 1
f2.1 <- f2.1a + inset_element(f2.1b, 0.05, 0.6, 0.4, 0.95)
legend2.1 <- get_legend(f2.1a + theme(legend.position="bottom",
                                      legend.spacing.x = unit(0.1, "cm"),
                                      legend.key.size = unit(0.5, "cm"),
                                      legend.text = element_text(size = 10)))
Fig2.1 <- cowplot::plot_grid(f2.1, legend2.1, ncol=1, rel_heights = c(1, 0.2))
Fig2.1

## -Transition Fig. 2.2- ----

# -2.2a- Time and trend

# Note the replacement of your_data and your_model with your actual data and model names
unique_wavelengths <- unique(group_transition$lightWavelength)    # replace data here
pred_list <- lapply(unique_wavelengths, generate_prediction_data_gamm, data = group_transition, gamm_object = group_transition.model1)  # replace data and model here
pred_data_all <- do.call(rbind, pred_list)
pred_data_all$CI_low <- pred_data_all$fit - 1.96 * pred_data_all$se.fit
pred_data_all$CI_high <- pred_data_all$fit + 1.96 * pred_data_all$se.fit
pred_data_all$lightWavelength <- factor(pred_data_all$lightWavelength, levels = unique_wavelengths)

# Plotting graphs, including raw data points and predicted confidence intervals
f2.2a <- ggplot() +
  geom_point(data = group_transition, aes(x = day, y = value, group = lightWavelength, fill = lightWavelength, color = lightWavelength), size = 1, alpha = 0.1) +   #  replace data here
  geom_ribbon(data = pred_data_all, aes(x = day, ymin = CI_low, ymax = CI_high, fill = lightWavelength), alpha = 0.2, show.legend = FALSE) +
  geom_line(data = pred_data_all, aes(x = day, y = fit, color = lightWavelength)) +
  scale_color_manual(values = wavelength_colors) +
  scale_fill_manual(values = wavelength_colors) +
  labs(title = "Trend of Value Over Time by Light Wavelength",
       x = "Day", y = "Time in transition zone",  # replace new name for y
       color = "Light Wavelength", fill = "Light Wavelength") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.5, colour = "black")
  ) +
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE), fill = guide_legend(nrow = 1, byrow = TRUE))

f2.2a

# - f2.2b -

fixed_effects <- summary(group_transition.model1$lme)$tTable   # replace the model here

group_transition_effectSizeLight <- data.frame(
  lightWavelength = c("dark", "white", "365nm", "420nm", "470nm", "500nm", "520nm", "590nm", "620nm", "635nm", "660nm"),
  Estimate = c(0, -28.290933, -58.909867, -53.406467, -77.312867, -21.609400, -17.545667, -37.324333, -28.759867, -38.625000, -35.069667),
  Std.Error = c(1.391870,  1.968402,  1.968402,  1.968402,  1.968402,  1.968402,  1.968402,  1.968402,  1.968402,  1.968402,  1.968402)
)

# calculate CI
group_transition_effectSizeLight$LowerCI <- group_transition_effectSizeLight$Estimate - 1.96 * group_transition_effectSizeLight$Std.Error
group_transition_effectSizeLight$UpperCI <- group_transition_effectSizeLight$Estimate + 1.96 * group_transition_effectSizeLight$Std.Error

f2.2b <- ggplot(group_transition_effectSizeLight, aes(x = lightWavelength, y = Estimate, ymin = LowerCI, ymax = UpperCI, color = lightWavelength)) +
  geom_point(size = 0.5) +
  geom_linerange(aes(ymin = LowerCI, ymax = UpperCI), linewidth = 0.25) +
  # geom_point(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), size = 1) +  # 
  # geom_linerange(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), ymin = subset(distance_effectSizeLight, lightWavelength == "470nm")$LowerCI, ymax = subset(distance_effectSizeLight, lightWavelength == "470nm")$UpperCI, linewidth = 0.5) +  # 加粗470nm的线
  # geom_text(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), aes(label = "470nm"), vjust = 2, size = 3.5, nudge_y = -2) +  # 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # 
  scale_color_manual(values = wavelength_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "light wavelength effect size", x = "light wavelength", y = "Partial for light wavelength") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

f2.2b

# a and b to figure 1
f2.2 <- f2.2a + inset_element(f2.2b, 0.05, 0.05, 0.4, 0.4)
legend2.2 <- get_legend(f2.2a + theme(legend.position="bottom",
                                      legend.spacing.x = unit(0.1, "cm"),
                                      legend.key.size = unit(0.5, "cm"),
                                      legend.text = element_text(size = 10)))
Fig2.2 <- cowplot::plot_grid(f2.2, legend2.2, ncol=1, rel_heights = c(1, 0.2))
Fig2.2

## -Center Fig. 2.3- ----

# -2.3a- Time and trend 

# Note the replacement of your_data and your_model with your actual data and model names
unique_wavelengths <- unique(group_center$lightWavelength)    # replace data here
pred_list <- lapply(unique_wavelengths, generate_prediction_data_gamm, data = group_center, gamm_object = group_center.model1)  # replace data and model here
pred_data_all <- do.call(rbind, pred_list)
pred_data_all$CI_low <- pred_data_all$fit - 1.96 * pred_data_all$se.fit
pred_data_all$CI_high <- pred_data_all$fit + 1.96 * pred_data_all$se.fit
pred_data_all$lightWavelength <- factor(pred_data_all$lightWavelength, levels = unique_wavelengths)

# Plotting graphs, including raw data points and predicted confidence intervals
f2.3a <- ggplot() +
  geom_point(data = group_center, aes(x = day, y = value, group = lightWavelength, fill = lightWavelength, color = lightWavelength), size = 1, alpha = 0.1) +   #  replace data here
  geom_ribbon(data = pred_data_all, aes(x = day, ymin = CI_low, ymax = CI_high, fill = lightWavelength), alpha = 0.2, show.legend = FALSE) +
  geom_line(data = pred_data_all, aes(x = day, y = fit, color = lightWavelength)) +
  scale_color_manual(values = wavelength_colors) +
  scale_fill_manual(values = wavelength_colors) +
  labs(title = "Trend of Value Over Time by Light Wavelength",
       x = "Day", y = "Time in center zone",  # replace new name for y
       color = "Light Wavelength", fill = "Light Wavelength") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.5, colour = "black")
  ) +
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE), fill = guide_legend(nrow = 1, byrow = TRUE))

f2.3a

# - f2.2b -

fixed_effects <- summary(group_center.model1$lme)$tTable   # replace the model here

effect_size_df <- data.frame(
  lightWavelength = rownames(fixed_effects),
  Estimate = fixed_effects[, "Value"],  
  Std.Error = fixed_effects[, "Std.Error"]
)
effect_size_df$Estimate
effect_size_df$Std.Error

group_center_effectSizeLight <- data.frame(
  lightWavelength = c("dark", "white", "365nm", "420nm", "470nm", "500nm", "520nm", "590nm", "620nm", "635nm", "660nm"),
  Estimate = c(0, -32.4560667, -40.5436667, -62.1898000, -69.4140667, -26.5793333, -32.5346000, -21.9888000, -33.8693333, -19.9700000, -36.9471333),
  Std.Error = c(2.575697,  3.642586,  3.642586,  3.642586,  3.642586,  3.642586,  3.642586,  3.642586,  3.642586,  3.642586,  3.642586)
)

# calculate CI
group_center_effectSizeLight$LowerCI <- group_center_effectSizeLight$Estimate - 1.96 * group_center_effectSizeLight$Std.Error
group_center_effectSizeLight$UpperCI <- group_center_effectSizeLight$Estimate + 1.96 * group_center_effectSizeLight$Std.Error

f2.3b <- ggplot(group_center_effectSizeLight, aes(x = lightWavelength, y = Estimate, ymin = LowerCI, ymax = UpperCI, color = lightWavelength)) + # replave data here
  geom_point(size = 0.5) +
  geom_linerange(aes(ymin = LowerCI, ymax = UpperCI), linewidth = 0.25) +
  # geom_point(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), size = 1) +  #
  # geom_linerange(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), ymin = subset(distance_effectSizeLight, lightWavelength == "470nm")$LowerCI, ymax = subset(distance_effectSizeLight, lightWavelength == "470nm")$UpperCI, linewidth = 0.5) +  # 加粗470nm的线
  # geom_text(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), aes(label = "470nm"), vjust = 2, size = 3.5, nudge_y = -2) +  #
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # 
  scale_color_manual(values = wavelength_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "light wavelength effect size", x = "light wavelength", y = "Partial for light wavelength") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

f2.3b

# a and b to figure 1
f2.3 <- f2.3a + inset_element(f2.3b, 0.05, 0.05, 0.4, 0.4)
legend2.3 <- get_legend(f2.3a + theme(legend.position="bottom",
                                      legend.spacing.x = unit(0.1, "cm"),
                                      legend.key.size = unit(0.5, "cm"),
                                      legend.text = element_text(size = 10)))
Fig2.3 <- cowplot::plot_grid(f2.3, legend2.3, ncol=1, rel_heights = c(1, 0.2))
Fig2.3

## -Distance Fig. 2.4- ----  

# -2.4a- Time and trend 

# Note the replacement of your_data and your_model with your actual data and model names
unique_wavelengths <- unique(group_distance$lightWavelength)    # replace data here
pred_list <- lapply(unique_wavelengths, generate_prediction_data_gamm, data = group_distance, gamm_object = group_distance.model1)  # replace data and model here
pred_data_all <- do.call(rbind, pred_list)
pred_data_all$CI_low <- pred_data_all$fit - 1.96 * pred_data_all$se.fit
pred_data_all$CI_high <- pred_data_all$fit + 1.96 * pred_data_all$se.fit
pred_data_all$lightWavelength <- factor(pred_data_all$lightWavelength, levels = unique_wavelengths)

# Plotting graphs, including raw data points and predicted confidence intervals
f2.4a <- ggplot() +
  geom_point(data = group_distance, aes(x = day, y = value, group = lightWavelength, fill = lightWavelength, color = lightWavelength), size = 1, alpha = 0.1) +   #  replace data here
  geom_ribbon(data = pred_data_all, aes(x = day, ymin = CI_low, ymax = CI_high, fill = lightWavelength), alpha = 0.2, show.legend = FALSE) +
  geom_line(data = pred_data_all, aes(x = day, y = fit, color = lightWavelength)) +
  scale_color_manual(values = wavelength_colors) +
  scale_fill_manual(values = wavelength_colors) +
  labs(title = "Trend of Value Over Time by Light Wavelength",
       x = "Day", y = "Total movement distance",  # replace new name for y
       color = "Light Wavelength", fill = "Light Wavelength") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.5, colour = "black")
  ) +
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE), fill = guide_legend(nrow = 1, byrow = TRUE))

f2.4a

# - f2.2b -

fixed_effects <- summary(group_distance.model1$lme)$tTable   # replace the model here

effect_size_df <- data.frame(
  lightWavelength = rownames(fixed_effects),
  Estimate = fixed_effects[, "Value"],  # 
  Std.Error = fixed_effects[, "Std.Error"]
)
effect_size_df$Estimate
effect_size_df$Std.Error

group_distance_effectSizeLight <- data.frame(
  lightWavelength = c("dark", "white", "365nm", "420nm", "470nm", "500nm", "520nm", "590nm", "620nm", "635nm", "660nm"),
  Estimate = c(0, -356.78547,  -544.14460,    54.14933, -1070.08953,  1499.01200,   458.45720,    11.96307, -1109.62987,   -16.74113,   -11.87000),
  Std.Error = c(66.23468,  93.66998,  93.66998,  93.66998,  93.66998,  93.66998,  93.66998,  93.66998,  93.66998,  93.66998,  93.66998)
)

# calculate CI
group_distance_effectSizeLight$LowerCI <- group_distance_effectSizeLight$Estimate - 1.96 * group_distance_effectSizeLight$Std.Error
group_distance_effectSizeLight$UpperCI <- group_distance_effectSizeLight$Estimate + 1.96 * group_distance_effectSizeLight$Std.Error

f2.4b <- ggplot(group_distance_effectSizeLight, aes(x = lightWavelength, y = Estimate, ymin = LowerCI, ymax = UpperCI, color = lightWavelength)) + # replave data here
  geom_point(size = 0.5) +
  geom_linerange(aes(ymin = LowerCI, ymax = UpperCI), linewidth = 0.25) +
  # geom_point(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), size = 1) +  # 
  # geom_linerange(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), ymin = subset(distance_effectSizeLight, lightWavelength == "470nm")$LowerCI, ymax = subset(distance_effectSizeLight, lightWavelength == "470nm")$UpperCI, linewidth = 0.5) +  # 加粗470nm的线
  # geom_text(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), aes(label = "470nm"), vjust = 2, size = 3.5, nudge_y = -2) +  #
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # 
  scale_color_manual(values = wavelength_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "light wavelength effect size", x = "light wavelength", y = "Partial for light wavelength") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

f2.4b

# a and b to figure 1
f2.4 <- f2.4a + inset_element(f2.4b, 0.05, 0.64, 0.4, 0.99)
legend2.4 <- get_legend(f2.4a + theme(legend.position="bottom",
                                      legend.spacing.x = unit(0.1, "cm"),
                                      legend.key.size = unit(0.5, "cm"),
                                      legend.text = element_text(size = 10)))
Fig2.4 <- cowplot::plot_grid(f2.4, legend2.4, ncol=1, rel_heights = c(1, 0.2))
Fig2.4

## -Inter_individual_distance Fig. 2.4- ----  

# -2.5a- Time and trend 

# Note the replacement of your_data and your_model with your actual data and model names
unique_wavelengths <- unique(group_inter_individual_distance$lightWavelength)    # replace data here
pred_list <- lapply(unique_wavelengths, generate_prediction_data_gamm, data = group_inter_individual_distance, gamm_object = group_inter_individual_distance.model1)  # replace data and model here
pred_data_all <- do.call(rbind, pred_list)
pred_data_all$CI_low <- pred_data_all$fit - 1.96 * pred_data_all$se.fit
pred_data_all$CI_high <- pred_data_all$fit + 1.96 * pred_data_all$se.fit
pred_data_all$lightWavelength <- factor(pred_data_all$lightWavelength, levels = unique_wavelengths)

# Plotting graphs, including raw data points and predicted confidence intervals
f2.5a <- ggplot() +
  geom_point(data = group_inter_individual_distance, aes(x = day, y = value, group = lightWavelength, fill = lightWavelength, color = lightWavelength), size = 1, alpha = 0.1) +   #  replace data here
  geom_ribbon(data = pred_data_all, aes(x = day, ymin = CI_low, ymax = CI_high, fill = lightWavelength), alpha = 0.2, show.legend = FALSE) +
  geom_line(data = pred_data_all, aes(x = day, y = fit, color = lightWavelength)) +
  scale_color_manual(values = wavelength_colors) +
  scale_fill_manual(values = wavelength_colors) +
  labs(title = "Trend of Value Over Time by Light Wavelength",
       x = "Day", y = "Inter individual distance",  # replace new name for y
       color = "Light Wavelength", fill = "Light Wavelength") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.5, colour = "black")
  ) +
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE), fill = guide_legend(nrow = 1, byrow = TRUE))

f2.5a

# - f2.5b -

fixed_effects <- summary(group_inter_individual_distance.model1$lme)$tTable   # replace the model here

effect_size_df <- data.frame(
  lightWavelength = rownames(fixed_effects),
  Estimate = fixed_effects[, "Value"],  # 
  Std.Error = fixed_effects[, "Std.Error"]
)
effect_size_df$Estimate
effect_size_df$Std.Error

group_inter_individual_distance_effectSizeLight <- data.frame(
  lightWavelength = c("dark", "white", "365nm", "420nm", "470nm", "500nm", "520nm", "590nm", "620nm", "635nm", "660nm"),
  Estimate = c(0, -5.7437333, -5.8491333, -7.0866667, -6.4932667, -7.7922667, -4.4887333, -7.0268000, -6.3095333, -5.3044667, -6.5193333),
  Std.Error = c(0.4525893, 0.6400580, 0.6400580, 0.6400580, 0.6400580, 0.6400580, 0.6400580, 0.6400580, 0.6400580, 0.6400580, 0.6400580)
)

# calculate CI
group_inter_individual_distance_effectSizeLight$LowerCI <- group_inter_individual_distance_effectSizeLight$Estimate - 1.96 * group_inter_individual_distance_effectSizeLight$Std.Error
group_inter_individual_distance_effectSizeLight$UpperCI <- group_inter_individual_distance_effectSizeLight$Estimate + 1.96 * group_inter_individual_distance_effectSizeLight$Std.Error

f2.5b <- ggplot(group_inter_individual_distance_effectSizeLight, aes(x = lightWavelength, y = Estimate, ymin = LowerCI, ymax = UpperCI, color = lightWavelength)) + # replave data here
  geom_point(size = 0.5) +
  geom_linerange(aes(ymin = LowerCI, ymax = UpperCI), linewidth = 0.25) +
  # geom_point(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), size = 1) +  # 
  # geom_linerange(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), ymin = subset(distance_effectSizeLight, lightWavelength == "470nm")$LowerCI, ymax = subset(distance_effectSizeLight, lightWavelength == "470nm")$UpperCI, linewidth = 0.5) +  # 加粗470nm的线
  # geom_text(data = subset(distance_effectSizeLight, lightWavelength == "470nm"), aes(label = "470nm"), vjust = 2, size = 3.5, nudge_y = -2) +  # 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # 
  scale_color_manual(values = wavelength_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "light wavelength effect size", x = "light wavelength", y = "Partial for light wavelength") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

f2.5b

# a and b to figure 1
f2.5 <- f2.5a + inset_element(f2.5b, 0.6, 0.6, 0.95, 0.95)
legend2.5 <- get_legend(f2.5a + theme(legend.position="bottom",
                                      legend.spacing.x = unit(0.1, "cm"),
                                      legend.key.size = unit(0.5, "cm"),
                                      legend.text = element_text(size = 10)))
Fig2.5 <- cowplot::plot_grid(f2.5, legend2.5, ncol=1, rel_heights = c(1, 0.2))
Fig2.5

## save figures
ggsave("Fig2.1.png", plot = Fig2.1, device = "png", width = 12, height = 9)
ggsave("Fig2.2.png", plot = Fig2.2, device = "png", width = 12, height = 9)
ggsave("Fig2.3.png", plot = Fig2.3, device = "png", width = 12, height = 9)
ggsave("Fig2.4.png", plot = Fig2.4, device = "png", width = 12, height = 9)
ggsave("Fig2.5.png", plot = Fig2.5, device = "png", width = 12, height = 9)



## - Part 3: offspring -

# get mean and se
ALAN_summary_distance <- ALAN_offspring %>%
  group_by(lightWavelength) %>%
  summarise(mean_distance = mean(totalMovementDistance),
            se = sd(totalMovementDistance) / sqrt(n()))  

ALAN_summary_activity <- ALAN_offspring %>%
  group_by(lightWavelength) %>%
  summarise(mean_activity = mean(activity),
            se = sd(activity) / sqrt(n())) 


f3.1 <- ggplot(ALAN_summary_distance, aes(x = lightWavelength, y = mean_distance, fill = lightWavelength)) +
  # geom_point(data = ALAN_offspring, aes(x = lightWavelength, y = totalMovementDistance, color = lightWavelength), position = position_jitter(0.2), size = 2, shape = 20) + 
  geom_bar(stat="identity", color = "black", width = 0.9) +  
  geom_errorbar(aes(ymin = mean_distance - se, ymax = mean_distance + se), width = 0.2) + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size = 12)) +
  scale_fill_manual(values = wavelength_colors) +  
  labs(x = "light", y = "movement distance (cm)") + 
  guides(fill = guide_legend("light")) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,250))

print(f3.1)

f3.2 <- ggplot(ALAN_summary_activity, aes(x = lightWavelength, y = mean_activity, fill = lightWavelength)) +
  # geom_point(data = ALAN_offspring, aes(x = lightWavelength, y = activity, color = lightWavelength), position = position_jitter(0.2), size = 2, shape = 20) + 
  geom_bar(stat="identity", color = "black", width = 0.9) +  
  geom_errorbar(aes(ymin = mean_activity - se, ymax = mean_activity + se), width = 0.2) +  
  theme_classic() + 
  theme(axis.text.x = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size = 12)) +
  scale_fill_manual(values = wavelength_colors) + 
  labs(x = "light", y = "activity (s)") + 
  guides(fill = guide_legend("light")) +  
  scale_y_continuous(expand = c(0,0), limits = c(0,40))

print(f3.2)

ggsave("Fig3.1.png", plot = f3.1, device = "png", width = 12, height = 9)
ggsave("Fig3.2.png", plot = f3.2, device = "png", width = 12, height = 9)




###### MAXX Effect size calculations #####

# get mean and se and n for each wavelength, behavior type, and day
df_sum <- ALAN_individual %>%
  group_by(lightWavelength, behaviorType, day) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    n = n_distinct(individualAdultID),
    .groups = "drop"
  )

# separate control from ALAN wavelengths
control_df <- df_sum %>%
  filter(lightWavelength == "dark") %>%
  rename(
    mean_control = mean,
    sd_control = sd,
    n_control = n
  )

alan_df <- df_sum %>%
  filter(lightWavelength != "dark") %>%
  rename(
    wavelength = lightWavelength,
    mean_treatment = mean,
    sd_treatment = sd,
    n_treatment = n
  )

# join wavelength/day to match control/day
es_df <- alan_df %>%
  left_join(control_df, by = c("behaviorType", "day")) %>%
  mutate(
    yi = log(mean_treatment / mean_control),
    vi = (sd_treatment^2 / (n_treatment * mean_treatment^2)) +
      (sd_control^2 / (n_control * mean_control^2)),
    Study_ID = "Zebrafish_wavelength_study"
  )


# filter out unncessary variables
es_df <- es_df %>%
  filter(behaviorType != "Time_In_Transition_Zone" &
           behaviorType != "Time_In_Center_Zone")

library(readr)

write_csv(es_df, "es_df.csv")



log(4.4022/3.954)


(0.6888^2 / (90*4.4022^2)) + (0.3218370115^2/ (85*3.954022989))
