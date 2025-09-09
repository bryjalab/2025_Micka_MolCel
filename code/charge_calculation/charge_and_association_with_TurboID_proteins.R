library("ggplot2")

### Define sequences charge proximal to DEP
linkerWT <- "MSPSLSTITSTSSSITSSIPDTER"
plinkerWT <- "MXPXLXZIZXZXXXIZXXIPDZER"

linkerdel <-  "AWVSHTAAMTGTFPAYGMSPDTER"
plinkerdel <- "AWVXHZAAMZGZFPAYGMXPDZER"

linkerA <- "MSPALAAIAAAAAAIAAAIPDTER"
plinkerA <- "MXPALAAIAAAAAAIAAAIPDZER"

linkerE <- "MSPELEEIEEEEEEIEEEIPDTER"
plinkerE <- "MXPELEEIEEEEEEIEEEIPDZER"


#### calculate charge
chargeWT_7.2 <- netCharge(linkerWT, pH = 7.2, pKaSet = custom_pka)
chargepWT_7.2 <- netCharge(plinkerWT, pH = 7.2, pKaSet = custom_pka)

chargedel_7.2 <- netCharge(linkerdel, pH = 7.2, pKaSet = custom_pka)
chargepdel_7.2 <- netCharge(plinkerdel, pH = 7.2, pKaSet = custom_pka)

chargeA_7.2 <- netCharge(linkerA, pH = 7.2, pKaSet = custom_pka)
chargepA_7.2 <- netCharge(plinkerA, pH = 7.2, pKaSet = custom_pka)

chargeE_7.2 <- netCharge(linkerE, pH = 7.2, pKaSet = custom_pka)
chargepE_7.2 <- netCharge(plinkerE, pH = 7.2, pKaSet = custom_pka)


### FZDs
## FZD3 association Based on Charge
# Dataframe with charge and association
FZD3 <- data.frame(charge = c(chargepA_7.2, chargepdel_7.2, chargepE_7.2, chargepWT_7.2)*(-1),
                  FZD3 = c(10.6426098636964, 10.4119552020055, 5.93149393418557, 0),
                  color = c("#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))

# Fit a sigmoidal curve (logistic function) using nls()
fit_FZD3 <- nls(FZD3 ~ A / (1 + exp(k * (charge - x0))),
               data = FZD3,
               start = list(A = 10, k = 0.6, x0 = 13))

# Print the coefficients and check the fit
x0_value_FZD3 <- coef(fit_FZD3)["x0"]
summary(fit_FZD3)

# Get fitted values from the model
FZD3$fitted <- predict(fit_FZD3)

# Create a fine grid of charge values for a smoother curve
charge_grid_FZD3 <- seq(0, 30, length.out = 500)

# Predict the activity for each charge in the fine grid
predicted_activity_FZD3 <- predict(fit_FZD3, newdata = data.frame(charge = charge_grid_FZD3))

# Create a new dataframe for plotting the smooth curve
smooth_data_FZD3 <- data.frame(charge = charge_grid_FZD3, FZD3 = predicted_activity_FZD3)

# Plot original data and fitted sigmoidal curve
svg(filename = "FZD3_charge.svg", width = 3, height = 2)
ggplot(FZD3, aes(x = charge, y = FZD3)) +
  geom_line(data = smooth_data_FZD3, aes(x = charge, y = FZD3), color = "grey50", linewidth = 1)+
  geom_point(size = 4.5, fill = FZD3$color,  shape = 21) +  # Plot the original data points
  geom_vline(xintercept=x0_value_FZD3, linetype = "dotted", color = "black", linewidth = 1)+
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "30"))+
  coord_cartesian(ylim = c(-1, max(FZD3$FZD3) + 1))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)

dev.off()


## FZD5 association Based on Charge
# Dataframe with charge and association
FZD5 <- data.frame(charge = c(chargepA_7.2, chargepdel_7.2, chargepE_7.2, chargepWT_7.2)*(-1),
                   FZD5 = c(13.1139026641177, 10.0057236571395, 0.337684349957501, 0),
                   color = c("#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))

# Fit a sigmoidal curve (logistic function) using nls()
fit_FZD5 <- nls(FZD5 ~ A / (1 + exp(k * (charge - x0))),
                data = FZD5,
                start = list(A = 10, k = 0.6, x0 = 13))

# Print the coefficients and check the fit
x0_value_FZD5 <- coef(fit_FZD5)["x0"]
summary(fit_FZD5)

# Get fitted values from the model
FZD5$fitted <- predict(fit_FZD5)

# Create a fine grid of charge values for a smoother curve
charge_grid_FZD5 <- seq(0, 30, length.out = 500)

# Predict the activity for each charge in the fine grid
predicted_activity_FZD5 <- predict(fit_FZD5, newdata = data.frame(charge = charge_grid_FZD5))

# Create a new dataframe for plotting the smooth curve
smooth_data_FZD5 <- data.frame(charge = charge_grid_FZD5, FZD5 = predicted_activity_FZD5)

# Plot original data and fitted sigmoidal curve
svg(filename = "FZD5_charge.svg", width = 3, height = 2)
ggplot(FZD5, aes(x = charge, y = FZD5)) +
  geom_line(data = smooth_data_FZD5, aes(x = charge, y = FZD5), color = "grey50", linewidth = 1)+
  geom_point(size = 4.5, fill = FZD5$color,  shape = 21) +  # Plot the original data points
  geom_vline(xintercept=x0_value_FZD5, linetype = "dotted", color = "black", linewidth = 1)+
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25,30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-1, max(FZD5$FZD5) + 1))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()

## FZD6 association Based on Charge
# Dataframe with charge and association
FZD6 <- data.frame(charge = c(chargepA_7.2, chargepdel_7.2, chargepE_7.2, chargepWT_7.2)*(-1),
                   FZD6 = c(8.27573911415643, 8.13342458876843,6.40770743012332, 0),
                   color = c("#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))

# Fit a sigmoidal curve (logistic function) using nls()
fit_FZD6 <- nls(FZD6 ~ A / (1 + exp(k * (charge - x0))),
                data = FZD6,
                start = list(A = 10, k = 0.6, x0 = 13))

# Print the coefficients and check the fit
x0_value_FZD6 <- coef(fit_FZD6)["x0"]
summary(fit_FZD6)

# Get fitted values from the model
FZD6$fitted <- predict(fit_FZD6)

# Create a fine grid of charge values for a smoother curve
charge_grid_FZD6 <- seq(0, 30, length.out = 500)

# Predict the activity for each charge in the fine grid
predicted_activity_FZD6 <- predict(fit_FZD6, newdata = data.frame(charge = charge_grid_FZD6))

# Create a new dataframe for plotting the smooth curve
smooth_data_FZD6 <- data.frame(charge = charge_grid_FZD6, FZD6 = predicted_activity_FZD6)

# Plot original data and fitted sigmoidal curve
svg(filename = "FZD6_charge.svg", width = 3, height = 2)
ggplot(FZD6, aes(x = charge, y = FZD6)) +
  geom_line(data = smooth_data_FZD6, aes(x = charge, y = FZD6), color = "grey50", linewidth = 1)+
  geom_point(size = 4.5, fill = FZD6$color, shape = 21) +  # Plot the original data points
  geom_vline(xintercept=x0_value_FZD6, linetype = "dotted", color = "black", linewidth = 1)+
  labs(x = "Charge", y = "FZD6")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-1, max(FZD6$FZD6) + 1))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()

### Proteins with no enrichment 
## Axin association Based on Charge
# Dataframe with charge and association
Axin <- data.frame(charge = c(chargepA_7.2, chargepdel_7.2, chargepE_7.2, chargepWT_7.2)*(-1),
                   Axin1 = c(2.60211044690386, 3.82221595235317, 0.631077893440392, 0),
                   Axin2 = c(3.07712341376974, 1.20480803148565, 0.753491003788647, 0),
                   color = c("#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))


# Plot original data
svg(filename = "Axin_charge.svg", width = 3, height = 2)
ggplot(Axin, aes(x = charge, y = Axin1)) +
  geom_point(size = 4.5, fill = Axin$color, shape = 22) +
  geom_point(aes(x = charge, y = Axin2), size = 4.5, fill = Axin$color, shape = 24) +
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-1, max(FZD6$FZD6)+1))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()

## DVL association Based on Charge
# Dataframe with charge and association
DVL <- data.frame(charge = c(chargepA_7.2, chargepdel_7.2, chargepE_7.2, chargepWT_7.2)*(-1),
                   DVL1 = c(1.22723172932593, 1.5097069609382, 0.79337343110206, 0),
                   DVL2 = c(4.51254125187909, 7.63590389234427, 3.94449773638389, 0)*(-1),
                   DVL3 = c(1.54834570375174, 4.17334976811018, 0.00767525914558874, 0),
                   color = c("#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))


# Plot original data for DVL1 and DVL2
svg(filename = "DVL_charge.svg", width = 3, height = 2)
ggplot(DVL, aes(x = charge, y =DVL1)) +
  geom_point(size = 4.5, fill = DVL$color, shape = 22) +
  geom_point(aes(x = charge, y = DVL2), size = 4.5, fill = DVL$color, shape = 24) +
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30), 
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-8.6, 8.6))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()

# Plot original data for DVL3
svg(filename = "DVL3_charge.svg", width = 3, height = 2)
ggplot(DVL, aes(x = charge, y =DVL3)) +
  geom_point(size = 4.5, fill = DVL$color, shape = 23) + 
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-8.6, 8.6))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()

## VANGL association Based on Charge
# Dataframe with charge and association
VANGL <- data.frame(charge = c(chargepA_7.2, chargepdel_7.2, chargepE_7.2, chargepWT_7.2)*(-1),
                    VANGL1 = c(0.43262730265152, 0.115704482472179, 2.58043818909244, 0),
                    VANGL2 = c(0.358738430753415, 0.52372314195846, 3.55813158776423, 0),
                    color = c("#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))


# Plot original data
svg(filename = "VANGL_charge.svg", width = 3, height = 2)
ggplot(VANGL, aes(x = charge, y = VANGL1)) +
  geom_point(size = 4.5, fill = VANGL$color, shape = 22) +
  geom_point(aes(x = charge, y = VANGL2), size = 4.5, fill = VANGL$color, shape = 24) +
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30), 
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-8.6, 8.6))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()

## CK1 association Based on Charge
# Dataframe with charge and association
CSNK1 <- data.frame(charge = c(chargepA_7.2, chargepdel_7.2, chargepE_7.2, chargepWT_7.2)*(-1),
                    CSNK1D = c(1.42191912749513, 1.45253010382519, 2.43885121389078, 0)*(-1),
                    CSNK1E = c(1.18143317216923, 3.43934395674114, 3.45689224367542, 0)*(-1),
                    color = c("#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))


# Plot original data
svg(filename = "CK1_charge.svg", width = 3, height = 2)
ggplot(CSNK1, aes(x = charge, y = CSNK1D)) +
  geom_point(size = 4.5, fill = CSNK1$color, shape = 22) +
  geom_point(aes(x = charge, y = CSNK1E), size = 4.5, fill = CSNK1$color, shape = 24) +
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-8.6, 8.6))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()


### Proteins enriched in WT
## TNKS1 association Based on Charge
# Dataframe with charge and association
TNKS1 <- data.frame(charge = c(chargepA_7.2, chargepdel_7.2, chargepE_7.2, chargepWT_7.2)*(-1),
                    TNKS1 = c(4.69050452172476,5.76111889078931, 3.11976412193476, 0)*(-1),
                    color = c("#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))

# Plot original data
svg(filename = "TNKS1_charge.svg", width = 3, height = 2)
ggplot(TNKS1, aes(x = charge, y = TNKS1)) +
  geom_point(size = 4.5, fill = TNKS1$color, shape = 21) +
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(min(TNKS1$TNKS1)-1, 1))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()

### TNKS2 association Based on Charge
# Dataframe with charge and association
TNKS2 <- data.frame(charge = c(chargepA_7.2, chargepdel_7.2, chargepE_7.2, chargepWT_7.2)*(-1),
                    TNKS2 = c(11.0172395032917,11.2485153561147, 5.55637335716639, 0)*(-1),
                    color = c("#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))

# Plot original data
svg(filename = "TNKS2_charge.svg", width = 3, height = 2)
ggplot(TNKS2, aes(x = charge, y = TNKS2)) +
  geom_point(size = 4.5, fill = TNKS2$color, shape = 21) +
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(min(TNKS2$TNKS2)-1, 1))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()


## CCDC88A association Based on Charge
# Dataframe with charge and association
CCDC88A <- data.frame(charge = c(chargepA_7.2, chargepdel_7.2, chargepE_7.2, chargepWT_7.2)*(-1),
                      CCDC88A = c(14.3587780230887,13.3622249774127, 13.6875022394778, 0)*(-1),
                   color = c("#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))

# Plot original data
svg(filename = "CCDC88A_charge.svg", width = 3, height = 2)
ggplot(CCDC88A, aes(x = charge, y = CCDC88A)) +
  geom_point(size = 4.5, fill = CCDC88A$color, shape = 21) +
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(min(CCDC88A$CCDC88A)-1, 1))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()

## CCDC88C association Based on Charge
# Dataframe with charge and association
CCDC88C <- data.frame(charge = c(chargepA_7.2, chargepdel_7.2, chargepE_7.2, chargepWT_7.2)*(-1),
                      CCDC88C = c(16.039826572818,15.077377435508, 12.8378946064521, 0)*(-1),
                      color = c("#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))


# Plot original data
svg(filename = "CCDC88C_charge.svg", width = 3, height = 2)
ggplot(CCDC88C, aes(x = charge, y = CCDC88C)) +
  geom_point(size = 4.5, fill = CCDC88C$color, shape = 21) +
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(min(CCDC88C$CCDC88C)-1, 1))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()


### Proteins enriched in phospho-switch mutants
## ADGRA2 association Based on Charge
# dataframe with charge and association
ADGRA2 <- data.frame(charge = c(chargepA_7.2, chargepdel_7.2, chargepE_7.2, chargepWT_7.2)*(-1),
                     ADGRA2 = c(14.9904125338438,15.011234468656,10.6548122122413, 0),
                     color = c("#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))


# Plot original data
svg(filename = "ADGRA2_charge.svg", width = 3, height = 2)
ggplot(ADGRA2, aes(x = charge, y = ADGRA2)) +
  geom_point(size = 4.5, fill = FZD5$color,  shape = 21) +
   labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-1, max(ADGRA2$ADGRA2) + 1))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()

## RYK association Based on Charge
# dataframe with charge and association
RYK <- data.frame(charge = c(chargepA_7.2, chargepdel_7.2, chargepE_7.2, chargepWT_7.2)*(-1),
                  RYK = c(6.75999337805435,7.87150641248253,3.16870435635016, 0),
                     color = c("#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))


# Plot original data
svg(filename = "RYK_charge.svg", width = 3, height = 2)
ggplot(RYK, aes(x = charge, y = RYK)) +
  geom_point(size = 4.5, fill = FZD5$color,  shape = 21) +
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-1, max(RYK$RYK) + 1))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()


## ZNRF3 association Based on Charge
# dataframe with charge and association
ZNRF3 <- data.frame(charge = c(chargepA_7.2, chargepdel_7.2, chargepE_7.2, chargepWT_7.2)*(-1),
                    ZNRF3 = c(5.52822522341622,4.65693271284799,4.1705105535001, 0),
                     color = c("#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))

# Plot original data
svg(filename = "ZNRF3_charge.svg", width = 3, height = 2)
ggplot(ZNRF3, aes(x = charge, y = ZNRF3)) +
  geom_point(size = 4.5, fill = FZD5$color,  shape = 21) +
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-1, max(ZNRF3$ZNRF3) + 1))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()


## RNF43 association Based on Charge
# Dataframe with charge and association
RNF43 <- data.frame(charge =  c(chargepA_7.2, chargepdel_7.2, chargepE_7.2, chargepWT_7.2)*(-1),
                    RNF43 = c(8.06213011679581, 7.32992480878802,1.12416404651768, 0),
                    color = c("#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))

# Plot original data
svg(filename = "RNF43_charge.svg", width = 3, height = 2)
ggplot(RNF43, aes(x = charge, y = RNF43)) +
  geom_point(size = 4.5, fill = FZD5$color,  shape = 21) +
  labs(x = "Charge", y = "RNF43")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-1, max(RNF43$RNF43) + 1))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()
