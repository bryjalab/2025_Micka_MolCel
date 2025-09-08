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
chargeWT <- netCharge(linkerWT, pH = 6.5, pKaSet = custom_pka)
chargepWT <- netCharge(plinkerWT, pH = 6.5, pKaSet = custom_pka)

chargedel <- netCharge(linkerdel, pH = 6.5, pKaSet = custom_pka)
chargepdel <- netCharge(plinkerdel, pH = 6.5, pKaSet = custom_pka)

chargeA <- netCharge(linkerA, pH = 6.5, pKaSet = custom_pka)
chargepA <- netCharge(plinkerA, pH = 6.5, pKaSet = custom_pka)

chargeE <- netCharge(linkerE, pH = 6.5, pKaSet = custom_pka)
chargepE <- netCharge(plinkerE, pH = 6.5, pKaSet = custom_pka)


### plot change of charge
# dataframe
dcharge <- data.frame(mutant = rep(c("WT", "del", "A", "E"), each = 2),
                      state = rep(c("nonphos", "phos"), times = 4),
                      values = c(chargeWT, chargepWT,
                                 chargedel, chargepdel,
                                 chargeA, chargepA,
                                 chargeE, chargepE)
                      )

dcharge$mutant <- factor(dcharge$mutant, levels = c("WT", "del", "A", "E"))

# arrows dataframe
arrows_df <- data.frame(
  mutant = c("WT", "del", "A", "E"),
  y_start = c(chargeWT, chargedel, chargeA, chargeE),
  y_end   = c(chargepWT, chargepdel, chargepA, chargepE)
)
# plot
svg(filename = "change_charge.svg", width = 3, height = 2)
ggplot(dcharge, aes(x = mutant, y = values, group = mutant)) +
  geom_point(aes(color = state), size = 4.5, show.legend = FALSE) +
  scale_color_manual(values = c("nonphos" = "black", "phos" = "red")) +
  geom_segment(data = arrows_df,
               aes(x = mutant, xend = mutant, y = y_start, yend = y_end),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "grey50") +
  geom_line(color = "grey50") +
  coord_cartesian(ylim = c(min(dcharge$values) - 1, max(dcharge$values) + 1))+
  labs(x = "Variant", y = "Net charge") +
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.x = element_blank(),
        aspect.ratio = 2 / 3)

dev.off()


### DEP association Based on Charge
#   with charge and association as percentage of CSPmax (WT CSP)
## average values
DEP <- data.frame(charge = c(chargeWT,chargepA, chargepdel, chargepE, chargepWT)*(-1),
                       DEP = c(0,15.17094017,32.26993,62.1618434,100),
                       error = c(0,1.495726496, 5.199553568, 12.15329639,0),
                       color = c("black", "#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))

# Fit a sigmoidal curve (logistic function) using nls()
fit_DEP <- nls(DEP ~ 100 / (1 + exp(-k * (charge - x0))),
           data = DEP,
           start = list(k = 0.1, x0 = 13))


# Print the coefficients and check the fit
x0_value_DEP <- coef(fit_DEP)["x0"]
summary(fit_DEP)

# Get fitted values from the model
DEP$fitted <- predict(fit_DEP)

# Create a fine grid of charge values for a smoother curve
charge_grid_DEP <- seq(0, 30, length.out = 500)

# Predict the activity for each charge in the fine grid
predicted_CSP_DEP <- predict(fit_DEP, newdata = data.frame(charge = charge_grid_DEP))

# Create a new dataframe for plotting the smooth curve
smooth_data_DEP <- data.frame(charge = charge_grid_DEP, DEP = predicted_CSP_DEP)

# Plot original data and fitted sigmoidal curve
svg(filename = "DEP_charge.svg", width = 3, height = 2)
ggplot(DEP, aes(x = charge, y = DEP)) +
  geom_line(data = smooth_data_DEP, aes(x = charge, y = DEP), color = "grey50", linewidth = 1)+
  geom_point(size = 4.5, fill = DEP$color, shape = 21) +  # Plot the original data points
  geom_errorbar(aes(ymin = DEP - error, ymax = DEP + error), width =1)+
  geom_vline(xintercept=x0_value_DEP, linetype = "dotted", color = "black", linewidth = 1)+
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-3, max(DEP$DEP) + 3))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()

## T459 values
T459 <- data.frame(charge = c(chargeWT,chargepA, chargepdel, chargepE, chargepWT)*(-1),
                  T459 = c(0,16.66666667,40.66666667,56,100),
                  color = c("black", "#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))

# Fit a sigmoidal curve (logistic function) using nls()
fit_T459 <- nls(T459 ~ 100 / (1 + exp(-k * (charge - x0))),
               data = T459,
               start = list(k = 0.1, x0 = 13))


# Print the coefficients and check the fit
x0_value_T459 <- coef(fit_T459)["x0"]
summary(fit_T459)

# Get fitted values from the model
DEP$fitted <- predict(fit_T459)

# Create a fine grid of charge values for a smoother curve
charge_grid_T459 <- seq(0, 30, length.out = 500)

# Predict the activity for each charge in the fine grid
predicted_CSP_T459 <- predict(fit_T459, newdata = data.frame(charge = charge_grid_T459))

# Create a new dataframe for plotting the smooth curve
smooth_data_T459 <- data.frame(charge = charge_grid_T459, T459 = predicted_CSP_T459)

# Plot original data and fitted sigmoidal curve
svg(filename = "T459_charge.svg", width = 3, height = 2)
ggplot(T459, aes(x = charge, y = T459)) +
  geom_line(data = smooth_data_T459, aes(x = charge, y = T459), color = "grey50", linewidth = 1)+
  geom_point(size = 4.5, fill = T459$color, shape = 21) +  # Plot the original data points
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-3, max(T459$T459) + 3))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()

## T480 values
T480 <- data.frame(charge = c(chargeWT,chargepA, chargepdel, chargepE, chargepWT)*(-1),
                   T480 = c(0,13.67521368,23.07692308,44.01709402,100),
                   color = c("black", "#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))

# Fit a sigmoidal curve (logistic function) using nls()
fit_T480 <- nls(T480 ~ 100 / (1 + exp(-k * (charge - x0))),
                data = T480,
                start = list(k = 0.1, x0 = 13))


# Print the coefficients and check the fit
x0_value_T480 <- coef(fit_T480)["x0"]
summary(fit_T480)

# Get fitted values from the model
DEP$fitted <- predict(fit_T480)

# Create a fine grid of charge values for a smoother curve
charge_grid_T480 <- seq(0, 30, length.out = 500)

# Predict the activity for each charge in the fine grid
predicted_CSP_T480 <- predict(fit_T480, newdata = data.frame(charge = charge_grid_T480))

# Create a new dataframe for plotting the smooth curve
smooth_data_T480 <- data.frame(charge = charge_grid_T480, T480 = predicted_CSP_T480)

# Plot original data and fitted sigmoidal curve
svg(filename = "T480_charge.svg", width = 3, height = 2)
ggplot(T480, aes(x = charge, y = T480)) +
  geom_line(data = smooth_data_T480, aes(x = charge, y = T480), color = "grey50", linewidth = 1)+
  geom_point(size = 4.5, fill = T480$color, shape = 21) +  # Plot the original data points
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-3, max(T480$T480) + 3))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()

## K483 values
K483 <- data.frame(charge = c(chargeWT,chargepA, chargepdel, chargepE, chargepWT)*(-1),
                   K483 = c(0,75.58685446,34.27230047,71.78398829,100),
                   color = c("black", "#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))

# Fit a sigmoidal curve (logistic function) using nls()
fit_K483 <- nls(K483 ~ 100 / (1 + exp(-k * (charge - x0))),
                data = K483,
                start = list(k = 0.1, x0 = 13))


# Print the coefficients and check the fit
x0_value_K483 <- coef(fit_K483)["x0"]
summary(fit_K483)

# Get fitted values from the model
DEP$fitted <- predict(fit_K483)

# Create a fine grid of charge values for a smoother curve
charge_grid_K483 <- seq(0, 30, length.out = 500)

# Predict the activity for each charge in the fine grid
predicted_CSP_K483 <- predict(fit_K483, newdata = data.frame(charge = charge_grid_K483))

# Create a new dataframe for plotting the smooth curve
smooth_data_K483 <- data.frame(charge = charge_grid_K483, K483 = predicted_CSP_K483)

# Plot original data and fitted sigmoidal curve
svg(filename = "K483_charge.svg", width = 3, height = 2)
ggplot(K483, aes(x = charge, y = K483)) +
  geom_line(data = smooth_data_K483, aes(x = charge, y = K483), color = "grey50", linewidth = 1)+
  geom_point(size = 4.5, fill = K483$color, shape = 21) +  # Plot the original data points
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-3, max(K483$K483) + 3))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()

## S487 values
S487 <- data.frame(charge = c(chargeWT,chargepA, chargepdel, chargepE, chargepWT)*(-1),
                   S487 = c(0,56.17021277,31.06382979,76.84629128,100),
                   color = c("black", "#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))

# Fit a sigmoidal curve (logistic function) using nls()
fit_S487 <- nls(S487 ~ 100 / (1 + exp(-k * (charge - x0))),
                data = S487,
                start = list(k = 0.1, x0 = 13))


# Print the coefficients and check the fit
x0_value_S487 <- coef(fit_S487)["x0"]
summary(fit_S487)

# Get fitted values from the model
DEP$fitted <- predict(fit_S487)

# Create a fine grid of charge values for a smoother curve
charge_grid_S487 <- seq(0, 30, length.out = 500)

# Predict the activity for each charge in the fine grid
predicted_CSP_S487 <- predict(fit_S487, newdata = data.frame(charge = charge_grid_S487))

# Create a new dataframe for plotting the smooth curve
smooth_data_S487 <- data.frame(charge = charge_grid_S487, S487 = predicted_CSP_S487)

# Plot original data and fitted sigmoidal curve
svg(filename = "S487_charge.svg", width = 3, height = 2)
ggplot(S487, aes(x = charge, y = S487)) +
  geom_line(data = smooth_data_S487, aes(x = charge, y = S487), color = "grey50", linewidth = 1)+
  geom_point(size = 4.5, fill = S487$color, shape = 21) +  # Plot the original data points
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-3, max(S487$S487) + 3))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()

