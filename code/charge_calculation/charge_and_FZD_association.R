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


### FZDs association Based on Charge
## FZD3
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
    breaks = c(0, 5, 10, 15, 20, 25, 30),  # Specify the x-ticks you want to label
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "30"))+
  coord_cartesian(ylim = c(-1, max(FZD3$FZD3) + 1))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)

dev.off()


## FZD5
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
    breaks = c(0, 5, 10, 15, 20, 25,30),  # Specify the x-ticks you want to label
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-1, max(FZD5$FZD5) + 1))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()

## FZD6
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
    breaks = c(0, 5, 10, 15, 20, 25, 30),  # Specify the x-ticks you want to label
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-1, max(FZD6$FZD6) + 1))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()

### Axin association Based on Charge
# Dataframe with charge and association
Axin <- data.frame(charge = c(chargepA_7.2, chargepdel_7.2, chargepE_7.2, chargepWT_7.2)*(-1),
                   Axin1 = c(2.60211044690386, 3.82221595235317, 0.631077893440392, 0),
                   Axin2 = c(3.07712341376974, 1.20480803148565, 0.753491003788647, 0),
                   color = c("#e7298aff","#1b9e77ff","#66a61eff","#d95f02ff"))


# Plot original data
svg(filename = "Axin_charge.svg", width = 3, height = 2)
ggplot(Axin, aes(x = charge, y = Axin1)) +
  geom_point(size = 4.5, fill = Axin$color, shape = 22) +  # Plot Axin1 data points
  geom_point(aes(x = charge, y = Axin2), size = 4.5, fill = Axin$color, shape = 24) + # Plot Axin2 data points
  labs(x = "Charge")+
  scale_x_continuous(
    limits = c(0, 30),
    breaks = c(0, 5, 10, 15, 20, 25, 30),  # Specify the x-ticks you want to label
    labels = c("0", "-5", "-10", "-15", "-20", "-25", "-30"))+
  coord_cartesian(ylim = c(-1, max(FZD6$FZD6)+1))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.y = element_blank(),
        aspect.ratio = 2 / 3)
dev.off()