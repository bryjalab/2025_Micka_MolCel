library("ggplot2")
library("here")

## PDZ-IDR2-DEP sequence upload
# X is phospho serine
# Z is phospho threonine
PDZ_DEP <- read.table(here("input", "PDZ-DEP.txt"))
PDZ_DEP <- toString(PDZ_DEP)
pPDZ_DEP <- read.table(here("input", "phosphoPDZ-DEP.txt"))
pPDZ_DEP <- toString(pPDZ_DEP)


PDZ_DEP_deltaST <- read.table(here("input", "PDZ-DEP_deltaST.txt"))
PDZ_DEP_deltaST <- toString(PDZ_DEP_deltaST)
pPDZ_DEP_deltaST <- read.table(here("input", "phosphoPDZ-DEP_deltaST.txt"))
pPDZ_DEP_deltaST <- toString(pPDZ_DEP_deltaST)

Dsh_PDZ_DEP <- read.table(here("input", "DshPDZ-DEP.txt"))
Dsh_PDZ_DEP <- toString(Dsh_PDZ_DEP)
Dsh_pPDZ_DEP <- read.table(here("input", "DshphosphoPDZ-DEP.txt"))
Dsh_pPDZ_DEP <- toString(Dsh_pPDZ_DEP)


## calculate charge
# WT
pPDZ_DEP_charge <- chargeCalculationLocal(pPDZ_DEP, window = 9, pH = 6.5, pKaSet = custom_pka)
PDZ_DEP_charge <- chargeCalculationLocal(PDZ_DEP, window = 9, pH = 6.5, pKaSet = custom_pka)

# delta S/T
pPDZ_DEP_deltaST_charge <- chargeCalculationLocal(pPDZ_DEP_deltaST, window = 9, pH = 6.5, pKaSet = custom_pka)
PDZ_DEP_deltaST_charge <- chargeCalculationLocal(PDZ_DEP_deltaST, window = 9, pH = 6.5, pKaSet = custom_pka)

# Dsh
DshpPDZ_DEP_charge <- chargeCalculationLocal(Dsh_pPDZ_DEP, window = 9, pH = 6.5, pKaSet = custom_pka)
DshPDZ_DEP_charge <- chargeCalculationLocal(Dsh_PDZ_DEP, window = 9, pH = 6.5, pKaSet = custom_pka)

## plot
# WT
PDZ_DEP_charge %>% ggplot() +
  geom_hline(yintercept = 0, linetype = "solid", 
            color = "black", linewidth = 0.5, alpha = 1) +
  geom_line(aes(x = Position+242, y = windowCharge, color= windowCharge), linewidth = 0.5)+
  geom_line(aes(x = pPDZ_DEP_charge$Position+242, y = pPDZ_DEP_charge$windowCharge, color= pPDZ_DEP_charge$windowCharge), linewidth = 0.5)+
  scale_colour_gradient2(low = "red", mid =  "white", high = "blue", limits = c(-1,1), oob = scales::squish)+
  ylim(-1.6,1)+
  scale_x_continuous(breaks = c(250, 300, 350, 400, 450, 500))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "#cccccc"), axis.line = element_line(colour = "black"),
        axis.text.x.top = element_text(angle = 90, vjust = 0.5), legend.position="none") +
  labs(x = "PDZ-DEP sequence",
      y = "Local charge")
ggsave("charge.eps", width = 1000,height = 500,units = "px")


# deletion mutant
PDZ_DEP_deltaST_charge %>% ggplot() +
  geom_hline(yintercept = 0, linetype = "solid", 
             color = "black", linewidth = 0.5, alpha = 1) +
  geom_line(aes(x = Position+242, y = windowCharge, color= windowCharge), linewidth = 0.5)+
  geom_line(aes(x = pPDZ_DEP_deltaST_charge$Position+242, y = pPDZ_DEP_deltaST_charge$windowCharge, color= pPDZ_DEP_deltaST_charge$windowCharge), linewidth = 0.5)+
  scale_colour_gradient2(low = "red", mid =  "white", high = "blue", limits = c(-1,1), oob = scales::squish)+
  ylim(-1.6,1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "#cccccc"), axis.line = element_line(colour = "black"),
        axis.text.x.top = element_text(angle = 90, vjust = 0.5), legend.position="none") +
  labs(x = "PDZ-DEP sequence",
       y = "Local charge")
ggsave("charge_deltaST_nogap.eps", width = 1000,height = 500,units = "px")


# Dsh
DshPDZ_DEP_charge %>% ggplot() +
  geom_hline(yintercept = 0, linetype = "solid", 
             color = "black", linewidth = 0.5, alpha = 1) +
  geom_line(aes(x = Position+245, y = windowCharge, color= windowCharge), linewidth = 0.5)+
  geom_line(aes(x = DshpPDZ_DEP_charge$Position+245, y = DshpPDZ_DEP_charge$windowCharge, color= DshpPDZ_DEP_charge$windowCharge), linewidth = 0.5)+
  scale_colour_gradient2(low = "red", mid =  "white", high = "blue", limits = c(-1,1), oob = scales::squish)+
  ylim(-1.6,1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "#cccccc"), axis.line = element_line(colour = "black"),
        axis.text.x.top = element_text(angle = 90, vjust = 0.5), legend.position="none") +
  labs(x = "PDZ-DEP sequence",
       y = "Local charge")
ggsave("charge_Dsh.eps", width = 1000,height = 500,units = "px")
