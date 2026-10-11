# CHEM 124 Lab 1 - Aspartame
# Melanie, Austin, Jacob
# spectrophotometer: Bryan
# lambda max = 257 nm, cuvette = 1 cm

install.packages("outliers")
library(outliers)

mw <- 294.30    # molar mass of aspartame (g/mol)


############################################
### PART 1: my calibration curve ###
############################################

# stock concentration values and calculation
mass <- 58.6           # mg aspartame
volume <- 0.020        # L (20 mL flask)
stock <- mass / mw / volume     # in mM
stock

# concentrations of my standards (in M), multiplied by (stock/10) to account
# for us using 58.6g aspartame instead of 58.8
conc <- c(0.5, 1, 2, 4, 6, 8, 10) * stock / 10 / 1000
conc

# our absorbances at 257 nm
abs <- c(0.173, 0.232, 0.580, 1.163, 1.754, 2.241, 2.609)

# linear regression with all 7 points, y-intercept set to 0 as instructed
# in lab
fit1 <- lm(abs ~ 0 + conc)
summary(fit1)
# slope = molar absorptivity (M^-1 cm^-1)
# "Std. Error" next to the slope = standard error of the molar absorptivity
# "Residual standard error" = standard error of the regression
# "Multiple R-squared" = R^2

# set epsilon, the molar absorptivity coefficient
eps <- fit1$coefficients[1]
eps_se <- summary(fit1)$coefficients[1, 2]

# graph for our calibration curve
plot(conc * 1000, abs, pch = 19,
     xlab = "Aspartame concentration (mM)",
     ylab = "Absorbance at 257 nm",
     main = "Aspartame calibration curve")
abline(0, eps / 1000)
legend("bottomright", bty = "n", cex = 1.2,
       legend = c(as.expression(bquote(A == .(round(eps / 1000, 4)) * x)),
                  as.expression(bquote(epsilon == .(round(eps, 1)) %+-% .(round(eps_se, 1)) ~ M^-1 ~ cm^-1))))
# x = concentration in mM (matches the x axis)
# x axis is in mM so I divided the slope by 1000 to draw the line
# saved with Export -> Save as Image -> PNG


############################################
### PART 2: Skittles Mango Tangelo drink mix ###
############################################

# packet masses (g)
packets <- c(2.9851, 2.9259, 3.0667)
mean(packets)
sd(packets)
(mean(packets) - 3) / 3 * 100     # percent difference from 3 g

# absorbance of each packet at 3.125 mg mix / mL (P1, P2, P3)
mix_conc <- 3.125
packet_abs <- c(0.377, 1.924, 1.878)

# outlier tests on the 3 absorbances
# Q test by hand: Q = gap / range
sorted_abs <- sort(packet_abs)
Q <- (sorted_abs[2] - sorted_abs[1]) / (sorted_abs[3] - sorted_abs[1])
Q
# Q crit (n = 3): 0.941 at 90%, 0.970 at 95%. Q > Q crit means outlier

dixon.test(packet_abs)    # Q test in R
grubbs.test(packet_abs)   # Grubbs test
# p < 0.05 means the most extreme value is an outlier

# P1 (0.377) is an outlier, so kept P2 and P3
kept_abs <- packet_abs[c(2, 3)]
avg_abs <- mean(kept_abs)
avg_abs
sd(kept_abs)

# aspartame concentration (M) from A = eps * c
asp_M <- avg_abs / eps
asp_M

# change to mg/mL (mol/L * g/mol = g/L = mg/mL)
asp_mg_mL <- asp_M * mw

# mg aspartame per g of mix
asp_mg_g <- asp_mg_mL / mix_conc * 1000
asp_mg_g

# SD and 95% confidence interval from the 2 packets
packet_vals <- kept_abs / eps * mw / mix_conc * 1000
sd(packet_vals)
qt(0.975, 1) * sd(packet_vals) / sqrt(2)

# t-test against estimated aspartame content
# Skittles doesn't list aspartame amount, so using Wyler's Light estimate
# from Open Food Facts: 8.8% aspartame = 88 mg aspartame / g mix
stated_mg_g <- 0.088 * 1000
stated_mg_g
mean(packets) * 0.088 * 1000    # estimated mg aspartame per packet
t.test(packet_vals, mu = stated_mg_g)


############################################
### PART 3: spectrophotometer group analysis (Bryan) ###
############################################
# both groups used the same standard concentrations
# using target concentrations since I don't know the other group's stock mass

conc_nom <- c(0.5, 1, 2, 4, 6, 8, 10) / 1000    # M (dilution 6 -> stock)

bryan1 <- c(0.003, 0.056, 0.185, 0.82, 0.675, 1.146, 1.539)   # Tue Bryan/Gemma
bryan2 <- c(0.173, 0.232, 0.580, 1.163, 1.754, 2.241, 2.609)  # Thu me/Austin/Jacob

# put both groups in one table (each row = one group)
bryan <- rbind(bryan1, bryan2)
bryan

# calibration curve with all Bryan data, y-intercept set to 0
xb <- rep(conc_nom, times = 2)   # concentration for every point
yb <- as.vector(t(bryan))        # every absorbance in the same order
fit_bryan <- lm(yb ~ 0 + xb)
summary(fit_bryan)

eps_b <- fit_bryan$coefficients[1]
eps_b_se <- summary(fit_bryan)$coefficients[1, 2]

# mean and SD at each concentration
means <- colMeans(bryan)
sds <- apply(bryan, 2, sd)
means
sds

# graph 1: averages with error bars (+/- 1 SD)
plot(conc_nom * 1000, means, pch = 19, ylim = c(0, 3),
     xlab = "Aspartame concentration (mM)",
     ylab = "Absorbance at 257 nm",
     main = "Calibration curve, Bryan (mean +/- SD)")
arrows(conc_nom * 1000, means - sds, conc_nom * 1000, means + sds,
       angle = 90, code = 3, length = 0.05)
abline(0, eps_b / 1000)
legend("bottomright", bty = "n", cex = 1.2,
       legend = c(as.expression(bquote(A == .(round(eps_b / 1000, 4)) * x)),
                  as.expression(bquote(epsilon == .(round(eps_b, 1)) %+-% .(round(eps_b_se, 1)) ~ M^-1 ~ cm^-1))))
# x = concentration in mM (matches the x axis)

# graph 2: every single data point
plot(xb * 1000, yb, pch = 19, ylim = c(0, 3),
     col = rep(c("blue", "red"), each = 7),
     xlab = "Aspartame concentration (mM)",
     ylab = "Absorbance at 257 nm",
     main = "Calibration curve, Bryan (individual points)")
abline(0, eps_b / 1000)
legend("bottomright", bty = "n", cex = 1.2,
       legend = c(as.expression(bquote(A == .(round(eps_b / 1000, 4)) * x)),
                  as.expression(bquote(epsilon == .(round(eps_b, 1)) %+-% .(round(eps_b_se, 1)) ~ M^-1 ~ cm^-1))))
# x = concentration in mM (matches the x axis)
legend("topleft", c("Tuesday group", "Thursday group (mine)"),
       col = c("blue", "red"), pch = 19)

# Grubbs test on the residuals
# only 2 values per concentration, so instead I tested the residuals
# (distance from each point to the line) from all 14 points
res <- residuals(fit_bryan)
res
grubbs.test(res)
# p < 0.05 means the most extreme point is an outlier from the line

# re-calculate my drink mix with the Bryan curve
asp_M_b <- avg_abs / eps_b
asp_mg_g_b <- asp_M_b * mw / mix_conc * 1000
asp_mg_g_b

packet_vals_b <- kept_abs / eps_b * mw / mix_conc * 1000
sd(packet_vals_b)
qt(0.975, 1) * sd(packet_vals_b) / sqrt(2)    # 95% CI

# t-test with the Bryan curve against the same estimate
t.test(packet_vals_b, mu = stated_mg_g)