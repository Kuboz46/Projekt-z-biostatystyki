# Biblioteki

library(survival)
library(foreign)
library(rms)
library(survMisc)
library(reshape2)
library(survminer)
library(flexsurv)

# Wczytanie danych
setwd("C:/Users/Kuboz/Documents/Biostatystyka/Bie¿¹cy rok/Projekty/Projekt koñcowy")
df <- read.csv("pharynx.csv")
# View(df)
head(df)

nrow(df)
# 195

sum(is.na(df))
# 0

length(unique(df$CASE))
# 195

# Zmienn¹ "CASE" usuwam, poniewa¿ ka¿dy identyfikator chorego jest równy numerowi obserwacji.
df <- df[, -1]

unique(df$GRADE)
# 1 2 3 9

unique(df$INST)
# 2 5 4 6 3 1



# Porównania funkcji prze¿ycia dla obu grup

# SEX

unique(df$SEX)
# 2 1

KM_SEX <- survfit(Surv(TIME, STATUS) ~ SEX, data = df, conf.type = "none")

X11()
par(mfrow = c(1, 3))
plot(KM_SEX, col = c("blue", "red"), main = "P³eæ")
legend("topright", legend = c("mê¿czyzna", "kobieta"), lty=1, col = c("blue", "red"))

(test_SEX <- survdiff(Surv(TIME, STATUS) ~ SEX, data = df)) 
# P-wartoœæ: p= 0.4


# T_STAGE

unique(df$T_STAGE)
# 3 2 4 1

KM_T_STAGE <- survfit(Surv(TIME, STATUS) ~ T_STAGE, data = df, conf.type = "none")

plot(KM_T_STAGE, col = c("black", "red", "green", "blue"), main = "Wielkoœæ guza")
legend("topright", legend = c("1", "2", "3", "4"), lty=1, col = c("black", "red", "green", "blue"))

(test_T_STAGE <- survdiff(Surv(TIME, STATUS) ~ T_STAGE, data = df)) 
# P-wartoœæ: p= 0.01

# N_STAGE

unique(df$N_STAGE)
# 1 3 0 2

KM_N_STAGE <- survfit(Surv(TIME, STATUS) ~ N_STAGE, data = df, conf.type = "none")

plot(KM_N_STAGE, col = c("black", "red", "green", "blue"), main = "Przerzuty do wêz³ów ch³onnych")
legend("topright", legend = c("0", "1", "2", "3"), lty=1, col = c("black", "red", "green", "blue"))

(test_N_STAGE <- survdiff(Surv(TIME, STATUS) ~ N_STAGE, data = df)) 
# P-wartoœæ: p= 0.01


# Dopasowanie modelu Coxa

# Wiek
ph_pusty <- coxph(Surv(TIME, STATUS)~1, data = df)
reszty_mart <- resid(ph_pusty)
X11()
plot(df$AGE, reszty_mart, xlab = "Wiek w chwili diagnozy", ylab = "Wartoœæ reszty", main = "Reszty martynga³owe vs wiek w chwili diagnozy")
lines(lowess(df$AGE, reszty_mart,iter=0,f=0.5))
# Zatem mo¿na dodaæ wiek jako funkcjê liniow¹.


mod_PH_full <- coxph(Surv(TIME, STATUS) ~ as.factor(INST) + as.factor(SEX) + TX + GRADE + AGE + 
                       COND + as.factor(SITE) + T_STAGE + as.factor(N_STAGE), data = df)

print(mod_PH_full)

(test_full <- cox.zph(mod_PH_full))
# Ok

# reszty dewiancji
ph_full_dev_res <- residuals(mod_PH_full, type = "deviance")
fit_val_full <- mod_PH_full$linear.predictors
X11()
plot(fit_val_full, ph_full_dev_res, ylim=c(-3.5,3.5), ylab="",xlab="")
abline(h=0, col=2); abline(h=-3, col=2, lty=2); abline(h=3, col=2, lty=2)


X11()
ggcoxdiagnostics(mod_PH_full, type='deviance', linear.predictions = FALSE)

(length(which(ph_full_dev_res > 2)) + length(which(ph_full_dev_res < -2))) / nrow(df) * 100
# 6.153846 %



X11()
plot(test_full[5], df=4, nsmo=10, se=TRUE, xlab = "Czas", ylab = " ")
title("Beta(t) dla wieku")
abline(mod_PH_full$coeff[9], 0, lty=17, col = "red")




quantile(df$AGE, c(1/4, 1/2, 3/4))
# 25% 50% 75% 
# 52  60  68 

df2 <- df
df2$AGE2 <- rep(0, nrow(df2))
df2$AGE2[which(df2$AGE <= 52)] <- 1
df2$AGE2[which(df2$AGE > 52 & df2$AGE <= 60)] <- 2
df2$AGE2[which(df2$AGE > 60 & df2$AGE <= 68)] <- 3
df2$AGE2[which(df2$AGE > 68)] <- 4

head(df2)
df2 <- df2[, -5]

mod_PH_skat_AGE <- coxph(Surv(TIME, STATUS) ~ as.factor(INST) + as.factor(SEX) + TX + GRADE + AGE2 + 
                           COND + as.factor(SITE) + T_STAGE + as.factor(N_STAGE), data = df2)
summary(mod_PH_skat_AGE)

(test_skat_AGE <- cox.zph(mod_PH_skat_AGE))

X11()
plot(test_skat_AGE[5], df=4, nsmo=10, se=TRUE, xlab = "Czas", ylab = " ")
title("Beta(t) dla skategoryzowanego wieku")
abline(mod_PH_skat_AGE$coeff[9], 0, lty=17, col = "red")

AIC(mod_PH_skat_AGE)
# 1325.577

# Warstwowanie po skategoryzowanym wieku
mod_PH_warstw_AGE <- coxph(Surv(TIME, STATUS) ~ as.factor(INST) + as.factor(SEX) + TX + GRADE + strata(AGE2) + 
                           COND + as.factor(SITE) + T_STAGE + as.factor(N_STAGE), data = df2)
summary(mod_PH_warstw_AGE)

(test_warstw_AGE <- cox.zph(mod_PH_warstw_AGE))

AIC(mod_PH_warstw_AGE)
# 936.9776



# cluster

mod_cluster <- coxph(Surv(TIME, STATUS) ~ cluster(INST) + as.factor(SEX) + TX + GRADE + AGE + 
                       COND + as.factor(SITE) + T_STAGE + as.factor(N_STAGE), data = df)

summary(mod_cluster)

(test_cluster <- cox.zph(mod_cluster))
# Ok

# reszty dewiancji
ph_cluster_dev_res <- residuals(mod_cluster, type = "deviance")
fit_val_cluster <- mod_cluster$linear.predictors
X11()
plot(fit_val_cluster, ph_cluster_dev_res, ylim=c(-3.5,3.5), ylab="",xlab="")
abline(h=0, col=2); abline(h=-3, col=2, lty=2); abline(h=3, col=2, lty=2)


X11()
plot(test_cluster[4], df=4, nsmo=10, se=TRUE, xlab = "Czas", ylab = " ")
title("Beta(t) dla wieku")
abline(mod_cluster$coeff[4], 0, lty=17, col = "red")


X11()
ggcoxdiagnostics(mod_cluster, type='deviance', linear.predictions = FALSE)

(length(which(ph_cluster_dev_res > 2)) + length(which(ph_cluster_dev_res < -2))) / nrow(df) * 100
# 6.153846

#### SELEKCJA ZMIENNYCH TOP-DOWN ####
(top_down_ph <- step(mod_PH_full, scope=list(upper=~ ., lower=~1) , data = df,direction = "backward"))

(test_top_down1 <- cox.zph(top_down_ph, transform = "identity"))

AIC(top_down_ph)
# 1310.127

AIC(mod_PH_full)
# 1325.844

AIC(mod_cluster)
# 1317.762





# A model cluster z warstwowaniem po wieku?
mod_PH_warstw_AGE_cluster <- coxph(Surv(TIME, STATUS) ~ cluster(INST) + as.factor(SEX) + TX + GRADE + strata(AGE2) + 
                             COND + as.factor(SITE) + T_STAGE + as.factor(N_STAGE), data = df2)
summary(mod_PH_warstw_AGE_cluster)

(test_warstw_AGE_cluster <- cox.zph(mod_PH_warstw_AGE_cluster))

AIC(mod_PH_warstw_AGE_cluster)
# 928.796


# Dopasowanie modelu AFT

# Dopasowanie uogolnionego modelu F

ogolny_F <- flexsurvreg(Surv(TIME, STATUS) ~ as.factor(INST) + as.factor(SEX) + TX + GRADE + AGE + COND + 
                          as.factor(SITE) + T_STAGE + as.factor(N_STAGE), data = df, dist = "genf")
# Komunikat ostrzegawczy:
# W poleceniu 'flexsurvreg(Surv(TIME, STATUS) ~ as.factor(INST) + as.factor(SEX) + ':
# Optimisation has probably not converged to the maximum likelihood - Hessian is not positive definite. 

ogolny_F <- flexsurvreg(Surv(TIME, STATUS) ~ INST + SEX + TX + GRADE + AGE + 
                          COND + SITE + T_STAGE + N_STAGE, data = df, dist = "genf")
ogolny_F



# Log-logistyczny

ogolny_LL <- flexsurvreg(Surv(TIME, STATUS) ~ INST + SEX + TX + GRADE + AGE + 
                           COND + SITE + T_STAGE + N_STAGE, data = df ,dist="genf",
                         inits=c(3,.2,0,1,1,1,1,1,1),fixedpars = c(3,4))

(LL <- 1 - pchisq(2 * (ogolny_F$loglik - ogolny_LL$loglik), 2))
# 0.0002124637

# Log-normalny

ogolny_LN <- flexsurvreg(Surv(TIME, STATUS) ~ INST + SEX + TX + GRADE + AGE + 
                           COND + SITE + T_STAGE + N_STAGE, data = df, dist = "lnorm")
(LN <- 1 - pchisq(2 * (ogolny_F$loglik - ogolny_LN$loglik), 2))
# 6.544084e-05

# Weibull

ogolny_W <- flexsurvreg(Surv(TIME, STATUS) ~ INST + SEX + TX + GRADE + AGE + 
                          COND + SITE + T_STAGE + N_STAGE, data = df, dist = "weibull")
(W <- 1 - pchisq(2 * (ogolny_F$loglik - ogolny_W$loglik), 2))
# 2.498289e-09

wyniki <- data.frame("log-log" = LL, "log-norm" = LN, "Weibull" = W)

wyniki




# A rozklad uogolniony Gamma?
ogolny_GG <- flexsurvreg(Surv(TIME, STATUS) ~ INST + SEX + TX + GRADE + AGE + 
                           COND + SITE + T_STAGE + N_STAGE, data = df, dist = "gengamma")
(GG <- 1 - pchisq(2 * (ogolny_F$loglik - ogolny_GG$loglik), 2))
# 1



X11()
par(mfrow=c(1,6))
AFT.Wei1<-psm(Surv(TIME, STATUS) ~ INST + SEX + TX + GRADE + AGE + 
                COND + SITE + T_STAGE + N_STAGE,data=df,dist="weibull",y=TRUE)
AFT.LL1<-psm(Surv(TIME, STATUS) ~ INST + SEX + TX + GRADE + AGE + 
               COND + SITE + T_STAGE + N_STAGE, data=df, dist="loglogistic",y=TRUE)
AFT.LogN1<-psm(Surv(TIME, STATUS) ~ INST + SEX + TX + GRADE + AGE + 
                 COND + SITE + T_STAGE + N_STAGE, data=df, dist="lognormal",y=TRUE)
res.Weib1 <-resid(AFT.Wei1,type="cens")
survplot(npsurv(res.Weib1 ~1),conf="none", xlab="Rezydua",ylab="P-stwo prze¿ycia (Gumbela)")
lines(res.Weib1)
res.LL1 <-resid(AFT.LL1,type="cens")
survplot(npsurv(res.LL1 ~1),conf="none", xlab="Rezydua",ylab="P-stwo prze¿ycia (logistyczny)")
lines(res.LL1)
res.LogN1 <-resid(AFT.LogN1,type="cens")
survplot(npsurv(res.LogN1 ~1),conf="none", xlab="Rezydua",ylab="P-stwo prze¿ycia (normalny)")
lines(res.LogN1)



# Wykresy reszt Coxa-Snell

weibull_CS <- psm(Surv(TIME, STATUS) ~ INST + SEX + TX + GRADE + AGE + 
                 COND + SITE + T_STAGE + N_STAGE,data=df, dist="weibull", y=TRUE)
res_weibull <-resid(weibull_CS,type="cens") # standaryzowane reszty
rozklad.reszt_weib<-npsurv(res_weibull ~1)
coxsnellres_weib<-(-log(rozklad.reszt_weib$surv))#-rozklad.reszt$n.censor*log(2)
y<-(-log(1-pnorm(rozklad.reszt_weib$time)))
plot(y,coxsnellres_weib,pch=16,cex=0.7,xlab="Reszty Coxa-Snell",ylab="-ln S(reszty Cox-Snell)",ylim=c(0,5))
abline(0,1, col="red", lty=2)
title("Reszty Coxa-Snell - Weibull")


loglogistic_CS <- psm(Surv(TIME, STATUS) ~ INST + SEX + TX + GRADE + AGE + 
                      COND + SITE + T_STAGE + N_STAGE,data=df, dist="loglogistic", y=TRUE)
res_loglogistic <-resid(loglogistic_CS,type="cens") # standaryzowane reszty
rozklad.reszt_loglogistic<-npsurv(res_loglogistic ~1)
coxsnellres_loglogistic<-(-log(rozklad.reszt_loglogistic$surv))#-rozklad.reszt$n.censor*log(2)
y<-(-log(1-pnorm(rozklad.reszt_loglogistic$time)))
plot(y,coxsnellres_loglogistic,pch=16,cex=0.7,xlab="Reszty Coxa-Snell",ylab="-ln S(reszty Cox-Snell)",ylim=c(0,5))
abline(0,1, col="red", lty=2)
title("Reszty Coxa-Snell - Log-logistyczny")

lognormal_CS <- psm(Surv(TIME, STATUS) ~ INST + SEX + TX + GRADE + AGE + 
                      COND + SITE + T_STAGE + N_STAGE,data=df, dist="lognormal", y=TRUE)
res_lognormal <-resid(lognormal_CS,type="cens") # standaryzowane reszty
rozklad.reszt_lognormal<-npsurv(res_lognormal ~1)
coxsnellres_lognormal<-(-log(rozklad.reszt_lognormal$surv))#-rozklad.reszt$n.censor*log(2)
y<-(-log(1-pnorm(rozklad.reszt_lognormal$time)))
plot(y,coxsnellres_lognormal,pch=16,cex=0.7,xlab="Reszty Coxa-Snell",ylab="-ln S(reszty Cox-Snell)",ylim=c(0,5))
abline(0,1, col="red", lty=2)
title("Reszty Coxa-Snell - Log-normalny")


# Model ostateczny - diagnostyka
# Wybieram model z warstwowaniem po AGE i cluster po INST.

# reszty dewiancji
ph_warstw_AGE_cluster_dev_res <- residuals(mod_PH_warstw_AGE_cluster, type = "deviance")
fit_val_warstw_AGE_cluster <- mod_PH_warstw_AGE_cluster$linear.predictors
X11()
plot(fit_val_warstw_AGE_cluster, ph_warstw_AGE_cluster_dev_res, ylim=c(-3.5,3.5), ylab="",xlab="")
abline(h=0, col=2); abline(h=-3, col=2, lty=2); abline(h=3, col=2, lty=2)


X11()
ggcoxdiagnostics(mod_PH_warstw_AGE_cluster, type='deviance', linear.predictions = FALSE)

(length(which(ph_warstw_AGE_cluster_dev_res > 2)) + length(which(ph_warstw_AGE_cluster_dev_res < -2))) / nrow(df2) * 100
# 4.102564
