library(here)
library(readxl)
library(dplyr)
library(ggplot2)
library(pROC)
library(furniture)
library(tidyverse)
here()


###########
# GLM model develop and validate
###########
#
# dat_f = data

df_val <- dat_f %>% filter(year>2011)
df <- dat_f %>% filter(year<2012)

# vars <- c("age", "sex", "diabetes", "hist_pd", "hist_nd", 
#           "dialysis", "preop_state", "angina", "excard_arter",
#           "prev_surg", "endocar", "lv_mod", "lv_poor1", "interval_mi_surg", 
#           "hypertension", "urg_emerg", "urg_urg", "urg_salv", "thoracic_aorta",
#           "outcome","year", "month")

vars <- c("age", "sex", "diabetes", "hist_pd", "hist_nd", 
          "preop_state", "angina", "excard_arter",
          "prev_surg", "endocar", "lv_mod", "lv_poor1", "interval_mi_surg", 
          "urg_emerg", "urg_urg", "thoracic_aorta",
          "outcome","year", "month")
x <- df[vars]
x <- x %>% replace(is.na(.), 0)
dat <- na.omit(x)
y <- dat$outcome
x1 <- as.matrix(sapply(dat[,1:(dim(dat)[2]-3)], as.numeric))

# Run model
fit <- glm(y ~ x1, family = binomial(link = logit))
summary(fit)
#save coefs
fit_coef <- t(fit$coefficients)


# add intercept term to val data
x <- df_val[vars]
dat <- na.omit(x)
y <- dat$outcome
x1 <- as.matrix(sapply(dat[,1:(dim(dat)[2]-3)], as.numeric))
x1 <- cbind(1,x1)

# Multiply each row of data by coef vector
p_glm <- mapply(FUN = `*`, as.data.frame(x1), fit_coef)

# Sum over rows
glm_pred <- as.data.frame(rowSums(p_glm))
names(glm_pred)[1] <- "prob"

# Transfor to risk
glm_pred$risk <- exp(glm_pred$prob)/(1+exp(glm_pred$prob))

glm_pred <- cbind(glm_pred,dat$outcome,dat$year,dat$month, dat$year+dat$month/12-1/12,
                  (dat$year-2012)*12+dat$month)

names(glm_pred)[3] <- "outcome"
names(glm_pred)[4] <- "year"
names(glm_pred)[5] <- "month"
names(glm_pred)[6] <- "time"
names(glm_pred)[7] <- "time_m"


# Check response is between 0 and 1
summary(glm_pred$risk) # looks correct


####################
#### Validation ####
####################

# Group data by decile of predicted risk
glm_pred$grp100 <- ntile(glm_pred$risk, 100)

m100<- group_by(glm_pred, grp100) %>%
  summarise(m=mean(risk), pr = m*100)
m100 <- as.data.frame(m100)

ob_r100 <- group_by(glm_pred, grp100) %>%
  summarise(table(outcome)[2], or= ifelse(is.na(table(outcome)[2])==TRUE,0,table(outcome)[2]/length(outcome)*100))
t100_m <- cbind(ob_r100, m100) 

## Create calibration plots for the 100 quantiles
t100_all <- t100_m[-1]
cal_all <- ggplot(t100_all, aes(x=pr, y=or)) + geom_point(aes(y = or, x = pr)) +
  geom_abline(aes(intercept=0, slope=1)) + xlab("Predicted risk") + ylab("Observed risk")
cal_all


# c-slope
mod1 <- glm(outcome ~ prob, family="binomial",x=TRUE,y=TRUE, data= glm_pred)
#mod1 # 
mod1_ci <- confint(mod1)

# CITL
mod2 <- glm(outcome~offset(prob),family="binomial", data = glm_pred)
mod2_ci <- confint(mod2)

# Discrimination
g <- roc(outcome ~ risk, data = glm_pred, direction = "<")
plot(g)
g
g$auc # AUC is 0.9073
gci <- ci(roc(outcome ~ risk, data = glm_pred, direction = "<"))

# E/O
prop.table(table(glm_pred$outcome))
O <- prop.table(table(glm_pred$outcome))[2]*length(glm_pred$outcome)
E <- mean(glm_pred$risk)*length(glm_pred$outcome)
o_e_u <- ((sqrt(O)+1.96*0.5)^2)/E
o_e_l <- ((sqrt(O)-1.96*0.5)^2)/E

bs <- 1/length(glm_pred[,1]) * sum((glm_pred$risk - glm_pred$outcome)^2) 

glm_val_alldat <- c("12-18", mod2$coefficients, mod2_ci[1],
                    mod2_ci[2], mod1$coefficients[2],
                    mod1_ci[2,1], mod1_ci[2,2],
                    O/E, o_e_u, o_e_l, bs, g$auc, gci[1],gci[3])


glm_val_results <- c()
for(i in 2012:2019){
  dat <- glm_pred %>% filter(year==i)
  mod1 <- glm(outcome ~ prob, family="binomial",x=TRUE,y=TRUE, data= dat)
  c_slope <- mod1$coefficients[2]
  mod1_ci <- confint(mod1)
  mod2 <- glm(outcome~offset(prob),family="binomial", data = dat)
  citl <- mod2$coefficients
  mod2_ci <- confint(mod2)
  g <- roc(outcome ~ risk, data = dat, direction = "<")
  disc <- g$auc
  gci <- ci(roc(outcome ~ risk, data = dat, direction = "<"))
  O <- prop.table(table(dat$outcome))[2]*length(dat$outcome)
  E <- mean(dat$risk)*length(dat$outcome)
  o_e <- O/E
  o_e_u <- ((sqrt(O)+1.96*0.5)^2)/E
  o_e_l <- ((sqrt(O)-1.96*0.5)^2)/E
  bs <- 1/length(dat[,1]) * sum((dat$risk - dat$outcome)^2) 
  glm_val_results <- rbind(glm_val_results, c(i,citl, mod2_ci[1], mod2_ci[2],
                                              c_slope, mod1_ci[2,1], mod1_ci[2,2],
                                              o_e, o_e_u, o_e_l, bs, disc, gci[1],gci[3]))
}

glm_val_results <- rbind(glm_val_results, glm_val_alldat)
glm_val_results <- as.data.frame(glm_val_results)
names(glm_val_results)[9] <- "2"
names(glm_val_results)[10] <- "3"
old_names <- names(glm_val_results)
new_names <- c("year", "citl", "citl_li", "citl_ui", "c_slope", "c_slope_li", "c_slope_ui",
               "o_e", "o_e_ui", "o_e_li", "brier", "disc", "disc_li", "disc_ui")
glm_val_results <- glm_val_results %>% data.table::setnames(old = old_names, new = new_names)
# view(glm_val_results)


gdf <- glm_val_results[-9,]
gdf[] <- lapply(gdf, function(x) as.numeric(as.character(x)))
gdf$mod <- "glm"

#########
# GLM yearly update
#######

x <- df[vars]
x <- x %>% replace(is.na(.), 0)
dat <- na.omit(x)
y <- dat$outcome
x1 <- as.matrix(sapply(dat[,1:(dim(dat)[2]-3)], as.numeric))

# Run model
fit <- glm(y ~ x1, family = binomial(link = logit))
summary(fit)
#save coefs
fit_coef <- t(fit$coefficients)


#


# add intercept term to val data
x <- df_val[vars]
x <- x %>% replace(is.na(.), 0)
dat <- na.omit(x)
y <- dat$outcome
x1 <- as.matrix(sapply(dat[,1:(dim(dat)[2]-3)], as.numeric))
x1 <- cbind(1,x1)

# Multiply each row of data by coef vector
p_glm <- mapply(FUN = `*`, as.data.frame(x1), fit_coef)

# Sum over rows to get LP
lp <- rowSums(p_glm)
df_val$lp11 <- lp

df_val12 <- df_val %>% filter(year==2012)
#df_val13 <- df_val %>% filter(year==2013)
#df_val14 <- df_val %>% filter(year==2014)
#df_val15 <- df_val %>% filter(year==2015)
#df_val16 <- df_val %>% filter(year==2016)
#df_val17 <- df_val %>% filter(year==2017)
#df_val18 <- df_val %>% filter(year==2018)


# in 2012 data lp11 is lp
#df_val12$lp <- df_val12$lp11


fit <- glm(outcome ~ lp11, family = binomial(link = logit), data = df_val12)
#save coefs
fit_coef12 <- t(fit$coefficients)


x1 <- cbind(1,df_val$lp11)
p_glm <- mapply(FUN = `*`, as.data.frame(x1), fit_coef12)
lp <- rowSums(p_glm)
df_val$lp12 <- lp

# now update using 2013 data and predict LP in all val data
df_val13 <- df_val %>% filter(year==2013)

fit <- glm(outcome ~ lp12, family = binomial(link = logit), data = df_val13)
#save coefs
fit_coef13 <- t(fit$coefficients)

x1 <- cbind(1,df_val$lp12)
p_glm <- mapply(FUN = `*`, as.data.frame(x1), fit_coef13)
lp <- rowSums(p_glm)
df_val$lp13 <- lp


# now update using 2014 data and predict LP in all val data
df_val14 <- df_val %>% filter(year==2014)

fit <- glm(outcome ~ lp13, family = binomial(link = logit), data = df_val14)
#save coefs
fit_coef14 <- t(fit$coefficients)

x1 <- cbind(1,df_val$lp13)
p_glm <- mapply(FUN = `*`, as.data.frame(x1), fit_coef14)
lp <- rowSums(p_glm)
df_val$lp14 <- lp


# now update using 2015 data and predict LP in all val data
df_val15 <- df_val %>% filter(year==2015)

fit <- glm(outcome ~ lp14, family = binomial(link = logit), data = df_val15)
#save coefs
fit_coef15 <- t(fit$coefficients)

x1 <- cbind(1,df_val$lp14)
p_glm <- mapply(FUN = `*`, as.data.frame(x1), fit_coef15)
lp <- rowSums(p_glm)
df_val$lp15 <- lp


# now update using 2016 data and predict LP in all val data
df_val16 <- df_val %>% filter(year==2016)

fit <- glm(outcome ~ lp15, family = binomial(link = logit), data = df_val16)
#save coefs
fit_coef16 <- t(fit$coefficients)

x1 <- cbind(1,df_val$lp15)
p_glm <- mapply(FUN = `*`, as.data.frame(x1), fit_coef16)
lp <- rowSums(p_glm)
df_val$lp16 <- lp


# now update using 2017 data and predict LP in all val data
df_val17 <- df_val %>% filter(year==2017)

fit <- glm(outcome ~ lp16, family = binomial(link = logit), data = df_val17)
#save coefs
fit_coef17 <- t(fit$coefficients)

x1 <- cbind(1,df_val$lp16)
p_glm <- mapply(FUN = `*`, as.data.frame(x1), fit_coef17)
lp <- rowSums(p_glm)
df_val$lp17 <- lp


# now update using 2018 data and predict LP in all val data
df_val18 <- df_val %>% filter(year==2018)

fit <- glm(outcome ~ lp17, family = binomial(link = logit), data = df_val18)
#save coefs
fit_coef18 <- t(fit$coefficients)

x1 <- cbind(1,df_val$lp17)
p_glm <- mapply(FUN = `*`, as.data.frame(x1), fit_coef18)
lp <- rowSums(p_glm)
df_val$lp18 <- lp

# now update using 2019 data and predict LP in all val data
df_val19 <- df_val %>% filter(year==2019)

fit <- glm(outcome ~ lp18, family = binomial(link = logit), data = df_val19)
#save coefs
fit_coef19 <- t(fit$coefficients)

x1 <- cbind(1,df_val$lp18)
p_glm <- mapply(FUN = `*`, as.data.frame(x1), fit_coef19)
lp <- rowSums(p_glm)
df_val$lp19 <- lp

names(df_val)


df_val$prob <- ifelse(df_val$year==2012, df_val$lp11, 99)
df_val$prob <- ifelse(df_val$year==2013, df_val$lp12, df_val$prob)
df_val$prob <- ifelse(df_val$year==2014, df_val$lp13, df_val$prob)
df_val$prob <- ifelse(df_val$year==2015, df_val$lp14, df_val$prob)
df_val$prob <- ifelse(df_val$year==2016, df_val$lp15, df_val$prob)
df_val$prob <- ifelse(df_val$year==2017, df_val$lp16, df_val$prob)
df_val$prob <- ifelse(df_val$year==2018, df_val$lp17, df_val$prob)
df_val$prob <- ifelse(df_val$year==2019, df_val$lp18, df_val$prob)


glm_pred_up <- as.data.frame(df_val$prob)
names(glm_pred_up)[1] <- "prob"

# Transfor to risk
glm_pred_up$risk <- exp(glm_pred_up$prob)/(1+exp(glm_pred_up$prob))

glm_pred_up <- cbind(glm_pred_up,dat$outcome,dat$year,dat$month, dat$year+dat$month/12-1/12,
                     (dat$year-2012)*12+dat$month)

names(glm_pred_up)[3] <- "outcome"
names(glm_pred_up)[4] <- "year"
names(glm_pred_up)[5] <- "month"
names(glm_pred_up)[6] <- "time"
names(glm_pred_up)[7] <- "time_m"


# Check response is between 0 and 1
summary(glm_pred_up$risk) # looks correct


# Group data by decile of predicted risk
glm_pred_up$grp100 <- ntile(glm_pred_up$risk, 100)

m100<- group_by(glm_pred_up, grp100) %>%
  summarise(m=mean(risk), pr = m*100)
m100 <- as.data.frame(m100)

ob_r100 <- group_by(glm_pred_up, grp100) %>%
  summarise(table(outcome)[2], or= ifelse(is.na(table(outcome)[2])==TRUE,0,table(outcome)[2]/length(outcome)*100))
t100_m <- cbind(ob_r100, m100) 

## Create calibration plots for the 100 quantiles
t100_all <- t100_m[-1]
cal_all <- ggplot(t100_all, aes(x=pr, y=or)) + geom_point(aes(y = or, x = pr)) +
  geom_abline(aes(intercept=0, slope=1)) + xlab("Predicted risk") + ylab("Observed risk")
cal_all


# c-slope
mod1 <- glm(outcome ~ prob, family="binomial",x=TRUE,y=TRUE, data= glm_pred_up)
#mod1 # 
mod1_ci <- confint(mod1)

# CITL
mod2 <- glm(outcome~offset(prob),family="binomial", data = glm_pred_up)
mod2_ci <- confint(mod2)

# Discrimination
g <- roc(outcome ~ risk, data = glm_pred_up, direction = "<")
plot(g)
g
g$auc # AUC is 0.9073
gci <- ci(roc(outcome ~ risk, data = glm_pred_up, direction = "<"))

# E/O
O <- prop.table(table(glm_pred_up$outcome))[2]*length(glm_pred_up$outcome)
E <- mean(glm_pred_up$risk)*length(glm_pred_up$outcome)
o_e_u <- ((sqrt(O)+1.96*0.5)^2)/E
o_e_l <- ((sqrt(O)-1.96*0.5)^2)/E

bs <- 1/length(glm_pred_up[,1]) * sum((glm_pred_up$risk - glm_pred_up$outcome)^2) 

glm_val_alldat <- c("12-18", mod2$coefficients, mod2_ci[1],
                    mod2_ci[2], mod1$coefficients[2],
                    mod1_ci[2,1], mod1_ci[2,2],
                    O/E, o_e_u, o_e_l, bs, g$auc, gci[1],gci[3])

# Loop generating validation results for each year separately 
glm_val_results <- c()
for(i in 2012:2019){
  dat <- glm_pred_up %>% filter(year==i)
  mod1 <- glm(outcome ~ prob, family="binomial",x=TRUE,y=TRUE, data= dat)
  c_slope <- mod1$coefficients[2]
  mod1_ci <- confint(mod1)
  mod2 <- glm(outcome~offset(prob),family="binomial", data = dat)
  citl <- mod2$coefficients
  mod2_ci <- confint(mod2)
  g <- roc(outcome ~ risk, data = dat, direction = "<")
  disc <- g$auc
  gci <- ci(roc(outcome ~ risk, data = dat, direction = "<"))
  O <- prop.table(table(dat$outcome))[2]*length(dat$outcome)
  E <- mean(dat$risk)*length(dat$outcome)
  o_e <- O/E
  o_e_u <- ((sqrt(O)+1.96*0.5)^2)/E
  o_e_l <- ((sqrt(O)-1.96*0.5)^2)/E
  bs <- 1/length(dat[,1]) * sum((dat$risk - dat$outcome)^2) 
  glm_val_results <- rbind(glm_val_results, c(i,citl, mod2_ci[1], mod2_ci[2],
                                              c_slope, mod1_ci[2,1], mod1_ci[2,2],
                                              o_e, o_e_u, o_e_l, bs, disc, gci[1],gci[3]))
}

glm_val_results <- rbind(glm_val_results, glm_val_alldat)
glm_val_results <- as.data.frame(glm_val_results)
names(glm_val_results)[9] <- "2"
names(glm_val_results)[10] <- "3"
old_names <- names(glm_val_results)
new_names <- c("year", "citl", "citl_li", "citl_ui", "c_slope", "c_slope_li", "c_slope_ui",
               "o_e", "o_e_ui", "o_e_li", "brier", "disc", "disc_li", "disc_ui")
glm_val_results <- glm_val_results %>% data.table::setnames(old = old_names, new = new_names)


gudf <- glm_val_results[-9,]
gudf[] <- lapply(gudf, function(x) as.numeric(as.character(x)))
gudf$year <- gudf$year + 0.1
gudf$mod <- "glm_up"


#######
# VC mod
#######


x <- df[vars]
dat <- na.omit(x)
tm <- (dat$year-2009)*12+dat$month
y <- dat$outcome
x1 <- as.matrix(sapply(dat[,1:(dim(dat)[2]-3)], as.numeric))

# Run model
fit <- glm(y ~ x1 + tm, family = binomial(link = logit))
summary(fit)
#save coefs
fit_coef <- t(fit$coefficients)

# add intercept term to val data
x <- df_val[vars]
dat <- na.omit(x)
#tm1 <- (dat$year-2012)*12+dat$month

tm1 <- rep(36, length(dat$year))
y <- dat$outcome
x1 <- as.matrix(sapply(dat[,1:(dim(dat)[2]-3)], as.numeric))
x1 <- cbind(1,x1,tm1)

# Multiply each row of data by coef vector
p_glm <- mapply(FUN = `*`, as.data.frame(x1), fit_coef)

# Sum over rows
glm_pred <- as.data.frame(rowSums(p_glm))
names(glm_pred)[1] <- "prob"

# Transfor to risk
glm_pred$risk <- exp(glm_pred$prob)/(1+exp(glm_pred$prob))

glm_pred <- cbind(glm_pred,dat$outcome,dat$year,dat$month, dat$year+dat$month/12-1/12,
                  (dat$year-2012)*12+dat$month)

names(glm_pred)[3] <- "outcome"
names(glm_pred)[4] <- "year"
names(glm_pred)[5] <- "month"
names(glm_pred)[6] <- "time"
names(glm_pred)[7] <- "time_m"


# Check response is between 0 and 1
summary(glm_pred$risk) # looks correct


####################
#### Validation ####
####################

# Group data by decile of predicted risk
glm_pred$grp100 <- ntile(glm_pred$risk, 100)

m100<- group_by(glm_pred, grp100) %>%
  summarise(m=mean(risk), pr = m*100)
m100 <- as.data.frame(m100)

ob_r100 <- group_by(glm_pred, grp100) %>%
  summarise(table(outcome)[2], or= ifelse(is.na(table(outcome)[2])==TRUE,0,table(outcome)[2]/length(outcome)*100))
t100_m <- cbind(ob_r100, m100) 

## Create calibration plots for the 100 quantiles
t100_all <- t100_m[-1]
cal_all <- ggplot(t100_all, aes(x=pr, y=or)) + geom_point(aes(y = or, x = pr)) +
  geom_abline(aes(intercept=0, slope=1)) + xlab("Predicted risk") + ylab("Observed risk")
cal_all


# c-slope
mod1 <- glm(outcome ~ prob, family="binomial",x=TRUE,y=TRUE, data= glm_pred)
#mod1 # 
mod1_ci <- confint(mod1)

# CITL
mod2 <- glm(outcome~offset(prob),family="binomial", data = glm_pred)
mod2_ci <- confint(mod2)

# Discrimination
g <- roc(outcome ~ risk, data = glm_pred, direction = "<")
plot(g)
g
g$auc # AUC is 0.9073
gci <- ci(roc(outcome ~ risk, data = glm_pred, direction = "<"))

# E/O
O <- prop.table(table(glm_pred$outcome))[2]*length(glm_pred$outcome)
E <- mean(glm_pred$risk)*length(glm_pred$outcome)
o_e_u <- ((sqrt(O)+1.96*0.5)^2)/E
o_e_l <- ((sqrt(O)-1.96*0.5)^2)/E

bs <- 1/length(glm_pred[,1]) * sum((glm_pred$risk - glm_pred$outcome)^2) 

glm_val_alldat <- c("12-18", mod2$coefficients, mod2_ci[1],
                    mod2_ci[2], mod1$coefficients[2],
                    mod1_ci[2,1], mod1_ci[2,2],
                    O/E, o_e_u, o_e_l, bs, g$auc, gci[1],gci[3])

# Loop generating validation results for each year separately 
glm_val_results <- c()
for(i in 2012:2019){
  dat <- glm_pred %>% filter(year==i)
  mod1 <- glm(outcome ~ prob, family="binomial",x=TRUE,y=TRUE, data= dat)
  c_slope <- mod1$coefficients[2]
  mod1_ci <- confint(mod1)
  mod2 <- glm(outcome~offset(prob),family="binomial", data = dat)
  citl <- mod2$coefficients
  mod2_ci <- confint(mod2)
  g <- roc(outcome ~ risk, data = dat, direction = "<")
  disc <- g$auc
  gci <- ci(roc(outcome ~ risk, data = dat, direction = "<"))
  O <- prop.table(table(dat$outcome))[2]*length(dat$outcome)
  E <- mean(dat$risk)*length(dat$outcome)
  o_e <- O/E
  o_e_u <- ((sqrt(O)+1.96*0.5)^2)/E
  o_e_l <- ((sqrt(O)-1.96*0.5)^2)/E
  bs <- 1/length(dat[,1]) * sum((dat$risk - dat$outcome)^2) 
  glm_val_results <- rbind(glm_val_results, c(i,citl, mod2_ci[1], mod2_ci[2],
                                              c_slope, mod1_ci[2,1], mod1_ci[2,2],
                                              o_e, o_e_u, o_e_l, bs, disc, gci[1],gci[3]))
}

glm_val_results <- rbind(glm_val_results, glm_val_alldat)
glm_val_results <- as.data.frame(glm_val_results)
names(glm_val_results)[9] <- "2"
names(glm_val_results)[10] <- "3"
old_names <- names(glm_val_results)
new_names <- c("year", "citl", "citl_li", "citl_ui", "c_slope", "c_slope_li", "c_slope_ui",
               "o_e", "o_e_ui", "o_e_li", "brier", "disc", "disc_li", "disc_ui")
glm_val_results <- glm_val_results %>% data.table::setnames(old = old_names, new = new_names)
# view(glm_val_results)


vcdf <- glm_val_results[-9,]
vcdf[] <- lapply(vcdf, function(x) as.numeric(as.character(x)))
vcdf$year <- vcdf$year + 0.2
vcdf$mod <- "VC"



##########
# DM Validation #
##########
### Start with DM with forgetting


# load coef data
load(file = "b_coefs_DM_3021_0.9997.rda")

dm <- as.data.frame(betas_dm)
dim(dm)
names(dm)
# drop last coef obs
dm <- dm[-(dim(dm)[1]),]

dm_val <- dm[(dim(dm)[1]-dim(df_val)[1]+1):(dim(dm)[1]),-(dim(dm)[2]-1)]
# DM val data is now same dimention as dm betas

# generate val data with variables in DM and add intercept term to val data
x <- df_val[vars]
x <- x %>% replace(is.na(.), 0)
dat <- na.omit(x)
y <- dat$outcome
x1 <- as.data.frame(sapply(dat[,1:(dim(dat)[2]-3)], as.numeric))
x1 <- cbind(1,x1)


# dm_val and val data now match.
# multiply and then row sum to get LP
dm_val <- dm_val * x1

# Sum over rows
dm_val <- as.data.frame(rowSums(dm_val))
names(dm_val)[1] <- "prob"

# Transfor to risk
dm_val$risk <- exp(dm_val$prob)/(1+exp(dm_val$prob))

dm_val <- cbind(dm_val,dat$outcome,dat$year,dat$month, dat$year+dat$month/12-1/12,
                (dat$year-2012)*12+dat$month)

names(dm_val)[3] <- "outcome"
names(dm_val)[4] <- "year"
names(dm_val)[5] <- "month"
names(dm_val)[6] <- "time"
names(dm_val)[7] <- "time_m"


# Check response is between 0 and 1
summary(dm_val$risk) # looks correct


dm_val <- dm_val %>% filter(risk!="NaN")

# c-slope
mod1 <- glm(outcome ~ prob, family="binomial",x=TRUE,y=TRUE, data= dm_val)
mod1_ci <- confint(mod1)
#mod1 # 

# CITL
mod2 <- glm(outcome~offset(prob),family="binomial", data = dm_val)
mod2_ci <- confint(mod2)

# Discrimination
g <- roc(outcome ~ risk, data = dm_val, direction = "<")
plot(g)
g$auc # AUC is 0.9073
gci <- ci(roc(outcome ~ risk, data = dm_val, direction = "<"))

# E/O
prop.table(table(dm_val$outcome))
O <- prop.table(table(glm_pred$outcome))[2]*length(glm_pred$outcome)
E <- mean(glm_pred$risk)*length(glm_pred$outcome)
o_e_u <- ((sqrt(O)+1.96*0.5)^2)/E
o_e_l <- ((sqrt(O)-1.96*0.5)^2)/E


bs <- 1/length(dm_val[,1]) * sum((dm_val$risk - dm_val$outcome)^2) 

dm_val_alldat <- c("12-18", mod2$coefficients, mod2_ci[1],
                   mod2_ci[2], mod1$coefficients[2],
                   mod1_ci[2,1], mod1_ci[2,2],
                   O/E, o_e_u, o_e_l, bs, g$auc, gci[1],gci[3])

# Loop generating validation results for each year separately 
dm_val_results <- c()
for(i in 2012:2019){
  dat <- dm_val %>% filter(year==i)
  mod1 <- glm(outcome ~ prob, family="binomial",x=TRUE,y=TRUE, data= dat)
  c_slope <- mod1$coefficients[2]
  mod1_ci <- confint(mod1)
  mod2 <- glm(outcome~offset(prob),family="binomial", data = dat)
  citl <- mod2$coefficients
  mod2_ci <- confint(mod2)
  g <- roc(outcome ~ risk, data = dat, direction = "<")
  disc <- g$auc
  gci <- ci(roc(outcome ~ risk, data = dat, direction = "<"))
  O <- prop.table(table(dat$outcome))[2]*length(dat$outcome)
  E <- mean(dat$risk)*length(dat$outcome)
  o_e <- O/E
  o_e_u <- ((sqrt(O)+1.96*0.5)^2)/E
  o_e_l <- ((sqrt(O)-1.96*0.5)^2)/E
  bs <- 1/length(dat[,1]) * sum((dat$risk - dat$outcome)^2) 
  dm_val_results <- rbind(dm_val_results, c(i,citl, mod2_ci[1], mod2_ci[2],
                                            c_slope, mod1_ci[2,1], mod1_ci[2,2],
                                            o_e, o_e_u, o_e_l, bs, disc, gci[1],gci[3]))
}

dm_val_results <- rbind(dm_val_results, dm_val_alldat)
dm_val_results <- as.data.frame(dm_val_results)
names(dm_val_results)[9] <- "2"
names(dm_val_results)[10] <- "3"
old_names <- names(dm_val_results)
new_names <- c("year", "citl", "citl_li", "citl_ui", "c_slope", "c_slope_li", "c_slope_ui",
               "o_e", "o_e_ui", "o_e_li", "brier", "disc", "disc_li", "disc_ui")
dm_val_results <- dm_val_results %>% data.table::setnames(old = old_names, new = new_names)
# view(dm_val_results)

ddf <- dm_val_results[-9,]
ddf[] <- lapply(ddf, function(x) as.numeric(as.character(x)))
ddf$year <- ddf$year -0.1
ddf$mod <- "dm"

#### Graph results ##########################################################

df <- rbind(gdf,gudf,ddf, vcdf)

ggplot(df, aes(x=year, y=citl, color=mod)) + geom_point(shape=3,size=1) +
  geom_segment(aes(x = year, y = citl_li, xend = year, yend = citl_ui, color=mod), data = df) + 
  theme_light() + theme(legend.position = "bottom", legend.title = element_blank())

df_long <- gather(df, pm, pm_num ,c(citl, c_slope, disc, o_e), factor_key=TRUE)
lower <- gather(df, pm, lower ,c(citl_li, c_slope_li, disc_li, o_e_li), factor_key=TRUE)
upper <- gather(df, pm, upper ,c(citl_ui, c_slope_ui, disc_ui, o_e_ui), factor_key=TRUE)
df_long$lower <- lower$lower
df_long$upper <- upper$upper


# Plot CITL, c_slope and discrimination
g <- ggplot(df_long, aes(x=year, y=pm_num, color=mod)) + geom_point(shape=3,size=1) +
  geom_segment(aes(x = year, y = lower, xend = year, yend = upper, color=mod), data = df_long)

g.labs <- c("Calibration-in-the-large", "Calibration slope", "Discrimination", "Observed-expected ratio")
names(g.labs) <- c("citl","c_slope", "disc", "o_e")

g + facet_wrap(~ pm, labeller = labeller(pm = g.labs), scales = "free")  + 
  scale_color_manual(labels = c("Bayesian model", "Logistic model", "yearly updated", "Varying coefficient model"),
                     values = c("royalblue3", "green3", "coral1", "mediumorchid3")) +
  labs(x = "Year", y = "Performance value", color = "Model") + theme_light() + 
  scale_x_continuous(breaks=seq(2012,2019,1), labels=seq(2012,2019,1)) + geom_hline()

g.labs <- c("Calibration-in-the-large", "Calibration slope", "Discrimination", "Observed-expected ratio")
names(g.labs) <- c("citl","c_slope", "disc", "o_e")

g + facet_wrap(~ pm, labeller = labeller(pm = g.labs), scales = "free")  + 
  scale_color_manual(labels = c("Bayesian model", "Logistic model", "Yearly updated model", "Varying coefficient model"),
                     values = c("royalblue3", "green3",  "grey45", "coral1")) +
  labs(x = "Year", y = "Performance value", color = "Model") + theme_light() + 
  scale_x_continuous(breaks=seq(2012,2019,1), labels=seq(2012,2019,1)) + theme(legend.position = "bottom")

g <- ggplot(df_long, aes(x=year, y=pm_num, color=mod)) + geom_line(size=0.5) +
  geom_line(aes(x = year, y = lower, color=mod), data = df_long) +
  geom_line(aes(x = year, y = upper, color=mod), data = df_long)

g + facet_wrap(~ pm, labeller = label_context, scales = "free")  + 
  scale_color_manual(labels = c("DM", "GLM", "GLM_UP", "VC"),
                     values = c("royalblue3", "green3", "coral1", "mediumorchid1")) +
  labs(x = "Year", y = "Performance value", color = "Model") + theme_light( )+ 
  scale_x_continuous(breaks=seq(0,75,12), labels=seq(2012,2018,1))


df_long <- gather(df, pm, pm_num ,c(citl, c_slope, disc,o_e, brier), factor_key=TRUE)

g <- ggplot(df_long, aes(x=year, y=pm_num, color=mod)) + geom_line(size=0.5) 

g + facet_wrap(~ pm, labeller = label_context, scales = "free")  + 
  scale_color_manual(labels = c("DM", "GLM", "GLM_UP", "VC"),
                     values = c("royalblue3", "green3", "coral1", "mediumorchid1")) +
  labs(x = "Year", y = "Performance value", color = "Model") + theme_light( )+ 
  scale_x_continuous(breaks=seq(0,75,12), labels=seq(2012,2018,1))




####
df <- rbind(gdf,gudf,ddf)

df_long <- gather(df, pm, pm_num ,c(citl, c_slope, disc, o_e), factor_key=TRUE)
lower <- gather(df, pm, lower ,c(citl_li, c_slope_li, disc_li, o_e_li), factor_key=TRUE)
upper <- gather(df, pm, upper ,c(citl_ui, c_slope_ui, disc_ui, o_e_ui), factor_key=TRUE)
df_long$lower <- lower$lower
df_long$upper <- upper$upper


# Plot CITL, c_slope and discrimination
g <- ggplot(df_long, aes(x=year, y=pm_num, color=mod)) + geom_point(shape=3,size=1) +
  geom_segment(aes(x = year, y = lower, xend = year, yend = upper, color=mod), data = df_long)

g.labs <- c("Calibration-in-the-large", "Calibration slope", "Discrimination", "Observed-expected ratio")
names(g.labs) <- c("citl","c_slope", "disc", "o_e")

g + facet_wrap(~ pm, labeller = labeller(pm = g.labs), scales = "free")  + 
  scale_color_manual(labels = c("DM", "GLM", "GLM_UP"),
                     values = c("royalblue3", "green3", "coral1")) +
  labs(x = "Year", y = "Performance value", color = "Model") + theme_light() + 
  scale_x_continuous(breaks=seq(2012,2018,1), labels=seq(2012,2018,1))

g.labs <- c("Calibration-in-the-large", "Calibration slope", "Discrimination", "Observed-expected ratio")
names(g.labs) <- c("citl","c_slope", "disc", "o_e")

g + facet_wrap(~ pm, labeller = labeller(pm = g.labs), scales = "free")  + 
  scale_color_manual(labels = c("Dynamic model", "GLM", "GLM yearly update"),
                     values = c("royalblue3", "green3", "coral1")) +
  labs(x = "Year", y = "Performance value", color = "Model") + theme_light() + 
  scale_x_continuous(breaks=seq(2012,2018,1), labels=seq(2012,2018,1))
#
