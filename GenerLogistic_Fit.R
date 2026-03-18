# exploration of Density Dependence (DD)  -West African Giraffes
library(tidyverse)
library(ggplot2)
library(ggtext)
library(dplyr)
#install.packages("easynls")
library(easynls)

dat <- read.csv('./data/WestAfrGir.csv')

#dat_GR <- dat %>%          #pipe operator/pipeline - passes the dat data frame through a series of transformations
  #mutate(N_lag = lag(N),     # mutate - adds new column N_lag = lag(N) (shifts N column down by 1 row), column GR1 and column GR
        #GR1 = N/N_lag,         # calculating growth rate lambda
         #GR = log(N/N_lag)) %>%     # taking ln to get r grwth rate (in R it is not ln, but log)
  #filter(! is.na(GR))                 # removing rows with NA value
#dat_GR

dat_gomp <- dat %>%
  mutate(N_lead = lead(N),
         logN = log(N),
         logN_lead = log(N_lead),
         r = log(N_lead / N)) %>%
  filter(!is.na(r))



# fit a simple linear model
# modDD <- lm(GR~ N_lag, data = dat_GR)
# plot(modDD)
# summary(modDD)
# modDD$coefficients  ## to know where we get the intercept and slope
# for our plot

# Fit the Gompertz model
summary(model_gomp)
plot(model_gomp)
model_gomp$coefficients

#logN - logN_lead ~ Rmax*(1 - (N_lead / K)^B)
# N ~ N_lead*exp((Rmax*(1 - (N_lead/K)^B)))
ricker_model <- nls(
  formula = N ~ N_lead*exp((Rmax*(1 - (N_lead/K)^B))),
  data = dat_gomp,
  start = list(Rmax = 1.2, K = 1000, B = 1),
  algorithm = "port", 
  lower = c(Rmax = 0, K = 300, B = 0),
  upper = c(Rmax = 5, K = 2000, B = 2),
  nls.control(maxiter = 2000)
)

summary(ricker_model)
fitted(ricker_model)

# explore the fit
newdat = data.frame(N_lead = seq(100, 2000, by = 100))
newdat$predN <- predict(gompertz_model, newdata = newdat)

newdat$r_pred <- log(newdat$N_lead / newdat$predN)
ggplot(newdat, aes(x = N_lead, y = r_pred)) +
  geom_point() +
  xlab('Population size (N<sub>t</sub>)') +
  ylab('Population growth rate, r') +
  theme_bw() +
  theme(axis.title.x = ggtext::element_markdown(),
        panel.grid.major = element_blank())
  


## some test with SSricker function
require(ggplot2)
set.seed(123)
x <- 1:30
y <- 30 * x * exp(-0.3 * x) + rnorm(30, 0, 0.25)
dat <- data.frame(x = x, y = y)
fit <- nls(y ~ SSricker(x, a, b), data = dat)
## plot
ggplot(data = dat, aes(x = x, y = y)) + 
  geom_point() + 
  geom_line(aes(y = fitted(fit)))

dat1 <- dat  %>% 
  mutate(r = log(x/y))

ggplot(aes(dat1), x = x, y = r) + geom_point()

gompertz_model <- nls(
  formula = N ~ (exp(Rmax)*N_lead)/(1 + exp(Rmax*(1 - (N_lead/K)^B))),
  data = dat_gomp,
  start = list(Rmax = 1.2, K = 1000, B = 1),
  algorithm = "port", 
  lower = c(Rmax = 0, K = 300, B = 0),
  upper = c(Rmax = 5, K = 2000, B = 2),
  nls.control(maxiter = 2000)
)

plot(gompertz_model)

# Extract coefficients
a <- coef(model_gomp)[1]
b <- -coef(model_gomp)[2]  # Gompertz slope is negative
r <- b
log_K <- a / b
K <- exp(log_K)

cat("Estimated intrinsic growth rate r =", round(r, 3), "\n")
cat("Estimated carrying capacity K =", round(K, 2), "\n")





dat <- read.csv('./data/WestAfrGir.csv')

gomp_model = nlsfit(dat, model = 10, start = c(a = 1200, b = 2, c = 0.1))
gomp_model

nlsplot(dat, model = 10, start = c(a = 1200, b = 2, c = 0.1),
        xlab = 'Abundance', ylab = 'GR', position = 1)