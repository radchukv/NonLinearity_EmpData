
## first exploring the behaviour of sigmoid function as specified currently, i.e. to give as output rnage of values between -2 and 2 (i.e. if limate affects traits)
fun_sigm <- function(beta, Trait_value, interc){
  Dem_rate <- (1/(1+exp(-5*(interc + beta*Trait_value))) - 0.5)*4  ## check if sd_temp is needed?
  return(Dem_rate)
}

fun_sigm_cl <- function(beta, Trait_value, interc){
  Dem_rate <- (1/(1+exp(-(interc + beta*Trait_value))))   ## check if sd_temp is needed?
  return(Dem_rate)
}

fun_sigm_cl(beta = 1, Trait_value = 0, interc = 0)
## adjusting the explored par-rs to highlight the range tested in the model and the one
## identified in the empirical sTrait
Trait_v <- seq(-2, 2, by = 0.1)

interc <- seq(-5, 5, by = 1)

beta <- c(seq(1, 1, by = 1))

DemRate <- unlist(lapply(Trait_v, FUN = function(x){fun_sigm_cl(beta = 1, interc = 0, Trait_value = x)}))
plot(Trait_v, DemRate)


## and now trying to apply the function across diff. intercept and beta values
dat <- data.frame(interc = numeric(0), demRate = numeric(0), beta = numeric(0),
                  Trait_v = numeric(0))
for(i in interc){
  for(j in beta){
    newDat <- data.frame(interc = rep(i, length(Trait_v)), beta = rep(j, length(Trait_v)),
                         Trait_v = Trait_v)
    newDat$demRate <- unlist(lapply(Trait_v, FUN = function(x){fun_sigm_cl(beta = j, interc = i, Trait_value = x)}))
    dat <- rbind(dat, newDat)
  }
}


subs <- subset(dat, interc == 1)
## to colour the panels reflecting the par-r values inferred from the empiricla data
## add a df
dat_col <- data.frame(interc = rep(c(-1, 0, 1, 2), 2),
                      beta = rep(c(-0.1, 0.1),
                                 each = length(c(-1, 0, 1, 2))),
                      col = 'blue', demRate = 0, Trait_v = 0)

pdf('./output/output_nonL/SupplFig_ParSpace_Sigmoid_betaRows_IntercColumns_classicFormula.pdf', height = 8, width = 6)
ggplot(dat, aes(x = Trait_v, y = demRate)) +
  geom_line() +
  # geom_rect(data = dat_col, aes(fill = col),
  #           xmin = -Inf, xmax = Inf,
  #           ymin = -Inf, ymax = Inf,
  #           fill = 'grey', alpha = 0.4) +
  facet_wrap(interc) +
  theme_bw() + ylab('Demographic rate') +
  xlab('Trait') +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white'))

dev.off()

## so, the function does behave well, but under low betas it resembles linear.. and that is
## what we get with the data
ggplot(subs, aes(x = Trait_v, y = demRate)) +
  geom_line()




## For probs of surivval we in fact do not need these scaling (i.e. multiple by 5, -0.5) because we end up with probabilities, and they rnage between 0 and 1
fun_sigm <- function(beta, Trait_value, interc){
  Dem_rate <- (1/(1+exp(-(interc + beta*Trait_value))))
  return(Dem_rate)
}

fun_sigm(beta = 1, Trait_value = 0, interc = 0)
## adjusting the explored par-rs to highlight the range tested in the model and the one
## identified in the empirical sTrait
Trait_v <- seq(-2, 2, by = 0.1)

interc <- seq(0, 5, by = 1)

beta <- c(seq(-5, -1, by = 1), -0.1, 0, 0.1, seq(1, 5, by = 1))

DemRate <- unlist(lapply(Trait_v, FUN = function(x){fun_sigm(beta = 1, interc = 0, Trait_value = x)}))
plot(Trait_v, DemRate)


## and now trying to apply the function across diff. intercept and beta values
dat <- data.frame(interc = numeric(0), demRate = numeric(0), beta = numeric(0),
                  Trait_v = numeric(0))
for(i in interc){
  for(j in beta){
    newDat <- data.frame(interc = rep(i, length(Trait_v)), beta = rep(j, length(Trait_v)),
                         Trait_v = Trait_v)
    newDat$demRate <- unlist(lapply(Trait_v, FUN = function(x){fun_sigm(beta = j, interc = i, Trait_value = x)}))
    dat <- rbind(dat, newDat)
  }
}

## to colour the panels reflecting the par-r values inferred from the empiricla data
## add a df
dat_col <- data.frame(interc = rep(c(-1, 0, 1, 2), 2),
                      beta = rep(c(-0.1, 0.1),
                                 each = length(c(-1, 0, 1, 2))),
                      col = 'blue', demRate = 0, Trait_v = 0)

pdf('./output/output_nonL/SupplFig_ParSpace_Sigmoid_betaRows_IntercColumns_correctedFormula.pdf', height = 8, width = 6)
ggplot(dat, aes(x = Trait_v, y = demRate)) +
  geom_line() +
  geom_rect(data = dat_col, aes(fill = col),
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf,
            fill = 'grey', alpha = 0.4) +
  facet_grid(beta ~ interc) +
  theme_bw() + ylab('Demographic rate') +
  xlab('Trait') +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white'))

dev.off()


# testing an extension of function for fecundity


## For probs of surivval we in fact do not need these scaling (i.e. multiple by 5, -0.5) because we end up with probabilities, and they rnage between 0 and 1
fun_fecund <- function(beta, Trait_value, interc){
  Dem_rate <- (1/(1+exp(-5*(interc + beta*Trait_value))) - 0.5)*20
 if (Dem_rate < 0){
    Dem_rate = 0
  }
  return(Dem_rate)
}

# splitting the function by trait value (> or < 0)
fun_fecund1 <- function(beta, Trait_value, interc){
  if (Trait_value >= 0){
  Dem_rate <- (1/(1+exp(-5*(interc + beta*Trait_value))) - 0.5)*20 }
  else {
    Dem_rate <- (1/(1+exp(-5*(interc + beta*Trait_value))) - 0.5)*4
  }
  if (Dem_rate < 0){
    Dem_rate = 0
  }
  return(Dem_rate)
}


fun_fecund1(beta = 1, Trait_value = -0.5, interc = 0)
## adjusting the explored par-rs to highlight the range tested in the model and the one
## identified in the empirical sTrait
Trait_v <- seq(-2, 2, by = 0.1)

interc <- seq(0, 5, by = 1)

beta <- c(seq(-5, -1, by = 1), -0.1, 0, 0.1, seq(1, 5, by = 1))

DemRate <- unlist(lapply(Trait_v, FUN = function(x){fun_fecund1(beta = 1, interc = 0, Trait_value = x)}))
plot(Trait_v, DemRate)


## and now trying to apply the function across diff. intercept and beta values
dat <- data.frame(interc = numeric(0), demRate = numeric(0), beta = numeric(0),
                  Trait_v = numeric(0))
for(i in interc){
  for(j in beta){
    newDat <- data.frame(interc = rep(i, length(Trait_v)), beta = rep(j, length(Trait_v)),
                         Trait_v = Trait_v)
    newDat$demRate <- unlist(lapply(Trait_v, FUN = function(x){fun_fecund(beta = j, interc = i, Trait_value = x)}))
    dat <- rbind(dat, newDat)
  }
}

## to colour the panels reflecting the par-r values inferred from the empiricla data
## add a df
dat_col <- data.frame(interc = rep(c(-1, 0, 1, 2), 2),
                      beta = rep(c(-0.1, 0.1),
                                 each = length(c(-1, 0, 1, 2))),
                      col = 'blue', demRate = 0, Trait_v = 0)

pdf('./output/output_nonL/SupplFig_ParSpace_Sigmoid_betaRows_IntercColumns_Fecundity.pdf', height = 8, width = 6)
ggplot(dat, aes(x = Trait_v, y = demRate)) +
  geom_line() +
  geom_rect(data = dat_col, aes(fill = col),
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf,
            fill = 'grey', alpha = 0.4) +
  facet_grid(beta ~ interc) +
  theme_bw() + ylab('Demographic rate') +
  xlab('Trait') +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white'))

dev.off()


# now that I am at it, let me also test Poisson...
fun_pois <- function(beta, Trait_value, interc){
  Dem_rate <- exp(interc + beta*Trait_value)
  return(Dem_rate)
}

fun_pois(beta = 1, Trait_value = 0, interc = 0)
## adjusting the explored par-rs to highlight the range tested in the model and the one
## identified in the empirical sTrait
Trait_v <- seq(-2, 2, by = 0.1)

interc <- seq(-5, 5, by = 1)

beta <- c(seq(1, 1, by = 1))

DemRate <- unlist(lapply(Trait_v, FUN = function(x){fun_pois(beta = 1, interc = 0, Trait_value = x)}))
plot(Trait_v, DemRate)


## and now trying to apply the function across diff. intercept and beta values
dat <- data.frame(interc = numeric(0), demRate = numeric(0), beta = numeric(0),
                  Trait_v = numeric(0))

interc <- seq(1, 5, by = 1)

beta <- c(seq(-2, 2, by = 1))

for(i in interc){
  for(j in beta){
    newDat <- data.frame(interc = rep(i, length(Trait_v)), beta = rep(j, length(Trait_v)),
                         Trait_v = Trait_v)
    newDat$demRate <- unlist(lapply(Trait_v, FUN = function(x){fun_pois(beta = j, interc = i, Trait_value = x)}))
    dat <- rbind(dat, newDat)
  }
}

## to colour the panels reflecting the par-r values inferred from the empiricla data
# ## add a df
# dat_col <- data.frame(interc = rep(c(0:5), 2),
#                       beta = rep(c(-1.0, -0.5,  0.0,  0.5,  1.0),
#                                  each = length(c(0:5))),
#                       col = 'blue', demRate = 0, Trait_v = 0)

pdf('./output/output_nonL/Test_ParSpace_Poisson_betaRows_IntercColumns.pdf', height = 8, width = 6)
ggplot(dat, aes(x = Trait_v, y = demRate)) +
  geom_line() +
  # geom_rect(data = dat_col, aes(fill = col),
  #           xmin = -Inf, xmax = Inf,
  #           ymin = -Inf, ymax = Inf,
  #           fill = 'grey', alpha = 0.4) +
  facet_grid(beta ~ interc, scales = 'free') +
  theme_bw() + ylab('Fecundity') +
  xlab('Trait') +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white'))

dev.off()


## so, based on this test, I do not think I need to fit a Poisson -
# the relations observed from the data do not seem oto follow this relation


