
## first exploring the behaviour of sigmoid function as specified currently, i.e. to give as output rnage of values between -2 and 2 (i.e. if limate affects traits)
fun_sigm <- function(beta, Trait_value, interc){
  Dem_rate <- (1/(1+exp(-5*(interc + beta*Trait_value))) - 0.5)*4
  return(Dem_rate)
}

fun_sigm(beta = 1, Trait_value = 0, interc = 0)
## adjusting the explored par-rs to highlight the range tested in the model and the one
## identified in the empirical sTrait
Trait_v <- seq(-2, 2, by = 0.1)

interc <- seq(-2, 2, by = 1)

beta <- c(seq(-5, 5, by = 1), -0.1, 0.1)

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

pdf('./output/output_nonL/SupplFig_ParSpace_Sigmoid_betaRows_IntercColumns_innerForm.pdf', height = 8, width = 6)
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

## so, the function does behave well, but under low betas it resembles linear.. and that is
## what we get with the data



## For probs of surivval we in fact do not need these scaling (i.e. multiple by 5, -0.5) because we end up with probabilities, and they rnage between 0 and 1
fun_sigm <- function(beta, Trait_value, interc){
  Dem_rate <- (1/(1+exp(interc + beta*Trait_value)))
  return(Dem_rate)
}

fun_sigm(beta = 1, Trait_value = 0, interc = 0)
## adjusting the explored par-rs to highlight the range tested in the model and the one
## identified in the empirical sTrait
Trait_v <- seq(-2, 2, by = 0.1)

interc <- seq(-2, 2, by = 1)

beta <- seq(-5, 5, by = 1)

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


