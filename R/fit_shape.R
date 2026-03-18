## fitting and comparing the fits of different shapes to the relation

# surv_rate  is a Boolean to select whether the relations (sigmoidal and quadratic)
# will be fitted on the scale of the predictor (thi sis the case for survival where probabilities are
# responses; and when it is FALSE (default), then the response is used without any link function, i.e.
# the version of the sigmoid will be used that scales the values back to -2 and 2 range
# and the quadratic will be implemented to response directly

##  a function to assess the fit of diff. shapes to the data
fit_shape <- function(data = all_trans,
                      x = 'Temp_z',
                      y = 'Trait_z',
                      ID = 1,
                      Thresh = 2,
                      surv_rate = FALSE,
                      out_folder = './output/output_nonL/shapes/'){
  sub_data <- droplevels(data[data$ID ==  ID, ])
sub_data <- sub_data[! (is.na(sub_data[[x]]) | is.na(sub_data[[y]])), ]

  lin_mod_form <- paste0(y, '~ interc + beta *', x)

  if(surv_rate){
    sigm_mod_form <- paste0(y, '~ 1/(1+exp((interc + beta*', x, ')))')
    quad_mod_form <- paste0(y, '~ 1/(1+exp((interc + beta *', x, '+ beta2 * ', x, '^2)))')
  }else{
    sigm_mod_form <- paste0(y, '~ (1/(1+exp(-5*(interc + beta*', x, '))) - 0.5)*4')
    quad_mod_form <- paste0(y, '~ interc + beta *', x, '+ beta2 * ', x, '^2')
  }
  null_mod_form <- paste0(y, '~ interc')

  ## have to catch the errors while fitting

  tt.error.null <- tryCatch(NullMod <- nlsLM(null_mod_form,
                                              start = list(interc = 1),# upper = c(10),
                                              data = sub_data, algorithm  = "LM",control = list(maxiter = 200)),
                              error=function(e) e)
  if(is(tt.error.null,"error")){
    warning(cat('error when fitting null model \n',
                tt.error.null[1]$message, 'for study',
                unique(sub_data$ID), '\n'))
    for (i in c(-10:10)){
        tt.error.null_test <- tryCatch(null_test <- nlsLM(null_mod_form,
                                                               # start = list(interc = i),
                                                                data = sub_data, algorithm  = "LM",control = list(maxiter = 200)),
                                          error=function(e) e)
        if(is(tt.error.null_test,"error")){
          warning(cat('error when fitting null model \n',
                      tt.error.null_test[1]$message, 'for study',
                      unique(sub_data$ID), '\n'))
        }
        if(! is(tt.error.null_test,"error")){
          message("successfully fitted null model")
          tt.error.null <- tt.error.null_test
          NullMod <- null_test
        }
    }
  } ## trying to find the fitting init par-rs for those fits when sigmoid starts with the value

  tt.error.linRel <- tryCatch(linRel <- nlsLM(lin_mod_form,
                                              start = list(interc = 0, beta = 1), upper = c(10, 10),
                                              data = sub_data, algorithm  = "LM",control = list(maxiter = 200)),
                              error=function(e) e)
  if(is(tt.error.linRel,"error")){
    warning(cat('error when fitting linear relation \n',
                tt.error.linRel[1]$message, 'for study',
                unique(sub_data$ID), '\n'))
  }

  if(! is(tt.error.linRel,"error")){
    message("successfully fitted linear relation")
    if(is.infinite(MuMIn::AICc(linRel))){
      warning(cat('infinite AIC for linear relation \n'))
    }
  }
  tt.error.quadRel <- tryCatch(quadRel <- nlsLM(quad_mod_form,
                                                start = list(interc = 0, beta = 1, beta2 = 1), upper = c(10, 10, 10),
                                                data = sub_data, algorithm  = "LM",control = list(maxiter = 200)),
                               error=function(e) e)
  if(is(tt.error.quadRel,"error")){
    warning(cat('error when fitting quadratic relation \n',
                tt.error.quadRel[1]$message, 'for study',
                unique(sub_data$ID), '\n'))
    for (i in c(-15:15)){
      for (j in c(-2:2)){
          for (k in c(-2:2)){
        tt.error.quadRel_test <- tryCatch(quadRel_test <- nlsLM(quad_mod_form,
                                                                start = list(interc = i, beta = j, beta2 = k),
                                                                data = sub_data, algorithm  = "LM",control = list(maxiter = 200)),
                                          error=function(e) e)
        if(is(tt.error.quadRel_test,"error")){
          warning(cat('error when fitting sigmoid relation \n',
                      tt.error.quadRel_test[1]$message, 'for study',
                      unique(sub_data$ID), '\n'))
        }
        if(! is(tt.error.quadRel_test,"error")){
          message("successfully fitted quadratic relation")
          tt.error.quadRel <- tt.error.quadRel_test
          quadRel <- quadRel_test
        }
      }
      }
    }
  } ## trying to find the fitting init par-rs for those fits when quadratic starts with the values that lead to a singular gradient matrix

  if(! is(tt.error.quadRel,"error")) {
    message("successfully fitted quadratic relation")
    if(is.infinite(MuMIn::AICc(quadRel))){
      warning(cat('infinite AIC for quadratic relation \n'))
    }
  }

  tt.error.sigmRel <- tryCatch(sigmRel <- nlsLM(sigm_mod_form,
                                                start = list(interc = 0, beta = 1),
                                                data = sub_data, algorithm  = "LM",control = list(maxiter = 200)),
                               error=function(e) e)
  if(is(tt.error.sigmRel,"error")){
    warning(cat('error when fitting sigmoid relation \n',
                tt.error.sigmRel[1]$message, 'for study',
                unique(sub_data$ID), '\n'))
    for (i in c(-15:15)){
      for (j in c(-2:2)){
        tt.error.sigmRel_test <- tryCatch(sigmRel_test <- nlsLM(sigm_mod_form,
                                                                start = list(interc = i, beta = j),
                                                                data = sub_data, algorithm  = "LM",control = list(maxiter = 200)),
                                          error=function(e) e)
        if(is(tt.error.sigmRel_test,"error")){
          warning(cat('error when fitting sigmoid relation \n',
                      tt.error.sigmRel[1]$message, 'for study',
                      unique(sub_data$ID), '\n'))
        }
        if(! is(tt.error.sigmRel_test,"error")){
          message("successfully fitted sigmoid relation")
          tt.error.sigmRel <- tt.error.sigmRel_test
          sigmRel <- sigmRel_test
        }
      }
    }
  } ## trying to find the fitting init par-rs for those fits when sigmoid starts with the values that lead to a singular gradient matrix

  if(! is(tt.error.sigmRel,"error")) {
    message("successfully fitted sigmoid relation")
    if(is.infinite(MuMIn::AICc(sigmRel))){
      warning(cat('infinite AIC for sigmoid relation \n'))
    }
  }

  ## saving a plot as a pdf
  if(! (is(tt.error.linRel,"error") | is(tt.error.quadRel,"error") |
        is(tt.error.sigmRel,"error") )){
    pdf(paste0(out_folder, ID, '_', unique(sub_data$Study_Authors),
               '_', unique(sub_data$Species), '_', x,'_', y, '.pdf'))
    plot(sub_data[[x]], sub_data[[y]],
         ylab = y,
         xlab = x)
    points(sub_data[[x]], fitted(linRel), pch = 19, col = 'red')
    points(sub_data[[x]], fitted(sigmRel), pch = 21, col = 'grey')
    points(sub_data[[x]], fitted(quadRel), pch = 19, col = 'blue')
    legend('bottomright', c('linear', 'sigmoid', 'quadratic'),
           col = c('red', 'grey', 'blue'), pch = c(19, 21, 19))
    dev.off()

  }

 # extracting the fits
  if(! (is(tt.error.linRel,"error") | is(tt.error.quadRel,"error") |
        is(tt.error.sigmRel,"error") )){
    AIC_lin <- MuMIn::AICc(linRel)
    AIC_quad <- MuMIn::AICc(quadRel)
    AIC_sigm <- MuMIn::AICc(sigmRel)
    AIC_null <- MuMIn::AICc(NullMod)
    if(! (is.infinite(AIC_lin) | is.infinite(AIC_quad) | is.infinite(AIC_sigm) | is.infinite(AIC_null))){

      res <- data.frame('ID' = unique(sub_data$ID),
                        'AIC_null' = AIC_null,
                        'AIC_Lin' = AIC_lin,
                        'AIC_Quad' = AIC_quad,
                        'AIC_Sigm' = AIC_sigm,
                        'Int_Lin' = coef(linRel)[names(coef(linRel)) == 'interc'],
                        'Beta_Lin' = coef(linRel)[names(coef(linRel)) == 'beta'],
                        'Int_Quad' = coef(quadRel)[names(coef(quadRel)) == 'interc'],
                        'Beta_Quad' = coef(quadRel)[names(coef(quadRel)) == 'beta'],
                        'Beta2_Quad' = coef(quadRel)[names(coef(quadRel)) == 'beta2'],
                        'Int_Sigm' = coef(sigmRel)[names(coef(sigmRel)) == 'interc'],
                        'Beta_Sigm' = coef(sigmRel)[names(coef(sigmRel)) == 'beta'])
      vect_AIC <- c(AIC_lin, AIC_quad, AIC_sigm)
      names(vect_AIC) <- c('linear', 'quadratic', 'sigmoid')
      ord <- order(c(AIC_lin, AIC_quad, AIC_sigm))
      min_AIC <- min(vect_AIC)
      Delta1 <- vect_AIC[ord][2] - min_AIC
      Delta2 <- vect_AIC[ord][3] - min_AIC


      res$minAIC <- min_AIC
      res$mod_minAIC <- names(which(vect_AIC == min_AIC))
      res$Delta1_AIC <- Delta1
      res$Delta2_AIC <- Delta2
      res$Delta_Best_null <- min_AIC - AIC_null
      if(Delta1 > Thresh){res$Selected <- names(vect_AIC[ord][1])
      } else {
        if(Delta1 <= Thresh & Delta2 <= Thresh){
          if(length(intersect(names(vect_AIC[ord][c(1,2,3)]),
                              c('linear', 'sigmoid', 'quadratic'))) == 3){
            res$Selected <- 'linear/sigmoid'
          }
        } else {
          if(Delta1 <= Thresh & Delta2 > Thresh){
            if(length(intersect(names(vect_AIC[ord][c(1,2)]),
                                c('linear', 'sigmoid'))) == 2){
              res$Selected <- 'linear/sigmoid'
            }
            if(length(intersect(names(vect_AIC[ord][c(1,2)]),
                                c('linear', 'quadratic'))) == 2){
              res$Selected <- 'linear'
            }
            if(length(intersect(names(vect_AIC[ord][c(1,2)]),
                                c('sigmoid', 'quadratic'))) == 2){
              res$Selected <- 'sigmoid'
            }
          }
        }
      }
      return(res)
    } else {
      message(paste('infinite AIC for at least one of the models'))
    }
  } else {
    message(paste('error in convergence of at least one of the models'))
  }
}
