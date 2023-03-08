## fitting and comparing the fits of different shapes to the relation

##  a function to assess the fit of diff. shapes to the data
fit_shape <- function(data = all_trans,
                      x = 'Temp_z',
                      y = 'Trait_z',
                      ID = 1,
                      Thresh = 2,
                      out_folder = './output/output_nonL/shapes/'){
  sub_data <- droplevels(data[data$ID ==  ID, ])

  lin_mod_form <- paste0(y, '~ interc + beta *', x)
  quad_mod_form <- paste0(y, '~ interc + beta *', x, '+ beta2 * ', x, '^2')
  sigm_mod_form <- paste0(y, '~ (1/(1+exp(-5*(interc + beta*', x, '))) - 0.5)*4')
  ## have to catch the errors while fitting
  tt.error.linRel <- tryCatch(linRel <- nlsLM(lin_mod_form,
                                              start = list(interc = 0, beta = 1), upper = c(10, 10),
                                              data = sub_data, algorithm  = "LM",control = list(maxiter = 200)),
                              error=function(e) e)
  if(is(tt.error.linRel,"error")){
    warning(cat('error when fitting linear relation \n',
                tt.error.linRel[1]$message, '\n'))
  }

  if(! is(tt.error.linRel,"error")){
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
                tt.error.quadRel[1]$message, '\n'))
  }
  if(! is(tt.error.quadRel,"error")) {
    if(is.infinite(MuMIn::AICc(quadRel))){
      warning(cat('infinite AIC for quadratic relation \n'))
    }
  }

  tt.error.sigmRel <- tryCatch(sigmRel <- nlsLM(sigm_mod_form,
                                                start = list(interc = -5, beta = -5), upper = c(10, 10),
                                                data = sub_data, algorithm  = "LM",control = list(maxiter = 200)),
                               error=function(e) e)
  if(is(tt.error.sigmRel,"error")){
    warning(cat('error when fitting sigmoid relation \n',
                tt.error.sigmRel[1]$message, '\n'))
  }
  if(! is(tt.error.sigmRel,"error")) {
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

  if(! (is(tt.error.linRel,"error") | is(tt.error.quadRel,"error") |
        is(tt.error.sigmRel,"error") )){
    AIC_lin <- MuMIn::AICc(linRel)
    AIC_quad <- MuMIn::AICc(quadRel)
    AIC_sigm <- MuMIn::AICc(sigmRel)
    if(! (is.infinite(AIC_lin) | is.infinite(AIC_quad) | is.infinite(AIC_sigm))){

      res <- data.frame('ID' = unique(sub_data$ID),
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
      if(Delta1 > Thresh){res$Selected <- names(vect_AIC[ord][1])
      } else {
        if(Delta1 <= Thresh & Delta2 <= Thresh){
          # res$Selected <- 'linear'
          ## this still has to be adjusted accordingly
          if(length(intersect(names(vect_AIC[ord][c(1,2,3)]),
                              c('linear', 'sigmoid', 'quadratic'))) == 3){
            res$Selected <- 'linear/sigmoid'
          }
          # if(length(intersect(names(vect_AIC[ord][c(1,2)]),
          #                     c('linear', 'quadratic'))) == 2){
          #   res$Selected <- 'linear'
          # }
          # if(length(intersect(names(vect_AIC[ord][c(1,2)]),
          #                     c('sigmoid', 'quadratic'))) == 2){
          #   res$Selected <- 'sigmoid'
          # }
          #
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
