################################################################################
#### Séries chronologiques TP3: fonctions utiles
#### LSTAT2170
################################################################################
# Comp.Sarima: Comparaison au moyen du critère AIC de différents mod?les S-ARIMA ajust?s ? une s?rie temporelle, avec ?limination de la tendance et de la saisonnalit? sans estimation (m?thode des diff?rentes it?r?es).
# Coef.p : Test de signification sur les coefficients du mod?le
# OneAhead : Somme des carr? des ?carts pour une pr?diction ?
# 	un pas de temps, pour 20% de la longueur de la s?rie.
################################################################################


################################################################################
Comp.Sarima = function(donnees, d, saison, D, p.max, q.max, P.max, Q.max){
#===========================================================================
# OBJET:
# Comparaison au moyen du crit?re AIC de diff?rents mod?les S-ARIMA ajust?s ? une s?rie temporelle, avec ?limination de la tendance et de la saisonnalit? sans estimation (m?thode des diff?rentes it?r?es).
#============================================================================
# SPECIFICATIONS
#----------------------------------------------------------------------------
# PRE:
# donnees : s?rie temporelle univari?e - objet-R  de type ts()
# d: nombres de diff?rentiations de lag 1, pour l'?limination de la tendance
# saison : p?riodicit? de la saisonnalit?
# D: nombre de diff?rentiations de lag 'saison' pour l'?limination de la saisonnalit?
# p.max, q.max, P.max, Q.max : ordres maximums consid?r?s pour l'ajustement du mod?le S-ARIMA de la forme mod?le (p,d,q)x(P,D,Q)_saison
#-----------------------------------------------------------------------------
# POST:
# Ecriture dans la fen?tre de commande des meilleurs mod?les s?lectionn?s sur base d'un crit?re AIC (le nombre de mod?les affich?s est ?gal ? 10% du nombre de mod?les test?s).
# Pour chacun de ces mod?les, on donne le nombre de param?tres impliqu?s et l'?cart du AIC correcpondant au AIC minimal atteint
#-----------------------------------------------------------------------------
# RMQ:
#(!) si p.max, q.max, P.max, Q.max sont grands, le temps de calcul peut ?tre important.
#=============================================================================
# (c): C. Timmermans, UCL - stat2170, 01/11/07
#============================================================================

# initialisation d'une table ? 4D reprennant les AIC des mod?les test?s
AIC.table <-  array(dim = c(p.max + 1, q.max + 1 , P.max + 1, Q.max + 1))

# construction de la table des AIC
for(i in 0:p.max){
  for(j in 0:q.max){
    for(k in 0:P.max){
      for(l in 0:Q.max){
          tmp <- arima(donnees, order = c(i, d, j), seasonal = list(order = c(k, D, l), period = saison))
          AIC.table[i + 1, j + 1, k + 1, l + 1] <- tmp$aic
      }
    }
  }
}

# Selection du minimum des AIC  et construction de AIC.table2 = AIC.table -min(AIC)  -- accroit la lisibilit? des r?sultats --
AICmin <- min(AIC.table)
AIC.table2 <- AIC.table - AICmin

# Quantile 10% de la distribution des AIC.table2
AICQ10 <- quantile(AIC.table2, 0.3, na.rm = TRUE)

#Pour les 10% meilleurs mod?les ajust?s, affichage des r?sultats.
for(i in 0:p.max){
  for(j in 0:q.max){
    for(k in 0:P.max){
      for(l in 0:Q.max){
        if(AIC.table2[i + 1, j + 1, k + 1, l + 1] < AICQ10){
            cat('modele (p,d,q)x(P,D,Q)_saison : ', i, d, j, 'x', k, D, l, '_', saison, ':  nb param: ', i+j+k+l,
                '    AIC:', AIC.table2[i+1,j+1,k+1,l+1], '\n')
        }
      }
    }
  }
}
}
################################################################################
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
################################################################################
## Coef.p()
## functions for time series course
# November 7th, 2006
# (c) Hilmar Boehm, Universite Catholique de Louvain, boehm@stat.ucl.ac.be
#--------------------------------------------------------------------------
# Significance test for coefficients...
# coeff.vect: vector of coefficients to be tested
# var.vect: vector of variances for coefficients
# output: vector of p-values for two-sided univariate tests
# H0: coef=0, H1: coef !=0
#---------------------------------------------------------------------------
coef.p <- function (coeff.vect, var.vect){
if (length(coeff.vect) != length(var.vect)) output = NaN
else {
	sigma.vect = sqrt(var.vect)
	n = length(sigma.vect)
	output <- sapply(1:n, function(i) 2*pnorm(abs(coeff.vect[i]) / sigma.vect[i], lower.tail = F))
	}
  output
}

#-------------------------------------------------------------------------------
################################################################################

################################################################################
#-------------------------------------------------------------------------------
# OneAhead()
## functions for time series course
# November 7th, 2006
# (c) Hilmar Boehm, Universite Catholique de Louvain, boehm@stat.ucl.ac.be
# modified: Catherine Timmermans
# December, 2006
#--------------------------------------------------------------------------
# compute one step prediction error for SARIMA time series
# ts: time series, order: vector (p,d,q)
# seasonal: list (order: vector (P,D,Q), period:numeric)
#--------------------------------------------------------------------------
OneAhead <- function(ts1, order, seasonal = list(order = c(0,0,0), period = 0)){
  n <- length(ts1)
  n80 <- floor(0.8 * n)
  n20 <- n - n80
  tmp <- numeric(n)
  for(i in n80:n){
    ts1.part <- ts1[1:(i-1)]
    tmp.model <- arima(ts1.part, order = order, seasonal = seasonal)
    tmp[i] <- predict(tmp.model, n.ahead = 1)$pred[1]
  }
  error <-  sum(((tmp - ts1)[n80:n])^2) / n20
  tspred <- c(ts1[1:(n80 - 1)], tmp[n80:n])
  return(list(tspred = tspred, error = error))
}
#-------------------------------------------------------------------------------
################################################################################

Fpe.ar <- function(serie, modele){
  modele$var.pred * (length(serie) + modele$order) / (length(serie) - modele$order)
}

Aic.ar <- function(serie, modele){
  log(modele$var.pred) + 2 * modele$order / length(serie)
}

Bic.ar <- function(serie,modele){
  log(modele$var.pred) + (modele$order) / length(serie) * log(length(serie))
}

Aic.arima <- function(serie, modele){
  log(modele$sigma2) + 2 * (length(modele$coef) - 1) / length(serie)
}

Bic.arima <- function(serie, modele){
  log(modele$sigma2) + (length(modele$coef) - 1) * log(length(serie)) / length(serie)
}




