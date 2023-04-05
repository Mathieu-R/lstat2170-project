select.sarima <- function(ts, d, D, max.p, max.q, max.P, max.Q, season) {
	sarima_table <- NULL

	for (i in 0:max.p) {
		for (j in 0:max.q) {
			for (k in 0:max.P) {
				for (l in 0:max.Q) {
					sarima_model <- Arima(
						ts, 
						order = c(i, d, j),
						seasonal = list(order = c(k, D, l), period = season),
						method = "ML"
					)

					sarima_table <- rbind(sarima_table, compute.sarima_summary(sarima_model, i, d, j, k, D, l, season))
				}
			}
		}
	}

	return(sarima_table)
}

compute.sarima_summary <- function(sarima_model, i, d, j, k, D, l, s) {
	sarima_summary <- data.frame(
		Order = paste("(", i, d, j, ") x (", k, D, l, ")_", s),
		AIC = round(sarima_model$aic, 3),
		BIC = round(sarima_model$bic, 3),
		LogLik = round(sarima_model$loglik, 3),
		Parameters = length(sarima_model$coef)
	)

	return(sarima_summary)
}

coefficients.significance <- function(sarima_model) {
	coefficients <- sarima_model$coef
	variance <- diag(sarima_model$var.coef)
	std_dev <- sqrt(variance)

	n <- length(coefficients)

	significance <- sapply(1:n, function(i) {
		2 * pnorm(abs(coefficients[i]) / std_dev[i], mean = 0, sd = 1, lower.tail = FALSE)
	})

	return(significance)
}

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

predictive.power <- function(ts, order, seasonal) {
	n <- length(ts)
	
	ts_attributes = tsp(ts)
	ts_start = ts_attributes[1]
	ts_end = ceiling(ts_attributes[2])
	ts_freq = ts_attributes[3]

	delta <- ts_end - ts_start
	delta20_year <- delta * 0.2
	delta20_month <- ts_freq * delta20_year

	first80 <- ((n - delta20_month) + 1)
	last20 <- n - first80

	model_pred <- numeric(n)

	for (i in first80:n) {
		ts_slice <- ts[1:(i - 1)]
		model_slice <- Arima(
			ts_slice,
			order = order,
			seasonal = seasonal,
			method = "ML"
		)

		model_pred[i] <- predict(model_slice, n.ahead = 1)$pred[1]
	}

	pred_mse <- sum(((model_pred - ts)[first80:n]^2)) / last20
	ts_pred <- ts(model_pred[first80:n], start = ts_end - delta20_year, frequency = ts_freq)

	return(list(
		data = data.frame(
			fitted = ts_pred,
			time = zoo::as.Date(ts_pred)
		),
		mse = pred_mse
	))
}