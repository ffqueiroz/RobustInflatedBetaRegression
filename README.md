# RobustInflatedBetaRegression

Function to fit the robust estimators for the zero-or-one inflated beta regression model.
The function used to fit is <tt>fit.betainflated()</tt>.

```{r usage functions, eval=FALSE}
fit.betainflated(y, S, X, Z, alphaoptimal = TRUE, alpha0 = 0, weights = FALSE,
                 linkmu="logit", linkphi="log", linkvartheta = "logit")
```
