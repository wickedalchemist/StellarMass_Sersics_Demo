# Table of Contents
1. [Demo Summary](README.md#demo-summary)
 
# Demo Summary
In this snippet I show how to integrate a sersic profile, given the coefficients of best-fit. 
Using those results, I convert the calculated total magnitude within a specified radius to a stellar mass using the stellar sysnthesis population models generated with [EZGAL](https://www.baryons.com/ezgal).
Finally, I pull this information together with two other samples of galaxy clusters with stellar mass measurements and look at the stellar mass as a fuction of total cluster mass, M_500.
I fit this relation: log(stellar mass) = alpha*log(m500)+beta with an orthogonal distance regression (ODR), which takes into accoutn uncertainties in both axes.
