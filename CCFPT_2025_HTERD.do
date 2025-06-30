/**************************************************************************
##	Project: 
##    ``Heterogeneous Treatment Effects in Regression Discontinuity Designs''
##	  by Calonico, Cattaneo, Farrell, Palomba, and Titiunik
##
##	Purpose: 
##    Replicate and extend the heterogeneity analysis of AER 2022, 112(2): 442–493
##    by Akhtari, Moreira, and Trucco, specifically heterogeneity by income
##    as reported in the original Table A.21.
##    https://doi.org/10.3886/E150323V1
##
**************************************************************************/

clear all
net install rdhte, from(https://raw.githubusercontent.com/rdpackages/rdhte/main/stata) replace
net install rdrobust, from(https://raw.githubusercontent.com/rdpackages/rdrobust/master/stata) replace


**************************************************************************
** Read in data and set up variables

import delimited "AMT_2022_AER.csv"

* Set up main variables
gen Y = expthisschl_lessthan2_dpb
gen X = -1*px
gen cluster_var = cod_municipio
gen W = medianincome2000
qui sum W, d
replace W = `r(p99)' if W>=`r(p99)' & W<.
replace W = W/100
gen W2 = W^2

**************************************************************************
** Discretize income by different quantiles

*also save quantile midpoints for plotting

qui sum W
local W_min = `r(min)'
local W_max = `r(max)'

* Median
pctile tmp = W if tag==1, nq(2)
egen W_med = cut(W), at(`W_min' `r(r1)' `W_max') label
replace W_med=1 if W==`W_max'
gen tmp1 = `W_min' in 1
replace tmp1 = tmp[_n-1] in 2
replace tmp1 = `W_max' in 3
gen W_med_midpoints = tmp1[_n] + (tmp1[_n+1] - tmp1[_n])/2
drop tmp tmp1

* Quartiles
pctile tmp = W if tag==1, nq(4)
egen W_qrt = cut(W), at(`W_min' `r(r1)' `r(r2)' `r(r3)' `W_max') label
replace W_qrt=3 if W==`W_max'
gen tmp1 = `W_min' in 1
replace tmp1 = tmp[_n-1] in 2/4
replace tmp1 = `W_max' in 5
gen W_qrt_midpoints = tmp1[_n] + (tmp1[_n+1] - tmp1[_n])/2
drop tmp tmp1

* Deciles
pctile tmp = W if tag==1, nq(10)
egen W_dec = cut(W), at(`W_min' `r(r1)' `r(r2)' `r(r3)' `r(r4)' `r(r5)' `r(r6)' `r(r7)' `r(r8)' `r(r9)' `W_max') label
replace W_dec=9 if W==`W_max'
gen tmp1 = `W_min' in 1
replace tmp1 = tmp[_n-1] in 2/10
replace tmp1 = `W_max' in 11
gen W_dec_midpoints = tmp1[_n] + (tmp1[_n+1] - tmp1[_n])/2
drop tmp tmp1



**************************************************************************
** Estimation and inference

** WARNING: The inference results below will not exactly match what is 
* reported in the the paper because HC3 is not available with clustering 
* in Stata. Here we use HC2 (HC1 is also available). In R the default is HC3.


* Standard RD treatment effect (no heterogeneity)
rdhte Y X, vce(hc2 cluster_var)
local rd_ate_estimate = e(tau_hat)[1,1]


* HTE - Binary
rdhte Y X, covs_hte(i.W_med) vce(hc2 cluster_var)

*save for plotting
gen med_est = .
gen med_ci_l = .
gen med_ci_r = .
forvalues i = 1/2 {
    replace med_est = e(tau_hat)[`i', 1] in `i'
    replace med_ci_l = e(tau_ci_lb)[`i', 1] in `i'
    replace med_ci_r = e(tau_ci_ub)[`i', 1] in `i'
}


* HTE - Quartiles
rdhte Y X, covs_hte(i.W_qrt) vce(hc2 cluster_var)

*save for plotting
	gen qrt_est = .
	gen qrt_ci_l = .
	gen qrt_ci_r = .
	forvalues i = 1/4 {
		replace qrt_est = e(tau_hat)[`i', 1] in `i'
		replace qrt_ci_l = e(tau_ci_lb)[`i', 1] in `i'
		replace qrt_ci_r = e(tau_ci_ub)[`i', 1] in `i'
	}


* HTE - Deciles
rdhte Y X, covs_hte(i.W_dec) vce(hc2 cluster_var)

*save for plotting
	gen dec_est = .
	gen dec_ci_l = .
	gen dec_ci_r = .
	forvalues i = 1/10 {
		replace dec_est = e(tau_hat)[`i', 1] in `i'
		replace dec_ci_l = e(tau_ci_lb)[`i', 1] in `i'
		replace dec_ci_r = e(tau_ci_ub)[`i', 1] in `i'
	}

	
* HTE - Continuous, Linear in Income
rdhte Y X, covs_hte(W) vce(hc2 cluster_var)
*save for plotting
	scalar b0_lin = _b[T]
	scalar b1_lin = _b[T#c.W]

*same result using c. syntax
rdhte Y X, covs_hte(c.W) vce(hc2 cluster_var)


* HTE - Continuous, Quadratic in Income
rdhte Y X, covs_hte(W W2) vce(hc2 cluster_var)
*save for plotting
	scalar b0_quad = _b[T]
	scalar b1_quad = _b[T#c.W]
	scalar b2_quad = _b[T#c.W2]

*same result using c. syntax
rdhte Y X, covs_hte(c.W##c.W) vce(hc2 cluster_var)

*note that single # gives only the quadratic term
rdhte Y X, covs_hte(c.W#c.W) vce(hc2 cluster_var)
rdhte Y X, covs_hte(W2) vce(hc2 cluster_var)










**************************************************************************
** Figure 1 - with linear & quadratic fits

* Create grid for fitted lines
gen W_grid = `W_min' + (`W_max' - `W_min') * (_n - 1) / 499 if _n <= 500

* Generate linear and quadratic predictions
gen linear_fit = b0_lin + b1_lin * W_grid if _n <= 500
gen quad_fit = b0_quad + b1_quad * W_grid + b2_quad * W_grid^2 if _n <= 500

* for axis labels
local b0_lin_round = round(b0_lin,0.001)
local rd_ate_round = round(`rd_ate_estimate',0.001)

* Median
twoway ///
	(line linear_fit W_grid if _n <= 500, lcolor(black)) ///
	(line quad_fit W_grid if _n <= 500, lpattern(dash) lwidth(medthick) lcolor(black)) ///
	(scatter med_est W_med_midpoints, msymbol(O) msize(small) color(black)) ///
	(rcap med_ci_l med_ci_r W_med_midpoints, lcolor(black)) ///
	, ///
	legend(off) ///
	yline(`rd_ate_estimate', lpattern(dot) lwidth(medthick) lcolor(black)) ///
	ylabel(-0.75 0 0.75 `rd_ate_round' `b0_lin_round', angle(horizontal) labsize(small)) ///
	xlabel(0(0.5)2.5, labsize(small)) ///
	xtitle("Income") ///
	ytitle("Estimate")


* Quartiles
twoway ///
	(line linear_fit W_grid if _n <= 500, lcolor(black)) ///
	(line quad_fit W_grid if _n <= 500, lpattern(dash) lwidth(medthick) lcolor(black)) ///
	(scatter qrt_est W_qrt_midpoints, msymbol(O) msize(small) color(black)) ///
	(rcap qrt_ci_l qrt_ci_r W_qrt_midpoints, lcolor(black)) ///
	, ///
	legend(off) ///
	yline(`rd_ate_estimate', lpattern(dot) lwidth(medthick) lcolor(black)) ///
	ylabel(-0.75 0 0.75 `rd_ate_round' `b0_lin_round', angle(horizontal) labsize(small)) ///
	xlabel(0(0.5)2.5, labsize(small)) ///
	xtitle("Income") ///
	ytitle("Estimate")


* Deciles
twoway ///
	(line linear_fit W_grid if _n <= 500, lcolor(black)) ///
	(line quad_fit W_grid if _n <= 500, lpattern(dash) lwidth(medthick) lcolor(black)) ///
	(scatter dec_est W_dec_midpoints, msymbol(O) msize(small) color(black)) ///
	(rcap dec_ci_l dec_ci_r W_dec_midpoints, lcolor(black)) ///
	, ///
	legend(off) ///
	yline(`rd_ate_estimate', lpattern(dot) lwidth(medthick) lcolor(black)) ///
	ylabel(-0.75 0 0.75 `rd_ate_round' `b0_lin_round', angle(horizontal) labsize(small)) ///
	xlabel(0(0.5)2.5, labsize(small)) ///
	xtitle("Income") ///
	ytitle("Estimate")
  
  
  

