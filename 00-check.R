library(gravity)

# bvu ---------------------------------------------------------------------

bvu(y="flow", dist="distw", x=c("rta", "contig", "comcur"), 
    inc_o="gdp_o", inc_d="gdp_d", vce_robust=T, data=gravity_no_zeros)

data(gravity_no_zeros)
# choose exemplarily 10 biggest countries for check data
countries_chosen <- names(sort(table(gravity_no_zeros$iso_o), decreasing = TRUE)[1:10])
grav_small <- gravity_no_zeros[gravity_no_zeros$iso_o %in% countries_chosen,]
bvu(y="flow", dist="distw", x=c("rta"), inc_o="gdp_o", inc_d="gdp_d", vce_robust=TRUE, data=grav_small)
