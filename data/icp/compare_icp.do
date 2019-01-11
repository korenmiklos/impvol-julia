clear all
import delimited using "icp2005-GD.ZS.csv", varnames(1) case(preserve)
ren CountryCode iso
ren ClassificationCode SNA
ren Value expenditure_share

* private service consumption, government collective consumption, construction
gen byte service = inlist(SNA,"1106","1107","1108","1109","1110","1111","1112","14","1502")
* exclue correction terms from denominator
gen byte in_total = inlist(SNA,"11","14","15")

egen service_expenditure = sum(service*expenditure_share), by(iso)
egen total_expenditure = sum(in_total*expenditure_share), by(iso)

gen service_expenditure_share = service_expenditure/total_expenditure*100

keep iso CountryName service_expenditure_share
drop if missing(iso,CountryName)

* only one row per country
egen tag = tag(iso)
keep if tag
drop tag
* drop aggregage regions
drop if inlist(iso,"AFR","ASP","CIS","OEE","SAM","WAS","WLD")

tempfile icp
save `icp'

clear
import delimited using "../sectoral_value_added.csv", varnames(1) case(preserve)
keep if year==2005

egen total_va = rsum(va_*)
gen service_va_share = va_services/total_va*100

* match on names
ren country CountryName
replace CountryName="Belgium" if CountryName=="Belgium and Luxembourg"
replace CountryName="Korea, Rep." if CountryName=="South Korea"
merge 1:1 CountryName using `icp'

egen ROW_share = mean(cond(_m==2,service_expenditure_share,.))
replace service_expenditure_share = ROW_share if CountryName=="ROW"
drop if _m==2

label var service_va_share "Services share in value added"
label var service_expenditure_share "Services share in final expenditure"

scatter service_va_share service_expenditure_share, scheme(538w) mlabel(CountryName)
graph export service_shares.pdf, replace
