clear all
save `1', replace emptyok

local N 25
forval n=1/`N' {
	import delimited `1'`n'.csv, clear
	gen year = _n+1971
	gen country = `n'
	reshape long x, i(year) j(sector)
	ren x `1'
	append using `1'
	saveold `1', replace
}
