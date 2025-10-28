capture program drop aldsim
program define aldsim, rclass
syntax [anything] , cn(numlist integer max=1)  pn(numlist integer max=1) pi(numlist max=1) [ci(numlist max=1) nc(numlist integer max=1) totn(numlist max=1) ///
 slp(numlist max=1) slpsd(numlist max=1) slpintcor(numlist max=1) ///
 int(numlist max=1) intsd(numlist max=1) ///
 cid(numlist max=1)  csd(numlist max=1) cidper(numlist max=1)  csdper(numlist max=1)  cidratio(numlist max=1) csdratio(numlist max=1) csdtotal cidtotal fixedcohort  cidrange(numlist max=1 >=0 <=1) csdrange(numlist max=1 >=0 <=1) ///
 agestart(numlist max=1) agesd(numlist max=1) agelb(numlist max=1) ageub(numlist max=1)  ///
 resid(numlist max=1) ///
 alpha(numlist max=1) alphasd(numlist max=1) nltype(string)   ageinflect(numlist max=1) ageinflectshift(numlist max=1) ///
 ATTRition(numlist max=1 >=0 <1 ) gamma(numlist max=1) attrtype(string) attrstatic ///
 period(numlist max=1) periodshift(numlist max=1) ///
 effsize(numlist max=1) effsizecon gcr(numlist max=1) graph print]

qui {


if "`totn'" == "" & "`nc'"=="" {
	dis as error "Either nc or totn option must be specified"
	exit
}


if "`totn'" != "" {
		local nc=floor(`totn'/`cn')
}

if "`ci'"== "" {
	local ci=`pi'
	noi dis as text "ci set equal to pi (value=`pi')"
}

*********************************
***Specification of Effect Sizes with Defaults
*********************************

	*Slope and Effect Size Specified
	if "`slpsd'" == "" & "`effsize'" != "" & "`slp'" != "" {
		local slpsd=`slp'/`effsize'
	}
	*Slope and Slope SD Specified
	if "`slpsd'" != "" & "`effsize'" == "" & "`slp'" != "" {
		local effsize=`slp'/`slpsd'
	}
	*Slope SD  and Effect Size Specified
	if "`slpsd'" != "" & "`effsize'" != "" & "`slp'" == "" {
		local slp = `effsize'*`slpsd'
	}
	*Effect Size alone Specified
	if "`slpsd'" == "" & "`effsize'" != "" & "`slp'" == "" {
		local slp=1
		local slpsd = `slp'/`effsize'
	}
	*Slope alone Specified
	if "`slpsd'" == "" & "`effsize'" == "" & "`slp'" != "" {
		local effsize=0.5
		local slpsd = `slp'/`effsize'
	}
	*Slope SD alone Specified
	if "`slpsd'" != "" & "`effsize'" == "" & "`slp'" == "" {
		local effsize=0.5
		local slp = `effsize'*`slpsd'
	}
	*None Specified
	if "`slpsd'" == "" & "`effsize'" == "" & "`slp'" == "" {
		local effsize=0.5
		local slp=1
		local slpsd = `slp'/`effsize'
	}	
	*ALL Specified
	if "`slpsd'" != "" & "`effsize'" != "" & "`slp'" != "" {
		local temp=`slp'/`effsize'
		if  `temp' != `slpsd' & `slpsd' !=0  & `temp' !=.{
			noi dis as error "Slope SD changed to `temp' to match effect size"
			local slpsd = `slp'/`effsize'	
		}
		
	}	


*********************************
***Specification of Intercept Means/SDs
*********************************	

if "`int'" == "" & "`intsd'" == "" {
	local int=0
	local intsd=`slpsd'

}
if "`int'" != "" & "`intsd'" == "" {
	local intsd=`slpsd'
}
if "`int'" == "" & "`intsd'" != "" {
	local int=0
}

********************************************
***Specification of Cohort Slope Differences
********************************************

*Neither Slope Differences or Slope Ratio Specified
if "`csdratio'" == "" &  "`csd'" == "" {
	*local csdratio=0
	local csd=0
}
*Slope Ratio Only Supplied--And not desired to be the total difference from first to last cohort
if "`csdratio'" != "" &  "`csd'" == "" & "`csdtotal'"=="" {
	local csd=`slpsd'*`csdratio'
}
*Slope Ratio Only Supplied-restrict total difference
if "`csdratio'" != "" &  "`csd'" == "" & "`csdtotal'"!="" {
	local csd=(`slpsd'*`csdratio')*(`cn'-1)
}
*Slope Differences Only Supplied-And not desired to be the total difference from first to last cohort
if "`csdratio'" == "" &  "`csd'" != "" & "`csdtotal'"=="" {
	local csdratio=`csd'/`slpsd'
}
*Percentage Growth change
if "`csdper'" != "" {
	local csd=`slp'*`csdper'
}
*Slope Differences Only Supplied-restrict total difference
if "`csdratio'" == "" &  "`csd'" != "" & "`csdtotal'"!="" {
	if `cn'==1 {
		local csd=`csd'/(`cn')
	}
	if `cn'!=1 {
		local csd=`csd'/(`cn'-1)
	}
	local csdratio=`csd'/`slpsd'
}


*Both Slope Differences and Ratio Supplied
if "`csdratio'" != "" &  "`csd'" != "" {
	local temp =`csd'/`csdratio'
	if `temp' != `slpsd' & `temp' !=. {
		noi dis as error "Slope SD changed to `temp' to match CSD info"
		local slpsd=`csd'/`csdratio'
		if "`csdtotal'" != "" {
			noi dis as error "NOTE: CSDTOTAL option cannot be used when both CSD and CSD-Ratio are supplied"
		}
	}
	
	local temp=`effsize'*`slpsd'
	if `temp' != `slp' & `temp' !=.{
		noi dis as error "Slope changed to `temp' to match Effect Size and Slope SD Info"
		local slp=`effsize'*`slpsd'
		if "`csdtotal'" != "" {
			noi dis as error "NOTE: CSDTOTAL option cannot be used when both CSD and CSD-Ratio are supplied"
		}		
	}
	
}


********************************************
***Specification of Cohort Intercept Differences
********************************************

*Neither Intercept Differences or Slope Ratio Specified
if "`cidratio'" == "" &  "`cid'" == "" {
	*local cidratio=0
	local cid=0
}
*Intercept Ratio Only Supplied
if "`cidratio'" != "" &  "`cid'" == "" & "`cidtotal'"=="" {
	local cid=`intsd'*`intratio'
}
*Intercept Ratio Only Supplied-restrict total difference
if "`cidratio'" != "" &  "`cid'" == "" & "`cidtotal'"!="" {
	local cid=(`intsd'*`intratio')*(`cn'-1)
}
*Intercept Differences Only Supplied
if "`cidratio'" == "" &  "`cid'" != "" & "`cidtotal'"=="" {
	local cidratio=`cid'/`intsd'
}
*Percentage Growth change
if "`cidper'" != "" {
		local cid=`int'*`cidper'
}
*Intercept Differences Only Supplied-restrict total difference
if "`cidratio'" == "" &  "`cid'" != "" & "`cidtotal'"!="" {
	if `cn'==1 {
		local cid=`cid'/(`cn')
	}
	if `cn'!=1 {
		local cid=`cid'/(`cn'-1)
	}
	local cidratio=`cid'/`intsd'
}
*Both Intercept Differences and Ratio Supplied
if "`cidratio'" != "" &  "`cid'" != "" {
	local temp =`cid'/`cidratio'
	if `temp' != `intsd' & `temp' !=. {
		noi dis as error "Intercept SD changed to `temp' to match CID info"
		local intsd =`cid'/`cidratio'
		if "`cidtotal'" != "" {
			noi dis as error "NOTE: CIDTOTAL option cannot be used when both CID and CID-Ratio are supplied"
		}	
	}
	
}



******************************************************************
***Computation of Residuals Based on GCR. 
*******GCR used to determine residual strength
******************************************************************
	/*Formula for: gcr=(`interceptsd'^2 + (`agestart'^2*`slopesd'^2) + (2*`agestart'*`slpintrho'*`slopesd'*`interceptsd')) 
				/ (`interceptsd'^2 + (`agestart'^2*`slopesd'^2) + (2*`agestart'*`slpintrho'*`slopesd'*`interceptsd') + `resid'^2)
	*/
	
	if "`slpintcor'" == "" {
		local slpintcor=0
	}
	*If using Random Uniform for age distribution
	if "`agestart'" == "" & "`ageub'" != "" & "`agelb'" != "" & "`agesd'" == "" {
		local agestart=(`ageub'-`agelb')/2
	}
	if "`agestart'" == "" & "`ageub'" == "" & "`agelb'" == "" & "`agesd'" == "" {
		noi dis as error "Must specify starting age or age range information. See 'help aldsim' for more information."
		exit
	}
	
	*For linear models
	if "`alpha'"=="" {
	
		if "`resid'" == "" & "`gcr'" == "" {	
			local gcr=0.8
			local resid=sqrt(((`intsd'^2 + (`agestart'^2*`slpsd'^2) + (2*`agestart'*`slpintcor'*`slpsd'*`intsd'))/`gcr') - (`intsd'^2 + (`agestart'^2*`slpsd'^2) + (2*`agestart'*`slpintcor'*`slpsd'*`intsd'))) 
		}
		if "`resid'" != "" & "`gcr'" == "" {	
			local gcr=(`intsd'^2 + (`agestart'^2*`slpsd'^2) + (2*`agestart'*`slpintcor'*`slpsd'*`intsd')) / (`intsd'^2 + (`agestart'^2*`slpsd'^2) + (2*`agestart'*`slpintcor'*`slpsd'*`intsd') + `resid'^2)
		}

		if "`resid'" == "" & "`gcr'" != "" {	
			local resid=sqrt(((`intsd'^2 + (`agestart'^2*`slpsd'^2) + (2*`agestart'*`slpintcor'*`slpsd'*`intsd'))/`gcr') - (`intsd'^2 + (`agestart'^2*`slpsd'^2) + (2*`agestart'*`slpintcor'*`slpsd'*`intsd'))) 
		}
		if "`resid'" != "" & "`gcr'" != "" {	
			local temp=sqrt(((`intsd'^2 + (`agestart'^2*`slpsd'^2) + (2*`agestart'*`slpintcor'*`slpsd'*`intsd'))/`gcr') - (`intsd'^2 + (`agestart'^2*`slpsd'^2) + (2*`agestart'*`slpintcor'*`slpsd'*`intsd'))) 
			if `resid' != `temp' & `temp' != . {
				noi dis as error "Residual changed to `temp' to match GCR info"
				local resid=sqrt(((`intsd'^2 + (`agestart'^2*`slpsd'^2) + (2*`agestart'*`slpintcor'*`slpsd'*`intsd'))/`gcr') - (`intsd'^2 + (`agestart'^2*`slpsd'^2) + (2*`agestart'*`slpintcor'*`slpsd'*`intsd'))) 

			}
		}	
	}

	*For nonlinear models
	if "`alpha'"!="" {
	
		*Inflection Point for Non-Linear Models.
		if "`ageinflect'" == "" & "`ageinflectshift'" == "" {
			local agenlstart= `agestart' - (`agestart' + (`pi'*(`pn'-1) + `ci'*(`cn'-1))/2)
		}
		if "`ageinflect'" == "" & "`ageinflectshift'" != "" {
			local agenlstart= `agestart' - (`agestart' + (`pi'*(`pn'-1) + `ci'*(`cn'-1))/2 + (`j'*`ageinflectshift'))
		}
		if "`ageinflect'" != "" & "`ageinflectshift'" == "" {
			local agenlstart=`agestart'-`ageinflect'
		}
		if "`ageinflect'" != "" & "`ageinflectshift'" != "" {
			local agenlstart=`agestart'-(`ageinflect'+ (`j'*`ageinflectshift'))
		}	
	
	
		if "`resid'" == "" & "`gcr'" == "" {	
			local gcr=0.8
			local resid=sqrt(((`intsd'^2 + (`agenlstart'^2*`slpsd'^2) + (2*`agenlstart'*`slpintcor'*`slpsd'*`intsd'))/`gcr') - (`intsd'^2 + (`agenlstart'^2*`slpsd'^2) + (2*`agenlstart'*`slpintcor'*`slpsd'*`intsd'))) 
		}
		if "`resid'" != "" & "`gcr'" == "" {	
			local gcr=(`intsd'^2 + (`agenlstart'^2*`slpsd'^2) + (2*`agenlstart'*`slpintcor'*`slpsd'*`intsd')) / (`intsd'^2 + (`agenlstart'^2*`slpsd'^2) + (2*`agenlstart'*`slpintcor'*`slpsd'*`intsd') + `resid'^2)
		}

		if "`resid'" == "" & "`gcr'" != "" {	
			local resid=sqrt(((`intsd'^2 + (`agenlstart'^2*`slpsd'^2) + (2*`agenlstart'*`slpintcor'*`slpsd'*`intsd'))/`gcr') - (`intsd'^2 + (`agenlstart'^2*`slpsd'^2) + (2*`agenlstart'*`slpintcor'*`slpsd'*`intsd'))) 
		}
		if "`resid'" != "" & "`gcr'" != "" {	
			local temp=sqrt(((`intsd'^2 + (`agenlstart'^2*`slpsd'^2) + (2*`agenlstart'*`slpintcor'*`slpsd'*`intsd'))/`gcr') - (`intsd'^2 + (`agenlstart'^2*`slpsd'^2) + (2*`agenlstart'*`slpintcor'*`slpsd'*`intsd'))) 
			if `resid' != `temp' & `temp' != . {
				noi dis as error "Residual changed to `temp' to match GCR info"
				local resid=sqrt(((`intsd'^2 + (`agenlstart'^2*`slpsd'^2) + (2*`agenlstart'*`slpintcor'*`slpsd'*`intsd'))/`gcr') - (`intsd'^2 + (`agenlstart'^2*`slpsd'^2) + (2*`agenlstart'*`slpintcor'*`slpsd'*`intsd'))) 

			}
		}	
	}	


******************************************************************************************
*Attrition

*If attrition is present, check to ensure other dependencies are met. Create larger output matrix
if "`attrition'" != "" {	
	capture assert "`gamma'" != ""
	if _rc!= 0 {
		noi dis as error "{it:gamma} option must be specified when incorporating attrition. Type 'help aldsim' for more information"
		exit
	} 
	capture assert "`attrtype'" !=""
	if _rc!= 0 {
		noi dis as error "{it:attrtype} option must be specified when incorporating attrition. Type 'help aldsim' for more information"
		exit
	}
	
	
	if "`attrtype'" == "scd" {
		*Number of Observations in equivalent SCD - Assumed Balanced Design
		local M=ceil((`pi'*(`pn'-1) + `ci'*(`cn'-1))/`pi' + 1)	
		
		clear
		set obs `M'
		gen measurement=_n
		gen t_p=(measurement-1)/(`M'-1)
		gen prob=(1-`attrition')^(t_p^`gamma') - (1-`attrition')^((t_p[_n+1])^`gamma')
		replace prob=1-`attrition' if prob==.
		gen hicut=sum(prob)
		drop if measurement > `pn'
		replace hicut=1 if measurement==`pn'
		gen locut=hicut[_n-1]
		replace locut=0 in 1
	}
	if "`attrtype'" == "ald" {
	
		clear
		set obs `pn'
		gen measurement=_n
		gen t_p=(measurement-1)/(`pn'-1)
		gen prob=(1-`attrition')^(t_p^`gamma') - (1-`attrition')^((t_p[_n+1])^`gamma')
		replace prob=1-`attrition' if prob==.
		gen hicut=sum(prob)	
		gen locut=hicut[_n-1]
		replace locut=0 in 1	
	}
	*Assign Probability Cut-Points for attrition 
	forvalues i=1(1)`pn' {
		local hicut`i'=hicut in `i'
		local locut`i'=locut in `i'
	}		
}

*********************************************************************************
*Simulation within Each Cohort
forvalues x=1(1)`cn' {
	local j=`x'-1
	
	
	********************************************************************
	*Determine Values for Slope and Intercept for each Cohort
	********************************************************************
	
	if "`fixedcohort'" != "" {
		local cid`x'=(`j'*`cid')
		local ival=`int'+(`j'*`cid')
		
		
		local csd`x'=(`j'*`csd')
		local sval=`slp' + (`j'*`csd')		
	}
	else {
		if "`cidrange'" != "" {
			local cid_lb=`cid'*(1-`cidrange')
			local cid_ub=`cid'*(1+`cidrange')
			local cid`x'=(`cid_lb'+runiform()*(`cid_ub'-`cid_lb'))*`j'
		}
		if "`cidrange'" == "" {
			local cid`x'=abs(rnormal(0,`cid'))*`j'
		}
		if "`csdrange'" != "" {
			local csd_lb=`csd'*(1-`csdrange')
			local csd_ub=`csd'*(1+`csdrange')
			local csd`x'=(`csd_lb'+runiform()*(`csd_ub'-`csd_lb'))*`j'
		}			
		if "`csdrange'" == "" {
			local csd`x'=abs(rnormal(0,`csd'))*`j'
		}

		local ival=`int'+`cid`x''
		local sval=`slp' +`csd`x''
	}

	*If Cohort Slope EffectSize is to be Constant in the presence of CSD
	if "`effsizecon'" != "" {
		local slpsd = `sval'/`effsize'
	} 


	matrix C=(1, `slpintcor' \ `slpintcor', 1)
	drawnorm intercept slope, n(`nc') means(`ival' `sval') sds(`intsd' `slpsd')  corr(C)  clear 
	
	
	
	*********************************************
	*Determine Starting Age (Age distribution at first period)
	*********************************************
	
	*If Truncated Normal Distribution Specified
	if "`agelb'" != "" & "`ageub'" != "" & "`agesd'" != "" {
		if `x'==1 {
			noi dis as text "Using Truncated Normal Distribution for Age"
		}
		local mval=`agestart' + (`j'*`ci')
		local mval_lb=`agelb'+ (`j'*`ci')
		local mval_ub=`ageub'+ (`j'*`ci')
	
		rtnormal, mean(`mval') sd(`agesd') lower(`mval_lb') upper(`mval_ub') gen(age1)
	}
	*If Uniform Distribution Specified
	if "`agelb'" != "" & "`ageub'" != "" & "`agesd'" == "" {
		if `x'==1 {
			noi dis as text "Using Uniform Distribution for Age"
		}
		local mval_lb=`agelb'+ (`j'*`ci')
		local mval_ub=`ageub'+ (`j'*`ci')
		
		gen age1=runiform(`mval_lb',`mval_ub') 
	}
	*If Normal Distribution Specified	
	if "`agelb'" == "" & "`ageub'" == "" & "`agesd'" != "" {
		if `x'==1 {
			noi dis as text  "Using Normal Distribution for Age"
		}
		local mval=`agestart' + (`j'*`ci')
		gen age1=rnormal(`mval', `agesd')
	}
	
	*Errors in age specification
	if "`agelb'" == "" & "`ageub'" == "" & "`agesd'" == "" {
		if `x'==1 {
			noi dis as error "Must specify age SD and/or upper/lower bounds"
			continue, break
		}
	}
	if "`agelb'" != "" & "`ageub'" == "" {
		if `x'==1 {
			noi dis as error "Must specify age upper bound when using lower bound "
			continue, break
		}
	}	
	if "`agelb'" == "" & "`ageub'" != "" {
		if `x'==1 {
			noi dis as error "Must specify age lower bound when using upper bound "
			continue, break
		}	
	}		
	
	
	
	********************************************************************
	*Generate Age Distributions for other periods 
	********************************************************************	
	
	forvalues i=2(1)`pn' {
		local z=`i'-1
		gen age`i'=age1+(`z'*`pi')
	}
	
	
	*Remove values when attrition is present
	if "`attrition'" != "" {
		
		if "`attrstatic'" != "" {
			gen rand=_n/_N
		}
		else {
			gen rand=runiform()
		}
		forvalues i=1(1)`pn' {
			forvalues y=1(1)`pn' {
				local z=-1*`y' + `pn' +1
				if `i' < `z' {
					replace age`z' = .  if rand >= `locut`i'' & rand < `hicut`i''
				}
			}
		}
		drop rand	
	}
	
	
	********************************************************************
	*Generate Outcomes
	********************************************************************	
	gen id=_n+(`j'*`nc')
	if "`alphasd'" != "" {
			gen alpha=rnormal(`alpha', `alphasd')
	}
	
	reshape long age, i(id) j(period)
		
	*Linear	
	if "`alpha'"=="" {
		if `x'==1 {
			noi dis as result "Linear Model Assumed"
		}
		if "`period'" == "" {
			gen y=intercept+(slope*(age))+rnormal(0, `resid')
		}
		if "`period'" != "" {
			if "`periodshift'" == "" {
				gen y=intercept+(slope*(age)) + `period'*(period-1) + rnormal(0, `resid')
			}
			if "`periodshift'" != "" {
				gen y=intercept+(slope*(age)) + (`period'+(`j'*`periodshift'))*(period-1) +  + rnormal(0, `resid')
			}
		
		}
	}
	*Nonliear
	if "`alpha'"!="" {
		if `x'==1 {
			noi dis as result "Non-Linear Model Assumed"
		}
		*Inflection Point for Non-Linear Models.
		if "`ageinflect'" == "" & "`ageinflectshift'" == "" {
			gen age_c= age - (`agestart' + (`pi'*(`pn'-1) + `ci'*(`cn'-1))/2)
		}
		if "`ageinflect'" == "" & "`ageinflectshift'" != "" {
			gen age_c= age - (`agestart' + (`pi'*(`pn'-1) + `ci'*(`cn'-1))/2 + (`j'*`ageinflectshift'))
		}
		if "`ageinflect'" != "" & "`ageinflectshift'" == "" {
			gen age_c=age-`ageinflect'
		}
		if "`ageinflect'" != "" & "`ageinflectshift'" != "" {
			gen age_c=age-(`ageinflect'+ (`j'*`ageinflectshift'))
		}	
		
		*Additive Models (alphasd=0)
		if "`alphasd'" == "" {	
			if "`nltype'" == "exp" {
				gen y=intercept+slope*(1-exp(-`alpha'*age_c))+rnormal(0, `resid')     /*Nonlinear Exponential Additive Model*/
			}
			if "`nltype'" == "log" | "`nltype'" == ""  {
				gen y=intercept+(slope/(1+exp(-`alpha'*age_c)))+rnormal(0, `resid')   /*Logistic Normal Model*/
			}
			if "`nltype'" == "gom" {
				gen y=intercept+(slope*exp(-exp(-`alpha'*age_c))) +rnormal(0, `resid') /*Gomepertz*/
			}
		}
		*Multiplicative Models (alphasd!=0)
		if "`alphasd'" != "" {
			
			if "`nltype'" == "exp" {
				gen y=intercept+slope*(1-exp(-alpha*age_c))+rnormal(0, `resid')     /*Nonlinear Exponential Additive Model*/
			}
			if "`nltype'" == "log" | "`nltype'" == ""  {
				gen y=intercept+(slope/(1+exp(-alpha*age_c)))+rnormal(0, `resid')   /*Logistic Normal Model*/
			}
			if "`nltype'" == "gom" {
				gen y=intercept+(slope*exp(-exp(-alpha*age_c))) +rnormal(0, `resid') /*Gomepertz*/
			}
			drop alpha
		}	
		*drop age_c
	}	
	
	
	drop slope intercept	
	gen cohort=`x'
	tempfile temp`x'
	save `temp`x'', replace 
		
}/*foreach cohort*/
	
	
*Combine together	
capture use `temp1', clear
if _rc != 0 {
	clear
	exit
}
forvalues i=2(1)`cn' {
	append using `temp`i''
}	
	
	
*If graph desired
if "`graph'" != "" {
	local lines
	forvalues i=1(1)`cn' {

		local j=`i'-1
		local lines `lines' (line y age if cohort==`i', connect(ascending)) 
	}
	twoway `lines' , legend(off)

}
*

*If output display is desired
		if "`print'" != "" {

			
				noisily di as text _col(5) "_________________________________________________________________________________________________________________________________________________________________________"  ///
				_newline  as text _col(5) "---Design Elements---" _col(30) "---Slope---"	 _col(45) "-Intercept-"	_col(60) "-Effect Size-"	 _col(79) "--Between Cohort Characteristics--"			///
				_newline as text _col(5) "Nc"  _col(12) "Cn, Ci, Pn, Pi"	 ///
				_col(31) "b" _col(37) "SD"   ///
				_col(46) "b" _col(52) "SD" ///
				_col(62) "Eff" _col(68) "GCR" ///
				_col(79) "CID" _col(85) "CID-Ratio" _col(98) "CSD"  _col(104) "CSD-Ratio" ///
				_newline as result _col(5) %-6.0fc `nc'  as result _col(12) `cn' as result _col(16) `ci' as result _col(20) `pn' as result _col(24) `pi' ///
				 as result _col(31) %-5.1fc `slp' as result _col(37) %-5.1fc `slpsd'   ///
				as result _col(46) %-5.1fc `int' as result _col(52) %-5.1fc `intsd'   ///
				as result _col(62) %-5.1fc `effsize' as result _col(68) %-5.1fc `gcr' ///
				as result _col(79) %-5.1fc `cid' as result _col(88) %-5.1fc `cidratio' as result _col(98) %-5.1fc `csd' as result _col(107) %-5.1fc `csdratio'  
		
		}
		
return scalar  csdratio=`csdratio'



}/*qui*/
end
	
	
