* 1. Parameters ----
clear
set more off
set linesize 180
capture cd "D:/"
capture cd "/Volumes/USB DRIVE/"

global nurse_fld	UKHLS 2-3, Nurse Assessment/stata/stata11
global main_fld		UKHLS 1-8 & BHPS 1-18, Special Licence/stata/stata13_se
global act_fld		Projects/UKHLS Work-Life Histories/Data
global curr_fld 	Projects/Allostatic Load Scarring Paper/Data

	* Unemployment time range parameters
global age_low = 16
global age_high = 24
global fte_high = 5


* 2. Parental NS-SEC ----
* BHPS
tempfile temp
foreach i of numlist 1 8 9 10 11 12 13 14 15 16 17 18{
	local j: word `i' of `c(alpha)'
	use pidp b`j'_panssec8_dv b`j'_paju ///
		using "${main_fld}/bhps_w`i'/b`j'_indresp_protect.dta", clear
	rename b`j'_panssec8_dv FatherNSSEC8_14
	rename b`j'_paju FatherNotWorking_14
	gen wave = `i'
	capture append using "`temp'"
	save "`temp'", replace
}

recode FatherNSSEC8_14 (min/-1=.)
replace FatherNSSEC8_14 = 7 + FatherNotWorking_14 ///
	if inrange(FatherNotWorking_14, 2, 4)	
by pidp (wave), sort: gen ob = _n if !missing(FatherNSSEC8_14)
by pidp (wave), sort: egen first_ob = min(ob)
by pidp (wave), sort: replace FatherNSSEC8_14 = FatherNSSEC8_14[first_ob]

keep pidp FatherNSSEC8_14
duplicates drop
drop if missing(FatherNSSEC8_14)

gen FatherNSSEC5_14 = 0 if inrange(FatherNSSEC8_14, 1, 3)
replace FatherNSSEC5_14 = FatherNSSEC8_14 - 3 if inrange(FatherNSSEC8_14, 4, 6)
replace FatherNSSEC5_14 = 4 if inrange(FatherNSSEC8_14, 7, 8)
replace FatherNSSEC5_14 = FatherNSSEC8_14 - 4 if inrange(FatherNSSEC8_14, 9, 11)

gen FatherNSSEC3_14 = 0 if inrange(FatherNSSEC8_14, 1, 3)
replace FatherNSSEC3_14 = 1 if inrange(FatherNSSEC8_14, 4, 5)
replace FatherNSSEC3_14 = 2 if inrange(FatherNSSEC8_14, 6, 8)
replace FatherNSSEC3_14 = FatherNSSEC8_14 - 6 if inrange(FatherNSSEC8_14, 9, 11)
drop FatherNSSEC8_14
save "${curr_fld}/family_nssec", replace

* UKHLS
tempfile temp
forval i = 1/8{
	local j: word `i' of `c(alpha)'
	use pidp `j'_paju `j'_panssec5_dv ///
		using "${main_fld}/ukhls_w`i'/`j'_indresp_protect.dta", clear
	rename `j'_panssec5_dv FatherNSSEC5_14
	rename `j'_paju FatherNotWorking_14
	gen wave = `i'
	capture append using "`temp'"
	save "`temp'", replace
}

recode FatherNSSEC5_14 (min/-1=.)
replace FatherNSSEC5_14 = FatherNSSEC5_14 - 1
replace FatherNSSEC5_14 = 3 + FatherNotWorking_14 ///
	if inrange(FatherNotWorking_14, 2, 4)
by pidp (wave), sort: gen ob = _n if !missing(FatherNSSEC5_14)
by pidp (wave), sort: egen first_ob = min(ob)
by pidp (wave), sort: replace FatherNSSEC5_14 = FatherNSSEC5_14[first_ob]

keep pidp FatherNSSEC5_14
duplicates drop	
drop if missing(FatherNSSEC5_14)
label values FatherNSSEC5_14
	
gen FatherNSSEC3_14 = 0 if FatherNSSEC5_14==0
replace FatherNSSEC3_14 = 1 if inrange(FatherNSSEC5_14, 1, 2)
replace FatherNSSEC3_14 = 2 if inrange(FatherNSSEC5_14, 3, 4)
replace FatherNSSEC3_14 = FatherNSSEC5_14 - 2 if inrange(FatherNSSEC5_14, 5, 7)
merge 1:1 pidp using "${curr_fld}/family_nssec", nogen update

compress
save "${curr_fld}/family_nssec", replace


* 2. Collect Main Variables ----
	* a. Load Data
global xwave 	dobm_dv doby_dv ukborn lvag16 agelh lvag14 racel_dv yr2uk4 ///
				maedqf paedqf
global xns		scghq1_dv hiqual_dv jbstat hiqual_dv jbnssec5_dv ///
				sf1 health marstat wave strata psu
global bind 	b_hidp b_scghq1_dv b_hhorig b_sampst b_smever b_smnow ///
				b_ncigs b_smcigs b_sceverdrnk b_scfalcdrnk b_scalcl7d
global cind		c_hidp c_scghq1_dv
global dind		d_scghq1_dv
global orig		$xwave $xns $bind $cind $dind
global bhh		b_fihhmngrs_dv b_tenure_dv
global chh		c_fihhmngrs_dv c_tenure_dv

tempfile Temp
use pidp $xwave pa* using "${main_fld}/ukhls_wx/xwavedat_protect", clear
merge 1:1 pidp using "${nurse_fld}/xindresp_ns", ///
	nogen keepusing($xns)
preserve
	use pidp $bind using "${main_fld}/ukhls_w2/b_indresp_protect", clear
	save "`Temp'", replace
restore
merge 1:1 pidp using "`Temp'",  nogen
preserve
	use pidp $cind using "${main_fld}/ukhls_w3/c_indresp_protect", clear
	save "`Temp'", replace
restore
merge 1:1 pidp using "`Temp'",  nogen
preserve
	use pidp $dind using "${main_fld}/ukhls_w4/d_indresp_protect", clear
	save "`Temp'", replace
restore
merge 1:1 pidp using "`Temp'",  nogen
preserve 
	use pidp bm_lvag* using "${main_fld}/bhps_w13/bm_indresp_protect", clear
	save "`Temp'", replace
restore
merge 1:1 pidp using "`Temp'",  nogen
merge m:1 b_hidp ///
	using "${main_fld}/ukhls_w2/b_hhresp_protect", ///
	nogen keepusing($bhh)
merge m:1 c_hidp ///
	using "${main_fld}/ukhls_w3/c_hhresp_protect", ///
	nogen keepusing($chh)
keep pidp $orig $bhh $chh bm_lvag*
drop if missing(pidp)

cls
qui labelbook
label drop `r(notused)'
numlabel, add

recode $orig *tenure_dv bm_lvag*  (min/-1=.)

	* b. Clean variables
sum psu
gen PSU=psu
replace PSU = `r(max)' + _n if missing(PSU)
gen Strata=strata
gen hhorig=b_hhorig
label values hhorig b_hhorig
gen sampst=b_sampst
label values sampst b_sampst

gen Birth_MY=ym(doby_dv,dobm_dv)
gen Foreign=ukborn==5 if !missing(ukborn)
gen AgetoUK=(ym(yr2uk4,6)-Birth_MY)/12
replace AgetoUK=0 if Foreign==0

gen ParentsHH_14=0 if lvag16==1 | inrange(lvag14,1,4) ///
						| bm_lvag16==1 | inrange(bm_lvag14_bh,1,4) 
replace ParentsHH_14=1 if inrange(lvag14,5,6) | inrange(bm_lvag14_bh,5,6)
replace ParentsHH_14=2 if inrange(lvag14,7,97) | inrange(bm_lvag14_bh,7,8)

merge 1:1 pidp using "${curr_fld}/family_nssec", nogen

gen Ethnicity = 1-(inrange(racel_dv,1,4)) if !missing(racel_dv)

gen FatherEdu=0 if inlist(paedqf,1,2)
replace FatherEdu=paedqf-2 if inrange(paedqf,3,5)
replace FatherEdu=1 if paedqf==97
gen MotherEdu=0 if inlist(maedqf,1,2)
replace MotherEdu=maedqf-2 if inrange(maedqf,3,5)
replace MotherEdu=1 if maedqf==97
gen ParentEdu=max(FatherEdu,MotherEdu)

gen Education=hiqual_dv-1 if inrange(hiqual_dv,1,5)
replace Education=4 if hiqual_dv==9

gen GHQ_Likert = cond(wave == 2, c_scghq1_dv, d_scghq1_dv) if inlist(wave, 2, 3)
recode GHQ_Likert (min/-1 = .)

gen Income=b_fihhmngrs_dv if wave==2 
replace Income=c_fihhmngrs_dv if wave==3
replace Income=0.01 if Income==0
replace Income=log(Income)
qui sum Income, d
replace Income=. if !inrange(Income,r(p1),r(p99))

gen tenure=cond(wave==2,b_tenure,c_tenure)
gen Tenure=0 if tenure==1
replace Tenure=1 if tenure==2
replace Tenure=2 if inrange(tenure,5,8)
replace Tenure=3 if inlist(tenure,3,4)
drop tenure

gen Status=0 if inlist(jbstat,2,9,10)
replace Status=1 if jbstat==1
replace Status=2 if jbstat==3
replace Status=3 if jbstat==4
replace Status=4 if inlist(jbstat,5,6,7,8,97)
gen NSSEC5=jbnssec5_dv-1 if inrange(jbnssec5_dv,1,5)
replace NSSEC5=3+Status if inrange(Status,2,4)
drop Status

gen Alcohol_W2=0 if b_sceverdrnk==2 | b_scfalcdrnk==8 | (b_scfalcdrnk==9 & b_scalcl7d==2)
replace Alcohol_W2=1 if  inrange(b_scfalcdrnk,6,7)
replace Alcohol_W2=2 if  b_scfalcdrnk==5
replace Alcohol_W2=3 if  b_scfalcdrnk==4
replace Alcohol_W2=4 if  b_scfalcdrnk==3
replace Alcohol_W2=5 if  inrange(b_scfalcdrnk, 1, 2)

gen Smoke_W2=0 if b_smever==2 
replace Smoke_W2=1 if inrange(b_smcigs,2,3) | b_ncigs==0
replace Smoke_W2=2 if b_smcigs==1
replace Smoke_W2=3 if inrange(b_ncigs,1,10)
replace Smoke_W2=4 if b_ncigs>10 & !missing(b_ncigs)
	
drop $orig $bhh $chh bm_lvag*
compress

label define Binary 0 "No" 1 "Yes"
label define Ethnicity 0 "White" 1 "Non-White"
label define Foreign 0 "UK-Born" 1 "Foreign-Born"
label define ParentsHH_14 0 "Two Parents" 1 "Single Parent" 2 "Other"
label define ParentEdu 0 "No Quals" 1 "School/Other" 2 "FE" 3 "HE"
label define Education 0 "Degree" 1 "Other HE" 2 "FE" 3 "GCSE" 4 "None/Other"
label define Alcohol_W2 0 "Not in last 12 months" 1 "Every few months" 2 "Every month" ///
	3 "Every week" 4 "3-4 times a week" 5 "5+ times a week"
label define NSSEC5 0 "Higher" 1 "Intermediate" 2 "Small Employers" 3 "Lower Supervisory" ///
					4 "Routine/Manual" 5 "Unemployed" 6 "Retired" 7 "Inactive"
label define ParentNSSEC5 0 "Higher" 1 "Intermediate" 2 "Small Employers" 3 "Lower Supervisory" ///
					4 "Routine/Manual" 5 "Not Working" 6 "Deceased" 7 "Not Present"
label define ParentNSSEC3 0 "Higher" 1 "Intermediate" 2 "Routine" ///
					3 "Not Working" 4 "Deceased" 5 "Not Present"
label define Smoke_W2 0 "Never Smoker" 1 "Occassional" 2 "Former" 3 "1-10" 4 "> 10"
label define Tenure 0 "Own" 1 "Mortgage" 2 "Private Rent" 3 "Social Rent"

label values Ethnicity Ethnicity 
label values Foreign Foreign
label values ParentsHH_14* ParentsHH_14
label values FatherEdu MotherEdu ParentEdu ParentEdu
label values Education Education
label values Alcohol_W2 Alcohol_W2
label values NSSEC5 NSSEC5
label values *NSSEC5_14 ParentNSSEC5
label values *NSSEC3_14 ParentNSSEC3
label values Smoke_W2 Smoke_W2
label values Tenure Tenure

save "${curr_fld}/collected_variables", replace


* 3. Unemployment Duration ----
capture program drop prog_unem
program define prog_unem
	args var lb ub
	replace Start_MY = max(Start_MY, `lb')
	replace End_MY = min(End_MY, `ub')
	drop if Start_MY >= End_MY

	gen Unem_Duration = End_MY - Start_MY if Status == 3
	by pidp (Spell), sort: egen `var' = max(Unem_Duration)
	replace `var' = `var'>=6 & !missing(`var')
	by pidp (Spell), sort: replace `var' = . ///
		if Start_MY[1] > `lb' & `var' == 0
	gen Missing = Status ==.m & End_MY - Start_MY >= 6
	by pidp (Spell), sort: egen Any_Missing = max(Missing)	
	replace `var' = . if Any_Missing == 1 & `var' == 0
	label define Unem 0 "<6 Months Unemployment" 1 "6+ Months Unemployment"
	label values `var' Unem	
	drop Missing Any_Missing Unem_Duration
end

* a. Birth_MY
use "${act_fld}/Merged Dataset" if !missing(Birth_MY), clear
prog_unem Unem_Age "Birth_MY + ${age_low}*12" "Birth_MY + `=${age_high}+1'*12"
keep pidp Unem_Age
duplicates drop
compress
save "${curr_fld}/unemployment_duration", replace


* 4. Allostatic Load ----
* Nurse Assessment
use "${nurse_fld}/xindresp_ns", clear
recode * (min/-1=.)

gen Pulse=ompulval
gen Systolic=omsysval + 10*inlist(hyper2om,2,3) if !missing(hyper2om)
replace Systolic = omsysval + 10*(inlist(hyper1, 2, 3)) ///
    if !missing(hyper1) & missing(Systolic)
gen Diastolic=omdiaval + 5*inlist(hyper2om,2,3) if !missing(hyper2om)
replace Diastolic = omdiaval + 10*(inlist(hyper1, 2, 3)) ///
    if !missing(hyper1) & missing(Diastolic)
gen WtH_Ratio=wstval/htval

gen Nurse_Weight_X=indnsub_xw

gen Gender = (nsex==2) if inrange(nsex,1,2)
label define Gender 0 "Male" 1 "Female"
label values Gender Gender
    
gen Weight=wtval if wtok==1 & wtval>=0
gen Age = confage
gen IntDate_MY=ym(nurdayy,nurdaym) if nurdaym>0 & nurdayy>0
rename wave Wave
label define Wave 2 "Wave 2" 3 "Wave 3"
label values Wave Wave

global nurse     Pulse Systolic Diastolic WtH_Ratio	
keep pidp Gender Age IntDate_MY Wave Nurse_Weight_X Weight $nurse
tempfile Temp
save "`Temp'", replace
    
    * Blood Biomarkers
    // Use Information from Too Small to Detect (-32 LINEAR RANGE, -31 BELOW DETECTION LIMIT)
use "${nurse_fld}/xlabblood_ns", clear
recode * (min/-1=.)
merge 1:1 pidp using "`Temp'", ///
    keepusing(Weight Gender Age)
gen Blood_Weight_X=indbdub_xw
gen Blood_Sample = inlist(_merge, 1, 3)
drop _merge

    *Lowest
tab1 ecre
rename igfi Insulin
gen Creatinine=(140-Age)*Weight*cond(Gender==1,1.04,1.23)/ecre ///
    if inlist(Gender,0,1)
rename dheas DHEAS
 
    * Highest
tab1 hba1c
rename cfib ClaussFib
rename hscrp CReactive
gen HDL_Ratio=chol/hdl
rename trig Trig
rename hba1c HbA1c


    * Merge
global blood    Insulin Creatinine DHEAS ClaussFib ///
                CReactive HDL_Ratio Trig HbA1c
keep pidp $blood Blood_Weight_X Blood_Sample
merge 1:1 pidp using "`Temp'", nogen
order pidp Blood_Weight_X IntDate_MY Gender Age
compress

preserve
	keep pidp Blood_Sample Blood_Weight_X Age  $nurse $blood
	order pidp Blood_Sample Blood_Weight_X Age  $nurse $blood
	save "${curr_fld}/df_sample", replace
restore
drop Blood_Sample

    * Drop those with likely recent infection
count if CReactive>10 & !missing(CReactive)
di in red "There are `r(N)' cases with C-Reactive>10 mg/L" _newline
drop if CReactive>10 & !missing(CReactive)

	* Drop individuals outside age range
keep if inrange(Age, ${age_high}+1, 64)


    * Create Allostatic Index
global biomarkers $nurse $blood

global lowest	Insulin Creatinine DHEAS
global highest: list global(biomarkers) - global(lowest)

egen miss_bio = rowmiss(${biomarkers})

di in red "Cut-offs for high risk category, male and female"
qui foreach var of global biomarkers{
    local ptile=cond(strpos("$lowest", "`var'") > 0, 25, 75)	
    local op=cond(strpos("$lowest", "`var'") > 0, "<=", ">=")
//  local weight=cond(strpos("$blood","`var'")>0,"Blood_Weight_X","Nurse_Weight_X")
	local weight Blood_Weight_X

    gen `var'_Quartile=.
    local risk	`var':
    forval i=0/1{
        _pctile `var' [pweight=`weight'] if Gender==`i' & miss_bio == 0, percentiles(`ptile')
        replace `var'_Quartile=(`var'`op'r(r1)) if !missing(`var') & Gender==`i'
        local xx=round(r(r1),0.001)
        local xx = substr("`xx'",1,strpos("`xx'",".")+3)
        local risk `risk' `xx'
    }
    di in red "`risk'"
}
egen Allostatic_Index = rowtotal(*Quartile)
replace Allostatic_Index = . if miss_bio != 0
    
    * Missing Patterns
misstable pattern $biomarkers if !missing(Blood_Weight_X), asis

	* Log-Transform
sktest $biomarkers
sum $biomarkers
foreach var of global biomarkers {
	replace `var'=ln(`var')
	
// 	local weight=cond(strpos("$blood","`var'")>0, "Blood_Weight_X", "Nurse_Weight_X")
	local weight Blood_Weight_X
	forval i = 0/1{
		mean `var' [pweight = `weight'] if Gender == `i' & miss_bio == 0
		estat sd
		replace `var' = (`var' - r(mean)[1,1])/ r(sd)[1,1] if Gender == `i'
	}
	label values `var' 
}
foreach var of global lowest{
	replace `var' = -`var'
}
egen Allostatic_Z = rowtotal(${biomarkers})
replace Allostatic_Z = . if miss_bio != 0

forval i = 0/1{
	mean Allostatic_Z [pweight = Blood_Weight_X] if Gender == `i'
	estat sd
	replace Allostatic_Z = (Allostatic_Z - r(mean)[1,1])/ r(sd)[1,1] if Gender == `i'
}

drop miss_bio

    * Format Dataset
drop Weight
save "${curr_fld}/allostatic_load", replace


* 5. Combine Datasets ----
use "${curr_fld}/allostatic_load", clear

* Merge with other data
merge 1:1 pidp using "${curr_fld}/unemployment_duration", keep(match master) nogen 
merge 1:1 pidp using "${curr_fld}/collected_variables", keep(match master) nogen

* Choose Sample
drop if Blood_Weight == 0 | missing(Blood_Weight_X)

* Final Formatting
drop AgetoUK hhorig sampst Birth_MY MotherEdu ///
	IntDate_MY Nurse_Weight_X Strata FatherNSSEC3_14
compress
save "${curr_fld}/df_analysis", replace
