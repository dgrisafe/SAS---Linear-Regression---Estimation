*==================================================================*
Program: Linear Regression Exploratory Data Analysis
Purpose: Exploratory data analysis
Date: 06/04/18
Note1: if comments are in ALL CAPS it is a place to manually enter code
*==================================================================*;

/* Fresh start when running SAS code */

	*clear log and output;
	dm 'log;clear;output;clear;';

	*clear the Results window within your program;
	dm 'odsresults; clear';

	*avoid getting too many labels error;
	ods graphics on / LABELMAX=2000;


*==================================================================*
Initialize dataset and look at all variables
*==================================================================*;

libname datset 'C:\Users\Dom\Documents\Screening Linear\datasets';
libname home 'C:\Users\Dom\Documents\Screening Linear';

*initialize dataset;
data stuff;
	*PATH TO APPROPRIATE DATASET;
	set datset.chol_313;
run;

*look at for all variables dataset;
proc contents;
run;


*<<<<<<<<<<<<<<<==================================================================*
 <<<<<<<<<<<<<<<	DEFINE GLOBAL MACROS:
 <<<<<<<<<<<<<<<		id variables
 <<<<<<<<<<<<<<<		outcome
 <<<<<<<<<<<<<<<		continuous predictors
 <<<<<<<<<<<<<<<		categorical predictors
 <<<<<<<<<<<<<<<		all variables
*<<<<<<<<<<<<<<<==================================================================*;

%let idvar		= ID;
%let outco		= SBP;
%let expos		= CHOL;
%let contcov	= AGE BMI DBP HDL HT LDL SKIN TG VLDL WT;
%let catcov		= FEMALE PROBAND;
*all variables, except id variables;
%let allvar 	= &outco &expos &contcov &catcov;


*==================================================================*
Preliminary data formatting:

	Define formats

	Format and label variables

	Remove unrealistic observations:
		Age range 0 to 122 - Guiness world record
		Weight max 635 kg = 1400 lbs - Guiness world record
		Highest recorded blood pressure in young adults preforming heavy, dynamic weight lifting 370/360 mmHg/mmHg - Narloch (1995) et. al.;
*==================================================================*;
title;

*DEFINE CODED FORMATS;
proc format;
	value ynf 0='No' 1='Yes';
	value FMLF 0='MALE' 1='FEMALE';
run;

data stuff;
	set stuff;

	*ASSIGN CODED FORMATS;
	format proband ynf. female FMLF.;

	*ASSIGN VARIABLE NAMES;
	label CHOL = 'Total cholesterol (dg/mL)';
	label AGE = 'Age (years)';
	label BMI = 'Body mass index (kg/m^2)';
	label DBP = 'Diastolic blood pressure (mmHg)';
	label HDL = 'High density lipoprotein (mg/dL)';
	label HT = 'Height (in)';
	label LDL = 'Low density lipoprotein (mg/dL)';
	label SKIN = 'Blood glucose levels (mg/dL)';
	label TG = 'Triglycerides (mg/dL)';
	label VLDL = 'Very low density lipoprotein (mg/dL)';
	label WT = 'Weight (lbs)';
	Label FEMALE = 'Gender (Female)';
	label PROBAND = 'Proband?';
	label SBP = 'Systolic blood pressure (mmHg)';

	*REMOVE UNREALISTIC DATA OBSERVATIONS;
	
		*Guiness book world record oldest person 122 yo;
		if AGE < 0 or AGE > 122 then AGE = .;

		*Guiness book world record heaviest person 635 kg = 1400 lbs;
		if WT > 1400 then WT = .;

		*Highest recorded blood pressure in young adults preforming heavy, dynamic weight lifting: 370/360 Narloch (1995) et. al.;
		if SBP > 370 then SBP = .;
		if DBP > 360 then DBP = .;
	
	*make temporary id to use in merging proc means output;
	tempid=1;

run;
proc print data=stuff (obs=10);
run;


*==================================================================*
Check for and delete duplicates
*==================================================================*;
title;

proc sort data = stuff nodupkey dupout = dups;
	*identify variables to check for duplicates;
 	by &idvar &allvar;
run;
proc print data=dups (obs=20);
	title 'First 20 duplicate observations';
run;


*====================================================================================================================================*
*====================================================================================================================================*
		EXPLORATORY DATA ANALYSIS:
*====================================================================================================================================*
*====================================================================================================================================*;

*==================================================================*
Table 1

Continuous outcome and covariates normality assessment & summary statistics:
	n
	means
	standard deviations
	ranges

Categorical covariates:
	n
	frequencies
	percentages
	cumulative frequencies
	number missing values
*==================================================================*;
title;

*view means and summary statistics of Continuous outcome, exposure, and covariates;
proc means data=stuff maxdec=2;
	var &outco &expos &contcov;
	*class female; *good if there's a group to compare by;
	title 'Continuous outcome, expos and covariates means, std dev, ranges';
	title2 'Use for Table 1';
run;

*view frequencies of Categorical varibles;
proc freq;
	tables &catcov;
	title 'Categorical variables frequencies';
	title2 'Use for Table 1';
run;


*==================================================================*
Continuous outcome and covariates normality assessment & summary statistics:
	histogram, sideways
	boxplot
	QQ-plot
	histogram, proper orientation

Categorical covariates:
	boxplots
*==================================================================*;
title;


*assess linearity of outcome to see if need to transform;
proc univariate normal plots;
	var &outco;
	histogram;
	title "&outco (outcome) need to be transformed to satisfy linearity?";
	title2 "Look at extreme 5 observations to assess for impossible values";
	title3 "Normality? Does mean = median = mode?";
run;


*histograms of all (including Continuous) variables;
*tests and plots to assess nomality;
%macro univar(var);
	proc univariate normal plots;
		var &var;
		histogram;
		title "&var univariate analysis";
		title2 "Look at extreme 5 observations to assess for impossible values";
		title3 "Independent variables do NOT have to satisfy Normality (mean = median = mode)";
	run;
%mend univar;
%univar(&expos)
%univar(AGE);
%univar(BMI);
%univar(DBP);
%univar(HDL);
%univar(HT);
%univar(LDL);
%univar(SKIN);
%univar(TG);
%univar(VLDL);
%univar(WT);
%univar(FEMALE);
%univar(PROBAND);

*boxplots of Categorical variables;
%macro boxplot(outcount, catvar);
	proc sgplot;
		vbox &outcount / category=&catvar;
		title "Boxplots of &outcount by &catvar (Categorical variable)";
	proc sgplot;
		histogram &outcount / group=&catvar transparency=0.75;
		density &outcount / type=kernel group=&catvar;
		title "Comparative histograms of &outcount by &catvar (Categorical variable)";
	run;
%mend boxplot;
%boxplot(&outco,female)
%boxplot(&outco,proband)


*====================================================================================================================================*
*====================================================================================================================================*
		UNIVARIATE ANALYSIS:
*====================================================================================================================================*
*====================================================================================================================================*;

*==================================================================*
Center continuous variables on their means
*==================================================================*;
title;

*make means of all variables;
proc means data = stuff mean std n;
	by tempid; *all subjects have tempid=1;
	var &allvar;
	output out = stuffmeans mean= / autoname;
	title "Means of all variables for centering";
run;

*center all continuous independent variables on their means;
data stuff;
	merge stuff stuffmeans;
	by tempid;

	*CENTER VALUES BY MANUALLY INPUTTING VARIABLES AND THEIR MEANS FROM PROC MEANS ABOVE;
	CHOL_Cent 		= CHOL 		- CHOL_mean;
	AGE_Cent 		= AGE 		- AGE_mean;
	BMI_Cent 		= BMI 		- BMI_mean;
	DBP_Cent 		= DBP 		- DBP_mean;
	HDL_Cent 		= HDL 		- HDL_mean;
	HT_Cent 		= HT		- HT_mean;
	LDL_Cent 		= LDL 		- LDL_mean;
	SKIN_Cent 		= SKIN 		- SKIN_mean;
	TG_Cent 		= TG 		- TG_mean;
	VLDL_Cent 		= VLDL 		- VLDL_mean;
	WT_Cent 		= WT 		- WT_mean;

	*ASSIGN CENTERED VARIABLE NAMES;
	label CHOL_Cent = 'Total cholesterol (dg/mL), centered';
	label AGE_Cent = 'Age (years), centered';
	label BMI_Cent = 'Body mass index (kg/m^2), centered';
	label DBP_Cent = 'Diastolic blood pressure (mmHg), centered';
	label HDL_Cent = 'High density lipoprotein (mg/dL), centered';
	label HT_Cent = 'Height (in), centered';
	label LDL_Cent = 'Low density lipoprotein (mg/dL), centered';
	label SKIN_Cent = 'Blood glucose levels (mg/dL, centered)';
	label TG_Cent = 'Triglycerides (mg/dL), centered';
	label VLDL_Cent = 'Very low density lipoprotein (mg/dL), centered';
	label WT_Cent = 'Weight (lbs), centered';

	*TRANSFORMATIONS OF OUTCOME - IF NEEDED AFTER LOOKING AT:
		1) SCATTER OF OUTCOME VS PRIMARY EXPOSURE
		2) RESID VS. PRED PLOT AFTER FINAL MODEL);

		*Tukey's Ladder - outcome^lambda;
		*outco_log = log(&outco); 	*lambda = 0;
		*outco_sqrt = sqrt(&outco); *lambda = 1/2 (sqrt);
		*outco_2 = &outco*&outco;	*lambda = 2;

		*Note: In ESTIMATION, only consider transformations of the outcome after looking at the scatter of outcome vs main exposure 
		(not covariates). Covariates are not transformed in modeling!;

	/*add single spline term of primary exposure (at 100 in this example) - p.6-8;
	if &expos ne . then
		do;
			if &expos > 100 the expos_sp1 = &expos - 100;
			else expos_sp1 = 0;
		end;
	*/

*check to make sure dataset correct;
*proc print data=stuff (obs=10);
run;

*look at univariate analysis of transformed variable;
*%univar(&outco);
*%univar(OUTCO_LOG);

*<<<<<<<<<<<<<<<==================================================================*
 <<<<<<<<<<<<<<<	DEFINE GLOBAL MACROS FOR CENTERED VARIABLES:
 <<<<<<<<<<<<<<<		outcome, NOT centered
 <<<<<<<<<<<<<<<		exposure of interest, centered
 <<<<<<<<<<<<<<<		continuous covariates, centered
 <<<<<<<<<<<<<<<		categorical covariates, NOT centered
 <<<<<<<<<<<<<<<		all variables, centered plus outcome
*<<<<<<<<<<<<<<<==================================================================*;

%let outTrans	= &OUTCO; *change if transforming outcome above;
%let expCent	= CHOL_Cent;
%let contCent	= AGE_Cent BMI_Cent DBP_Cent HDL_Cent HT_Cent LDL_Cent SKIN_Cent TG_Cent VLDL_Cent WT_Cent;
*all variables (w/centered) except id variables;
%let allCent	= &outTrans &expCent &contcent &catcov;


*==================================================================*
Scatter plots of univariate relationships, w/macro
*==================================================================*;
title;

*scatter plots macro;
%macro scatthat(x,y,dataset);
	*simple scatter plot;
	proc sgplot data=&dataset;
		scatter x=&x y=&y / datalabel = &idvar;
		reg x=&x y=&y;
		title "Scatter plot of &y on &x";
	run;
%mend scatthat;

*main univariate relationship;
%scatthat(&expCent,&outTrans,stuff)
%scatthat(AGE_Cent,&outTrans,stuff)
%scatthat(BMI_Cent,&outTrans,stuff)
%scatthat(DBP_Cent,&outTrans,stuff)
%scatthat(HDL_Cent,&outTrans,stuff)
%scatthat(HT_Cent,&outTrans,stuff)
%scatthat(LDL_Cent,&outTrans,stuff)
%scatthat(SKIN_Cent,&outTrans,stuff)
%scatthat(TG_Cent,&outTrans,stuff)
%scatthat(VLDL_Cent,&outTrans,stuff)
%scatthat(WT_Cent,&outTrans,stuff)
%scatthat(FEMALE,&outTrans,stuff)
%scatthat(PROBAND,&outTrans,stuff)


*==================================================================*
Correlations between all independent variables and continuous outcome
	Look at top row only for interesting variables
*==================================================================*;
title;

*Pearsons correlations of all variables;
%macro corrthat(vars);
	proc corr data=stuff best = 6; *shows the top 6 correlations for each variable;
	  var &vars;
	  title "Pearson correlations: Percent variability shared = R^2 * 100%";
	  title2 "> 0.60-0.70 b/w exposure or covariates suggests collinearity";
	  *By squaring the correlation and then multiplying by 100%, 
	  you can determine what percentage of the variability is shared.  
	  Let’s round 0.59678 to be 0.6, which when squared would be .36, multiplied by 100 would be 36%.  
	  Hence read shares about 36% of its variability with write.;
	run;
%mend corrthat;
%corrthat(&allCent)

*EXPLAIN IF ANY VARIABLES WILL BE REMOVED DUE TO HIGH CORRELATION (I.E. THEY ARE COLLINEAR)
	REMOVE WEIGHT BECAUSE IT IS COLLINEAR WITH BMI (R = 0.867)
	REMOVE LDL BECAUSE IT IS COLLINEAR WITH CHOLESTEROL (R = 0.972)
	REMOVE TG BECAUSE IT IS COLLINEAR WITH VLDL (r = 0.911)


*==================================================================*
t-test for categorical independent variables and continuous outcome
	Should be the same p-value as for the univariate lin reg
*==================================================================*;
title;

%macro ttestthat(cat,out);
	proc ttest data=stuff;
		class &cat;
		var &out;
		title "T-test of &out and &cat";
	run;
%mend;

*PICK CATEGORICAL VARIABLE AND OUTCOME VARIABLE;
%ttestthat(FEMALE,&outTrans)
%ttestthat(PROBAND,&outTrans)


*==================================================================*
One-way ANOVA for categorical variables and outcome
	Useful if the exposure of interest or covariate is multilevel categorical
		If varible is unordered categorical, should get same p-value as for univariate lin reg
*==================================================================*;
title;

%macro anovathat(cat,out);
	proc glm data=stuff;
		class &cat;
		model &out=&cat;
		means &cat;
		title "1-Way ANOVA of &out and &cat";
	run;
	quit;
%mend anovathat;

%anovathat(FEMALE,&outTrans)
%anovathat(PROBAND,&outTrans)


*==================================================================*
Univariate linear regressions of outcome on each independent variable 
*==================================================================*;
title;

*produces ALL OUTPUT of proc glm;
%macro unireglong(out,var);
	proc glm data=stuff;
		model &out=&var / clparm;
		title "Univariate linear regression of &out on &var";
	run;
	quit;
%mend unireglong;
%unireglong(&outTrans,&expCent)
%unireglong(&outTrans,AGE_Cent)
%unireglong(&outTrans,BMI_Cent)
%unireglong(&outTrans,DBP_Cent)
%unireglong(&outTrans,HDL_Cent)
%unireglong(&outTrans,HT_Cent)
%unireglong(&outTrans,SKIN_Cent)
%unireglong(&outTrans,TG_Cent)
%unireglong(&outTrans,VLDL_Cent)
%unireglong(&outTrans,FEMALE)
%unireglong(&outTrans,PROBAND)

*shows ONLY PARAMETER ESTIMATES of each regression;
%macro uniregshort(out,var);
	ods select ParameterEstimates;
	proc glm data=stuff;
		model &out=&var / clparm;
		title "Univariate linear regression of &out on &var";
		title2 "Parameter estimates only for Table 2";
	run;
	quit;
	ods select all;
%mend uniregshort;
%uniregshort(&outTrans,&expCent)
%uniregshort(&outTrans,AGE_Cent)
%uniregshort(&outTrans,BMI_Cent)
%uniregshort(&outTrans,DBP_Cent)
%uniregshort(&outTrans,HDL_Cent)
%uniregshort(&outTrans,HT_Cent)
%uniregshort(&outTrans,SKIN_Cent)
%uniregshort(&outTrans,TG_Cent)
%uniregshort(&outTrans,VLDL_Cent)
%uniregshort(&outTrans,FEMALE)
%uniregshort(&outTrans,PROBAND)


*====================================================================================================================================*
*====================================================================================================================================*
		MULTIVARIATE ANALYSIS:
*====================================================================================================================================*
*====================================================================================================================================*;

*==================================================================*
Assess Confounding
	Check if adding each independent variable into model of exposure on outcome
	Changes the univariate beta estimate by 10-15% 
*==================================================================*;
title;

*univariate linear regression of outcome on exposure of interest;
%uniregshort(&outTrans,&expCent)

*multivariate regression of outcome on exposure of interest, adjusted for single covariate;
*shows ONLY PARAMETER ESTIMATES of each regression;
%macro confoundregshort(out,var,adj);
	ods select ParameterEstimates;
	proc glm data=stuff;
		model &out=&var &adj / clparm;
		title "Univariate linear regression of &out on &var, adjusted for &adj";
		title2 "Parameter estimates only for Table 3";
	run;
	quit;
	ods select all;
%mend confoundregshort;

*unadjusted, base model;
%confoundregshort(&outco,&expCent)

*model adjusted for each covariate;
%confoundregshort(&outTrans,&expCent,AGE_cent)
%confoundregshort(&outTrans,&expCent,BMI)
%confoundregshort(&outTrans,&expCent,DBP_Cent)
%confoundregshort(&outTrans,&expCent,HDL_Cent)
%confoundregshort(&outTrans,&expCent,HT_Cent)
%confoundregshort(&outTrans,&expCent,SKIN_Cent)
%confoundregshort(&outTrans,&expCent,VLDL_Cent)
%confoundregshort(&outTrans,&expCent,FEMALE)
%confoundregshort(&outTrans,&expCent,PROBAND)

*look at Table 3 to assess which variables are potential confounders
 and may be used in the preliminary main effects model;


*<<<<<<<<<<<<<<<==================================================================*
 <<<<<<<<<<<<<<<	DEFINE GLOBAL MACROS FOR VARIABLES IN PRELIMINARY MAIN EFFECTS (PME) MODEL:
 <<<<<<<<<<<<<<<		outcome, no change
 <<<<<<<<<<<<<<<		exposure of interest, no change
 <<<<<<<<<<<<<<<		all covariates used in the PME model
 <<<<<<<<<<<<<<<		all variables used in PME model
*<<<<<<<<<<<<<<<==================================================================*;

*all covariates in preliminary main effects model;
%let PMEcov	= AGE_Cent FEMALE BMI_Cent DBP_Cent VLDL_Cent;
	*NOTE: Even though Age or Gender may not be confounders by this definition, 
	 they are both included in the model because the literature suggests these variables are important 
	 in this exposure-disease association;
	*NOTE2: Weight and LDL were removed
	 because they were collinear with BMI and LDL, respectively;


*all variables (w/centered) except id variables;
%let allPME	= &outTrans &expCent &PMEcov;


*==================================================================*
Preliminary main effects model
	Outcome on primary exposure, adjusted for each confounder
	Model refinement: LINE assumptions
*==================================================================*;
title;

*macro to produce the Prelminary Main Effects model;
%macro premaineffects(out,var,adj);
	proc glm data=stuff;
		model &out=&var &adj / clparm;
		title "Preliminary main effects model of &out on &var,";
		title2 "Adjusted for &adj";
		title3 "Do NOT formally assess LINE assumptions until after final model selected";
		output out=premaineffects p=pred r=resid;
	run;
	quit;
	proc sgplot data=premaineffects;
		scatter x=pred y=resid;
		loess  x=pred y=resid / clm smooth=0.4;
		refline 0;
		title "Quick check of Linearity:"; 
		title2 "Quick check of Homoscedasticity: Are Variances (of residuals) Equal?";
		title3 "Residuals vs. Predicted &out to look for pattern (Loess) or unequal variance";
	proc univariate data=premaineffects normal plots;
		var resid;
		title "Are Residuals Normally Distributed?";
		title2 "Quick check of Normality of &out on &var,";
		title3 "Adjusted for &adj";
	run;
%mend premaineffects;

*ADD IN ALL THE COVARIATES THAT SHOWED CONFOUNDING;
%premaineffects(&outTrans,&expCent,&PMEcov)

*QUICKLY CHECK INTERPRETATIONS OF LINE ASSUMPTIONS AND TOLERANCE (COLLINEARITY) DIAGNOSTIC
	1) LINEARITY IS NOT SATISFIED, LOOKING AT ONE POINT TOWARDS THE EXTREME PREDICTED VALUE
	2) INDEPENDENCE IS ASSUMED BECAUSE EACH SUBJECT IS THEIR OWN. THIS MAY BE AN INAPPROPRIATE ASSUMPTION BECAUSE MANY SUBJECTS ARE RELATED
	3) NORMAILITY IS ASSUMED BY LOOKING AT THE GRAPHS OF THE RESIDUALS AND THE QQ PLOT
	4) EQUAL VARIANCES/HOMOSCEDASTICITY IS ASSUMED BECAUSE THE RESIDUALS SEEM TO BE SPREAD RELATIVELY EVENLY THROUGHOUT PREDICTED VALUES
;

*==================================================================*
Assess Effect Modification
	Check if adding each interaction term to preliminary multivariate model
	Produces a statistically significant beta estimate (for interaction term)
*==================================================================*;
title;

*preliminary multivariate regression of outcome on exposure of interest, adjusted for all confounders;
*shows ONLY PARAMETER ESTIMATES of each regression;
%macro interactionRegShort(out,var,adj,ixn);
	ods select ParameterEstimates;
	proc glm data=stuff;
		model &out=&var &adj  &ixn / clparm;
		title1 "Assessing interaction term &ixn";
		title2 "Multivariate linear regression of &out on &var, adjusted for &adj";
		title3 "Parameter estimates only for Table 4";
	run;
	quit;
	ods select all;
%mend interactionRegShort;
%interactionRegShort(&outTrans,&expCent,&PMEcov,&expCent*AGE_CENT)
%interactionRegShort(&outTrans,&expCent,&PMEcov,&expCent*BMI_CENT)
%interactionRegShort(&outTrans,&expCent,&PMEcov,&expCent*DBP_CENT)
%interactionRegShort(&outTrans,&expCent,&PMEcov,&expCent*VLDL_CENT)
%interactionRegShort(&outTrans,&expCent,&PMEcov,&expCent*FEMALE)

*look at Table 4 to assess which interaction terms are significant
 and may be used in the Preliminary Final model;

*<<<<<<<<<<<<<<<	SELECT ALL SIGNIFICANT INTERACTION TERMS TO POTENTIALLY INCLUDE IN PREMINIARY FINAL MODEL;
%let ixnAllSig = &expCent*AGE_CENT &expCent*BMI_CENT &expCent*DBP_CENT &expCent*VLDL_CENT;

*Model that contains all significant interaction terms;
*check to see if all the interaction terms maintain significance;
%interactionRegShort(&outTrans,&expCent,&PMEcov,&ixnAllSig)

*remove &expCent*BMI_CENT term (p = 0.799);
%interactionRegShort(&outTrans,&expCent,&PMEcov,&expCent*AGE_CENT &expCent*DBP_CENT &expCent*VLDL_CENT)

*remove &expCentVLDL*_CENT term (p = 0.488);
%interactionRegShort(&outTrans,&expCent,&PMEcov,&expCent*AGE_CENT &expCent*DBP_CENT)


*<<<<<<<<<<<<<<<==================================================================*
 <<<<<<<<<<<<<<<	DEFINING GLOBAL MACROS FOR IXN VARIABLES IN PRELIMINARY FINAL (PF) MODEL:
 <<<<<<<<<<<<<<<		outcome, no change
 <<<<<<<<<<<<<<<		exposure of interest, no change
 <<<<<<<<<<<<<<<		all covariates used in the PF model, no change
 <<<<<<<<<<<<<<<		all significant interaction terms used in PF model
 <<<<<<<<<<<<<<<		all variables used in PF, transformed
*<<<<<<<<<<<<<<<==================================================================*;

*SELECT FINAL INTERACTION TERMS THAT REMAIN SIGNIFICANT IN MODEL;
%let ixnFin = &expCent*AGE_Cent &expCent*DBP_CENT;

*all variables (w/centered) except id variables;
%let allPF	= &outTrans &expCent &PMEcov &ixnFin;


*==================================================================*
Preliminary Final Model
	Assess LINE assumptions
	Model Diagnostics: Tolerance, Cooks D, Leverage
*==================================================================*;
title;


*macro to produce the Prelminary Main Effects model;
%macro prelimFinalModelReg(out,var,adj,ixn,ident);
	proc glm data=stuff plot=diagnostics(label);
		id &ident;
		model &out=&var &adj &ixn / clparm tolerance;
		title "Preliminary final model of &out on &var,";
		title2 "Adjusted for &adj, w/interaction terms &ixnFin";
		title3 "Tolerance < 0.1 indicates substantial collinearity";
		output out=prefinalmodel p=Predicted r=Residual student=StudentResid cookd=CooksD 
				h=Leverage rstudent=JackknifeResid dffits=DFFITs covratio=CovRatio;
	run;
	quit;
	proc sgplot data=prefinalmodel;
		scatter x=predicted y=residual;
		loess  x=predicted y=residual / clm smooth=0.4;
		refline 0;
		title "Assess Linearity:"; 
		title2 "Assess Homoscedasticity: Are Variances (of residuals) Equal?";
		title3 "Residuals vs. Predicted &out to look for pattern (Loess) or unequal variance";
	proc univariate data=prefinalmodel normal plots;
		var residual;
		title "Are Residuals Normally Distributed?";
		title2 "Assess Normality of &out on &var,";
		title3 "Adjusted for &adj, w/interaction terms &ixnFin";
	run;
%mend prelimFinalModelReg;

*Preliminary Final Model that contains all interaction terms that maintain significance;
%prelimFinalModelReg(&outTrans,&expCent,&PMEcov,&ixnFin,&idvar)

*WRITE OUT INTERPRETATIONS OF LINE ASSUMPTIONS AND TOLERANCE (COLLINEARITY) DIAGNOSTIC
	1) LINEARITY IS NOT SATISFIED, LOOKING AT ONE POINT TOWARDS THE EXTREME PREDICTED VALUE
	2) INDEPENDENCE IS ASSUMED BECAUSE EACH SUBJECT IS THEIR OWN. THIS MAY BE AN INAPPROPRIATE ASSUMPTION BECAUSE MANY SUBJECTS ARE RELATED
	3) NORMAILITY IS ASSUMED BY LOOKING AT THE GRAPHS OF THE RESIDUALS AND THE QQ PLOT
	4) EQUAL VARIANCES/HOMOSCEDASTICITY IS ASSUMED BECAUSE THE RESIDUALS SEEM TO BE SPREAD RELATIVELY EVENLY THROUGHOUT PREDICTED VALUES

*look graphically for outliers;
	*diagnostics;
		%scatthat(leverage,predicted,prefinalmodel)
		%scatthat(cooksd,predicted,prefinalmodel)
		%scatthat(DFFITS,predicted,prefinalmodel)
	*scatter plot of residuals vs predicted;
		%scatthat(residual,predicted,prefinalmodel)
		%scatthat(studentResid,predicted,prefinalmodel)
		%scatthat(jackknifeResid,predicted,prefinalmodel)

*assess most extreme diagnostic values;
%macro extremeDiag(diag);
	data diagnosticPFM;
		set prefinalmodel;
		keep id &outTrans &expCent predicted residual &diag;
	proc sort data=diagnosticPFM;
		by descending &diag;
	proc print data=diagnosticPFM(obs = 25);
		title "&diag: 25 most extreme observations";
		title2 "n = # observations ; k = # model parameters";
		title3 "Leverage > 2(k+1)/n	; CooksD > 1 ; DFFITs > 2*sqrt(k/n)";
run;
%mend;
*diagnostics;
	%extremeDiag(leverage)
	%extremeDiag(cooksd)
	%extremeDiag(DFFITS)
*residuals;
*	%extremeDiag(residual);
*	%extremeDiag(studentResid);
*	%extremeDiag(jackknifeResid);

*extreme predicted values;
%extremeDiag(predicted)

*look for collinearity in the variables in the preliminary final model;
%corrthat(&outTrans &expCent &PMEcov)


*==================================================================*
Sensitivity analysis

*==================================================================*;
title;

*macro to produce the Prelminary Main Effects model;
*PUT THE ID NUMBERS OF THE OUTLIERS IN THE GREEN VALUES FOR 42 AND 115;
%macro SensitivityAnalysis(out,var,adj,ixn,ident);
	proc glm data=stuff(where=(&ident ne 42 and &ident ne 115)) plot=diagnostics(label);
		id &ident;
		model &out=&var &adj &ixn / clparm;
		title "Sensitivity analysis of Final Model (&out on &var), w/adjusting vars and ixns";
		title2 "Eliminating potentially influential points";
		output out=SensitivityAnalysis p=Predicted r=Residual student=StudentResid cookd=CooksD 
				h=Leverage rstudent=JackknifeResid dffits=DFFITs covratio=CovRatio;
	run;
	quit;
	proc sgplot data=SensitivityAnalysis;
		scatter x=predicted y=residual;
		loess  x=predicted y=residual / clm smooth=0.4;
		refline 0;
		title "Assess Linearity:"; 
		title2 "Assess Homoscedasticity: Are Variances (of residuals) Equal?";
		title3 "Residuals vs. Predicted &out to look for pattern (Loess) or unequal variance";
	proc univariate data=SensitivityAnalysis normal plots;
		var residual;
		title "Are Residuals Normally Distributed?";
		title2 "Assess Normality of &out on &var,";
		title3 "Adjusted for &adj, w/interaction terms &ixn";
	run;
%mend SensitivityAnalysis;

*Sensitivity Analysis that contains all interaction terms that maintain significance;
%SensitivityAnalysis(&outTrans,&expCent,&PMEcov,&ixnFin,&idvar)

*NOTES TO SELF:
4) Write up plan for prediction modeling

*LESSONS
3) After assessing PME, if Residuals vs Prections plot doesnt look good, try looking at residuals vs each independent variable separately
	Allows you to parse out which exposure/covariate is causing the violation in the overall model.
	If one of the x variables looks guilty, you can try to either
		a)transform outcome - good for increasing spread with increasing x variable
		b)add another term for the guilty x variable (x^2 or x^3) - HARD!!! (p.524)
4) Effect modification, look for all significant interactions individually,
	Include all significant interactions in the model, and reassess each term for significance
	Keep all interaction terms in the model if they maintain significance (alpha < 0.05) or marginal significance (alpha ~0.05)
5) Check LINE assumptions
6) Check model diagnostics
7) Call it a day
;

