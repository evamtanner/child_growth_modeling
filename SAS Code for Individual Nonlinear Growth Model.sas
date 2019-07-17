/*******************************************************************************************************
SAS Code to Implement Double Logistic Model for Individual Characterizing Child Growth Trajectories
** Starting values specified for infant and child growth
** Different model confidgurations & starting values necessary depending on available growth data
*******************************************************************************************************/

/*Create more data points for plotting predicted values*/
data one_obs;
	set growth_data_long; /*data in long format*/
	if age_mo = 0; /*select each child once*/
	keep id; 
	run;
data for_plot;
	set one_obs;
	do age_mo = 0.5 to 66.5 by 2; /*unobserved ages to predict weight*/
		output;
	end; 
	keep id age_mo;
	run;
data nlin_1; /*label since may need multiple runs*/
	set growth_data_long for_plot;
	run;
proc sort data=nlin_1; 
by id age_mo; 
run;

/*Double logistic model*/
proc nlin data=nlin_1 noitprint maxiter=500 outest=betas_1 /*obtain SSE*/;
	by id; 
	/*Specify starting values in kg and months*/;
	parms a1 = 6 to 12 by 0.5 /*1st weight plateau between 6-12 kg*/  
		  b1 = 0.1 to 1 by 0.1 /*infant velocity should be near 9 kg per 12 months*/ 
		  del1 = 1 to 11 by 1 /*infant peak velocity before 1 year*/
          b2 = 0.1 to 1 by 0.1 /*slower mid-childhood velocity*/  
		  del2 = 10 to 100 by 20; /*unknown mid-childhood velocity timing*/  
   	H=40; /*fixed mid-childhood weight plateau to reduce #estimated parameters*/
   	t=age_mo; 
	/*double logistic model*/
   	mu = a1/(1+exp(-b1*(t-del1))) + (H-a1)/(1+exp(-b2*(t-del2)));
	/*first derivative for velocity*/
	dmu_dt = -a1*(1+exp(-b1*(t-del1)))**(-2)*exp(-b1*(t-del1))*(-b1) + 
        	   -(H-a1)*(1+exp(-b2*(t-del2)))**(-2)*exp(-b2*(t-del2))*(-b2); 
	model weight_kg = mu;
    id dmu_dt; /*reconstruct weight velocities*/
   	output out=pred_1 p=pred; /*weight predictions*/
	ods output ParameterEstimates=param_1 /*growth parameter estimates*/
     		   ConvergenceStatus=conv_1; /*log of convergence info*/
   	run; quit;
data sse_1; /*=500 iterations so delete to free computer space*/
	set betas_1; 
	if _TYPE_ ^= 'FINAL' then delete; 
	keep id _SSE_; 
	run;

/*Plot Individual Trajectories*/
proc sgplot data=pred_1;
   by id;
   scatter y=weight_kg x=age_mo / markerattrs=(symbol=circlefilled color=black size=10);
   series y=pred x=age_mo / lineattrs=(pattern=solid thickness=3 color=black);
   series y=dmu_dt x=age_mo / y2axis lineattrs=(pattern=shortdash thickness=3 color=black);
   yaxis minor label="Weight (kg)";
   y2axis minor label="Weight Velocity (kg/month)";
   xaxis minor label="Age (months)";
   run; quit;	
