
//UTILITIES

function vSum(arr) {
	/**
	Calculates the sum of elements in an array.
	**/
	var sum = 0.0;
	for (var i = 0; i < arr.length; i++) {
	  sum += arr[i];
	}
	return sum;
}

function vLog(arr) {
	/**
	Takes in an array of numbers and returns an array of the natural log of
	each number.
	**/
	for (var i = 0; i < arr.length; i++) {
		arr[i] = Math.log(arr[i]);
	}
	return arr;
}

function argMax(array) {
	/**
	Retrieve the array key corresponding to the largest element in the array.
	**/
	return array.map((x, i) => [x, i]).reduce((r, a) => (a[0] > r[0] ? a : r))[1];
}

function arange(start, end, step=1){
  	var rangeArray = [];
  	if(start < end){
	    for(var i=start; i<end; ){
	    	rangeArray.push(i);
	    	i = i + step;
	    }
	}
	return rangeArray;
}


function fisher_combined_pvalue(pvalues) {
	/** 
	Find the p-value for Fisher's combined test statistic

    Parameters
    ----------
    pvalues : array_like
        Array of p-values to combine

    Returns
    -------
    float
    	p-value for Fisher's combined test statistic

	**/
	if (pvalues.includes(0)) { 
	// check empty vals
		return 0;
	} else {
		obs = -2 * vSum(vLog(pvalues));
	}
	return 1 - stat_functions.chi2Cdf(df=2*pvalues.length, x=obs); //What is rmin in incGamma?
}





function maximize_fisher_combined_pvalue(N, overall_margin, pvalue_funs, precise=true, plausible_lambda_range=null){
    /**
    Find the smallest Fisher's combined statistic for p-values obtained 
    by testing two null hypotheses at level alpha using data X=(X1, X2) 

    Parameters
    ----------
    N : array_like
        Array of stratum sizes
    overall_margin : int
        the difference in votes for reported winner and loser in the population
    pvalue_funs : array_like
        functions for computing p-values. The observed statistics/sample and known parameters should be
        plugged in already. The function should take the lambda allocation AS INPUT and output a p-value.
    precise : bool
        Optional, should we refine the maximum found by minimize_scalar? Default is True
    plausible_lambda_range : array-like
        lower and upper limits to search over lambda. Optional, but will speed up the search
    
    Returns
    -------
    dict with 
    
    float
        maximum combined p-value
    float
        minimum value of Fisher's combined test statistic
    float
        lambda, the parameter that minimizes the Fisher's combined statistic/maximizes the combined p-value
    **/

    //assert len(N)==2
    //assert len(pvalue_funs)==2

    if (plausible_lambda_range == null) {
    	lambda_upper = int(Math.min([2*N[0]/overall_margin, 1+2*N[1]/overall_margin])) + 1;
    	lambda_lower = int(Math.max([-2*N[0]/overall_margin, 1-2*N[1]/overall_margin]));
    }
    else {
    	lambda_lower = plausible_lambda_range[0];
    	lambda_upper = plausible_lambda_range[1];
	}

	fisher_pvalues = [];
    cvr_pvalues = [];
    test_lambdas = arange(lambda_lower, lambda_upper+1, 0.5); // fix
    for (var lam of test_lambdas) {
    	pvalue = Math.min([1, pvalue_funs[0](lam)]);
    	if (pvalue1 < 0.01) {
    		fisher_pvalues.push(0);
    	} else {
    		pvalue2 = Math.min([1, pvalue_funs[1](1-lam)]);
    		fisher_pvalues.push(fisher_combined_pvalue([pvalue1,pvalue2]));
    	}
    }
    pvalue = Math.max(fisher_pvalues);
    alloc_lambda = test_lambdas[argMax(fisher_pvalues)];

    if (precise) {
        fisher_pvalues = []
        test_lambdas = arange(alloc_lambda-0.5, alloc_lambda+0.5, 0.1); //fix
    	for (var lam of test_lambdas) {
            pvalue1 = Math.min([1, pvalue_funs[0](lam)]);
            pvalue2 = Math.min([1, pvalue_funs[1](1-lam)]);
            fisher_pvalues.push(fisher_combined_pvalue([pvalue1, pvalue2]));
     
	        if (Math.max(fisher_pvalues) > pvalue) {
	            pvalue = Math.max(fisher_pvalues);
	            alloc_lambda = test_lambdas[argMax(fisher_pvalues)];
	        }
	    }
    }
    var dict = {
    		'max_pvalue' : pvalue,
            'min_chisq' : stat_functions.chi2Inv(p=1 - pvalue, df=4), // assuming chi2inv is the ppf?
            'allocation lambda' : alloc_lambda
            };
    return dict;
}

