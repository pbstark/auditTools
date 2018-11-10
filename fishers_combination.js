
//UTILITIES

function getRandomSubarray(arr, size) {
    var shuffled = arr.slice(0), i = arr.length, temp, index;
    while (i--) {
        index = Math.floor((i + 1) * Math.random());
        temp = shuffled[index];
        shuffled[index] = shuffled[i];
        shuffled[i] = temp;
    }
    return shuffled.slice(0, size);
}

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

function vMultiply(constant, arr) {
    var return_arr = []
    for (var i = 0; i < arr.length; i++) {
	  return_arr.push(constant * arr[i]);
	}
	return return_arr;
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

function argMin(array) {
	/**
	Retrieve the array key corresponding to the largest element in the array.
	**/
	return array.map((x, i) => [x, i]).reduce((r, a) => (a[0] < r[0] ? a : r))[1];
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


function create_modulus(n1, n2, n_w2, n_l2, N1, V_wl, gamma) {
    /**
    Find the smallest Fisher's combined statistic for p-values obtained 	    The modulus of continuity for the Fisher's combined p-value.
    by testing two null hypotheses at level alpha using data X=(X1, X2) 	    This function returns the modulus of continuity, as a function of
    the distance between two lambda values.

    n1 : int
        sample size in the ballot comparison stratum
    n2 : int
        sample size in the ballot polling stratum
    n_w2 : int
        votes for the reported winner in the ballot polling sample
    n_l2 : int
        votes for the reported loser in the ballot polling sample
    N1 : int
        total number of votes in the ballot comparison stratum
    V_wl : int
        margin (in votes) between w and l in the whole contest
    gamma : float
        gamma from the ballot comparison audit
    **/

    var Wn = n_w2;
    var Ln = n_l2;
    var Un = n2-n_w2-n_l2

    function return_func(delta) {

        var return_arr = []

        var arr1 = vMultiply(2*Wn, vLog(1 + V_wl*delta/2));
        var arr2 = vMultiply(2*Ln, vLog(1 + V_wl*delta/2));
        var arr3 = vMultiply(2*Un, vLog(1 + 2*V_wl*delta));
        var arr4 = vMultiply(2*n1, vLog(1 + V_wl*delta/(2*N1*gamma)));

        for (var i = 0; i < arr1.length; i++) {
          return_arr[i] = arr1[i] + arr2[i] + arr3[i] + arr4[i];
        }

        return return_arr;


    }

    return return_func;

}


function maximize_fisher_combined_pvalue(N_w1, N_l1, N1, N_w2, N_l2, N2,
    pvalue_funs, stepsize=0.05, modulus=None, alpha=0.05, feasible_lambda_range=None) {
    /**
    Grid search to find the maximum P-value.

    Find the smallest Fisher's combined statistic for P-values obtained
    by testing two null hypotheses at level alpha using data X=(X1, X2).
    Parameters
    ----------
    N_w1 : int
        votes for the reported winner in the ballot comparison stratum
    N_l1 : int
        votes for the reported loser in the ballot comparison stratum
    N1 : int
        total number of votes in the ballot comparison stratum
    N_w2 : int
        votes for the reported winner in the ballot polling stratum
    N_l2 : int
        votes for the reported loser in the ballot polling stratum
    N2 : int
        total number of votes in the ballot polling stratum
    pvalue_funs : array_like
        functions for computing p-values. The observed statistics/sample and known parameters should be
        plugged in already. The function should take the lambda allocation AS INPUT and output a p-value.
    stepsize : float
        size of the grid for searching over lambda. Default is 0.05
    modulus : function
        the modulus of continuity of the Fisher's combination function.
        This should be created using `create_modulus`.
        Optional (Default is None), but increases the precision of the grid search.
    alpha : float
        Risk limit. Default is 0.05.
    feasible_lambda_range : array-like
        lower and upper limits to search over lambda.
        Optional, but a smaller interval will speed up the search.

    Returns
    -------
    dict with

    max_pvalue: float
        maximum combined p-value
    min_chisq: float
        minimum value of Fisher's combined test statistic
    allocation lambda : float
        the parameter that minimizes the Fisher's combined statistic/maximizes the combined p-value
    refined : bool
        was the grid search refined after the first pass?
    stepsize : float
        the final grid step size used
    tol : float
        if refined is True, this is an upper bound on potential approximation error of min_chisq
    **/



    if (plausible_lambda_range == null) {
    	feasible_lambda_range = calculate_lambda_range(N_w1, N_l1, N1, N_w2, N_l2, N2)
    }

    var lambda_lower = feasible_lambda_range[0];
    var lambda_upper = feasible_lambda_range[1];

    var test_lambdas = arange(lambda_lower, lambda_upper+stepsize, stepsize);
    if test_lambdas.length < 5 {
        stepsize = (lambda_upper + 1 - lambda_lower)/5;
        test_lambdas = arange(lambda_lower, lambda_upper+stepsize, stepsize);
    }


	var fisher_pvalues = [];
    for (var lam of test_lambdas) {
    	pvalue1 = Math.min([1, pvalue_funs[0](lam)]);
    	pvalue2 = Math.min([1, pvalue_funs[1](1-lam)]);
    	fisher_pvalues.push(fisher_combined_pvalue([pvalue1,pvalue2]));
    }

    var pvalue = Math.max(fisher_pvalues);
    var alloc_lambda = test_lambdas[argMax(fisher_pvalues)];


    // If p-value is over the risk limit, then there's no need to refine the
    // maximization. We have a lower bound on the maximum.
    if (pvalue > alpha || modulus is null) {
        return {'max_pvalue' : pvalue,
                'min_chisq' : stat_functions.chi2Pdf(1 - pvalue, df=4),
                'allocation lambda' : alloc_lambda,
                'tol' : None,
                'stepsize' : stepsize,
                'refined' : False
                }
    }


    // Use modulus of continuity for the Fisher combination function to check
    // how close this is to the true max
    var fisher_fun_obs = stat_functions.chi2Pdf(1-pvalue, df=4);
    var fisher_fun_alpha = stat_functions.chi2Pdf(1-alpha, df=4);
    var dist = Math.abs(fisher_fun_obs - fisher_fun_alpha)
    var mod = modulus(stepsize)

    if mod <= dist {
        return {'max_pvalue' : pvalue,
                'min_chisq' : fisher_fun_obs,
                'allocation lambda' : alloc_lambda,
                'stepsize' : stepsize,
                'tol' : mod,
                'refined' : False
                }
    } else {
        var lambda_lower = alloc_lambda - 2*stepsize
        var lambda_upper = alloc_lambda + 2*stepsize
        var refined = maximize_fisher_combined_pvalue(N_w1, N_l1, N1, N_w2, N_l2, N2,
            pvalue_funs, stepsize=stepsize/10, modulus=modulus, alpha=alpha,
            feasible_lambda_range=(lambda_lower, lambda_upper))
        refined['refined'] = True
        return refined
    }

};

function simulate_fisher_combined_audit(N_w1, N_l1, N1, N_w2, N_l2, N2, n1, n2, alpha,
    reps=10000, verbose=false, plausible_lambda_range=null) {


    var margin = (N_w1+N_w2)-(N_l1+N_l2);
    var N1 = N_w1+N_l1;
    var N2 = N_w2+N_l2;
    var Vwl = (N_w1 + N_w2) - (N_l1 + N_l2)
    var pop2 = new Array(N2).fill(0);


    for (var i = 0; i < N_w2; i++) {
        pop2[i] = 1;
    };


    for (var i = N_w2 + N_l2; i < N2; i++) {
        pop2[i] = NaN;
    };


    function cvr_pvalue(alloc) {
        return ballot_comparison_pvalue(n=n1, gamma=1.03905, o1=0, u1=0, o2=0, u2=0, reported_margin=margin, N=N1, null_lambda=alloc)
    };


    var fisher_pvalues = new Array(reps).fill(0);




    for (var i = 0; i < reps; i++) {
        if (verbose) {
            console.log(i)
        }

        var sam = getRandomSubarray(pop2, n2);
        var nw2 = 0;
        var nl2 = 0;
        for (var i = 0; i < sam.length; i++ ){
            if (sam[i] == 1) {
                nw2 += 1;
            } else if (sam[i] == 1) {
                nl2 += 1;
            }
        };
        var mod = create_modulus(n1, n2, nw2, nl2, N1, Vwl, 1.03905)


        function nocvr_pvalue(alloc) {
            return ballot_polling_sprt(sam, N2, alpha, N_w2, N_l2, (N_w2-N_l2) - alloc*margin)['pvalue']
        };

        fisher_pvalues[i] = maximize_fisher_combined_pvalue(N_w1, N_l1,
                               N1, N_w2, N_l2, N2,
                               [cvr_pvalue, nocvr_pvalue],
                               mod,
                               plausible_lambda_range)['max_pvalue']
    };

    var sum = 0;

    for(var i = 0; i < fisher_pvalues.length; i++ ){
        if (fisher_pvalues[i] <= alpha) {
            sum += 1
        }
    };

    var avg = sum/fisher_pvalues.length;
    return avg
};


function calculate_lambda_range(N_w1, N_l1, N1, N_w2, N_l2, N2) {
    var V1 = N_w1 - N_l1;
    var V2 = N_w2 - N_l2;
    var V = V1+V2;
    var lb = argMin([2*N1/V, 1+2*N2/V,1-(N2+V2)/V])
    var ub = argMax([-2*N1/V, 1-2*N2/V,  1+(N2-V2)/V])
    return (lb, ub)
};

function bound_fisher_fun(N_w1, N_l1, N1, N_w2, N_l2, N2,
                     pvalue_funs, plausible_lambda_range=null, stepsize=0.5) {

        if plausible_lambda_range is null:
            plausible_lambda_range = calculate_lambda_range(N_w1, N_l1, N1, N_w2, N_l2, N2)
        var lambda_lower = plausible_lambda_range[0]
        var lambda_upper = plausible_lambda_range[1]

        var cvr_pvalue = pvalue_funs[0]
        var nocvr_pvalue = pvalue_funs[1]
        var cvr_pvalues = []
        var nocvr_pvalues = []
        var lambdas = []

        var lam = lambda_lower

        while (lam< lambda_upper + 1) {
            var pvalue1 = argMin([1, cvr_pvalue(lam)])
            var pvalue2 = argMin([1, nocvr_pvalue(1-lam)])
            cvr_pvalues.push(pvalue1)
            nocvr_pvalues.push(pvalue2)
            lambdas.push(lam)
            lam += stepsize
        }

        lower_bounds = []
        upper_bounds = []

        for (var i = 0; i < cvr_pvalues - 1; i++) {
            lower_bounds.push(fisher_combined_pvalue([cvr_pvalues[i+1], nocvr_pvalues[i]]))
            upper_bounds.push(fisher_combined_pvalue([cvr_pvalues[i], nocvr_pvalues[i+1]]))
        }

        sample_points = []

        for (var i = 0; i < cvr_pvalues - 1; i++) {
            sample_points.push(fisher_combined_pvalue([cvr_pvalues[i], nocvr_pvalues[i]]))
        }

        return {'sample_points' : sample_points,
                'upper_bounds' : upper_bounds,
                'lower_bounds' : lower_bounds,
                'grid' : lambdas
                }

}
