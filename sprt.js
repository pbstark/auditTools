function sprt(sample, popsize, alpha, Vw, Vl, null_margin=0) {
	// Set parameters
	var beta, lower, upper, n, Wn, Ln, Un, decision
	beta = 0;
	lower = beta / (1 - alpha);
	upper = (1 - beta) / alpha;
	n = sample.length;
	Wn = 0;
	Ln = 0;
	for (i = 0; i < n; i++) {
		if (sample[i] == 1) {
			Wn += 1;
		} else if (sample[i] == 0) {
			Ln += 1;
		}
	}
	Un = n - Wn - Ln;
	decision = "None";


	// Set up likelihood for null and alternative hypotheses
	var Vw, Vl, Vu
	Vw = int(Vw);
	Vl = int(Vl);
	Vu = int(popsize - Vw - Vl);
	if (!(Vw >= Wn || Vl >= Ln || Vu >= Un)) {
		document.write("Alternative hypothesis isn't consistent with sample.")
	}


	var alt_W, alt_L, alt_U, alt_logLR;
	alt_W = 0;
	alt_L = 0;
	alt_U = 0;
	for (i = 0; i < Wn; i++) {
		alt_W += Math.log(Vw - i);
	}
	for (i = 0; i < Ln; i++) {
		alt_L += Math.log(Vl - i);
	}
	for (i = 0; i < Un; i++) {
		alt_U += Math.log(Vu - i);
	}
	alt_logLR = alt_W + alt_L + alt_U;



	function null_logLR(Nw) {
		var winners, losers, others;
		for (i = 0; i < Wn; i++) {
			winners += Math.log(Nw - i);
		}
		for (i = 0; i < Ln; i++) {
			losers += Math.log(Nw - null_margin - i);
		}
		for (i = 0; i < Un; i++) {
			others += Math.log(popsize - (2 * Nw) + null_margin - i);
		}
		return winners + losers + others;
	}


	var upper_Nw_limit = Math.ceil((popsize - Un + null_margin) / 2) - 1;
	var lower_Nw_limit = Math.max(Wn, Ln+null_margin);

	if (lower_Nw_limit > upper_Nw_limit) {
		var temp = lower_Nw_limit;
		lower_Nw_limit = upper_Nw_limit;
		upper_Nw_limit = temp;
	}

	var intervals = [lower_Nw_limit, upper_Nw_limit];

	function f_upper(x) {
		var winners, losers;
		winners = 0;
		losers = 0;
		for (i = 0; i < Wn; i++) {
			winners += Math.log(x - i);
		}
		for (i = 0; i < Ln; i++) {
			losers += Math.log(x - null_margin - i);
		}
		return winners + losers;
	}
	function f_lower(x) {
		var others = 0;
		for (i = 0; i < Un; i++) {
			others += Math.log(popsize - (2 * x) + null_margin - i);
		}
		return others;
	}

	f_eval = [null_logLR(intervals[0]), null_logLR(intervals[1])];
	f_upper_eval = [f_upper(intervals[0]), f_upper(intervals[1])];
	f_lower_eval = [f_lower(intervals[0]), f_lower(intervals[1])];
	f_ub = [n * Math.log(popsize)];

	function argMax(array) {
	/**
	Retrieve the array key corresponding to the largest element in the array.
	**/
		return array.map((x, i) => [x, i]).reduce((r, a) => (a[0] > r[0] ? a : r))[1];
	}

	while (max(f_ub) > max(f_eval)) {
		var max_index = argMax(f_ub);
		var midpoint_Nw = int((intervals[max_index] + intervals[max_index + 1] / 2))
		if (intervals.indexOf(midpoint_Nw) == -1) {
			break;
		}
		intervals.insert(max_index + 1, midpoint_Nw);
		f_eval.insert(max_index + 1, f_upper(midpoint_Nw));
		f_lower.insert(max_index + 1, f_lower(midpoint_Nw));
		f_ub = f_upper_eval.slice(1, f_upper_eval.length) + f_lower_eval[0: f_lower_eval.length - 1]
	}

	fmax = argMax(f_eval);
	nuisance_param = intervals[argMax(f_eval)];
	number_invalid = popsize - (2 * nuisance_param) + null_margin;

	if (nuisance_param < 0 || nuisance_param > popsize) {
		return {
			'decision' : 'Number invalid is incompatible with the null.',
			'lower_threshold' : lower,
			'upper_threshold' : upper,
			'LR' : Infinity,
			'pvalue' : 0,
			'sample_proportion' : (Wn / n, Ln / n, Un / n),
			'Nu_used' : number_invalid,
			'Nw_used' : nuisance_param
		};
	} 

	if (nuisance_param < Wn || (nuisance_param - null_margin) < Ln ||
		number_invalid < Un) {
		return {
			'decision' : 'Null is impossible given the sample.',
			'lower_threshold' : lower,
			'upper_threshold' : upper,
			'LR' : Infinity,
			'pvalue' : 0,
			'sample_proportion' : (Wn / n, Ln / n, Un / n),
			'Nu_used' : number_invalid,
			'Nw_used' : nuisance_param
		};
	}

	var logLR = alt_logLR - fmax;
	var LR = Math.exp(logLR);

	if (LR <= lower) {
		decision = 0;
	} else if (LR >= upper) {
		decision = 1;
	}

	return {
		'decision' : decision,
		'lower_threshold' : lower,
		'upper_threshold' : upper,
		'LR' : LR,
		'pvalue' : min(1, 1/LR),
		'sample_proportion' : (Wn / n, Ln / n, Un / n),
		'Nu_used' : number_invalid,
		'Nw_used' : nuisance_param
	};
}
