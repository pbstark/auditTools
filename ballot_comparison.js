function ballot_comparison_pvalue(n, gamma, o1, u1, o2, u2, reported_margin, N, null_lambda) {
    /*
    Compute the p-value for a ballot comparison audit using Kaplan-Markov
    
    Parameters
    ----------
    n : int
        sample size
    gamma : float
        value > 1 to inflate the error bound, to avoid requiring full hand count for a single 2-vote overstatement
    o1 : int
        number of ballots that overstate any 
        margin by one vote but no margin by two votes
    u1 : int
        number of ballots that understate any margin by 
        exactly one vote, and every margin by at least one vote
    o2 : int
        number of ballots that overstate any margin by two votes
    u2 : int
        number of ballots that understate every margin by two votes
    reported_margin : float
        the smallest reported margin *in votes* between a winning
        and losing candidate for the contest as a whole, including any other strata
    N : int
        number of votes cast in the stratum
    null_lambda : float
        fraction of the overall margin (in votes) to test for in the stratum. If the overall margin is reported_margin,
        test that the overstatement in this stratum does not exceed null_lambda*reported_margin
    Returns
    -------
    pvalue
    */
    null_lambda = (null_lambda !== undefined) ? null_lambda : 1; // set default value of null_lambda to 1...
    var U_s = 2*N/reported_margin;
    log_pvalue = n*Math.log(1 - null_lambda/(gamma*U_s)) - 
                    o1*Math.log(1 - 1/(2*gamma)) - 
                    o2*Math.log(1 - 1/gamma) - 
                    u1*Math.log(1 + 1/(2*gamma)) - 
                    u2*Math.log(1 + 1/gamma);
    pvalue = Math.exp(log_pvalue);
    return Math.min(pvalue, 1);
}


function findNmin_ballot_comparison(alpha, gamma, o1, u1, o2, u2,
                                reported_margin, N, null_lambda) {

    /*
    Compute the smallest sample size for which a ballot comparison 
    audit, using Kaplan-Markov, with the given statistics could stop
    
    Parameters
    ----------
    alpha : float
        risk limit
    gamma : float
        value > 1 to inflate the error bound, to avoid requiring full hand count for a single 2-vote overstatement
    o1 : int
        number of ballots that overstate any 
        margin by one vote but no margin by two votes
    u1 : int
        number of ballots that understate any margin by 
        exactly one vote, and every margin by at least one vote
    o2 : int
        number of ballots that overstate any margin by two votes
    u2 : int
        number of ballots that understate every margin by two votes
    reported_margin : float
        the smallest reported margin *in votes* between a winning
        and losing candidate in the contest as a whole, including any other strata
    N : int
        number of votes cast in the stratum 
    null_lambda : float
        fraction of the overall margin (in votes) to test for in the stratum. If the overall margin is reported_margin,
        test that the overstatement in this stratum does not exceed null_lambda*reported_margin
        
    Returns
    -------
    n
    */
    null_lambda = (null_lambda !== undefined) ? null_lambda : 1; // set default value of null_lambda to 1...
    U_s = 2*N/reported_margin;
    val = -gamma*U_s/null_lambda * (Math.log(alpha) +
                o1*Math.log(1 - 1/(2*gamma)) + 
                o2*Math.log(1 - 1/gamma) + 
                u1*Math.log(1 + 1/(2*gamma)) + 
                u2*Math.log(1 + 1/gamma) );
    val2 = o1+o2+u1+u2;
    return Math.max(Math.floor(val)+1, val2);
}


function findNmin_ballot_comparison_rates(alpha, gamma, r1, s1, r2, s2,
                                reported_margin, N, null_lambda) {

    /**
    Compute the smallest sample size for which a ballot comparison 
    audit, using Kaplan-Markov, with the given statistics could stop
    
    Parameters
    ----------
    alpha : float
        risk limit
    gamma : float
        value > 1 to inflate the error bound, to avoid requiring full hand count for a single 2-vote overstatement
    r1 : int
        hypothesized rate of ballots that overstate any 
        margin by one vote but no margin by two votes
    s1 : int
        hypothesizedrate of ballots that understate any margin by 
        exactly one vote, and every margin by at least one vote
    r2 : int
        hypothesizedrate of ballots that overstate any margin by two votes
    s2 : int
        hypothesizedrate of ballots that understate every margin by two votes
    reported_margin : float
        the smallest reported margin *in votes* between a winning
        and losing candidate in the contest as a whole, including any other strata
    N : int
        number of votes cast in the stratum
    null_lambda : float
        fraction of the overall margin (in votes) to test for in the stratum. If the overall margin is reported_margin,
        test that the overstatement in this stratum does not exceed null_lambda*reported_margin
        
    Returns
    -------
    n
    **/
    null_lambda = (null_lambda !== undefined) ? null_lambda : 1; // set default value of null_lambda to 1...
    U_s = 2*N/reported_margin;
    denom = (Math.log(1 - null_lambda/(U_s*gamma)) -
                r1*Math.log(1 - 1/(2*gamma))- 
                r2*Math.log(1 - 1/gamma) - 
                s1*Math.log(1 + 1/(2*gamma)) - 
                s2*Math.log(1 + 1/gamma) );
    return denom < 0 ? Math.ceil(Math.log(alpha)/denom) : Number.NaN;
}