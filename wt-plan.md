## Tasks for Election Audit Web Tools
### 24 September 2018

### Utilities

1. Implement the randint function from cryptorandom in javascript.

1. Implement the fykd sampling method and sample_by_index from cryptorandom in javascript.

1. Make a Javascript library including the SHA256 generator, randint, sample_by_index, and fykd

### Ballot polling

1. Make a version of ballotPollTool that uses sampling without replacement.
    + Rely on the new Javascript library for PRNGs and sampling
    + Implement risk calculations for ballot-polling without replacement (BBPwoR)
    + Develop formulae for estimating required sample size for BBPwoR
    
1. Make a new ballotPollTool page that implements both sampling with replacement and sampling without
replacement

### Ballot-level comparisons

1. Make a version of auditTool that uses sampling without replacement.
    + Rely on the new Javascript library for PRNGs and sampling
    + Implement risk calculations for ballot-level comparison audits without replacement
    + Develop formulae for estimating required sample size 
    
1. Make a new auditTool page that implements both sampling with replacement and sampling without
replacement

### SUITE

1. Implement SUITE for two strata in Javascript
    + Implement risk calculation for BBPwoR against a fixed threshold
    + Implement risk calculation for ballot-level comparison audits without replacement
    + Implement Fisher's combining function and corresponding p-value calculation 
(SticiGui has the chi-square distribution already--can borrow from that)
        - Implement piecewise constant upper and lower bounds on the p-value
        _ Implement search for maximum p-value across allocations of misfit   
    + Implement rules for initial sample sizes
    + Implement escalation rule
    + Extend the auditTools page to accommodate two ballot manifests
    + Rely on the new JavaScript PRNG library for sampling from the two strata independently
    + Design the UI
    