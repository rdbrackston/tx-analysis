# tx-analysis
Code for analysis of transcriptional burst kinetics from single cell data. The included code infers parameters for a zero-inflated negative binomial model given single cell mRNA-FISH data. A Bayesian inference is implemented using a Metropolis-Hastings MCMC algorithm. Maximum a posteriori parameter estimates and 95% credible interbvals are then evaluated from the MCMC chains.

`Run_ZINB.jl` runs the MCMC on the files within the folder which must be specified within the script.
`Stats.jl` extracts the MAP and CI from the MCMC chains.
`Utils.jl` contains various functions used within the other scripts.
