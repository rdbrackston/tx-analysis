# tx-analysis
Code for analysis of transcriptional burst kinetics from single cell data. The included code infers parameters for a zero-inflated negative binomial model given single cell mRNA-FISH data. A Bayesian inference is implemented using a Metropolis-Hastings MCMC algorithm. Maximum *a posteriori* (MAP) parameter estimates and 95% credible intervals (CI) can then be evaluated from the MCMC chains.

`Run_ZINB.jl` runs the MCMC inference on the files within a specified folder, generating a chain of inferred parameter samples.

`Stats.jl` extracts the MAP and CI from the MCMC chains.

`Utils.jl` contains various functions used within the other scripts.

## Requirements and Installation

The code is written in the programming language Julia which can be downloaded from [https://julialang.org/downloads/](https://julialang.org/downloads/) and installed on Windows, OSX or Linux. The code in this repository has been tested on Julia 1.3.0 but should run on all versions of Julia >1.0.

With Julia installed, one now needs to download the code and install the required Julia packages. First clone this repository, navigate into the project folder and launch Julia.
```bash
git clone https://github.com/rdbrackston/tx-analysis.git
cd tx-analysis
julia
```
From within Julia, enter the Pkg REPL-mode by pressing the `]` key. Then activate and instantiate the project to download and install the required packages.
```bash
(v1.3) pkg> activate .
(tx-analysis) pkg> instantiate
```

## Running the Demonstration

To run the code on the included demo data, first install Julia and activate the project as directed above. With an instance of julia initiated within the repsoitory folder, first run the MCMC inference (takes around 1 minute).
```bash
julia> include("Run_ZINB.jl")
```
This will run the inference algorithm for the parameters of a zero-inflated negative binomial distribution. The algorithm will operate on the file "PspF_High.csv". The resulting MCMC chain will be saved to a file within a folder "Chains", while a comparison of the distributions will be saved as a .pdf in a folder "Distributions".

To subsequently run analysis on the saved MCMC chain (may take several minutes).
```bash
julia> include("Stats.jl")
```
This will extract the maximum *a posteriori* estimates of the parameter values and the 95% credible intervals. These results are plotted and placed in a folder "MCMC", and also written to a file "Statistics.txt".

## Instructions for General Use

The process described above can be run on any single cell data conforming to certain restrictions. Each data file should be in .csv format with each entry corresponding to the copy number of mRNA within a single cell. Many such files can be placed within a single folder and the analysis will be run on each in turn. This folder must be specified within "Run_ZINB.jl".

## Related Publication

The code here accompanies the following paper:

C. Engl, G. Jovanovic, R. D. Brackston, I. Kotta-Loizou & M. Buck. The route to transcription initiation determines the mode of transcriptional bursting in bacteria. 2020, *To appear in Nat. Comms.*
