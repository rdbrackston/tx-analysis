# Script to run through the files and fit the zero-inflated negative binomial distrubution

using Plots, Distributions, DelimitedFiles, Base.Printf
import SignalUtils, TxModels, StatsUtils, PltUtils
gr()

# Set some parameters
spread = 0.95
folder = "Chains/"

Files = readdir(folder)
Stats = Array{Any,2}(undef,length(Files)+1,13)
Stats[1,:] = ["Name" "p" "w" "r" "K" "r_min" "r_max" "p_min" "p_max" "w_min" "w_max" "K_min" "K_max"]

for (ii,file) in enumerate(Files)

    # Load and plot the distribution data
    chain = readdlm(folder*file)
    
	# Extract parameters and plot
    r = TxModels.find_MAP(chain,idx=1)
    p = TxModels.find_MAP(chain,idx=2)
    w = TxModels.find_MAP(chain,idx=3)
    K = TxModels.find_MAP(chain,idx=4)
    int_r = TxModels.credibleintervals(chain,idx=1, spread=spread)
    int_p = TxModels.credibleintervals(chain,idx=2, spread=spread)
    int_w = TxModels.credibleintervals(chain,idx=3, spread=spread)
    int_K = TxModels.credibleintervals(chain,idx=4, spread=spread)

    # Plot and save
    plt = TxModels.plot_chain(chain, [r,p,w,K], [int_r,int_p,int_w,int_K])
    Plots.pdf("MCMC/"*file)

    # Save data
    println(file)
    Stats[ii+1,:] = [file p w r K r-int_r[1] int_r[2]-r p-int_p[1] int_p[2]-p w-int_w[1] int_w[2]-w K-int_K[1] int_K[2]-K]

end

DelimitedFiles.writedlm("Statistics.txt", Stats)