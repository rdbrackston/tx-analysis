# Script to run through the files and fit the zero-inflated negative binomial distrubution

using Plots, Distributions, DelimitedFiles, Base.Printf
import SignalUtils, TxModels, StatsUtils, PltUtils
gr()

# Set some parameters
Lchain = 1000000
burn = 300000
thin = 100
folder = "/Users/rowanbrackston/Box Sync/06-Projects/07-PspNoise/01-Suf_Data/"
restart = "OxyR"

"""
Function to remove the spurious zeros that arise.
"""
function rmv_zeros(chain)

    l = length(chain[:,1])
    cnt = 0
    for ii=1:l
        if iszero(chain[ii,3])
            cnt+=1
        end
    end

    chain_red = Array{Float64,2}(undef, l-cnt,3)
    idx = 1
    for ii=1:l
        if !iszero(chain[ii,3])
            chain_red[idx,:] = chain[ii,:]
            idx+=1
        end
    end
    return chain_red
    
end

Files = readdir(folder)
filter!(f->occursin(".csv",f), Files)

for (ii,file) in enumerate(Files)
    if occursin(".csv",file) && occursin(restart*".csv",file)

        name = replace(file, ".csv"=>"")

    	# Load and plot the distribution data
        data = TxModels.load_data(file,folder, 151)
        x,y = SignalUtils.genpdf(Integer.(round.(data)))

        # Fit the zero-inflated negatiove binomial model
        lFunc = p->StatsUtils.log_likelihood_zi_negbinom(data, p)
        guess = [1.0,0.5,0.5]

        # Run the MCMC, then remove spurious zeros from the chain
        priors = [Truncated(Normal(0,20),0,Inf), Uniform(0.0,1.0),Uniform(0.0,1.0)]
        if !isequal(restart,:none)
            @assert isequal(restart,name) "Loaded chain and current file must be the same."
            chain = readdlm("Chains/"*restart)[:,1:3]
            chain = TxModels.mcmc_metropolis(chain, lFunc, Lchain; prior=priors, propStd=0.03,scaleProp=false, burn=burn,step=thin);
        else
            chain = TxModels.mcmc_metropolis(guess, lFunc, Lchain; prior=priors, propStd=0.03,scaleProp=false, burn=burn,step=thin);
		end
        chain_red = rmv_zeros(chain)

        # Chain for K/ν
        tmp = (1.0.-chain_red[:,2])./chain_red[:,2]
        chain_red = [chain_red tmp]
        
		# Extract parameters and plot
        r = TxModels.find_MAP(chain_red,idx=1)
        p = TxModels.find_MAP(chain_red,idx=2)
        w = TxModels.find_MAP(chain_red,idx=3)
        K = TxModels.find_MAP(chain_red,idx=4)

        # Plot and save
        plt = plot(plot(title=name), plot(legend=false), size=(600,800),layout=(2,1))
        bar!(plt[1], x,y, line=0,ylabel="Probability, P(n)", label="Data")
        plot!(plt[2], x,cumsum(y), xlabel="Copy number, n",ylabel="Cumulative probability, C(n)", line=:steppost)
        m = MixtureModel([DiscreteUniform(0,0), NegativeBinomial(r,p)], [w,1-w])
        plot!(plt[1], x, x->pdf(m,x), label="Model")
        plot!(plt[2], 0:maximum(x),n->cdf(m,n))

        # Save data
        println(name)
        writedlm("Chains/"*name,chain_red)
        PltUtils.rend_pdf("Distributions/"*name)
    end
end
