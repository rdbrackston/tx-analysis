# Script to run through the files and fit the zero-inflated negative binomial distrubution

using Plots, Distributions, DelimitedFiles, Base.Printf
import SignalUtils, TxModels, StatsUtils, PltUtils
gr()

# Set some parameters
Lchain = 400000
burn = 100000
thin = 100
folder = "/Users/rowanbrackston/Box Sync/06-Projects/07-PspNoise/01-Psp_Data/"

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
# Make a stats variable to be writtten to a file
Stats = Array{Any,2}(undef,length(Files)+1,16)
Stats[1,:] = ["Name" "Mean" "Variance" "No. cells" "p" "w" "r" "K" "r_min" "r_max" "p_min" "p_max" "w_min" "w_max" "K_min" "K_max"]

for (ii,file) in enumerate(Files)
    if occursin(".csv",file)

        name = replace(file, ".csv"=>"")

    	# Load and plot the distribution data
        data = TxModels.load_data(file,folder, 151)
        x,y = SignalUtils.genpdf(Integer.(round.(data)))

        # Fit the zero-inflated negatiove binomial model
        lFunc = p->StatsUtils.log_likelihood_zi_negbinom(data, p)
        # d = StatsUtils.mle_negbinom(filter(x->x>0, data))
        # wGuess = max((y[1]-pdf(d,0))/(1-pdf(d,0)),0.0)
        # guess = [params(d)[1],params(d)[2],wGuess]
        # println(guess)
        guess = [1.0,0.5,0.5]

        # Run the MCMC, then remove spurious zeros from the chain
        priors = [Truncated(Normal(0,20),0,Inf), Uniform(0.0,1.0),Uniform(0.0,1.0)]
		chain = TxModels.mcmc_metropolis(guess, lFunc, Lchain; prior=priors, propStd=0.03,scaleProp=false, burn=burn,step=thin);
		chain_red = rmv_zeros(chain)

        # Chain for K/ν
        tmp = (1.0.-chain_red[:,2])./chain_red[:,2]
        chain_red = [chain_red tmp]
        
		# Extract parameters and plot
        r = TxModels.find_MAP(chain_red,idx=1)
        p = TxModels.find_MAP(chain_red,idx=2)
        w = TxModels.find_MAP(chain_red,idx=3)
        K = TxModels.find_MAP(chain_red,idx=4)
        int_r = TxModels.credibleintervals(chain_red,idx=1)
        int_p = TxModels.credibleintervals(chain_red,idx=2)
        int_w = TxModels.credibleintervals(chain_red,idx=3)
        int_K = TxModels.credibleintervals(chain_red,idx=4)

        # Plot and save
        plt = TxModels.plot_chain(chain_red, [r,p,w,K], [int_r,int_p,int_w,int_K])
        Plots.pdf("MCMC/"*name)
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
        Stats[ii+1,:] = [name mean(data) var(data) length(data) p w r K r-int_r[1] int_r[2]-r p-int_p[1] int_p[2]-p w-int_w[1] int_w[2]-w K-int_K[1] int_K[2]-K]
    end
end

DelimitedFiles.writedlm("Statistics.txt", Stats)