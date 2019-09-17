# Script to run through the files and fit the zero-inflated negative binomial distrubution

using Plots, Distributions, DelimitedFiles, Base.Printf
import Utils
gr()

# Set some MCMC parameters
Lchain = 400000
burn = 100000
thin = 100
folder = "/Users/rowanbrackston/Box Sync/06-Projects/07-PspNoise/01-Psp_Data/"
# restart = "OxyR"

Files = readdir(folder)
filter!(f->occursin(".csv",f), Files)

for (ii,file) in enumerate(Files)
    if occursin(".csv",file)

        name = replace(file, ".csv"=>"")

    	# Load and plot the distribution data
        data = Utils.load_data(file,folder, 151)
        x,y = Utils.genpdf(Integer.(round.(data)))

        # Specify likelihood function and initial guess
        lFunc = p->Utils.log_likelihood_zi_negbinom(data, p)
        guess = [1.0,0.5,0.5]

        # Run the MCMC
        priors = [Truncated(Normal(0,20),0,Inf), Uniform(0.0,1.0),Uniform(0.0,1.0)]
		chain = Utils.mcmc_metropolis(guess, lFunc, Lchain; prior=priors, propStd=0.03,scaleProp=false, burn=burn,step=thin);
		chain_red = Utils.rmv_zeros(chain)

        # Chain for K/ν
        tmp = (1.0.-chain_red[:,2])./chain_red[:,2]
        chain_red = [chain_red tmp]
        
		# Extract parameters and plot
        r = Utils.find_MAP(chain_red,idx=1)
        p = Utils.find_MAP(chain_red,idx=2)
        w = Utils.find_MAP(chain_red,idx=3)
        K = Utils.find_MAP(chain_red,idx=4)

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
        Plots.pdf("Distributions/"*name)
    end
end
