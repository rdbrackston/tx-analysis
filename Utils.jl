module Utils

using Distributions, Discretizers, Optim, CSV, Plots
import Printf
import KernelDensity; const KDE = KernelDensity


"""
Function to perform the MCMC metropolis algorithm.
"""
function mcmc_metropolis(x0::AbstractArray, logPfunc::Function, Lchain::Integer;
                         propStd::Union{AbstractFloat,AbstractArray}=0.1,
                         step::Integer=500, burn::Integer=500, printFreq::Integer=10000,
                         prior=:none, verbose::Bool=true, scaleProp::Bool=true)

	if length(size(x0)) > 1 # restart from old chain
        if verbose; println("Restarting from old chain"); end
		xOld = x0[end,:]
		n = length(xOld)
	else
	    xOld = x0
	    n = length(x0)
	end
    chain = zeros(Float64, Lchain,n)
    chain[1,:] = xOld
    acc = 0

    if prior == :none
    	logpOld = logPfunc(xOld)
    else
        logpOld = logPfunc(xOld)
        for  (ip,prr) in enumerate(prior)
        	logpOld += log(Distributions.pdf(prr,xOld[ip]))
        end
    end

    for ii=2:Lchain
        if scaleProp
            proposal = MvNormal(propStd.*xOld)
        else
            if isa(propStd, Array)
                proposal = MvNormal(propStd)
            else
                proposal = MvNormal(propStd.*ones(size(xOld)))
            end
        end
        xNew = xOld + Distributions.rand(proposal)

        # Reject if any x are negative
        if any(x->x<0, xNew); continue; end

        # Evaluate new log-likelihood
        if prior == :none
	    	logpNew = logPfunc(xNew)
	    else
	        logpNew = 0.0
	        for  (ip,prr) in enumerate(prior)
	        	logpNew += log(Distributions.pdf(prr,xNew[ip]))
	        end
	        if !isinf(logpNew); logpNew += logPfunc(xNew); end
	    end
        a = exp(logpNew - logpOld)

        if Base.rand() < a
            xOld = xNew
            logpOld = logpNew
            acc += 1
        end
        chain[ii,:] = xOld

        if ii % printFreq == 0
            if verbose
                Printf.@printf("Completed iteration %i out of %i. \n", ii,Lchain)
            end
        end
    end

    if verbose
        Printf.@printf("Acceptance ratio of %.2f (%i out of %i).\n", acc/Lchain,acc,Lchain)
    end

	if length(size(x0)) > 1
		chainRed = [x0; chain[step:step:end,:]] # Append new chain to old
	else
		chainRed = chain[burn:step:end,:]
	end
	return chainRed

end


"""
Function to find the MAP from the posterior distribution.
Finds independently for each parameter.
"""
function find_MAP(chain; idx=1)

    if length(size(chain))>1
        data = chain[:,idx]
    else
        data = copy(chain)
    end

    # Perform the KDE
    kdeObj = kde_wrpr(data)

    # Find MAP via crude optimisation
    func = x->KDE.pdf(kdeObj,x[1])
    x = collect(range(min(data...),stop=max(data...),length=10000))
    y = map(func,x)
    MAP = x[y.==maximum(y)][1]

    return MAP

end


"""
Function to find the MAP from the posterior distribution.
Finds independently for each parameter.
"""
function credibleintervals(chain; idx=1, spread=0.682)

    if length(size(chain))>1
        data = chain[:,idx]
    else
        data = copy(chain)
    end

    # Perform the KDE
    kdeObj = kde_wrpr(data)
    func = z->KDE.pdf(kdeObj,z)
    x = collect(range(min(data...),stop=max(data...),length=10000))
    dx = x[2]-x[1]
    y = map(func,x)

    # Move out from MAP in both directions
    # Stop once intervals contain "spread" of probability
    MAP = find_MAP(data)
    xMin = last(x[x.<MAP])
    xMax = first(x[x.>MAP])
    tot = 0.5*dx*(KDE.pdf(kdeObj,xMin)+KDE.pdf(kdeObj,xMax))
    flag = false
    while tot<spread
        xMinN = xMin-dx
        xMaxN = xMax+dx
        if xMinN<0.0
            flag = true
            break
        end
        if KDE.pdf(kdeObj,xMinN)>KDE.pdf(kdeObj,xMaxN)
            tot += 0.5*dx*(KDE.pdf(kdeObj,xMin)+KDE.pdf(kdeObj,xMinN))
            xMin = xMinN
        else
            tot += 0.5*dx*(KDE.pdf(kdeObj,xMax)+KDE.pdf(kdeObj,xMaxN))
            xMax = xMaxN
        end
    end

    if flag
        sort!(data)
        idx = Int(round(spread*length(data)))
        xMin = 0.0
        xMax = data[idx]
    end

    return xMin,xMax

end


"""
Wrapper around the KDE function to ensure that all calls to KDE use the same settings
"""
function kde_wrpr(data)
    KDE.kde(data, boundary=(0.0,2*max(data...)))
end


"""
Function to plot the full chains and posterior distributions.
"""
function plot_chain(chain, MAP=:none, Ints=:none)

	nVar = size(chain)[2]
	plots = Array{Any,1}(undef,nVar*2)

    # Set MAP to the mean if not provided
    if isequal(MAP,:none)
        MAP = mean(chain,dims=1)
        lbl = "Mean"
    else
        lbl = "MAP"
    end


	for ii=1:nVar
		plots[ii] = plot(chain[:,ii], legend=false, title=Printf.@sprintf("Parameter %i",ii));

		smpls = chain[:,ii]
		kdeObj = kde_wrpr(smpls)
		x = collect(range(Base.minimum(smpls),stop=Base.maximum(smpls),length=100))
		y = map(z->KDE.pdf(kdeObj,z),x)
		if ii ==1
			tmp = plot(x,y, label="Posterior", legend=:top);
			plot!(tmp, [MAP[ii],MAP[ii]],[0,Base.maximum(y)], label=lbl, line=(2,:black));
		else
			tmp = plot(x,y, legend=false);
			plot!(tmp, [MAP[ii],MAP[ii]],[0,Base.maximum(y)], line=(2,:black));
		end
        if !isequal(Ints,:none)
            plot!(tmp, [Ints[ii][1],Ints[ii][1]],[0,Base.maximum(y)], line=(2,:black,:dash), label="");
            plot!(tmp, [Ints[ii][2],Ints[ii][2]],[0,Base.maximum(y)], line=(2,:black,:dash), label="");
        end
		plots[ii+nVar] = tmp;

	end

	if nVar == 2
		plt = plot(plots[1],plots[2],plots[3],plots[4],
			       layout=(2,nVar), size=(800,800))
	elseif nVar == 3
		plt = plot(plots[1],plots[2],plots[3],plots[4],plots[5],plots[6],
			       layout=(2,nVar), size=(1200,800))
	elseif nVar == 4
		plt = plot(plots[1],plots[2],plots[3],plots[4],plots[5],plots[6],plots[7],plots[8],
			       layout=(2,nVar), size=(1600,800))
	elseif nVar == 5
		plt = plot(plots[1],plots[2],plots[3],plots[4],plots[5],plots[6],plots[7],plots[8],plots[9],plots[10],
			       layout=(2,nVar), size=(2000,800))
	elseif nVar == 6
		plt = plot(plots[1],plots[2],plots[3],plots[4],plots[5],plots[6],plots[7],plots[8],plots[9],plots[10],plots[11],plots[12],
			       layout=(2,nVar), size=(2400,800))
	else
		Printf.@printf("nVar of %i not currently catered for.", nVar)
		plt = plot()
	end

	return plt

end


"""
Function to load in the experimental data, filtering out missings/NaNs and
applying upper cut-off in terms of standard deviations from the mean.
"""
function load_data(File::String, Folder::String, cutOff::Number=Inf, format=:Int)

    if occursin(".csv",File)
        rnaData = CSV.read(Folder*File, datarow=1)[:,1]
    else
	    rnaData = CSV.read(Folder*File*".csv", datarow=1)[:,1]
    end

	filter!(x -> !isnan(x), rnaData)

    fltData = filter(x->x<cutOff, rnaData)

    if format==:Int
	    return Integer.(round.(fltData))
    else
        return fltData
    end

end


"""
1D integer version of genpdf, intended for use when the underlying data is discrete.
"""
function genpdf(data::Array{Int64,1}, nbin=:auto::Union{Int,Symbol})

    # Set bins to all the integers between lo and hi
    lo, hi = extrema(data)
    edges = collect(lo-1:hi) .+ 0.5

    disc = LinearDiscretizer(edges)
    counts = get_discretization_counts(disc, data)

    x = edges[1:end-1]+0.5*diff(edges)
    y = (counts./(diff(edges)*float(size(data)[1]))) # Normalise
    return (x,y)

end


"""
Function to remove the spurious zeros that can arise during MCMC.
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


"""
Return log-likelihood of getting data X from zero-inflated negative binomial with prms
"""
function log_likelihood_zi_negbinom(X, prms, verbose=false)

    if verbose
        println(prms)
    end

    m = MixtureModel([DiscreteUniform(0,0), NegativeBinomial(prms[1:2]...)], [prms[3],1-prms[3]])
    sum(log.([pdf(m,X[ii]) for ii=1:length(X)]))
end


"""
Obtain confidence intervals for the statistic specified in statFunc using bootstrapping
"""
function bootstrap_intervals(data, statFunc::Function, N=10000)

	statistic = statFunc(data)
	nSamp = length(data)

	statVals = Array{Float64}(undef, N)
	for ii=1:N
		# Resample from data with replacement
		idxs = Int.(round.(rand(nSamp).*nSamp .+0.5))
		samps = data[idxs]
		statVals[ii] = statFunc(samps)
	end

	ints = std(statVals)
	return statistic-ints,statistic,statistic+ints

end

end # module
