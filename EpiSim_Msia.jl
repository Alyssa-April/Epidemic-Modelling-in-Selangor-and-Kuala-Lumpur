module EpiSim_Msia

    import Pkg
    Pkg.activate("new_environment")
    Pkg.add("FFMPEG")
    Pkg.add("Graphs")
    using Random

    export epistep!, episim_twocpt_samerate_MCOCMCO, episim_nochangepoint,
    plotqnt_normal, plotqnt_log, homophily_epistep!,
    episim_homophily_twocpt_samerate_MCOCMCO, episim_homophily_nochangepoint

    using ProgressBars, CSV, Plots, Statistics, Graphs

    function epistep!(state, net, β, σ, γ, ξ) 
        # perform a single step of the epidemic propagation process (one day)
        # loop over the nodes, updating the infection states
        # state is a vector and net is a contact network with commensurate dimensions
        # beta, sigma, gamma and xi are the epidemic parameters
        # (transmission rate, latency period, recovery rate, reinfection rate)
        for v in Graphs.vertices(net) # looping through nodes 
            if state[v] == 1 # look at current state of the node, if it is susceptible, check its neighbors
                for n in all_neighbors(net, v) # check all neighbors of the node v
                    if state[n] == 3 && rand(Float64) .< β # if neighbour is infectious, and pr < transmission rate
                        state[v] = 2   # S becomes E when meeting an I node
                    elseif state[n] == 13 && rand(Float64) .< β
                        state[v] = 2   # S becomes E when meeting a newly I node
                    elseif state[n] == 14 && rand(Float64) .< β
                        state[v] = 2   # S becomes E when meeting a newly R node
                    end
                end
            elseif state[v] == 2 && rand(Float64) .< σ
                state[v] = 13 # exposed becomes infectious (state 13 is newly infectious)
            elseif state[v] == 3 && rand(Float64) .< γ
                state[v] = 14 # infectious becomes removed (state 14 is newly removed)
            elseif state[v] == 4 && rand(Float64) .< ξ
                state[v] = 1
            end
        end
        # update the newly infected and removed nodes
        # a newly removed node cannot enter state S within the same day, it can enter S after one day
        state[state.>10] = state[state.>10] .- 10 
    end

    function homophily_epistep!(state, net, ethnic_group, β_intra, β_inter, σ, γ, ξ)
        # perform epistep for homophily network, where
        # the transmission rate between groups is smaller than within groups
        for v in Graphs.vertices(net)
            if state[v] == 1
                for n in all_neighbors(net, v)
                    
                    # Determine the transmission rate based on group relationship
                    β = ethnic_group[v] == ethnic_group[n] ? β_intra : β_inter
                    
                    if state[n] == 3 && rand(Float64) .< β
                        state[v] = 2   # S becomes E when meeting an I node
                    elseif state[n] == 13 && rand(Float64) .< β
                        state[v] = 2   # S becomes E when meeting a newly I node
                    elseif state[n] == 14 && rand(Float64) .< β
                        state[v] = 2   # S becomes E when meeting a newly R node
                    end
                end
            elseif state[v] == 2 && rand(Float64) .< σ
                state[v] = 13 # exposed becomes infectious (state 13 is newly infectious)
            elseif state[v] == 3 && rand(Float64) .< γ
                state[v] = 14 # infectious becomes removed (state 14 is newly removed)
            elseif state[v] == 4 && rand(Float64) .< ξ
                state[v] = 1
            end
        end
        # update the newly infected and removed nodes
        state[state.>10] = state[state.>10] .- 10
    end
    
    # use for BA only
    function episim_nochangepoint(net, epiparam, ndays::Int64 = 127, 
        nsims::Int64 = 1000, seed::Int64 = 1234)
        # no changepoint model (scale-free network with no MCO implementation)

        # epidemic parameters
        pop = epiparam["pop"]         # population
        β = epiparam["β_fixed"]       # transmission rate throughout
        σ = epiparam["σ"]             # latency period
        γ = epiparam["γ"]             # recovery rate throughout
        ξ = epiparam["ξ"]             # reinfection rate
        nseeds = epiparam["nseeds"]   # number of infected nodes to start with

        # initialise and get set to track the various tallies
        st = Array{UInt64,2}(undef,ndays,nsims)
        ex = Array{UInt64,2}(undef,ndays,nsims)
        fe = Array{UInt64,2}(undef,ndays,nsims)
        rm = Array{UInt64,2}(undef,ndays,nsims)

        iter = ProgressBar(1:nsims) 

        # main iteratarion loop nsims simulations
        Threads.@threads for j in iter
            # Set random seed for consistency
            Random.seed!(seed) 

            # reinitialise the state vector
            state = Array{Int8,2}(undef,1,pop) 
            state[1:pop].=1; # start with the whole pop in S
            state[rand((1:pop),(1,nseeds))].=3; # nseeds are infected

            i = 1
            # loop for a single epidemic time course 
            while i <= (ndays)

                # updating the SEIR states
                epistep!(state, net, β, σ, γ, ξ)

                # count the respective totals
                st[i,j] = count(x->x==1, state)
                ex[i,j] = count(x->x==2, state)
                fe[i,j] = count(x->x==3, state)
                rm[i,j] = count(x->x==4, state)

                i += 1  # done one more day (one time step)

            end
        end

        return st, ex, fe, rm

    end 

    # use for BA only
    function episim_homophily_nochangepoint(net, epiparam, ethnic_group, 
        ndays::Int64 = 127, nsims::Int64 = 1000, seed::Int64 = 1234)
        # no changepoint model (scale-free network with no MCO implementation)

        # epidemic parameters
        pop = epiparam["pop"]         # population
        β_intra = epiparam["β_intra"] # transmission rate for within groups
        β_inter = epiparam["β_inter"] # transmission rate for between groups
        σ = epiparam["σ"]             # latency period
        γ = epiparam["γ"]             # recovery rate throughout
        ξ = epiparam["ξ"]             # reinfection rate
        nseeds = epiparam["nseeds"]   # number of infected nodes to start with

        # initialise and get set to track the various tallies
        st = Array{UInt64,2}(undef,ndays,nsims)
        ex = Array{UInt64,2}(undef,ndays,nsims)
        fe = Array{UInt64,2}(undef,ndays,nsims)
        rm = Array{UInt64,2}(undef,ndays,nsims)

        iter = ProgressBar(1:nsims) 

        # main iteration loop nsims simulations
        Threads.@threads for j in iter 
            # Set random seed for consistency 
            Random.seed!(seed) 

            # reinitialise the state vector
            state = Array{Int8,2}(undef,1,pop)
            state[1:pop].=1; # start with the whole pop in S
            state[rand((1:pop),(1,nseeds))].=3; # nseeds are infected

            i = 1
            # loop for a single epidemic time course 
            while i <= (ndays)

                # updating the SEIR states
                homophily_epistep!(state, net, ethnic_group, β_intra, β_inter, σ, γ, ξ)

                # count the respective totals
                st[i,j] = count(x->x==1, state)
                ex[i,j] = count(x->x==2, state)
                fe[i,j] = count(x->x==3, state)
                rm[i,j] = count(x->x==4, state)

                i += 1  # done one more day (one time step)

            end
        end

        return st, ex, fe, rm

    end 

    # use for BA, switch to MCO, then switch to CMCO until end (B4 and after MCO same rates)
    function episim_twocpt_samerate_MCOCMCO(net1, net2, net3, epiparam, MCOstart::Int64 = 44, CMCOstart::Int64 = 91,
        postdays::Int64 = 36, nsims::Int64 = 1000, seed::Int64 = 1234)
        # epidemic parameters
        pop = epiparam["pop"]         # population 
        β = epiparam["β_fixed"]       # fixed transmission rate throughout (no homophily, only one rate)
        σ = epiparam["σ"]             # fixed latency period
        γ = epiparam["γ"]             # fixed recovery rate 
        ξ = epiparam["ξ"]             # fixed reinfection rate
        nseeds = epiparam["nseeds"]   # number of infected nodes to start with (2 for this paper)

        #initialise and get set to track the various tallies
        ndays = CMCOstart + postdays # 127 days overall
        st = Array{UInt64,2}(undef,ndays,nsims)
        ex = Array{UInt64,2}(undef,ndays,nsims)
        fe = Array{UInt64,2}(undef,ndays,nsims)
        rm = Array{UInt64,2}(undef,ndays,nsims)

        iter = ProgressBar(1:nsims) # progress bar

        # main iteratarion loop nsims simulations
        Threads.@threads for j in iter
            # Set random seed for consistency
            Random.seed!(seed) 

            # start with first network before MCO
            net = net1

            # reinitialise the state vector
            state = Array{Int8,2}(undef,1,pop) 
            state[1:pop].=1;
            state[rand((1:pop),(1,nseeds))].=3; # number of infected cases to start with

            i = 1
            # loop for a single epidemic time course
            while i <= ndays

                # switch infection model at the start of the MCO 
                # once achieved, the transmission model switches to net2
                if i == MCOstart   
                    net = net2 # change structure to MCO network                  
                end

                # switch infection model for the start of CMCO 
                # once achieved, the transmission model switches to net3
                if i == CMCOstart
                    net = net3
                end

                # updating the infection states
                epistep!(state, net, β, σ, γ, ξ)

                # count the respective totals
                st[i,j] = count(x->x==1, state)
                ex[i,j] = count(x->x==2, state)
                fe[i,j] = count(x->x==3, state)
                rm[i,j] = count(x->x==4, state)

                i += 1  # done one more day

            end
        end

        return st, ex, fe, rm

    end # end of episim_twocpt_samerate_MCOCMCO

    # use for BA, MCO, then CMCO until end (B4 and after MCO same rates)
    function episim_homophily_twocpt_samerate_MCOCMCO(net1, net2, net3, epiparam, 
        ethnic_group, MCOstart::Int64 = 44, CMCOstart::Int64 = 91,
        postdays::Int64 = 36, nsims::Int64 = 1000, seed::Int64 = 1234)
        # epidemic parameters
        pop = epiparam["pop"]         # population 
        β_intra = epiparam["β_intra"] # fixed transmission rate for within group
        β_inter = epiparam["β_inter"] # fixed transmission rate for between group  
        σ = epiparam["σ"]             # fixed latency period
        γ = epiparam["γ"]             # fixed recovery rate 
        ξ = epiparam["ξ"]             # reinfection rate
        nseeds = epiparam["nseeds"]   # number of infected nodes to start with

        #initialise and get set to track the various tallies
        ndays = CMCOstart + postdays
        st = Array{UInt64,2}(undef,ndays,nsims)
        ex = Array{UInt64,2}(undef,ndays,nsims)
        fe = Array{UInt64,2}(undef,ndays,nsims)
        rm = Array{UInt64,2}(undef,ndays,nsims)

        iter = ProgressBar(1:nsims) # progress bar

        # main iteratarion loop nsims simulations
        Threads.@threads for j in iter
            # Set random seed for consistency
            Random.seed!(seed) 

            # start with first network before MCO
            net = net1

            # reinitialise the state vector
            state = Array{Int8,2}(undef,1,pop) 
            state[1:pop].=1;
            state[rand((1:pop),(1,nseeds))].=3; # number of infected cases to start with

            i = 1
            # loop for a single epidemic time course (one day)
            while i <= ndays

                # switch infection model at the start of the MCO 
                # once achieved, the transmission model switches to net2
                if i == MCOstart   
                    net = net2 # change structure to MCO network                  
                end

                # switch infection model for the start of CMCO 
                # once achieved, the transmission model switches to net3
                if i == CMCOstart
                    net = net3
                end

                # updating the infection states
                homophily_epistep!(state, net, ethnic_group, β_intra, β_inter, σ, γ, ξ)

                # count the respective totals
                st[i,j] = count(x->x==1, state)
                ex[i,j] = count(x->x==2, state)
                fe[i,j] = count(x->x==3, state)
                rm[i,j] = count(x->x==4, state)

                i += 1  # done one more day 

            end
        end

        return st, ex, fe, rm

    end # end of episim_homophily_twocpt_samerate_MCOCMCO


    function plotqnt_normal(x, y, col, labl, qnt::Float64 = 0.45) # we use the median 
        nt,ny = size(y)
        low = Array{Float64,1}(undef,nt)
        mid = Array{Float64,1}(undef,nt)
        hig = Array{Float64,1}(undef,nt)  
        for i in 1:nt
            low[i],mid[i],hig[i] = quantile(y[i,:],[0.5-qnt, 0.5, 0.5+qnt]) # to add ribbon
        end
        plot!(x, mid, grid = false, ribbon = (mid-low,hig-mid), fillalpha = 0.25,
        lw = 1.5, seriescolor = col, label = labl, #xticks=(x, string.(x)), #xrotation = 90, 
        xtickfontsize = 8)
        #return low, mid, hig   
    end

    function plotqnt_log(x, y, col, labl, qnt::Float64 = 0.45) 
        nt, ny = size(y)
        low = Array{Float64,1}(undef,nt) 
        mid = Array{Float64,1}(undef,nt)
        hig = Array{Float64,1}(undef,nt)
        for i in 1:nt
            low[i],mid[i],hig[i] = quantile(y[i,:],[0.5-qnt, 0.5, 0.5+qnt])
        end
        plot!(x, mid, grid = false, ribbon = (mid-low,hig-mid), fillalpha = 0.25,
        lw = 1.5, seriescolor = col, label = labl, #xticks=(x, string.(x)), xrotation = 90, 
        xtickfontsize = 8, yaxis = :log) # log scale for y axis to enable better visualisation
    end
end 
