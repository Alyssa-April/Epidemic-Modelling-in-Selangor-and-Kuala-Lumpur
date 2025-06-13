module ContactModels_Msia

    export fullmixing, limitpop, limitpop_sbm, limitmix, limitmix_for_sbm, isolation, socialdist, rewire!
    export generate_sbm, remap!, rewire_sbm!, socialdist_sbm
    using ProgressBars, Graphs, Statistics
    using SparseArrays, Random, Combinatorics

    
    # function to generate SBM so that we can start with an SBM model instead of a BA model
    # this function randomly assigns edges based on the number of expected edges within and between groups
    function generate_sbm(number_per_gp::Vector{Int}, pr_mat::Matrix{Float64}, seed = 1234)
        Random.seed!(seed)

        # total number of nodes in the different ethnicity groups
        total_nodes = sum(number_per_gp)

        # set the index for each node in the network based on their ethnicity (1111...,22222...)
        node_ethnicity = Vector{Int}(undef, total_nodes)
        start_index = 1
        for (i, count) in enumerate(number_per_gp)
            node_ethnicity[start_index:start_index + count - 1] .= i
            start_index += count
        end 

        # obtain the start and end indices for each group
        # this ensures that the edges are materialised between the correct pair of nodes
        ethnic_start = cumsum([1; number_per_gp[1:end-1]]) # 1, 1+number in group 1
        ethnic_end = cumsum(number_per_gp) # number in group 1, number in group 1 + number in group 2
        # for example if number_per_gp = [10, 20], the first group starts at 1 and ends at 10, the second group starts at 11 and ends at 30

        # create a simple graph based on the total number of nodes
        # this is a placeholder and will be replaced with the actual adjacency matrix later
        # this avoids storing sparse matrices in memory and allows for more efficient memory usage
        net = SimpleGraph(total_nodes)
        
        total_number_of_expected_edges = 0  

        # loop through each group and assign edges based on the expected number of edges
        for i in 1:length(number_per_gp)
            for j in i:length(number_per_gp)
                # get the start and end indices for the two groups to ensure correct group assignment
                ni = number_per_gp[i]
                nj = number_per_gp[j]
                pr_ij = pr_mat[i, j]

                # expected number of edges for within and between groups
                expected_edges = i == j ? round(Int, pr_ij * ni * (ni - 1) / 2) : round(Int, pr_ij * ni * nj)

                total_number_of_expected_edges += expected_edges

                # loop through the nodes in the two groups and assign edges based on the expected number of edges
                for counter in 1:expected_edges
                    # randomly select two nodes from the two groups 
                    u = rand(ethnic_start[i]:ethnic_end[i])
                    v = rand(ethnic_start[j]:ethnic_end[j])
                    # ensure that the selected nodes are not the same to avoid self-loops
                    if u != v
                        add_edge!(net, u, v)
                    end
                    print(counter)
                end
            end 
        end
        println("The expected number of edges is $total_number_of_expected_edges")
        println("The actual number of edges is $(ne(net))") 
        return net 
    end


    function fullmixing(pop)
        # model a population of pop fully mixed nodes
        # a BA scale-free network model, mean degree 4
        if pop < 5
            println("Creating a complete graph for a population of $pop")
            net = complete_graph(pop)
            return net
        end
        println("Creating a Barabasi-Albert graph for a population of $pop")
        net = barabasi_albert(pop, 3, 2)
        # need to add six extra edges
        verts = vertices(net)
        for i in 1:6
            edg = Edge(rand(verts),rand(verts))
            while edg in edges(net)
                edg = Edge(rand(verts),rand(verts))
            end
            add_edge!(net,edg)
            println("Added edge $i/6")
        end
        println("Completed adding six extra edges to the network")
        return net
    end

# "Redistribution of pruned edges in scale-free BA networks continues transmission paths, 
#  while in SBM networks representing multiracial contact structure, such redistribution 
#  disrupts community boundaries and artificially inflates between-group mixing."


    function limitpop(net, gth)
        # prune contact patterns to limit gatherings on network to no more than gth people. 
        # edges are redistributed to preserve mean degree
        nlost = 0
        nedges = ne(net)
        println("Starting edge pruning to limit gatherings to $gth people")

        # prune edges from hubs
        @Threads.threads for vert in vertices(net)
            neigh = neighbors(net, vert)
            while length(neigh) > gth
                rem_edge!(net, vert, rand(neigh))
                nlost += 1
                println("Pruned an edge from vertex $vert (nlost: $nlost)")
                neigh = neighbors(net, vert)
            end
        end
        println("Total pruned edges: $nlost")

        # redistribute them elsewhere 
        vert = vertices(net)
        @Threads.threads for i in 1:nlost # will complain if there are not enough edges to add back in
            verts = vert[degree(net) .< (gth - 1)]
            if length(verts) < 2
                break
            end
            v1 = rand(verts)
            v2 = v1
            while v2 == v1
                v2 = rand(verts)
            end
            edg = Edge(v1,v2)
            add_edge!(net,edg)
            println("Redistributed edge $i/$nlost between $v1 and $v2")
        end
        println("Completed redistributing edges")  

        # maybe a few short, add those too
        iteration = 0
    	while ne(net) < nedges
            verts = vert[degree(net) .< (gth - 1)]
            if length(verts) < 2
                break
            end
    	    edg = Edge(rand(verts),rand(verts))
    	    while edg in edges(net)
                edg = Edge(rand(verts),rand(verts))
    	    end
    	    add_edge!(net,edg)
            iteration += 1
            println("Restored missing edge in iteration $iteration")
        end
        println("Final edge count restored to original: $nedges")
        return net
    end

    function limitpop_sbm(net, gth, groups::Dict{Int, Vector{Int}})
        # prune contact patterns to limit gatherings on network to no more than gth people. 
        # edges are redistributed to preserve mean degree
        nlost = 0
        nedges = ne(net)
        println("Starting edge pruning to limit gatherings to $gth people")

        # prune edges from hubs
        @Threads.threads for vert in vertices(net)
            neigh = neighbors(net, vert)
            while length(neigh) > gth
                rem_edge!(net, vert, rand(neigh))
                nlost += 1
                println("Pruned an edge from vertex $vert (nlost: $nlost)")
                neigh = neighbors(net, vert)
            end
        end
        println("Total pruned edges: $nlost")

        # redistribute them elsewhere 
        vert = vertices(net)
        group_map = Dict(v => gid for (gid, nodes) in groups for v in nodes) # for each node v, get its group id
        @Threads.threads for i in 1:round(Int, (0.694*nlost)) 
            verts = [v for v in vert if degree(net, v) < (gth - 1) && group_map[v] == 1] # only consider group 1 nodes for 69.4% of the lost edges
            if length(verts) < 2
                break
            end
            v1 = rand(verts)
            v2 = v1
            while v2 == v1
                v2 = rand(verts)
            end
            edg = Edge(v1,v2)
            add_edge!(net,edg)
            println("Redistributed edge $i/$nlost between $v1 and $v2")
        end
        @Threads.threads for i in 1:round(Int, (0.306*nlost)) 
            verts = [v for v in vert if degree(net, v) < (gth - 1) && group_map[v] == 2] # only consider group 2 nodes for 30.6% of the lost edges
            if length(verts) < 2
                break
            end
            v1 = rand(verts)
            v2 = v1
            while v2 == v1
                v2 = rand(verts)
            end
            edg = Edge(v1,v2)
            add_edge!(net,edg)
            println("Redistributed edge $i/$nlost between $v1 and $v2")
        end
        println("Completed redistributing edges")

        # maybe a few short, add those too 
        iteration = 0
        while ne(net) < nedges
            edge_added = false # track if an edge was added in this iteration
            for group_id in [1, 2]
                verts = [v for v in vert if degree(net, v) < (gth - 1) && group_map[v] == group_id]
                if length(verts) < 2
                    continue
                end
                vpair = shuffle(collect(combinations(verts, 2)))
                for (v1, v2) in vpair
                    if !has_edge(net, v1, v2) 
                        add_edge!(net, v1, v2)
                        iteration += 1
                        edge_added = true
                        println("Restored missing edge in iteration $iteration")
                        break
                    end
                end
                if ne(net) >= nedges
                    break
                end
            end
            if !edge_added
                println("No more valid edges to add. Final edge count: $(ne(net)). Target: $nedges.")
                break
            end
        end
        # Count inter-group edges
        inter_edges = 0
        for e in edges(net)
            v1, v2 = src(e), dst(e)
            if group_map[v1] != group_map[v2]
                inter_edges += 1
            end
        end
        println("Number of inter-group edges: $inter_edges")
        return net
    end

    function limitmix(pop,gth)
        # model a population of pop fully mixed nodes
        # a BA scale-free network model, mean degree 4
        # No more than gth people (put a cap on mass gatherings)
        println("Starting limitmix for population $pop with limit $gth")
        return limitpop(fullmixing(pop), gth)
    end

    function limitmix_for_sbm(net, gth, groups::Dict{Int, Vector{Int}})
        # model a population of pop fully mixed nodes
        # a BA scale-free network model, mean degree 4
        # No more than gth people
        println("Starting limitmix for SBM network with limit $gth")
        return limitpop_sbm(net, gth, groups)
    end

    function isolation(pop)
        # model a population of pop (must be a square number) nodes
        # on a 2D grid, mean degree 4
        println("Starting lattice model for a population of $pop")
        spop = Int(floor(sqrt(pop)))
        if abs(spop^2-pop) > 0
            spop1 = Int(floor(sqrt(pop)))  
            spop2 = Int(ceil(sqrt(pop)))
            if spop2*spop1 > spop
                spop2 = spop2-1
            end
            println("population not square ", pop-spop1*spop2, " node(s) isolated")
            net = SimpleGraph(Graphs.grid([spop1, spop2], periodic=true))
            add_vertices!(net,pop-spop1*spop2)
        else
            println("Population is a perfect square. Creating a grid of size $spop x $spop.")
            net = SimpleGraph(Graphs.grid([spop, spop], periodic=true))
        end
        println("Completed lattice model setup.")
        return net
    end


    function socialdist(pop, cmpl)
        # model a population of pop (must be a square number) nodes
        # on a 2D grid, with social distancing compliance cmpl, mean degree 4
        println("Starting social distancing model for a population of $pop with compliance $cmpl")
        q = 1-cmpl^(1/4)  # as elsewhere, this whole thing is assuming a mean degree of 4
        println("Rewiring probability, $q")
        net = isolation(pop)
        println("Network created. Applying rewiring for social distancing.")
        rewire!(net,q)
        println("Completed rewiring. Social distancing model setup finished.")
        return net
    end


    function rewire!(net, p)
        # rewire edges of a graph with probability p
        vert = vertices(net)
        
        for edg in collect(edges(net))
            if rand(1)[1] < p
                u = src(edg)
                # filter out the current node and nodes that are already connected
                candidates = filter(x -> x != u && !has_edge(net, u, x), vert) # ensures no self-loops and multiedges
                if !isempty(candidates)
                    v_new = rand(candidates)
                    rem_edge!(net, edg) # remove edge only if rewiring happens, if not leave the edge as is
                    add_edge!(net, u, v_new)
                end
            end
        end
        return net
    end


    # remap node indices in each community to match the global graph
    function remap!(g::SimpleGraph, new_ids::Vector{Int})
        mapping = Dict(i => new_ids[i] for i in 1:length(new_ids))
        g2 = SimpleGraph(length(new_ids))
        for e in edges(g)
            u = mapping[src(e)]
            v = mapping[dst(e)]
            add_edge!(g2, u, v)
        end
        return g2
    end

    function rewire_sbm!(net::SimpleGraph, comp::Float64, groups::Dict{Int, Vector{Int}})
        verts = collect(vertices(net))
        group_map = Dict(v => gid for (gid, nodes) in groups for v in nodes) # for each node v, get its group id

        for e in collect(edges(net)) 
            if rand(1)[1] < comp
                u = src(e)
                g = group_map[u]

                # obtain the candidates for rewiring
                candidates = rand(1)[1] < 0.694 ? 
                    filter(x -> x != u && !has_edge(net, u, x), groups[g]) :
                    filter(x -> group_map[x] != g && !has_edge(net, u, x), verts)

                if !isempty(candidates)
                    v_new = rand(candidates)
                    rem_edge!(net, e) # remove edge only if rewiring happens, if not leave the edge as is
                    add_edge!(net, u, v_new)
                end
            end
        end  
    end 

    # build a lattice model with community structure and
    # rewire edges (based on compliance levels during MCO)
    # 69.4% of the time, edges are rewired within the same group, and the rest between groups
    function socialdist_sbm(groups::Dict{Int, Vector{Int}}, cmpl::Float64)
        net = SimpleGraph(sum(length, values(groups))) # create an empty graph with the total number of nodes in both groups

        # build community-wise lattice networks 
        for (gid, nodes) in groups
            subpopulation = length(nodes)
            subgraph = isolation(subpopulation)
            g_sub = remap!(subgraph, nodes) # remap local node indices to match the global graph
            for e in edges(g_sub)
                add_edge!(net, src(e), dst(e)) 
            end 
        end

        # apply rewiring based on compliance
        q = 1 - cmpl^(1/4)
        rewire_sbm!(net, q, groups)
        return net
    end

end  


