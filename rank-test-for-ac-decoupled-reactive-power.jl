# Copyright 2021 Mark Jeeninga 
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Email: mark<dot>jeeninga<at>gmail<dot>com



# This code considers all test cases in MATPOWER 7.0.
# A DC power grid is derived from each test case by considering the decoupled reactive power flow.
# This means that we take the branch inductances in each case as the line resistances in the DC grid.
# In this analogy, PQ buses correspond to constant-power loads and PV buses correspond to voltages sources.
#
# The goal of this code is to answer the following question:
# Does there exist a connected component in the induced subgraph of PQ buses with at least two nodes
# such that the Laplacian submatrix that describes the inductive connections between these PQ buses and the 
# neighboring PV buses has full column rank?
# If so, then any increase in the edges within this subgraph could violate the power flow feasibility for
# some choice of power demands at the loads, and hence Braess' paradox for power flow feasibilty can occur.

using LinearAlgebra
using SparseArrays

# This codes makes use of the LightGraphs and PowerModels packages, and can be added using the Julia package manager.
# https://github.com/sbromberger/LightGraphs.jl
# https://lanl-ansi.github.io/PowerModels.jl/stable/
using LightGraphs
import PowerModels

# A list of all the MATPOWER test cases, available at
# https://matpower.org
Matpower_cases = [
    "data/case118.m",
    "data/case1354pegase.m",
    "data/case13659pegase.m",
    "data/case14.m",
    "data/case141.m",
    "data/case145.m",
    "data/case18.m",
    "data/case1888rte.m",
    "data/case1951rte.m",
    "data/case22.m",
    "data/case2383wp.m",
    "data/case24_ieee_rts.m",
    "data/case2736sp.m",
    "data/case2737sop.m",
    "data/case2746wop.m",
    "data/case2746wp.m",
    "data/case2848rte.m",
    "data/case2868rte.m",
    "data/case2869pegase.m",
    "data/case30.m",
    "data/case300.m",
    "data/case3012wp.m",
    "data/case30pwl.m",
    "data/case30Q.m",
    "data/case3120sp.m",
    "data/case3375wp.m",
    "data/case33bw.m",
    "data/case39.m",
    "data/case4_dist.m",
    "data/case4gs.m",
    "data/case5.m",
    "data/case57.m",
    "data/case6468rte.m",
    "data/case6470rte.m",
    "data/case6495rte.m",
    "data/case6515rte.m",
    "data/case69.m",
    "data/case6ww.m",
    "data/case85.m",
    "data/case89pegase.m",
    "data/case9.m",
    "data/case9241pegase.m",
    "data/case9Q.m",
    "data/case9target.m",
    "data/case_ACTIVSg10k.m",
    "data/case_ACTIVSg200.m",
    "data/case_ACTIVSg2000.m",
    "data/case_ACTIVSg25k.m",
    "data/case_ACTIVSg500.m",
    "data/case_ACTIVSg70k.m",
    "data/case_ieee30.m",
    "data/case_RTS_GMLC.m",
    "data/case_SyntheticUSA.m"
]

# Suppress the warnings by PowerModels when loading the MATPOWER test cases
PowerModels.silence()

# Perform the test for each MATPOWER case (mpc)
for mpc in Matpower_cases

    # Load the test case into PowerModels
    network = PowerModels.parse_file(mpc);

    # Determine the number of buses and branches
    nr_buses = length(network["bus"])
    nr_branches = length(network["branch"])

    # Collect the keys of the nodes and branches.
    # We need to do this since PowerModels does not use sequential numbers 
    # to identify the nodes, but instead labels them with a string.
    # We use the indices in these arrays to identify the nodes and branches.
    bus_keys = collect(keys(network["bus"]))
    branch_keys = collect(keys(network["branch"]))

    # Define the graph and the weighted laplacian matrix associated to the power grid
    graph = SimpleGraph(nr_buses);
    Laplacian = spzeros(nr_buses, nr_buses);

    # Add the edges to the graph and fill the Laplacian matrix based on the branch data
    for (branchKey, branch) in network["branch"]
        f_bus = string(branch["f_bus"]);
        t_bus = string(branch["t_bus"]);

        # Since PowerModels keeps track of the bus names using strings, we have to
        # look up to the index of each node that the string is referring to
        i = findfirst(x -> x == f_bus, bus_keys)
        j = findfirst(x -> x == t_bus, bus_keys)

        add_edge!(graph, i, j);

        # The branch susceptence
        # Since we consider lossless AC power grids we ignore the resistance branch["br_r"]
        w = 1/branch["br_x"];

        # Add the weight to the Laplacian matrix
        Laplacian[i,j] -= w;
        Laplacian[j,i] -= w;
        Laplacian[i,i] += w;
        Laplacian[j,j] += w;
    end

    # Filter the buses based on their type
    PQ_buses = [];
    PV_buses = [];
    reference_buses = [];
    isolated_buses = [];

    # Iterate over the buses and sort them based on their type (PQ = 1, PV = 2, reference = 3, isolated = 4)
    for key_id in 1:length(bus_keys)
        # key_id is our indexing of the buses.
        # We convert this to the indexing of PowerModels by using bus_keys.
        # In this way we identify the branch in the PowerModels data and its bus type.
        key = bus_keys[key_id];
        bus = network["bus"][key];
        bus_type = bus["bus_type"];

        # Add the bus to the array corresponding to its type
        if bus_type == 1
            push!(PQ_buses, key_id)
        elseif bus_type == 2
            push!(PV_buses, key_id)
        elseif bus_type == 3
            push!(reference_buses, key_id)
        elseif bus_type == 4
            push!(isolated_buses, key_id)
        end
    end

    # Make sure there are PV buses in the network, and skip this case otherwise.
    if length(PV_buses) == 0
        println(mpc, " has no PV buses and was skipped.")
        continue;
    end

    # Partition the subgraph induced by the PQ buses into connected components.
    PQ_subgraph, PQ_vector_mapping = induced_subgraph(graph, [PQ_buses...])
    PQ_connected_components = connected_components(PQ_subgraph);

    # pass is false unless we find otherwise.
    pass = false;

    # For each connected component we consider we include the neighbors which are PV buses.
    for component in PQ_connected_components
        # The indices of the nodes in this component that correspond to the buses in the original graph.
        original_indices = PQ_vector_mapping[component]

        # We are only interested in the connected components that have at least two PQ buses.
        if length(component) >= 2
            # We list the neighbors in the original graph of the PQ buses in this connected component.
            PQ_component_neighbors = [];
            for i in original_indices
                PQ_component_neighbors = [PQ_component_neighbors..., neighbors(graph, i)...]
            end
            unique!(PQ_component_neighbors);

            # We are only interested in cases where the neighbors that are PV buses
            # since this list of neighbors also includes the PQ buses.
            intersect!(PQ_component_neighbors, PV_buses);

            # We skip this subgraph if there is only one PV neighbor, since in this case braess paradox does not occurs
            if length(PQ_component_neighbors) < 2
                continue
            end

            # We want to verify the claim that this submatrix of the Laplacian matrix has full column rank.
            my_rank = rank(Laplacian[original_indices, PQ_component_neighbors]);
            my_goal = length(PQ_component_neighbors);

            # Check if the rank is equal to the goal.
            # If this is true then we break the loop over the PQ components and move to the next MATPOWER case.
            if my_rank == my_goal
                pass = true
                break
            end
        end
    end

    # Print wheter the MATPOWER case passed the test or not
    if pass
        println(mpc, "\tpassed the test.");
    else
        println(mpc, "\tfailed the test.");
    end

end
