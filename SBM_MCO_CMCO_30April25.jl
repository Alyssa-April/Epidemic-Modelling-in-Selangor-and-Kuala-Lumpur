#________________________________________________________________
##### Include the relevant modules #####
#________________________________________________________________
include("ContactModels_Msia.jl") 
include("EpiSim_Msia.jl")

import Pkg
Pkg.activate("new_environment")
Pkg.add("FFMPEG")

# load all the relevant packages
using GraphPlot
using Compose
using Colors
using Graphs
using SimpleWeightedGraphs
using GraphRecipes
using Plots
using Statistics
using NetworkLayout
using Serialization
using Combinatorics

# this part is to plot the real data by KKM
using DataFrames, CSV, Dates 

#________________________________________________________________________
##### Extract the active cases data from Selangor and Kuala Lumpur based 
# on a certain timeframe #####
#_________________________________________________________________________

# Load the csv file into vscode
data = CSV.read("KKM_cases_bystate.csv", DataFrame)

# Filter for Selangor and W.P. Kuala Lumpur only
selangor_kl = filter(row -> row.state in ["Selangor", "W.P. Kuala Lumpur"], data)

# Group by date and sum the columns to get an aggregate
sel_kl_data = combine(groupby(selangor_kl, :date), 
[:cases_new, :cases_import, :cases_recovered, :cases_active] .=> sum)

# Convert the date column to Date type
sel_kl_data.date = Date.(sel_kl_data.date, "m/d/yyyy")

# Define the start and end dates for filtering (use these range of dates)
start_date = Date("2020-02-04") 
end_date = Date("2020-06-09")

# Filter the DataFrame to include rows within the specified date range
final_data = filter(row -> start_date <= row.date <= end_date, sel_kl_data)

# to view the first few rows of the data
first(final_data, 4)
last(final_data.cases_active_sum, 4)

# Extract dates and the column of interest (e.g., total_new_vax for total new vaccinations)
dates = final_data.date
new_cases = final_data.cases_active_sum
recovered_cases = final_data.cases_recovered_sum 
 
#__________________________________________________________________________________________________
##### Set the parameters for SEIRS epistep function #####
#__________________________________________________________________________________________________

# set the parameters for the scale-free network 
epiparam = Dict()
transmission_rates = [0.063]  # the final transmission rate used in the paper
epiparam["β_intra"] = 0.063 # transmission rate for intra-community transmission in SBM baseline model
epiparam["β_inter"] = 0.0315 # transmission rate for inter-community transmission in SBM baseline model
epiparam["σ"] = 1/5 # 5 day latency period
epiparam["γ"] = 1/11 # 11 day infectious period (recovery rate)
epiparam["ξ"] = 0.0345 # reinfection rate
epiparam["nseeds"] = 2 # only 2 active case in Selangor and KL on 4 February 2020 (infected)

nsims = 1000 # number of simulations

#______________________________________________________________________________________________
##### BA, MCO, CMCO (4 Feb 2020 to 9 June 2020) - 127 days #####
#_______________________________________________________________________________________________

ndays = 127


epiparam["pop"] = 1024 # population size (change size according to simulation)
#epiparam["pop"] = 10000
#epiparam["pop"] = 100489

using Random

function assign_ethnic_group_to_nodes(num_nodes, proportions)
    ethnic_groups = []
    for (i, prop) in enumerate(proportions)
        append!(ethnic_groups, fill(i, round(Int, num_nodes * prop)))
    end
    shuffle!(ethnic_groups)
    return ethnic_groups
end

# Assign ethnic groups

num_nodes = 1024
#num_nodes = 10000
#num_nodes = 100489
proportions = [0.694, 0.306] # Majority (bumiputera) and Minority (Chinese, Indian, Others) proportions
ethnic_group = assign_ethnic_group_to_nodes(num_nodes, proportions)

# check the proportions of each ethnicity
unique_ethinicity = unique(ethnic_group)
counts = Dict(ethnic => count(x -> x == ethnic, ethnic_group) for ethnic in unique_ethinicity)

# number of nodes in the Majority and Minority groups respectively
number_per_gp = [711, 313] # change proportions according to number of nodes
#number_per_gp = [6940, 3060]
#number_per_gp = [69739, 30750]

# the probabilities of edge formation between and within ethnic groups
pr_mat = [
    0.071 0.012;
    0.012 0.156
]


groups = Dict(
    1 => collect(1:711),
    2 => collect(712:1024) 
)


#=
groups = Dict(
    1 => collect(1:6940),
    2 => collect(6941:10000)  
)


groups = Dict(
    1 => collect(1:69739),
    2 => collect(69740:100489) 
)
=#


# use SBM then limit ten for CMCO (cap on mass gatherings)
SBM = ContactModels_Msia.generate_sbm(number_per_gp, pr_mat, 1234)
serialize("networks_SBMBASELINE_A_1k.dat", (SBM))

SBM_CMCO = ContactModels_Msia.limitmix_for_sbm(SBM, 10, groups)
serialize("networks_SBMBASELINE_B_1k.dat", (SBM_CMCO))


MCOCMCO_95 = ContactModels_Msia.socialdist_sbm(groups, 0.95)
serialize("networks_SBMBASELINE_C_1k.dat", (MCOCMCO_95))
MCOCMCO_90 = ContactModels_Msia.socialdist_sbm(groups, 0.9)
serialize("networks_SBMBASELINE_D_1k.dat", (MCOCMCO_90))
MCOCMCO_85 = ContactModels_Msia.socialdist_sbm(groups, 0.85)
serialize("networks_SBMBASELINE_E_1k.dat", (MCOCMCO_85))
MCOCMCO_80 = ContactModels_Msia.socialdist_sbm(groups, 0.8)
serialize("networks_SBMBASELINE_F_1k.dat", (MCOCMCO_80))
MCOCMCO_70 = ContactModels_Msia.socialdist_sbm(groups, 0.7)
serialize("networks_SBMBASELINE_G_1k.dat", (MCOCMCO_70))
MCOCMCO_60 = ContactModels_Msia.socialdist_sbm(groups, 0.6)
serialize("networks_SBMBASELINE_H_1k.dat", (MCOCMCO_60))
MCOCMCO_50 = ContactModels_Msia.socialdist_sbm(groups, 0.5)
serialize("networks_SBMBASELINE_I_1k.dat", (MCOCMCO_50))
MCOCMCO_40 = ContactModels_Msia.socialdist_sbm(groups, 0.4)
serialize("networks_SBMBASELINE_J_1k.dat", (MCOCMCO_40))
MCOCMCO_30 = ContactModels_Msia.socialdist_sbm(groups, 0.3)
serialize("networks_SBMBASELINE_K_1k.dat", (MCOCMCO_30))
MCOCMCO_20 = ContactModels_Msia.socialdist_sbm(groups, 0.2)
serialize("networks_SBMBASELINE_L_1k.dat", (MCOCMCO_20))
MCOCMCO_10 = ContactModels_Msia.socialdist_sbm(groups, 0.1)
serialize("networks_SBMBASELINE_M_1k.dat", (MCOCMCO_10))


# Save to a file
serialize("networks_SBMBASELINE_1k.dat", (SBM, SBM_CMCO, MCOCMCO_95, MCOCMCO_90, MCOCMCO_85, 
MCOCMCO_80, MCOCMCO_70, MCOCMCO_60, MCOCMCO_50, MCOCMCO_40, MCOCMCO_30, MCOCMCO_20, MCOCMCO_10))


# Load from the file
(SBM, SBM_CMCO, MCOCMCO_95, MCOCMCO_90, MCOCMCO_85, 
MCOCMCO_80, MCOCMCO_70, MCOCMCO_60, MCOCMCO_50, MCOCMCO_40,
MCOCMCO_30, MCOCMCO_20, MCOCMCO_10) = deserialize("networks_SBMBASELINE_1k.dat")


# -------------------------------------------------------------------------------------
#### Plot the lineplots for SBM Baseline with 1024/10000/100489 population

results_St_SBMBASELINE_1k = Dict()
results_Et_SBMBASELINE_1k = Dict()
results_It_SBMBASELINE_1k = Dict()
results_Rt_SBMBASELINE_1k = Dict()
results_SBM_SBMBASELINE_1k = Dict()

for β in transmission_rates
    epiparam["β_intra"] = β
    epiparam["β_inter"] = β/2

    St95,Et95,It95,Rt95 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(SBM, MCOCMCO_95, SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St90,Et90,It90,Rt90 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(SBM, MCOCMCO_90, SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St85,Et85,It85,Rt85 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(SBM, MCOCMCO_85, SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St80,Et80,It80,Rt80 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(SBM, MCOCMCO_80, SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St70,Et70,It70,Rt70 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(SBM, MCOCMCO_70, SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St60,Et60,It60,Rt60 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(SBM, MCOCMCO_60, SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St50,Et50,It50,Rt50 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(SBM, MCOCMCO_50, SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St40,Et40,It40,Rt40 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(SBM, MCOCMCO_40, SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St30,Et30,It30,Rt30 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(SBM, MCOCMCO_30, SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St20,Et20,It20,Rt20 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(SBM, MCOCMCO_20, SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St10,Et10,It10,Rt10 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(SBM, MCOCMCO_10, SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)

    St,Et,It,Rt = EpiSim_Msia.episim_homophily_nochangepoint(SBM, epiparam, ethnic_group, ndays, nsims)

    sbm_vals = Dict("SBM_St" => St, "SBM_Et" => Et, "SBM_It" => It, "SBM_Rt" => Rt)
    results_SBM_SBMBASELINE_1k[(β)] = sbm_vals

    St_vals = Dict("St_95" => St95, "St_90" => St90, "St_85" => St85, "St_80" => St80, "St_70" => St70, "St_60" => St60,
    "St_50" => St50, "St_40" => St40, "St_30" => St30, "St_20" => St20, "St_10" => St10)

    Et_vals = Dict("Et_95" => Et95, "Et_90" => Et90, "Et_85" => Et85, "Et_80" => Et80, "Et_70" => Et70, "Et_60" => Et60,
    "Et_50" => Et50, "Et_40" => Et40, "Et_30" => Et30, "Et_20" => Et20, "Et_10" => Et10)

    It_vals = Dict("It_95" => It95, "It_90" => It90, "It_85" => It85, "It_80" => It80, "It_70" => It70, "It_60" => It60,
    "It_50" => It50, "It_40" => It40, "It_30" => It30, "It_20" => It20, "It_10" => It10)

    Rt_vals = Dict("Rt_95" => Rt95, "Rt_90" => Rt90, "Rt_85" => Rt85, "Rt_80" => Rt80, "Rt_70" => Rt70, "Rt_60" => Rt60,
    "Rt_50" => Rt50, "Rt_40" => Rt40, "Rt_30" => Rt30, "Rt_20" => Rt20, "Rt_10" => Rt10)

    for compliance in [95, 90, 85, 80, 70, 60, 50, 40, 30, 20, 10]
        results_St_SBMBASELINE_1k[(β, compliance)] = St_vals["St_$(compliance)"]
        results_Et_SBMBASELINE_1k[(β, compliance)] = Et_vals["Et_$(compliance)"]
        results_It_SBMBASELINE_1k[(β, compliance)] = It_vals["It_$(compliance)"]
        results_Rt_SBMBASELINE_1k[(β, compliance)] = Rt_vals["Rt_$(compliance)"]
    end

    plot(title = "SBM Baseline Model - 1024 Nodes", 
    ylabel = "Number of Active Cases", legend = :topleft,
    legendfontpointsize = 5, titlefontsize = 10)
    EpiSim_Msia.plotqnt_normal(1:127, It95, :grey, "95%", 0.0)
    EpiSim_Msia.plotqnt_normal(1:127, It90, :yellow, "90%", 0.0)
    EpiSim_Msia.plotqnt_normal(1:127, It85, :blue, "85%", 0.0)
    EpiSim_Msia.plotqnt_normal(1:127, It80, :purple, "80%", 0.0)
    EpiSim_Msia.plotqnt_normal(1:127, It70, :green, "70%", 0.0)
    EpiSim_Msia.plotqnt_normal(1:127, It60, :orange, "60%", 0.0)
    EpiSim_Msia.plotqnt_normal(1:127, It50, :pink, "50%", 0.0)
    EpiSim_Msia.plotqnt_normal(1:127, It40, :black, "40%", 0.0)
    EpiSim_Msia.plotqnt_normal(1:127, It30, :lightblue, "30%", 0.0)
    EpiSim_Msia.plotqnt_normal(1:127, It20, :lightgreen, "20%", 0.0)
    EpiSim_Msia.plotqnt_normal(1:127, It10, :mistyrose, "10%", 0.0)

    #EpiSim_Msia.plotqnt_normal(1:127, It, :red, "BA", 0.0)

    RealData_NormalScale_SBMBASELINE = plot!(1:127, new_cases,
        xlabel = "Days", ylabel = "Number of Active Cases", 
        label = "Real Data", linecolor = :brown, linewidth = 1.5,
        xtickfontsize = 8, grid = false)

    savefig(RealData_NormalScale_SBMBASELINE, "SBM Baseline Model - 1024 Nodes (Normal).png")

    plot(title = "SBM Baseline Model - 1024 Nodes", 
    ylabel = "Number of Active Cases", legend = :topleft, yaxis = :log,
    legendfontpointsize = 5, titlefontsize = 10)
    EpiSim_Msia.plotqnt_log(1:127, It95, :grey, "95%", 0.0)
    EpiSim_Msia.plotqnt_log(1:127, It90, :yellow, "90%", 0.0)
    EpiSim_Msia.plotqnt_log(1:127, It85, :blue, "85%", 0.0)
    EpiSim_Msia.plotqnt_log(1:127, It80, :purple, "80%", 0.0)
    EpiSim_Msia.plotqnt_log(1:127, It70, :green, "70%", 0.0)
    EpiSim_Msia.plotqnt_log(1:127, It60, :orange, "60%", 0.0)
    EpiSim_Msia.plotqnt_log(1:127, It50, :pink, "50%", 0.0)
    EpiSim_Msia.plotqnt_log(1:127, It40, :black, "40%", 0.0)
    EpiSim_Msia.plotqnt_log(1:127, It30, :lightblue, "30%", 0.0)
    EpiSim_Msia.plotqnt_log(1:127, It20, :lightgreen, "20%", 0.0)
    EpiSim_Msia.plotqnt_log(1:127, It10, :mistyrose, "10%", 0.0)

    #EpiSim_Msia.plotqnt_log(1:127, It, :red, "BA", 0.0)

    RealData_LogScale_SBMBASELINE = plot!(1:127, new_cases,
        xlabel = "Days", ylabel = "Number of Active Cases", 
        label = "Real Data", linecolor = :brown, linewidth = 1.5,
        xtickfontsize = 8, grid = false, yaxis = :log)

    savefig(RealData_LogScale_SBMBASELINE, "SBM Baseline Model - 1024 Nodes (Log).png")
end

# Define file names to save each dictionary
file_St_SBMBASELINE_1k = "results_St_SBMBASELINE_1k.jls"
file_Et_SBMBASELINE_1k = "results_Et_SBMBASELINE_1k.jls"
file_It_SBMBASELINE_1k = "results_It_SBMBASELINE_1k.jls"
file_Rt_SBMBASELINE_1k = "results_Rt_SBMBASELINE_1k.jls"
file_SBM_SBMBASELINE_1k = "results_SBM_SBMBASELINE_1k.jls"

# Save the dictionaries to files (Can reload results for future use)
serialize(file_St_SBMBASELINE_1k, results_St_SBMBASELINE_1k)
serialize(file_Et_SBMBASELINE_1k, results_Et_SBMBASELINE_1k)
serialize(file_It_SBMBASELINE_1k, results_It_SBMBASELINE_1k)
serialize(file_Rt_SBMBASELINE_1k, results_Rt_SBMBASELINE_1k)
serialize(file_SBM_SBMBASELINE_1k, results_SBM_SBMBASELINE_1k)

#___________________________________________________________________________________________________

# Load the dictionaries from files
results_St_SBMBASELINE_1k = deserialize("results_St_SBMBASELINE_1k.jls")
results_Et_SBMBASELINE_1k = deserialize("results_Et_SBMBASELINE_1k.jls")
results_It_SBMBASELINE_1k = deserialize("results_It_SBMBASELINE_1k.jls")
results_Rt_SBMBASELINE_1k = deserialize("results_Rt_SBMBASELINE_1k.jls")
results_SBM_SBMBASELINE_1k = deserialize("results_SBM_SBMBASELINE_1k.jls")

#___________________________________________________________________________________________________

#### Malaysia's Compliance rate

# chosen tranmission rate 
chosen_transmission_rate = 0.063

# Filter results_It_SBMBASELINE_10k for the chosen transmission rate.
chosen_results_It = Dict(key => value for (key, value) in results_It_SBMBASELINE_10k 
if key[1] == chosen_transmission_rate)

# Real data for comparison
real_data = new_cases  

# Calculate MSE and MAE for each compliance rate
# The smallest errors show the compliance rate closest to the real data
mse_values = Dict()
mae_values = Dict()

for (comp_rate, sim_cases) in chosen_results_It
    log_real_data = log.(real_data) # use log number of active cases for better comparison
    log_sim_cases = log.(sim_cases)
    
    # Mean Squared Error (NSE)
    mse = mean((log_real_data .- log_sim_cases).^2)
    
    # Mean Absolute Error (MAE)
    mae = mean(abs.(log_real_data .- log_sim_cases))
    
    mse_values[comp_rate] = mse
    mae_values[comp_rate] = mae
end

for comp in keys(mse_values)
    println(comp)
end

# Find the compliance rate that is closest to the real_data
min_mse_comp = minimum([(mse_values[comp], comp) for comp in keys(mse_values)])
min_mae_comp = minimum([(mae_values[comp], comp) for comp in keys(mae_values)])

println("Compliance rate closest to real data based on MSE: $(min_mse_comp[2][2])% with MSE value $(min_mse_comp[1])")
println("Compliance rate closest to real data based on MAE: $(min_mae_comp[2][2])% with MAE value $(min_mae_comp[1])")

#___________________________________________________________________________________________________
#___________________________________________________________________________________________________

# Quantifying the effectiveness of control measures by comparing the sum of simulated active cases
# throughout all 127 time steps for a (transmission rate of 0.063 and MCO compliance rate of 85%) 
# to the scale-free network.
# Also calculate at day 90

# This is the simulation of active cases for the SBM model (Baseline with no control measures at all)
sim_SBM_cases_no_control = results_SBM_SBMBASELINE_10k[(0.063)]["SBM_It"] 

# min_mae_compliance[2][2] represents the chosen compliance rate for the chosen transmission rate
chosen_comp_rate = min_mae_comp[2][2]

# Filter results for 6.3% transmission rate at 60% compliance
chosen_results_It85 = Dict(key => value for (key, value) in results_It_SBMBASELINE_10k 
if key == (chosen_transmission_rate, chosen_comp_rate))

# Calculate Percentage Reduction
sim_cases_85 = chosen_results_It85[chosen_transmission_rate, chosen_comp_rate]

# Plot the SBM and 85% compliance line first 
plot(title = "Active Cases Based on Simulations (With and Without Control Measures)", 
ylabel = "Number of Active Cases", xlabel = "Days", legend = :topleft, yaxis = :log,
legendfontpointsize = 7, titlefontsize = 10)

EpiSim_Msia.plotqnt_log(1:127, sim_cases_85, :blue, "With Control Measures", 0.0)

sbm_85_plot = EpiSim_Msia.plotqnt_log(1:127, sim_SBM_cases_no_control, :red, "Without Control Measures", 0.0)

savefig(sbm_85_plot, "w_and_wo_control_85_plot.png")

# Since there are 1000 simulations, compute the average across simulations for each time step
median_sim_It85 = median(sim_cases_85, dims = 2)  # Averages across columns (simulations)
median_sbm = median(sim_SBM_cases_no_control, dims = 2)

percentage_reduction = (1 - sum(median_sim_It85) / sum(median_sbm)) * 100

println("Percentage Reduction in cases: $(percentage_reduction)%")

# Find the reduction at time step 90 (when MCO ends)

# number of cases at time step 90
cases_sim_It85_day90 = median_sim_It85[90]
cases_sbm_day90 = median_sbm[90]

percentage_reduction_day90 = (1 - cases_sim_It85_day90 / cases_sbm_day90) * 100

println("Percentage Reduction in cases on day 90: $(percentage_reduction_day90)%")


#___________________________________________________________________________________________________
#___________________________________________________________________________________________________

#### Peak Calculation

# to set the position of the control measure labels
function relative(p::Plots.Subplot, rx, ry) 
    xlims = Plots.xlims(p)
    ylims = Plots.ylims(p)
    return xlims[1] + rx * (xlims[2]-xlims[1]), ylims[1] + ry * (ylims[2] - ylims[1])
 end

# Get the number of time steps (x-axis) 
time_steps = 1:size(median_sim_It85, 1) # Assuming one row per day

# Find peaks
peak_85 = Int(maximum(median_sim_It85))
peak_85_day = argmax(median_sim_It85)[1] # obtain the peak day

peak_sbm = Int(maximum(median_sbm))
peak_sbm_day = argmax(median_sbm)[1] 

# Plotting
pl = plot(time_steps, median_sim_It85, label = "With Control Measures", grid = false,
color = :blue, lw = 2, xlabel = "Days", ylabel = "Number of Active Cases",
xtickfontsize = 8, legendfontsize = 7)
plot!(time_steps, median_sbm, label = "Without Control Measures", color = :red, lw = 2)


# Annotate peaks
annotate!(sp = 1, relative(pl[1], 0.06, 0.4)..., Plots.text("Peak: $peak_85 cases\nDay: $peak_85_day", :left, 10, :blue))
annotate!(sp = 1, relative(pl[1], 0.06, 0.6)..., Plots.text("Peak: $peak_sbm cases\nDay: $peak_sbm_day", :left, 10, :red))

# Add title for plot
title!("Active Cases Based on Simulations (With and Without Control Measures)", titlefont = 9)

# Save the plot
savefig("case_peaks.png")

# Calculate the reduction in peak cases
reduction_percentage = ((peak_sbm - peak_85) / peak_sbm) * 100 

#________________________________________________________________

# Plot the SBM Network

using GraphMakie, CairoMakie

G = MCOCMCO_80

# set the colours for the Majority and Minority ethnic groups
ethnic_colors = Dict(1 => :dodgerblue, 2 => :orange)  

# Create the node_colors vector based on ethnic group
node_colours = [ethnic_colors[ethnic_group[i]] for i in 1:nv(G)]

# Plot the graph
fig, ax, p = GraphMakie.graphplot(G, node_size = 6, node_color = node_colours)

# Save and display the plot
save("Plot SBM Baseline Model - 1024 Nodes (MCOCMCO_80).png", fig)

fig