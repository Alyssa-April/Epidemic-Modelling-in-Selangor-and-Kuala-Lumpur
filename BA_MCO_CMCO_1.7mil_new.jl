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
transmission_rates = [0.063] # the final transmission rate used in the paper for BA baseline model
epiparam["β_intra"] = 0.063 # transmission rate for intra-community transmission in modified networks based on SBM
epiparam["β_inter"] = 0.0315 # transmission rate for inter-community transmission in modified networks based on SBM
epiparam["σ"] = 1/5 # 5 day latency period
epiparam["γ"] = 1/11 # 11 day infectious period (recovery rate)
epiparam["ξ"] = 0.0345 # reinfection rate
epiparam["nseeds"] = 2 # only 2 active case in Selangor and KL on 4 February 2020 (infected)
epiparam["pop"] = 1747684 # population Kuala Lumpur rounded up to be a square number

nsims = 1000 # number of simulations 

#______________________________________________________________________________________________
##### BA, MCO, CMCO (4 Feb 2020 to 9 June 2020) - 127 days #####
#_______________________________________________________________________________________________

ndays = 127 # number of time steps

BA = ContactModels_Msia.fullmixing(1747684)
serialize("networks_BABASELINE_A_1M.dat", (BA))
CMCO = ContactModels_Msia.limitmix(1747684, 10) # use BA then limit ten for CMCO (cap on mass gatherings)
serialize("networks_BABASELINE_B_1M.dat", (CMCO))

BA_MCOCMCO_95 = ContactModels_Msia.socialdist(1747684, 0.95)
serialize("networks_BABASELINE_C_1M.dat", (BA_MCOCMCO_95)) 
BA_MCOCMCO_90 = ContactModels_Msia.socialdist(1747684, 0.9)
serialize("networks_BABASELINE_D_1M.dat", (BA_MCOCMCO_90))
BA_MCOCMCO_85 = ContactModels_Msia.socialdist(1747684, 0.85)
serialize("networks_BABASELINE_E_1M.dat", (BA_MCOCMCO_85))
BA_MCOCMCO_80 = ContactModels_Msia.socialdist(1747684, 0.8)
serialize("networks_BABASELINE_F_1M.dat", (BA_MCOCMCO_80))
BA_MCOCMCO_70 = ContactModels_Msia.socialdist(1747684, 0.7)
serialize("networks_BABASELINE_G_1M.dat", (BA_MCOCMCO_70))
BA_MCOCMCO_60 = ContactModels_Msia.socialdist(1747684, 0.6)
serialize("networks_BABASELINE_H_1M.dat", (BA_MCOCMCO_60))
BA_MCOCMCO_50 = ContactModels_Msia.socialdist(1747684, 0.5)
serialize("networks_BABASELINE_I_1M.dat", (BA_MCOCMCO_50))
BA_MCOCMCO_40 = ContactModels_Msia.socialdist(1747684, 0.4)
serialize("networks_BABASELINE_J_1M.dat", (BA_MCOCMCO_40))
BA_MCOCMCO_30 = ContactModels_Msia.socialdist(1747684, 0.3)
serialize("networks_BABASELINE_K_1M.dat", (BA_MCOCMCO_30))
BA_MCOCMCO_20 = ContactModels_Msia.socialdist(1747684, 0.2)
serialize("networks_BABASELINE_L_1M.dat", (BA_MCOCMCO_20))
BA_MCOCMCO_10 = ContactModels_Msia.socialdist(1747684, 0.1)
serialize("networks_BABASELINE_M_1M.dat", (BA_MCOCMCO_10))


# if need to deserialize the networks
(BA) = deserialize("networks_BABASELINE_A_1M.dat")
(CMCO) = deserialize("networks_BABASELINE_B_1M.dat")
(BA_MCOCMCO_95) = deserialize("networks_BABASELINE_C_1M.dat")
(BA_MCOCMCO_90) = deserialize("networks_BABASELINE_D_1M.dat")
(BA_MCOCMCO_85) = deserialize("networks_BABASELINE_E_1M.dat")
(BA_MCOCMCO_80) = deserialize("networks_BABASELINE_F_1M.dat")
(BA_MCOCMCO_70) = deserialize("networks_BABASELINE_G_1M.dat")
(BA_MCOCMCO_60) = deserialize("networks_BABASELINE_H_1M.dat")
(BA_MCOCMCO_50) = deserialize("networks_BABASELINE_I_1M.dat")
(BA_MCOCMCO_40) = deserialize("networks_BABASELINE_J_1M.dat")
(BA_MCOCMCO_30) = deserialize("networks_BABASELINE_K_1M.dat")
(BA_MCOCMCO_20) = deserialize("networks_BABASELINE_L_1M.dat")
(BA_MCOCMCO_10) = deserialize("networks_BABASELINE_M_1M.dat")

# Save to a file
serialize("networks_BABASELINE_1M.dat", (BA, CMCO, BA_MCOCMCO_95, BA_MCOCMCO_90, BA_MCOCMCO_85, 
BA_MCOCMCO_80, BA_MCOCMCO_70, BA_MCOCMCO_60, BA_MCOCMCO_50, BA_MCOCMCO_40, BA_MCOCMCO_30, 
BA_MCOCMCO_20, BA_MCOCMCO_10))

# Load from the file
(BA, CMCO, BA_MCOCMCO_95, BA_MCOCMCO_90, BA_MCOCMCO_85, BA_MCOCMCO_80, BA_MCOCMCO_70, 
BA_MCOCMCO_60, BA_MCOCMCO_50, BA_MCOCMCO_40, BA_MCOCMCO_30, 
BA_MCOCMCO_20, BA_MCOCMCO_10) = deserialize("networks_BABASELINE_1M.dat")

#___________________________________________________________________________________________________

#### Plot the lineplots for BA Baseline with 1747684 population

results_St_BABASELINE_1M = Dict()
results_Et_BABASELINE_1M = Dict()
results_It_BABASELINE_1M = Dict()
results_Rt_BABASELINE_1M = Dict()
results_BA_BABASELINE_1M = Dict()

for β in transmission_rates
    epiparam["β_fixed"] = β # only one transmission rate, no community structure

    St95,Et95,It95,Rt95 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_95, CMCO, epiparam, 44, 91, 36, nsims)
    St90,Et90,It90,Rt90 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_90, CMCO, epiparam, 44, 91, 36, nsims)
    St85,Et85,It85,Rt85 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_85, CMCO, epiparam, 44, 91, 36, nsims)
    St80,Et80,It80,Rt80 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_80, CMCO, epiparam, 44, 91, 36, nsims)
    St70,Et70,It70,Rt70 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_70, CMCO, epiparam, 44, 91, 36, nsims)
    St60,Et60,It60,Rt60 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_60, CMCO, epiparam, 44, 91, 36, nsims)
    St50,Et50,It50,Rt50 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_50, CMCO, epiparam, 44, 91, 36, nsims)
    St40,Et40,It40,Rt40 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_40, CMCO, epiparam, 44, 91, 36, nsims)
    St30,Et30,It30,Rt30 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_30, CMCO, epiparam, 44, 91, 36, nsims)
    St20,Et20,It20,Rt20 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_20, CMCO, epiparam, 44, 91, 36, nsims)
    St10,Et10,It10,Rt10 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_10, CMCO, epiparam, 44, 91, 36, nsims)

    St,Et,It,Rt = EpiSim_Msia.episim_nochangepoint(BA, epiparam, ndays, nsims)

    BA_vals = Dict("BA_St" => St, "BA_Et" => Et, "BA_It" => It, "BA_Rt" => Rt)
    results_BA_BABASELINE_1M[(β)] = BA_vals

    St_vals = Dict("St_95" => St95, "St_90" => St90, "St_85" => St85, "St_80" => St80, "St_70" => St70, "St_60" => St60,
    "St_50" => St50, "St_40" => St40, "St_30" => St30, "St_20" => St20, "St_10" => St10)

    Et_vals = Dict("Et_95" => Et95, "Et_90" => Et90, "Et_85" => Et85, "Et_80" => Et80, "Et_70" => Et70, "Et_60" => Et60,
    "Et_50" => Et50, "Et_40" => Et40, "Et_30" => Et30, "Et_20" => Et20, "Et_10" => Et10)

    It_vals = Dict("It_95" => It95, "It_90" => It90, "It_85" => It85, "It_80" => It80, "It_70" => It70, "It_60" => It60,
    "It_50" => It50, "It_40" => It40, "It_30" => It30, "It_20" => It20, "It_10" => It10)

    Rt_vals = Dict("Rt_95" => Rt95, "Rt_90" => Rt90, "Rt_85" => Rt85, "Rt_80" => Rt80, "Rt_70" => Rt70, "Rt_60" => Rt60,
    "Rt_50" => Rt50, "Rt_40" => Rt40, "Rt_30" => Rt30, "Rt_20" => Rt20, "Rt_10" => Rt10)

    for compliance in [95, 90, 85, 80, 70, 60, 50, 40, 30, 20, 10]
        results_St_BABASELINE_1M[(β, compliance)] = St_vals["St_$(compliance)"]
        results_Et_BABASELINE_1M[(β, compliance)] = Et_vals["Et_$(compliance)"]
        results_It_BABASELINE_1M[(β, compliance)] = It_vals["It_$(compliance)"]
        results_Rt_BABASELINE_1M[(β, compliance)] = Rt_vals["Rt_$(compliance)"]
    end

    # plot on a normal scale
    plot(title = "BA Baseline Model - 1747684 Nodes", 
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

    #EpiSim_Msia.plotqnt_normal(1:127, It, :red, "BA", 0.0) # no control measures at all (not included due to large magnitude)

    RealData_NormalScale_BABASELINE = plot!(1:127, new_cases,
        xlabel = "Days", ylabel = "Number of Active Cases", 
        label = "Real Data", linecolor = :brown, linewidth = 1.5,
        xtickfontsize = 8, grid = false)

    savefig(RealData_NormalScale_BABASELINE, "BA Baseline Model - 1747684 Nodes (Normal).png") # save plot

    # plot on a log scale
    plot(title = "BA Baseline Model - 1747684 Nodes", 
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

    RealData_LogScale_BABASELINE = plot!(1:127, new_cases,
        xlabel = "Days", ylabel = "Number of Active Cases", 
        label = "Real Data", linecolor = :brown, linewidth = 1.5,
        xtickfontsize = 8, grid = false, yaxis = :log)

    savefig(RealData_LogScale_BABASELINE, "BA Baseline Model - 1747684 Nodes (Log).png")
end

# Define file names to save each dictionary
file_St_BABASELINE_1M = "results_St_BABASELINE_1M.jls"
file_Et_BABASELINE_1M = "results_Et_BABASELINE_1M.jls"
file_It_BABASELINE_1M = "results_It_BABASELINE_1M.jls"
file_Rt_BABASELINE_1M = "results_Rt_BABASELINE_1M.jls"
file_BA_BABASELINE_1M = "results_BA_BABASELINE_1M.jls"

# Save the dictionaries to files (Can reload results for future use)
serialize(file_St_BABASELINE_1M, results_St_BABASELINE_1M)
serialize(file_Et_BABASELINE_1M, results_Et_BABASELINE_1M)
serialize(file_It_BABASELINE_1M, results_It_BABASELINE_1M)
serialize(file_Rt_BABASELINE_1M, results_Rt_BABASELINE_1M)
serialize(file_BA_BABASELINE_1M, results_BA_BABASELINE_1M)

# Load the dictionaries from files
results_St_BABASELINE_1M = deserialize("results_St_BABASELINE_1M.jls")
results_Et_BABASELINE_1M = deserialize("results_Et_BABASELINE_1M.jls")
results_It_BABASELINE_1M = deserialize("results_It_BABASELINE_1M.jls")
results_Rt_BABASELINE_1M = deserialize("results_Rt_BABASELINE_1M.jls")
results_BA_BABASELINE_1M = deserialize("results_BA_BABASELINE_1M.jls")

#___________________________________________________________________________________________________

#### Population 1024, 10000, 100489 (Change the numbers for different population sizes)
# The populations have to be square numbers for lattice model convenience

# epiparam["pop"] = 1024
epiparam["pop"] = 10000
#epiparam["pop"] = 100489
nsims = 1000 # number of simulations 

# use BA then limit ten for CMCO (cap on mass gatherings)
BA = ContactModels_Msia.fullmixing(10000)
serialize("networks_BABASELINE_A_10k.dat", (BA))
CMCO = ContactModels_Msia.limitmix(10000, 10)
serialize("networks_BABASELINE_B_10k.dat", (CMCO))

BA_MCOCMCO_95 = ContactModels_Msia.socialdist(10000, 0.95)
serialize("networks_BABASELINE_C_10k.dat", (BA_MCOCMCO_95))
BA_MCOCMCO_90 = ContactModels_Msia.socialdist(10000, 0.9)
serialize("networks_BABASELINE_D_10k.dat", (BA_MCOCMCO_90))
BA_MCOCMCO_85 = ContactModels_Msia.socialdist(10000, 0.85)
serialize("networks_BABASELINE_E_10k.dat", (BA_MCOCMCO_85))
BA_MCOCMCO_80 = ContactModels_Msia.socialdist(10000, 0.8)
serialize("networks_BABASELINE_F_10k.dat", (BA_MCOCMCO_80))
BA_MCOCMCO_70 = ContactModels_Msia.socialdist(10000, 0.7)
serialize("networks_BABASELINE_G_10k.dat", (BA_MCOCMCO_70))
BA_MCOCMCO_60 = ContactModels_Msia.socialdist(10000, 0.6)
serialize("networks_BABASELINE_H_10k.dat", (BA_MCOCMCO_60))
BA_MCOCMCO_50 = ContactModels_Msia.socialdist(10000, 0.5)
serialize("networks_BABASELINE_I_10k.dat", (BA_MCOCMCO_50))
BA_MCOCMCO_40 = ContactModels_Msia.socialdist(10000, 0.4)
serialize("networks_BABASELINE_J_10k.dat", (BA_MCOCMCO_40))
BA_MCOCMCO_30 = ContactModels_Msia.socialdist(10000, 0.3)
serialize("networks_BABASELINE_K_10k.dat", (BA_MCOCMCO_30))
BA_MCOCMCO_20 = ContactModels_Msia.socialdist(10000, 0.2)
serialize("networks_BABASELINE_L_10k.dat", (BA_MCOCMCO_20))
BA_MCOCMCO_10 = ContactModels_Msia.socialdist(10000, 0.1)
serialize("networks_BABASELINE_M_10k.dat", (BA_MCOCMCO_10))


# Save to a file
serialize("networks_BABASELINE_10k.dat", (BA, CMCO, BA_MCOCMCO_95, BA_MCOCMCO_90, BA_MCOCMCO_85, 
BA_MCOCMCO_80, BA_MCOCMCO_70, BA_MCOCMCO_60, BA_MCOCMCO_50, BA_MCOCMCO_40, BA_MCOCMCO_30, 
BA_MCOCMCO_20, BA_MCOCMCO_10))

# Load from the file
(BA, CMCO, BA_MCOCMCO_95, BA_MCOCMCO_90, BA_MCOCMCO_85, BA_MCOCMCO_80, BA_MCOCMCO_70, 
BA_MCOCMCO_60, BA_MCOCMCO_50, BA_MCOCMCO_40, BA_MCOCMCO_30, 
BA_MCOCMCO_20, BA_MCOCMCO_10) = deserialize("networks_BABASELINE_10k.dat")

# -------------------------------------------------------------------------------------

#### Plot the lineplots for BA Baseline with 1024/10000/100489 population

results_St_BABASELINE_10k = Dict()
results_Et_BABASELINE_10k = Dict()
results_It_BABASELINE_10k = Dict()
results_Rt_BABASELINE_10k = Dict()
results_BA_BABASELINE_10k = Dict() 

for β in transmission_rates
    epiparam["β_fixed"] = β

    St95,Et95,It95,Rt95 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_95, CMCO, epiparam, 44, 91, 36, nsims)
    St90,Et90,It90,Rt90 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_90, CMCO, epiparam, 44, 91, 36, nsims)
    St85,Et85,It85,Rt85 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_85, CMCO, epiparam, 44, 91, 36, nsims)
    St80,Et80,It80,Rt80 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_80, CMCO, epiparam, 44, 91, 36, nsims)
    St70,Et70,It70,Rt70 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_70, CMCO, epiparam, 44, 91, 36, nsims)
    St60,Et60,It60,Rt60 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_60, CMCO, epiparam, 44, 91, 36, nsims)
    St50,Et50,It50,Rt50 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_50, CMCO, epiparam, 44, 91, 36, nsims)
    St40,Et40,It40,Rt40 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_40, CMCO, epiparam, 44, 91, 36, nsims)
    St30,Et30,It30,Rt30 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_30, CMCO, epiparam, 44, 91, 36, nsims)
    St20,Et20,It20,Rt20 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_20, CMCO, epiparam, 44, 91, 36, nsims)
    St10,Et10,It10,Rt10 = EpiSim_Msia.episim_twocpt_samerate_MCOCMCO(BA, BA_MCOCMCO_10, CMCO, epiparam, 44, 91, 36, nsims)

    St,Et,It,Rt = EpiSim_Msia.episim_nochangepoint(BA, epiparam, ndays, nsims)

    BA_vals = Dict("BA_St" => St, "BA_Et" => Et, "BA_It" => It, "BA_Rt" => Rt)
    results_BA_BABASELINE_10k[(β)] = BA_vals

    St_vals = Dict("St_95" => St95, "St_90" => St90, "St_85" => St85, "St_80" => St80, "St_70" => St70, "St_60" => St60,
    "St_50" => St50, "St_40" => St40, "St_30" => St30, "St_20" => St20, "St_10" => St10)

    Et_vals = Dict("Et_95" => Et95, "Et_90" => Et90, "Et_85" => Et85, "Et_80" => Et80, "Et_70" => Et70, "Et_60" => Et60,
    "Et_50" => Et50, "Et_40" => Et40, "Et_30" => Et30, "Et_20" => Et20, "Et_10" => Et10)

    It_vals = Dict("It_95" => It95, "It_90" => It90, "It_85" => It85, "It_80" => It80, "It_70" => It70, "It_60" => It60,
    "It_50" => It50, "It_40" => It40, "It_30" => It30, "It_20" => It20, "It_10" => It10)

    Rt_vals = Dict("Rt_95" => Rt95, "Rt_90" => Rt90, "Rt_85" => Rt85, "Rt_80" => Rt80, "Rt_70" => Rt70, "Rt_60" => Rt60,
    "Rt_50" => Rt50, "Rt_40" => Rt40, "Rt_30" => Rt30, "Rt_20" => Rt20, "Rt_10" => Rt10)

    for compliance in [95, 90, 85, 80, 70, 60, 50, 40, 30, 20, 10]
        results_St_BABASELINE_10k[(β, compliance)] = St_vals["St_$(compliance)"]
        results_Et_BABASELINE_10k[(β, compliance)] = Et_vals["Et_$(compliance)"]
        results_It_BABASELINE_10k[(β, compliance)] = It_vals["It_$(compliance)"]
        results_Rt_BABASELINE_10k[(β, compliance)] = Rt_vals["Rt_$(compliance)"]
    end

    # plot on a normal scale
    plot(title = "BA Baseline Model - 10000 Nodes", 
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

    RealData_NormalScale_BABASELINE = plot!(1:127, new_cases,
        xlabel = "Days", ylabel = "Number of Active Cases", 
        label = "Real Data", linecolor = :brown, linewidth = 1.5,
        xtickfontsize = 8, grid = false)

    savefig(RealData_NormalScale_BABASELINE, "BA Baseline Model - 10000 Nodes (Normal).png")

    # plot on a log scale
    plot(title = "BA Baseline Model - 10000 Nodes", 
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

    # EpiSim_Msia.plotqnt_log(1:127, It, :red, "BA", 0.0)

    RealData_LogScale_BABASELINE = plot!(1:127, new_cases,
        xlabel = "Days", ylabel = "Number of Active Cases", 
        label = "Real Data", linecolor = :brown, linewidth = 1.5,
        xtickfontsize = 8, grid = false, yaxis = :log)

    savefig(RealData_LogScale_BABASELINE, "BA Baseline Model - 10000 Nodes (Log).png")
end

# Define file names to save each dictionary
file_St_BABASELINE_10k = "results_St_BABASELINE_10k.jls"
file_Et_BABASELINE_10k = "results_Et_BABASELINE_10k.jls"
file_It_BABASELINE_10k = "results_It_BABASELINE_10k.jls"
file_Rt_BABASELINE_10k = "results_Rt_BABASELINE_10k.jls"
file_BA_BABASELINE_10k = "results_BA_BABASELINE_10k.jls"

# Save the dictionaries to files (Can reload results for future use)
serialize(file_St_BABASELINE_10k, results_St_BABASELINE_10k)
serialize(file_Et_BABASELINE_10k, results_Et_BABASELINE_10k)
serialize(file_It_BABASELINE_10k, results_It_BABASELINE_10k)
serialize(file_Rt_BABASELINE_10k, results_Rt_BABASELINE_10k)
serialize(file_BA_BABASELINE_10k, results_BA_BABASELINE_10k)

#___________________________________________________________________________________________________

# Load the dictionaries from files
results_St_BABASELINE_10k = deserialize("results_St_BABASELINE_10k.jls")
results_Et_BABASELINE_10k = deserialize("results_Et_BABASELINE_10k.jls")
results_It_BABASELINE_10k = deserialize("results_It_BABASELINE_10k.jls")
results_Rt_BABASELINE_10k = deserialize("results_Rt_BABASELINE_10k.jls")
results_BA_BABASELINE_10k = deserialize("results_BA_BABASELINE_10k.jls")

#_______________________________________________________________________________________________

#=

# Network Modification via Scaled SBM Rewiring

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
#num_nodes = 1747684
proportions = [0.694, 0.306] # Majority (bumiputera) and Minority (Chinese, Indian, Others) proportions
ethnic_group = assign_ethnic_group_to_nodes(num_nodes, proportions)

# check the proportions of each ethnicity
unique_ethinicity = unique(ethnic_group)
counts = Dict(ethnic => count(x -> x == ethnic, ethnic_group) for ethnic in unique_ethinicity)

ethnic_group_dict = Dict(i => ethnic_group[i] for i in 1:length(ethnic_group))

# the probabilities of edge formation between and within ethnic groups
ethnic_prs = [
    0.071 0.012;
    0.012 0.156
]

# number of nodes in the Majority and Minority groups respectively
number_per_gp = [711, 313] # change proportions according to number of nodes
#number_per_gp = [6940, 3060]
#number_per_gp = [69739, 30750]
#number_per_gp = [1212893, 534791]

using SparseArrays

# implement sbm on the orginal networks, set the ethnicity group of each node, 
# and count the expected edges between and within groups
# if there happens to be a surplus or deficit of edges, we will swap edges 
# between groups to bring the network closer to the expected edges
function modify_net_using_sbm(net::SimpleGraph, ethnic_group::Vector{Any}, ethnic_prs::Matrix{Float64}, 
    number_per_gp::Vector{Int}, rewired_edges::Vector{Tuple{Int, Int}})
    
    net = deepcopy(net) # make a copy of the original network to avoid modifying it directly
    
    n = nv(net) # obtain the number of nodes in the original network
    m = ne(net) # obtain the number of edges in the original network

    rewired_nodes = unique(vcat([e[1] for e in rewired_edges], [e[2] for e in rewired_edges]))

    # compute the expected number of edges for each ethnicity group
    # the number of edges also needs to be scaled to adjust to the number of edges in the original network

    expected_edges_before_scaling = 0

    expected_counts = Dict{Tuple{Int, Int}, Int}()
    for i in 1:2, j in i:2 # obtain expected total number of edges within groups 1,1 and 2,2 and between groups 1,2
        if i == j
            expected_edges_before_scaling += ethnic_prs[i,j] * number_per_gp[i] * (number_per_gp[i] - 1) /2
        else
            expected_edges_before_scaling += ethnic_prs[i,j] * number_per_gp[i] * number_per_gp[j]
        end
    end 

    println("Expected total edges before scaling: ", round(expected_edges_before_scaling, digits = 3)) # three decimal places

    # scaling factor is to make sure the number of edges in the modified network is the same as the original network
    # this is important because the number of edges in the original network is not equal to the expected number of edges
    scaling_factor = m / expected_edges_before_scaling
    println("Scaling factor: ", scaling_factor)

    # scale the ethnicity group probabilities
    ethnic_prs_scaled = copy(ethnic_prs)
    for i in 1:2, j in 1:2
        ethnic_prs_scaled[i, j] *= scaling_factor
    end

    # obtain the number of expected edges for within and between ethnicity groups
    expected_edges = Dict{Tuple{Int, Int}, Int}()
    for i in 1:2, j in i:2
        if i == j
            expected = round(Int, ethnic_prs_scaled[i, j] * number_per_gp[i] * (number_per_gp[i] - 1) / 2)
        else
            expected = round(Int, ethnic_prs_scaled[i, j] * number_per_gp[i] * number_per_gp[j])
        end
        expected_edges[(i, j)] = expected
    end

    println("\nExpected edges after scaling: ")
    for ((i, j), count) in expected_edges
        println("Group ($i, $j): $count")
    end

    # count the actual number of edges in the original network
    actual_edges = Dict{Tuple{Int, Int}, Int}()
    for ((i,j), _) in expected_edges
        actual_edges[(i,j)] = 0
    end
    
    # loop through all the edges in the original network and count the number of edges for each ethnicity group
    for edge in edges(net)
        i, j = src(edge), dst(edge)
        p, q = ethnic_group[i], ethnic_group[j]
        key = p <= q ? (p, q) : (q, p) # here 1,2 and 2,1 are the same, regarded as 1,2 (no double counting)
        actual_edges[key] += 1
    end

    println("\nCurrent edges in the original network: ")
    for ((i, j), count) in actual_edges
        println("Group ($i, $j): $count")
    end

    # now we look at the deficits and surpluses of edges in the original network for each grouping
    # edges from groups with surplus edges will be swapped to groups with deficit edges (only one end is swapped)
    # this swapping brings the original network closer to the expected edges within and between groups, 
    # according to the SBM structure
    surpluses = Dict{Tuple{Int, Int}, Int}()
    deficits = Dict{Tuple{Int, Int}, Int}()    

    # determine if groups 1,1 / 2,2 / 1,2 are surplus or deficit, and how many edges need to be added or removed
    for (key, expected) in expected_edges
        actual = actual_edges[key]
        if actual > expected
            surpluses[key] = actual - expected
        elseif actual < expected
            deficits[key] = expected - actual
        end
    end

    println("\nSurplus edges (the amount that needs to be removed by swapping): ")
    for (key, count) in surpluses
        println("Group $key: remove $count edges")
    end

    println("\nDeficit edges (the amount that needs to be added by swapping): ")
    for (key, count) in deficits
        println("Group $key: add $count edges")
    end

    # think of swapping like rewiring current edges 
    iterations = 0 # how many times through the rewiring loop? ensures the number of attempts does not exceed the max attempts
    successful_swaps = 0 # how many times a swap was made (there could be 100 attempts but 50 successful swaps)
    max_attempts = 10000 # maximum number of attempts to swap edge

    # make sure there are still surpluses and deficits to swap, and it does not go over 
    # the max attempts set to avoid infinite loop
    while !isempty(surpluses) && !isempty(deficits) && iterations < max_attempts

        # randomly select a surplus edge 
        surplus_keys = collect(keys(surpluses))
        surplus_key = rand(surplus_keys)

        # find all the surplus edges from this group pair
        surplus_edges = [e for e in edges(net) if begin
            i, j = src(e), dst(e)
            p, q = ethnic_group[i], ethnic_group[j]
            k = p <= q ? (p, q) : (q, p) # here 1,2 and 2,1 are the same, regarded as 1,2 (no double counting)
            k == surplus_key
        end]

        if isempty(surplus_edges)
            delete!(surpluses, surplus_key) # remove the group pair from surpluses if there are no edges in this group
            continue
        end

        e = rand(surplus_edges) # randomly select a surplus edge from the surplus edges
        i, j = src(e), dst(e) # get the nodes of the surplus edge
        p, q = ethnic_group[i], ethnic_group[j] # get the ethnic groups of the nodes of the surplus edge

        # allow for flexibility by allowing any end of the edge to be rewired to a deficit group
        ends = [(i, p, j, q), (j, q, i, p)] # (fixed node, fixed ethnicity, rewiring node, rewiring ethnicity)
        rewired = false # flag to indicate if rewiring was successful

        for (fixed_node, fixed_ethnicity, rewiring_node, rewiring_ethnicity) in ends
            for (target_key, _) in deficits

                # check if the fixed node is not involved in the current deficit pair, making it irrelevant and can be skipped
                if fixed_ethnicity != target_key[1] && fixed_ethnicity != target_key[2]
                    continue # skip this iteration if the fixed node is not in the current deficit pair
                end

                if target_key == (fixed_ethnicity, fixed_ethnicity) # so a deficit happens within group
                    candidate_nodes = [node for node in rewired_nodes if
                    ethnic_group[node] == fixed_ethnicity && 
                    node != fixed_node && # avoid self-loops
                    !has_edge(net, fixed_node, node) # avoid existing edges to avoid multiple edges
                    ]
                    print("Candidate nodes for rewiring within group $target_key: ", candidate_nodes, "\n")
                else
                    # a between group edge is deficit
                    other_ethnicity = target_key[1] == fixed_ethnicity ? target_key[2] : target_key[1]
                    candidate_nodes = [node for node in rewired_nodes if
                    ethnic_group[node] == other_ethnicity &&
                    node != fixed_node && # avoid self-loops
                    !has_edge(net, fixed_node, node) # avoid existing edges to avoid multiple edges
                    ]
                    print("Candidate nodes for rewiring between groups $target_key: ", candidate_nodes, "\n")
                end

                if isempty(candidate_nodes)
                    continue # no candidate nodes available for rewiring, so skip 
                end

                new_rewiring_partner = rand(candidate_nodes) # randomly select a rewiring partner from the candidate nodes
                new_ethnicities = (ethnic_group[fixed_node], ethnic_group[new_rewiring_partner]) 
                new_key = new_ethnicities[1] <= new_ethnicities[2] ? (new_ethnicities[1], new_ethnicities[2]) : (new_ethnicities[2], new_ethnicities[1])

                # do the rewiring
                rem_edge!(net, i, j)
                add_edge!(net, fixed_node, new_rewiring_partner) 

                actual_edges[surplus_key] -= 1
                actual_edges[new_key] += 1

                surpluses[surplus_key] -= 1
                if surpluses[surplus_key] == 0
                    delete!(surpluses, surplus_key) # remove the group pair from surpluses if there are no edges left in this group  
                end

                deficits[target_key] -= 1
                if deficits[target_key] == 0
                    delete!(deficits, target_key)
                end

                # println("Rewired edge ($i, $j) to ($fixed_node, $new_rewiring_partner), moved from group $surplus_key to group $new_key")

                successful_swaps += 1 # increment the total number of successful swaps made
                rewired = true # set the flag to true to indicate rewiring was successful
                break # stop looking for deficits for this edge in ends
            end

            if rewired
                break # stop looping through ends since reiwiring was already done for one end
            end
        end

        iterations += 1 # increment the total number of attempts made
    end

    println("\nTotal successful swaps done: $successful_swaps")
    if !isempty(deficits)
        println("Unfortunately, not all deficits were satisfied. Max Attempts were reached.")
    else
        println("All deficits were satisfied.")
    end

    return net # return the modified network 
end
                    
# modify the networks from earlier according to the SBM structure SBM
modified_SBM_BA = modify_net_using_sbm(BA, ethnic_group, ethnic_prs, number_per_gp)
modified_SBM_CMCO = modify_net_using_sbm(CMCO, ethnic_group, ethnic_prs, number_per_gp)
modified_SBM_BA_MCOCMCO_95 = modify_net_using_sbm(BA_MCOCMCO_95, ethnic_group, ethnic_prs, number_per_gp, rewired_edges_95)
modified_SBM_BA_MCOCMCO_90 = modify_net_using_sbm(BA_MCOCMCO_90, ethnic_group, ethnic_prs, number_per_gp)
modified_SBM_BA_MCOCMCO_85 = modify_net_using_sbm(BA_MCOCMCO_85, ethnic_group, ethnic_prs, number_per_gp)
modified_SBM_BA_MCOCMCO_80 = modify_net_using_sbm(BA_MCOCMCO_80, ethnic_group, ethnic_prs, number_per_gp)
modified_SBM_BA_MCOCMCO_70 = modify_net_using_sbm(BA_MCOCMCO_70, ethnic_group, ethnic_prs, number_per_gp)
modified_SBM_BA_MCOCMCO_60 = modify_net_using_sbm(BA_MCOCMCO_60, ethnic_group, ethnic_prs, number_per_gp)
modified_SBM_BA_MCOCMCO_50 = modify_net_using_sbm(BA_MCOCMCO_50, ethnic_group, ethnic_prs, number_per_gp)
modified_SBM_BA_MCOCMCO_40 = modify_net_using_sbm(BA_MCOCMCO_40, ethnic_group, ethnic_prs, number_per_gp)
modified_SBM_BA_MCOCMCO_30 = modify_net_using_sbm(BA_MCOCMCO_30, ethnic_group, ethnic_prs, number_per_gp)
modified_SBM_BA_MCOCMCO_20 = modify_net_using_sbm(BA_MCOCMCO_20, ethnic_group, ethnic_prs, number_per_gp)
modified_SBM_BA_MCOCMCO_10 = modify_net_using_sbm(BA_MCOCMCO_10, ethnic_group, ethnic_prs, number_per_gp)


serialize("modified_SBM_networks_1k.dat", (modified_SBM_BA, modified_SBM_CMCO, modified_SBM_BA_MCOCMCO_95, 
modified_SBM_BA_MCOCMCO_90, modified_SBM_BA_MCOCMCO_85, modified_SBM_BA_MCOCMCO_80, modified_SBM_BA_MCOCMCO_70, 
modified_SBM_BA_MCOCMCO_60, modified_SBM_BA_MCOCMCO_50, modified_SBM_BA_MCOCMCO_40, modified_SBM_BA_MCOCMCO_30, 
modified_SBM_BA_MCOCMCO_20, modified_SBM_BA_MCOCMCO_10))


(modified_SBM_BA, modified_SBM_CMCO, modified_SBM_BA_MCOCMCO_95, 
modified_SBM_BA_MCOCMCO_90, modified_SBM_BA_MCOCMCO_85, modified_SBM_BA_MCOCMCO_80, modified_SBM_BA_MCOCMCO_70, 
modified_SBM_BA_MCOCMCO_60, modified_SBM_BA_MCOCMCO_50, modified_SBM_BA_MCOCMCO_40, modified_SBM_BA_MCOCMCO_30, 
modified_SBM_BA_MCOCMCO_20, modified_SBM_BA_MCOCMCO_10) = deserialize("modified_SBM_networks_1k.dat")


results_St_MODIFIEDSBM_1k = Dict()
results_Et_MODIFIEDSBM_1k = Dict()
results_It_MODIFIEDSBM_1k = Dict()
results_Rt_MODIFIEDSBM_1k = Dict()
results_BA_MODIFIEDSBM_1k = Dict()

nsims = 1000

for β in transmission_rates
    epiparam["β_intra"] = β
    epiparam["β_inter"] = β/2

    St95,Et95,It95,Rt95 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(modified_SBM_BA, modified_SBM_BA_MCOCMCO_95, modified_SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St90,Et90,It90,Rt90 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(modified_SBM_BA, modified_SBM_BA_MCOCMCO_90, modified_SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St85,Et85,It85,Rt85 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(modified_SBM_BA, modified_SBM_BA_MCOCMCO_85, modified_SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St80,Et80,It80,Rt80 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(modified_SBM_BA, modified_SBM_BA_MCOCMCO_80, modified_SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St70,Et70,It70,Rt70 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(modified_SBM_BA, modified_SBM_BA_MCOCMCO_70, modified_SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St60,Et60,It60,Rt60 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(modified_SBM_BA, modified_SBM_BA_MCOCMCO_60, modified_SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St50,Et50,It50,Rt50 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(modified_SBM_BA, modified_SBM_BA_MCOCMCO_50, modified_SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St40,Et40,It40,Rt40 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(modified_SBM_BA, modified_SBM_BA_MCOCMCO_40, modified_SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St30,Et30,It30,Rt30 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(modified_SBM_BA, modified_SBM_BA_MCOCMCO_30, modified_SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St20,Et20,It20,Rt20 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(modified_SBM_BA, modified_SBM_BA_MCOCMCO_20, modified_SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)
    St10,Et10,It10,Rt10 = EpiSim_Msia.episim_homophily_twocpt_samerate_MCOCMCO(modified_SBM_BA, modified_SBM_BA_MCOCMCO_10, modified_SBM_CMCO, epiparam, ethnic_group, 44, 91, 36, nsims)

    St,Et,It,Rt = EpiSim_Msia.episim_homophily_nochangepoint(modified_SBM_BA, epiparam, ethnic_group, ndays, nsims)

    BA_vals = Dict("BA_St" => St, "BA_Et" => Et, "BA_It" => It, "BA_Rt" => Rt)
    results_BA_MODIFIEDSBM_1k[(β)] = BA_vals

    St_vals = Dict("St_95" => St95, "St_90" => St90, "St_85" => St85, "St_80" => St80, "St_70" => St70, "St_60" => St60,
    "St_50" => St50, "St_40" => St40, "St_30" => St30, "St_20" => St20, "St_10" => St10)

    Et_vals = Dict("Et_95" => Et95, "Et_90" => Et90, "Et_85" => Et85, "Et_80" => Et80, "Et_70" => Et70, "Et_60" => Et60,
    "Et_50" => Et50, "Et_40" => Et40, "Et_30" => Et30, "Et_20" => Et20, "Et_10" => Et10)

    It_vals = Dict("It_95" => It95, "It_90" => It90, "It_85" => It85, "It_80" => It80, "It_70" => It70, "It_60" => It60,
    "It_50" => It50, "It_40" => It40, "It_30" => It30, "It_20" => It20, "It_10" => It10)

    Rt_vals = Dict("Rt_95" => Rt95, "Rt_90" => Rt90, "Rt_85" => Rt85, "Rt_80" => Rt80, "Rt_70" => Rt70, "Rt_60" => Rt60,
    "Rt_50" => Rt50, "Rt_40" => Rt40, "Rt_30" => Rt30, "Rt_20" => Rt20, "Rt_10" => Rt10)

    for compliance in [95, 90, 85, 80, 70, 60, 50, 40, 30, 20, 10]
        results_St_MODIFIEDSBM_1k[(β, compliance)] = St_vals["St_$(compliance)"]
        results_Et_MODIFIEDSBM_1k[(β, compliance)] = Et_vals["Et_$(compliance)"]
        results_It_MODIFIEDSBM_1k[(β, compliance)] = It_vals["It_$(compliance)"]
        results_Rt_MODIFIEDSBM_1k[(β, compliance)] = Rt_vals["Rt_$(compliance)"]
    end

    Plots.plot(title = "Modified Networks via SBM - 1024 Nodes", 
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

    RealData_NormalScale_MODIFIEDSBM = Plots.plot!(1:127, new_cases,
        xlabel = "Days", ylabel = "Number of Active Cases", 
        label = "Real Data", linecolor = :brown, linewidth = 1.5,
        xtickfontsize = 8, grid = false)

    savefig(RealData_NormalScale_MODIFIEDSBM, "Modified Networks via SBM - 1024 Nodes (Normal).png")

    Plots.plot(title = "Modified Networks via SBM - 1024 Nodes", 
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

    RealData_LogScale_MODIFIEDSBM = Plots.plot!(1:127, new_cases,
        xlabel = "Days", ylabel = "Number of Active Cases", 
        label = "Real Data", linecolor = :brown, linewidth = 1.5,
        xtickfontsize = 8, grid = false, yaxis = :log)

    savefig(RealData_LogScale_MODIFIEDSBM, "Modified Networks via SBM - 1024 Nodes (Log).png")
end

# Define file names to save each dictionary
file_St_MODIFIEDSBM_1k = "results_St_MODIFIEDSBM_1k.jls"
file_Et_MODIFIEDSBM_1k = "results_Et_MODIFIEDSBM_1k.jls"
file_It_MODIFIEDSBM_1k = "results_It_MODIFIEDSBM_1k.jls"
file_Rt_MODIFIEDSBM_1k = "results_Rt_MODIFIEDSBM_1k.jls"
file_BA_MODIFIEDSBM_1k = "results_BA_MODIFIEDSBM_1k.jls"

# Save the dictionaries to files (Can reload results for future use)
serialize(file_St_MODIFIEDSBM_1k, results_St_MODIFIEDSBM_1k)
serialize(file_Et_MODIFIEDSBM_1k, results_Et_MODIFIEDSBM_1k)
serialize(file_It_MODIFIEDSBM_1k, results_It_MODIFIEDSBM_1k)
serialize(file_Rt_MODIFIEDSBM_1k, results_Rt_MODIFIEDSBM_1k)
serialize(file_BA_MODIFIEDSBM_1k, results_BA_MODIFIEDSBM_1k)

#___________________________________________________________________________________________________

# Load the dictionaries from files
results_St_MODIFIEDSBM_1k = deserialize("results_St_MODIFIEDSBM_1k.jls")
results_Et_MODIFIEDSBM_1k = deserialize("results_Et_MODIFIEDSBM_1k.jls")
results_It_MODIFIEDSBM_1k = deserialize("results_It_MODIFIEDSBM_1k.jls")
results_Rt_MODIFIEDSBM_1k = deserialize("results_Rt_MODIFIEDSBM_1k.jls")
results_BA_MODIFIEDSBM_1k = deserialize("results_BA_MODIFIEDSBM_1k.jls")

#___________________________________________________________________________________________________

=#

### Plot networks to view
 
using GraphMakie, CairoMakie

G = BA

node_colours = [:blue for _ in 1:nv(G)] # set the node colours

# Plot the graph
fig, ax, p = GraphMakie.graphplot(G, node_size = 6, node_color = node_colours)

# Save the figure as png file
save("Plot BA Baseline Model - 1024 Nodes (BA).png", fig) 

# Display the plot
fig

#________________________________________________________________

# Plot the graph with colors based on ethnic group (Modified SBM)

# set the colours for the Majority and Minority ethnic groups
ethnic_colors = Dict(1 => :dodgerblue, 2 => :orange)  

# Create the node_colors vector based on ethnic group
node_colours = [ethnic_colors[ethnic_group[i]] for i in 1:nv(G)]

# Plot the graph
fig, ax, p = GraphMakie.graphplot(G, node_size = 6, node_color = node_colours)

# Save and display the plot
save("Plot Modified SBM Model - 1024 Nodes.png", fig)

fig

#________________________________________________________________

### Done for BA Baseline Model and Modified SBM Model




