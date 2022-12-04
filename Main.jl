cd(dirname(@__FILE__))
using Pkg
Pkg.activate("")
using PowerModels
using Ipopt
using JuMP
using JSON
using Gurobi

logger_config!("debug")

json_str = readline(stdin)
ding0_grid = ARGS[1]
results_path = ARGS[2]
method = ARGS[3]
silence_moi = ARGS[4].=="True"

function optimize_edisgo()                                        
  # read in data and create multinetwork
  data_edisgo = parse_json(json_str)
  data_edisgo_mn = make_multinetwork(data_edisgo)

  if method == "soc" # Second order cone
    # Set solver attributes
    gurobi = optimizer_with_attributes(Gurobi.Optimizer, MOI.Silent() => silence_moi, "Presolve" => 1,  "QCPDual" =>1, "BarQCPConvTol" => 1e-4, "BarConvTol" => 1e-6, "BarHomogeneous"=> 1)#,"FeasibilityTol"=>1e-4, "NumericFocus"=> 1)
    # Solve SOC model
    result_soc = solve_mn_opf_bf_flex(data_edisgo_mn, SOCBFPowerModelEdisgo, gurobi)
    # Check if SOC constraint is tight
    exactness = check_SOC_equality(result, data_edisgo)
    open(joinpath(results_path, ding0_grid*"_SOC_tightness.json"), "w") do f
        write(f, JSON.json(exactness))
    end
    update_data!(data_edisgo_mn, result_soc["solution"])
    set_ac_bf_start_values!(data_edisgo_mn["nw"]["1"])
    result_nc_ws = solve_mn_opf_bf_flex(data_edisgo_mn, NCBFPowerModelEdisgo, ipopt)
  elseif method == "nc" # Non-Convex
    # Set solver attributes
    ipopt = optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => silence_moi, "sb" => "yes")#, "tol"=>1e-4)
    # Solve NC model
    result = solve_mn_opf_bf_flex(data_edisgo_mn, NCBFPowerModelEdisgo, ipopt)
    update_data!(data_edisgo_mn, result["solution"])
  end

  # Update network data with optimization results and print to stdout
  print(JSON.json(data_edisgo_mn))
end

if abspath(PROGRAM_FILE) == @__FILE__
  optimize_edisgo()
end