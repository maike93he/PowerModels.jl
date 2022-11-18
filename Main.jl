cd(dirname(@__FILE__))
using Pkg
Pkg.activate("")
import Pkg; Pkg.add("Gurobi")
using PowerModels
using Ipopt
using JuMP
using JSON
using Gurobi

const ipopt = optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => false, "sb" => "yes", "tol"=>1e-6)
const gurobi = optimizer_with_attributes(Gurobi.Optimizer, MOI.Silent() => false, "Presolve" => 1, "FeasibilityTol"=>1e-6, "QCPDual" =>1, "BarQCPConvTol" => 1e-12)#, "BarHomogeneous"=> 1)

grid_path = ARGS[1]
ding0_grid = ARGS[2]
results_path = ARGS[3]
method = ARGS[4]
# read in data and create multinetwork
data_edisgo = parse_json(grid_path)
data_edisgo_mn = make_multinetwork(data_edisgo)
println("multinetwork created")

# solve opf with soc or nc method
if method == "soc"
  # Second order cone
  result = solve_mn_opf_bf_flex(data_edisgo_mn, SOCBFPowerModelEdisgo, gurobi)
elseif method == "pm"
  # Second order cone
  result = solve_mn_opf_bf_strg(data_edisgo_mn, SOCBFPowerModel, ipopt)
elseif method == "nc"
  # Non-Convex
  result = solve_mn_opf_flex(data_edisgo_mn, ACRPowerModel, ipopt)
end

test = check_SOC_equality(result, data_edisgo)
stringdata = JSON.json(test)
open(results_path*"/"*ding0_grid*"_SOC_eq.json", "w") do f
  write(f, stringdata)
end

# save result in json file
#json_string = JSON.json(result)
update_data!(data_edisgo_mn, result["solution"])
# save result in json file
json_string = JSON.json(data_edisgo_mn)

open(results_path*"/"*ding0_grid*".json","w") do f
  JSON.print(f, json_string)
end