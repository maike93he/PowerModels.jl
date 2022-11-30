cd(dirname(@__FILE__))
using Pkg
Pkg.activate("")
using PowerModels
using Ipopt
using JuMP
using JSON
using Gurobi

const ipopt = optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => true, "sb" => "yes",
                                        "tol"=>1e-4)
const gurobi = optimizer_with_attributes(Gurobi.Optimizer, MOI.Silent() => false, "Presolve" => 1,
                                        "FeasibilityTol"=>1e-4, "QCPDual" =>1,
                                         "BarQCPConvTol" => 1e-12, "BarHomogeneous"=> 1)

json_str = readline(stdin)
ding0_grid = ARGS[1]
results_path = ARGS[2]
method = ARGS[3]

# read in data and create multinetwork
data_edisgo = parse_json(json_str)
data_edisgo_mn = make_multinetwork(data_edisgo)

# solve opf with soc or nc method
if method == "soc"
  # Second order cone
  result = solve_mn_opf_bf_flex(data_edisgo_mn, SOCBFPowerModelEdisgo, gurobi)
elseif method == "nc"
  # Non-Convex
  result = solve_mn_opf_bf_flex(data_edisgo_mn, NCBFPowerModelEdisgo, ipopt)
end

exactness = check_SOC_equality(result, data_edisgo)
stringdata = JSON.json(exactness)
print(stringdata)
println(" ")
update_data!(data_edisgo_mn, result["solution"])

json_string = JSON.json(data_edisgo_mn)
print(json_string)