""
function solve_opf_bf(file, model_type::Type{T}, optimizer; kwargs...) where T <: AbstractBFModel
    return solve_model(file, model_type, optimizer, build_opf_bf; kwargs...)
end

""
function solve_mn_opf_bf_strg(file, model_type::Type{T}, optimizer; kwargs...) where T <: AbstractBFModel
    return solve_model(file, model_type, optimizer, build_mn_opf_bf_strg; multinetwork=true, kwargs...)
end

""
function solve_mn_opf_bf_flex(file, model_type::Type{T}, optimizer; kwargs...) where T <: AbstractBFModel
    return solve_model(file, model_type, optimizer, build_mn_opf_bf_flex; multinetwork=true, kwargs...)
end

""
function build_opf_bf(pm::AbstractPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_branch_current(pm)
    variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost(pm)

    constraint_model_current(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_power_losses(pm, i)
        constraint_voltage_magnitude_difference(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end
end

"Build multinetwork branch flow storage OPF"
function build_mn_opf_bf_strg(pm::AbstractPowerModel)
    for (n, network) in nws(pm)
        variable_bus_voltage(pm, nw=n)
        variable_gen_power(pm, nw=n)
        variable_storage_power_mi(pm, nw=n)
        variable_branch_power(pm, nw=n)
        variable_branch_current(pm, nw=n)
        variable_dcline_power(pm, nw=n)

        constraint_model_current(pm, nw=n)

        for i in ids(pm, :ref_buses, nw=n)
            constraint_theta_ref(pm, i, nw=n)
        end

        for i in ids(pm, :bus, nw=n)
            constraint_power_balance(pm, i, nw=n)
        end

        for i in ids(pm, :storage, nw=n)
            constraint_storage_complementarity_mi(pm, i, nw=n)
            constraint_storage_losses(pm, i, nw=n)
            constraint_storage_thermal_limit(pm, i, nw=n)
        end

        for i in ids(pm, :branch, nw=n)
            constraint_power_losses(pm, i, nw=n)
            constraint_voltage_magnitude_difference(pm, i, nw=n)

            constraint_voltage_angle_difference(pm, i, nw=n)

            constraint_thermal_limit_from(pm, i, nw=n)
            constraint_thermal_limit_to(pm, i, nw=n)
        end

        for i in ids(pm, :dcline, nw=n)
            constraint_dcline_power_losses(pm, i, nw=n)
        end
    end

    network_ids = sort(collect(nw_ids(pm)))

    n_1 = network_ids[1]
    for i in ids(pm, :storage, nw=n_1)
        constraint_storage_state(pm, i, nw=n_1)
    end

    for n_2 in network_ids[2:end]
        for i in ids(pm, :storage, nw=n_2)
            constraint_storage_state(pm, i, n_1, n_2)
        end
        n_1 = n_2
    end

    objective_min_fuel_and_flow_cost(pm)
end


"Build multinetwork branch flow OPF with multiple flexibilities"
function build_mn_opf_bf_flex(pm::AbstractBFModelEdisgo)
    for (n, network) in nws(pm)
        # VARIABLES
        if ref(pm, 1, :opt_version) in(1, 2, 3, 4)
            variable_branch_power_radial(pm, nw=n)  # Eq. ():  branch power <= rate_a (s_nom)
            if ref(pm, 1, :opt_version) in(1, 3)
                variable_branch_current(pm, nw=n, bounded=false)
            else
                variable_branch_current(pm, nw=n)  # Eq. ()
            end
            variable_bus_voltage(pm, nw=n)  # Eq. (29)
            variable_gen_power_curt(pm, nw=n)  #  Eq. (20)
            variable_battery_storage_power(pm, nw=n)  # Eq. (21), (22)
            variable_heat_storage(pm, nw=n)  # Eq. (22)
            variable_cp_power(pm, nw=n)  #  Eq. (23), (24)
            variable_heat_pump_power(pm, nw=n)  # Eq. (25)
            variable_dsm_storage_power(pm, nw=n)  # Eq. (26), (27)
            variable_slack_gen(pm, nw=n)  # Eq. (28)
            variable_slack_grid_restrictions(pm, nw=n)
            variable_slack_HV_requirements(pm, nw=n)
        else
            throw(ArgumentError("OPF version $(ref(pm, 1, :opt_version)) is not implemented! Choose between version 1 to 4."))
        end
        
        # CONSTRAINTS
        for i in ids(pm, :bus, nw=n)  
            constraint_power_balance_bf(pm, i, nw=n) # Eq. (2)-(5)
        end
        for i in ids(pm, :branch, nw=n)
            constraint_voltage_magnitude_difference_radial(pm, i, nw=n) # Eq. (6)
        end
        constraint_model_current(pm, nw=n)  # Eq. (7) as SOC
        for i in ids(pm, :heatpumps, nw=n)  
            constraint_hp_operation(pm, i, n) # Eq. (14)
        end

        for i in ids(pm, :HV_requirements, nw=n)  
            constraint_HV_requirements(pm, i, n) # Eq. (15)-(19)
        end

    end

    # CONSTRAINTS
    network_ids = sort(collect(nw_ids(pm)))
    for kind in ["storage", "heat_storage", "dsm"]
        n_1 = network_ids[1]
        for i in ids(pm, Symbol(kind), nw=n_1)
            constraint_store_state(pm, i, nw=n_1, kind=kind)  # Eq. (8), (10)
        end

        for n_2 in network_ids[2:end]
            for i in ids(pm, Symbol(kind), nw=n_2)
                constraint_store_state(pm, i, n_1, n_2, kind) # Eq. (9), (11)
            end
            n_1 = n_2
        end
    end

    n_1 = network_ids[1]
    for i in ids(pm, :electromobility, nw=n_1)
        constraint_cp_state_initial(pm, n_1, i)  # Eq. (12)
    end

    for n_2 in network_ids[2:end]
        for i in ids(pm, :electromobility, nw=n_2)
            constraint_cp_state(pm, n_1, n_2, i) # Eq. (13)
        end
        n_1 = n_2
    end

    # OBJECTIVE FUNCTION
    if ref(pm, 1, :opt_version) in(1,3)
        objective_min_losses(pm)  # Eq. (1)
        if (ref(pm, 1, :opt_version) == 1)
            #objective_min_hv_slacks(pm)
            # Set multiple objectives
            # https://www.gurobi.com/documentation/9.1/refman/specifying_multiple_object.html
        end
    elseif ref(pm, 1, :opt_version) in(2,4)
        objective_min_losses_slacks(pm)  # Eq. (1)
        if (ref(pm, 1, :opt_version) == 2)
            #objective_min_hv_slacks(pm)
        end
    end
end
