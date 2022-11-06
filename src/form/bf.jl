# this file contains (balanced) convexified DistFlow formulation, in W space

""
function variable_buspair_current_magnitude_sqr(pm::AbstractBFModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    branch = ref(pm, nw, :branch)

    ccm = var(pm, nw)[:ccm] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :branch)], base_name="$(nw)_ccm",
        start = comp_start_value(branch[i], "ccm_start")
    )

    if bounded
        bus = ref(pm, nw, :bus)
        for (i, b) in branch
            rate_a = Inf
            if haskey(b, "rate_a")
                rate_a = b["rate_a"]
            end
            ub = ((rate_a*b["tap"])/(bus[b["f_bus"]]["vmin"]))^2

            JuMP.set_lower_bound(ccm[i], 0.0)
            if !isinf(ub)
                JuMP.set_upper_bound(ccm[i], ub)
            end
        end
    end

    report && sol_component_value(pm, nw, :branch, :ccm, ids(pm, nw, :branch), ccm)
end

""
function variable_branch_current(pm::AbstractBFModel; kwargs...)
    variable_buspair_current_magnitude_sqr(pm; kwargs...)
end

""
function variable_bus_voltage(pm::AbstractBFModel; kwargs...)
    variable_bus_voltage_magnitude_sqr(pm; kwargs...)
end

""
# function variable_slack_gen(pm::SOCBFPowerModelEdisgo; kwargs...)
#     pgs = var(pm, nw)[:pgs] = JuMP.@variable(pm.model,
#         [i in ids(pm, nw, :gen)], base_name="$(nw)_pgs"
#     )
#     report && sol_component_value(pm, nw, :gen, :pgs, ids(pm, nw, :gen), pgs)
# end


"""
Defines branch flow model power flow equations
"""
function constraint_power_losses(pm::AbstractBFModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to, tm)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)
    w_fr = var(pm, n, :w, f_bus)
    w_to = var(pm, n, :w, t_bus)
    ccm =  var(pm, n, :ccm, i)

    ym_sh_sqr = g_sh_fr^2 + b_sh_fr^2

    JuMP.@constraint(pm.model, p_fr + p_to == r*(ccm + ym_sh_sqr*(w_fr/tm^2) - 2*(g_sh_fr*p_fr - b_sh_fr*q_fr)) + g_sh_fr*(w_fr/tm^2) + g_sh_to*w_to)
    JuMP.@constraint(pm.model, q_fr + q_to == x*(ccm + ym_sh_sqr*(w_fr/tm^2) - 2*(g_sh_fr*p_fr - b_sh_fr*q_fr)) - b_sh_fr*(w_fr/tm^2) - b_sh_to*w_to)
end


"""
Defines voltage drop over a branch, linking from and to side voltage magnitude
"""
function constraint_voltage_magnitude_difference(pm::AbstractBFModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    w_fr = var(pm, n, :w, f_bus)
    w_to = var(pm, n, :w, t_bus)
    ccm =  var(pm, n, :ccm, i)

    ym_sh_sqr = g_sh_fr^2 + b_sh_fr^2

    JuMP.@constraint(pm.model, (1+2*(r*g_sh_fr - x*b_sh_fr))*(w_fr/tm^2) - w_to ==  2*(r*p_fr + x*q_fr) - (r^2 + x^2)*(ccm + ym_sh_sqr*(w_fr/tm^2) - 2*(g_sh_fr*p_fr - b_sh_fr*q_fr)))
end


function constraint_voltage_magnitude_difference(pm::SOCBFPowerModelEdisgo, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, tm)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    w_fr = var(pm, n, :w, f_bus)
    w_to = var(pm, n, :w, t_bus)
    ccm =  var(pm, n, :ccm, i)

    JuMP.@constraint(pm.model, ((w_fr/tm^2)  - w_to ==  2*(r*p_fr + x*q_fr) - (r^2 + x^2)*ccm))
end

"""
Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude
"""
function constraint_model_current(pm::AbstractBFQPModel, n::Int)
    _check_missing_keys(var(pm, n), [:p,:q,:w,:ccm], typeof(pm))

    p  = var(pm, n, :p)
    q  = var(pm, n, :q)
    w  = var(pm, n, :w)
    ccm = var(pm, n, :ccm)

    for (i,branch) in ref(pm, n, :branch)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        f_idx = (i, f_bus, t_bus)
        tm = branch["tap"]

        JuMP.@constraint(pm.model, p[f_idx]^2 + q[f_idx]^2 <= (w[f_bus]/tm^2)*ccm[i])
    end
end


"""
Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude
"""
function constraint_model_current(pm::AbstractBFConicModel, n::Int)
    _check_missing_keys(var(pm, n), [:p,:q,:w,:ccm], typeof(pm))

    p  = var(pm, n, :p)
    q  = var(pm, n, :q)
    w  = var(pm, n, :w)
    ccm = var(pm, n, :ccm)

    for (i,branch) in ref(pm, n, :branch)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        f_idx = (i, f_bus, t_bus)
        tm = branch["tap"]

        JuMP.@constraint(pm.model, [w[f_bus]/tm^2, ccm[i]/2, p[f_idx], q[f_idx]] in JuMP.RotatedSecondOrderCone())
    end
end


function constraint_voltage_angle_difference(pm::AbstractBFModel, n::Int, f_idx, angmin, angmax)
    i, f_bus, t_bus = f_idx
    t_idx = (i, t_bus, f_bus)

    branch = ref(pm, n, :branch, i)
    tm = branch["tap"]
    g_fr = branch["g_fr"]
    g_to = branch["g_to"]
    b_fr = branch["b_fr"]
    b_to = branch["b_to"]

    tr, ti = calc_branch_t(branch)

    r = branch["br_r"]
    x = branch["br_x"]

    # getting the variables
    w_fr = var(pm, n, :w, f_bus)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)

    tzr = r*tr + x*ti
    tzi = r*ti - x*tr

    JuMP.@constraint(pm.model,
        tan(angmin)*((tr + tzr*g_fr + tzi*b_fr)*(w_fr/tm^2) - tzr*p_fr + tzi*q_fr)
                 <= ((ti + tzi*g_fr - tzr*b_fr)*(w_fr/tm^2) - tzi*p_fr - tzr*q_fr)
        )
    JuMP.@constraint(pm.model,
        tan(angmax)*((tr + tzr*g_fr + tzi*b_fr)*(w_fr/tm^2) - tzr*p_fr + tzi*q_fr)
                 >= ((ti + tzi*g_fr - tzr*b_fr)*(w_fr/tm^2) - tzi*p_fr - tzr*q_fr)
        )
end

"""
Defines linear branch flow model power flow equations
"""
function constraint_power_losses(pm::AbstractBFAModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to, tm)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)
    w_fr = var(pm, n, :w, f_bus)
    w_to = var(pm, n, :w, t_bus)

    JuMP.@constraint(pm.model, p_fr + p_to ==  g_sh_fr*(w_fr/tm^2) + g_sh_to*w_to)
    JuMP.@constraint(pm.model, q_fr + q_to == -b_sh_fr*(w_fr/tm^2) - b_sh_to*w_to)
end


"""
Defines voltage drop over a branch, linking from and to side voltage magnitude
"""
function constraint_voltage_magnitude_difference(pm::AbstractBFAModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    w_fr = var(pm, n, :w, f_bus)
    w_to = var(pm, n, :w, t_bus)

    JuMP.@constraint(pm.model, (w_fr/tm^2) - w_to ==  2*(r*p_fr + x*q_fr))
end


""
function constraint_model_current(pm::AbstractBFAModel, n::Int)

end

""
function variable_buspair_current_magnitude_sqr(pm::AbstractBFAModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)

end


"Neglects the active and reactive loss terms associated with the squared current magnitude."
function constraint_storage_losses(pm::AbstractBFAModel, n::Int, i, bus, r, x, p_loss, q_loss; conductors=[1])
    ps = var(pm, n, :ps, i)
    qs = var(pm, n, :qs, i)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)
    qsc = var(pm, n, :qsc, i)


    JuMP.@constraint(pm.model,
        sum(ps[c] for c in conductors) + (sd - sc)
        ==
        p_loss
    )

    JuMP.@constraint(pm.model,
        sum(qs[c] for c in conductors)
        ==
        qsc + q_loss
    )
end


""
function constraint_power_balance(pm::SOCBFPowerModelEdisgo, n::Int, i, bus_gens_nd, bus_arcs_to, bus_arcs_from, bus_lines_to, bus_storage, bus_pg, bus_qg, bus_pg_nd, bus_qg_nd, bus_pd, bus_qd, bus_gs, bus_bs, branch_r, branch_x)
    w    = var(pm, n, :w, i)
    pt   = get(var(pm, n),  :p, Dict()); _check_var_keys(pt, bus_arcs_to, "active power", "branch")
    qt   = get(var(pm, n),  :q, Dict()); _check_var_keys(qt, bus_arcs_to, "reactive power", "branch")
    pf   = get(var(pm, n),  :p, Dict()); _check_var_keys(pf, bus_arcs_from, "active power", "branch")
    qf   = get(var(pm, n),  :q, Dict()); _check_var_keys(qf, bus_arcs_from, "reactive power", "branch")
    ps   = get(var(pm, n),  :ps, Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, n),  :qs, Dict()); _check_var_keys(qs, bus_storage, "reactive power", "storage")
    pgc  = get(var(pm, n),  :pgc, Dict()); _check_var_keys(pgc, bus_gens_nd, "active power", "curtailment")
    qgc  = get(var(pm, n),  :qgc, Dict()); _check_var_keys(qgc, bus_gens_nd, "reactive power", "curtailment")
    ccm  = get(var(pm, n),  :ccm, Dict()); _check_var_keys(ccm, bus_lines_to, "active power", "branch")

    # TODO: add Slack_gen, ccm * R / ccm * X

    cstr_p = JuMP.@constraint(pm.model,
        sum(pt[a] for a in bus_arcs_to)
        ==
        sum(pf[a] for a in bus_arcs_from)
        + sum(ccm[a] * branch_r[a] for a in bus_lines_to) # *R
        - sum(pg for pg in values(bus_pg))
        - sum(pg for pg in values(bus_pg_nd))
        + sum(ps[s] for s in bus_storage)
        + sum(pd for pd in values(bus_pd))
        #+ sum(gs for gs in values(bus_gs))*w   # TODO: + oder -?
        + sum(pgc[g] for g in bus_gens_nd)  # TODO: + oder -?
    )
    cstr_q = JuMP.@constraint(pm.model,
        sum(qt[a] for a in bus_arcs_to)
        ==
        sum(qf[a] for a in bus_arcs_from)
        + sum(ccm[a] * branch_x[a] for a in bus_lines_to) # *X
        - sum(qg for qg in values(bus_qg))
        - sum(qg for qg in values(bus_qg_nd))
        + sum(qs[s] for s in bus_storage)
        + sum(qd for qd in values(bus_qd))
        #- sum(bs for bs in values(bus_bs))*w   # TODO: + oder -?
        - sum(qgc[g] for g in bus_gens_nd)  # TODO: + oder -?
    )

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end