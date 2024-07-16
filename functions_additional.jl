function λ_2_staple_residual_full(λ_2::Array{Float64,1},C::Array{Float64,1},p_B_vec::Array{Float64,1},p_M_vec::Array{Float64,1}
    ,ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64
    ,ϕ_B::Float64,c̄_S::Float64,Q_S::Float64,z::Array{Float64,1},
    ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,z_index::Int64)

    x_S_tmp = ((1 + τ_S) * p_x./(λ_2 * ζ * ϕ_S *z[z_index] )).^(1 / (ζ - 1));
    q_s_tmp = ϕ_S * x_S_tmp.^ζ *z[z_index];
    P_S_tmp = (λ_2.^(1 - ϵ)*ψ_S^ϵ .+ p_B_vec.^(1 -ϵ)*ψ_B^ϵ .+ p_M_vec.^(1 -ϵ)*ψ_M^ϵ).^(1 / (1 - ϵ));
    c_s_tmp = c̄_S .+ λ_2.^( -ϵ).*ψ_S.^ϵ.*C.*P_S_tmp.^ϵ
    return -(q_s_tmp - c_s_tmp).^2
end
function coeff_λ_2_s_approx_staples_full(coeff_λ_2_s::Array{Float64,2},C::Array{Float64,1},fspace_C_fine::Dict{Symbol,Any},p_B::Float64,p_M::Float64)
    cons_vec_local = repeat(C,outer = [ 1, 3])
    cons_vec_local[:,2] .= p_B
    cons_vec_local[:,3] .= p_M
    Phi_C_w_prices = funbase(fspace_joint_C, cons_vec_local);
    return diag(Phi_C_w_prices * coeff_λ_2_s);
end
function λ_2_staple_constrained(C::Array{Float64,1},θ::Array{Float64,1},
    fspace_C_fine::Dict{Symbol,Any},s::Array{Float64,2},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,
        p_B::Float64,p_M::Float64,ϕ_B::Float64,c̄_S::Float64,Q_S::Float64,
        ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,κ::Float64,ns::Int64,TC_S_c3_constrained::Array{Float64,1},
        c_S_c3_constrained::Array{Float64,1},a_min::Float64)
    condition = max.((((c_S_c3_constrained .- c̄_S)./ψ_S.^ϵ./C).^((1 - ϵ)/ϵ) .- ψ_S.^ϵ),a_min);
    λ_2_constrained_tmp = ((p_B^(1 - ϵ) *ψ_B^ϵ + p_M^(1 - ϵ) *ψ_M^ϵ ) ./condition).^(1 /(1 - ϵ));
    λ_S_constrained_tmp = (λ_2_constrained_tmp .* ζ .* ϕ_S .* θ)./ ((TC_S_c3_constrained ).^(1 - ζ) .* ((1 + τ_S) * p_x).^ζ) .- 1;
    #λ_S_constrained_tmp1 = (((c_S_c3_constrained .- c̄_S).^((1 - ϵ)/ϵ)*ψ_S.^(ϵ-1).*C.^(( ϵ - 1)/ϵ) .-  ψ_S.^ϵ).*c_S_c3_constrained.^((1- ζ)*(1 -ϵ )/ζ
    #)./((p_B^(1 - ϵ) *ψ_B^ϵ + p_M^(1 - ϵ) *ψ_M^ϵ ) .* ((1 + τ_S) * p_x).^( ϵ - 1))).^(1/( ϵ - 1)) .* ζ .* (ϕ_S .* θ).^(1/ζ) .- 1
    return λ_2_constrained_tmp,λ_S_constrained_tmp
end
function staples_c3_objects(C::Array{Float64,1},ϕ_S::Float64,τ_S::Float64,p_x::Float64,
    p_B::Float64,p_M::Float64,ϕ_B::Float64,τ_B::Float64,Q_S::Float64,ϵ::Float64,ψ_S::Float64,
    ψ_B::Float64,ψ_M::Float64,ζ::Float64,θ::Array{Float64,1},coeff_λ_2_s::Array{Float64,2},fspace_C_fine::Dict{Symbol,Any},
    s::Array{Float64,2},κ::Float64,c̄_S::Float64,C_max_staple::Array{Float64,1},C_min_staple::Array{Float64,1},ns::Int64,C_max_staple_constrained::Array{Float64,1},
    C_min_staple_constrained::Array{Float64,1},TC_S_c3_constrained::Array{Float64,1},
    x_S_c3_constrained::Array{Float64,1},q_S_c3_constrained::Array{Float64,1},c_S_c3_constrained::Array{Float64,1},a_min::Float64)
    TC_mat_S_c3 = zeros(ns,2);
    λ_2_mat_S_c3 = copy(TC_mat_S_c3);
    λ_S_mat_S_c3 = copy(TC_mat_S_c3);
    x_S_mat_S_c3  = copy(TC_mat_S_c3);
    q_S_mat_S_c3  = copy(TC_mat_S_c3);
    c_S_mat_S_c3  = copy(TC_mat_S_c3);
    feasibility_mat_S_c3 = copy(TC_mat_S_c3);
    λ_2_tmp = coeff_λ_2_s_approx_staples(coeff_λ_2_s,C,fspace_C_fine);
    λ_2_mat_S_c3[:,1] = λ_2_tmp;
    ( λ_2_mat_S_c3[:,2],λ_S_mat_S_c3[:,2]) = λ_2_staple_constrained(C,θ,fspace_C_fine,s,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,κ,ns,TC_S_c3_constrained,c_S_c3_constrained,a_min);
    TC_mat_S_c3[:,2] = TC_S_c3_constrained;
    x_S_mat_S_c3[:,2] = TC_S_c3_constrained./ ((1 + τ_S) * p_x)
    x_S_tmp = ((1 + τ_S) * p_x./(λ_2_tmp .* ζ .* ϕ_S .*θ )).^(1 / (ζ - 1));
    x_S_mat_S_c3[:,1] =x_S_tmp;
    TC_mat_S_c3[:,1]  = (1 + τ_S) * p_x * x_S_tmp;
    TC_true, constrained_staple_index = findmin(TC_mat_S_c3; dims=2);
    x_S_tmp = x_S_mat_S_c3[constrained_staple_index];
    q_S_tmp = ϕ_S * x_S_tmp.^ζ .*θ;
    c_s_tmp = q_S_tmp;
    λ_2_tmp = λ_2_mat_S_c3[constrained_staple_index];
    λ_S_tmp = λ_S_mat_S_c3[constrained_staple_index];
    shadow_price = min.(max.(1,λ_2_tmp),(1 + Q_S));
    P_S_c3_tmp = (shadow_price.^(1 - ϵ)*ψ_S^ϵ .+ p_B^(1 -ϵ)*ψ_B^ϵ .+ p_M^(1 -ϵ)*ψ_M^ϵ).^(1 / (1 - ϵ));
    c_B_tmp = p_B^( -ϵ)*ψ_B^ϵ * C .* P_S_c3_tmp.^ϵ;
    c_M_tmp = p_M^( -ϵ)*ψ_M^ϵ * C .* P_S_c3_tmp.^ϵ;
    Y_S_tmp = p_B * c_B_tmp + c_M_tmp * p_M + (1 + τ_S) * p_x * x_S_tmp;
    feasibility_mat_S_c3[:,1] =  (C_max_staple.<C) +  ((C_min_staple).>C);
    feasibility_mat_S_c3[:,2] = (C_max_staple_constrained.<C) +  ((C_min_staple_constrained).>C);
    feasibility = feasibility_mat_S_c3[constrained_staple_index] + (λ_2_tmp.>(1 + Q_S)) + (λ_2_tmp.<1) + (λ_S_tmp.<0.0 );
    return Y_S_tmp,P_S_c3_tmp,c_s_tmp,c_B_tmp,c_M_tmp,x_S_tmp,shadow_price,constrained_staple_index,feasibility
end




function coeff_creator(s::Array{Float64,2},ns::Int64,
    z::Array{Float64,1},z_W::Array{Float64,1},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,
    ϕ_B::Float64,τ_B::Float64,ρ::Float64,c̄_S::Float64,a_min::Float64,a_max::Float64,γ::Float64,n::Array{Int64,1},κ::Float64,Q_S::Float64,
    ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,agrid_fine::Array{Float64,1},fspace_C_fine::Dict{Symbol,Any},agrid::Array{Float64,1},
    tol::Float64 = 1e-8)

    s =Baseline_parameter.s
    ns=Baseline_parameter.ns
    ns_fine = Baseline_parameter.ns_fine #
    n_fine = Baseline_parameter.n_fine #
    z=Baseline_parameter.z
    z_W=Baseline_parameter.z_W
    ϕ_S = Baseline_parameter.ϕ_S
    ζ = Baseline_parameter.ζ
    τ_S = Baseline_parameter.τ_S
    p_x = Baseline_parameter.p_x
    ϕ_B = Baseline_parameter.ϕ_B
    τ_B = Baseline_parameter.τ_B
    ρ = Baseline_parameter.ρ
    c̄_S = Baseline_parameter.c̄_S
    a_min = Baseline_parameter.a_min
    a_max = Baseline_parameter.a_max
    γ = Baseline_parameter.γ
    n = Baseline_parameter.n
    κ = Baseline_parameter.κ
    Q_S = Baseline_parameter.Q_S
    ϵ = Baseline_parameter.ϵ
    ψ_S = Baseline_parameter.ψ_S
    ψ_B = Baseline_parameter.ψ_B
    ψ_M = Baseline_parameter.ψ_M
    agrid_fine = Baseline_parameter.agrid_fine
    fspace_C_fine = Baseline_parameter.fspace_C_fine
    agrid = Baseline_parameter.agrid
    ones_agrid_fine = 1 .* agrid_fine

    no_labor_shock = size(z_W)[1];
    no_prod_shock = convert(Int64,n[2]/no_labor_shock);
    θ = z[((convert(Array{Int64,1},s[:,2]).-1) .% no_prod_shock .+1)];
    labor_prod = z_W[convert(Array{Int64,1},floor.((s[:,2]/ no_prod_shock .-0.01) .+ 1))];
    ones_tmp = ones(ns);
    p_B_max = 4.0
    p_B_min = 0.1
    p_B_grid_no = 20;
    p_M_max = 0.5
    p_M_min = 0.05
    p_M_grid_no = 20;
    # p_B_max = 0.9210526315789473
    # p_B_min = 0.9210526315789473
    # p_B_grid_no = 1;
    # p_M_max = 0.07368421052631577
    # p_M_min = 0.07368421052631577
    # p_M_grid_no = 1;


    p_B_grid = range(p_B_min,p_B_max,length=p_B_grid_no);
    p_M_grid = range(p_M_min,p_M_max,length=p_M_grid_no);
    

    fspace_p_B = fundef((:spli, p_B_grid, 0,1));
    fspace_p_M = fundef((:spli, p_M_grid, 0,1));
    fspace_joint_C = fundef((:spli, Baseline_parameter.C_grid_fine, 0,1),(:spli, p_B_grid, 0,1),(:spli, p_M_grid, 0,1));# joint
    cons_state_vec    = funnode(fspace_joint_C)[1];
    n_cons_state = size(cons_state_vec);

    Phi_C_fine = funbase(fspace_C_fine, cons_state_vec[:,1]);
    Phi_p_B = funbase(fspace_p_B, cons_state_vec[:,2]);
    Phi_p_M = funbase(fspace_p_M, cons_state_vec[:,3]);

    Phi_λ = funbase(fspace_joint_C, cons_state_vec)  

    #Phi_λ1 = row_kron(Phi_p_M,row_kron(Phi_p_B,Phi_C_fine)); # ALWAYS DO IT IN REVERSE ORDER of the states!

    ones_λ = ones(n_cons_state[1])

    # Hopefully only the c_S=q_S case is necessary to be handled here:

    # Stores the consumption solutions at the double grid.
    # Each column is for different θ, row is for different C. Solution must be for each C,p_B,p_M combination
    ones_λ = ones(n_cons_state[1]);
    ones_tmp = ones(ns);
    ones_tmp_fine = ones(ns_fine);
    coeff_λ_2_s = zeros(n_cons_state[1],ns);
    Value_λ_2_s = zeros(n_cons_state[1],ns);
    V_temp_staple_residual = zeros(n_cons_state[1],ns);
    exit_flag_mat = zeros(n_cons_state[1],ns);
    TC_S_c3_constrained = zeros(ns);
    x_S_c3_constrained = zeros(ns);
    q_S_c3_constrained = zeros(ns);
    c_S_c3_constrained = zeros(ns);
    exit_flag_mat_fine = zeros(n_cons_state[1],ns_fine);
    TC_S_c3_constrained_fine = zeros(ns_fine);
    x_S_c3_constrained_fine = zeros(ns_fine);
    q_S_c3_constrained_fine = zeros(ns_fine);
    c_S_c3_constrained_fine = zeros(ns_fine);
    C_tmp = cons_state_vec[:,1];
    p_B_vec = cons_state_vec[:,2];
    p_M_vec = cons_state_vec[:,3];

    for z_index= 1:(convert(Int64,n[2]/no_labor_shock))

        (Value_λ_2_s_tmp,V_temp) = goldenx(λ_2_staple_residual_full,ones_λ,(1 + Q_S) *ones_λ,tol,C_tmp,p_B_vec,p_M_vec,ϕ_S,ζ,τ_S,p_x,ϕ_B,c̄_S,Q_S,z,ϵ,ψ_S,ψ_B,ψ_M,z_index);
        θ_index_start = (z_index-1)*n[1] + 1;
        θ_index_end = z_index*n[1];
        labor_shock_same_index= convert(Int64,ns/no_labor_shock);
        for l_index=0:(no_labor_shock-1)
            Value_λ_2_s[:,(labor_shock_same_index*l_index+ θ_index_start):(labor_shock_same_index*l_index+θ_index_end)].=Value_λ_2_s_tmp;
            V_temp_staple_residual[:,(labor_shock_same_index*l_index+ θ_index_start):(labor_shock_same_index*l_index+θ_index_end)].=V_temp;
            coeff_λ_2_s[:,(labor_shock_same_index*l_index+ θ_index_start):(labor_shock_same_index*l_index+θ_index_end)].=Phi_λ\Value_λ_2_s_tmp
        end
        for a_index= 1:n[1]
            # Constrained case: bounds
            TC_S_constrained_tmp = κ*agrid[a_index];
            x_S_constrained_tmp = TC_S_constrained_tmp/((1 + τ_S) * p_x)
            q_S_constrained_tmp = ϕ_S * x_S_constrained_tmp.^ζ .*z[z_index];
            condition=(((q_S_constrained_tmp .- c̄_S)./ψ_S.^ϵ./C_tmp).^((1 - ϵ)/ϵ) .- ψ_S.^ϵ);
            exitflag_tmp = 0* copy(ones_λ);
            exitflag_tmp[condition.<0].=-1;
            for l_index=0:(no_labor_shock-1)
                exit_flag_mat[:,(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))].=exitflag_tmp;
                TC_S_c3_constrained[(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))] = TC_S_constrained_tmp;
                x_S_c3_constrained[(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))] = x_S_constrained_tmp;
                q_S_c3_constrained[(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))] = q_S_constrained_tmp;
                c_S_c3_constrained[(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))] = q_S_constrained_tmp;
            end
        end
        #Now do it for the finer grid:
        for a_index= 1:n_fine[1]
            # Constrained case: bounds
            TC_S_constrained_tmp_fine = κ*agrid_fine[a_index];
            x_S_constrained_tmp_fine = TC_S_constrained_tmp_fine/((1 + τ_S) * p_x)
            q_S_constrained_tmp_fine = ϕ_S * x_S_constrained_tmp_fine.^ζ .*z[z_index];
            condition=(((q_S_constrained_tmp_fine .- c̄_S)./ψ_S.^ϵ./C_tmp).^((1 - ϵ)/ϵ) .- ψ_S.^ϵ);
            exitflag_tmp_fine = 0* copy(ones_λ);
            exitflag_tmp_fine[condition.<0].=-1;
            for l_index=0:(no_labor_shock-1)
                exit_flag_mat_fine[:,(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))].=exitflag_tmp_fine;
                TC_S_c3_constrained_fine[(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))] = TC_S_constrained_tmp_fine;
                x_S_c3_constrained_fine[(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))] = x_S_constrained_tmp_fine;
                q_S_c3_constrained_fine[(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))] = q_S_constrained_tmp_fine;
                c_S_c3_constrained_fine[(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))] = q_S_constrained_tmp_fine;
            end
        end
        # for this simple case, this isnt needed as c_s_solve_mat = coeff_λ_2_s!!!
    end

    # Consumption bounds - do it for staples:
    C_grid_mat = repeat(C_tmp,1,ns);
    C_max_staple = maximum(C_grid_mat .* (V_temp_staple_residual.> -tol),dims = 1)[:] .+ a_min;
    C_max_staple_constrained = maximum(C_grid_mat .* (exit_flag_mat.== 0),dims = 1)[:];
    C_min_staple_constrained = a_min * ones_tmp;
    tmp2=zeros(Int64,ns);
    for ii=1:ns
        if isnothing(findfirst(myCondition, (C_grid_mat .* (V_temp_staple_residual.> -tol))[:,ii]))
            tmp2[ii]=1
        else
            tmp2[ii]=findfirst(myCondition, (C_grid_mat .* (V_temp_staple_residual.> -tol))[:,ii])
        end
    end
    C_min_staple = C_grid_mat[tmp2];
    C_grid_fine_mat = repeat(C_tmp,1,ns_fine);
    C_max_staple_constrained_fine = maximum(C_grid_fine_mat .* (exit_flag_mat_fine.== 0),dims = 1)[:];
    C_min_staple_constrained_fine = a_min * ones_tmp_fine;

    # Create functions, use staples_objects later.

    # Cash crop farmer
    #π_possible = zeros(ns,3);

    # #3c) q_S=c_S, hence transaction cost is not paid, but production is distorted, see functions 
    # coeff_λ_2_cashcrop_residual_unconstrained = zeros(n_cons_state[1],ns);
    # coeff_λ_2_cashcrop_residual_constrained = zeros(n_cons_state[1],ns);
    # V_temp_cashcrop_residual_unconstrained = zeros(n_cons_state[1],ns);
    # V_temp_cashcrop_residual_constrained = zeros(n_cons_state[1],ns);
    #     #    z_index = 1
    #     #    c_S_cashcrop_residual_unconstrained((c̄_S+tol) * ones_agrid_fine,ϕ_S,ζ,
    #     #    τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,z,ϵ,ψ_S,ψ_B,ψ_M,C_grid_fine,z_index,ones_agrid_fine,agrid_fine,tol,1)
    #     for z_index= 1:(convert(Int64,n[2]/no_labor_shock))
    #          (coeff_λ_2_cashcrop_residual_unconstrained_tmp,V_temp) = goldenx(λ_2_residual_unconstrained,  ones_agrid_fine
    #          ,(1 + Q_S) *ones_agrid_fine ,tol,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,z,ϵ,ψ_S,ψ_B,ψ_M,C_grid_fine,z_index,ρ);
    
    #         θ_index_start = (z_index-1)*n[1] + 1;
    #         θ_index_end = z_index*n[1];
    #         labor_shock_same_index= convert(Int64,ns/no_labor_shock);
    #         for l_index=0:(no_labor_shock-1)
    #             coeff_λ_2_cashcrop_residual_unconstrained[:,(labor_shock_same_index*l_index+ θ_index_start):(labor_shock_same_index*l_index+θ_index_end)].=coeff_λ_2_cashcrop_residual_unconstrained_tmp;
    #             V_temp_cashcrop_residual_unconstrained[:,(labor_shock_same_index*l_index+ θ_index_start):(labor_shock_same_index*l_index+θ_index_end)].=V_temp;
    #         end
    #         #unconstrained ends, constrained begins:
    #         for a_index= 1:n[1]
    #             (coeff_λ_2_cashcrop_residual_constrained_tmp,V_temp) = goldenx(λ_2_residual_constrained,  ones_agrid_fine
    #             ,(1 + Q_S) *ones_agrid_fine ,tol,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,z,ϵ,ψ_S,ψ_B,ψ_M,C_grid_fine,z_index,agrid,a_index,κ,ρ);
    #             for l_index=0:(no_labor_shock-1)
    #                 coeff_λ_2_cashcrop_residual_constrained[:,(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))].=coeff_λ_2_cashcrop_residual_constrained_tmp;
    #                 V_temp_cashcrop_residual_constrained[:,(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))].=V_temp;
    #             end
    #         end
    #     end

# ind1=3890 #intr range 3890:4000
# ind2_1=1
# ind2_2 = 55
#     ind2_3 = 120
#     ind2_4 = 302



#     coeff_λ_2_s[3895:3920, 1:30]
#     C_tmp = cons_state_vec[3895:3920, 1]
#     p_B_tmp = cons_state_vec[3895:3920, 2]
#     p_M_tmp = cons_state_vec[3895:3920, 3]

coeff_λ_2_s
Phi_prime_a = funbase(fspace_a_co, future_asset);

end