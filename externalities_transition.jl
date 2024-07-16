
function Residual_transition_sequential_epsilon(price_trans_actual::Array{Float64,2},capital_trans::Array{Float64,1},
    distr_store::Array{Float64,2},T::Int64,parameter_end::Parameter_type,coeff_store::Array{Float64,3},τ_trans::Array{Float64,1},r::Float64,cons_level_substinence::Float64,epsilon_u::Float64,epsilon_r::Float64)
    # Extract the necessary local parameters:
    (δ,ζ,ρ,α,σ,β,ϵ,ψ_S,ψ_B,ψ_M,ϕ_S,ϕ_B,c̄_S,F_W,F_S,F_B,FM_W,FM_S,FM_B,Q_S,p_x,τ_S,τ_B,a_D,b_D,K_a,K_b,γ,A_W,
    ρ_S,ρ_SW,σ_S,ρ_W,σ_W,n,n_fine,agrid,agrid_fine,a_min,a_max,spliorder,fspace_a,fspace_a_fine,fspace_C_fine,C_grid_fine_no,C_grid_fine,s,ns,
    s_fine,ns_fine,z,z_W,Phi_z,Phi_z_fine,Phi,Phi_aug,P_kron,P_kron1,P_kron_fine,κ) = local_parameters(parameter_end);
    # Policy function
    (   a_prime_fine_store,future_occupation_fine_local_store,
    cons_fine_local_store,c_S_W_fine_store,c_B_W_fine_store,c_M_W_fine_store,Y_W_fine_policy_store,Y_S_fine_store,
    c_S_S_fine_store,c_B_S_fine_store,c_M_S_fine_store,q_S_S_fine_store,P_S_fine_store,
    x_S_S_fine_store,solve_staple_index_S_fine_store,λ_2_S_fine_store,future_asset_S_store,
    future_asset_C_store,c_S_B_fine_store,c_B_B_fine_store,c_M_B_fine_store,x_SC_fine_store,x_BC_fine_store,
    land_C_fine_store,λ_2_fine_store,P_B_fine_store,Y_B_fine_store,q_S_B_fine_store,
    q_B_B_fine_store,solve_cash_crop_index_B_fine_store,solve_staple_index_B_fine_store,TC_fine_store,residual_store,
    coeff_λ_2_cashcrop_residual_unconstrained_store,coeff_λ_2_cashcrop_residual_constrained_store,C_max_unconstrained_store,
    C_max_constrained_store,C_min_unconstrained_store,C_min_constrained_store,coeff_λ_2_s_store,C_max_staple_store,C_min_staple_store,
    C_max_staple_constrained_store,C_min_staple_constrained_store,TC_S_c3_constrained_store,x_S_c3_constrained_store,q_S_c3_constrained_store,
    c_S_c3_constrained_store) =local_var_creator_policy(ns,ns_fine,T,C_grid_fine_no);
    #Precheck the prices
    price_check_sum = 0;
    println("Backward loop")
    for t in (T+1):-1:2
        # Prices and other time dependent quantities:
        τ_S_loc = τ_trans[t];
        prices_loc = price_trans_actual[:,t-1];
        # Epsilon distortions
        cttilde=prices_loc[4]

        if (τ_S_loc==0 || prices_loc[3]==0.0)
            balanced_share = 0.0;
        else
            balanced_share = 1.0;
        end
        p_B,p_M,R,r,w,τ_W= price_reshaper_fixed_r(prices_loc[1:3],δ,ϵ,ψ_S,ψ_B,ψ_M,p_x,τ_S_loc,τ_B,α,balanced_share,r);
        price_check_tmp = p_B<p_B_min + p_B>p_B_max + p_M<p_M_min + p_M>p_M_max + τ_W<τ_W_min + τ_W>τ_W_max;
        # Price checks might be necessary - leave it for now.
        #0.01<p_B
        price_check_sum = price_check_tmp + price_check_sum;
        #price_check_sum = 0
        # Get the coefficients
        if price_check_sum>0.0
            print(t-1,",")
        else
            coefficients_next_tmp = coeff_store[:,:,t+1];
            coeff_next = coeff_store[:,:,t];
            (coeff_store[:,:,t],coeff_λ_2_cashcrop_residual_unconstrained_store[:,:,t-1],coeff_λ_2_cashcrop_residual_constrained_store[:,:,t-1],
            C_max_unconstrained_store[:,t-1],C_max_constrained_store[:,t-1],C_min_unconstrained_store[:,t-1],C_min_constrained_store[:,t-1],
            coeff_λ_2_s_store[:,:,t-1],C_max_staple_store[:,t-1],C_min_staple_store[:,t-1],C_max_staple_constrained_store[:,t-1],
            C_min_staple_constrained_store[:,t-1],TC_S_c3_constrained_store[:,t-1],x_S_c3_constrained_store[:,t-1],q_S_c3_constrained_store[:,t-1],
            c_S_c3_constrained_store[:,t-1] )= Residual_transition_backward_epsilon(s,ns,
            z,z_W,ϕ_S,ζ,τ_S_loc,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n,κ,Q_S,ϵ,ψ_S,ψ_B,ψ_M,agrid_fine,
            fspace_C_fine,agrid,coefficients_next_tmp,coeff_next,C_grid_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,P_kron,Phi,Phi_z,β,fspace_a,σ,cons_level_substinence, epsilon_u, epsilon_r,cttilde)
        end
        print(t-1,",")
    end

    if price_check_sum>0
        println("Guess out of bounds")
        residual_store = 1000*ones(6,T);
    else
        #Forward loop
        println("Forward loop")
        for t in 2:(T+1)
            # Prices and other time dependent quantities:
            τ_S_loc = τ_trans[t];
            prices_loc = price_trans_actual[:,t-1];
            # Epsilon distortions
            cttilde=prices_loc[4]

            if (τ_S_loc==0 || prices_loc[3]==0.0)
                balanced_share = 0.0;
            else
                balanced_share = 1.0;
            end
            p_B,p_M,R,r,w,τ_W= price_reshaper_fixed_r(prices_loc[1:3],δ,ϵ,ψ_S,ψ_B,ψ_M,p_x,τ_S_loc,τ_B,α,balanced_share,r);
             (residual_store[:,t-1],distr_store[:,t]) = Residual_transition_forward_epsilon(coeff_store[:,:,t],distr_store[:,t-1],s_fine,ns_fine,
                 z,z_W,ϕ_S,ζ,τ_S_loc,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n_fine,κ,Q_S,
                 ϵ,ψ_S,ψ_B,ψ_M,coeff_λ_2_cashcrop_residual_unconstrained_store[:,:,t-1],
                 coeff_λ_2_cashcrop_residual_constrained_store[:,:,t-1],
                 C_max_unconstrained_store[:,t-1] ,C_max_constrained_store[:,t-1],C_min_unconstrained_store[:,t-1],
                 C_min_constrained_store[:,t-1],coeff_λ_2_s_store[:,:,t-1],agrid_fine,fspace_C_fine,C_max_staple_store[:,t-1],
                 C_min_staple_store[:,t-1],C_max_staple_constrained_store[:,t-1],
                 C_min_staple_constrained_store[:,t-1],TC_S_c3_constrained_store[:,t-1],
                 x_S_c3_constrained_store[:,t-1],q_S_c3_constrained_store[:,t-1],c_S_c3_constrained_store[:,t-1],capital_trans[t],
                 balanced_share,τ_W,C_grid_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,P_kron_fine,Phi_z_fine,β,fspace_a,σ,P_kron1,fspace_a_fine,a_D,b_D,R,δ,α,cons_level_substinence, epsilon_u, epsilon_r,cttilde);
            print( t-1,",")
        end
    end
    return residual_store,distr_store,coeff_store
end


function Residual_transition_forward_epsilon(coeff::Array{Float64,2},distr_previous::Array{Float64,1},s_fine::Array{Float64,2},ns_fine::Int64,
    z::Array{Float64,1},z_W::Array{Float64,1},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,
    p_B::Float64,p_M::Float64,ϕ_B::Float64,τ_B::Float64,ρ::Float64,w::Float64,r::Float64,
    c̄_S::Float64,a_min::Float64,a_max::Float64,γ::Float64,n_fine::Array{Int64,1},κ::Float64,Q_S::Float64,
    ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,
    coeff_λ_2_cashcrop_residual_unconstrained::Array{Float64,2},
    coeff_λ_2_cashcrop_residual_constrained::Array{Float64,2},
    C_max_unconstrained::Array{Float64,1} ,C_max_constrained::Array{Float64,1},C_min_unconstrained::Array{Float64,1},
    C_min_constrained::Array{Float64,1},
    coeff_λ_2_s::Array{Float64,2},agrid_fine::Array{Float64,1},fspace_C_fine::Dict{Symbol,Any},C_max_staple::Array{Float64,1},
    C_min_staple::Array{Float64,1},C_max_staple_constrained::Array{Float64,1},
    C_min_staple_constrained::Array{Float64,1},TC_S_c3_constrained::Array{Float64,1},
    x_S_c3_constrained::Array{Float64,1},q_S_c3_constrained::Array{Float64,1},c_S_c3_constrained::Array{Float64,1},foreign_supply_capital::Float64,
    balanced_share::Float64,τ_W::Float64,C_grid_fine::Array{Float64,1},F_W::Float64,F_S::Float64,F_B::Float64,FM_W::Float64,FM_S::Float64,FM_B::Float64,
    P_kron_fine::SparseMatrixCSC{Float64, Int64},
    Phi_z_fine::SparseMatrixCSC{Float64, Int64},β::Float64,
    fspace_a::Dict{Symbol,Any},σ::Float64,P_kron1::SparseMatrixCSC{Float64, Int64},fspace_a_fine::Dict{Symbol,Any},a_D::Float64,b_D::Float64,
    R::Float64,δ::Float64,α::Float64,cons_level_substinence::Float64, epsilon_u::Float64, epsilon_r::Float64,cttilde::Float64,tol::Float64 = 1e-8)
   
# coeff = coeff_store[:,:,t]; 
# distr_previous = distr_store[:,t]; 
# τ_S = τ_S_loc;
# coeff_λ_2_cashcrop_residual_unconstrained = coeff_λ_2_cashcrop_residual_unconstrained_store[:,:,t-1]
# coeff_λ_2_cashcrop_residual_constrained = coeff_λ_2_cashcrop_residual_constrained_store[:,:,t-1]
# C_max_unconstrained = C_max_unconstrained_store[:,t-1]
# C_max_constrained = C_max_constrained_store[:,t-1]
# C_min_unconstrained = C_min_unconstrained_store[:,t-1]
# C_min_constrained = C_min_constrained_store[:,t-1]
# coeff_λ_2_s = coeff_λ_2_s_store[:,:,t-1]
# C_max_staple = C_max_staple_store[:,t-1]
# C_min_staple = C_min_staple_store[:,t-1]
# C_max_staple_constrained = C_max_staple_constrained_store[:,t-1]
# C_min_staple_constrained = C_min_staple_constrained_store[:,t-1]
# TC_S_c3_constrained = TC_S_c3_constrained_store[:,t-1]
# x_S_c3_constrained = x_S_c3_constrained_store[:,t-1]
# q_S_c3_constrained = q_S_c3_constrained_store[:,t-1]
# c_S_c3_constrained = c_S_c3_constrained_store[:,t-1]
# foreign_supply_capital = capital_trans[t]

    ctilde = copy(cons_level_substinence)
    ϕ_S_loc = ϕ_S * exp(-epsilon_r * (cttilde - 0.2))
    ϕ_B_loc = ϕ_B * exp(-epsilon_r * (cttilde - 0.2))
    # Aggregates
        (θ_fine,labor_prod_fine,tol,P_W_fine,Y_W_fine,coeff_λ_2_cashcrop_residual_unconstrained_fine,coeff_λ_2_cashcrop_residual_constrained_fine,
        x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
        coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,labor_allocated_interior_c3a_fine,
        λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,q_S_c2_fine,q_B_c1_fine,q_S_c3a_fine,
        q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,Y_B_c3b_fine, c_S_mat_fine,c_B_mat_fine,
        c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,q_B_mat_fine,land_B_mat_fine, λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,feasibility_mat_fine,C_max_mat_fine,
        C_min_mat_fine,q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,
        Y_S_potential_fine,TC_mat_fine,C_max_staple_fine,C_min_staple_fine,C_max_staple_constrained_fine,
        C_min_staple_constrained_fine,TC_S_c3_constrained_fine,x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine, x_S_mat_3c_fine,x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine) =  income_creator_no_approx_ext(s_fine,ns_fine,
        z,z_W,ϕ_S_loc,ζ,τ_S,p_x,p_B,p_M,ϕ_B_loc,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n_fine,κ,Q_S,ϵ,ψ_S,ψ_B,ψ_M,coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,
        C_max_unconstrained ,C_max_constrained,C_min_unconstrained,C_min_constrained, coeff_λ_2_s,C_grid_fine,fspace_C_fine,C_max_staple,C_min_staple,C_max_staple_constrained,
    C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained, epsilon_u, cttilde, ctilde);


       min_C_applied_fine,max_C_applied_fine = bounds_consumption(P_W_fine,Y_W_fine,s_fine,r,ρ,w,
       coeff_λ_2_cashcrop_residual_unconstrained_fine,coeff_λ_2_cashcrop_residual_constrained_fine,θ_fine,
       fspace_C_fine,ϕ_S_loc,ζ,τ_S,p_x,p_B,p_M,ϕ_B_loc,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,
       a_min,a_max,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
       coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,
       labor_allocated_interior_c3a_fine,λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,
       x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,q_S_c2_fine,
       q_B_c1_fine,q_S_c3a_fine,q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,
       x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,Y_B_c3b_fine,
       c_S_mat_fine,c_B_mat_fine,c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,q_B_mat_fine,land_B_mat_fine,
       λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,feasibility_mat_fine,C_max_mat_fine,C_min_mat_fine,
       q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,
       x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat_fine,C_max_staple_fine,
       C_min_staple_fine,C_max_staple_constrained_fine,C_min_staple_constrained_fine,TC_S_c3_constrained_fine,x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine,
       x_S_mat_3c_fine,x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine);
       Q_trans = spzeros(ns_fine*3,ns_fine*3);
       (Q_trans_prime,cons_fine_local,a_prime_fine_local,future_occupation_fine_local,V_saved_local) = Q_transition(coeff,ns_fine,P_kron_fine,min_C_applied_fine,max_C_applied_fine,Phi_z_fine,
           β,fspace_a,σ,P_W_fine,Y_W_fine,s_fine,r,ρ,w,coeff_λ_2_cashcrop_residual_unconstrained_fine,coeff_λ_2_cashcrop_residual_constrained_fine,θ_fine,
           fspace_C_fine,ϕ_S_loc,ζ,τ_S,p_x,p_B,p_M,ϕ_B_loc,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,κ,tol,a_min,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
           coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,labor_allocated_interior_c3a_fine,λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,
           x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,q_S_c2_fine,q_B_c1_fine,q_S_c3a_fine,q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,
           x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,Y_B_c3b_fine,c_S_mat_fine,c_B_mat_fine,c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,
           q_B_mat_fine,land_B_mat_fine,λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,feasibility_mat_fine,C_max_mat_fine,C_min_mat_fine,q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,
           c_M_staples_fine,P_S_staples_fine,x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat_fine,a_max,P_kron1,Q_trans,
           C_max_staple_fine,C_min_staple_fine,C_max_staple_constrained_fine,
           C_min_staple_constrained_fine,TC_S_c3_constrained_fine,x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine,fspace_a_fine,
           x_S_mat_3c_fine,x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine);

       distr_current = Q_trans_prime*distr_previous;
        # Start acquiring objects for the equilibrium
        # Distribution of past occupations
        worker_past_dist = distr_previous[(ns_fine *0 + 1):(ns_fine *1)];
        staple_past_dist = distr_previous[(ns_fine *1 + 1):(ns_fine *2)];
        cash_crop_past_dist = distr_previous[(ns_fine *2 + 1):(ns_fine *3)];

        # Policy functions:
        c_S_W_fine= zeros(ns_fine,3);
        c_B_W_fine= zeros(ns_fine,3);
        c_M_W_fine= zeros(ns_fine,3);
        Y_W_fine_policy = zeros(ns_fine,3);
        Y_S_fine = zeros(ns_fine,3);
        c_S_S_fine= zeros(ns_fine,3);
        c_B_S_fine= zeros(ns_fine,3);
        c_M_S_fine= zeros(ns_fine,3);
        q_S_S_fine= zeros(ns_fine,3);
        P_S_fine= zeros(ns_fine,3);
        x_S_S_fine= zeros(ns_fine,3);
        solve_staple_index_S_fine= zeros(Int64,ns_fine,3);
        λ_2_S_fine= zeros(ns_fine,3);
        future_asset_S= zeros(ns_fine,3);
        future_asset_C= zeros(ns_fine,3);
        c_S_B_fine= zeros(ns_fine,3);
        c_B_B_fine= zeros(ns_fine,3);
        c_M_B_fine= zeros(ns_fine,3);
        x_SC_fine= zeros(ns_fine,3);
        x_BC_fine= zeros(ns_fine,3);
        land_C_fine= zeros(ns_fine,3);
        λ_2_fine= zeros(ns_fine,3);
        P_B_fine= zeros(ns_fine,3);
        Y_B_fine= zeros(ns_fine,3);
        q_S_B_fine= zeros(ns_fine,3);
        q_B_B_fine= zeros(ns_fine,3);
        solve_cash_crop_index_B_fine= zeros(Int64,ns_fine,3);
        solve_staple_index_B_fine= zeros(Int64,ns_fine,3);
        TC_fine= zeros(ns_fine,3);
        for j = 1:3
            for jj =1:3
            if jj == 1
                Y_W_fine_policy[:,j] =  P_W_fine *cons_fine_local[:,j]  + Y_W_fine  .+ w*FM_W .+ w*F_W * (j != 1);
            end
                if jj == 2
                    (future_asset_S[:,j],Y_S_fine[:,j],c_S_S_fine[:,j],c_B_S_fine[:,j],c_M_S_fine[:,j],q_S_S_fine[:,j],
                    P_S_fine[:,j],x_S_S_fine[:,j],solve_staple_index_S_fine[:,j],λ_2_S_fine[:,j]) = policy_function_creator(
                    cons_fine_local[:,j],jj,j,P_W_fine,Y_W_fine,s_fine,r,ρ,w,coeff_λ_2_cashcrop_residual_unconstrained_fine,
                    coeff_λ_2_cashcrop_residual_constrained_fine,θ_fine,fspace_C_fine,ϕ_S_loc,ζ,τ_S,p_x,p_B,
                    p_M,ϕ_B_loc,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,a_min,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
                    coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,
                    labor_allocated_interior_c3a_fine,λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,
                    q_S_c2_fine,q_B_c1_fine,q_S_c3a_fine,q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,
                    Y_B_c3b_fine,c_S_mat_fine,c_B_mat_fine,c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,q_B_mat_fine,land_B_mat_fine,λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,
                    feasibility_mat_fine,C_max_mat_fine,C_min_mat_fine,q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,x_S_staples_fine,
                    λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat_fine,C_max_staple_fine,C_min_staple_fine,C_max_staple_constrained_fine,
                    C_min_staple_constrained_fine,TC_S_c3_constrained_fine,x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine,
                    x_S_mat_3c_fine,x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine);
                elseif jj == 3
                    (future_asset_C[:,j],c_S_B_fine[:,j],c_B_B_fine[:,j],c_M_B_fine[:,j],x_SC_fine[:,j],x_BC_fine[:,j],
                    land_C_fine[:,j],λ_2_fine[:,j],P_B_fine[:,j],Y_B_fine[:,j],q_S_B_fine[:,j],q_B_B_fine[:,j],solve_cash_crop_index_B_fine[:,j]
                    ,solve_staple_index_B_fine[:,j],TC_fine[:,j]) = policy_function_creator(
                    cons_fine_local[:,j],jj,j,P_W_fine,Y_W_fine,s_fine,r,ρ,w,coeff_λ_2_cashcrop_residual_unconstrained_fine,
                    coeff_λ_2_cashcrop_residual_constrained_fine,θ_fine,fspace_C_fine,ϕ_S_loc,ζ,τ_S,p_x,p_B,
                    p_M,ϕ_B_loc,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,a_min,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
                    coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,
                    labor_allocated_interior_c3a_fine,λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,
                    q_S_c2_fine,q_B_c1_fine,q_S_c3a_fine,q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,
                    Y_B_c3b_fine,c_S_mat_fine,c_B_mat_fine,c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,q_B_mat_fine,land_B_mat_fine,λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,
                    feasibility_mat_fine,C_max_mat_fine,C_min_mat_fine,q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,x_S_staples_fine,
                    λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat_fine,C_max_staple_fine,C_min_staple_fine,C_max_staple_constrained_fine,
                    C_min_staple_constrained_fine,TC_S_c3_constrained_fine,x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine,
                    x_S_mat_3c_fine,x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine);
                end
            end
        end
        matcheck= zeros(ns_fine,12)
        matcheck[:,1:3]= future_occupation_fine_local;
        matcheck[:,4:6] = future_asset_S;#s_fine
        matcheck[:,7:9] = future_asset_C;#s_fine
        #matcheck[:,6:8] = cons_fine_local
        #matcheck[:,9:11] = solve_cash_crop_index_B_fine
        # Distribution of current occupations
        current_distr_store = zeros(3 * ns_fine);

        past_distr_store = zeros(ns_fine,3);

        stay_workers = worker_past_dist.*(future_occupation_fine_local[:,1].==1);
        exit_staple_to_work = staple_past_dist.*(future_occupation_fine_local[:,2].==1);
        exit_cashcrop_to_work = cash_crop_past_dist.*(future_occupation_fine_local[:,3].==1);
        current_workers = stay_workers + exit_staple_to_work + exit_cashcrop_to_work;
        # Worker policy functions are easy:
        c_S_W_fine[:,1] = (c̄_S .+ (1 + Q_S)^(-ϵ)*ψ_S.^ϵ.*cons_fine_local[:,1].*P_W_fine.^ϵ);
        c_S_W_fine[:,2] = (c̄_S .+ (1 + Q_S)^(-ϵ)*ψ_S.^ϵ.*cons_fine_local[:,2].*P_W_fine.^ϵ);
        c_S_W_fine[:,3] = (c̄_S .+ (1 + Q_S)^(-ϵ)*ψ_S.^ϵ.*cons_fine_local[:,3].*P_W_fine.^ϵ) ;
        c_B_W_fine[:,1] = (p_B^(-ϵ)*ψ_B.^ϵ.*cons_fine_local[:,1].*P_W_fine.^ϵ);
        c_B_W_fine[:,2] = (p_B^(-ϵ)*ψ_B.^ϵ.*cons_fine_local[:,2].*P_W_fine.^ϵ);
        c_B_W_fine[:,3] = (p_B^(-ϵ)*ψ_B.^ϵ.*cons_fine_local[:,3].*P_W_fine.^ϵ) ;
        c_M_W_fine[:,1] = (p_M^(-ϵ)*ψ_M.^ϵ.*cons_fine_local[:,1].*P_W_fine.^ϵ);
        c_M_W_fine[:,2] = (p_M^(-ϵ)*ψ_M.^ϵ.*cons_fine_local[:,2].*P_W_fine.^ϵ);
        c_M_W_fine[:,3] = (p_M^(-ϵ)*ψ_M.^ϵ.*cons_fine_local[:,3].*P_W_fine.^ϵ) ;

        c_S_worker_sum = sum( c_S_W_fine[:,1].*stay_workers + c_S_W_fine[:,2].*exit_staple_to_work +
        c_S_W_fine[:,3] .*exit_cashcrop_to_work );
        c_B_worker_sum = sum(c_B_W_fine[:,1] .*stay_workers +c_B_W_fine[:,2] .*exit_staple_to_work +
        c_B_W_fine[:,3] .*exit_cashcrop_to_work );
        c_M_worker_sum = sum(c_M_W_fine[:,1] .*stay_workers +c_M_W_fine[:,2] .*exit_staple_to_work +
        c_M_W_fine[:,3] .*exit_cashcrop_to_work );
        urban_labor_supply_sum = sum(current_workers.*labor_prod_fine);
        transaction_cost_worker_sum = c_S_worker_sum * Q_S;

        entrants_staple_from_workers = worker_past_dist.*(future_occupation_fine_local[:,1].==2);
        incumbents_staple = staple_past_dist.*(future_occupation_fine_local[:,2].==2);
        exit_cashcrop_to_staple = cash_crop_past_dist.*(future_occupation_fine_local[:,3].==2);
        current_staple = entrants_staple_from_workers + incumbents_staple + exit_cashcrop_to_staple;
        # Staple policy function now depends on consumption, hence the decision is always different!


        c_S_staple_sum = sum(c_S_S_fine[:,1] .*entrants_staple_from_workers + c_S_S_fine[:,2] .*incumbents_staple
        + c_S_S_fine[:,3] .*exit_cashcrop_to_staple);
        c_B_staple_sum = sum(c_B_S_fine[:,1] .*entrants_staple_from_workers + c_B_S_fine[:,2] .*incumbents_staple
        + c_B_S_fine[:,3] .*exit_cashcrop_to_staple);
        c_M_staple_sum = sum(c_M_S_fine[:,1] .*entrants_staple_from_workers + c_M_S_fine[:,2] .*incumbents_staple
        + c_M_S_fine[:,3] .*exit_cashcrop_to_staple);
        q_S_staple_sum = sum(q_S_S_fine[:,1] .*entrants_staple_from_workers + q_S_S_fine[:,2] .*incumbents_staple
        + q_S_S_fine[:,3] .*exit_cashcrop_to_staple);
        x_S_staple_sum = sum(x_S_S_fine[:,1] .*entrants_staple_from_workers + x_S_S_fine[:,2] .*incumbents_staple
        + x_S_S_fine[:,3] .*exit_cashcrop_to_staple);

        transaction_cost_staple_sum = Q_S * sum(max.(c_S_S_fine[:,1] - q_S_S_fine[:,1],0) .*entrants_staple_from_workers +
            max.(c_S_S_fine[:,2] - q_S_S_fine[:,2],0) .*incumbents_staple +
            max.(c_S_S_fine[:,3] - q_S_S_fine[:,3],0) .*exit_cashcrop_to_staple);

        entrants_cashcrop_from_workers = worker_past_dist.*(future_occupation_fine_local[:,1].==3);
        entrants_from_staple_to_cashcrop= staple_past_dist.*(future_occupation_fine_local[:,2].==3);
        incumbents_cashcrop = cash_crop_past_dist.*(future_occupation_fine_local[:,3].==3);
        current_cashcrop = entrants_cashcrop_from_workers + incumbents_cashcrop + entrants_from_staple_to_cashcrop;
        # Cash crop producer policy functions now depends on consumption, hence the decision is always different!

        c_S_cashcrop_sum = sum(c_S_B_fine[:,1] .*entrants_cashcrop_from_workers + c_S_B_fine[:,3] .*incumbents_cashcrop
        + c_S_B_fine[:,2] .*entrants_from_staple_to_cashcrop);
        c_B_Cashcrop_sum = sum(c_B_B_fine[:,1] .*entrants_cashcrop_from_workers + c_B_B_fine[:,3] .*incumbents_cashcrop
        + c_B_B_fine[:,2] .*entrants_from_staple_to_cashcrop);
        c_M_cashcrop_sum = sum(c_M_B_fine[:,1] .*entrants_cashcrop_from_workers + c_M_B_fine[:,3] .*incumbents_cashcrop
        + c_M_B_fine[:,2] .*entrants_from_staple_to_cashcrop);
        q_S_cashcrop_sum = sum(q_S_B_fine[:,1] .*entrants_cashcrop_from_workers + q_S_B_fine[:,3] .*incumbents_cashcrop
        + q_S_B_fine[:,2] .*entrants_from_staple_to_cashcrop);
        q_B_cashcrop_sum = sum(q_B_B_fine[:,1] .*entrants_cashcrop_from_workers + q_B_B_fine[:,3] .*incumbents_cashcrop
        + q_B_B_fine[:,2] .*entrants_from_staple_to_cashcrop);
        x_S_cashcrop_sum = sum(x_SC_fine[:,1] .*entrants_cashcrop_from_workers + x_SC_fine[:,3] .*incumbents_cashcrop
        + x_SC_fine[:,2] .*entrants_from_staple_to_cashcrop);
        x_B_cashcrop_sum = sum(x_BC_fine[:,1] .*entrants_cashcrop_from_workers + x_BC_fine[:,3] .*incumbents_cashcrop
        + x_BC_fine[:,2] .*entrants_from_staple_to_cashcrop);

        transaction_cost_cashcrop_sum = Q_S * sum(max.(c_S_B_fine[:,1] - q_S_B_fine[:,1],0) .*entrants_cashcrop_from_workers +
            max.(c_S_B_fine[:,2] - q_S_B_fine[:,2],0) .*entrants_from_staple_to_cashcrop +
            max.(c_S_B_fine[:,3] - q_S_B_fine[:,3],0) .*incumbents_cashcrop);

        current_distr_store[(ns_fine *0 + 1):(ns_fine *1)] = current_workers;
        current_distr_store[(ns_fine *1 + 1):(ns_fine *2)] = current_staple;
        current_distr_store[(ns_fine *2 + 1):(ns_fine *3)] = current_cashcrop;


        past_distr_store[:,1]=distr_previous[(ns_fine *0 + 1):(ns_fine *1)];
        past_distr_store[:,2]=distr_previous[(ns_fine *1 + 1):(ns_fine *2)];
        past_distr_store[:,3]=distr_previous[(ns_fine *2 + 1):(ns_fine *3)];

        current_worker_pop = sum(current_workers);
        current_staple_pop = sum(current_staple);
        current_cashcrop_pop = sum(current_cashcrop);
        # Entry cost accounting - measured in labor, convert with wages later
        entry_costs_to_workers =  F_W*(sum(exit_staple_to_work) + sum(exit_cashcrop_to_work));
        entry_costs_to_staples = F_S*(sum(entrants_staple_from_workers)); #Not needed: + sum(exit_cashcrop_to_work));
        entry_costs_to_cashcrops = F_B*(sum(entrants_cashcrop_from_workers) + sum(entrants_from_staple_to_cashcrop));
        total_entry_cost = entry_costs_to_workers + entry_costs_to_staples + entry_costs_to_cashcrops;
        # Maintenance cost accounting
        maintenance_costs_for_workers =  FM_W*(current_worker_pop);
        maintenance_costs_to_staples = FM_S*(current_staple_pop); #Not needed: + sum(exit_cashcrop_to_work));
        maintenance_costs_to_cashcrops = FM_B*(current_cashcrop_pop);
        total_maintenance_cost = maintenance_costs_for_workers + maintenance_costs_to_staples + maintenance_costs_to_cashcrops;

        # Savings and debt
        worker_bond_holding = s_fine[:,1] .*current_workers;
        worker_bond_holding_sum = sum(worker_bond_holding);
        staple_bond_holding = s_fine[:,1] .*current_staple;
        staple_bond_holding_sum = sum(staple_bond_holding[staple_bond_holding .> 0.0]);
        staple_debt_holding_sum = sum(staple_bond_holding[staple_bond_holding .< 0.0]);
        cashcrop_bond_holding = s_fine[:,1] .*current_cashcrop;
        cashcrop_bond_holding_sum = sum(cashcrop_bond_holding[cashcrop_bond_holding .> 0.0]);
        cashcrop_debt_holding_sum = sum(cashcrop_bond_holding[cashcrop_bond_holding .< 0.0]);
        asset_supply = worker_bond_holding_sum + (staple_bond_holding_sum + staple_debt_holding_sum ) + (cashcrop_bond_holding_sum + cashcrop_debt_holding_sum);


        # Manufacturing firm:
        labor_used = urban_labor_supply_sum - (total_maintenance_cost + total_entry_cost);

        # Govt spending side as only the manufacturing firm pays taxes --- # KAROL FIXED THIS
        input_staple= x_S_staple_sum + x_S_cashcrop_sum; # change current_staple to current_staple2 ? etc.
        input_cashcrop= copy(x_B_cashcrop_sum);
        foreign_demand_cash_crop = a_D*p_B^b_D#foreign_demand/p_B;
        #println(input_staple)
        #println(input_cashcrop)
        
        # Alternatively, govt borrowing would make a ton of sense here...
        Government_expenditure = balanced_share * p_x * (τ_S * input_staple + τ_B * input_cashcrop);
        Import_value = p_x * ( input_staple + input_cashcrop);
        Export_value = p_B * foreign_demand_cash_crop;
        
        #current_account = (Government_expenditure + Net_Foreign_Factor_Income - Import_value+Export_value)
        current_account_residual = Export_value - Import_value;
        #Net_Foreign_Factor_Income = -current_account_residual;
        #foreign_supply_capital = -Net_Foreign_Factor_Income/R;
        capital_used = copy(asset_supply) + foreign_supply_capital;
        Net_Foreign_Factor_Income = R*(copy(asset_supply) - capital_used)

        # if capital_used<0 || labor_used<0
        #     println("ERROR? Capital_used: ", capital_used, " Labor_used: ", labor_used)
        # end
        #
        # KLratio=α/(1-α)*((1+τ_W)*w)/r
        # r_new = -δ + α*p_M * KLratio^(α-1) #(max(0.01,capital_used)/max(0.01,labor_used))^(α-1)
        # w_new = (1 - α)/(1 + τ_W)*p_M* KLratio^α #(max(0.01,capital_used)/max(0.01,labor_used))^α # this is when we calculate wages paid by the manufacturing firm
        τ_W_new = - balanced_share * p_x * (τ_S * input_staple + τ_B * input_cashcrop) / (labor_used * w);
        K_L_ratio = maximum([capital_used./labor_used, tol]);
        r_new = -δ + α*p_M *K_L_ratio^(α-1);
        w_new = (1 - α)/(1 + τ_W)*p_M* K_L_ratio^α; # this is when we calculate wages paid by the manufacturing firm
        #labor_demand=((1 - α)/(1 + τ_W)/w_new*p_M*capital_used^α)^(1/α);

        prod_staple = q_S_staple_sum + q_S_cashcrop_sum;
        prod_cashcrop = q_B_cashcrop_sum;
        prod_manuf = sum(max(tol,capital_used).^α.*max(tol,labor_used).^(1-α));
        #foreign_demand = p_x * ( input_staple + input_cashcrop);
        
        #p_M_new =
        #τ_W*(labor_used * w)

                # Production -- KAROL ADDED THIS
        residual_goods = zeros(6);
        #residual[1] = (cash_crop_cons + p_x/p_B*input_staple + p_x/p_B*input_cashcrop - prod_cashcrop)/(cash_crop_cons + p_x/p_B*input_staple + p_x/p_B*input_cashcrop + prod_cashcrop); # -- KAROL ADDED THIS
        #residual[1] = 0#(cash_crop_cons + p_x*input_staple + p_x*input_cashcrop - prod_cashcrop)/(cash_crop_cons + p_x*input_staple + p_x*input_cashcrop + prod_cashcrop); # -- KAROL ADDED THIS
        #residual[1] = (cash_crop_cons  - prod_cashcrop)/(cash_crop_cons + prod_cashcrop); # -- KAROL ADDED THIS
        residual_goods[1] = (c_B_worker_sum + c_B_staple_sum + c_B_Cashcrop_sum+ foreign_demand_cash_crop  - prod_cashcrop)/(
                            c_B_worker_sum + c_B_staple_sum + c_B_Cashcrop_sum+ foreign_demand_cash_crop  + prod_cashcrop);

        #println(c_B_worker_sum + c_B_staple_sum + c_B_Cashcrop_sum)
        #println(foreign_demand_cash_crop)
        #println(prod_cashcrop)

        residual_goods[2] = (c_M_worker_sum + c_M_staple_sum + c_M_cashcrop_sum - prod_manuf)/(c_M_worker_sum + c_M_staple_sum + c_M_cashcrop_sum + prod_manuf);

        #residual[2] = (manuf_cons + input_staple + input_cashcrop - prod_manuf)/(manuf_cons + input_staple + input_cashcrop + prod_manuf); # -- KAROL ADDED THIS
        residual_goods[3] = 0 #No longer needed for the calibration as the foreign_supply_capital offsets capital markets(r - r_new)/(r+r_new);


        #residual[4] = (labor_demand - labor_used)/(labor_demand + labor_used);
        #staple_mkt_clr = (sum(c_current_S) + p_x*input_staple - prod_staple)/(sum(c_current_S) + prod_staple + p_x*input_staple) # -- KAROL ADDED THIS
        #staple_mkt_clr = (sum(c_current_S) - prod_staple)/(sum(c_current_S) + prod_staple) # -- KAROL ADDED THIS
        residual_goods[5] = (c_S_worker_sum + c_S_staple_sum + c_S_cashcrop_sum +
            transaction_cost_worker_sum + transaction_cost_staple_sum + transaction_cost_cashcrop_sum - prod_staple)/(
            c_S_worker_sum + c_S_staple_sum + c_S_cashcrop_sum +
                transaction_cost_worker_sum + transaction_cost_staple_sum + transaction_cost_cashcrop_sum + prod_staple);# This is the main residual, since other markets are tempered with
        if balanced_share >0.0
            residual_goods[6] = (τ_W - τ_W_new)/(τ_W + τ_W_new)
        else
            residual_goods[6] = 0.0
        end
        undernutritioned_workers = sum((c_S_W_fine[:, 1] .< cons_level_substinence) .* stay_workers
        + (c_S_W_fine[:, 2] .< cons_level_substinence) .* exit_staple_to_work +
        (c_S_W_fine[:, 3] .< cons_level_substinence) .* exit_cashcrop_to_work)
        undernutritioned_staple_farmer = sum((c_S_S_fine[:, 1] .< cons_level_substinence) .* entrants_staple_from_workers +
                    (c_S_S_fine[:, 2] .< cons_level_substinence) .* incumbents_staple +
                    (c_S_S_fine[:, 3] .< cons_level_substinence) .* exit_cashcrop_to_staple)
        undernutritioned_cashcrop_farmer = sum((c_S_B_fine[:, 1] .< cons_level_substinence) .* entrants_cashcrop_from_workers
                        + (c_S_B_fine[:, 2] .< cons_level_substinence) .* entrants_from_staple_to_cashcrop +
                        (c_S_B_fine[:, 3] .< cons_level_substinence) .* incumbents_cashcrop)

        cttilde_new = undernutritioned_workers + undernutritioned_staple_farmer + undernutritioned_cashcrop_farmer
        residual_goods[4] = (cttilde_new - cttilde)/(cttilde_new+cttilde);#(w - w_new)/(w+w_new);

        # if current_staple_pop ==0
        #     residual[5] = 100.0;
        # end

        #rural_pop_only_staples_printable=(sum(entrants_cashcrop_from_workers.* (q_B_B_fine[:,1].==0.0)) + sum(incumbents_cashcrop.*(q_B_B_fine[:,3].==0.0))
        #+ sum( entrants_from_staple_to_cashcrop.*(q_B_B_fine[:,2].==0.0)) + current_staple_pop);
        #rural_pop_only_Bashcrop_printable=(sum(entrants_cashcrop_from_workers.* (land_C_fine[:,1].==1.0)) + sum(incumbents_cashcrop.*(land_C_fine[:,3].==1.0))
        #+ sum( entrants_from_staple_to_cashcrop.*(land_C_fine[:,2].==1.0)));
        #println("Pops: Worker:",current_worker_pop," Staple:",current_staple_pop," Cash crop:",current_cashcrop_pop," Only producing staples", rural_pop_only_staples_printable
        #," Only producing cashcrop", rural_pop_only_Bashcrop_printable)
    return residual_goods,distr_current
end

function Residual_transition_backward_epsilon(s::Array{Float64,2},ns::Int64,
    z::Array{Float64,1},z_W::Array{Float64,1},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,
    p_B::Float64,p_M::Float64,ϕ_B::Float64,τ_B::Float64,ρ::Float64,w::Float64,r::Float64,
    c̄_S::Float64,a_min::Float64,a_max::Float64,γ::Float64,n::Array{Int64,1},κ::Float64,Q_S::Float64,
    ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,agrid_fine::Array{Float64,1},
    fspace_C_fine::Dict{Symbol,Any},agrid::Array{Float64,1},
    coefficients_next_tmp::Array{Float64,2},coeff_next::Array{Float64,2},C_grid_fine::Array{Float64,1},
    F_W::Float64,F_S::Float64,F_B::Float64,FM_W::Float64,FM_S::Float64,FM_B::Float64,P_kron::SparseMatrixCSC{Float64, Int64},Phi::SparseMatrixCSC{Float64, Int64},
    Phi_z::SparseMatrixCSC{Float64, Int64},β::Float64,
    fspace_a::Dict{Symbol,Any},σ::Float64,cons_level_substinence::Float64, epsilon_u::Float64, epsilon_r::Float64,cttilde::Float64,tol::Float64 = 1e-8)

    ctilde = copy(cons_level_substinence)
    ϕ_S_loc = ϕ_S * exp(-epsilon_r * (cttilde - 0.2))
    ϕ_B_loc = ϕ_B * exp(-epsilon_r * (cttilde - 0.2))
    (θ,labor_prod,tol,P_W,Y_W,coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
    coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,λ_B_interior_c3a,x_SC_interior_c3a,x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,
    q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b, c_S_mat,c_B_mat,
    c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat, λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,q_S_staples,c_S_staples,c_B_staples,
    c_M_staples,P_S_staples,x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,C_max_unconstrained ,C_max_constrained,C_min_unconstrained,
    C_min_constrained,TC_mat,C_max_staple,C_min_staple,C_max_staple_constrained,C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,
    q_S_c3_constrained,c_S_c3_constrained,cbar_violated, x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c) = income_creator_ext(s,ns,z,z_W,ϕ_S_loc,ζ,τ_S,p_x,p_B,p_M,ϕ_B_loc,τ_B,ρ,w,r,
        c̄_S,a_min,a_max,γ,n,κ,Q_S,ϵ,ψ_S,ψ_B,ψ_M,C_grid_fine,fspace_C_fine,agrid, epsilon_u, cttilde, ctilde);
    #println("cbar: ", cbar_violated)


    min_C_applied,max_C_applied = bounds_consumption(P_W,Y_W,s,r,ρ,w,
        coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,θ,
        fspace_C_fine,ϕ_S_loc,ζ,τ_S,p_x,p_B,p_M,ϕ_B_loc,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns,κ,tol,a_min,a_max,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
        coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,
        λ_B_interior_c3a,x_SC_interior_c3a,x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,
        x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b, c_S_mat,c_B_mat,
        c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat, λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,q_S_staples,c_S_staples,c_B_staples,
        c_M_staples,P_S_staples,x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,F_W,F_S,F_B,
        FM_W,FM_S,FM_B,TC_mat,C_max_staple,C_min_staple,C_max_staple_constrained,C_min_staple_constrained,TC_S_c3_constrained,
        x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained, x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c);

        x_tmp = zeros(ns,3);
        V_tmp = zeros(ns,3);
        conv = 10.0

    (coeff, conv, conv_ind)  = Bellman_iteration(coefficients_next_tmp,coeff_next,x_tmp,V_tmp,P_kron,Phi,min_C_applied,max_C_applied,Phi_z,β,
        fspace_a,σ,P_W,Y_W,s,r,ρ,w,coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,θ,
        fspace_C_fine,ϕ_S_loc,ζ,τ_S,p_x,p_B,p_M,ϕ_B_loc,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns,κ,tol,a_min,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
        coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,
        λ_B_interior_c3a,x_SC_interior_c3a,x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,
        x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b, c_S_mat,c_B_mat,
        c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat, λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,q_S_staples,c_S_staples,c_B_staples,
        c_M_staples,P_S_staples,x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat,a_max,C_max_staple,C_min_staple,C_max_staple_constrained,
        C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained,x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c);
    return (coeff,coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,C_max_unconstrained,C_max_constrained,
    C_min_unconstrained,C_min_constrained,coeff_λ_2_s,C_max_staple,C_min_staple,C_max_staple_constrained,C_min_staple_constrained,
    TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained)
end


function Residual_transition_forward_epsilon_detailed(coeff::Array{Float64,2},distr_previous::Array{Float64,1},s_fine::Array{Float64,2},ns_fine::Int64,
    z::Array{Float64,1},z_W::Array{Float64,1},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,
    p_B::Float64,p_M::Float64,ϕ_B::Float64,τ_B::Float64,ρ::Float64,w::Float64,r::Float64,
    c̄_S::Float64,a_min::Float64,a_max::Float64,γ::Float64,n_fine::Array{Int64,1},κ::Float64,Q_S::Float64,
    ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,
    coeff_λ_2_cashcrop_residual_unconstrained::Array{Float64,2},
    coeff_λ_2_cashcrop_residual_constrained::Array{Float64,2},
    C_max_unconstrained::Array{Float64,1} ,C_max_constrained::Array{Float64,1},C_min_unconstrained::Array{Float64,1},
    C_min_constrained::Array{Float64,1},
    coeff_λ_2_s::Array{Float64,2},agrid_fine::Array{Float64,1},fspace_C_fine::Dict{Symbol,Any},C_max_staple::Array{Float64,1},
    C_min_staple::Array{Float64,1},C_max_staple_constrained::Array{Float64,1},
    C_min_staple_constrained::Array{Float64,1},TC_S_c3_constrained::Array{Float64,1},
    x_S_c3_constrained::Array{Float64,1},q_S_c3_constrained::Array{Float64,1},c_S_c3_constrained::Array{Float64,1},foreign_supply_capital::Float64,
    balanced_share::Float64,τ_W::Float64,C_grid_fine::Array{Float64,1},F_W::Float64,F_S::Float64,F_B::Float64,FM_W::Float64,FM_S::Float64,FM_B::Float64,
    P_kron_fine::SparseMatrixCSC{Float64, Int64},
    Phi_z_fine::SparseMatrixCSC{Float64, Int64},β::Float64,
    fspace_a::Dict{Symbol,Any},σ::Float64,P_kron1::SparseMatrixCSC{Float64, Int64},fspace_a_fine::Dict{Symbol,Any},a_D::Float64,b_D::Float64,
    R::Float64,δ::Float64,α::Float64,cons_level_substinence::Float64, epsilon_u::Float64, epsilon_r::Float64,cttilde::Float64,tol::Float64 = 1e-8)
   
    # coeff = coeff_store[:,:,t]; 
    # distr_previous = distr_store[:,t]; 
    # τ_S = τ_S_loc;
    # coeff_λ_2_cashcrop_residual_unconstrained = coeff_λ_2_cashcrop_residual_unconstrained_store[:,:,t-1]
    # coeff_λ_2_cashcrop_residual_constrained = coeff_λ_2_cashcrop_residual_constrained_store[:,:,t-1]
    # C_max_unconstrained = C_max_unconstrained_store[:,t-1]
    # C_max_constrained = C_max_constrained_store[:,t-1]
    # C_min_unconstrained = C_min_unconstrained_store[:,t-1]
    # C_min_constrained = C_min_constrained_store[:,t-1]
    # coeff_λ_2_s = coeff_λ_2_s_store[:,:,t-1]
    # C_max_staple = C_max_staple_store[:,t-1]
    # C_min_staple = C_min_staple_store[:,t-1]
    # C_max_staple_constrained = C_max_staple_constrained_store[:,t-1]
    # C_min_staple_constrained = C_min_staple_constrained_store[:,t-1]
    # TC_S_c3_constrained = TC_S_c3_constrained_store[:,t-1]
    # x_S_c3_constrained = x_S_c3_constrained_store[:,t-1]
    # q_S_c3_constrained = q_S_c3_constrained_store[:,t-1]
    # c_S_c3_constrained = c_S_c3_constrained_store[:,t-1]
    # foreign_supply_capital = capital_trans[t]
    
    ctilde = copy(cons_level_substinence)
    ϕ_S_loc = ϕ_S * exp(-epsilon_r * (cttilde - 0.2))
    ϕ_B_loc = ϕ_B * exp(-epsilon_r * (cttilde - 0.2))
    # Aggregates
        (θ_fine,labor_prod_fine,tol,P_W_fine,Y_W_fine,coeff_λ_2_cashcrop_residual_unconstrained_fine,coeff_λ_2_cashcrop_residual_constrained_fine,
        x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
        coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,labor_allocated_interior_c3a_fine,
        λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,q_S_c2_fine,q_B_c1_fine,q_S_c3a_fine,
        q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,Y_B_c3b_fine, c_S_mat_fine,c_B_mat_fine,
        c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,q_B_mat_fine,land_B_mat_fine, λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,feasibility_mat_fine,C_max_mat_fine,
        C_min_mat_fine,q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,
        Y_S_potential_fine,TC_mat_fine,C_max_staple_fine,C_min_staple_fine,C_max_staple_constrained_fine,
        C_min_staple_constrained_fine,TC_S_c3_constrained_fine,x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine, x_S_mat_3c_fine,x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine) =  income_creator_no_approx_ext(s_fine,ns_fine,
        z,z_W,ϕ_S_loc,ζ,τ_S,p_x,p_B,p_M,ϕ_B_loc,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n_fine,κ,Q_S,ϵ,ψ_S,ψ_B,ψ_M,coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,
        C_max_unconstrained ,C_max_constrained,C_min_unconstrained,C_min_constrained, coeff_λ_2_s,C_grid_fine,fspace_C_fine,C_max_staple,C_min_staple,C_max_staple_constrained,
    C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained, epsilon_u, cttilde, ctilde);


        min_C_applied_fine,max_C_applied_fine = bounds_consumption(P_W_fine,Y_W_fine,s_fine,r,ρ,w,
        coeff_λ_2_cashcrop_residual_unconstrained_fine,coeff_λ_2_cashcrop_residual_constrained_fine,θ_fine,
        fspace_C_fine,ϕ_S_loc,ζ,τ_S,p_x,p_B,p_M,ϕ_B_loc,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,
        a_min,a_max,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
        coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,
        labor_allocated_interior_c3a_fine,λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,
        x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,q_S_c2_fine,
        q_B_c1_fine,q_S_c3a_fine,q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,
        x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,Y_B_c3b_fine,
        c_S_mat_fine,c_B_mat_fine,c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,q_B_mat_fine,land_B_mat_fine,
        λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,feasibility_mat_fine,C_max_mat_fine,C_min_mat_fine,
        q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,
        x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat_fine,C_max_staple_fine,
        C_min_staple_fine,C_max_staple_constrained_fine,C_min_staple_constrained_fine,TC_S_c3_constrained_fine,x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine,
        x_S_mat_3c_fine,x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine);
        Q_trans = spzeros(ns_fine*3,ns_fine*3);
        (Q_trans_prime,cons_fine_local,a_prime_fine_local,future_occupation_fine_local,V_saved_local) = Q_transition(coeff,ns_fine,P_kron_fine,min_C_applied_fine,max_C_applied_fine,Phi_z_fine,
            β,fspace_a,σ,P_W_fine,Y_W_fine,s_fine,r,ρ,w,coeff_λ_2_cashcrop_residual_unconstrained_fine,coeff_λ_2_cashcrop_residual_constrained_fine,θ_fine,
            fspace_C_fine,ϕ_S_loc,ζ,τ_S,p_x,p_B,p_M,ϕ_B_loc,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,κ,tol,a_min,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
            coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,labor_allocated_interior_c3a_fine,λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,
            x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,q_S_c2_fine,q_B_c1_fine,q_S_c3a_fine,q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,
            x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,Y_B_c3b_fine,c_S_mat_fine,c_B_mat_fine,c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,
            q_B_mat_fine,land_B_mat_fine,λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,feasibility_mat_fine,C_max_mat_fine,C_min_mat_fine,q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,
            c_M_staples_fine,P_S_staples_fine,x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat_fine,a_max,P_kron1,Q_trans,
            C_max_staple_fine,C_min_staple_fine,C_max_staple_constrained_fine,
            C_min_staple_constrained_fine,TC_S_c3_constrained_fine,x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine,fspace_a_fine,
            x_S_mat_3c_fine,x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine);

        distr_current = Q_trans_prime*distr_previous;
        # Start acquiring objects for the equilibrium
        # Distribution of past occupations
        worker_past_dist = distr_previous[(ns_fine *0 + 1):(ns_fine *1)];
        staple_past_dist = distr_previous[(ns_fine *1 + 1):(ns_fine *2)];
        cash_crop_past_dist = distr_previous[(ns_fine *2 + 1):(ns_fine *3)];

        # Policy functions:
        c_S_W_fine= zeros(ns_fine,3);
        c_B_W_fine= zeros(ns_fine,3);
        c_M_W_fine= zeros(ns_fine,3);
        Y_W_fine_policy = zeros(ns_fine,3);
        Y_S_fine = zeros(ns_fine,3);
        c_S_S_fine= zeros(ns_fine,3);
        c_B_S_fine= zeros(ns_fine,3);
        c_M_S_fine= zeros(ns_fine,3);
        q_S_S_fine= zeros(ns_fine,3);
        P_S_fine= zeros(ns_fine,3);
        x_S_S_fine= zeros(ns_fine,3);
        solve_staple_index_S_fine= zeros(Int64,ns_fine,3);
        λ_2_S_fine= zeros(ns_fine,3);
        future_asset_S= zeros(ns_fine,3);
        future_asset_C= zeros(ns_fine,3);
        c_S_B_fine= zeros(ns_fine,3);
        c_B_B_fine= zeros(ns_fine,3);
        c_M_B_fine= zeros(ns_fine,3);
        x_SC_fine= zeros(ns_fine,3);
        x_BC_fine= zeros(ns_fine,3);
        land_C_fine= zeros(ns_fine,3);
        λ_2_fine= zeros(ns_fine,3);
        P_B_fine= zeros(ns_fine,3);
        Y_B_fine= zeros(ns_fine,3);
        q_S_B_fine= zeros(ns_fine,3);
        q_B_B_fine= zeros(ns_fine,3);
        solve_cash_crop_index_B_fine= zeros(Int64,ns_fine,3);
        solve_staple_index_B_fine= zeros(Int64,ns_fine,3);
        TC_fine= zeros(ns_fine,3);
        for j = 1:3
            for jj =1:3
            if jj == 1
                Y_W_fine_policy[:,j] =  P_W_fine *cons_fine_local[:,j]  + Y_W_fine  .+ w*FM_W .+ w*F_W * (j != 1);
            end
                if jj == 2
                    (future_asset_S[:,j],Y_S_fine[:,j],c_S_S_fine[:,j],c_B_S_fine[:,j],c_M_S_fine[:,j],q_S_S_fine[:,j],
                    P_S_fine[:,j],x_S_S_fine[:,j],solve_staple_index_S_fine[:,j],λ_2_S_fine[:,j]) = policy_function_creator(
                    cons_fine_local[:,j],jj,j,P_W_fine,Y_W_fine,s_fine,r,ρ,w,coeff_λ_2_cashcrop_residual_unconstrained_fine,
                    coeff_λ_2_cashcrop_residual_constrained_fine,θ_fine,fspace_C_fine,ϕ_S_loc,ζ,τ_S,p_x,p_B,
                    p_M,ϕ_B_loc,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,a_min,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
                    coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,
                    labor_allocated_interior_c3a_fine,λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,
                    q_S_c2_fine,q_B_c1_fine,q_S_c3a_fine,q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,
                    Y_B_c3b_fine,c_S_mat_fine,c_B_mat_fine,c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,q_B_mat_fine,land_B_mat_fine,λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,
                    feasibility_mat_fine,C_max_mat_fine,C_min_mat_fine,q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,x_S_staples_fine,
                    λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat_fine,C_max_staple_fine,C_min_staple_fine,C_max_staple_constrained_fine,
                    C_min_staple_constrained_fine,TC_S_c3_constrained_fine,x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine,
                    x_S_mat_3c_fine,x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine);
                elseif jj == 3
                    (future_asset_C[:,j],c_S_B_fine[:,j],c_B_B_fine[:,j],c_M_B_fine[:,j],x_SC_fine[:,j],x_BC_fine[:,j],
                    land_C_fine[:,j],λ_2_fine[:,j],P_B_fine[:,j],Y_B_fine[:,j],q_S_B_fine[:,j],q_B_B_fine[:,j],solve_cash_crop_index_B_fine[:,j]
                    ,solve_staple_index_B_fine[:,j],TC_fine[:,j]) = policy_function_creator(
                    cons_fine_local[:,j],jj,j,P_W_fine,Y_W_fine,s_fine,r,ρ,w,coeff_λ_2_cashcrop_residual_unconstrained_fine,
                    coeff_λ_2_cashcrop_residual_constrained_fine,θ_fine,fspace_C_fine,ϕ_S_loc,ζ,τ_S,p_x,p_B,
                    p_M,ϕ_B_loc,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,a_min,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
                    coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,
                    labor_allocated_interior_c3a_fine,λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,
                    q_S_c2_fine,q_B_c1_fine,q_S_c3a_fine,q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,
                    Y_B_c3b_fine,c_S_mat_fine,c_B_mat_fine,c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,q_B_mat_fine,land_B_mat_fine,λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,
                    feasibility_mat_fine,C_max_mat_fine,C_min_mat_fine,q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,x_S_staples_fine,
                    λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat_fine,C_max_staple_fine,C_min_staple_fine,C_max_staple_constrained_fine,
                    C_min_staple_constrained_fine,TC_S_c3_constrained_fine,x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine,
                    x_S_mat_3c_fine,x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine);
                end
            end
        end
        matcheck= zeros(ns_fine,12)
        matcheck[:,1:3]= future_occupation_fine_local;
        matcheck[:,4:6] = future_asset_S;#s_fine
        matcheck[:,7:9] = future_asset_C;#s_fine
        #matcheck[:,6:8] = cons_fine_local
        #matcheck[:,9:11] = solve_cash_crop_index_B_fine
        # Distribution of current occupations
        current_distr_store = zeros(3 * ns_fine);

        past_distr_store = zeros(ns_fine,3);

        stay_workers = worker_past_dist.*(future_occupation_fine_local[:,1].==1);
        exit_staple_to_work = staple_past_dist.*(future_occupation_fine_local[:,2].==1);
        exit_cashcrop_to_work = cash_crop_past_dist.*(future_occupation_fine_local[:,3].==1);
        current_workers = stay_workers + exit_staple_to_work + exit_cashcrop_to_work;
        # Worker policy functions are easy:
        c_S_W_fine[:,1] = (c̄_S .+ (1 + Q_S)^(-ϵ)*ψ_S.^ϵ.*cons_fine_local[:,1].*P_W_fine.^ϵ);
        c_S_W_fine[:,2] = (c̄_S .+ (1 + Q_S)^(-ϵ)*ψ_S.^ϵ.*cons_fine_local[:,2].*P_W_fine.^ϵ);
        c_S_W_fine[:,3] = (c̄_S .+ (1 + Q_S)^(-ϵ)*ψ_S.^ϵ.*cons_fine_local[:,3].*P_W_fine.^ϵ) ;
        c_B_W_fine[:,1] = (p_B^(-ϵ)*ψ_B.^ϵ.*cons_fine_local[:,1].*P_W_fine.^ϵ);
        c_B_W_fine[:,2] = (p_B^(-ϵ)*ψ_B.^ϵ.*cons_fine_local[:,2].*P_W_fine.^ϵ);
        c_B_W_fine[:,3] = (p_B^(-ϵ)*ψ_B.^ϵ.*cons_fine_local[:,3].*P_W_fine.^ϵ) ;
        c_M_W_fine[:,1] = (p_M^(-ϵ)*ψ_M.^ϵ.*cons_fine_local[:,1].*P_W_fine.^ϵ);
        c_M_W_fine[:,2] = (p_M^(-ϵ)*ψ_M.^ϵ.*cons_fine_local[:,2].*P_W_fine.^ϵ);
        c_M_W_fine[:,3] = (p_M^(-ϵ)*ψ_M.^ϵ.*cons_fine_local[:,3].*P_W_fine.^ϵ) ;

        c_S_worker_sum = sum( c_S_W_fine[:,1].*stay_workers + c_S_W_fine[:,2].*exit_staple_to_work +
        c_S_W_fine[:,3] .*exit_cashcrop_to_work );
        c_B_worker_sum = sum(c_B_W_fine[:,1] .*stay_workers +c_B_W_fine[:,2] .*exit_staple_to_work +
        c_B_W_fine[:,3] .*exit_cashcrop_to_work );
        c_M_worker_sum = sum(c_M_W_fine[:,1] .*stay_workers +c_M_W_fine[:,2] .*exit_staple_to_work +
        c_M_W_fine[:,3] .*exit_cashcrop_to_work );
        urban_labor_supply_sum = sum(current_workers.*labor_prod_fine);
        transaction_cost_worker_sum = c_S_worker_sum * Q_S;

        entrants_staple_from_workers = worker_past_dist.*(future_occupation_fine_local[:,1].==2);
        incumbents_staple = staple_past_dist.*(future_occupation_fine_local[:,2].==2);
        exit_cashcrop_to_staple = cash_crop_past_dist.*(future_occupation_fine_local[:,3].==2);
        current_staple = entrants_staple_from_workers + incumbents_staple + exit_cashcrop_to_staple;
        # Staple policy function now depends on consumption, hence the decision is always different!


        c_S_staple_sum = sum(c_S_S_fine[:,1] .*entrants_staple_from_workers + c_S_S_fine[:,2] .*incumbents_staple
        + c_S_S_fine[:,3] .*exit_cashcrop_to_staple);
        c_B_staple_sum = sum(c_B_S_fine[:,1] .*entrants_staple_from_workers + c_B_S_fine[:,2] .*incumbents_staple
        + c_B_S_fine[:,3] .*exit_cashcrop_to_staple);
        c_M_staple_sum = sum(c_M_S_fine[:,1] .*entrants_staple_from_workers + c_M_S_fine[:,2] .*incumbents_staple
        + c_M_S_fine[:,3] .*exit_cashcrop_to_staple);
        q_S_staple_sum = sum(q_S_S_fine[:,1] .*entrants_staple_from_workers + q_S_S_fine[:,2] .*incumbents_staple
        + q_S_S_fine[:,3] .*exit_cashcrop_to_staple);
        x_S_staple_sum = sum(x_S_S_fine[:,1] .*entrants_staple_from_workers + x_S_S_fine[:,2] .*incumbents_staple
        + x_S_S_fine[:,3] .*exit_cashcrop_to_staple);

        transaction_cost_staple_sum = Q_S * sum(max.(c_S_S_fine[:,1] - q_S_S_fine[:,1],0) .*entrants_staple_from_workers +
            max.(c_S_S_fine[:,2] - q_S_S_fine[:,2],0) .*incumbents_staple +
            max.(c_S_S_fine[:,3] - q_S_S_fine[:,3],0) .*exit_cashcrop_to_staple);

        entrants_cashcrop_from_workers = worker_past_dist.*(future_occupation_fine_local[:,1].==3);
        entrants_from_staple_to_cashcrop= staple_past_dist.*(future_occupation_fine_local[:,2].==3);
        incumbents_cashcrop = cash_crop_past_dist.*(future_occupation_fine_local[:,3].==3);
        current_cashcrop = entrants_cashcrop_from_workers + incumbents_cashcrop + entrants_from_staple_to_cashcrop;
        # Cash crop producer policy functions now depends on consumption, hence the decision is always different!

        c_S_cashcrop_sum = sum(c_S_B_fine[:,1] .*entrants_cashcrop_from_workers + c_S_B_fine[:,3] .*incumbents_cashcrop
        + c_S_B_fine[:,2] .*entrants_from_staple_to_cashcrop);
        c_B_Cashcrop_sum = sum(c_B_B_fine[:,1] .*entrants_cashcrop_from_workers + c_B_B_fine[:,3] .*incumbents_cashcrop
        + c_B_B_fine[:,2] .*entrants_from_staple_to_cashcrop);
        c_M_cashcrop_sum = sum(c_M_B_fine[:,1] .*entrants_cashcrop_from_workers + c_M_B_fine[:,3] .*incumbents_cashcrop
        + c_M_B_fine[:,2] .*entrants_from_staple_to_cashcrop);
        q_S_cashcrop_sum = sum(q_S_B_fine[:,1] .*entrants_cashcrop_from_workers + q_S_B_fine[:,3] .*incumbents_cashcrop
        + q_S_B_fine[:,2] .*entrants_from_staple_to_cashcrop);
        q_B_cashcrop_sum = sum(q_B_B_fine[:,1] .*entrants_cashcrop_from_workers + q_B_B_fine[:,3] .*incumbents_cashcrop
        + q_B_B_fine[:,2] .*entrants_from_staple_to_cashcrop);
        x_S_cashcrop_sum = sum(x_SC_fine[:,1] .*entrants_cashcrop_from_workers + x_SC_fine[:,3] .*incumbents_cashcrop
        + x_SC_fine[:,2] .*entrants_from_staple_to_cashcrop);
        x_B_cashcrop_sum = sum(x_BC_fine[:,1] .*entrants_cashcrop_from_workers + x_BC_fine[:,3] .*incumbents_cashcrop
        + x_BC_fine[:,2] .*entrants_from_staple_to_cashcrop);

        transaction_cost_cashcrop_sum = Q_S * sum(max.(c_S_B_fine[:,1] - q_S_B_fine[:,1],0) .*entrants_cashcrop_from_workers +
            max.(c_S_B_fine[:,2] - q_S_B_fine[:,2],0) .*entrants_from_staple_to_cashcrop +
            max.(c_S_B_fine[:,3] - q_S_B_fine[:,3],0) .*incumbents_cashcrop);

        current_distr_store[(ns_fine *0 + 1):(ns_fine *1)] = current_workers;
        current_distr_store[(ns_fine *1 + 1):(ns_fine *2)] = current_staple;
        current_distr_store[(ns_fine *2 + 1):(ns_fine *3)] = current_cashcrop;


        past_distr_store[:,1]=distr_previous[(ns_fine *0 + 1):(ns_fine *1)];
        past_distr_store[:,2]=distr_previous[(ns_fine *1 + 1):(ns_fine *2)];
        past_distr_store[:,3]=distr_previous[(ns_fine *2 + 1):(ns_fine *3)];

        current_worker_pop = sum(current_workers);
        current_staple_pop = sum(current_staple);
        current_cashcrop_pop = sum(current_cashcrop);
        # Entry cost accounting - measured in labor, convert with wages later
        entry_costs_to_workers =  F_W*(sum(exit_staple_to_work) + sum(exit_cashcrop_to_work));
        entry_costs_to_staples = F_S*(sum(entrants_staple_from_workers)); #Not needed: + sum(exit_cashcrop_to_work));
        entry_costs_to_cashcrops = F_B*(sum(entrants_cashcrop_from_workers) + sum(entrants_from_staple_to_cashcrop));
        total_entry_cost = entry_costs_to_workers + entry_costs_to_staples + entry_costs_to_cashcrops;
        # Maintenance cost accounting
        maintenance_costs_for_workers =  FM_W*(current_worker_pop);
        maintenance_costs_to_staples = FM_S*(current_staple_pop); #Not needed: + sum(exit_cashcrop_to_work));
        maintenance_costs_to_cashcrops = FM_B*(current_cashcrop_pop);
        total_maintenance_cost = maintenance_costs_for_workers + maintenance_costs_to_staples + maintenance_costs_to_cashcrops;

        # Savings and debt
        worker_bond_holding = s_fine[:,1] .*current_workers;
        worker_bond_holding_sum = sum(worker_bond_holding);
        staple_bond_holding = s_fine[:,1] .*current_staple;
        staple_bond_holding_sum = sum(staple_bond_holding[staple_bond_holding .> 0.0]);
        staple_debt_holding_sum = sum(staple_bond_holding[staple_bond_holding .< 0.0]);
        cashcrop_bond_holding = s_fine[:,1] .*current_cashcrop;
        cashcrop_bond_holding_sum = sum(cashcrop_bond_holding[cashcrop_bond_holding .> 0.0]);
        cashcrop_debt_holding_sum = sum(cashcrop_bond_holding[cashcrop_bond_holding .< 0.0]);
        asset_supply = worker_bond_holding_sum + (staple_bond_holding_sum + staple_debt_holding_sum ) + (cashcrop_bond_holding_sum + cashcrop_debt_holding_sum);


        # Manufacturing firm:
        labor_used = urban_labor_supply_sum - (total_maintenance_cost + total_entry_cost);

        # Govt spending side as only the manufacturing firm pays taxes --- # KAROL FIXED THIS
        input_staple= x_S_staple_sum + x_S_cashcrop_sum; # change current_staple to current_staple2 ? etc.
        input_cashcrop= copy(x_B_cashcrop_sum);
        foreign_demand_cash_crop = a_D*p_B^b_D#foreign_demand/p_B;
        #println(input_staple)
        #println(input_cashcrop)
        
        # Alternatively, govt borrowing would make a ton of sense here...
        Government_expenditure = balanced_share * p_x * (τ_S * input_staple + τ_B * input_cashcrop);
        Import_value = p_x * ( input_staple + input_cashcrop);
        Export_value = p_B * foreign_demand_cash_crop;
        
        #current_account = (Government_expenditure + Net_Foreign_Factor_Income - Import_value+Export_value)
        current_account_residual = Export_value - Import_value;
        #Net_Foreign_Factor_Income = -current_account_residual;
        #foreign_supply_capital = -Net_Foreign_Factor_Income/R;
        capital_used = copy(asset_supply) + foreign_supply_capital;
        Net_Foreign_Factor_Income = R*(copy(asset_supply) - capital_used)

        # if capital_used<0 || labor_used<0
        #     println("ERROR? Capital_used: ", capital_used, " Labor_used: ", labor_used)
        # end
        #
        # KLratio=α/(1-α)*((1+τ_W)*w)/r
        # r_new = -δ + α*p_M * KLratio^(α-1) #(max(0.01,capital_used)/max(0.01,labor_used))^(α-1)
        # w_new = (1 - α)/(1 + τ_W)*p_M* KLratio^α #(max(0.01,capital_used)/max(0.01,labor_used))^α # this is when we calculate wages paid by the manufacturing firm
        τ_W_new = - balanced_share * p_x * (τ_S * input_staple + τ_B * input_cashcrop) / (labor_used * w);
        K_L_ratio = maximum([capital_used./labor_used, tol]);
        r_new = -δ + α*p_M *K_L_ratio^(α-1);
        w_new = (1 - α)/(1 + τ_W)*p_M* K_L_ratio^α; # this is when we calculate wages paid by the manufacturing firm
        #labor_demand=((1 - α)/(1 + τ_W)/w_new*p_M*capital_used^α)^(1/α);

        prod_staple = q_S_staple_sum + q_S_cashcrop_sum;
        prod_cashcrop = q_B_cashcrop_sum;
        prod_manuf = sum(max(tol,capital_used).^α.*max(tol,labor_used).^(1-α));
        #foreign_demand = p_x * ( input_staple + input_cashcrop);
        
        #p_M_new =
        #τ_W*(labor_used * w)

                # Production -- KAROL ADDED THIS
        residual_goods = zeros(6);
        #residual[1] = (cash_crop_cons + p_x/p_B*input_staple + p_x/p_B*input_cashcrop - prod_cashcrop)/(cash_crop_cons + p_x/p_B*input_staple + p_x/p_B*input_cashcrop + prod_cashcrop); # -- KAROL ADDED THIS
        #residual[1] = 0#(cash_crop_cons + p_x*input_staple + p_x*input_cashcrop - prod_cashcrop)/(cash_crop_cons + p_x*input_staple + p_x*input_cashcrop + prod_cashcrop); # -- KAROL ADDED THIS
        #residual[1] = (cash_crop_cons  - prod_cashcrop)/(cash_crop_cons + prod_cashcrop); # -- KAROL ADDED THIS
        residual_goods[1] = (c_B_worker_sum + c_B_staple_sum + c_B_Cashcrop_sum+ foreign_demand_cash_crop  - prod_cashcrop)/(
                            c_B_worker_sum + c_B_staple_sum + c_B_Cashcrop_sum+ foreign_demand_cash_crop  + prod_cashcrop);

        #println(c_B_worker_sum + c_B_staple_sum + c_B_Cashcrop_sum)
        #println(foreign_demand_cash_crop)
        #println(prod_cashcrop)

        residual_goods[2] = (c_M_worker_sum + c_M_staple_sum + c_M_cashcrop_sum - prod_manuf)/(c_M_worker_sum + c_M_staple_sum + c_M_cashcrop_sum + prod_manuf);

        #residual[2] = (manuf_cons + input_staple + input_cashcrop - prod_manuf)/(manuf_cons + input_staple + input_cashcrop + prod_manuf); # -- KAROL ADDED THIS
        residual_goods[3] = 0 #No longer needed for the calibration as the foreign_supply_capital offsets capital markets(r - r_new)/(r+r_new);


        #residual[4] = (labor_demand - labor_used)/(labor_demand + labor_used);
        #staple_mkt_clr = (sum(c_current_S) + p_x*input_staple - prod_staple)/(sum(c_current_S) + prod_staple + p_x*input_staple) # -- KAROL ADDED THIS
        #staple_mkt_clr = (sum(c_current_S) - prod_staple)/(sum(c_current_S) + prod_staple) # -- KAROL ADDED THIS
        residual_goods[5] = (c_S_worker_sum + c_S_staple_sum + c_S_cashcrop_sum +
            transaction_cost_worker_sum + transaction_cost_staple_sum + transaction_cost_cashcrop_sum - prod_staple)/(
            c_S_worker_sum + c_S_staple_sum + c_S_cashcrop_sum +
                transaction_cost_worker_sum + transaction_cost_staple_sum + transaction_cost_cashcrop_sum + prod_staple);# This is the main residual, since other markets are tempered with
        if balanced_share >0.0
            residual_goods[6] = (τ_W - τ_W_new)/(τ_W + τ_W_new)
        else
            residual_goods[6] = 0.0
        end
        undernutritioned_workers = sum((c_S_W_fine[:, 1] .< cons_level_substinence) .* stay_workers
        + (c_S_W_fine[:, 2] .< cons_level_substinence) .* exit_staple_to_work +
        (c_S_W_fine[:, 3] .< cons_level_substinence) .* exit_cashcrop_to_work)
        undernutritioned_staple_farmer = sum((c_S_S_fine[:, 1] .< cons_level_substinence) .* entrants_staple_from_workers +
                    (c_S_S_fine[:, 2] .< cons_level_substinence) .* incumbents_staple +
                    (c_S_S_fine[:, 3] .< cons_level_substinence) .* exit_cashcrop_to_staple)
        undernutritioned_cashcrop_farmer = sum((c_S_B_fine[:, 1] .< cons_level_substinence) .* entrants_cashcrop_from_workers
                        + (c_S_B_fine[:, 2] .< cons_level_substinence) .* entrants_from_staple_to_cashcrop +
                        (c_S_B_fine[:, 3] .< cons_level_substinence) .* incumbents_cashcrop)

        cttilde_new = undernutritioned_workers + undernutritioned_staple_farmer + undernutritioned_cashcrop_farmer
        residual_goods[4] = (cttilde_new - cttilde)/(cttilde_new+cttilde);#(w - w_new)/(w+w_new);

        # if current_staple_pop ==0
        #     residual[5] = 100.0;
        # end

        #rural_pop_only_staples_printable=(sum(entrants_cashcrop_from_workers.* (q_B_B_fine[:,1].==0.0)) + sum(incumbents_cashcrop.*(q_B_B_fine[:,3].==0.0))
        #+ sum( entrants_from_staple_to_cashcrop.*(q_B_B_fine[:,2].==0.0)) + current_staple_pop);
        #rural_pop_only_Bashcrop_printable=(sum(entrants_cashcrop_from_workers.* (land_C_fine[:,1].==1.0)) + sum(incumbents_cashcrop.*(land_C_fine[:,3].==1.0))
        #+ sum( entrants_from_staple_to_cashcrop.*(land_C_fine[:,2].==1.0)));
        #println("Pops: Worker:",current_worker_pop," Staple:",current_staple_pop," Cash crop:",current_cashcrop_pop," Only producing staples", rural_pop_only_staples_printable
        #," Only producing cashcrop", rural_pop_only_Bashcrop_printable)
        mean_land_share_to_staples_among_cc_model =(sum(entrants_cashcrop_from_workers .* (land_C_fine[:,1].>0.0).* (1 .- land_C_fine[:,1])) +
        sum(incumbents_cashcrop.*(land_C_fine[:,3].>0.0).*(1 .- land_C_fine[:,3]))
            + sum( entrants_from_staple_to_cashcrop.*(land_C_fine[:,2].>0.0).*(1 .- land_C_fine[:,2])))/current_cashcrop_pop;
        avg_land_to_cashcrop = 1 - mean_land_share_to_staples_among_cc_model
        relative_land_to_cashcrop = avg_land_to_cashcrop *current_cashcrop_pop/(current_cashcrop_pop + current_staple_pop);
        relative_land_to_staples = 1 - relative_land_to_cashcrop;
        avg_land_to_staples = 1 - avg_land_to_cashcrop;
        staple_productivity = prod_staple/(avg_land_to_staples *current_cashcrop_pop  + current_staple_pop);
        cashcrop_productivity = prod_cashcrop/(avg_land_to_cashcrop *current_cashcrop_pop);
        manuf_productivity = prod_manuf/(current_worker_pop);
        aggregate_consumption = sum(cons_fine_local[:,1] .* distr_previous[1:ns_fine] +
        cons_fine_local[:,2] .* distr_previous[(ns_fine + 1):(2*ns_fine)] +
        cons_fine_local[:,3] .* distr_previous[(2*ns_fine + 1):(3*ns_fine)])

        mean_land_share_staples =   (mean_land_share_to_staples_among_cc_model*current_cashcrop_pop
        + current_staple_pop)/(1 - current_worker_pop ) 

        undernutritioned_workers = sum( (c_S_W_fine[:,1].<cons_level_substinence).*stay_workers
        + (c_S_W_fine[:,2].<cons_level_substinence).*exit_staple_to_work +
       (c_S_W_fine[:,3].<cons_level_substinence) .*exit_cashcrop_to_work );
       undernutritioned_staple_farmer = sum( (c_S_S_fine[:,1].<cons_level_substinence).*entrants_staple_from_workers +
       (c_S_S_fine[:,2].<cons_level_substinence).*incumbents_staple +
       (c_S_S_fine[:,3].<cons_level_substinence) .*exit_cashcrop_to_staple );
       undernutritioned_cashcrop_farmer = sum( (c_S_B_fine[:,1].<cons_level_substinence).*entrants_cashcrop_from_workers
        + (c_S_B_fine[:,2].<cons_level_substinence).*entrants_from_staple_to_cashcrop +
       (c_S_B_fine[:,3].<cons_level_substinence) .*incumbents_cashcrop );
       undernourished = undernutritioned_workers + undernutritioned_staple_farmer +undernutritioned_cashcrop_farmer
    fertilizer_use =    (input_staple + input_cashcrop);
    worker_pop_effective=(urban_labor_supply_sum - total_maintenance_cost )/urban_labor_supply_sum;
    #YL_manuf=(p_M*prod_manuf - w*labor_used - (r+ δ)*capital_used)/(current_worker_pop*worker_pop_effective);
    #YL_manuf=(p_M*prod_manuf - w*labor_used)/(current_worker_pop*worker_pop_effective);
    YL_manuf=(p_M*prod_manuf - w*total_entry_cost)/(current_worker_pop*worker_pop_effective);
    YL_agr=(prod_staple + p_B*prod_cashcrop - p_x*fertilizer_use - (total_maintenance_cost)*w)/(1-current_worker_pop*worker_pop_effective);
    APG=YL_manuf/YL_agr;

    # Return on "land"
    APland_cashcrop_B_val = zeros(ns_fine*3);
    APland_cashcrop_B_val[1:ns_fine] = q_B_B_fine[:,1].*p_B./land_C_fine[:,1];
    APland_cashcrop_B_val[(ns_fine+1):(2*ns_fine)] = q_B_B_fine[:,2].*p_B./land_C_fine[:,2];
    APland_cashcrop_B_val[(2*ns_fine+1):(3*ns_fine)] = q_B_B_fine[:,3].*p_B./land_C_fine[:,3]; 

    APland_cashcrop_S_val = zeros(ns_fine*3);
    APland_cashcrop_S_val[1:ns_fine] = q_S_B_fine[:,1]./(1 .-land_C_fine[:,1]);
    APland_cashcrop_S_val[(ns_fine+1):(2*ns_fine)] = q_S_B_fine[:,2]./(1 .-land_C_fine[:,2]);
    APland_cashcrop_S_val[(2*ns_fine+1):(3*ns_fine)] = q_S_B_fine[:,3]./(1 .-land_C_fine[:,3]); 

    APland_staples_S_val = zeros(ns_fine*3);
    APland_staples_S_val[1:ns_fine] = q_S_S_fine[:,1];
    APland_staples_S_val[(ns_fine+1):(2*ns_fine)] = q_S_S_fine[:,2];
    APland_staples_S_val[(2*ns_fine+1):(3*ns_fine)] = q_S_S_fine[:,3]; 

    APland_cashcrop_val = zeros(ns_fine*3);
    APland_cashcrop_val[1:ns_fine] = (q_S_B_fine[:,1] + p_B * q_B_B_fine[:,1]);
    APland_cashcrop_val[(ns_fine+1):(2*ns_fine)] = (q_S_B_fine[:,2] + p_B * q_B_B_fine[:,2]);
    APland_cashcrop_val[(2*ns_fine+1):(3*ns_fine)] = (q_S_B_fine[:,3] + p_B * q_B_B_fine[:,3]);

    APland_cashcrop_B_val[isnan.(APland_cashcrop_B_val)].=0.0001;
    APland_staples_S_val[isnan.(APland_staples_S_val)].=0.0001;
    APland_cashcrop_S_val[isnan.(APland_cashcrop_S_val)].=0.0001; 
    # sum([entrants_cashcrop_from_workers
    #  entrants_from_staple_to_cashcrop
    #   incumbents_cashcrop] .* isnan.(APland_cashcrop_B_val))
    #   sum(isnan.(APland_cashcrop_B_val))
        # Construct the means:
    APland_mean_cashcrop_B_log = sum(log.(APland_cashcrop_B_val[1:ns_fine]) .* entrants_cashcrop_from_workers/current_cashcrop_pop +
    log.(APland_cashcrop_B_val[(ns_fine+1):(2*ns_fine)]) .* entrants_from_staple_to_cashcrop/current_cashcrop_pop +
    log.(APland_cashcrop_B_val[(2*ns_fine+1):(3*ns_fine)]) .* incumbents_cashcrop/current_cashcrop_pop)

    APland_sq_cashcrop_B_log = sum(log.(APland_cashcrop_B_val[1:ns_fine]).^2 .* entrants_cashcrop_from_workers/current_cashcrop_pop +
    log.(APland_cashcrop_B_val[(ns_fine+1):(2*ns_fine)]).^2 .* entrants_from_staple_to_cashcrop/current_cashcrop_pop +
    log.(APland_cashcrop_B_val[(2*ns_fine+1):(3*ns_fine)]).^2 .* incumbents_cashcrop/current_cashcrop_pop)

    var_APland_cashcrop_B = APland_sq_cashcrop_B_log-APland_mean_cashcrop_B_log.^2;

    APland_mean_cashcrop_S_log = sum(log.(APland_cashcrop_S_val[1:ns_fine]) .* entrants_cashcrop_from_workers/current_cashcrop_pop +
    log.(APland_cashcrop_S_val[(ns_fine+1):(2*ns_fine)]) .* entrants_from_staple_to_cashcrop/current_cashcrop_pop +
    log.(APland_cashcrop_S_val[(2*ns_fine+1):(3*ns_fine)]) .* incumbents_cashcrop/current_cashcrop_pop)

    APland_sq_cashcrop_S_log = sum(log.(APland_cashcrop_S_val[1:ns_fine]).^2 .* entrants_cashcrop_from_workers/current_cashcrop_pop +
    log.(APland_cashcrop_S_val[(ns_fine+1):(2*ns_fine)]).^2 .* entrants_from_staple_to_cashcrop/current_cashcrop_pop +
    log.(APland_cashcrop_S_val[(2*ns_fine+1):(3*ns_fine)]).^2 .* incumbents_cashcrop/current_cashcrop_pop)

    var_APland_cashcrop_S = APland_sq_cashcrop_S_log-APland_mean_cashcrop_S_log.^2;

    APland_mean_cashcrop_log = sum(log.((APland_cashcrop_S_val[1:ns_fine].*x_SC_fine[:,1] + p_B *APland_cashcrop_B_val[1:ns_fine].*x_BC_fine[:,1])./(x_SC_fine[:,1] + x_BC_fine[:,1])) .* entrants_cashcrop_from_workers/current_cashcrop_pop +
    log.((APland_cashcrop_S_val[(ns_fine+1):(2*ns_fine)].*x_SC_fine[:,2] + p_B *APland_cashcrop_B_val[(ns_fine+1):(2*ns_fine)].*x_BC_fine[:,2])./(x_SC_fine[:,2] + x_BC_fine[:,2])) .* entrants_from_staple_to_cashcrop/current_cashcrop_pop +
    log.((APland_cashcrop_S_val[(2*ns_fine+1):(3*ns_fine)].*x_SC_fine[:,3] + p_B *APland_cashcrop_B_val[(2*ns_fine+1):(3*ns_fine)].*x_BC_fine[:,3])./(x_SC_fine[:,3] + x_BC_fine[:,3])).* incumbents_cashcrop/current_cashcrop_pop)
    
    APland_sq_cashcrop_log = sum(log.((APland_cashcrop_S_val[1:ns_fine].*x_SC_fine[:,1] + p_B *APland_cashcrop_B_val[1:ns_fine].*x_BC_fine[:,1])./(x_SC_fine[:,1] + x_BC_fine[:,1])).^2 .* entrants_cashcrop_from_workers/current_cashcrop_pop +
    log.((APland_cashcrop_S_val[(ns_fine+1):(2*ns_fine)].*x_SC_fine[:,2] + p_B *APland_cashcrop_B_val[(ns_fine+1):(2*ns_fine)].*x_BC_fine[:,2])./(x_SC_fine[:,2] + x_BC_fine[:,2])).^2 .* entrants_from_staple_to_cashcrop/current_cashcrop_pop +
    log.((APland_cashcrop_S_val[(2*ns_fine+1):(3*ns_fine)].*x_SC_fine[:,3] + p_B *APland_cashcrop_B_val[(2*ns_fine+1):(3*ns_fine)].*x_BC_fine[:,3])./(x_SC_fine[:,3] + x_BC_fine[:,3])).^2 .* incumbents_cashcrop/current_cashcrop_pop)
    
    var_APland_cashcrop = APland_sq_cashcrop_log-APland_mean_cashcrop_log.^2;

    APland_mean_staples_S_log = sum(log.(APland_staples_S_val[1:ns_fine]) .* entrants_staple_from_workers/current_staple_pop +
    log.(APland_staples_S_val[(ns_fine+1):(2*ns_fine)]) .* incumbents_staple/current_staple_pop +
    log.(APland_staples_S_val[(2*ns_fine+1):(3*ns_fine)]) .* exit_cashcrop_to_staple/current_staple_pop)

    APland_sq_staples_S_log = sum(log.(APland_staples_S_val[1:ns_fine]).^2 .* entrants_staple_from_workers/current_staple_pop +
    log.(APland_staples_S_val[(ns_fine+1):(2*ns_fine)]).^2 .* incumbents_staple/current_staple_pop +
    log.(APland_staples_S_val[(2*ns_fine+1):(3*ns_fine)]).^2 .* exit_cashcrop_to_staple/current_staple_pop)

    var_APland_staples_S = APland_sq_staples_S_log-APland_mean_staples_S_log.^2;

    # For the entire agricultural sector:

    APland_mean_log = sum(log.(APland_staples_S_val[1:ns_fine]) .* entrants_staple_from_workers/(current_staple_pop + current_cashcrop_pop) +
    log.(APland_staples_S_val[(ns_fine+1):(2*ns_fine)]) .* incumbents_staple/(current_staple_pop + current_cashcrop_pop) +
    log.(APland_staples_S_val[(2*ns_fine+1):(3*ns_fine)]) .* exit_cashcrop_to_staple/(current_staple_pop + current_cashcrop_pop) +
    log.((APland_cashcrop_S_val[1:ns_fine].*x_SC_fine[:,1] + p_B *APland_cashcrop_B_val[1:ns_fine].*x_BC_fine[:,1])./(x_SC_fine[:,1] + x_BC_fine[:,1])) .* entrants_cashcrop_from_workers/(current_staple_pop + current_cashcrop_pop) +
    log.((APland_cashcrop_S_val[(ns_fine+1):(2*ns_fine)].*x_SC_fine[:,2] + p_B *APland_cashcrop_B_val[(ns_fine+1):(2*ns_fine)].*x_BC_fine[:,2])./(x_SC_fine[:,2] + x_BC_fine[:,2])) .* entrants_from_staple_to_cashcrop/(current_staple_pop + current_cashcrop_pop) +
    log.((APland_cashcrop_S_val[(2*ns_fine+1):(3*ns_fine)].*x_SC_fine[:,3] + p_B *APland_cashcrop_B_val[(2*ns_fine+1):(3*ns_fine)].*x_BC_fine[:,3])./(x_SC_fine[:,3] + x_BC_fine[:,3])) .* incumbents_cashcrop/(current_staple_pop + current_cashcrop_pop))

    APland_sq_log = sum(log.(APland_staples_S_val[1:ns_fine]).^2 .* entrants_staple_from_workers/(current_staple_pop + current_cashcrop_pop) +
    log.(APland_staples_S_val[(ns_fine+1):(2*ns_fine)]).^2 .* incumbents_staple/(current_staple_pop + current_cashcrop_pop) +
    log.(APland_staples_S_val[(2*ns_fine+1):(3*ns_fine)]).^2 .* exit_cashcrop_to_staple/(current_staple_pop + current_cashcrop_pop) +
    log.((APland_cashcrop_S_val[1:ns_fine].*x_SC_fine[:,1] + p_B *APland_cashcrop_B_val[1:ns_fine].*x_BC_fine[:,1])./(x_SC_fine[:,1] + x_BC_fine[:,1])).^2 .* entrants_cashcrop_from_workers/(current_staple_pop + current_cashcrop_pop) +
    log.((APland_cashcrop_S_val[(ns_fine+1):(2*ns_fine)].*x_SC_fine[:,2] + p_B *APland_cashcrop_B_val[(ns_fine+1):(2*ns_fine)].*x_BC_fine[:,2])./(x_SC_fine[:,2] + x_BC_fine[:,2])).^2 .* entrants_from_staple_to_cashcrop/(current_staple_pop + current_cashcrop_pop) +
    log.((APland_cashcrop_S_val[(2*ns_fine+1):(3*ns_fine)].*x_SC_fine[:,3] + p_B *APland_cashcrop_B_val[(2*ns_fine+1):(3*ns_fine)].*x_BC_fine[:,3])./(x_SC_fine[:,3] + x_BC_fine[:,3])).^2 .* incumbents_cashcrop/(current_staple_pop + current_cashcrop_pop))

    var_APland = APland_sq_log-APland_mean_log.^2;
    ### Return on fertilizers
    # Note that MPX = ζ * APX
    # MPX
    MPX_cashcrop_B_val = zeros(ns_fine*3);
    MPX_cashcrop_B_val[1:ns_fine] = ϕ_B_loc * θ_fine .* land_C_fine[:,1].^ρ .* ζ.*x_BC_fine[:,1].^(ζ-1);
    MPX_cashcrop_B_val[(ns_fine+1):(2*ns_fine)] = ϕ_B_loc * θ_fine .* land_C_fine[:,2].^ρ .* ζ.*x_BC_fine[:,2].^(ζ-1);
    MPX_cashcrop_B_val[(2*ns_fine+1):(3*ns_fine)] = ϕ_B_loc * θ_fine .* land_C_fine[:,3].^ρ .* ζ.*x_BC_fine[:,3].^(ζ-1); 

    MPX_cashcrop_S_val = zeros(ns_fine*3);
    MPX_cashcrop_S_val[1:ns_fine] = ϕ_S_loc * θ_fine .* (1 .- land_C_fine[:,1]).^ρ .* ζ.*x_SC_fine[:,1].^(ζ-1);
    MPX_cashcrop_S_val[(ns_fine+1):(2*ns_fine)] = ϕ_S_loc * θ_fine .* (1 .- land_C_fine[:,2]).^ρ .* ζ.*x_SC_fine[:,2].^(ζ-1);
    MPX_cashcrop_S_val[(2*ns_fine+1):(3*ns_fine)] = ϕ_S_loc * θ_fine .* (1 .- land_C_fine[:,3]).^ρ .* ζ.*x_SC_fine[:,3].^(ζ-1); 

    MPX_staples_S_val = zeros(ns_fine*3);
    MPX_staples_S_val[1:ns_fine] = ϕ_S_loc * θ_fine .* ζ.*x_S_S_fine[:,1].^(ζ-1);
    MPX_staples_S_val[(ns_fine+1):(2*ns_fine)] = ϕ_S_loc * θ_fine .* ζ.*x_S_S_fine[:,2].^(ζ-1);
    MPX_staples_S_val[(2*ns_fine+1):(3*ns_fine)] = ϕ_S_loc * θ_fine .* ζ.*x_S_S_fine[:,3].^(ζ-1);
    MPX_cashcrop_B_val[isnan.(MPX_cashcrop_B_val)].=0.0001;
    MPX_staples_S_val[isnan.(MPX_staples_S_val)].=0.0001;
    MPX_cashcrop_S_val[isnan.(MPX_cashcrop_S_val)].=0.0001; 
    # APX
    APX_cashcrop_B_val = zeros(ns_fine*3);
    APX_cashcrop_B_val[1:ns_fine] = q_B_B_fine[:,1]./x_BC_fine[:,1];
    APX_cashcrop_B_val[(ns_fine+1):(2*ns_fine)] = q_B_B_fine[:,2]./x_BC_fine[:,2];
    APX_cashcrop_B_val[(2*ns_fine+1):(3*ns_fine)] = q_B_B_fine[:,3]./x_BC_fine[:,3]; 

    APX_cashcrop_S_val = zeros(ns_fine*3);
    APX_cashcrop_S_val[1:ns_fine] = q_S_B_fine[:,1]./x_SC_fine[:,1];
    APX_cashcrop_S_val[(ns_fine+1):(2*ns_fine)] = q_S_B_fine[:,2]./x_SC_fine[:,2];
    APX_cashcrop_S_val[(2*ns_fine+1):(3*ns_fine)] = q_S_B_fine[:,3]./x_SC_fine[:,3]; 

    APX_staples_S_val = zeros(ns_fine*3);
    APX_staples_S_val[1:ns_fine] = q_S_S_fine[:,1] ./x_S_S_fine[:,1];
    APX_staples_S_val[(ns_fine+1):(2*ns_fine)] = q_S_S_fine[:,2] ./x_S_S_fine[:,2];
    APX_staples_S_val[(2*ns_fine+1):(3*ns_fine)] = q_S_S_fine[:,3] ./x_S_S_fine[:,3]; 

    APX_cashcrop_val = zeros(ns_fine*3);
    APX_cashcrop_val[1:ns_fine] = (q_S_B_fine[:,1] + p_B * q_B_B_fine[:,1])./(x_SC_fine[:,1] + x_BC_fine[:,1]);
    APX_cashcrop_val[(ns_fine+1):(2*ns_fine)] = (q_S_B_fine[:,2] + p_B * q_B_B_fine[:,2])./(x_SC_fine[:,2] + x_BC_fine[:,2]);
    APX_cashcrop_val[(2*ns_fine+1):(3*ns_fine)] = (q_S_B_fine[:,3] + p_B * q_B_B_fine[:,3])./(x_SC_fine[:,3] + x_BC_fine[:,3]);
    
    # Construct the means:
    MPX_mean_cashcrop_B_log = sum(log.(MPX_cashcrop_B_val[1:ns_fine]) .* entrants_cashcrop_from_workers/current_cashcrop_pop +
    log.(MPX_cashcrop_B_val[(ns_fine+1):(2*ns_fine)]) .* entrants_from_staple_to_cashcrop/current_cashcrop_pop +
    log.(MPX_cashcrop_B_val[(2*ns_fine+1):(3*ns_fine)]) .* incumbents_cashcrop/current_cashcrop_pop)

    MPX_sq_cashcrop_B_log = sum(log.(MPX_cashcrop_B_val[1:ns_fine]).^2 .* entrants_cashcrop_from_workers/current_cashcrop_pop +
    log.(MPX_cashcrop_B_val[(ns_fine+1):(2*ns_fine)]).^2 .* entrants_from_staple_to_cashcrop/current_cashcrop_pop +
    log.(MPX_cashcrop_B_val[(2*ns_fine+1):(3*ns_fine)]).^2 .* incumbents_cashcrop/current_cashcrop_pop)

    var_MPX_cashcrop_B = MPX_sq_cashcrop_B_log-MPX_mean_cashcrop_B_log.^2;

    MPX_mean_cashcrop_S_log = sum(log.(MPX_cashcrop_S_val[1:ns_fine]) .* entrants_cashcrop_from_workers/current_cashcrop_pop +
    log.(MPX_cashcrop_S_val[(ns_fine+1):(2*ns_fine)]) .* entrants_from_staple_to_cashcrop/current_cashcrop_pop +
    log.(MPX_cashcrop_S_val[(2*ns_fine+1):(3*ns_fine)]) .* incumbents_cashcrop/current_cashcrop_pop)

    MPX_sq_cashcrop_S_log = sum(log.(MPX_cashcrop_S_val[1:ns_fine]).^2 .* entrants_cashcrop_from_workers/current_cashcrop_pop +
    log.(MPX_cashcrop_S_val[(ns_fine+1):(2*ns_fine)]).^2 .* entrants_from_staple_to_cashcrop/current_cashcrop_pop +
    log.(MPX_cashcrop_S_val[(2*ns_fine+1):(3*ns_fine)]).^2 .* incumbents_cashcrop/current_cashcrop_pop)

    var_MPX_cashcrop_S = MPX_sq_cashcrop_S_log-MPX_mean_cashcrop_S_log.^2;

    MPX_mean_cashcrop_log = sum(log.((MPX_cashcrop_S_val[1:ns_fine].*x_SC_fine[:,1] + p_B *MPX_cashcrop_B_val[1:ns_fine].*x_BC_fine[:,1])./(x_SC_fine[:,1] + x_BC_fine[:,1])) .* entrants_cashcrop_from_workers/current_cashcrop_pop +
    log.((MPX_cashcrop_S_val[(ns_fine+1):(2*ns_fine)].*x_SC_fine[:,2] + p_B *MPX_cashcrop_B_val[(ns_fine+1):(2*ns_fine)].*x_BC_fine[:,2])./(x_SC_fine[:,2] + x_BC_fine[:,2])) .* entrants_from_staple_to_cashcrop/current_cashcrop_pop +
    log.((MPX_cashcrop_S_val[(2*ns_fine+1):(3*ns_fine)].*x_SC_fine[:,3] + p_B *MPX_cashcrop_B_val[(2*ns_fine+1):(3*ns_fine)].*x_BC_fine[:,3])./(x_SC_fine[:,3] + x_BC_fine[:,3])).* incumbents_cashcrop/current_cashcrop_pop)
    
    MPX_sq_cashcrop_log = sum(log.((MPX_cashcrop_S_val[1:ns_fine].*x_SC_fine[:,1] + p_B *MPX_cashcrop_B_val[1:ns_fine].*x_BC_fine[:,1])./(x_SC_fine[:,1] + x_BC_fine[:,1])).^2 .* entrants_cashcrop_from_workers/current_cashcrop_pop +
    log.((MPX_cashcrop_S_val[(ns_fine+1):(2*ns_fine)].*x_SC_fine[:,2] + p_B *MPX_cashcrop_B_val[(ns_fine+1):(2*ns_fine)].*x_BC_fine[:,2])./(x_SC_fine[:,2] + x_BC_fine[:,2])).^2 .* entrants_from_staple_to_cashcrop/current_cashcrop_pop +
    log.((MPX_cashcrop_S_val[(2*ns_fine+1):(3*ns_fine)].*x_SC_fine[:,3] + p_B *MPX_cashcrop_B_val[(2*ns_fine+1):(3*ns_fine)].*x_BC_fine[:,3])./(x_SC_fine[:,3] + x_BC_fine[:,3])).^2 .* incumbents_cashcrop/current_cashcrop_pop)
    
    var_MPX_cashcrop = MPX_sq_cashcrop_log-MPX_mean_cashcrop_log.^2;

    MPX_mean_staples_S_log = sum(log.(MPX_staples_S_val[1:ns_fine]) .* entrants_staple_from_workers/current_staple_pop +
    log.(MPX_staples_S_val[(ns_fine+1):(2*ns_fine)]) .* incumbents_staple/current_staple_pop +
    log.(MPX_staples_S_val[(2*ns_fine+1):(3*ns_fine)]) .* exit_cashcrop_to_staple/current_staple_pop)

    MPX_sq_staples_S_log = sum(log.(MPX_staples_S_val[1:ns_fine]).^2 .* entrants_staple_from_workers/current_staple_pop +
    log.(MPX_staples_S_val[(ns_fine+1):(2*ns_fine)]).^2 .* incumbents_staple/current_staple_pop +
    log.(MPX_staples_S_val[(2*ns_fine+1):(3*ns_fine)]).^2 .* exit_cashcrop_to_staple/current_staple_pop)

    var_MPX_staples_S = MPX_sq_staples_S_log-MPX_mean_staples_S_log.^2;

    # For the entire agricultural sector:

    MPX_mean_log = sum(log.(MPX_staples_S_val[1:ns_fine]) .* entrants_staple_from_workers/(current_staple_pop + current_cashcrop_pop) +
    log.(MPX_staples_S_val[(ns_fine+1):(2*ns_fine)]) .* incumbents_staple/(current_staple_pop + current_cashcrop_pop) +
    log.(MPX_staples_S_val[(2*ns_fine+1):(3*ns_fine)]) .* exit_cashcrop_to_staple/(current_staple_pop + current_cashcrop_pop) +
    log.((MPX_cashcrop_S_val[1:ns_fine].*x_SC_fine[:,1] + p_B *MPX_cashcrop_B_val[1:ns_fine].*x_BC_fine[:,1])./(x_SC_fine[:,1] + x_BC_fine[:,1])) .* entrants_cashcrop_from_workers/(current_staple_pop + current_cashcrop_pop) +
    log.((MPX_cashcrop_S_val[(ns_fine+1):(2*ns_fine)].*x_SC_fine[:,2] + p_B *MPX_cashcrop_B_val[(ns_fine+1):(2*ns_fine)].*x_BC_fine[:,2])./(x_SC_fine[:,2] + x_BC_fine[:,2])) .* entrants_from_staple_to_cashcrop/(current_staple_pop + current_cashcrop_pop) +
    log.((MPX_cashcrop_S_val[(2*ns_fine+1):(3*ns_fine)].*x_SC_fine[:,3] + p_B *MPX_cashcrop_B_val[(2*ns_fine+1):(3*ns_fine)].*x_BC_fine[:,3])./(x_SC_fine[:,3] + x_BC_fine[:,3])) .* incumbents_cashcrop/(current_staple_pop + current_cashcrop_pop))

    MPX_sq_log = sum(log.(MPX_staples_S_val[1:ns_fine]).^2 .* entrants_staple_from_workers/(current_staple_pop + current_cashcrop_pop) +
    log.(MPX_staples_S_val[(ns_fine+1):(2*ns_fine)]).^2 .* incumbents_staple/(current_staple_pop + current_cashcrop_pop) +
    log.(MPX_staples_S_val[(2*ns_fine+1):(3*ns_fine)]).^2 .* exit_cashcrop_to_staple/(current_staple_pop + current_cashcrop_pop) +
    log.((MPX_cashcrop_S_val[1:ns_fine].*x_SC_fine[:,1] + p_B *MPX_cashcrop_B_val[1:ns_fine].*x_BC_fine[:,1])./(x_SC_fine[:,1] + x_BC_fine[:,1])).^2 .* entrants_cashcrop_from_workers/(current_staple_pop + current_cashcrop_pop) +
    log.((MPX_cashcrop_S_val[(ns_fine+1):(2*ns_fine)].*x_SC_fine[:,2] + p_B *MPX_cashcrop_B_val[(ns_fine+1):(2*ns_fine)].*x_BC_fine[:,2])./(x_SC_fine[:,2] + x_BC_fine[:,2])).^2 .* entrants_from_staple_to_cashcrop/(current_staple_pop + current_cashcrop_pop) +
    log.((MPX_cashcrop_S_val[(2*ns_fine+1):(3*ns_fine)].*x_SC_fine[:,3] + p_B *MPX_cashcrop_B_val[(2*ns_fine+1):(3*ns_fine)].*x_BC_fine[:,3])./(x_SC_fine[:,3] + x_BC_fine[:,3])).^2 .* incumbents_cashcrop/(current_staple_pop + current_cashcrop_pop))

    var_MPX = MPX_sq_log-MPX_mean_log.^2;
    avg_labor_prod_rural = sum(labor_prod_fine .*(current_staple + current_cashcrop )/(current_cashcrop_pop + current_staple_pop))
    avg_labor_prod_urban = sum(labor_prod_fine .*(current_workers )/current_worker_pop)
    avg_agri_prod_rural = sum(θ_fine .*(current_staple + current_cashcrop )/(current_cashcrop_pop + current_staple_pop))
    avg_agri_prod_urban = sum(θ_fine .*(current_workers )/current_worker_pop)

    marketable_surplus_staple_sum =  sum(max.(q_S_S_fine[:,1] - c_S_S_fine[:,1],0) .*entrants_staple_from_workers +
    max.(q_S_S_fine[:,2] - c_S_S_fine[:,2],0) .*incumbents_staple +
    max.(q_S_S_fine[:,3] - c_S_S_fine[:,3],0) .*exit_cashcrop_to_staple);
    marketable_surplus_cashcrop_sum =  sum(max.(q_S_B_fine[:,1] - c_S_B_fine[:,1],0) .*entrants_cashcrop_from_workers +
        max.(q_S_B_fine[:,2] - c_S_B_fine[:,2],0) .*entrants_from_staple_to_cashcrop +
        max.(q_S_B_fine[:,3] - c_S_B_fine[:,3],0) .*incumbents_cashcrop);
    # marketable_surplus_cashcrop_sum =  sum(max.(q_S_B_fine[:,1],0) .*entrants_cashcrop_from_workers +
    #     max.(q_S_B_fine[:,2],0) .*entrants_from_staple_to_cashcrop +
    #     max.(q_S_B_fine[:,3],0) .*incumbents_cashcrop);
    marketable_agr_surplus=marketable_surplus_staple_sum+marketable_surplus_cashcrop_sum + p_B * prod_cashcrop
    nominal_GDP =  (marketable_agr_surplus + p_B * prod_cashcrop + p_M * prod_manuf);
    current_account_residual = current_account_residual /nominal_GDP 
    V_saved_reshaped = reshape(V_saved_local,ns_fine*3)
    return (residual_goods,distr_current,prod_staple,prod_cashcrop,prod_manuf,asset_supply,current_worker_pop,
    current_staple_pop,current_cashcrop_pop,current_account_residual,staple_productivity,cashcrop_productivity,
    manuf_productivity,aggregate_consumption,relative_land_to_cashcrop,mean_land_share_staples,undernourished,
    fertilizer_use,APG,var_APland,var_MPX,avg_labor_prod_rural,avg_labor_prod_urban,avg_agri_prod_rural,avg_agri_prod_urban,
    V_saved_reshaped,a_prime_fine_local)
end


function Residual_transition_sequential_epsilon_detailed(price_trans_actual::Array{Float64,2},capital_trans::Array{Float64,1},
    distr_store::Array{Float64,2},T::Int64,parameter_end::Parameter_type,coeff_store::Array{Float64,3},τ_trans::Array{Float64,1},r::Float64,V_saved_store::Array{Float64,2},cons_level_substinence::Float64,epsilon_u::Float64,epsilon_r::Float64)
    # Extract the necessary local parameters:
    (δ,ζ,ρ,α,σ,β,ϵ,ψ_S,ψ_B,ψ_M,ϕ_S,ϕ_B,c̄_S,F_W,F_S,F_B,FM_W,FM_S,FM_B,Q_S,p_x,τ_S,τ_B,a_D,b_D,K_a,K_b,γ,A_W,
    ρ_S,ρ_SW,σ_S,ρ_W,σ_W,n,n_fine,agrid,agrid_fine,a_min,a_max,spliorder,fspace_a,fspace_a_fine,fspace_C_fine,C_grid_fine_no,C_grid_fine,s,ns,
    s_fine,ns_fine,z,z_W,Phi_z,Phi_z_fine,Phi,Phi_aug,P_kron,P_kron1,P_kron_fine,κ) = local_parameters(parameter_end);
    # Policy function
    (   a_prime_fine_store,future_occupation_fine_local_store,
    cons_fine_local_store,c_S_W_fine_store,c_B_W_fine_store,c_M_W_fine_store,Y_W_fine_policy_store,Y_S_fine_store,
    c_S_S_fine_store,c_B_S_fine_store,c_M_S_fine_store,q_S_S_fine_store,P_S_fine_store,
    x_S_S_fine_store,solve_staple_index_S_fine_store,λ_2_S_fine_store,future_asset_S_store,
    future_asset_C_store,c_S_B_fine_store,c_B_B_fine_store,c_M_B_fine_store,x_SC_fine_store,x_BC_fine_store,
    land_C_fine_store,λ_2_fine_store,P_B_fine_store,Y_B_fine_store,q_S_B_fine_store,
    q_B_B_fine_store,solve_cash_crop_index_B_fine_store,solve_staple_index_B_fine_store,TC_fine_store,residual_store,
    coeff_λ_2_cashcrop_residual_unconstrained_store,coeff_λ_2_cashcrop_residual_constrained_store,C_max_unconstrained_store,
    C_max_constrained_store,C_min_unconstrained_store,C_min_constrained_store,coeff_λ_2_s_store,C_max_staple_store,C_min_staple_store,
    C_max_staple_constrained_store,C_min_staple_constrained_store,TC_S_c3_constrained_store,x_S_c3_constrained_store,q_S_c3_constrained_store,
    c_S_c3_constrained_store) =local_var_creator_policy(ns,ns_fine,T,C_grid_fine_no);
    (prod_staple_store,prod_cashcrop_store,prod_manuf_store,asset_supply_store,current_worker_pop_store,
    current_staple_pop_store,current_cashcrop_pop_store,current_account_residual_store,staple_productivity_store,
    cashcrop_productivity_store,manuf_productivity_store,aggregate_consumption_store,relative_land_to_cashcrop_store,
    mean_land_share_staples_store,undernourished_store,fertilizer_use_store,APG_store,var_APland_store,var_MPX_store,
    avg_labor_prod_rural_store,avg_labor_prod_urban_store,avg_agri_prod_rural_store,avg_agri_prod_urban_store) =local_var_creator_policy_details(ns,ns_fine,T,C_grid_fine_no);
    #Precheck the prices
    price_check_sum = 0;
    println("Backward loop")
    for t in (T+1):-1:2
        # Prices and other time dependent quantities:
        τ_S_loc = τ_trans[t];
        prices_loc = price_trans_actual[:,t-1];
        # Epsilon distortions
        cttilde=prices_loc[4]

        if (τ_S_loc==0 || prices_loc[3]==0.0)
            balanced_share = 0.0;
        else
            balanced_share = 1.0;
        end
        p_B,p_M,R,r,w,τ_W= price_reshaper_fixed_r(prices_loc[1:3],δ,ϵ,ψ_S,ψ_B,ψ_M,p_x,τ_S_loc,τ_B,α,balanced_share,r);
        price_check_tmp = p_B<p_B_min + p_B>p_B_max + p_M<p_M_min + p_M>p_M_max + τ_W<τ_W_min + τ_W>τ_W_max;
        # Price checks might be necessary - leave it for now.
        #0.01<p_B
        price_check_sum = price_check_tmp + price_check_sum;
        #price_check_sum = 0
        # Get the coefficients
        if price_check_sum>0.0
            print(t-1,",")
        else
            coefficients_next_tmp = coeff_store[:,:,t+1];
            coeff_next = coeff_store[:,:,t];
            (coeff_store[:,:,t],coeff_λ_2_cashcrop_residual_unconstrained_store[:,:,t-1],coeff_λ_2_cashcrop_residual_constrained_store[:,:,t-1],
            C_max_unconstrained_store[:,t-1],C_max_constrained_store[:,t-1],C_min_unconstrained_store[:,t-1],C_min_constrained_store[:,t-1],
            coeff_λ_2_s_store[:,:,t-1],C_max_staple_store[:,t-1],C_min_staple_store[:,t-1],C_max_staple_constrained_store[:,t-1],
            C_min_staple_constrained_store[:,t-1],TC_S_c3_constrained_store[:,t-1],x_S_c3_constrained_store[:,t-1],q_S_c3_constrained_store[:,t-1],
            c_S_c3_constrained_store[:,t-1] )= Residual_transition_backward_epsilon(s,ns,
            z,z_W,ϕ_S,ζ,τ_S_loc,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n,κ,Q_S,ϵ,ψ_S,ψ_B,ψ_M,agrid_fine,
            fspace_C_fine,agrid,coefficients_next_tmp,coeff_next,C_grid_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,P_kron,Phi,Phi_z,β,fspace_a,σ,cons_level_substinence, epsilon_u, epsilon_r,cttilde)
        end
        print(t-1,",")
    end

    if price_check_sum>0
        println("Guess out of bounds")
        residual_store = 1000*ones(6,T);
    else
        #Forward loop
        println("Forward loop")
        for t in 2:(T+1)
            # Prices and other time dependent quantities:
            τ_S_loc = τ_trans[t];
            prices_loc = price_trans_actual[:,t-1];
            # Epsilon distortions
            cttilde=prices_loc[4]

            if (τ_S_loc==0 || prices_loc[3]==0.0)
                balanced_share = 0.0;
            else
                balanced_share = 1.0;
            end
            p_B,p_M,R,r,w,τ_W= price_reshaper_fixed_r(prices_loc[1:3],δ,ϵ,ψ_S,ψ_B,ψ_M,p_x,τ_S_loc,τ_B,α,balanced_share,r);
            (residual_store[:,t-1],distr_store[:,t],prod_staple_store[t-1],prod_cashcrop_store[t-1],prod_manuf_store[t-1],asset_supply_store[t-1],
            current_worker_pop_store[t-1],current_staple_pop_store[t-1],current_cashcrop_pop_store[t-1],
            current_account_residual_store[t-1],staple_productivity_store[t-1],cashcrop_productivity_store[t-1],manuf_productivity_store[t-1],aggregate_consumption_store[t-1],relative_land_to_cashcrop_store[t-1],
            mean_land_share_staples_store[t-1],undernourished_store[t-1],
            fertilizer_use_store[t-1],APG_store[t-1],var_APland_store[t-1],var_MPX_store[t-1],avg_labor_prod_rural_store[t-1],avg_labor_prod_urban_store[t-1],avg_agri_prod_rural_store[t-1],avg_agri_prod_urban_store[t-1],
            V_saved_store[:,t],a_prime_fine_store[:,:,t-1]        ) = Residual_transition_forward_epsilon_detailed(coeff_store[:,:,t],distr_store[:,t-1],s_fine,ns_fine,
            z,z_W,ϕ_S,ζ,τ_S_loc,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n_fine,κ,Q_S,
            ϵ,ψ_S,ψ_B,ψ_M,coeff_λ_2_cashcrop_residual_unconstrained_store[:,:,t-1],
            coeff_λ_2_cashcrop_residual_constrained_store[:,:,t-1],
            C_max_unconstrained_store[:,t-1] ,C_max_constrained_store[:,t-1],C_min_unconstrained_store[:,t-1],
            C_min_constrained_store[:,t-1],coeff_λ_2_s_store[:,:,t-1],agrid_fine,fspace_C_fine,C_max_staple_store[:,t-1],
            C_min_staple_store[:,t-1],C_max_staple_constrained_store[:,t-1],
            C_min_staple_constrained_store[:,t-1],TC_S_c3_constrained_store[:,t-1],
            x_S_c3_constrained_store[:,t-1],q_S_c3_constrained_store[:,t-1],c_S_c3_constrained_store[:,t-1],capital_trans[t],
            balanced_share,τ_W,C_grid_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,P_kron_fine,Phi_z_fine,β,fspace_a,σ,P_kron1,fspace_a_fine,a_D,b_D,R,δ,α,cons_level_substinence, epsilon_u, epsilon_r,cttilde);
            print( t-1,",")
        end
    end
    return (residual_store,distr_store,coeff_store,prod_staple_store,prod_cashcrop_store,prod_manuf_store,asset_supply_store,current_worker_pop_store,current_staple_pop_store,current_cashcrop_pop_store,
        current_account_residual_store,staple_productivity_store,cashcrop_productivity_store,manuf_productivity_store,aggregate_consumption_store,relative_land_to_cashcrop_store,
        mean_land_share_staples_store,undernourished_store,fertilizer_use_store,APG_store,var_APland_store,var_MPX_store,avg_labor_prod_rural_store,avg_labor_prod_urban_store,avg_agri_prod_rural_store,
        avg_agri_prod_urban_store,V_saved_store,a_prime_fine_store)
end
function Residual_transition_iterative_epsilon(vec_price_trans::Array{Float64,1},capital_trans::Array{Float64,1},distr_store::Array{Float64,2},
    T::Int64,parameter_end::Parameter_type,coeff_store::Array{Float64,3},τ_trans::Array{Float64,1},cons_level_substinence::Float64,epsilon_u::Float64,epsilon_r::Float64)
    price_trans_actual = reshape(vec_price_trans,4,T);
    # price_trans_actual = copy(price_trans_next)
    residual_store,distr_store,coeff_store =  Residual_transition_sequential_epsilon(price_trans_actual,capital_trans,distr_store,T,parameter_end,coeff_store,τ_trans,moments[1],cons_level_substinence,epsilon_u,epsilon_r)
    maxval, maxval_ind = findmax(abs.(residual_store));
    println("Max residual: ",maxval," in time period: ",maxval_ind[2]," for market: ",maxval_ind[1])
    max_T = maxval_ind[2];
    price_trans_next = zeros(4,T);
    price_trans_next_maxT = copy(price_trans_actual);
    price_trans_next_single_update_plus = copy(price_trans_actual);
    res_cc = residual_store[1,:];
    maxval_cc, maxval_ind_cc = findmax(abs.(res_cc));
    res_tmp = res_cc.>0.2;
    res_cc[res_tmp] .= 0.2;
    res_tmp = res_cc.<-0.2;
    res_cc[res_tmp] .= -0.2;
    
    res_manu = residual_store[2,:];
    maxval_manu, maxval_ind_manu = findmax(abs.(res_manu));
    res_tmp = res_manu.>0.2;
    res_manu[res_tmp] .= 0.2;
    res_tmp = res_manu.<-0.2;
    res_manu[res_tmp] .= -0.2;
    res_tax = residual_store[6,:];
    maxval_tax, maxval_ind_tax = findmax(abs.(res_tax));
    res_tmp = res_tax.>0.2;
    res_tax[res_tmp] .= 0.2;
    res_tmp = res_tax.<-0.2;
    res_tax[res_tmp] .= -0.2;
    res_staple = residual_store[5,:];
    res_tmp = res_staple.>0.2;
    res_staple[res_tmp] .= 0.2;
    res_tmp = res_staple.<-0.2;
    res_staple[res_tmp] .= -0.2;

    res_cttilde = residual_store[4,:];
    res_tmp = res_cttilde.>0.2;
    res_cttilde[res_tmp] .= 0.2;
    res_tmp = res_cttilde.<-0.2;
    res_cttilde[res_tmp] .= -0.2;


    cross_terms = zeros(4,5); # This is where you can meddle with cross terms, for now, only contemporenous prices of the same markets matter
    cross_terms[1,1] =1;
    cross_terms[2,2] =1;
    cross_terms[3,3] =-1;
    cross_terms[1,4] =-1.0;
    cross_terms[2,4] =-1.0;
    cross_terms[3,4] = 1.0;
    cross_terms[4,5] =1.0;
    # 4th column contains the staple market, which I add with a negative sign

    price_trans_next[1,:] = (ones(T)+cross_terms[1,1]*res_cc+cross_terms[1,2]*res_manu +
    cross_terms[1,3]*res_tax + cross_terms[1,4] * res_staple + cross_terms[1,5] * res_cttilde).*price_trans_actual[1,:];
    price_trans_next[2,:] =  (ones(T)+cross_terms[2,1]*res_cc+cross_terms[2,2]*res_manu +
    cross_terms[2,3]*res_tax+ cross_terms[2,4] * res_staple + cross_terms[2,5] * res_cttilde).*price_trans_actual[2,:];
    price_trans_next[3,:] =  (ones(T)+cross_terms[3,1]*res_cc+cross_terms[3,2]*res_manu +
    cross_terms[3,3]*res_tax+ cross_terms[3,4] * res_staple + cross_terms[3,5] * res_cttilde).*price_trans_actual[3,:];
    price_trans_next[4,:] =  (ones(T)+cross_terms[4,1]*res_cc+cross_terms[4,2]*res_manu +
    cross_terms[4,3]*res_tax+ cross_terms[4,4] * res_staple + cross_terms[4,5] * res_cttilde).*price_trans_actual[4,:];
    residual_adj = zeros(6,T);
    residual_adj[1,:] = res_cc;
    residual_adj[2,:] = res_manu;
    residual_adj[4,:] = res_cttilde;
    residual_adj[5,:] = res_tax;
    price_trans_next_maxT[:,1:max_T] = price_trans_next[:,1:max_T];
    if maxval_cc== maxval
        price_trans_next_single_update_plus[1,maxval_ind_cc] = price_trans_next[1,maxval_ind_cc];
        #print(1)
    elseif maxval_manu== maxval
        price_trans_next_single_update_plus[2,maxval_ind_manu] = price_trans_next[2,maxval_ind_manu];
        #print(2)
    elseif maxval_tax== maxval
        price_trans_next_single_update_plus[3,maxval_ind_tax] = price_trans_next[3,maxval_ind_tax];
        #print(3)
    else
        price_trans_next_single_update_plus[:,max_T] = price_trans_next[:,max_T];
        #print(4)
    end
    vec_price_trans_next_maxT = reshape(price_trans_next_maxT,4*T);
    vec_price_trans_next = reshape(price_trans_next,4*T);
    vec_price_trans_next_single_update_plus = reshape(price_trans_next_single_update_plus,4*T);
    return vec_price_trans_next,vec_price_trans_next_maxT,vec_price_trans_next_single_update_plus,residual_store,distr_store,coeff_store
end