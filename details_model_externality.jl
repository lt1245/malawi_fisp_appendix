function details_model_externality(prices::Vector{Float64}, parameters_tmp::Parameter_type, out::Int64, moments::Vector{Float64}, balanced_share::Float64, foreign_supply_capital::Float64, cons_level_substinence::Float64, epsilon_u::Float64, epsilon_r::Float64)

    ctilde = copy(cons_level_substinence)

    if size(prices,1)<3 && balanced_share>0.0 
        error("Prices input wrong sized")
    end

    if balanced_share>0.0 && (prices[1]<0.0 || prices[2]<0.0 || prices[3]<0.0)

        if prices[1]<0.0
            residual_goods = [prices[1]^2*10000 + 1.00,prices[2]^2*1 + 1.00,
                            prices[3]^2*1 + 1.00,1.0^2*1 + 1.00,1.0, 100.00];
        elseif prices[2]<0.0
            residual_goods = [prices[1]^2*1 + 1.00,prices[2]^2*10000 + 1.00,
                            prices[3]^2*1 + 1.00,1.0^2*1 + 1.00,1.0, 100.00];
        elseif prices[3]<0.0
            residual_goods = [prices[1]^2*1 + 1.00,prices[2]^2*1 + 1.00,
                            prices[3]^2*10000 + 1.00,1.0^2*1 + 1.00,1.0, 100.00];
        end

        return residual_goods

    elseif balanced_share==0.0  && (prices[1]<0.0 || prices[2]<0.0)

        if prices[1]<0.0
            residual_goods = [prices[1]^2*10000 + 100.00,prices[2]^2*1 + 100.00,
                            1.00,1.00,1.0];
        elseif prices[2]<0.0
            residual_goods = [prices[1]^2*1 + 100.00,prices[2]^2*10000 + 100.00,
                            1.00,1.00,1.0];
        end
        return residual_goods

    else
        r = moments[1];
        L = moments[2];
        K_Y_ratio = moments[3];
        Y_agr_Y_tot = moments[4];
        exported_cash_crop_ratio = moments[5];
        fraction_staple_producers_without_surplus = moments[6];
        G_Y_ratio = moments[7];
        RCT_moment1_value = moments[8]; #RCT stuff
        RCT_moment2_share = moments[9];
        exp_ratio=moments[10];
        RU_migration_rate=moments[11];
        fraction_only_Bashcrops=moments[12];
        mean_land_share_to_staples_among_cc=moments[13];
        rural_pop_only_staples=moments[14];
        urban_rural_inc_ratio=moments[15];
        urban_rural_wealth_ratio=moments[16];
        urban_rural_consumption_ratio=moments[17];
        Top1_share_wealth_rural=moments[18];
        Top1_share_income_rural=moments[19];
        Top1_share_consumption_rural=moments[20];
        Top10_share_wealth_rural=moments[21];
        Top10_share_income_rural=moments[22];
        Top10_share_consumption_rural=moments[23];
        Top1_share_wealth_urban=moments[24];
        Top1_share_income_urban=moments[25];
        Top1_share_consumption_urban=moments[26];
        Top10_share_wealth_urban=moments[27];
        Top10_share_income_urban=moments[28];
        Top10_share_consumption_urban=moments[29];
        RCT_moment3_share=moments[30];
        RCT_moment4_increase=moments[31];
        fraction_borrowing_data =moments[32];
        fraction_cashcrop_suboptimal = moments[33]; #from Brune et al EDCC 2016 Footnote 9.
        APG_data= moments[34];
        # var_log_cons_rural=moments[24];
        # var_log_cons_urban=moments[25];
        # var_log_inc_rural=moments[26];
        # var_log_inc_urban=moments[27];
        # var_log_wealth_rural=moments[28];
        # var_log_wealth_urban=moments[29];
        exitflag = 0;

        if balanced_share>0.0
            cttilde=prices[4]
        else
            cttilde=prices[3]
        end

        # Initialization
        (δ,ζ,ρ,α,σ,β,ϵ,ψ_S,ψ_B,ψ_M,ϕ_S,ϕ_B,c̄_S,F_W,F_S,F_B,FM_W,FM_S,FM_B,Q_S,p_x,τ_S,τ_B,a_D,b_D,K_a,K_b,γ,A_W,
        ρ_S,ρ_SW,σ_S,ρ_W,σ_W,n,n_fine,agrid,agrid_fine,a_min,a_max,spliorder,fspace_a,fspace_a_fine,fspace_C_fine,C_grid_fine_no,C_grid_fine,s,ns,
            s_fine, ns_fine, z, z_W, Phi_z, Phi_z_fine, Phi, Phi_aug, P_kron, P_kron1, P_kron_fine, κ) = local_parameters_ext(parameters_tmp, epsilon_r,cttilde, ctilde)

        p_B,p_M,R,r,w,τ_W= price_reshaper_fixed_r(prices,δ,ϵ,ψ_S,ψ_B,ψ_M,p_x,τ_S,τ_B,α,balanced_share,r);



        #println("Trying p_B: ", p_B, " p_M: ", p_M, " r: ", r, " w: ", w)

        (θ,labor_prod,tol,P_W,Y_W,coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
                coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,λ_B_interior_c3a,x_SC_interior_c3a,x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,
                q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b, c_S_mat,c_B_mat,
                c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat, λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,q_S_staples,c_S_staples,c_B_staples,
                c_M_staples,P_S_staples,x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,C_max_unconstrained ,C_max_constrained,C_min_unconstrained,
                C_min_constrained,TC_mat,C_max_staple,C_min_staple,C_max_staple_constrained,C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,
                q_S_c3_constrained,c_S_c3_constrained,cbar_violated, x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c) = income_creator_ext(s,ns,z,z_W,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,
            c̄_S, a_min, a_max, γ, n, κ, Q_S, ϵ, ψ_S, ψ_B, ψ_M, C_grid_fine, fspace_C_fine, agrid, epsilon_u, cttilde, ctilde)
                #println("cbar: ", cbar_violated)

        #if cbar_violated==0
        min_C_applied,max_C_applied = bounds_consumption(P_W,Y_W,s,r,ρ,w,
            coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,θ,
                fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns,κ,tol,a_min,a_max,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
                    coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,
                    λ_B_interior_c3a,x_SC_interior_c3a,x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,
                    x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b, c_S_mat,c_B_mat,
                    c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat, λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,q_S_staples,c_S_staples,c_B_staples,
                    c_M_staples,P_S_staples,x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,F_W,F_S,F_B,
                    FM_W,FM_S,FM_B,TC_mat,C_max_staple,C_min_staple,C_max_staple_constrained,C_min_staple_constrained,TC_S_c3_constrained,
                    x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained, x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c);




        coeff = zeros(ns,3);
        coeff_next = zeros(ns,3);
        x_tmp = zeros(ns,3);
        V_tmp = zeros(ns,3);
        V_next_stacked = zeros(ns*3,1);
        iterate1 = 0;
        conv = 10.0
        Phi_prime_tmp = Array{SparseMatrixCSC{Float64,Int64}}(undef, 3);
        D_deriv_tmp_block = Array{SparseMatrixCSC{Float64,Int64}}(undef, 3,1);
        Q_trans = spzeros(ns_fine*3,ns_fine*3);
        for i = 1:3
        #while conv> 10^(-7)
        (coeff[:], conv, conv_ind)  = Bellman_iteration(coeff,coeff_next,x_tmp,V_tmp,P_kron,Phi,min_C_applied,max_C_applied,Phi_z,β,
            fspace_a,σ,P_W,Y_W,s,r,ρ,w,coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,θ,
            fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns,κ,tol,a_min,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
            coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,
            λ_B_interior_c3a,x_SC_interior_c3a,x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,
            x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b, c_S_mat,c_B_mat,
            c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat, λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,q_S_staples,c_S_staples,c_B_staples,
            c_M_staples,P_S_staples,x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat,a_max,C_max_staple,C_min_staple,C_max_staple_constrained,
            C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained,x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c);
            #println("Conv criterion:", conv, "Max error index:", conv_ind)
        end
        while conv> 10^(-7)
        (coeff[:], conv,iterate1,exitflag,conv_ind) =     Newton_iteration(coeff,x_tmp,V_tmp,P_kron,Phi,min_C_applied,max_C_applied,Phi_z,β,
                fspace_a,σ,P_W,Y_W,s,r,ρ,w,coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,θ,
                fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns,κ,tol,a_min,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
                coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,
                λ_B_interior_c3a,x_SC_interior_c3a,x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,
                x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b, c_S_mat,c_B_mat,
                c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat, λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,q_S_staples,c_S_staples,c_B_staples,
                c_M_staples,P_S_staples,x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat,a_max,
                V_next_stacked,iterate1,Phi_prime_tmp,D_deriv_tmp_block,Phi_aug,P_kron1,exitflag,C_max_staple,C_min_staple,C_max_staple_constrained,
                C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained,x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c);
            if iterate1>20
                println("Conv criterion:", conv, "Max error index:", conv_ind, prices)
                conv = 0;
                exitflag=4;
            end
                #println("Conv criterion:", conv, "Max error index:", conv_ind,coeff[conv_ind[2][1],conv_ind[2][2]])
        end
        #end
        if balanced_share>0.0
            residual_goods = [p_B^2*100 + 1.00,p_M^2*100 + 1.00,
            r^2*100 + 1.00,w^2*100 + 1.00,100,τ_W^2*100 + 1.00,ctilde^2*100 + 1.00]; # Placeholder if everything is fine, this is overwritten
            tmp_var=ones(7);
        else
            residual_goods = [p_B^2*100 + 1.00,p_M^2*100 + 1.00,
                r^2 * 100 + 1.00, w^2 * 100 + 1.00, 100, ctilde^2 * 100 + 1.00] # Placeholder if everything is fine, this is overwritten
            tmp_var=ones(6);
        end

        if exitflag==4
            residual_goods = 100.0*tmp_var;
        elseif cbar_violated==1
            residual_goods = 200.0*tmp_var;
        else
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
                    z,z_W,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n_fine,κ,Q_S,ϵ,ψ_S,ψ_B,ψ_M,coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,
                    C_max_unconstrained ,C_max_constrained,C_min_unconstrained,C_min_constrained, coeff_λ_2_s,C_grid_fine,fspace_C_fine,C_max_staple,C_min_staple,C_max_staple_constrained,
                C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained, epsilon_u, cttilde, ctilde);


            min_C_applied_fine,max_C_applied_fine = bounds_consumption(P_W_fine,Y_W_fine,s_fine,r,ρ,w,
            coeff_λ_2_cashcrop_residual_unconstrained_fine,coeff_λ_2_cashcrop_residual_constrained_fine,θ_fine,
            fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,
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

            (Q_trans_prime,cons_fine_local,a_prime_fine_local,future_occupation_fine_local,V_saved_local) = Q_transition(coeff,ns_fine,P_kron_fine,min_C_applied_fine,max_C_applied_fine,Phi_z_fine,
                β,fspace_a,σ,P_W_fine,Y_W_fine,s_fine,r,ρ,w,coeff_λ_2_cashcrop_residual_unconstrained_fine,coeff_λ_2_cashcrop_residual_constrained_fine,θ_fine,
                fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,κ,tol,a_min,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
                coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,labor_allocated_interior_c3a_fine,λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,
                x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,q_S_c2_fine,q_B_c1_fine,q_S_c3a_fine,q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,
                x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,Y_B_c3b_fine,c_S_mat_fine,c_B_mat_fine,c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,
                q_B_mat_fine,land_B_mat_fine,λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,feasibility_mat_fine,C_max_mat_fine,C_min_mat_fine,q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,
                c_M_staples_fine,P_S_staples_fine,x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat_fine,a_max,P_kron1,Q_trans,
                C_max_staple_fine,C_min_staple_fine,C_max_staple_constrained_fine,
                C_min_staple_constrained_fine,TC_S_c3_constrained_fine,x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine,fspace_a_fine,
                x_S_mat_3c_fine,x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine);

            exitflag = predict_irreducibility(future_occupation_fine_local,exitflag);


            stat_distr,exitflag_tmp = stationary_distribution(Q_trans_prime,ns_fine,exitflag);

            if exitflag_tmp ==5
                if balanced_share>0.0
                    residual_goods = [p_B^2*100 + 1.00,p_M^2*100 + 1.00,
                    r^2*100 + 1.00,w^2*100 + 1.00,100,τ_W^2*100 + 1.00]; # Placeholder if everything is fine, this is overwritten
                else
                    residual_goods = [p_B^2*100 + 1.00,p_M^2*100 + 1.00,
                    r^2*100 + 1.00,w^2*100 + 1.00,100]; # Placeholder if everything is fine, this is overwritten
                end
            else
                # Start acquiring objects for the equilibrium
                # Distribution of past occupations
                worker_past_dist = stat_distr[(ns_fine *0 + 1):(ns_fine *1)];
                staple_past_dist = stat_distr[(ns_fine *1 + 1):(ns_fine *2)];
                cash_crop_past_dist = stat_distr[(ns_fine *2 + 1):(ns_fine *3)];

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
                        Y_W_fine_policy[:,j] =  P_W *cons_fine_local[:,j]  + Y_W_fine  .+ w*FM_W .+ w*F_W * (j != 1);
                    end
                        if jj == 2
                            (future_asset_S[:,j],Y_S_fine[:,j],c_S_S_fine[:,j],c_B_S_fine[:,j],c_M_S_fine[:,j],q_S_S_fine[:,j],
                            P_S_fine[:,j],x_S_S_fine[:,j],solve_staple_index_S_fine[:,j],λ_2_S_fine[:,j]) = policy_function_creator(
                            cons_fine_local[:,j],jj,j,P_W_fine,Y_W_fine,s_fine,r,ρ,w,coeff_λ_2_cashcrop_residual_unconstrained_fine,
                            coeff_λ_2_cashcrop_residual_constrained_fine,θ_fine,fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,
                            p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,a_min,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
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
                            coeff_λ_2_cashcrop_residual_constrained_fine,θ_fine,fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,
                            p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,a_min,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
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


                past_distr_store[:,1]=stat_distr[(ns_fine *0 + 1):(ns_fine *1)];
                past_distr_store[:,2]=stat_distr[(ns_fine *1 + 1):(ns_fine *2)];
                past_distr_store[:,3]=stat_distr[(ns_fine *2 + 1):(ns_fine *3)];

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


                if balanced_share>0.0
                    τ_W_new = - balanced_share * p_x * (τ_S * input_staple + τ_B * input_cashcrop) / (labor_used * w);
                end


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


                #residual[1] = (cash_crop_cons + p_x/p_B*input_staple + p_x/p_B*input_cashcrop - prod_cashcrop)/(cash_crop_cons + p_x/p_B*input_staple + p_x/p_B*input_cashcrop + prod_cashcrop); # -- KAROL ADDED THIS
                #residual[1] = 0#(cash_crop_cons + p_x*input_staple + p_x*input_cashcrop - prod_cashcrop)/(cash_crop_cons + p_x*input_staple + p_x*input_cashcrop + prod_cashcrop); # -- KAROL ADDED THIS
                #residual[1] = (cash_crop_cons  - prod_cashcrop)/(cash_crop_cons + prod_cashcrop); # -- KAROL ADDED THIS
                residual_goods[1] = (c_B_worker_sum + c_B_staple_sum + c_B_Cashcrop_sum+ foreign_demand_cash_crop  - prod_cashcrop)/(
                                    c_B_worker_sum + c_B_staple_sum + c_B_Cashcrop_sum+ foreign_demand_cash_crop  + prod_cashcrop);
                if current_cashcrop_pop <1e-3
                    residual_goods[1] = 100.0;
                end
                #println(c_B_worker_sum + c_B_staple_sum + c_B_Cashcrop_sum)
                #println(foreign_demand_cash_crop)
                #println(prod_cashcrop)

                residual_goods[2] = (c_M_worker_sum + c_M_staple_sum + c_M_cashcrop_sum - prod_manuf)/(c_M_worker_sum + c_M_staple_sum + c_M_cashcrop_sum + prod_manuf);

                #residual[2] = (manuf_cons + input_staple + input_cashcrop - prod_manuf)/(manuf_cons + input_staple + input_cashcrop + prod_manuf); # -- KAROL ADDED THIS
                residual_goods[3] = 0 #No longer needed for the calibration as the foreign_supply_capital offsets capital markets(r - r_new)/(r+r_new);
                residual_goods[4] = 0;#(w - w_new)/(w+w_new);
                if current_worker_pop <1e-3
                    residual_goods[4] = -100.0;
                end
                #residual[4] = (labor_demand - labor_used)/(labor_demand + labor_used);
                #staple_mkt_clr = (sum(c_current_S) + p_x*input_staple - prod_staple)/(sum(c_current_S) + prod_staple + p_x*input_staple) # -- KAROL ADDED THIS
                #staple_mkt_clr = (sum(c_current_S) - prod_staple)/(sum(c_current_S) + prod_staple) # -- KAROL ADDED THIS
                residual_goods[5] = (c_S_worker_sum + c_S_staple_sum + c_S_cashcrop_sum +
                    transaction_cost_worker_sum + transaction_cost_staple_sum + transaction_cost_cashcrop_sum - prod_staple)/(
                    c_S_worker_sum + c_S_staple_sum + c_S_cashcrop_sum +
                        transaction_cost_worker_sum + transaction_cost_staple_sum + transaction_cost_cashcrop_sum + prod_staple);# This is the main residual, since other markets are tempered with
                if prod_staple <1e-3
                    residual_goods[5] = 100.0;
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

                if balanced_share>0.0
                    residual_goods[6] =  (τ_W - τ_W_new)/(τ_W + τ_W_new)

                    residual_goods[7] = (cttilde_new - cttilde)/(cttilde_new+cttilde)
                else
                    residual_goods[6] = (cttilde_new - cttilde) / (cttilde_new + cttilde)
                end

               

                # if current_staple_pop ==0
                #     residual[5] = 100.0;
                # end

                rural_pop_only_staples_printable=(sum(entrants_cashcrop_from_workers.* (q_B_B_fine[:,1].==0.0)) + sum(incumbents_cashcrop.*(q_B_B_fine[:,3].==0.0))
                + sum( entrants_from_staple_to_cashcrop.*(q_B_B_fine[:,2].==0.0)) + current_staple_pop);
                rural_pop_only_Bashcrop_printable=(sum(entrants_cashcrop_from_workers.* (land_C_fine[:,1].==1.0)) + sum(incumbents_cashcrop.*(land_C_fine[:,3].==1.0))
                + sum( entrants_from_staple_to_cashcrop.*(land_C_fine[:,2].==1.0)));
                #println("Pops: Worker:",current_worker_pop," Staple:",current_staple_pop," Cash crop:",current_cashcrop_pop," Only producing staples", rural_pop_only_staples_printable
                #," Only producing cashcrop", rural_pop_only_Bashcrop_printable)
            end
        end
    #end

    # if current_staple_pop ==0
    #     residual[5] = 100.0;
    # end

    rural_pop_only_staples_printable=(sum(entrants_cashcrop_from_workers.* (q_B_B_fine[:,1].==0.0)) + sum(incumbents_cashcrop.*(q_B_B_fine[:,3].==0.0))
    + sum( entrants_from_staple_to_cashcrop.*(q_B_B_fine[:,2].==0.0)) + current_staple_pop);
    rural_pop_only_Bashcrop_printable=(sum(entrants_cashcrop_from_workers.* (land_C_fine[:,1].==1.0)) + sum(incumbents_cashcrop.*(land_C_fine[:,3].==1.0))
    + sum( entrants_from_staple_to_cashcrop.*(land_C_fine[:,2].==1.0)));
    #println("Pops: Worker:",current_worker_pop," Staple:",current_staple_pop," Cash crop:",current_cashcrop_pop," Only producing staples", rural_pop_only_staples_printable
    #," Only producing cashcrop", rural_pop_only_Bashcrop_printable)
            #Government_expenditure = balanced_share * p_x * (τ_S * input_staple + τ_B * input_cashcrop);
    #Net_Foreign_Factor_Income = R*(copy(asset_supply) - capital_used)
    current_account = (Government_expenditure + Net_Foreign_Factor_Income - Import_value+Export_value)
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
    # moment[5]
    # Targeting subsistence  - maybe focus on consumption to income ratios?
    fraction_model = sum(entrants_cashcrop_from_workers .*((c_S_B_fine[:,1] .- q_S_B_fine[:,1].>=0.0) .& (q_S_B_fine[:,1].>0.0)) +
                        incumbents_cashcrop.*((c_S_B_fine[:,3] .- q_S_B_fine[:,3].>=0.0) .& (q_S_B_fine[:,3].>0.0)) +
                        entrants_from_staple_to_cashcrop.*((c_S_B_fine[:,2] .- q_S_B_fine[:,2].>=0.0) .& (q_S_B_fine[:,2].>0.0)) +
                        entrants_staple_from_workers.*(c_S_S_fine[:,1] .- q_S_S_fine[:,1].>=0.0) +
                        incumbents_staple.*(c_S_S_fine[:,2] .- q_S_S_fine[:,2].>=0.0) +
                        exit_cashcrop_to_staple.*(c_S_S_fine[:,3] .- q_S_S_fine[:,3].>=0.0)) / (current_staple_pop +
                            sum(entrants_cashcrop_from_workers .*(q_S_B_fine[:,1].>0.0) +
                                incumbents_cashcrop.*(q_S_B_fine[:,3].>0.0) +
                                entrants_from_staple_to_cashcrop.*(q_S_B_fine[:,2].>0.0)));

    program_spending = -p_x * (τ_S * input_staple + τ_B * input_cashcrop)/nominal_GDP;


        # share_urban_unemployed = sum((current_workers .* (labor_prod_fine .== parameters_tmp.l_z_low)) ./current_worker_pop);
        # residual[7] = (share_urban_unemployed - Urban_unemp_rate)/(share_urban_unemployed + Urban_unemp_rate);
    # moment[8]

        #for res[9] with RCT stuff here....
    # RCT stuff - but if I am here, lets get other poverty measures too.
    # Modified: Laszlo 29/08/2022
    total_nominal_consumption = (p_B * (c_B_worker_sum + c_B_staple_sum + c_B_Cashcrop_sum) +
            (c_S_worker_sum + c_S_staple_sum + c_S_cashcrop_sum + transaction_cost_staple_sum + transaction_cost_cashcrop_sum + transaction_cost_worker_sum) +
        p_M * (c_M_worker_sum + c_M_staple_sum + c_M_cashcrop_sum));
    #total_nominal_consumption = (p_B * (c_B_worker_sum + c_B_staple_sum + c_B_Cashcrop_sum) +
    #    (1 + Q_S) * (c_S_worker_sum + c_S_staple_sum + c_S_cashcrop_sum) +
    #    p_M * (c_M_worker_sum + c_M_staple_sum + c_M_cashcrop_sum));
    mean_cons_tmp = total_nominal_consumption/ 1.0;

    income_worker_stay_workers = labor_prod_fine.*w;# .- w*FM_W #-Y_W_fine_policy[:,1];
    income_worker_staple_to_work = labor_prod_fine.*w;# .- w*FM_W .- w*F_W #-Y_W_fine_policy[:,2];
    income_worker_cashcrop_to_work = labor_prod_fine.*w;# .- w*FM_W .- w*F_W #-Y_W_fine_policy[:,3];
    # Worker policy functions: .*stay_workers .*exit_staple_to_work  .*exit_cashcrop_to_work
    #Staples
    c_S_worker_stay_workers = (c̄_S .+ (1 + Q_S)^(-ϵ)*ψ_S.^ϵ.*cons_fine_local[:,1].*P_W_fine.^ϵ)
    c_S_worker_staple_to_work = (c̄_S .+ (1 + Q_S)^(-ϵ)*ψ_S.^ϵ.*cons_fine_local[:,2].*P_W_fine.^ϵ)
    c_S_worker_cashcrop_to_work = (c̄_S .+ (1 + Q_S)^(-ϵ)*ψ_S.^ϵ.*cons_fine_local[:,3].*P_W_fine.^ϵ);
    #Cashcrop
    c_B_worker_stay_workers = (p_B^(-ϵ)*ψ_B.^ϵ.*cons_fine_local[:,1].*P_W_fine.^ϵ)
    c_B_worker_staple_to_work = (p_B^(-ϵ)*ψ_B.^ϵ.*cons_fine_local[:,2].*P_W_fine.^ϵ)
    c_B_worker_cashcrop_to_work = (p_B^(-ϵ)*ψ_B.^ϵ.*cons_fine_local[:,3].*P_W_fine.^ϵ);
    #Manufacturing
    c_M_worker_stay_workers = (p_M^(-ϵ)*ψ_M.^ϵ.*cons_fine_local[:,1].*P_W_fine.^ϵ)
    c_M_worker_staple_to_work = (p_M^(-ϵ)*ψ_M.^ϵ.*cons_fine_local[:,2].*P_W_fine.^ϵ)
    c_M_worker_cashcrop_to_work = (p_M^(-ϵ)*ψ_M.^ϵ.*cons_fine_local[:,3].*P_W_fine.^ϵ);
    #Worker # consumption
    worker_consumption_stay_workers = (1 + Q_S) *c_S_worker_stay_workers  + p_B*c_B_worker_stay_workers + p_M*c_M_worker_stay_workers;
    worker_consumption_staple_to_work = (1 + Q_S) *c_S_worker_staple_to_work  + p_B*c_B_worker_staple_to_work + p_M*c_M_worker_staple_to_work
    worker_consumption_cashcrop_to_work =  (1 + Q_S) *c_S_worker_cashcrop_to_work + p_B*c_B_worker_cashcrop_to_work + p_M*c_M_worker_cashcrop_to_work;
    # This is not cons_fine_local! We have to subtract .- c̄_S from c_S_worker to have equivalence
    #income:
    worker_gross_income_stay_workers =copy(income_worker_stay_workers)# (income_worker_stay_workers + worker_consumption_stay_workers);
    worker_gross_income_staple_to_work = copy(income_worker_staple_to_work)#(income_worker_staple_to_work + worker_consumption_staple_to_work);
    worker_gross_income_cashcrop_to_work =copy(income_worker_cashcrop_to_work)# (income_worker_cashcrop_to_work + worker_consumption_cashcrop_to_work);
    income_worker_mean = sum(worker_gross_income_stay_workers .* stay_workers
                            +worker_gross_income_staple_to_work .* exit_staple_to_work
                            +worker_gross_income_cashcrop_to_work .* exit_cashcrop_to_work)/current_worker_pop;

    # Staple farmers policy functions: .*stay_workers .*exit_staple_to_work  .*exit_cashcrop_to_work
    #Staple consumption

    test_mat_staple = c_S_S_fine .- q_S_S_fine .>= 0
    staple_market_consumption_staple_from_workers = c_S_S_fine[:,1] + Q_S .* test_mat_staple[:,1] .* (c_S_S_fine[:,1] - q_S_S_fine[:,1]) + p_B*c_B_S_fine[:,1] + p_M*c_M_S_fine[:,1]; # was: (1 + Q_S) *c_S_S_fine[:,1] + p_B*c_B_S_fine[:,1] + p_M*c_M_S_fine[:,1];
    staple_market_consumption_incumbents_staple = c_S_S_fine[:,2] + Q_S .* test_mat_staple[:,2] .* (c_S_S_fine[:,2] - q_S_S_fine[:,2]) + p_B*c_B_S_fine[:,2] + p_M*c_M_S_fine[:,2]; #was: (1 + Q_S) *c_S_S_fine[:,2] + p_B*c_B_S_fine[:,2] + p_M*c_M_S_fine[:,2];
    staple_market_consumption_cashcrop_to_staple = c_S_S_fine[:,3] + Q_S .* test_mat_staple[:,3] .* (c_S_S_fine[:,3] - q_S_S_fine[:,3]) + p_B*c_B_S_fine[:,3] + p_M*c_M_S_fine[:,3]; #was: (1 + Q_S) *c_S_S_fine[:,3] + p_B*c_B_S_fine[:,3] + p_M*c_M_S_fine[:,3];

    # income
    #I think fertilizer exp is already calculated in Y, and we correctly subtract consumption (i.e. we want to have pure profits net of input expenditures and fixed costs):
    revenue_staples = max.(q_S_S_fine - c_S_S_fine,0);
    staple_revenue_income_staple_from_workers = copy(revenue_staples[:,1]); #.- w*FM_S .- w*F_S);
    staple_revenue_income_incumbents_staple = copy(revenue_staples[:,2]); #.- w*FM_S);
    staple_revenue_income_cashcrop_to_staple =copy(revenue_staples[:,3]); #.- w*FM_S .- w*F_S);
    income_staple_mean = sum(staple_revenue_income_staple_from_workers .* entrants_staple_from_workers
    +staple_revenue_income_incumbents_staple.* incumbents_staple
    +staple_revenue_income_cashcrop_to_staple .* exit_cashcrop_to_staple)/current_staple_pop;
    #Cash crop

    test_mat_cashcrop = c_S_B_fine .- q_S_B_fine .>= 0 #similar changes as above
    cashcrop_market_consumption_cashcrop_from_workers = c_S_B_fine[:,1] + Q_S .* test_mat_cashcrop[:,1] .* (c_S_B_fine[:,1] - q_S_B_fine[:,1])  + p_B*c_B_B_fine[:,1] + p_M*c_M_B_fine[:,1]; # was: (1 + Q_S) *c_S_B_fine[:,1] + p_B*c_B_B_fine[:,1] + p_M*c_M_B_fine[:,1];
    cashcrop_market_consumption_staple_to_cashcrop = c_S_B_fine[:,2] + Q_S .* test_mat_cashcrop[:,2] .* (c_S_B_fine[:,2] - q_S_B_fine[:,2]) + p_B*c_B_B_fine[:,2] + p_M*c_M_B_fine[:,2];#(1 + Q_S) *c_S_B_fine[:,2] + p_B*c_B_B_fine[:,2] + p_M*c_M_B_fine[:,2];
    cashcrop_market_consumption_incumbents_cashcrop = c_S_B_fine[:,3] + Q_S .* test_mat_cashcrop[:,3] .* (c_S_B_fine[:,3] - q_S_B_fine[:,3]) + p_B*c_B_B_fine[:,3] + p_M*c_M_B_fine[:,3];#(1 + Q_S) *c_S_B_fine[:,3] + p_B*c_B_B_fine[:,3] + p_M*c_M_B_fine[:,3];

    #income:
    revenue_cashcrop = max.(q_S_B_fine - c_S_B_fine,0) + p_B*max.(q_B_B_fine - c_B_B_fine,0) # +
    cashcrop_revenue_income_cashcrop_from_workers =copy(revenue_cashcrop[:,1]);# (-Y_B_fine[:,1] + cashcrop_consumption_cashcrop_from_workers);# .- w*FM_B .- w*F_B);
    cashcrop_revenue_income_staple_to_cashcrop = copy(revenue_cashcrop[:,2]);#(-Y_B_fine[:,2] + cashcrop_consumption_staple_to_cashcrop);# .- w*FM_B .- w*F_B);
    cashcrop_revenue_income_incumbents_cashcrop = copy(revenue_cashcrop[:,3]);#(-Y_B_fine[:,3] + cashcrop_consumption_incumbents_cashcrop);# .- w*FM_B);

    income_cashcrop_mean = sum(cashcrop_revenue_income_cashcrop_from_workers .* entrants_cashcrop_from_workers
    +cashcrop_revenue_income_staple_to_cashcrop.* entrants_from_staple_to_cashcrop
    +cashcrop_revenue_income_incumbents_cashcrop .* incumbents_cashcrop)/current_cashcrop_pop;

    # Get the distribution for the urban
    current_distr_store_urban = copy(past_distr_store);
    current_distr_store_urban[1:ns_fine] = (current_distr_store_urban[1:ns_fine].*(future_occupation_fine_local[:,1] .== 1.0 ));
    current_distr_store_urban[(ns_fine+1):(2*ns_fine)] = (current_distr_store_urban[(ns_fine+1):(2*ns_fine)].*(future_occupation_fine_local[:,2] .== 1.0 ));
    current_distr_store_urban[(2*ns_fine + 1):(3*ns_fine)] =(current_distr_store_urban[(2*ns_fine + 1):(3*ns_fine)].*(future_occupation_fine_local[:,3] .== 1.0 ));
    current_distr_store_urban = current_distr_store_urban./current_worker_pop;
    # Get the distribution for the rural
    current_distr_store_rural = copy(past_distr_store);
    current_distr_store_rural[1:ns_fine] = (current_distr_store_rural[1:ns_fine].*(future_occupation_fine_local[:,1] .> 1.0 ));
    current_distr_store_rural[(ns_fine+1):(2*ns_fine)] = (current_distr_store_rural[(ns_fine+1):(2*ns_fine)].*(future_occupation_fine_local[:,2] .> 1.0 ));
    current_distr_store_rural[(2*ns_fine + 1):(3*ns_fine)] =(current_distr_store_rural[(2*ns_fine + 1):(3*ns_fine)].*(future_occupation_fine_local[:,3] .> 1.0 ));
    current_distr_store_rural = current_distr_store_rural./(1.0 - current_worker_pop);

    # Consumption ordered:
    consumption_policy = zeros(3*ns_fine);
    # consumption_policy_urban = zeros(ns_fine);
    # consumption_policy_rural = zeros(2*ns_fine);
    consumption_policy[1:ns_fine] = (worker_consumption_stay_workers.*(future_occupation_fine_local[:,1] .== 1.0 ) +
    staple_market_consumption_staple_from_workers.*(future_occupation_fine_local[:,1].== 2.0 ) +
    cashcrop_market_consumption_cashcrop_from_workers.*(future_occupation_fine_local[:,1].== 3.0 ));
    consumption_policy[(ns_fine+1):(2*ns_fine)] = (worker_consumption_staple_to_work.*(future_occupation_fine_local[:,2] .== 1.0 ) +
    staple_market_consumption_incumbents_staple.*(future_occupation_fine_local[:,2].== 2.0 ) +
    cashcrop_market_consumption_staple_to_cashcrop.*(future_occupation_fine_local[:,2].== 3.0 ));
    consumption_policy[(2*ns_fine + 1):(3*ns_fine)] =(worker_consumption_cashcrop_to_work.*(future_occupation_fine_local[:,3] .== 1.0 ) +
    staple_market_consumption_cashcrop_to_staple.*(future_occupation_fine_local[:,3].== 2.0 ) +
    cashcrop_market_consumption_incumbents_cashcrop.*(future_occupation_fine_local[:,3].== 3.0 ));
    consumption_sort_index = sortperm(consumption_policy);
    consumption_sorted = consumption_policy[consumption_sort_index];
    current_distr_consumption_sorted = current_distr_store[consumption_sort_index];
    cumu_consumption_distr =  cumsum(current_distr_consumption_sorted,dims = 1);
    # Consumption of urban ordered:
    consumption_policy_urban = zeros(3*ns_fine);
    # consumption_policy_urban = zeros(ns_fine);
    # consumption_policy_rural = zeros(2*ns_fine);
    consumption_policy_urban[1:ns_fine] = (worker_consumption_stay_workers.*(future_occupation_fine_local[:,1] .== 1.0 ));
    consumption_policy_urban[(ns_fine+1):(2*ns_fine)] = (worker_consumption_staple_to_work.*(future_occupation_fine_local[:,2] .== 1.0 ));
    consumption_policy_urban[(2*ns_fine + 1):(3*ns_fine)] =(worker_consumption_cashcrop_to_work.*(future_occupation_fine_local[:,3] .== 1.0 ));
    consumption_sort_index_urban = sortperm(consumption_policy_urban);
    consumption_sorted_urban = consumption_policy_urban[consumption_sort_index_urban];
    current_distr_consumption_sorted_urban = current_distr_store_urban[consumption_sort_index_urban];
    cumu_consumption_distr_urban =  cumsum(current_distr_consumption_sorted_urban,dims = 1);
    # Consumption of rural ordered:
    consumption_policy_rural = zeros(3*ns_fine);
    # consumption_policy_rural = zeros(ns_fine);
    # consumption_policy_rural = zeros(2*ns_fine);
    consumption_policy_rural[1:ns_fine] = (staple_market_consumption_staple_from_workers.*(future_occupation_fine_local[:,1].== 2.0 ) +
    cashcrop_market_consumption_cashcrop_from_workers.*(future_occupation_fine_local[:,1].== 3.0 ));
    consumption_policy_rural[(ns_fine+1):(2*ns_fine)] = (staple_market_consumption_incumbents_staple.*(future_occupation_fine_local[:,2].== 2.0 ) +
    cashcrop_market_consumption_staple_to_cashcrop.*(future_occupation_fine_local[:,2].== 3.0 ));
    consumption_policy_rural[(2*ns_fine + 1):(3*ns_fine)] =(staple_market_consumption_cashcrop_to_staple.*(future_occupation_fine_local[:,3].== 2.0 ) +
    cashcrop_market_consumption_incumbents_cashcrop.*(future_occupation_fine_local[:,3].== 3.0 ));
    consumption_sort_index_rural = sortperm(consumption_policy_rural);
    consumption_sorted_rural = consumption_policy_rural[consumption_sort_index_rural];
    current_distr_consumption_sorted_rural = current_distr_store_rural[consumption_sort_index_rural];
    cumu_consumption_distr_rural =  cumsum(current_distr_consumption_sorted_rural,dims = 1);

    mean_cons = sum(consumption_sorted.*current_distr_consumption_sorted)
    mean_cons_urban = sum(consumption_sorted_urban.*current_distr_consumption_sorted_urban)
    mean_cons_rural = sum(consumption_sorted_rural.*current_distr_consumption_sorted_rural)

    # Income ordered:
    income_policy = zeros(3*ns_fine);
    income_policy_urban = zeros(3*ns_fine);
    income_policy_rural = zeros(3*ns_fine);
    income_policy[1:ns_fine] = (worker_gross_income_stay_workers.*(future_occupation_fine_local[:,1] .== 1.0 ) +
    staple_revenue_income_staple_from_workers.*(future_occupation_fine_local[:,1].== 2.0 ) +
    cashcrop_revenue_income_cashcrop_from_workers.*(future_occupation_fine_local[:,1].== 3.0 ));
    income_policy[(ns_fine+1):(2*ns_fine)] = (worker_gross_income_staple_to_work.*(future_occupation_fine_local[:,2] .== 1.0 ) +
    staple_revenue_income_incumbents_staple.*(future_occupation_fine_local[:,2].== 2.0 ) +
    cashcrop_revenue_income_staple_to_cashcrop.*(future_occupation_fine_local[:,2].== 3.0 ));
    income_policy[(2*ns_fine + 1):(3*ns_fine)] =(worker_gross_income_cashcrop_to_work.*(future_occupation_fine_local[:,3] .== 1.0 ) +
    staple_revenue_income_cashcrop_to_staple.*(future_occupation_fine_local[:,3].== 2.0 ) +
    cashcrop_revenue_income_incumbents_cashcrop.*(future_occupation_fine_local[:,3].== 3.0 ));
    income_policy_urban[1:ns_fine] = (worker_gross_income_stay_workers.*(future_occupation_fine_local[:,1] .== 1.0 ));
    income_policy_urban[(ns_fine+1):(2*ns_fine)] = (worker_gross_income_staple_to_work.*(future_occupation_fine_local[:,2] .== 1.0 ));
    income_policy_urban[(2*ns_fine + 1):(3*ns_fine)] =(worker_gross_income_cashcrop_to_work.*(future_occupation_fine_local[:,3] .== 1.0 ));


    income_policy_rural[1:ns_fine] = (staple_revenue_income_staple_from_workers.*(future_occupation_fine_local[:,1].== 2.0 ) +
    cashcrop_revenue_income_cashcrop_from_workers.*(future_occupation_fine_local[:,1].== 3.0 ));
    income_policy_rural[(ns_fine+1):(2*ns_fine)] = (staple_revenue_income_incumbents_staple.*(future_occupation_fine_local[:,2].== 2.0 ) +
        cashcrop_revenue_income_staple_to_cashcrop.*(future_occupation_fine_local[:,2].== 3.0 ));
    income_policy_rural[(2*ns_fine + 1):(3*ns_fine)] =(staple_revenue_income_cashcrop_to_staple.*(future_occupation_fine_local[:,3].== 2.0 ) +
    cashcrop_revenue_income_incumbents_cashcrop.*(future_occupation_fine_local[:,3].== 3.0 ));

    #worker_gross_income_policy = (worker_gross_income_stay_workers.*(future_occupation_fine_local[:,1] .== 1.0 ) +
    #)
    income_sort_index = sortperm(income_policy);
    income_sorted = income_policy[income_sort_index];
    current_distr_income_sorted = current_distr_store[income_sort_index];
    cumu_income_distr =  cumsum(current_distr_income_sorted,dims = 1);
    income_sort_index_urban = sortperm(income_policy_urban);
    income_sorted_urban = income_policy_urban[income_sort_index_urban];
    current_distr_income_sorted_urban = current_distr_store_urban[income_sort_index_urban];
    cumu_income_distr_urban =  cumsum(current_distr_income_sorted_urban,dims = 1)
    income_sort_index_rural = sortperm(income_policy_rural);
    income_sorted_rural = income_policy_rural[income_sort_index_rural];
    current_distr_income_sorted_rural = current_distr_store_rural[income_sort_index_rural];
    cumu_income_distr_rural =  cumsum(current_distr_income_sorted_rural,dims = 1)
    income_mean = sum(income_sorted.*current_distr_income_sorted)
    income_mean_urban = sum(income_sorted_urban.*current_distr_income_sorted_urban)
    income_mean_rural = sum(income_sorted_rural.*current_distr_income_sorted_rural)
    #Second moments:
    #income_worker_2ndmom_mean = sum(Y_W_fine_policy[:,1].^2 .* stay_workers + Y_W_fine_policy[:,2].^2 .* exit_staple_to_work + Y_W_fine_policy[:,3].^2 .* exit_cashcrop_to_work);
    #income_staple_2ndmom_mean = sum(Y_S_fine[:,1].^2 .* entrants_staple_from_workers + Y_S_fine[:,2].^2 .* incumbents_staple + Y_S_fine[:,3].^2 .* exit_cashcrop_to_staple);
    #income_cashcrop_2ndmom_mean = sum(Y_B_fine[:,1].^2 .* entrants_cashcrop_from_workers + Y_B_fine[:,2].^2 .* entrants_from_staple_to_cashcrop + Y_B_fine[:,3].^2 .* incumbents_cashcrop);

    #Standard deviation
    #income_worker_std = (income_worker_2ndmom_mean - income_worker_mean.^2).^(1/2)
    #income_staple_std = (income_staple_2ndmom_mean - income_staple_mean.^2).^(1/2)
    #income_cashcrop_std = (income_cashcrop_2ndmom_mean - income_cashcrop_mean.^2).^(1/2)

    #### INC/CONS/WEALTH INEQUALITY MEASURES
    a_grid_size = size(agrid_fine)[1];
    z_grid_size = size(z)[1]*size(z_W)[1];
    worker_past_dist_mat = reshape(worker_past_dist,a_grid_size,z_grid_size);
    staple_past_dist_mat = reshape(staple_past_dist,a_grid_size,z_grid_size);
    cash_crop_past_dist_mat = reshape(cash_crop_past_dist,a_grid_size,z_grid_size);

    mean_wealth = sum(worker_past_dist.*s_fine[:,1] + staple_past_dist.*s_fine[:,1] +
    cash_crop_past_dist.*s_fine[:,1]);
    sd_wealth = (sum(worker_past_dist.*s_fine[:,1].^2 + staple_past_dist.*s_fine[:,1].^2 +
    cash_crop_past_dist.*s_fine[:,1].^2) - mean_wealth^2)^(1/2);

    dist0 = worker_past_dist_mat + staple_past_dist_mat + cash_crop_past_dist_mat;
    dist1 = sum(dist0,dims = 2);
    cumu_wealth = cumsum(dist1,dims = 1);
    p10_index_wealth = findfirst(cumu_wealth.>0.1);
    p50_index_wealth = findfirst(cumu_wealth.>0.5);
    p90_index_wealth = findfirst(cumu_wealth.>0.9);
    p99_index_wealth = findfirst(cumu_wealth.>0.99);
    p10_wealth_tmp = sum(dist1[p10_index_wealth[1]:end].*agrid_fine[p10_index_wealth[1]:end])/mean_wealth;
    p50_wealth_tmp = sum(dist1[p50_index_wealth[1]:end].*agrid_fine[p50_index_wealth[1]:end])/mean_wealth;
    p90_wealth_tmp = sum(dist1[p90_index_wealth[1]:end].*agrid_fine[p90_index_wealth[1]:end])/mean_wealth;
    p99_wealth_tmp = sum(dist1[p99_index_wealth[1]:end].*agrid_fine[p99_index_wealth[1]:end])/mean_wealth;

    wealth_of_workers = sum(worker_past_dist.*s_fine[:,1])/mean_wealth;
    wealth_of_staples = sum(staple_past_dist.*s_fine[:,1])/mean_wealth;
    wealth_of_cashcrop = sum(cash_crop_past_dist.*s_fine[:,1])/mean_wealth;


    worker_current_dist_mat = reshape(current_workers,a_grid_size,z_grid_size);
    staple_current_dist_mat = reshape(current_staple,a_grid_size,z_grid_size);
    cash_crop_current_dist_mat = reshape(current_cashcrop,a_grid_size,z_grid_size);
    #
    #+ staple_current_dist_mat + cash_crop_current_dist_mat
    dist0_urban = worker_current_dist_mat./current_worker_pop;
    dist1_urban = sum(dist0_urban,dims = 2);
    mean_wealth_urban = sum(agrid_fine.*dist1_urban);# - (entry_costs_to_workers + maintenance_costs_for_workers)/(current_worker_pop);
    cumu_wealth_urban = cumsum(dist1_urban,dims = 1);
    p10_index_wealth_urban = findfirst(cumu_wealth_urban.>0.1);
    p50_index_wealth_urban = findfirst(cumu_wealth_urban.>0.5);
    p90_index_wealth_urban = findfirst(cumu_wealth_urban.>0.9);
    p99_index_wealth_urban = findfirst(cumu_wealth_urban.>0.99);
    if isnothing(p10_index_wealth_urban)==0
        p10_wealth_urban = sum(dist1_urban[p10_index_wealth_urban[1]:end].*agrid_fine[p10_index_wealth_urban[1]:end])/mean_wealth_urban;
        p50_wealth_urban = sum(dist1_urban[p50_index_wealth_urban[1]:end].*agrid_fine[p50_index_wealth_urban[1]:end])/mean_wealth_urban;
        p90_wealth_urban = sum(dist1_urban[p90_index_wealth_urban[1]:end].*agrid_fine[p90_index_wealth_urban[1]:end])/mean_wealth_urban;
        p99_wealth_urban = sum(dist1_urban[p99_index_wealth_urban[1]:end].*agrid_fine[p99_index_wealth_urban[1]:end])/mean_wealth_urban;
    else
        p10_wealth_urban=0.0
        p50_wealth_urban=0.0
        p90_wealth_urban=0.0
        p99_wealth_urban=0.0
    end

    dist0_rural = (staple_current_dist_mat + cash_crop_current_dist_mat)/(current_staple_pop + current_cashcrop_pop);
    dist1_rural = sum(dist0_rural,dims = 2);
    mean_wealth_rural = sum(agrid_fine.*dist1_rural);# - (entry_costs_to_staples + entry_costs_to_cashcrops + maintenance_costs_to_staples +
                                                    #    maintenance_costs_to_cashcrops)/(current_staple_pop + current_cashcrop_pop);
    cumu_wealth_rural = cumsum(dist1_rural,dims = 1);
    p10_index_wealth_rural = findfirst(cumu_wealth_rural.>0.1);
    p50_index_wealth_rural = findfirst(cumu_wealth_rural.>0.5);
    p90_index_wealth_rural = findfirst(cumu_wealth_rural.>0.9);
    p99_index_wealth_rural = findfirst(cumu_wealth_rural.>0.99);
    if isnothing(p10_index_wealth_rural)==0
        p10_wealth_rural = sum(dist1_rural[p10_index_wealth_rural[1]:end].*agrid_fine[p10_index_wealth_rural[1]:end])/mean_wealth_rural;
        p50_wealth_rural = sum(dist1_rural[p50_index_wealth_rural[1]:end].*agrid_fine[p50_index_wealth_rural[1]:end])/mean_wealth_rural;
        p90_wealth_rural = sum(dist1_rural[p90_index_wealth_rural[1]:end].*agrid_fine[p90_index_wealth_rural[1]:end])/mean_wealth_rural;
        p99_wealth_rural = sum(dist1_rural[p99_index_wealth_rural[1]:end].*agrid_fine[p99_index_wealth_rural[1]:end])/mean_wealth_rural;
    else
        p10_wealth_urban=0.0
        p50_wealth_urban=0.0
        p90_wealth_urban=0.0
        p99_wealth_urban=0.0
    end

    p10_index_income = findfirst(cumu_income_distr.>0.1);
    p50_index_income = findfirst(cumu_income_distr.>0.5);
    p90_index_income = findfirst(cumu_income_distr.>0.9);
    p99_index_income = findfirst(cumu_income_distr.>0.99);
    p10_income_tmp = sum(current_distr_income_sorted[p10_index_income[1]:end].*income_sorted[p10_index_income[1]:end])/income_mean;
    p50_income_tmp = sum(current_distr_income_sorted[p50_index_income[1]:end].*income_sorted[p50_index_income[1]:end])/income_mean;
    p90_income_tmp = sum(current_distr_income_sorted[p90_index_income[1]:end].*income_sorted[p90_index_income[1]:end])/income_mean;
    p99_income_tmp = sum(current_distr_income_sorted[p99_index_income[1]:end].*income_sorted[p99_index_income[1]:end])/income_mean;

    p10_index_income_urban = findfirst(cumu_income_distr_urban.>0.1);
    p50_index_income_urban = findfirst(cumu_income_distr_urban.>0.5);
    p90_index_income_urban = findfirst(cumu_income_distr_urban.>0.9);
    p99_index_income_urban = findfirst(cumu_income_distr_urban.>0.99);
    p10_income_tmp_urban = sum(current_distr_income_sorted_urban[p10_index_income_urban[1]:end].*income_sorted_urban[p10_index_income_urban[1]:end])/income_mean_urban;
    p50_income_tmp_urban = sum(current_distr_income_sorted_urban[p50_index_income_urban[1]:end].*income_sorted_urban[p50_index_income_urban[1]:end])/income_mean_urban;
    p90_income_tmp_urban = sum(current_distr_income_sorted_urban[p90_index_income_urban[1]:end].*income_sorted_urban[p90_index_income_urban[1]:end])/income_mean_urban;
    p99_income_tmp_urban = sum(current_distr_income_sorted_urban[p99_index_income_urban[1]:end].*income_sorted_urban[p99_index_income_urban[1]:end])/income_mean_urban;

    p10_index_income_rural = findfirst(cumu_income_distr_rural.>0.1);
    p50_index_income_rural = findfirst(cumu_income_distr_rural.>0.5);
    p90_index_income_rural = findfirst(cumu_income_distr_rural.>0.9);
    p99_index_income_rural = findfirst(cumu_income_distr_rural.>0.99);
    p10_income_tmp_rural = sum(current_distr_income_sorted_rural[p10_index_income_rural[1]:end].*income_sorted_rural[p10_index_income_rural[1]:end])/income_mean_rural;
    p50_income_tmp_rural = sum(current_distr_income_sorted_rural[p50_index_income_rural[1]:end].*income_sorted_rural[p50_index_income_rural[1]:end])/income_mean_rural;
    p90_income_tmp_rural = sum(current_distr_income_sorted_rural[p90_index_income_rural[1]:end].*income_sorted_rural[p90_index_income_rural[1]:end])/income_mean_rural;
    p99_income_tmp_rural = sum(current_distr_income_sorted_rural[p99_index_income_rural[1]:end].*income_sorted_rural[p99_index_income_rural[1]:end])/income_mean_rural;


    p10_index_cons = findfirst(cumu_consumption_distr.>0.1);
    p50_index_cons = findfirst(cumu_consumption_distr.>0.5);
    p90_index_cons = findfirst(cumu_consumption_distr.>0.9);
    p99_index_cons = findfirst(cumu_consumption_distr.>0.99);
    p10_cons_tmp = sum(current_distr_consumption_sorted[p10_index_cons[1]:end].*consumption_sorted[p10_index_cons[1]:end])/mean_cons;
    p50_cons_tmp = sum(current_distr_consumption_sorted[p50_index_cons[1]:end].*consumption_sorted[p50_index_cons[1]:end])/mean_cons;
    p90_cons_tmp = sum(current_distr_consumption_sorted[p90_index_cons[1]:end].*consumption_sorted[p90_index_cons[1]:end])/mean_cons;
    p99_cons_tmp = sum(current_distr_consumption_sorted[p99_index_cons[1]:end].*consumption_sorted[p99_index_cons[1]:end])/mean_cons;
    #Urban cons
    p10_index_cons_urban = findfirst(cumu_consumption_distr_urban.>0.1);
    p50_index_cons_urban = findfirst(cumu_consumption_distr_urban.>0.5);
    p90_index_cons_urban = findfirst(cumu_consumption_distr_urban.>0.9);
    p99_index_cons_urban = findfirst(cumu_consumption_distr_urban.>0.99);
    p10_cons_tmp_urban = sum(current_distr_consumption_sorted_urban[p10_index_cons_urban[1]:end].*consumption_sorted_urban[p10_index_cons_urban[1]:end])/mean_cons_urban;
    p50_cons_tmp_urban = sum(current_distr_consumption_sorted_urban[p50_index_cons_urban[1]:end].*consumption_sorted_urban[p50_index_cons_urban[1]:end])/mean_cons_urban;
    p90_cons_tmp_urban = sum(current_distr_consumption_sorted_urban[p90_index_cons_urban[1]:end].*consumption_sorted_urban[p90_index_cons_urban[1]:end])/mean_cons_urban;
    p99_cons_tmp_urban = sum(current_distr_consumption_sorted_urban[p99_index_cons_urban[1]:end].*consumption_sorted_urban[p99_index_cons_urban[1]:end])/mean_cons_urban;
    #Rural cons
    p10_index_cons_rural = findfirst(cumu_consumption_distr_rural.>0.1);
    p50_index_cons_rural = findfirst(cumu_consumption_distr_rural.>0.5);
    p90_index_cons_rural = findfirst(cumu_consumption_distr_rural.>0.9);
    p99_index_cons_rural = findfirst(cumu_consumption_distr_rural.>0.99);
    p10_cons_tmp_rural = sum(current_distr_consumption_sorted_rural[p10_index_cons_rural[1]:end].*consumption_sorted_rural[p10_index_cons_rural[1]:end])/mean_cons_rural;
    p50_cons_tmp_rural = sum(current_distr_consumption_sorted_rural[p50_index_cons_rural[1]:end].*consumption_sorted_rural[p50_index_cons_rural[1]:end])/mean_cons_rural;
    p90_cons_tmp_rural = sum(current_distr_consumption_sorted_rural[p90_index_cons_rural[1]:end].*consumption_sorted_rural[p90_index_cons_rural[1]:end])/mean_cons_rural;
    p99_cons_tmp_rural = sum(current_distr_consumption_sorted_rural[p99_index_cons_rural[1]:end].*consumption_sorted_rural[p99_index_cons_rural[1]:end])/mean_cons_rural;

    # Now to the actual RCT evidence    p10_index_cons identifies the group we care about
    # Changed on 9/11/2022 to those who choose to live in rural
    # This questionable.
    RCT_factor = 0.25
    cons_rural_p10_level = consumption_sorted_rural[p10_index_cons_rural]
    s_fine_RCT_workers = copy(s_fine);
    s_fine_RCT_staples = copy(s_fine);
    s_fine_RCT_cashcrop = copy(s_fine);
    treated_past_worker = (consumption_policy_rural[(1):(ns_fine)].>0).* (consumption_policy_rural[(1):(ns_fine)].<=cons_rural_p10_level);
    treated_past_staples = (consumption_policy_rural[(ns_fine+1):(2*ns_fine)].>0).* (consumption_policy_rural[(ns_fine+1):(2*ns_fine)].<=cons_rural_p10_level);
    treated_past_cashcrop =(consumption_policy_rural[(2*ns_fine+1):(3*ns_fine)].>0).* (consumption_policy_rural[(2*ns_fine+1):(3*ns_fine)].<=cons_rural_p10_level);
    avg_cons_factor=RCT_factor *sum(worker_past_dist.*treated_past_worker.*consumption_policy_rural[1:ns_fine] +
                                    staple_past_dist.*treated_past_staples.*consumption_policy_rural[(ns_fine+1):(2*ns_fine)] +
                                    cash_crop_past_dist.*treated_past_cashcrop.*consumption_policy_rural[(2*ns_fine+1):(3*ns_fine)]) / sum(worker_past_dist.*treated_past_worker +
                                                                        staple_past_dist.*treated_past_staples + cash_crop_past_dist.*treated_past_cashcrop)
    #No need for the 1/P_W, because its already nominal consumption that should be added to wealth
    # s_fine_RCT_workers[:,1] = s_fine_RCT_workers[:,1] + P_W*RCT_factor * treated_past_worker.*consumption_policy[1:ns_fine];
    # s_fine_RCT_staples[:,1] = s_fine_RCT_staples[:,1] + P_W*RCT_factor * treated_past_staples.*consumption_policy[(ns_fine+1):(2*ns_fine)];
    # s_fine_RCT_cashcrop[:,1] = s_fine_RCT_cashcrop[:,1] + P_W*RCT_factor * treated_past_cashcrop.*consumption_policy[(2*ns_fine+1):(3*ns_fine)];
    s_fine_RCT_workers[:,1] = s_fine_RCT_workers[:,1] .+ treated_past_worker.*avg_cons_factor;
    s_fine_RCT_staples[:,1] = s_fine_RCT_staples[:,1] .+ treated_past_staples.*avg_cons_factor;
    s_fine_RCT_cashcrop[:,1] = s_fine_RCT_cashcrop[:,1] .+ treated_past_cashcrop.*avg_cons_factor;
    # Simplification is needed unfortunately - assume the transfer cannot be conditional on your previous occupation choice, therefore you get the maximum
    # consumption, irrespective of your occupation (not a far fetched theory). As we dont care about those who are workers,
    #the main difference is the cashcrop vs. staples are treated equal
    s_fine_RCT = copy(s_fine);
    s_fine_RCT_mat = zeros(ns_fine,3);
    s_fine_RCT_mat[:,1] =  s_fine_RCT_workers[:,1];
    s_fine_RCT_mat[:,2] =  s_fine_RCT_staples[:,1];
    s_fine_RCT_mat[:,3] =  s_fine_RCT_cashcrop[:,1];
    s_fine_RCT[:,1] = maximum(s_fine_RCT_mat,dims=2);
    # Change optimal consumption too, by assuming that the consumption to asset ratio stays the same
    #cons_fine_RCT = cons_fine_local.*min.((s_fine_RCT[:,1] ./s_fine[:,1]),1.0+RCT_factor);
    #cons_fine_RCT = copy(cons_fine_local);
    # Some of this is clearly too large - maximize to the 1.25*consumption

    # Evaluate the new approximation and objects with the new state vector:

        (θ_RCT,labor_prod_RCT,tol,P_W_RCT,Y_W_RCT,coeff_λ_2_cashcrop_residual_unconstrained_RCT,coeff_λ_2_cashcrop_residual_constrained_RCT,
        x_B_c1_RCT,π_B_only_B_c1_RCT,λ_B_only_B_c1_RCT,P_B_c1_RCT,Y_B_c1_RCT,
        coeff_λ_2_s_RCT,P_S_c1_RCT,P_S_c2_RCT,Y_S_c1_RCT,Y_S_c2_RCT,x_S_c1_RCT, x_S_c2_RCT,labor_allocated_interior_c3a_RCT,
        λ_B_interior_c3a_RCT,x_SC_interior_c3a_RCT,x_BC_interior_c3a_RCT,Y_B_c3a_RCT,P_B_c3a_RCT,P_B_c3b_RCT,q_S_c1_RCT,q_S_c2_RCT,q_B_c1_RCT,q_S_c3a_RCT,
        q_B_c3a_RCT,q_S_c3b_RCT,q_B_c3b_RCT,x_SC_interior_c3b_RCT,x_BC_interior_c3b_RCT,labor_allocated_interior_c3b_RCT,Y_B_c3b_RCT, c_S_mat_RCT,c_B_mat_RCT,
        c_M_mat_RCT,x_S_mat_RCT,x_B_mat_RCT,q_S_mat_RCT,q_B_mat_RCT,land_B_mat_RCT, λ_2_mat_RCT,P_B_mat_RCT,Y_B_mat_RCT,feasibility_mat_RCT,C_max_mat_RCT,
        C_min_mat_RCT,q_S_staples_RCT,c_S_staples_RCT,c_B_staples_RCT,c_M_staples_RCT,P_S_staples_RCT,x_S_staples_RCT,λ_2_S_staples_RCT,unfeasible_mat_RCT,
        Y_S_potential_RCT,TC_mat_RCT,C_max_staple_RCT,C_min_staple_RCT,C_max_staple_constrained_RCT,
        C_min_staple_constrained_RCT,TC_S_c3_constrained_RCT,x_S_c3_constrained_RCT,q_S_c3_constrained_RCT,c_S_c3_constrained_RCT,
        x_S_mat_3c_RCT,x_B_mat_3c_RCT,land_B_mat_3c_RCT,λ_2_mat_3c_RCT,TC_mat_3c_RCT) = income_creator_no_approx(s_fine_RCT,ns_fine,
        z,z_W,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n_fine,κ,Q_S,ϵ,ψ_S,ψ_B,ψ_M,coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,
        C_max_unconstrained ,C_max_constrained,C_min_unconstrained,C_min_constrained, coeff_λ_2_s,C_grid_fine,fspace_C_fine,C_max_staple,C_min_staple,C_max_staple_constrained,
    C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained);


    min_C_applied_RCT,max_C_applied_RCT =    bounds_consumption(P_W_RCT,Y_W_RCT,s_fine_RCT,r,ρ,w,
coeff_λ_2_cashcrop_residual_unconstrained_RCT,coeff_λ_2_cashcrop_residual_constrained_RCT,θ_RCT,
fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,
a_min,a_max,x_B_c1_RCT,π_B_only_B_c1_RCT,λ_B_only_B_c1_RCT,P_B_c1_RCT,Y_B_c1_RCT,
coeff_λ_2_s_RCT,P_S_c1_RCT,P_S_c2_RCT,Y_S_c1_RCT,Y_S_c2_RCT,x_S_c1_RCT, x_S_c2_RCT,
labor_allocated_interior_c3a_RCT,λ_B_interior_c3a_RCT,x_SC_interior_c3a_RCT,
x_BC_interior_c3a_RCT,Y_B_c3a_RCT,P_B_c3a_RCT,P_B_c3b_RCT,q_S_c1_RCT,q_S_c2_RCT,
q_B_c1_RCT,q_S_c3a_RCT,q_B_c3a_RCT,q_S_c3b_RCT,q_B_c3b_RCT,
x_SC_interior_c3b_RCT,x_BC_interior_c3b_RCT,labor_allocated_interior_c3b_RCT,Y_B_c3b_RCT,
c_S_mat_RCT,c_B_mat_RCT,c_M_mat_RCT,x_S_mat_RCT,x_B_mat_RCT,q_S_mat_RCT,q_B_mat_RCT,land_B_mat_RCT,
λ_2_mat_RCT,P_B_mat_RCT,Y_B_mat_RCT,feasibility_mat_RCT,C_max_mat_RCT,C_min_mat_RCT,
q_S_staples_RCT,c_S_staples_RCT,c_B_staples_RCT,c_M_staples_RCT,P_S_staples_RCT,
x_S_staples_RCT,λ_2_S_staples_RCT,unfeasible_mat_RCT,Y_S_potential_RCT,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat_RCT,C_max_staple_RCT,
C_min_staple_RCT,C_max_staple_constrained_RCT,C_min_staple_constrained_RCT,TC_S_c3_constrained_RCT,x_S_c3_constrained_RCT,q_S_c3_constrained_RCT,c_S_c3_constrained_RCT,
x_S_mat_3c_RCT,x_B_mat_3c_RCT,land_B_mat_3c_RCT,λ_2_mat_3c_RCT,TC_mat_3c_RCT);

    (cons_fine_RCT,a_prime_fine_RCT,future_occupation_fine_RCT) = RCT_opt_cons(coeff,ns_fine,P_kron_fine,min_C_applied_RCT,max_C_applied_RCT,Phi_z_fine,
            β,fspace_a,σ,P_W_RCT,Y_W_RCT,s_fine_RCT,r,ρ,w,coeff_λ_2_cashcrop_residual_unconstrained_RCT,coeff_λ_2_cashcrop_residual_constrained_RCT,θ_RCT,
            fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,κ,tol,a_min,x_B_c1_RCT,π_B_only_B_c1_RCT,λ_B_only_B_c1_RCT,P_B_c1_RCT,Y_B_c1_RCT,
            coeff_λ_2_s_RCT,P_S_c1_RCT,P_S_c2_RCT,Y_S_c1_RCT,Y_S_c2_RCT,x_S_c1_RCT, x_S_c2_RCT,labor_allocated_interior_c3a_RCT,λ_B_interior_c3a_RCT,x_SC_interior_c3a_RCT,
            x_BC_interior_c3a_RCT,Y_B_c3a_RCT,P_B_c3a_RCT,P_B_c3b_RCT,q_S_c1_RCT,q_S_c2_RCT,q_B_c1_RCT,q_S_c3a_RCT,q_B_c3a_RCT,q_S_c3b_RCT,q_B_c3b_RCT,
            x_SC_interior_c3b_RCT,x_BC_interior_c3b_RCT,labor_allocated_interior_c3b_RCT,Y_B_c3b_RCT,c_S_mat_RCT,c_B_mat_RCT,c_M_mat_RCT,x_S_mat_RCT,x_B_mat_RCT,q_S_mat_RCT,
            q_B_mat_RCT,land_B_mat_RCT,λ_2_mat_RCT,P_B_mat_RCT,Y_B_mat_RCT,feasibility_mat_RCT,C_max_mat_RCT,C_min_mat_RCT,q_S_staples_RCT,c_S_staples_RCT,c_B_staples_RCT,
            c_M_staples_RCT,P_S_staples_RCT,x_S_staples_RCT,λ_2_S_staples_RCT,unfeasible_mat_RCT,Y_S_potential_RCT,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat_RCT,a_max,P_kron1,Q_trans,
            C_max_staple_RCT,C_min_staple_RCT,C_max_staple_constrained_RCT,
            C_min_staple_constrained_RCT,TC_S_c3_constrained_RCT,x_S_c3_constrained_RCT,q_S_c3_constrained_RCT,c_S_c3_constrained_RCT,
            x_S_mat_3c_RCT,x_B_mat_3c_RCT,land_B_mat_3c_RCT,λ_2_mat_3c_RCT,TC_mat_3c_RCT);

    # Policy functions of all people with the new states:
    Y_W_fine_RCT = zeros(ns_fine,3);
    Y_S_fine_RCT = zeros(ns_fine,3);
    c_S_S_fine_RCT= zeros(ns_fine,3);
    c_B_S_fine_RCT= zeros(ns_fine,3);
    c_M_S_fine_RCT= zeros(ns_fine,3);
    q_S_S_fine_RCT= zeros(ns_fine,3);
    P_S_fine_RCT= zeros(ns_fine,3);
    x_S_S_fine_RCT= zeros(ns_fine,3);
    solve_staple_index_S_fine_RCT= zeros(Int64,ns_fine,3);
    λ_2_S_fine_RCT= zeros(ns_fine,3);
    future_asset_S_RCT= zeros(ns_fine,3);
    future_asset_C_RCT= zeros(ns_fine,3);
    c_S_B_fine_RCT= zeros(ns_fine,3);
    c_B_B_fine_RCT= zeros(ns_fine,3);
    c_M_B_fine_RCT= zeros(ns_fine,3);
    x_SC_fine_RCT= zeros(ns_fine,3);
    x_BC_fine_RCT= zeros(ns_fine,3);
    land_C_fine_RCT= zeros(ns_fine,3);
    λ_2_fine_RCT= zeros(ns_fine,3);
    P_B_fine_RCT= zeros(ns_fine,3);
    Y_B_fine_RCT= zeros(ns_fine,3);
    q_S_B_fine_RCT= zeros(ns_fine,3);
    q_B_B_fine_RCT= zeros(ns_fine,3);
    solve_cash_crop_index_B_fine_RCT= zeros(Int64,ns_fine,3);
    solve_staple_index_B_fine_RCT= zeros(Int64,ns_fine,3);
    TC_fine_RCT= zeros(ns_fine,3);
    for j = 1:3
        for jj =1:3
            if jj == 1
                Y_W_fine_RCT[:,j] =  P_W *cons_fine_RCT[:,j]  +Y_W_RCT  .+ w*FM_W .+ w*F_W * (j != 1);
            end
            if jj == 2
            (future_asset_S_RCT[:,j],Y_S_fine_RCT[:,j],c_S_S_fine_RCT[:,j],c_B_S_fine_RCT[:,j],c_M_S_fine_RCT[:,j],q_S_S_fine_RCT[:,j],
            P_S_fine_RCT[:,j],x_S_S_fine_RCT[:,j],solve_staple_index_S_fine_RCT[:,j],λ_2_S_fine_RCT[:,j]) = policy_function_creator(
                cons_fine_RCT[:,j],jj,j,P_W_RCT,Y_W_RCT,s_fine_RCT,r,ρ,w,coeff_λ_2_cashcrop_residual_unconstrained_RCT,
                coeff_λ_2_cashcrop_residual_constrained_RCT,θ_RCT,fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,
                p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,a_min,x_B_c1_RCT,π_B_only_B_c1_RCT,λ_B_only_B_c1_RCT,P_B_c1_RCT,Y_B_c1_RCT,
                coeff_λ_2_s_RCT,P_S_c1_RCT,P_S_c2_RCT,Y_S_c1_RCT,Y_S_c2_RCT,x_S_c1_RCT, x_S_c2_RCT,
                labor_allocated_interior_c3a_RCT,λ_B_interior_c3a_RCT,x_SC_interior_c3a_RCT,x_BC_interior_c3a_RCT,Y_B_c3a_RCT,P_B_c3a_RCT,P_B_c3b_RCT,q_S_c1_RCT,
                q_S_c2_RCT,q_B_c1_RCT,q_S_c3a_RCT,q_B_c3a_RCT,q_S_c3b_RCT,q_B_c3b_RCT,x_SC_interior_c3b_RCT,x_BC_interior_c3b_RCT,labor_allocated_interior_c3b_RCT,
                Y_B_c3b_RCT,c_S_mat_RCT,c_B_mat_RCT,c_M_mat_RCT,x_S_mat_RCT,x_B_mat_RCT,q_S_mat_RCT,q_B_mat_RCT,land_B_mat_RCT,λ_2_mat_RCT,P_B_mat_RCT,Y_B_mat_RCT,
                feasibility_mat_RCT,C_max_mat_RCT,C_min_mat_RCT,q_S_staples_RCT,c_S_staples_RCT,c_B_staples_RCT,c_M_staples_RCT,P_S_staples_RCT,x_S_staples_RCT,
                λ_2_S_staples_RCT,unfeasible_mat_RCT,Y_S_potential_RCT,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat_RCT,C_max_staple_RCT,C_min_staple_RCT,C_max_staple_constrained_RCT,
                C_min_staple_constrained_RCT,TC_S_c3_constrained_RCT,x_S_c3_constrained_RCT,q_S_c3_constrained_RCT,c_S_c3_constrained_RCT,
                x_S_mat_3c_RCT,x_B_mat_3c_RCT,land_B_mat_3c_RCT,λ_2_mat_3c_RCT,TC_mat_3c_RCT);
            elseif jj == 3
                    (future_asset_C_RCT[:,j],c_S_B_fine_RCT[:,j],c_B_B_fine_RCT[:,j],c_M_B_fine_RCT[:,j],x_SC_fine_RCT[:,j],x_BC_fine_RCT[:,j],
                    land_C_fine_RCT[:,j],λ_2_fine_RCT[:,j],P_B_fine_RCT[:,j],Y_B_fine_RCT[:,j],q_S_B_fine_RCT[:,j],q_B_B_fine_RCT[:,j],solve_cash_crop_index_B_fine_RCT[:,j]
                    ,solve_staple_index_B_fine_RCT[:,j],TC_fine_RCT[:,j]) = policy_function_creator(
                        cons_fine_RCT[:,j],jj,j,P_W_RCT,Y_W_RCT,s_fine_RCT,r,ρ,w,coeff_λ_2_cashcrop_residual_unconstrained_RCT,
                        coeff_λ_2_cashcrop_residual_constrained_RCT,θ_RCT,fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,
                        p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,a_min,x_B_c1_RCT,π_B_only_B_c1_RCT,λ_B_only_B_c1_RCT,P_B_c1_RCT,Y_B_c1_RCT,
                        coeff_λ_2_s_RCT,P_S_c1_RCT,P_S_c2_RCT,Y_S_c1_RCT,Y_S_c2_RCT,x_S_c1_RCT, x_S_c2_RCT,
                        labor_allocated_interior_c3a_RCT,λ_B_interior_c3a_RCT,x_SC_interior_c3a_RCT,x_BC_interior_c3a_RCT,Y_B_c3a_RCT,P_B_c3a_RCT,P_B_c3b_RCT,q_S_c1_RCT,
                        q_S_c2_RCT,q_B_c1_RCT,q_S_c3a_RCT,q_B_c3a_RCT,q_S_c3b_RCT,q_B_c3b_RCT,x_SC_interior_c3b_RCT,x_BC_interior_c3b_RCT,labor_allocated_interior_c3b_RCT,
                        Y_B_c3b_RCT,c_S_mat_RCT,c_B_mat_RCT,c_M_mat_RCT,x_S_mat_RCT,x_B_mat_RCT,q_S_mat_RCT,q_B_mat_RCT,land_B_mat_RCT,λ_2_mat_RCT,P_B_mat_RCT,Y_B_mat_RCT,
                        feasibility_mat_RCT,C_max_mat_RCT,C_min_mat_RCT,q_S_staples_RCT,c_S_staples_RCT,c_B_staples_RCT,c_M_staples_RCT,P_S_staples_RCT,x_S_staples_RCT,
                        λ_2_S_staples_RCT,unfeasible_mat_RCT,Y_S_potential_RCT,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat_RCT,C_max_staple_RCT,C_min_staple_RCT,C_max_staple_constrained_RCT,
                        C_min_staple_constrained_RCT,TC_S_c3_constrained_RCT,x_S_c3_constrained_RCT,q_S_c3_constrained_RCT,c_S_c3_constrained_RCT,
                        x_S_mat_3c_RCT,x_B_mat_3c_RCT,land_B_mat_3c_RCT,λ_2_mat_3c_RCT,TC_mat_3c_RCT)
                end
            end
    end
    # Now lets get some aggregates. Since no dynamic optimization, we do not have to care about nonrecipients reacting:
    entrants_staple_from_workers_RCT = (treated_past_worker.==1.0) .* entrants_staple_from_workers;
    incumbents_staple_RCT = (treated_past_staples .==1.0) .* incumbents_staple;
    exit_cashcrop_to_staple_RCT = (treated_past_cashcrop.==1.0) .* exit_cashcrop_to_staple;
    entrants_cashcrop_from_workers_RCT =  (treated_past_worker.==1.0) .*entrants_cashcrop_from_workers;
    entrants_from_staple_to_cashcrop_RCT = (treated_past_staples.==1.0) .* entrants_from_staple_to_cashcrop;
    incumbents_cashcrop_RCT = (treated_past_cashcrop.==1.0) .* incumbents_cashcrop;
    q_S_staple_sum_RCT = sum(q_S_S_fine_RCT[:,1] .*entrants_staple_from_workers_RCT + q_S_S_fine_RCT[:,2] .*incumbents_staple_RCT
    + q_S_S_fine_RCT[:,3] .*exit_cashcrop_to_staple_RCT);
    x_S_staple_sum_RCT = sum(x_S_S_fine_RCT[:,1] .*entrants_staple_from_workers_RCT + x_S_S_fine_RCT[:,2] .*incumbents_staple_RCT
    + x_S_S_fine_RCT[:,3] .*exit_cashcrop_to_staple_RCT);

    q_S_cashcrop_sum_RCT = sum(q_S_B_fine_RCT[:,1] .*entrants_cashcrop_from_workers_RCT + q_S_B_fine_RCT[:,3] .*incumbents_cashcrop_RCT
    + q_S_B_fine_RCT[:,2] .*entrants_from_staple_to_cashcrop_RCT);
    q_B_cashcrop_sum_RCT = sum(q_B_B_fine_RCT[:,1] .*entrants_cashcrop_from_workers_RCT + q_B_B_fine_RCT[:,3] .*incumbents_cashcrop_RCT
    + q_B_B_fine_RCT[:,2] .*entrants_from_staple_to_cashcrop_RCT);
    x_S_cashcrop_sum_RCT = sum(x_SC_fine_RCT[:,1] .*entrants_cashcrop_from_workers_RCT + x_SC_fine_RCT[:,3] .*incumbents_cashcrop_RCT
    + x_SC_fine_RCT[:,2] .*entrants_from_staple_to_cashcrop_RCT);
    x_B_cashcrop_sum_RCT = sum(x_BC_fine_RCT[:,1] .*entrants_cashcrop_from_workers_RCT + x_BC_fine_RCT[:,3] .*incumbents_cashcrop_RCT
    + x_BC_fine_RCT[:,2] .*entrants_from_staple_to_cashcrop_RCT);

    # Lets do only treated
    q_S_staple_sum_only_treated_RCT = sum(q_S_S_fine_RCT[:,1] .*entrants_staple_from_workers_RCT.* treated_past_worker
        + q_S_S_fine_RCT[:,2] .*incumbents_staple_RCT .*treated_past_staples
    + q_S_S_fine_RCT[:,3] .*exit_cashcrop_to_staple_RCT .*treated_past_cashcrop);
    x_S_staple_sum_only_treated_RCT = sum(x_S_S_fine_RCT[:,1] .*entrants_staple_from_workers_RCT.* treated_past_worker
        + x_S_S_fine_RCT[:,2] .*incumbents_staple_RCT .*treated_past_staples
    + x_S_S_fine_RCT[:,3] .*exit_cashcrop_to_staple_RCT .*treated_past_cashcrop);


    q_S_cashcrop_sum_only_treated_RCT = sum(q_S_B_fine_RCT[:,1] .*entrants_cashcrop_from_workers_RCT.* treated_past_worker
    + q_S_B_fine_RCT[:,3] .*incumbents_cashcrop_RCT.*treated_past_cashcrop
    + q_S_B_fine_RCT[:,2] .*entrants_from_staple_to_cashcrop_RCT.*treated_past_staples);
    q_B_cashcrop_sum_only_treated_RCT = sum(q_B_B_fine_RCT[:,1] .*entrants_cashcrop_from_workers_RCT.* treated_past_worker
    + q_B_B_fine_RCT[:,3] .*incumbents_cashcrop_RCT.*treated_past_cashcrop
    + q_B_B_fine_RCT[:,2] .*entrants_from_staple_to_cashcrop_RCT.*treated_past_staples);
    x_S_cashcrop_sum_only_treated_RCT = sum(x_SC_fine_RCT[:,1] .*entrants_cashcrop_from_workers_RCT .* treated_past_worker
    + x_SC_fine_RCT[:,3] .*incumbents_cashcrop_RCT.*treated_past_cashcrop
    + x_SC_fine_RCT[:,2] .*entrants_from_staple_to_cashcrop_RCT.*treated_past_staples);
    x_B_cashcrop_sum_only_treated_RCT = sum(x_BC_fine_RCT[:,1] .*entrants_cashcrop_from_workers_RCT.* treated_past_worker
    + x_BC_fine_RCT[:,3] .*incumbents_cashcrop_RCT.*treated_past_cashcrop
    + x_BC_fine_RCT[:,2] .*entrants_from_staple_to_cashcrop_RCT.*treated_past_staples);

    # Lets get treated group without intervention
    q_S_staple_sum_only_treated = sum(q_S_S_fine[:,1] .*entrants_staple_from_workers_RCT.* treated_past_worker
        + q_S_S_fine[:,2] .*incumbents_staple_RCT .*treated_past_staples
    + q_S_S_fine[:,3] .*exit_cashcrop_to_staple_RCT .*treated_past_cashcrop);
    x_S_staple_sum_only_treated = sum(x_S_S_fine[:,1] .*entrants_staple_from_workers_RCT.* treated_past_worker
        + x_S_S_fine[:,2] .*incumbents_staple_RCT .*treated_past_staples
    + x_S_S_fine[:,3] .*exit_cashcrop_to_staple_RCT .*treated_past_cashcrop);


    q_S_cashcrop_sum_only_treated = sum(q_S_B_fine[:,1] .*entrants_cashcrop_from_workers_RCT.* treated_past_worker
    + q_S_B_fine[:,3] .*incumbents_cashcrop_RCT.*treated_past_cashcrop
    + q_S_B_fine[:,2] .*entrants_from_staple_to_cashcrop_RCT.*treated_past_staples);
    q_B_cashcrop_sum_only_treated = sum(q_B_B_fine[:,1] .*entrants_cashcrop_from_workers_RCT.* treated_past_worker
    + q_B_B_fine[:,3] .*incumbents_cashcrop_RCT.*treated_past_cashcrop
    + q_B_B_fine[:,2] .*entrants_from_staple_to_cashcrop_RCT.*treated_past_staples);
    x_S_cashcrop_sum_only_treated = sum(x_SC_fine[:,1] .*entrants_cashcrop_from_workers_RCT .* treated_past_worker
    + x_SC_fine[:,3] .*incumbents_cashcrop_RCT.*treated_past_cashcrop
    + x_SC_fine[:,2] .*entrants_from_staple_to_cashcrop_RCT.*treated_past_staples);
    x_B_cashcrop_sum_only_treated = sum(x_BC_fine[:,1] .*entrants_cashcrop_from_workers_RCT.* treated_past_worker
    + x_BC_fine[:,3] .*incumbents_cashcrop_RCT.*treated_past_cashcrop
    + x_BC_fine[:,2] .*entrants_from_staple_to_cashcrop_RCT.*treated_past_staples);

    # I am not 100% about our treated group - we should somehow take into account their occupation choice?
    q_S_staple_sum_only_treated_RCT/q_S_staple_sum_only_treated
    x_S_staple_sum_only_treated_RCT/x_S_staple_sum_only_treated
    q_S_cashcrop_sum_only_treated_RCT/q_S_cashcrop_sum_only_treated
    q_B_cashcrop_sum_only_treated_RCT/q_B_cashcrop_sum_only_treated
    x_S_cashcrop_sum_only_treated_RCT/x_S_cashcrop_sum_only_treated
    x_B_cashcrop_sum_only_treated_RCT/x_B_cashcrop_sum_only_treated

    # Now get standard deviation of cross sectional standard deviation:
    # Second moments
    # Lets get treated group without intervention
    q_S_staple_sum_second_moment = sum(q_S_S_fine[:,1].^2 .*entrants_staple_from_workers_RCT
        + q_S_S_fine[:,2].^2 .*incumbents_staple_RCT
    + q_S_S_fine[:,3].^2 .*exit_cashcrop_to_staple_RCT);
    x_S_staple_sum_second_moment = sum(x_S_S_fine[:,1].^2 .*entrants_staple_from_workers_RCT
        + x_S_S_fine[:,2].^2 .*incumbents_staple_RCT
    + x_S_S_fine[:,3].^2 .*exit_cashcrop_to_staple_RCT);


    q_S_cashcrop_sum_second_moment = sum(q_S_B_fine[:,1].^2 .*entrants_cashcrop_from_workers_RCT
    + q_S_B_fine[:,3].^2 .*incumbents_cashcrop_RCT
    + q_S_B_fine[:,2].^2 .*entrants_from_staple_to_cashcrop_RCT);
    q_B_cashcrop_sum_second_moment = sum(q_B_B_fine[:,1].^2 .*entrants_cashcrop_from_workers_RCT
    + q_B_B_fine[:,3].^2 .*incumbents_cashcrop_RCT
    + q_B_B_fine[:,2].^2 .*entrants_from_staple_to_cashcrop_RCT);
    x_S_cashcrop_sum_second_moment = sum(x_SC_fine[:,1].^2 .*entrants_cashcrop_from_workers_RCT
    + x_SC_fine[:,3].^2 .*incumbents_cashcrop_RCT
    + x_SC_fine[:,2].^2 .*entrants_from_staple_to_cashcrop_RCT);
    x_B_cashcrop_sum_second_moment = sum(x_BC_fine[:,1].^2 .*entrants_cashcrop_from_workers_RCT
    + x_BC_fine[:,3].^2 .*incumbents_cashcrop_RCT
    + x_BC_fine[:,2].^2 .*entrants_from_staple_to_cashcrop_RCT);

    #First moment
    q_S_staple_sum_first_moment = sum(q_S_S_fine[:,1] .*entrants_staple_from_workers_RCT
        + q_S_S_fine[:,2] .*incumbents_staple_RCT
    + q_S_S_fine[:,3] .*exit_cashcrop_to_staple_RCT);
    x_S_staple_sum_first_moment = sum(x_S_S_fine[:,1] .*entrants_staple_from_workers_RCT
        + x_S_S_fine[:,2] .*incumbents_staple_RCT
    + x_S_S_fine[:,3] .*exit_cashcrop_to_staple_RCT);


    q_S_cashcrop_sum_first_moment = sum(q_S_B_fine[:,1].*entrants_cashcrop_from_workers_RCT
    + q_S_B_fine[:,3] .*incumbents_cashcrop_RCT
    + q_S_B_fine[:,2] .*entrants_from_staple_to_cashcrop_RCT);
    q_B_cashcrop_sum_first_moment = sum(q_B_B_fine[:,1].*entrants_cashcrop_from_workers_RCT
    + q_B_B_fine[:,3] .*incumbents_cashcrop_RCT
    + q_B_B_fine[:,2] .*entrants_from_staple_to_cashcrop_RCT);
    x_S_cashcrop_sum_first_moment = sum(x_SC_fine[:,1] .*entrants_cashcrop_from_workers_RCT
    + x_SC_fine[:,3] .*incumbents_cashcrop_RCT
    + x_SC_fine[:,2] .*entrants_from_staple_to_cashcrop_RCT);
    x_B_cashcrop_sum_first_moment = sum(x_BC_fine[:,1] .*entrants_cashcrop_from_workers_RCT
    + x_BC_fine[:,3] .*incumbents_cashcrop_RCT
    + x_BC_fine[:,2] .*entrants_from_staple_to_cashcrop_RCT);


    q_S_staple_sum_std = (q_S_staple_sum_second_moment - q_S_staple_sum_first_moment.^2)^(1/2)
    x_S_staple_sum_std = (x_S_staple_sum_second_moment - x_S_staple_sum_first_moment.^2)^(1/2)
    q_S_cashcrop_sum_std = (q_S_cashcrop_sum_second_moment - q_S_cashcrop_sum_first_moment.^2)^(1/2)
    q_B_cashcrop_sum_std = (q_B_cashcrop_sum_second_moment - q_B_cashcrop_sum_first_moment.^2)^(1/2)
    x_S_cashcrop_sum_std = (x_S_cashcrop_sum_second_moment - x_S_cashcrop_sum_first_moment.^2)^(1/2)
    x_B_cashcrop_sum_std = (x_B_cashcrop_sum_second_moment - x_B_cashcrop_sum_first_moment.^2)^(1/2)

    value_cashcrop_first_moment = q_S_staple_sum_only_treated + q_S_cashcrop_sum_only_treated + p_B * q_B_cashcrop_sum_only_treated;
    value_cashcrop_first_moment_RCT = q_S_staple_sum_only_treated_RCT + q_S_cashcrop_sum_only_treated_RCT + p_B * q_B_cashcrop_sum_only_treated_RCT;
    value_cashcrop_second_moment = q_S_cashcrop_sum_second_moment + p_B^2*q_B_cashcrop_sum_second_moment;
    #value_cashcrop_std = (value_cashcrop_second_moment -value_cashcrop_first_moment.^2)^(1/2);
    # tau_S = 0, program spending set to decrease Q_S or F_W
    #prod_value_improvement = (value_cashcrop_first_moment_RCT - value_cashcrop_first_moment)/ value_cashcrop_std;
    prod_value_improvement = (value_cashcrop_first_moment_RCT)/value_cashcrop_first_moment - 1;
    x_B_m10 = sum(x_BC_fine[:,1] .*entrants_cashcrop_from_workers + x_BC_fine[:,3] .*incumbents_cashcrop +
            x_BC_fine[:,2] .*entrants_from_staple_to_cashcrop);
    x_S_m10 = sum(x_SC_fine[:,1] .*entrants_cashcrop_from_workers + x_SC_fine[:,3] .*incumbents_cashcrop +
            x_SC_fine[:,2] .*entrants_from_staple_to_cashcrop + x_S_S_fine[:,1] .*entrants_staple_from_workers + x_S_S_fine[:,2] .*incumbents_staple +
            x_S_S_fine[:,3] .*exit_cashcrop_to_staple);
    exp_ratio_model=x_S_m10 / x_B_m10

    #    residual[10] = (exp_ratio_model - exp_ratio) / (exp_ratio_model + exp_ratio);
    mig_rate_model=sum((exit_staple_to_work + exit_cashcrop_to_work))/1; #/(current_cashcrop_pop + current_staple_pop);
    #Fixed how this is calculated in the data
    rural_pop_only_staples_model=(sum(entrants_cashcrop_from_workers.* (q_B_B_fine[:,1].==0.0)) + sum(incumbents_cashcrop.*(q_B_B_fine[:,3].==0.0))
    + sum( entrants_from_staple_to_cashcrop.*(q_B_B_fine[:,2].==0.0)) + current_staple_pop)/(current_cashcrop_pop + current_staple_pop);
    rural_pop_only_cashcrop_model=(sum(entrants_cashcrop_from_workers.* (land_C_fine[:,1].==1.0)) + sum(incumbents_cashcrop.*(land_C_fine[:,3].==1.0))
    + sum( entrants_from_staple_to_cashcrop.*(land_C_fine[:,2].==1.0)))/(current_cashcrop_pop + current_staple_pop);

    unconstrained_cashcrop_producer = sum(entrants_cashcrop_from_workers .*((TC_fine[:,1].<(κ*s_fine[:,1] .- tol) ).==1) +
    entrants_from_staple_to_cashcrop .*((TC_fine[:,2].<(κ*s_fine[:,1] .- tol)  ).==1) +
    incumbents_cashcrop .*((TC_fine[:,3].<(κ*s_fine[:,1] .- tol)) .==1) );
    fraction_cashcrop_suboptimal_model = 1 - unconstrained_cashcrop_producer/current_cashcrop_pop;

    mean_land_share_to_staples_among_cc_model =(sum(entrants_cashcrop_from_workers .* (land_C_fine[:,1].>0.0).* (1 .- land_C_fine[:,1])) +
        sum(incumbents_cashcrop.*(land_C_fine[:,3].>0.0).*(1 .- land_C_fine[:,3]))
    + sum( entrants_from_staple_to_cashcrop.*(land_C_fine[:,2].>0.0).*(1 .- land_C_fine[:,2])))/current_cashcrop_pop;

    avg_inc_rural = sum(staple_revenue_income_staple_from_workers .* entrants_staple_from_workers +
                staple_revenue_income_incumbents_staple .* incumbents_staple +
                staple_revenue_income_cashcrop_to_staple .* exit_cashcrop_to_staple +
                cashcrop_revenue_income_cashcrop_from_workers .* entrants_cashcrop_from_workers +
                cashcrop_revenue_income_staple_to_cashcrop .* entrants_from_staple_to_cashcrop +
                cashcrop_revenue_income_incumbents_cashcrop .* incumbents_cashcrop)/(current_cashcrop_pop+current_staple_pop);
    avg_inc_urban = sum(income_worker_stay_workers.*stay_workers .+ income_worker_staple_to_work.*exit_staple_to_work .+ income_worker_cashcrop_to_work.*exit_cashcrop_to_work)/current_worker_pop;

    urban_consumption_model = (p_B * (c_B_worker_sum) +
            (c_S_worker_sum + transaction_cost_worker_sum) +
        p_M * (c_M_worker_sum))/sum(current_workers);
    rural_consumption_model = (p_B * (c_B_staple_sum + c_B_Cashcrop_sum) +
            (c_S_staple_sum + c_S_cashcrop_sum + transaction_cost_staple_sum + transaction_cost_cashcrop_sum) +
        p_M * (c_M_staple_sum + c_M_cashcrop_sum))/sum((current_staple + current_cashcrop));
    # To test, I also print real consumption numbers - it potentially ignores migration too:
    urban_consumption_model_real = sum(worker_past_dist.*cons_fine_local[:,1])
    rural_consumption_model_real = sum(staple_past_dist.*cons_fine_local[:,2] + cash_crop_past_dist.*cons_fine_local[:,3])
    urban_rural_consumption_ratio_model = urban_consumption_model/rural_consumption_model
    urban_rural_consumption_ratio_model_real = urban_consumption_model_real/rural_consumption_model_real

    ### I THINK WE CAN DROP THIS MOMENT:
    fraction_selling_only_treated = sum(entrants_cashcrop_from_workers .*((c_S_B_fine[:,1] .- q_S_B_fine[:,1].<=0.0) .& (q_S_B_fine[:,1].>0.0)).* treated_past_worker +
                        incumbents_cashcrop.*((c_S_B_fine[:,3] .- q_S_B_fine[:,3].<=0.0) .& (q_S_B_fine[:,3].>0.0)).*treated_past_cashcrop +
                        entrants_from_staple_to_cashcrop.*((c_S_B_fine[:,2] .- q_S_B_fine[:,2].<=0.0) .& (q_S_B_fine[:,2].>0.0)).*treated_past_staples +
                        entrants_staple_from_workers.*(c_S_S_fine[:,1] .- q_S_S_fine[:,1].<=0.0).* treated_past_worker +
                        incumbents_staple.*(c_S_S_fine[:,2] .- q_S_S_fine[:,2].<=0.0).*treated_past_staples +
                        exit_cashcrop_to_staple.*(c_S_S_fine[:,3] .- q_S_S_fine[:,3].<=0.0).*treated_past_cashcrop);

    fraction_selling_only_treated_RCT = sum(entrants_cashcrop_from_workers .*((c_S_B_fine_RCT[:,1] .- q_S_B_fine_RCT[:,1].<=0.0) .& (q_S_B_fine_RCT[:,1].>0.0)).* treated_past_worker +
                        incumbents_cashcrop.*((c_S_B_fine_RCT[:,3] .- q_S_B_fine_RCT[:,3].<=0.0) .& (q_S_B_fine_RCT[:,3].>0.0)).*treated_past_cashcrop +
                        entrants_from_staple_to_cashcrop.*((c_S_B_fine_RCT[:,2] .- q_S_B_fine_RCT[:,2].<=0.0) .& (q_S_B_fine_RCT[:,2].>0.0)).*treated_past_staples +
                        entrants_staple_from_workers.*(c_S_S_fine_RCT[:,1] .- q_S_S_fine_RCT[:,1].<=0.0).* treated_past_worker +
                        incumbents_staple.*(c_S_S_fine_RCT[:,2] .- q_S_S_fine_RCT[:,2].<=0.0).*treated_past_staples +
                        exit_cashcrop_to_staple.*(c_S_S_fine_RCT[:,3] .- q_S_S_fine_RCT[:,3].<=0.0).*treated_past_cashcrop);

    share_selling_increase = (fraction_selling_only_treated_RCT)/fraction_selling_only_treated - 1;

    fraction_selling_any_only_treated = sum(entrants_cashcrop_from_workers .*((((c_S_B_fine[:,1] .- q_S_B_fine[:,1].<=0.0) + (c_B_B_fine[:,1] .- q_B_B_fine[:,1].<=0.0))).>0.0
    ) .* treated_past_worker +
                        incumbents_cashcrop.*(((c_S_B_fine[:,3] .- q_S_B_fine[:,3].<=0.0) + (c_B_B_fine[:,3] .- q_B_B_fine[:,3].<=0.0)).>0).*treated_past_cashcrop +
                        entrants_from_staple_to_cashcrop.*(((c_S_B_fine[:,2] .- q_S_B_fine[:,2].<=0.0) + (c_B_B_fine[:,2] .- q_B_B_fine[:,2].<=0.0)).>0).*treated_past_staples +
                        entrants_staple_from_workers.*((c_S_S_fine[:,1] .- q_S_S_fine[:,1].<=0.0).>0).* treated_past_worker +
                        incumbents_staple.*((c_S_S_fine[:,2] .- q_S_S_fine[:,2].<=0.0).>0).*treated_past_staples +
                        exit_cashcrop_to_staple.*((c_S_S_fine[:,3] .- q_S_S_fine[:,3].<=0.0).>0).*treated_past_cashcrop);

    fraction_selling_any_only_treated_RCT = sum(entrants_cashcrop_from_workers .*((((c_S_B_fine_RCT[:,1] .- q_S_B_fine_RCT[:,1].<=0.0) + (c_B_B_fine_RCT[:,1] .- q_B_B_fine_RCT[:,1].<=0.0))).>0.0
    ) .* treated_past_worker +
                        incumbents_cashcrop.*(((c_S_B_fine_RCT[:,3] .- q_S_B_fine_RCT[:,3].<=0.0) + (c_B_B_fine_RCT[:,3] .- q_B_B_fine_RCT[:,3].<=0.0)).>0).*treated_past_cashcrop +
                        entrants_from_staple_to_cashcrop.*(((c_S_B_fine_RCT[:,2] .- q_S_B_fine_RCT[:,2].<=0.0) + (c_B_B_fine_RCT[:,2] .- q_B_B_fine_RCT[:,2].<=0.0)).>0).*treated_past_staples +
                        entrants_staple_from_workers.*((c_S_S_fine_RCT[:,1] .- q_S_S_fine_RCT[:,1].<=0.0).>0).* treated_past_worker +
                        incumbents_staple.*((c_S_S_fine_RCT[:,2] .- q_S_S_fine_RCT[:,2].<=0.0).>0).*treated_past_staples +
                        exit_cashcrop_to_staple.*((c_S_S_fine_RCT[:,3] .- q_S_S_fine_RCT[:,3].<=0.0).>0).*treated_past_cashcrop);

    share_selling_any_increase = (fraction_selling_any_only_treated_RCT)/fraction_selling_any_only_treated - 1;
    fraction_borrowing_staples = sum(((x_S_S_fine[:,1]*p_x*(1+τ_S)).>s_fine[:,1]).*entrants_staple_from_workers + 
    ((x_S_S_fine[:,2]*p_x*(1+τ_S)).>s_fine[:,1]).*incumbents_staple + 
    ((x_S_S_fine[:,3]*p_x*(1+τ_S)).>s_fine[:,1]).*exit_cashcrop_to_staple ); 
    fraction_borrowing_cashcrop = sum((TC_fine[:,1].>s_fine[:,1]).*entrants_cashcrop_from_workers + 
    (TC_fine[:,2].>s_fine[:,1]).* entrants_from_staple_to_cashcrop + 
    (TC_fine[:,3].>s_fine[:,1]).*incumbents_cashcrop); 
    fraction_borrowing = (fraction_borrowing_staples + fraction_borrowing_cashcrop)

    worker_pop_effective=(urban_labor_supply_sum - total_maintenance_cost )/urban_labor_supply_sum;
    #YL_manuf=(p_M*prod_manuf - w*labor_used - (r+ δ)*capital_used)/(current_worker_pop*worker_pop_effective);
    #YL_manuf=(p_M*prod_manuf - w*labor_used)/(current_worker_pop*worker_pop_effective);
    YL_manuf=(p_M*prod_manuf - w*total_entry_cost)/(current_worker_pop*worker_pop_effective);
    YL_agr=(prod_staple + p_B*prod_cashcrop - p_x*(input_staple + input_cashcrop) - (total_maintenance_cost)*w)/(1-current_worker_pop*worker_pop_effective);
    APG=YL_manuf/YL_agr;


    transaction_cost_loss = (transaction_cost_worker_sum  + transaction_cost_staple_sum + transaction_cost_cashcrop_sum) /prod_staple;
    Phi_fine_aug = kron(Matrix(1.0I, 3, 3),row_kron(funbase(fspace_a,s_fine[:,1]),Phi_z_fine));
    welfare_val =  (Phi_fine_aug * vcat(vcat(coeff[:,1],coeff[:,2]),coeff[:,3]));
    marketable_agr_surplus_share = marketable_agr_surplus/nominal_GDP;
    exportshare_cashcrop = foreign_demand_cash_crop/prod_cashcrop;
    urban_rural_inc_ratio_model = avg_inc_urban/avg_inc_rural;
    urban_rural_wealth_ratio_model = mean_wealth_urban/mean_wealth_rural;
    avg_land_to_cashcrop = 1 - mean_land_share_to_staples_among_cc_model
    # Transaction cost - maintenance cost - migration cost comparison
    (transaction_cost_staple_sum + transaction_cost_cashcrop_sum)/ (prod_staple)
    (maintenance_costs_to_cashcrops*w)/ (p_B*prod_cashcrop)

    ### Return on fertilizers
    # Note that MPX = ζ * APX
    # MPX
    MPX_cashcrop_B_val = zeros(ns_fine*3);
    MPX_cashcrop_B_val[1:ns_fine] = ϕ_B * θ_fine .* land_C_fine[:,1].^ρ .* ζ.*x_BC_fine[:,1].^(ζ-1);
    MPX_cashcrop_B_val[(ns_fine+1):(2*ns_fine)] = ϕ_B * θ_fine .* land_C_fine[:,2].^ρ .* ζ.*x_BC_fine[:,2].^(ζ-1);
    MPX_cashcrop_B_val[(2*ns_fine+1):(3*ns_fine)] = ϕ_B * θ_fine .* land_C_fine[:,3].^ρ .* ζ.*x_BC_fine[:,3].^(ζ-1); 

    MPX_cashcrop_S_val = zeros(ns_fine*3);
    MPX_cashcrop_S_val[1:ns_fine] = ϕ_S * θ_fine .* (1 .- land_C_fine[:,1]).^ρ .* ζ.*x_SC_fine[:,1].^(ζ-1);
    MPX_cashcrop_S_val[(ns_fine+1):(2*ns_fine)] = ϕ_S * θ_fine .* (1 .- land_C_fine[:,2]).^ρ .* ζ.*x_SC_fine[:,2].^(ζ-1);
    MPX_cashcrop_S_val[(2*ns_fine+1):(3*ns_fine)] = ϕ_S * θ_fine .* (1 .- land_C_fine[:,3]).^ρ .* ζ.*x_SC_fine[:,3].^(ζ-1); 

    MPX_staples_S_val = zeros(ns_fine*3);
    MPX_staples_S_val[1:ns_fine] = ϕ_S * θ_fine .* ζ.*x_S_S_fine[:,1].^(ζ-1);
    MPX_staples_S_val[(ns_fine+1):(2*ns_fine)] = ϕ_S * θ_fine .* ζ.*x_S_S_fine[:,2].^(ζ-1);
    MPX_staples_S_val[(2*ns_fine+1):(3*ns_fine)] = ϕ_S * θ_fine .* ζ.*x_S_S_fine[:,3].^(ζ-1);
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

    second_moment_labor_prod_rural = sum(labor_prod_fine.^2 .*(current_staple + current_cashcrop )/(current_cashcrop_pop + current_staple_pop))
    second_moment_labor_prod_urban = sum(labor_prod_fine.^2 .*(current_workers )/current_worker_pop)
    second_moment_agri_prod_rural = sum(θ_fine.^2 .*(current_staple + current_cashcrop )/(current_cashcrop_pop + current_staple_pop))
    second_moment_agri_prod_urban = sum(θ_fine.^2 .*(current_workers )/current_worker_pop)

    std_labor_prod_rural = (second_moment_labor_prod_rural - avg_labor_prod_rural.^2).^(1/2);
    std_labor_prod_urban = (second_moment_labor_prod_urban - avg_labor_prod_urban.^2).^(1/2);
    std_agri_prod_rural = (second_moment_agri_prod_rural - avg_agri_prod_rural.^2).^(1/2);
    std_agri_prod_urban = (second_moment_agri_prod_urban- avg_agri_prod_urban.^2).^(1/2);

    coeff_var_labor_prod_rural = std_labor_prod_rural/avg_labor_prod_rural;
    coeff_var_labor_prod_urban = std_labor_prod_urban/avg_labor_prod_urban;
    coeff_var_agri_prod_rural = std_agri_prod_rural/avg_agri_prod_rural;
    coeff_var_agri_prod_urban = std_agri_prod_urban/avg_agri_prod_urban;

    # var_MPX_cashcrop=var(ϕ_B * θ_fine .* land_C_fine[:,1].^ρ .* ζ.*x_BC_fine[:,1].^(ζ-1) .* entrants_cashcrop_from_workers/current_cashcrop_pop +
    #     ϕ_B * θ_fine .* land_C_fine[:,2].^ρ .* ζ.*x_BC_fine[:,2].^(ζ-1) .* entrants_from_staple_to_cashcrop/current_cashcrop_pop +
    #     ϕ_B * θ_fine .* land_C_fine[:,3].^ρ .* ζ.*x_BC_fine[:,3].^(ζ-1) .* incumbents_cashcrop/current_cashcrop_pop +
    #     ϕ_S * θ_fine .* (1 .- land_C_fine[:,1]).^ρ .* ζ.*x_SC_fine[:,1].^(ζ-1) .* entrants_cashcrop_from_workers/current_cashcrop_pop +
    #     ϕ_S * θ_fine .* (1 .- land_C_fine[:,2]).^ρ .* ζ.*x_SC_fine[:,2].^(ζ-1) .* entrants_from_staple_to_cashcrop/current_cashcrop_pop +
    #     ϕ_S * θ_fine .* (1 .- land_C_fine[:,3]).^ρ .* ζ.*x_SC_fine[:,3].^(ζ-1) .* incumbents_cashcrop/current_cashcrop_pop);
    # var_MPX_staple=var(ϕ_S *θ_fine .*  ζ.*x_S_S_fine[:,1].^(ζ-1) .* entrants_cashcrop_from_workers/current_staple_pop +
    #     ϕ_S * θ_fine .* ζ.*x_S_S_fine[:,2].^(ζ-1) .* incumbents_staple/current_staple_pop +
    #     ϕ_S * θ_fine .* ζ.*x_S_S_fine[:,3].^(ζ-1) .* exit_cashcrop_to_staple/current_staple_pop);
    # var_MPL=var(labor_prod_fine .* stay_workers/current_worker_pop +
    #             labor_prod_fine .* exit_staple_to_work/current_worker_pop +
    #             labor_prod_fine .* exit_cashcrop_to_work/current_worker_pop)
    relative_land_to_cashcrop = avg_land_to_cashcrop *current_cashcrop_pop/(current_cashcrop_pop + current_staple_pop);
    relative_land_to_staples = 1 - relative_land_to_cashcrop;
    avg_land_to_staples = 1 - avg_land_to_cashcrop;
    staple_productivity = prod_staple/(avg_land_to_staples *current_cashcrop_pop  + current_staple_pop);
    cashcrop_productivity = prod_cashcrop/(avg_land_to_cashcrop *current_cashcrop_pop);
    manuf_productivity = prod_manuf/(current_worker_pop);

    constrained_cashcrop_producer = sum(entrants_cashcrop_from_workers .*(((TC_fine[:,1].<=(κ*s_fine[:,1] .+ 2*tol)) + (TC_fine[:,1].>=(κ*s_fine[:,1] .- 2*tol) )  ).==2) +
    entrants_from_staple_to_cashcrop .*(((TC_fine[:,2].<=(κ*s_fine[:,1] .+ 2*tol)) + (TC_fine[:,2].>=(κ*s_fine[:,1] .- 2*tol) )  ).==2) +
    incumbents_cashcrop .*(((TC_fine[:,3].<=(κ*s_fine[:,1] .+ 2*tol)) + (TC_fine[:,3].>=(κ*s_fine[:,1] .- 2*tol) )  ).==2));
    share_constrained_cashcrop = constrained_cashcrop_producer/current_cashcrop_pop;

    TC_staple_fine=p_x.*(1+τ_S).*x_S_S_fine;
    constrained_staple_producer = sum(entrants_staple_from_workers .*(((TC_staple_fine[:,1].<=(κ*s_fine[:,1] .+ 2*tol)) + (TC_staple_fine[:,1].>=(κ*s_fine[:,1] .- 2*tol) )  ).==2) +
    incumbents_staple .*(((TC_staple_fine[:,2].<=(κ*s_fine[:,1] .+ 2*tol)) + (TC_staple_fine[:,2].>=(κ*s_fine[:,1] .- 2*tol) )  ).==2) +
    exit_cashcrop_to_staple .*(((TC_staple_fine[:,3].<=(κ*s_fine[:,1] .+ 2*tol)) + (TC_staple_fine[:,3].>=(κ*s_fine[:,1] .- 2*tol) )  ).==2));
    share_constrained_staple = constrained_staple_producer/current_staple_pop;
    aggregate_consumption = sum(cons_fine_local[:,1] .* stat_distr[1:ns_fine] +
    cons_fine_local[:,2] .* stat_distr[(ns_fine + 1):(2*ns_fine)] +
    cons_fine_local[:,3] .* stat_distr[(2*ns_fine + 1):(3*ns_fine)])

    worker_pop_effective=(urban_labor_supply_sum - total_maintenance_cost )/urban_labor_supply_sum;
    #YL_manuf=(p_M*prod_manuf - w*labor_used - (r+ δ)*capital_used)/(current_worker_pop*worker_pop_effective);
    #YL_manuf=(p_M*prod_manuf - w*labor_used)/(current_worker_pop*worker_pop_effective);
    YL_manuf= 1/P_W* (p_M*prod_manuf - w*total_entry_cost)/(current_worker_pop*worker_pop_effective);
    YL_agr= 1/P_W*(prod_staple + p_B*prod_cashcrop - p_x*(input_staple + input_cashcrop) - (total_maintenance_cost)*w)/(1-current_worker_pop*worker_pop_effective);
    APG=YL_manuf/YL_agr;
    TFP =  YL_agr * (1-current_worker_pop*worker_pop_effective) * (ψ_S + ψ_M) + YL_manuf *  ψ_B* (current_worker_pop*worker_pop_effective)

    if out == 3
        # Experimental stuff
        future_occupation_fine_local_int = convert.(Int64,future_occupation_fine_local);
        producer_price = zeros(ns_fine*6);
        producer_distr = zeros(ns_fine*6);
        producer_price[1:ns_fine] = λ_2_fine[:,1];
        producer_distr[1:ns_fine] = entrants_cashcrop_from_workers;
        producer_price[ (ns_fine+1):(2*ns_fine)] = λ_2_fine[:,2];
        producer_distr[ (ns_fine+1):(2*ns_fine)] = entrants_from_staple_to_cashcrop;
        producer_price[(2*ns_fine+1):(3*ns_fine)] = λ_2_fine[:,3];
        producer_distr[(2*ns_fine+1):(3*ns_fine)] = incumbents_cashcrop;
        producer_price[ (3*ns_fine+1):(4*ns_fine)] = λ_2_S_fine[:,1];
        producer_distr[ (3*ns_fine+1):(4*ns_fine)] = entrants_staple_from_workers;
        producer_price[(4*ns_fine+1):(5*ns_fine)] = λ_2_S_fine[:,2];
        producer_distr[(4*ns_fine+1):(5*ns_fine)] = incumbents_staple;
        producer_price[(5*ns_fine+1):(6*ns_fine) ] = λ_2_S_fine[:,3];
        producer_distr[(5*ns_fine+1):(6*ns_fine) ] = exit_cashcrop_to_staple;
    
        ratio_price = zeros(ns_fine*6);
        ratio_price2 = zeros(ns_fine*6);
        ratio_distr = zeros(ns_fine*6);
        ratio_price[1:ns_fine] = ( q_S_B_fine[:, 1] + p_B * q_B_B_fine[:, 1] ) ./ (λ_2_fine[:, 1] .* q_S_B_fine[:,1] + p_B*q_B_B_fine[:,1] )
        #ratio_price2[1:ns_fine] = (λ_2_fine[:, 1] .* q_S_B_fine[:,1] + p_B*q_B_B_fine[:,1] ) ./ ((1+Q_S) .* q_S_B_fine[:,1] + p_B*q_B_B_fine[:,1] )
        ratio_price2[1:ns_fine] = (λ_2_fine[:, 1] .* q_S_B_fine[:,1] ) ./ ((1+Q_S) .* q_S_B_fine[:,1]  )
        ratio_distr[1:ns_fine] = entrants_cashcrop_from_workers;
        ratio_price[ (ns_fine+1):(2*ns_fine)] = ( q_S_B_fine[:, 2] + p_B * q_B_B_fine[:, 2] ) ./ (λ_2_fine[:, 2] .* q_S_B_fine[:,2] + p_B*q_B_B_fine[:,2] ) #λ_2_fine[:,2];
        #ratio_price2[ (ns_fine+1):(2*ns_fine)] = (λ_2_fine[:, 2] .* q_S_B_fine[:,2] + p_B*q_B_B_fine[:,2] ) ./ ((1+Q_S) .* q_S_B_fine[:,2] + p_B*q_B_B_fine[:,2] ) #λ_2_fine[:,2];
        ratio_price2[ (ns_fine+1):(2*ns_fine)] = (λ_2_fine[:, 2] .* q_S_B_fine[:,2] ) ./ ((1+Q_S) .* q_S_B_fine[:,2]) #λ_2_fine[:,2];
        ratio_distr[ (ns_fine+1):(2*ns_fine)] = entrants_from_staple_to_cashcrop;
        ratio_price[(2*ns_fine+1):(3*ns_fine)] = ( q_S_B_fine[:, 3] + p_B * q_B_B_fine[:, 3] ) ./ (λ_2_fine[:, 3] .* q_S_B_fine[:,3] + p_B*q_B_B_fine[:,3] ) #λ_2_fine[:,3];
        #ratio_price2[(2*ns_fine+1):(3*ns_fine)] = (λ_2_fine[:, 3] .* q_S_B_fine[:,3] + p_B*q_B_B_fine[:,3] ) ./ ((1+Q_S) .* q_S_B_fine[:,3] + p_B*q_B_B_fine[:,3] ) #λ_2_fine[:,3];
        ratio_price2[(2*ns_fine+1):(3*ns_fine)] = (λ_2_fine[:, 3] .* q_S_B_fine[:,3] ) ./ ((1+Q_S) .* q_S_B_fine[:,3]  )
        ratio_distr[(2*ns_fine+1):(3*ns_fine)] = incumbents_cashcrop;
        ratio_price[ (3*ns_fine+1):(4*ns_fine)] = ( q_S_S_fine[:, 1] ) ./ (λ_2_S_fine[:, 1] .* q_S_S_fine[:,1] ) #λ_2_S_fine[:,1];
        ratio_price2[ (3*ns_fine+1):(4*ns_fine)] = (λ_2_S_fine[:, 1] .* q_S_S_fine[:,1] ) ./ ((1+Q_S) .* q_S_S_fine[:,1] ) #λ_2_S_fine[:,1];
        ratio_distr[ (3*ns_fine+1):(4*ns_fine)] = entrants_staple_from_workers;
        ratio_price[(4*ns_fine+1):(5*ns_fine)] = ( q_S_S_fine[:, 2] ) ./ (λ_2_S_fine[:, 2] .* q_S_S_fine[:,2] ) #λ_2_S_fine[:,2];
        ratio_price2[(4*ns_fine+1):(5*ns_fine)] = (λ_2_S_fine[:, 2] .* q_S_S_fine[:,2] ) ./ ((1+Q_S) .* q_S_S_fine[:,2] ) #λ_2_S_fine[:,2];
        ratio_distr[(4*ns_fine+1):(5*ns_fine)] = incumbents_staple;
        ratio_price[(5*ns_fine+1):(6*ns_fine) ] = ( q_S_S_fine[:, 3] ) ./ (λ_2_S_fine[:, 3] .* q_S_S_fine[:,3] ) #λ_2_S_fine[:,3];
        ratio_price2[(5*ns_fine+1):(6*ns_fine) ] = (λ_2_S_fine[:, 3] .* q_S_S_fine[:,3] ) ./ ((1+Q_S) .* q_S_S_fine[:,3] )#λ_2_S_fine[:,3];
        ratio_distr[(5*ns_fine+1):(6*ns_fine) ] = exit_cashcrop_to_staple;
        #shadow_price_sort_index = sortperm(producer_price);
        #shadow_price_sorted = producer_price[shadow_price_sort_index];
        #producer_distr_sorted = producer_distr[shadow_price_sort_index];
        #cumu_producer_distr =  cumsum(producer_distr_sorted,dims = 1);
        #h = fit(Histogram,shadow_price_sorted, weights(producer_distr_sorted), nbins = 10 )
        #plot(h)
        #savefig("histogram2.pdf")
        # - correlation between:
        #- share of maize and log_value
        land_share_maize = ones(ns_fine*6);
        land_share_maize[(0*ns_fine+1):(1*ns_fine)] =   1 .- land_C_fine[:,1];
        land_share_maize[(1*ns_fine+1):(2*ns_fine)] =   1 .- land_C_fine[:,2];
        land_share_maize[(2*ns_fine+1):(3*ns_fine)] =   1 .- land_C_fine[:,3];
    
        revenues = ones(ns_fine*6);
        revenues[(0*ns_fine+1):(1*ns_fine)] =   revenue_cashcrop[:,1];
        revenues[(1*ns_fine+1):(2*ns_fine)] =   revenue_cashcrop[:,2];
        revenues[(2*ns_fine+1):(3*ns_fine)] =   revenue_cashcrop[:,3];
        revenues[(3*ns_fine+1):(4*ns_fine)] =   revenue_staples[:,1];
        revenues[(4*ns_fine+1):(5*ns_fine)] =   revenue_staples[:,2];
        revenues[(5*ns_fine+1):(6*ns_fine)] =   revenue_staples[:,3];
        # Clear observations without any revenues: 
        revenue_cleared_index = (revenues.!=0.0);
     
        
        values_harvest = ones(ns_fine*6);
        values_harvest[(0*ns_fine+1):(1*ns_fine)] =   q_S_B_fine[:,1] + p_B*q_B_B_fine[:,1];
        values_harvest[(1*ns_fine+1):(2*ns_fine)] =   q_S_B_fine[:,2] + p_B*q_B_B_fine[:,2];
        values_harvest[(2*ns_fine+1):(3*ns_fine)] =   q_S_B_fine[:,3] + p_B*q_B_B_fine[:,3];
        values_harvest[(3*ns_fine+1):(4*ns_fine)] =   q_S_S_fine[:,1];
        values_harvest[(4*ns_fine+1):(5*ns_fine)] =   q_S_S_fine[:,2];
        values_harvest[(5*ns_fine+1):(6*ns_fine)] =   q_S_S_fine[:,3];
    
    
        # Self consumed:
        self_consumed = 1 .- revenues./values_harvest;
    
        subsidized_dummy = ones(ns_fine*6);
        threshold_subsidized = 0.0555;
        subsidized_dummy[(0*ns_fine+1):(1*ns_fine)] =   threshold_subsidized .< (x_SC_fine[:,1] * p_x * (1+τ_S));
        subsidized_dummy[(1*ns_fine+1):(2*ns_fine)] =   threshold_subsidized .< (x_SC_fine[:,2] * p_x * (1+τ_S));
        subsidized_dummy[(2*ns_fine+1):(3*ns_fine)] =   threshold_subsidized .< (x_SC_fine[:,3] * p_x * (1+τ_S));
        subsidized_dummy[(3*ns_fine+1):(4*ns_fine)] =   threshold_subsidized .< (x_S_S_fine[:,1] * p_x * (1+τ_S));
        subsidized_dummy[(4*ns_fine+1):(5*ns_fine)] =   threshold_subsidized .< (x_S_S_fine[:,2] * p_x * (1+τ_S));
        subsidized_dummy[(5*ns_fine+1):(6*ns_fine)] =   threshold_subsidized .< (x_S_S_fine[:,3] * p_x * (1+τ_S));
    
        #Fertilizer use:
        fertilizer_use = ones(ns_fine*6);
        fertilizer_use[(0*ns_fine+1):(1*ns_fine)] =   x_SC_fine[:,1] + x_BC_fine[:,1];
        fertilizer_use[(1*ns_fine+1):(2*ns_fine)] =   x_SC_fine[:,2] + x_BC_fine[:,2];
        fertilizer_use[(2*ns_fine+1):(3*ns_fine)] =   x_SC_fine[:,3] + x_BC_fine[:,3];
        fertilizer_use[(3*ns_fine+1):(4*ns_fine)] =   x_S_S_fine[:,1];
        fertilizer_use[(4*ns_fine+1):(5*ns_fine)] =   x_S_S_fine[:,2];
        fertilizer_use[(5*ns_fine+1):(6*ns_fine)] =   x_S_S_fine[:,3];
    
        #Productivity of households:
        θ_fine_stacked = [θ_fine
        θ_fine
        θ_fine
        θ_fine
        θ_fine
        θ_fine]
        #Wealth of households:
        a_fine_stacked = [s_fine[:,1]
        s_fine[:,1]
        s_fine[:,1]
        s_fine[:,1]
        s_fine[:,1]
        s_fine[:,1]]
        # Get the sample of households
        prob_distr = [entrants_cashcrop_from_workers 
        entrants_from_staple_to_cashcrop 
        incumbents_cashcrop 
        entrants_staple_from_workers
        incumbents_staple
        exit_cashcrop_to_staple ]./(1- current_worker_pop)
    #    prob_distr_cleared = prob_distr[revenue_cleared_index];
    #    subsidized_dummy_cleared = subsidized_dummy[revenue_cleared_index];
    #    revenues_cleared = revenues[revenue_cleared_index];
    #    land_share_maize_cleared = land_share_maize[revenue_cleared_index];
    #    producer_price_cleared = producer_price[revenue_cleared_index];
    #    prob_distr_cleared = prob_distr_cleared./sum(prob_distr_cleared);
    #weights = Weights(prob_distr_cleared)
        weights = Weights(prob_distr)
    
        #N_sample =8544;#4279#5541 #
        N_sample_FISP =4879;#4279#5541 #
        N_sample_NoFISP = 3729#4279#5541 #
            N_sample = 500000#N_sample_NoFISP
    
            Atmp=[subsidized_dummy revenues land_share_maize producer_price fertilizer_use values_harvest self_consumed θ_fine_stacked a_fine_stacked ratio_price ratio_price2]
            indx = sample(axes(Atmp, 1), weights, N_sample)
            sampled_data=zeros(N_sample,11);
            sampled_data=Atmp[indx, :]
        
        # for ii=1:N_sample
        #     #sampled_data[ii,1]=sample(subsidized_dummy_cleared,weights)
        #     #sampled_data[ii,2]=sample(revenues_cleared,weights)
        #     #sampled_data[ii,3]=sample(land_share_maize_cleared,weights)
        #     #sampled_data[ii,4]=sample(producer_price_cleared,weights)
        #     sampled_data[ii,1]=sample(subsidized_dummy,weights)
        #     sampled_data[ii,2]=sample(revenues,weights)
        #     sampled_data[ii,3]=sample(land_share_maize,weights)
        #     sampled_data[ii,4]=sample(producer_price,weights)
        #     sampled_data[ii,5]=sample(fertilizer_use,weights)
        #     sampled_data[ii,6]=sample(values_harvest,weights)
        #     sampled_data[ii,7]=sample(self_consumed,weights)
        #     sampled_data[ii,8]=sample(θ_fine_stacked,weights)
        #     sampled_data[ii,9]=sample(a_fine_stacked,weights)
        #     sampled_data[ii,10]=sample(ratio_price,weights)
        #     #sampled_data[ii,5]=sample(repeat(kron(z,ones(n_fine[1])),size(l_z)[1]*3),weights) #rural productivy
        # end
    
        # Check the fraction of people "without subsidy":
        sum(sampled_data[:,1])/N_sample
    
        df_validate = DataFrame(share_maize=sampled_data[:, 3], value_prod=sampled_data[:, 6], share_selfconsumed=sampled_data[:, 7], subsidized_dummy = sampled_data[:,1], 
             inputs_used = sampled_data[:,5], producer_price_at_gate = sampled_data[:,4], price_ratio = sampled_data[:,10], price_ratio2 = sampled_data[:,11], theta = sampled_data[:,8], a = sampled_data[:,9])
    
        # if N_sample==N_sample_FISP
        #     df_validate[!, :ind_sub]=ones(N_sample)
        #     df_validate_FISP=copy(df_validate)
        # else
        #     df_validate[!, :ind_sub]=zeros(N_sample)
        #     df_validate_NoFISP=copy(df_validate)  
        # end
    
        ### Comment out and rUN THE ABOVE CHUNK TWICE FOR FISP AND NO FISP if you want all REGRESSIONS BELOW:
        df_validate_FISP = copy(df_validate)
        df_validate_FISP[!, :ind_sub] .= 1#ones(N_sample_FISP)
        # df_validate_NoFISP_PE[!, :ind_sub] = zeros(N_sample_NoFISP)
        # df_validate_NoFISP_GE[!, :ind_sub] = zeros(N_sample_NoFISP)
    
        # df_validateFISP_NoFISP_GE=[df_validate_FISP; df_validate_NoFISP_GE]
        # df_validateFISP_NoFISP_PE=[df_validate_FISP; df_validate_NoFISP_PE]
    
        # df_validateFISP_NoFISP_GE[!, :nind_sub] = df_validateFISP_NoFISP_GE.ind_sub / mean(df_validateFISP_NoFISP_GE.ind_sub)
        # df_validateFISP_NoFISP_PE[!, :nind_sub] = df_validateFISP_NoFISP_PE.ind_sub / mean(df_validateFISP_NoFISP_PE.ind_sub)
    
        df_validate_FISP[!, :nshare_maize] = df_validate_FISP.share_maize / mean(df_validate_FISP.share_maize)
        # df_validateFISP_NoFISP_GE[!, :nshare_maize] = df_validateFISP_NoFISP_GE.share_maize / mean(df_validateFISP_NoFISP_GE.share_maize)
        # df_validateFISP_NoFISP_PE[!, :nshare_maize] = df_validateFISP_NoFISP_PE.share_maize / mean(df_validateFISP_NoFISP_PE.share_maize)
    
        df_validate_FISP[!, :nvalue_prod] = df_validate_FISP.value_prod / mean(df_validate_FISP.value_prod)
        # df_validateFISP_NoFISP_GE[!, :nvalue_prod] = df_validateFISP_NoFISP_GE.value_prod / mean(df_validateFISP_NoFISP_GE.value_prod)
        # df_validateFISP_NoFISP_PE[!, :nvalue_prod] = df_validateFISP_NoFISP_PE.value_prod / mean(df_validateFISP_NoFISP_PE.value_prod)
    
        df_validate_FISP[!, :nshare_selfconsumed] = df_validate_FISP.share_selfconsumed / mean(df_validate_FISP.share_selfconsumed)
        # df_validateFISP_NoFISP_GE[!, :nshare_selfconsumed] = df_validateFISP_NoFISP_GE.share_selfconsumed / mean(df_validateFISP_NoFISP_GE.share_selfconsumed)
        # df_validateFISP_NoFISP_PE[!, :nshare_selfconsumed] = df_validateFISP_NoFISP_PE.share_selfconsumed / mean(df_validateFISP_NoFISP_PE.share_selfconsumed)
    
        df_validate_FISP[!, :ninputs_used] = df_validate_FISP.inputs_used / mean(df_validate_FISP.inputs_used)
        # df_validateFISP_NoFISP_GE[!, :ninputs_used] = df_validateFISP_NoFISP_GE.inputs_used / mean(df_validateFISP_NoFISP_GE.inputs_used)
        # df_validateFISP_NoFISP_PE[!, :ninputs_used] = df_validateFISP_NoFISP_PE.inputs_used / mean(df_validateFISP_NoFISP_PE.inputs_used)
      
        df_validate_FISP[!, :nprice_ratio] = df_validate_FISP.price_ratio / mean(df_validate_FISP.price_ratio)
        # df_validateFISP_NoFISP_GE[!, :nprice_ratio] = df_validateFISP_NoFISP_GE.price_ratio / mean(df_validateFISP_NoFISP_GE.price_ratio)
        # df_validateFISP_NoFISP_PE[!, :nprice_ratio] = df_validateFISP_NoFISP_PE.price_ratio / mean(df_validateFISP_NoFISP_PE.price_ratio)
    
        df_validate_FISP[!, :nprice_ratio2] = df_validate_FISP.price_ratio2 / mean(df_validate_FISP.price_ratio2)
        # df_validateFISP_NoFISP_GE[!, :nprice_ratio2] = df_validateFISP_NoFISP_GE.price_ratio2 / mean(df_validateFISP_NoFISP_GE.price_ratio2)
        # df_validateFISP_NoFISP_PE[!, :nprice_ratio2] = df_validateFISP_NoFISP_PE.price_ratio2 / mean(df_validateFISP_NoFISP_PE.price_ratio2)
    
        df_validate_FISP[!, :nproducer_price_at_gate] = df_validate_FISP.producer_price_at_gate / mean(df_validate_FISP.producer_price_at_gate)
        # df_validateFISP_NoFISP_GE[!, :nproducer_price_at_gate] = df_validateFISP_NoFISP_GE.producer_price_at_gate / mean(df_validateFISP_NoFISP_GE.producer_price_at_gate)
        # df_validateFISP_NoFISP_PE[!, :nproducer_price_at_gate] = df_validateFISP_NoFISP_PE.producer_price_at_gate / mean(df_validateFISP_NoFISP_PE.producer_price_at_gate)
    
        df_validate_FISP[!, :na] = df_validate_FISP.a / mean(df_validate_FISP.a)
    
    
        # reg1PE = lm(@formula(nshare_maize  ~ 1+ nind_sub), df_validateFISP_NoFISP_PE);
        # reg2PE = lm(@formula(nvalue_prod ~ 1+ nind_sub ), df_validateFISP_NoFISP_PE); 
        # reg3PE = lm(@formula(ninputs_used ~ 1 + nind_sub), df_validateFISP_NoFISP_PE);
        # reg4PE = lm(@formula(nshare_selfconsumed ~ 1 + nind_sub), df_validateFISP_NoFISP_PE);
        # reg5PE = lm(@formula(nproducer_price_at_gate ~ 1 + nind_sub), df_validateFISP_NoFISP_PE);
        # reg6PE = lm(@formula(nprice_ratio ~ 1 + nind_sub), df_validateFISP_NoFISP_PE);
        # reg7PE = lm(@formula(nprice_ratio2 ~ 1 + nind_sub), df_validateFISP_NoFISP_PE);
    
        # reg1GE = √;
        # reg2GE = lm(@formula(nvalue_prod ~ 1+ nind_sub ), df_validateFISP_NoFISP_GE); 
        # reg3GE = lm(@formula(ninputs_used ~ 1 + nind_sub), df_validateFISP_NoFISP_GE);
        # reg4GE = lm(@formula(nshare_selfconsumed ~ 1 + nind_sub), df_validateFISP_NoFISP_GE);
        # reg5GE = lm(@formula(nproducer_price_at_gate ~ 1 + nind_sub), df_validateFISP_NoFISP_GE);
        # reg6GE = lm(@formula(nprice_ratio ~ 1 + nind_sub), df_validateFISP_NoFISP_GE);
        # reg7GE = lm(@formula(nprice_ratio2 ~ 1 + nind_sub), df_validateFISP_NoFISP_GE);
    
        #regs1PE = lm(@formula(nshare_maize  ~ 1+ nshare_maize), df_validateFISP_NoFISP_PE);
        regs2only = lm(@formula(nvalue_prod ~ 1+ nshare_maize ), df_validate_FISP);
        regs3only = lm(@formula(ninputs_used ~ 1 + nshare_maize), df_validate_FISP);
        regs4only = lm(@formula(nshare_selfconsumed ~ 1 + nshare_maize), df_validate_FISP);
        regs5only = lm(@formula(nproducer_price_at_gate ~ 1 + nshare_maize), df_validate_FISP);
        regs6only = lm(@formula(nprice_ratio ~ 1 + nshare_maize), df_validate_FISP);
        regs7only = lm(@formula(nprice_ratio2 ~ 1 + nshare_maize), df_validate_FISP);
    
        #     regs2onlya = lm(@formula(nvalue_prod ~ 1+ nshare_maize  + a), df_validate_FISP);
        # regs3onlya = lm(@formula(ninputs_used ~ 1 + nshare_maize + a), df_validate_FISP);
        # regs4onlya = lm(@formula(nshare_selfconsumed ~ 1 + nshare_maize + a), df_validate_FISP);
        # regs5onlya = lm(@formula(nproducer_price_at_gate ~ 1 + nshare_maize + a), df_validate_FISP);
        # regs6onlya = lm(@formula(nprice_ratio ~ 1 + nshare_maize + a), df_validate_FISP);
        # regs7onlya = lm(@formula(nprice_ratio2 ~ 1 + nshare_maize + a), df_validate_FISP);
    
    
    
    
        #     reg1 = lm(@formula(share_maize  ~ 1+ value_prod + ind_sub), df_validate);
        # reg2 = lm(@formula(value_prod ~ 1+ share_maize + ind_sub ), df_validate); 
        # reg3 = lm(@formula(inputs_used ~ 1 + share_maize + value_prod + ind_sub), df_validate);
        # reg4 = lm(@formula(share_selfconsumed ~ 1 + share_maize + value_prod + ind_sub + inputs_used + household_land + household_wealth ), df_validate);
        # reg5 = lm(@formula(producer_price_at_gate ~ 1 + share_maize + value_prod + ind_sub), df_validate);
        # reg5 = lm(@formula(nprice_ratio ~ 1 + share_maize + value_prod + ind_sub), df_validate);
        # reg5 = lm(@formula(nprice_ratio2 ~ 1 + share_maize + value_prod + ind_sub), df_validate);
     
        # reg1 , r2(reg1)
        # reg2 , r2(reg2)
        # reg3 , r2(reg3)
        # reg4 , r2(reg4)
        # reg5 , r2(reg5)
        #lm(@formula(share_selfconsumed ~ 1 + share_maize + share_maize^2 + value_prod + value_prod^2 + FISP_recipient + inputs_used + household_land + household_wealth ), df_validate)
        #lm(@formula(share_selfconsumed ~ 1+ FISP_recipient ), df_validate)
        #lm(@formula(share_maize  ~ 1+ value_prod + FISP_recipient + household_land + household_wealth + ), df_validate);
        #lm(@formula(share_maize  ~ 1+ value_prod + value_prod^2 + FISP_recipient), df_validate)
        #sum(producer_price.*producer_distr./(current_staple_pop + current_cashcrop_pop))
        #sampled_data_FISP = copy(sampled_data)
        #sampled_data_FISP[:,1] .= 1
        #sampled_data_noFISP = copy(sampled_data)
        #sampled_data_noFISP[:,1] .= 0
        #sampled_data_combined = [sampled_data_FISP 
        #sampled_data_noFISP]
        end
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


        # sum((cons_fine_local .> consumption_threshold) .* [worker_past_dist staple_past_dist cash_crop_past_dist])
        # sum((income_policy .> income_threshold) .* [worker_past_dist 
        # staple_past_dist 
        # cash_crop_past_dist])
        # # Subsistence levels 
        # substinence_factor = 1.25
        # substinence_level = substinence_factor * c̄_S;
        # substinence_constrained= 
        # c_S_S_fine[:,1] .*entrants_staple_from_workers + c_S_S_fine[:,2] .*incumbents_staple
        # + c_S_S_fine[:,3] .*exit_cashcrop_to_staple
        # c_S_B_fine[:,1] .*entrants_cashcrop_from_workers + c_S_B_fine[:,3] .*incumbents_cashcrop
        #     + c_S_B_fine[:,2] .*entrants_from_staple_to_cashcrop
        # (c̄_S .+ (1 + Q_S)^(-ϵ)*ψ_S.^ϵ.*cons_fine_local[:,1].*P_W_fine.^ϵ) .*stay_workers +(c̄_S .+ (1 + Q_S)^(-ϵ)*ψ_S.^ϵ.*cons_fine_local[:,2].*P_W_fine.^ϵ) .*exit_staple_to_work +
        #     (c̄_S .+ (1 + Q_S)^(-ϵ)*ψ_S.^ϵ.*cons_fine_local[:,3].*P_W_fine.^ϵ) .*exit_cashcrop_to_work

    return (residual_goods, stat_distr, cons_fine_local, a_prime_fine_local,future_occupation_fine_local,x_S_S_fine,x_SC_fine,x_BC_fine, coeff,
        transaction_cost_loss,nominal_GDP,welfare_val,Import_value,Export_value,current_worker_pop,current_staple_pop,current_cashcrop_pop,marketable_agr_surplus_share,exportshare_cashcrop,
        fraction_model,program_spending,prod_value_improvement,share_selling_increase,exp_ratio_model,mig_rate_model,rural_pop_only_staples_model,rural_pop_only_cashcrop_model,
        mean_land_share_to_staples_among_cc_model,urban_rural_inc_ratio_model,urban_rural_wealth_ratio_model,urban_rural_consumption_ratio_model,
        p90_wealth_rural,p90_wealth_urban,p99_wealth_rural,p99_wealth_urban,p90_cons_tmp,p90_income_tmp,p99_cons_tmp,p99_income_tmp,
        staple_productivity,cashcrop_productivity,manuf_productivity,relative_land_to_staples,relative_land_to_cashcrop,share_constrained_cashcrop,var_MPX_staples_S,var_MPX_cashcrop_B,var_MPX_cashcrop_S,
        share_constrained_staple,APG,urban_rural_consumption_ratio_model_real,aggregate_consumption,worker_pop_effective,prod_manuf,total_entry_cost,prod_staple,
        prod_cashcrop,input_staple,input_cashcrop,total_maintenance_cost,current_account_residual,fraction_cashcrop_suboptimal_model,V_saved_local,avg_labor_prod_rural,avg_labor_prod_urban,
        avg_agri_prod_rural,avg_agri_prod_urban,var_MPX_cashcrop,var_MPX,TFP,YL_manuf,YL_agr,coeff_var_labor_prod_rural,coeff_var_labor_prod_urban,coeff_var_agri_prod_rural, coeff_var_agri_prod_urban,
        p90_wealth_tmp,p99_wealth_tmp,p90_cons_tmp_rural,p99_cons_tmp_rural,p90_cons_tmp_urban,p99_cons_tmp_urban,p90_income_tmp_rural,p99_income_tmp_rural,p90_income_tmp_urban,p99_income_tmp_urban,
        wealth_of_workers,wealth_of_staples,wealth_of_cashcrop,c_B_worker_sum,c_B_staple_sum,c_B_Cashcrop_sum,c_S_worker_sum,c_S_staple_sum ,c_S_cashcrop_sum,
        transaction_cost_staple_sum,transaction_cost_cashcrop_sum,transaction_cost_worker_sum,c_M_worker_sum,c_M_staple_sum,c_M_cashcrop_sum,
        MPX_mean_log, MPX_mean_staples_S_log,MPX_mean_cashcrop_log
        , APland_mean_log,APland_mean_cashcrop_log, APland_mean_staples_S_log,var_APland,var_APland_cashcrop,var_APland_staples_S,
        c_S_W_fine,c_B_W_fine,c_M_W_fine,c_S_S_fine,c_B_S_fine,c_M_S_fine,c_S_B_fine,c_B_B_fine,c_M_B_fine)
end
end
