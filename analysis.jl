#include("main.jl")
##########################################################################################################################################################
##
##          ANALYSIS: additional data construction
##
##########################################################################################################################################################


function mat_creator(y_value)
    y_mat = zeros(Baseline_parameter.n_fine[2],Baseline_parameter.n_fine[1]);
    for x_index = 1:Baseline_parameter.n_fine[2]
        for y_index=1:Baseline_parameter.n_fine[1]
            y_mat[x_index,y_index] = y_value[Baseline_parameter.n_fine[1]*(x_index-1) + y_index];
        end
    end
    return y_mat
end
function movmean(array::Array{Float64,1},window::Int64)
    array_size = size(array)[1];
    array_smooth = zeros(array_size);
    for i = 1:array_size
        if i<window
            array_smooth[i] = sum(array[1:i])/i
        else
            array_smooth[i] = sum(array[(i-window+1):i])/window
        end
    end
    return array_smooth
end
using Plots
#using PlotlyJS
#using LaTeXStrings
#;plotlyjs() using PlotlyJS
#using LaTeXStrings

R = r + Baseline_parameter.δ;
w_no_subsidy = prices_no_subsidy[2]*(1-Baseline_parameter.α)/(1+0)*(R/Baseline_parameter.α*1/prices_no_subsidy[2])^(Baseline_parameter.α/(Baseline_parameter.α-1));
w_subsidy_nb = prices_subsidy_nb[2]*(1-Baseline_parameter.α)/(1+0)*(R/Baseline_parameter.α*1/prices_subsidy_nb[2])^(Baseline_parameter.α/(Baseline_parameter.α-1));
w_subsidy_b = prices_subsidy_b[2]*(1-Baseline_parameter.α)/(1+prices_subsidy_b[3])*(R/Baseline_parameter.α*1/prices_subsidy_b[2])^(Baseline_parameter.α/(Baseline_parameter.α-1));

w_no_QS_no_subsidy = prices_no_QS_no_subsidy[2]*(1-Baseline_parameter.α)/(1+0)*(R/Baseline_parameter.α*1/prices_no_QS_no_subsidy[2])^(Baseline_parameter.α/(Baseline_parameter.α-1));
w_no_QS_subsidy_b = prices_no_QS_subsidy_b[2]*(1-Baseline_parameter.α)/(1+prices_no_QS_subsidy_b[3])*(R/Baseline_parameter.α*1/prices_no_QS_subsidy_b[2])^(Baseline_parameter.α/(Baseline_parameter.α-1));
w_inf = prices_inf[2]*(1-Baseline_parameter.α)/(1)*(R/Baseline_parameter.α*1/prices_inf[2])^(Baseline_parameter.α/(Baseline_parameter.α-1));
w_inf_sp = prices_inf_sp[2]*(1-Baseline_parameter.α)/(1+0)*(R/Baseline_parameter.α*1/prices_inf_sp[2])^(Baseline_parameter.α/(Baseline_parameter.α-1));

welfare_subsidy_partial = sum(stat_distr_no_subsidy  .* (exp.((welfare_val_subsidy_partial -welfare_val_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1; # Partial equilibrium
welfare_subsidy_nb = sum(stat_distr_no_subsidy  .* (exp.((welfare_val_subsidy_nb -welfare_val_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1; # Without balanced budget
welfare_subsidy_b = sum(stat_distr_no_subsidy  .* (exp.((welfare_val_subsidy_b -welfare_val_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
welfare_subsidy_b_nocbarQS = sum(stat_distr_no_cbarQS_no_subsidy  .* (exp.((welfare_val_no_cbarQS_subsidy_b -welfare_val_no_cbarQS_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
welfare_subsidy_b_nocbar = sum(stat_distr_no_cbar_no_subsidy  .* (exp.((welfare_val_no_cbar_subsidy_b -welfare_val_no_cbar_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
welfare_subsidy_b_noQS =sum(stat_distr_no_QS_no_subsidy  .* (exp.((welfare_val_no_QS_subsidy_b -welfare_val_no_QS_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
#welfare_subsidy_b_noQS_recab =sum(stat_distr_no_QS_recab_no_subsidy  .* (exp.((welfare_val_no_QS_recab_subsidy_b -welfare_val_no_QS_recab_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
welfare_inf = sum(stat_distr_no_subsidy .* (exp.(min.(welfare_val_inf - welfare_val_no_subsidy,10000) * (1.0 - Baseline_parameter.β)) .- 1)); 
welfare_inf_sp = sum(stat_distr_no_subsidy .* (exp.(min.(welfare_val_inf_sp - welfare_val_no_subsidy,10000) * (1.0 - Baseline_parameter.β)) .- 1)); 

V_saved_subsidy_partial_reshaped = reshape(V_saved_subsidy_partial,Baseline_parameter.ns_fine*3)
V_saved_subsidy_b_reshaped=reshape(V_saved_subsidy_b,Baseline_parameter.ns_fine*3)
V_saved_subsidy_nb_reshaped=reshape(V_saved_subsidy_nb,Baseline_parameter.ns_fine*3)
V_saved_no_cbarQS_subsidy_b_reshaped=reshape(V_saved_no_cbarQS_subsidy_b,Baseline_parameter.ns_fine*3)
V_saved_no_QS_subsidy_b_reshaped=reshape(V_saved_no_QS_subsidy_b,Baseline_parameter.ns_fine*3)
#V_saved_no_QS_recab_subsidy_b_reshaped=reshape(V_saved_no_QS_recab_subsidy_b,Baseline_parameter.ns_fine*3)
V_saved_no_cbar_subsidy_b_reshaped=reshape(V_saved_no_cbar_subsidy_b,Baseline_parameter.ns_fine*3)
V_saved_no_κ_subsidy_b_reshaped=reshape(V_saved_no_κ_subsidy_b,Baseline_parameter.ns_fine*3)
V_saved_no_F_W_subsidy_b_reshaped=reshape(V_saved_no_F_W_subsidy_b,Baseline_parameter.ns_fine*3)
V_saved_no_FM_B_subsidy_b_reshaped=reshape(V_saved_no_FM_B_subsidy_b,Baseline_parameter.ns_fine*3)

V_saved_no_subsidy_reshaped=reshape(V_saved_no_subsidy,Baseline_parameter.ns_fine*3)
V_saved_no_cbarQS_no_subsidy_reshaped=reshape(V_saved_no_cbarQS_no_subsidy,Baseline_parameter.ns_fine*3)
V_saved_no_QS_no_subsidy_reshaped=reshape(V_saved_no_QS_no_subsidy,Baseline_parameter.ns_fine*3)
#V_saved_no_QS_recab_no_subsidy_reshaped=reshape(V_saved_no_QS_recab_no_subsidy,Baseline_parameter.ns_fine*3)
V_saved_no_cbar_no_subsidy_reshaped=reshape(V_saved_no_cbar_no_subsidy,Baseline_parameter.ns_fine*3)
V_saved_no_κ_no_subsidy_reshaped=reshape(V_saved_no_κ_no_subsidy,Baseline_parameter.ns_fine*3)
V_saved_no_F_W_no_subsidy_reshaped=reshape(V_saved_no_F_W_no_subsidy,Baseline_parameter.ns_fine*3)
V_saved_no_FM_B_no_subsidy_reshaped=reshape(V_saved_no_FM_B_no_subsidy,Baseline_parameter.ns_fine*3)

V_saved_inf_reshaped=reshape(V_saved_inf,Baseline_parameter.ns_fine*3)
V_saved_inf_sp_reshaped=reshape(V_saved_inf_sp,Baseline_parameter.ns_fine*3)

welfare_subsidy_partial_real = sum(stat_distr_no_subsidy  .* (exp.((V_saved_subsidy_partial_reshaped -V_saved_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # Partial eq
welfare_subsidy_nb_real = sum(stat_distr_no_subsidy  .* (exp.((V_saved_subsidy_nb_reshaped -V_saved_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # Without balanced budget
welfare_subsidy_b_real = sum(stat_distr_no_subsidy  .* (exp.((V_saved_subsidy_b_reshaped -V_saved_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
welfare_subsidy_b_nocbarQS_real = sum(stat_distr_no_cbarQS_no_subsidy  .* (exp.((V_saved_no_cbarQS_subsidy_b_reshaped -V_saved_no_cbarQS_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
welfare_subsidy_b_nocbar_real = sum(stat_distr_no_cbar_no_subsidy  .* (exp.((V_saved_no_cbar_subsidy_b_reshaped -V_saved_no_cbar_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
#welfare_subsidy_b_noQS_real =sum(stat_distr_no_QS_no_subsidy  .* (exp.((V_saved_no_QS_subsidy_b_reshaped -V_saved_no_QS_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
#welfare_subsidy_b_noQS_recab_real =sum(stat_distr_no_QS_recab_no_subsidy  .* (exp.((V_saved_no_QS_recab_subsidy_b_reshaped -V_saved_no_QS_recab_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
welfare_subsidy_b_no_F_W_real = sum(stat_distr_no_F_W_no_subsidy  .* (exp.((V_saved_no_F_W_subsidy_b_reshaped -V_saved_no_F_W_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
welfare_subsidy_b_no_FM_B_real = sum(stat_distr_no_FM_B_no_subsidy  .* (exp.((V_saved_no_FM_B_subsidy_b_reshaped -V_saved_no_FM_B_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget

welfare_inf_real = sum(stat_distr_no_subsidy  .* (exp.((V_saved_inf_reshaped -V_saved_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
welfare_inf_sp_real = sum(stat_distr_no_subsidy  .* (exp.((V_saved_inf_sp_reshaped -V_saved_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget


welfare_subsidy_partial_real_alt = (1 - Baseline_parameter.β)* (sum(stat_distr_subsidy_partial .* V_saved_subsidy_partial_reshaped) - sum(stat_distr_no_subsidy.*V_saved_no_subsidy_reshaped)  ); # Partial equilibrium
welfare_subsidy_nb_real_alt = (1 - Baseline_parameter.β)* (sum(stat_distr_subsidy_nb .* V_saved_subsidy_nb_reshaped) - sum(stat_distr_no_subsidy.*V_saved_no_subsidy_reshaped)  ); # Without balanced budget
welfare_subsidy_b_real_alt = (1 - Baseline_parameter.β)* (sum(stat_distr_subsidy_b .* V_saved_subsidy_b_reshaped) - sum(stat_distr_no_subsidy.*V_saved_no_subsidy_reshaped)  ); # With balanced budget
welfare_subsidy_b_nocbarQS_real_alt = (1 - Baseline_parameter.β)* (sum(stat_distr_no_cbarQS_subsidy_b .* V_saved_no_cbarQS_subsidy_b_reshaped) - sum(stat_distr_no_cbarQS_no_subsidy.*V_saved_no_cbarQS_no_subsidy_reshaped)  ); # With balanced budget
welfare_subsidy_b_nocbar_real_alt = (1 - Baseline_parameter.β)* (sum(stat_distr_no_cbar_subsidy_b .* V_saved_no_cbar_subsidy_b_reshaped) - sum(stat_distr_no_cbar_no_subsidy.*V_saved_no_cbar_no_subsidy_reshaped)  );# With balanced budget
welfare_subsidy_b_noQS_real_alt =(1 - Baseline_parameter.β)* (sum(stat_distr_no_QS_subsidy_b .* V_saved_no_QS_subsidy_b_reshaped) - sum(stat_distr_no_QS_no_subsidy.*V_saved_no_QS_no_subsidy_reshaped)  ); # With balanced budget
welfare_subsidy_b_no_F_W_real_alt =(1 - Baseline_parameter.β)* (sum(stat_distr_no_F_W_subsidy_b .* V_saved_no_F_W_subsidy_b_reshaped) - sum(stat_distr_no_F_W_no_subsidy.*V_saved_no_F_W_no_subsidy_reshaped)  ); # With balanced budget
welfare_subsidy_b_no_FM_B_real_alt =(1 - Baseline_parameter.β)* (sum(stat_distr_no_FM_B_subsidy_b .* V_saved_no_FM_B_subsidy_b_reshaped) - sum(stat_distr_no_FM_B_no_subsidy.*V_saved_no_FM_B_no_subsidy_reshaped)  ); # With balanced budget



welfare_subsidy_partial_alt = (1 - Baseline_parameter.β)* (sum(stat_distr_subsidy_partial .* welfare_val_subsidy_partial) - sum(stat_distr_no_subsidy.*welfare_val_no_subsidy)  ); # Partial eq
welfare_subsidy_nb_alt = (1 - Baseline_parameter.β)* (sum(stat_distr_subsidy_nb .* welfare_val_subsidy_nb) - sum(stat_distr_no_subsidy.*welfare_val_no_subsidy)  ); # Without balanced budget
welfare_subsidy_b_alt = (1 - Baseline_parameter.β)* (sum(stat_distr_subsidy_b .* welfare_val_subsidy_b) - sum(stat_distr_no_subsidy.*welfare_val_no_subsidy)  ); # With balanced budget
welfare_subsidy_b_nocbarQS_alt = (1 - Baseline_parameter.β)* (sum(stat_distr_no_cbarQS_subsidy_b .* welfare_val_no_cbarQS_subsidy_b) - sum(stat_distr_no_cbarQS_no_subsidy.*welfare_val_no_cbarQS_no_subsidy)  ); # With balanced budget
welfare_subsidy_b_nocbar_alt = (1 - Baseline_parameter.β)* (sum(stat_distr_no_cbar_subsidy_b .* welfare_val_no_cbar_subsidy_b) - sum(stat_distr_no_cbar_no_subsidy.*welfare_val_no_cbar_no_subsidy)  );# With balanced budget
welfare_subsidy_b_noQS_alt =(1 - Baseline_parameter.β)* (sum(stat_distr_no_QS_subsidy_b .* welfare_val_no_QS_subsidy_b) - sum(stat_distr_no_QS_no_subsidy.*welfare_val_no_QS_no_subsidy)  ); # With balanced budget
welfare_subsidy_b_no_F_W_alt =(1 - Baseline_parameter.β)* (sum(stat_distr_no_F_W_subsidy_b .* welfare_val_no_F_W_subsidy_b) - sum(stat_distr_no_F_W_no_subsidy.*welfare_val_no_F_W_no_subsidy)  ); # With balanced budget
welfare_subsidy_b_no_FM_B_alt =(1 - Baseline_parameter.β)* (sum(stat_distr_no_FM_B_subsidy_b .* welfare_val_no_FM_B_subsidy_b) - sum(stat_distr_no_FM_B_no_subsidy.*welfare_val_no_FM_B_no_subsidy)  ); # With balanced budget

welfare_inf_real_alt = (1 - Baseline_parameter.β)* (sum(stat_distr_inf .* V_saved_inf_reshaped) - sum(stat_distr_no_subsidy.*V_saved_no_subsidy_reshaped)  ); 
welfare_inf_sp_real_alt = (1 - Baseline_parameter.β)* (sum(stat_distr_inf_sp .* V_saved_inf_sp_reshaped) - sum(stat_distr_no_subsidy.*V_saved_no_subsidy_reshaped)  ); 

welfare_inf_alt = (1 - Baseline_parameter.β)* (sum(stat_distr_inf .* welfare_val_inf) - sum(stat_distr_no_subsidy.*welfare_val_no_subsidy)  ); 
welfare_inf_sp_alt = (1 - Baseline_parameter.β)* (sum(stat_distr_inf_sp .* welfare_val_inf_sp) - sum(stat_distr_no_subsidy.*welfare_val_no_subsidy)  ); 

# Urban-rural welfare long run

welfare_subsidy_b_real_workers_lr = sum(stat_distr_no_subsidy[1:Baseline_parameter.ns_fine]./current_worker_pop_no_subsidy  .* (exp.((V_saved_subsidy_b_reshaped[1:Baseline_parameter.ns_fine] -V_saved_no_subsidy_reshaped[1:Baseline_parameter.ns_fine]) * (1.0 - Baseline_parameter.β) ) ))  - 1;
welfare_subsidy_b_real_farmers_lr = sum(stat_distr_no_subsidy[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]./(1.0 - current_worker_pop_no_subsidy)  .* (exp.((V_saved_subsidy_b_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)] -V_saved_no_subsidy_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]) * (1.0 - Baseline_parameter.β) ) ))  - 1;        
welfare_subsidy_b_real_alt_workers_lr = (1 - Baseline_parameter.β)* (sum(stat_distr_subsidy_b[1:Baseline_parameter.ns_fine]./current_worker_pop_subsidy_b .* V_saved_subsidy_b_reshaped[1:Baseline_parameter.ns_fine]) - sum(stat_distr_no_subsidy[1:Baseline_parameter.ns_fine]./current_worker_pop_no_subsidy .*V_saved_no_subsidy_reshaped[1:Baseline_parameter.ns_fine])  ); 
welfare_subsidy_b_real_alt_farmers_lr = (1 - Baseline_parameter.β)* (sum(stat_distr_subsidy_b[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]./(1.0 - current_worker_pop_subsidy_b) .* V_saved_subsidy_b_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]) - sum(stat_distr_no_subsidy[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]./(1.0 - current_worker_pop_no_subsidy).*V_saved_no_subsidy_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)])  ); 

welfare_subsidy_nb_real_workers_lr = sum(stat_distr_no_subsidy[1:Baseline_parameter.ns_fine]./current_worker_pop_no_subsidy  .* (exp.((V_saved_subsidy_nb_reshaped[1:Baseline_parameter.ns_fine] -V_saved_no_subsidy_reshaped[1:Baseline_parameter.ns_fine]) * (1.0 - Baseline_parameter.β) ) ))  - 1;
welfare_subsidy_nb_real_farmers_lr = sum(stat_distr_no_subsidy[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]./(1.0 - current_worker_pop_no_subsidy)  .* (exp.((V_saved_subsidy_nb_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)] -V_saved_no_subsidy_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]) * (1.0 - Baseline_parameter.β) ) ))  - 1;        
welfare_subsidy_nb_real_alt_workers_lr = (1 - Baseline_parameter.β)* (sum(stat_distr_subsidy_nb[1:Baseline_parameter.ns_fine]./current_worker_pop_subsidy_nb .* V_saved_subsidy_nb_reshaped[1:Baseline_parameter.ns_fine]) - sum(stat_distr_no_subsidy[1:Baseline_parameter.ns_fine]./current_worker_pop_no_subsidy .*V_saved_no_subsidy_reshaped[1:Baseline_parameter.ns_fine])  ); 
welfare_subsidy_nb_real_alt_farmers_lr = (1 - Baseline_parameter.β)* (sum(stat_distr_subsidy_nb[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]./(1.0 - current_worker_pop_subsidy_nb) .* V_saved_subsidy_nb_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]) - sum(stat_distr_no_subsidy[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]./(1.0 - current_worker_pop_no_subsidy).*V_saved_no_subsidy_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)])  ); 

welfare_inf_real_workers_lr = sum(stat_distr_no_subsidy[1:Baseline_parameter.ns_fine]./current_worker_pop_no_subsidy  .* (exp.((V_saved_inf_reshaped[1:Baseline_parameter.ns_fine] -V_saved_no_subsidy_reshaped[1:Baseline_parameter.ns_fine]) * (1.0 - Baseline_parameter.β) ) ))  - 1;
welfare_inf_real_farmers_lr = sum(stat_distr_no_subsidy[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]./(1.0 - current_worker_pop_no_subsidy)  .* (exp.((V_saved_inf_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)] -V_saved_no_subsidy_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]) * (1.0 - Baseline_parameter.β) ) ))  - 1;        
welfare_inf_real_alt_workers_lr = (1 - Baseline_parameter.β)* (sum(stat_distr_inf[1:Baseline_parameter.ns_fine]./current_worker_pop_inf .* V_saved_inf_reshaped[1:Baseline_parameter.ns_fine]) - sum(stat_distr_no_subsidy[1:Baseline_parameter.ns_fine]./current_worker_pop_no_subsidy .*V_saved_no_subsidy_reshaped[1:Baseline_parameter.ns_fine])  ); 
welfare_inf_real_alt_farmers_lr = (1 - Baseline_parameter.β)* (sum(stat_distr_inf[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]./(1.0 - current_worker_pop_inf) .* V_saved_inf_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]) - sum(stat_distr_no_subsidy[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]./(1.0 - current_worker_pop_no_subsidy).*V_saved_no_subsidy_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)])  ); 

welfare_inf_sp_real_workers_lr = sum(stat_distr_no_subsidy[1:Baseline_parameter.ns_fine]./current_worker_pop_no_subsidy  .* (exp.((V_saved_inf_sp_reshaped[1:Baseline_parameter.ns_fine] -V_saved_no_subsidy_reshaped[1:Baseline_parameter.ns_fine]) * (1.0 - Baseline_parameter.β) ) ))  - 1;
welfare_inf_sp_real_farmers_lr = sum(stat_distr_no_subsidy[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]./(1.0 - current_worker_pop_no_subsidy)  .* (exp.((V_saved_inf_sp_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)] -V_saved_no_subsidy_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]) * (1.0 - Baseline_parameter.β) ) ))  - 1;        
welfare_inf_sp_real_alt_workers_lr = (1 - Baseline_parameter.β)* (sum(stat_distr_inf_sp[1:Baseline_parameter.ns_fine]./current_worker_pop_inf_sp .* V_saved_inf_sp_reshaped[1:Baseline_parameter.ns_fine]) - sum(stat_distr_no_subsidy[1:Baseline_parameter.ns_fine]./current_worker_pop_no_subsidy .*V_saved_no_subsidy_reshaped[1:Baseline_parameter.ns_fine])  ); 
welfare_inf_sp_real_alt_farmers_lr = (1 - Baseline_parameter.β)* (sum(stat_distr_inf_sp[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]./(1.0 - current_worker_pop_inf_sp) .* V_saved_inf_sp_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]) - sum(stat_distr_no_subsidy[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]./(1.0 - current_worker_pop_no_subsidy).*V_saved_no_subsidy_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)])  ); 

welfare_subsidy_partial_real_workers_lr = sum(stat_distr_no_subsidy[1:Baseline_parameter.ns_fine]./current_worker_pop_no_subsidy  .* (exp.((V_saved_subsidy_partial_reshaped[1:Baseline_parameter.ns_fine] -V_saved_no_subsidy_reshaped[1:Baseline_parameter.ns_fine]) * (1.0 - Baseline_parameter.β) ) ))  - 1;
welfare_subsidy_partial_real_farmers_lr = sum(stat_distr_no_subsidy[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]./(1.0 - current_worker_pop_no_subsidy)  .* (exp.((V_saved_subsidy_partial_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)] -V_saved_no_subsidy_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]) * (1.0 - Baseline_parameter.β) ) ))  - 1;        
welfare_subsidy_partial_real_alt_workers_lr = (1 - Baseline_parameter.β)* (sum(stat_distr_subsidy_partial[1:Baseline_parameter.ns_fine]./current_worker_pop_subsidy_partial .* V_saved_subsidy_partial_reshaped[1:Baseline_parameter.ns_fine]) - sum(stat_distr_no_subsidy[1:Baseline_parameter.ns_fine]./current_worker_pop_no_subsidy .*V_saved_no_subsidy_reshaped[1:Baseline_parameter.ns_fine])  ); 
welfare_subsidy_partial_real_alt_farmers_lr = (1 - Baseline_parameter.β)* (sum(stat_distr_subsidy_partial[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]./(1.0 - current_worker_pop_subsidy_partial) .* V_saved_subsidy_partial_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]) - sum(stat_distr_no_subsidy[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]./(1.0 - current_worker_pop_no_subsidy).*V_saved_no_subsidy_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)])  ); 




welfare_trans_optimal_real_workers_lr = zeros(no_taus)
welfare_trans_optimal_real_farmers_lr = zeros(no_taus)
welfare_trans_optimal_real_alt_workers_lr = zeros(no_taus)
welfare_trans_optimal_real_alt_farmers_lr = zeros(no_taus)
for τ_index = 1:no_taus
    welfare_trans_optimal_real_workers_lr[τ_index] = sum(stat_distr_no_subsidy[1:Baseline_parameter.ns_fine]./current_worker_pop_no_subsidy  .* (exp.((V_saved_b_grid_reshaped[1:Baseline_parameter.ns_fine,τ_index] -V_saved_no_subsidy_reshaped[1:Baseline_parameter.ns_fine]) * (1.0 - Baseline_parameter.β) ) ))  - 1;
    welfare_trans_optimal_real_farmers_lr[τ_index] = sum(stat_distr_no_subsidy[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]./(1.0 - current_worker_pop_no_subsidy)  .* (exp.((V_saved_b_grid_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine),τ_index] -V_saved_no_subsidy_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]) * (1.0 - Baseline_parameter.β) ) ))  - 1;        
    welfare_trans_optimal_real_alt_workers_lr[τ_index] = (1 - Baseline_parameter.β)* (sum(stat_distr_subsidy_b_grid[1:Baseline_parameter.ns_fine,τ_index]./current_worker_pop_subsidy_b_grid[τ_index]  .* V_saved_b_grid_reshaped[1:Baseline_parameter.ns_fine,τ_index]) - sum(stat_distr_no_subsidy[1:Baseline_parameter.ns_fine]./current_worker_pop_no_subsidy .*V_saved_no_subsidy_reshaped[1:Baseline_parameter.ns_fine])  ); 
    welfare_trans_optimal_real_alt_farmers_lr[τ_index] = (1 - Baseline_parameter.β)* (sum(stat_distr_subsidy_b_grid[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine),τ_index]./(1.0 - current_worker_pop_subsidy_b_grid[τ_index]) .* V_saved_b_grid_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine),τ_index]) - sum(stat_distr_no_subsidy[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)]./(1.0 - current_worker_pop_no_subsidy).*V_saved_no_subsidy_reshaped[(Baseline_parameter.ns_fine+1):(3*Baseline_parameter.ns_fine)])  ); 
end


# Conversion to monetary values and poverty indeces

# Convert average consumption levels first, check poverty and only then convert to wealth 
#Rural avg cons: 1366 dollars table 2 in magalhaes-santaeulalia_15
P_level_naive = 0.847652638824993 # this is P_W_subsidy_b
rural_cons_in_model = P_level_naive *sum(stat_distr_subsidy_b[2561:end] .* reshape(cons_fine_local_subsidy_b[:,2:3],Baseline_parameter.ns_fine*2))/(1 - current_worker_pop_subsidy_b)

#Urban avg cons: 2912 dollars table 2 in magalhaes-santaeulalia_15
urban_cons_in_model = P_level_naive *sum(stat_distr_subsidy_b[1:2560] .* cons_fine_local_subsidy_b[:,1])/( current_worker_pop_subsidy_b)

cons_dollar_factor = (1366*0.8 + 0.2*2912)/(rural_cons_in_model*0.8 + urban_cons_in_model*0.2);

0.14581184120737242 # wages_subsidy_b
fixed_cost_converted = 0.14581184120737242* Baseline_parameter.F_W * cons_dollar_factor

income_threshold = (7.068* 365) / cons_dollar_factor;
consumption_threshold = income_threshold / (P_level_naive);
wealth_transformed_to_dollars_through_consumption = cons_dollar_factor * Baseline_parameter.s_fine[:,1]

poverty_level_subsidy_b = sum(([cons_fine_local_subsidy_b[:,1]
cons_fine_local_subsidy_b[:,2]
cons_fine_local_subsidy_b[:,3]] .< consumption_threshold) .* stat_distr_subsidy_b);
poverty_level_subsidy_nb = sum(([cons_fine_local_subsidy_nb[:,1]
cons_fine_local_subsidy_nb[:,2]
cons_fine_local_subsidy_nb[:,3]] .< consumption_threshold) .* stat_distr_subsidy_nb);
poverty_level_subsidy_partial = sum(([cons_fine_local_subsidy_partial[:,1]
cons_fine_local_subsidy_partial[:,2]
cons_fine_local_subsidy_partial[:,3]] .< consumption_threshold) .* stat_distr_subsidy_partial);
poverty_level_no_subsidy = sum(([cons_fine_local_no_subsidy[:,1]
cons_fine_local_no_subsidy[:,2]
cons_fine_local_no_subsidy[:,3]] .< consumption_threshold) .* stat_distr_no_subsidy);
poverty_level_inf = sum(([cons_fine_local_inf[:,1]
cons_fine_local_inf[:,2]
cons_fine_local_inf[:,3]] .< consumption_threshold) .* stat_distr_inf);
poverty_level_inf_sp = sum(([cons_fine_local_inf_sp[:,1]
cons_fine_local_inf_sp[:,2]
cons_fine_local_inf_sp[:,3]] .< consumption_threshold) .* stat_distr_inf_sp);

poverty_level_subsidy_b_grid = (sum(([cons_fine_local_subsidy_b_grid[:,1,:]
cons_fine_local_subsidy_b_grid[:,2,:]
cons_fine_local_subsidy_b_grid[:,3,:]] .< consumption_threshold) .* stat_distr_subsidy_b_grid, dims = 1))';

# Without subsidy but with infrastructure nb (not balance) nsp (no spillovers)
# # Construct APG changes with keeping the no_subsidy prices constant
# # First, construct the prices:
# p_B_no_subsidy,p_M_no_subsidy,R_no_subsidy,r_no_subsidy,w_no_subsidy,τ_W_no_subsidy= price_reshaper_fixed_r(prices_no_subsidy,No_subsidy_parameter.δ,No_subsidy_parameter.ϵ,No_subsidy_parameter.ψ_S,No_subsidy_parameter.ψ_B,No_subsidy_parameter.ψ_M,
# No_subsidy_parameter.p_x,No_subsidy_parameter.τ_S,No_subsidy_parameter.τ_B,No_subsidy_parameter.α,0.0,moments[1]);

# YL_manuf_subsidy_nb_cp=(p_M_no_subsidy*prod_manuf_subsidy_nb - w_no_subsidy*total_entry_cost_subsidy_nb)/(current_worker_pop_subsidy_nb*worker_pop_effective_subsidy_nb);
# YL_agr_subsidy_nb_cp=(prod_staple_subsidy_nb + p_B_no_subsidy*prod_cashcrop_subsidy_nb - Baseline_parameter.p_x*(input_staple_subsidy_nb + input_cashcrop_subsidy_nb) - (total_maintenance_cost_subsidy_nb
# )*w_no_subsidy)/(1-current_worker_pop_subsidy_nb*worker_pop_effective_subsidy_nb);
# APG_nb_cp=YL_manuf_subsidy_nb_cp/YL_agr_subsidy_nb_cp;

# YL_manuf_subsidy_b_cp=(p_M_no_subsidy*prod_manuf_subsidy_b - w_no_subsidy*total_entry_cost_subsidy_b)/(current_worker_pop_subsidy_b*worker_pop_effective_subsidy_b);
# YL_agr_subsidy_b_cp=(prod_staple_subsidy_b + p_B_no_subsidy*prod_cashcrop_subsidy_b - Baseline_parameter.p_x*(input_staple_subsidy_b + input_cashcrop_subsidy_b) - (total_maintenance_cost_subsidy_b
# )*w_no_subsidy)/(1-current_worker_pop_subsidy_b*worker_pop_effective_subsidy_b);
# APG_b_cp=YL_manuf_subsidy_b_cp/YL_agr_subsidy_b_cp;

# function lognormal(mean::Float64,variance::Float64)
#     lognormal_mean = exp.(mean + 1/2* variance);
#     lognormal_sd = (exp.(mean + 1/2* variance)* (exp(variance) - 1))^(1/2);
#     cv_var= lognormal_sd/lognormal_mean
#     return (lognormal_mean,lognormal_sd,cv_var)
# end
# (lognormal_mean_no_subsidy,lognormal_var_no_subsidy,cv_var_no_subsidy) = lognormal(MPX_mean_log_no_subsidy,var_MPX_no_subsidy)
# (lognormal_mean_subsidy_partial,lognormal_var_subsidy_partial,cv_var_subsidy_partial) =lognormal(MPX_mean_log_subsidy_partial,var_MPX_subsidy_partial)
# (lognormal_mean_subsidy_nb,lognormal_var_subsidy_nb,cv_var_subsidy_nb) =lognormal(MPX_mean_log_subsidy_nb,var_MPX_subsidy_nb)
# (lognormal_mean_subsidy_b,lognormal_var_subsidy_b,cv_var_subsidy_b) =lognormal(MPX_mean_log_subsidy_b,var_MPX_subsidy_b)
mean_land_share_staples_no_subsidy =  convert(Int64, round(100 * (mean_land_share_to_staples_among_cc_model_no_subsidy*current_cashcrop_pop_no_subsidy
 + current_staple_pop_no_subsidy)/(1 - current_worker_pop_no_subsidy ) ))
 mean_land_share_staples_subsidy_nb =  convert(Int64, round(100 * (mean_land_share_to_staples_among_cc_model_subsidy_nb*current_cashcrop_pop_subsidy_nb
 + current_staple_pop_subsidy_nb)/(1 - current_worker_pop_subsidy_nb ) ))
 mean_land_share_staples_subsidy_b =  convert(Int64, round(100 * (mean_land_share_to_staples_among_cc_model_subsidy_b*current_cashcrop_pop_subsidy_b
 + current_staple_pop_subsidy_b)/(1 - current_worker_pop_subsidy_b ) ))
 mean_land_share_staples_subsidy_partial =  convert(Int64, round(100 * (mean_land_share_to_staples_among_cc_model_subsidy_partial*current_cashcrop_pop_subsidy_partial
 + current_staple_pop_subsidy_partial)/(1 - current_worker_pop_subsidy_partial ) ))
 mean_land_share_staples_inf =  convert(Int64, round(100 * (mean_land_share_to_staples_among_cc_model_inf*current_cashcrop_pop_inf
 + current_staple_pop_inf)/(1 - current_worker_pop_inf ) ))
 mean_land_share_staples_inf_sp =  convert(Int64, round(100 * (mean_land_share_to_staples_among_cc_model_inf_sp*current_cashcrop_pop_inf_sp
 + current_staple_pop_inf_sp)/(1 - current_worker_pop_inf_sp ) ))

# Smoothing for plotting tauSgrid values for balanced budgets
aggregate_consumption_subsidy_b_grid_growth_smoothed = movmean((100*(aggregate_consumption_subsidy_b_grid./aggregate_consumption_subsidy_b_grid[no_taus]) .- 100),5)
aggregate_consumption_subsidy_b_grid_growth_smoothed = aggregate_consumption_subsidy_b_grid_growth_smoothed.-aggregate_consumption_subsidy_b_grid_growth_smoothed[no_taus];
welfare_subsidy_b_grid_real_smoothed = movmean(100*welfare_subsidy_b_grid_real,5)
welfare_subsidy_b_grid_real_smoothed = welfare_subsidy_b_grid_real_smoothed .- welfare_subsidy_b_grid_real_smoothed[no_taus];
welfare_subsidy_b_grid_real_alt_smoothed = movmean(100*welfare_subsidy_b_grid_real_alt,5)
welfare_subsidy_b_grid_real_alt_smoothed = welfare_subsidy_b_grid_real_alt_smoothed .- welfare_subsidy_b_grid_real_alt_smoothed[no_taus];

avg_agri_prod_rural_subsidy_b_grid_growth_smoothed = movmean((100*avg_agri_prod_rural_subsidy_b_grid./(avg_agri_prod_rural_subsidy_b_grid[no_taus]).-100),5)
avg_agri_prod_rural_subsidy_b_grid_growth_smoothed = avg_agri_prod_rural_subsidy_b_grid_growth_smoothed .- avg_agri_prod_rural_subsidy_b_grid_growth_smoothed[no_taus]
avg_labor_prod_urban_subsidy_b_grid_growth_smoothed = movmean((100*avg_labor_prod_urban_subsidy_b_grid./(avg_labor_prod_urban_subsidy_b_grid[no_taus]).-100),5)
avg_labor_prod_urban_subsidy_b_grid_growth_smoothed = avg_labor_prod_urban_subsidy_b_grid_growth_smoothed .- avg_labor_prod_urban_subsidy_b_grid_growth_smoothed[no_taus]
current_worker_pop_subsidy_b_grid_growth_smoothed = movmean((100*current_worker_pop_subsidy_b_grid./(current_worker_pop_subsidy_b_grid[no_taus]).-100),5)
current_worker_pop_subsidy_b_grid_growth_smoothed = current_worker_pop_subsidy_b_grid_growth_smoothed .- current_worker_pop_subsidy_b_grid_growth_smoothed[no_taus]

staple_productivity_subsidy_b_grid_smoothed = movmean((100*(staple_productivity_subsidy_b_grid./staple_productivity_subsidy_b_grid[no_taus]) .- 100),5)
staple_productivity_subsidy_b_grid_smoothed = staple_productivity_subsidy_b_grid_smoothed .- staple_productivity_subsidy_b_grid_smoothed[no_taus]

cashcrop_productivity_subsidy_b_grid_smoothed = movmean((100*(cashcrop_productivity_subsidy_b_grid./cashcrop_productivity_subsidy_b_grid[no_taus]) .- 100),5)
cashcrop_productivity_subsidy_b_grid_smoothed = cashcrop_productivity_subsidy_b_grid_smoothed .- cashcrop_productivity_subsidy_b_grid_smoothed[no_taus]

prod_staple_subsidy_b_grid_smoothed = movmean((100*(prod_staple_subsidy_b_grid./prod_staple_subsidy_b_grid[no_taus]) .- 100),5)
prod_staple_subsidy_b_grid_smoothed = prod_staple_subsidy_b_grid_smoothed .- prod_staple_subsidy_b_grid_smoothed[no_taus]

prod_cashcrop_subsidy_b_grid_smoothed = movmean((100*(prod_cashcrop_subsidy_b_grid./prod_cashcrop_subsidy_b_grid[no_taus]) .- 100),5)
prod_cashcrop_subsidy_b_grid_smoothed = prod_cashcrop_subsidy_b_grid_smoothed .- prod_cashcrop_subsidy_b_grid_smoothed[no_taus]


savings_subsidy_b_grid = zeros(no_taus);
for iterate = 1:no_taus
    savings_subsidy_b_grid[iterate] = sum(stat_distr_subsidy_b_grid[:, iterate] .* reshape(a_prime_fine_local_subsidy_b_grid[:,:, iterate],Baseline_parameter.ns_fine*3))
end

savings_subsidy_b_grid_smoothed = movmean((100*(savings_subsidy_b_grid./savings_subsidy_b_grid[no_taus]) .- 100),5)
savings_subsidy_b_grid_smoothed = savings_subsidy_b_grid_smoothed .- savings_subsidy_b_grid_smoothed[no_taus]

var_MPX_subsidy_b_grid_growth_smoothed = movmean((100*(var_MPX_subsidy_b_grid./var_MPX_subsidy_b_grid[no_taus]) .- 100),5)
var_MPX_subsidy_b_grid_growth_smoothed = var_MPX_subsidy_b_grid_growth_smoothed .- var_MPX_subsidy_b_grid_growth_smoothed[no_taus]
#a_prime_difference_subsidy_b_grid_growth_smoothed = movmean((100*a_prime_difference_subsidy_b_grid ),5)


c_S_worker_sum_subsidy_b_grid_growth_smoothed = movmean((100*(c_S_worker_sum_subsidy_b_grid./c_S_worker_sum_subsidy_b_grid[no_taus]) .- 100),5)
c_S_worker_sum_subsidy_b_grid_growth_smoothed = c_S_worker_sum_subsidy_b_grid_growth_smoothed .- c_S_worker_sum_subsidy_b_grid_growth_smoothed[no_taus]
c_S_staple_sum_subsidy_b_grid_growth_smoothed = movmean((100*(c_S_staple_sum_subsidy_b_grid./c_S_staple_sum_subsidy_b_grid[no_taus]) .- 100),5)
c_S_staple_sum_subsidy_b_grid_growth_smoothed = c_S_staple_sum_subsidy_b_grid_growth_smoothed .- c_S_staple_sum_subsidy_b_grid_growth_smoothed[no_taus]
c_S_cashcrop_sum_subsidy_b_grid_growth_smoothed = movmean((100*(c_S_cashcrop_sum_subsidy_b_grid./c_S_cashcrop_sum_subsidy_b_grid[no_taus]) .- 100),5)
c_S_cashcrop_sum_subsidy_b_grid_growth_smoothed = c_S_cashcrop_sum_subsidy_b_grid_growth_smoothed .- c_S_cashcrop_sum_subsidy_b_grid_growth_smoothed[no_taus]
c_M_worker_sum_subsidy_b_grid_growth_smoothed = movmean((100*(c_M_worker_sum_subsidy_b_grid./c_M_worker_sum_subsidy_b_grid[no_taus]) .- 100),5)
c_M_worker_sum_subsidy_b_grid_growth_smoothed = c_M_worker_sum_subsidy_b_grid_growth_smoothed .- c_M_worker_sum_subsidy_b_grid_growth_smoothed[no_taus]
c_M_staple_sum_subsidy_b_grid_growth_smoothed = movmean((100*(c_M_staple_sum_subsidy_b_grid./c_M_staple_sum_subsidy_b_grid[no_taus]) .- 100),5)
c_M_staple_sum_subsidy_b_grid_growth_smoothed = c_M_staple_sum_subsidy_b_grid_growth_smoothed .- c_M_staple_sum_subsidy_b_grid_growth_smoothed[no_taus]
c_M_cashcrop_sum_subsidy_b_grid_growth_smoothed = movmean((100*(c_M_cashcrop_sum_subsidy_b_grid./c_M_cashcrop_sum_subsidy_b_grid[no_taus]) .- 100),5)
c_M_cashcrop_sum_subsidy_b_grid_growth_smoothed = c_M_cashcrop_sum_subsidy_b_grid_growth_smoothed .- c_M_cashcrop_sum_subsidy_b_grid_growth_smoothed[no_taus]
c_B_worker_sum_subsidy_b_grid_growth_smoothed = movmean((100*(c_B_worker_sum_subsidy_b_grid./c_B_worker_sum_subsidy_b_grid[no_taus]) .- 100),5)
c_B_worker_sum_subsidy_b_grid_growth_smoothed = c_B_worker_sum_subsidy_b_grid_growth_smoothed .- c_B_worker_sum_subsidy_b_grid_growth_smoothed[no_taus]
c_B_staple_sum_subsidy_b_grid_growth_smoothed = movmean((100*(c_B_staple_sum_subsidy_b_grid./c_B_staple_sum_subsidy_b_grid[no_taus]) .- 100),5)
c_B_staple_sum_subsidy_b_grid_growth_smoothed = c_B_staple_sum_subsidy_b_grid_growth_smoothed .- c_B_staple_sum_subsidy_b_grid_growth_smoothed[no_taus]
c_B_cashcrop_sum_subsidy_b_grid_growth_smoothed = movmean((100*(c_B_cashcrop_sum_subsidy_b_grid./c_B_cashcrop_sum_subsidy_b_grid[no_taus]) .- 100),5)
c_B_cashcrop_sum_subsidy_b_grid_growth_smoothed = c_B_cashcrop_sum_subsidy_b_grid_growth_smoothed .- c_B_cashcrop_sum_subsidy_b_grid_growth_smoothed[no_taus]

APG_subsidy_b_grid_growth_smoothed = movmean((100*(APG_subsidy_b_grid./APG_subsidy_b_grid[no_taus]) .- 100),5)
APG_subsidy_b_grid_growth_smoothed = APG_subsidy_b_grid_growth_smoothed .- APG_subsidy_b_grid_growth_smoothed[no_taus]

fertilizer_use_subsidy_b_grid = Import_value_subsidy_b_grid/Baseline_parameter.p_x;
fertilizer_use_subsidy_b_grid_growth_smoothed = movmean((100*(fertilizer_use_subsidy_b_grid./fertilizer_use_subsidy_b_grid[no_taus]) .- 100),5)
fertilizer_use_subsidy_b_grid_growth_smoothed = fertilizer_use_subsidy_b_grid_growth_smoothed.-fertilizer_use_subsidy_b_grid_growth_smoothed[no_taus];

cashcrop_productivity_value_subsidy_b_grid_smoothed = movmean((100*(prices_subsidy_b_grid[1,:].*cashcrop_productivity_subsidy_b_grid./prices_subsidy_b_grid[1,no_taus]./cashcrop_productivity_subsidy_b_grid[no_taus]) .- 100),5)
cashcrop_productivity_value_subsidy_b_grid_smoothed = cashcrop_productivity_value_subsidy_b_grid_smoothed .- cashcrop_productivity_value_subsidy_b_grid_smoothed[no_taus]

# Cashcrop price
cashcrop_price_subsidy_b_grid_smoothed = movmean((100*(prices_subsidy_b_grid[1,:]./prices_subsidy_b_grid[1,no_taus]) .- 100),5)
cashcrop_price_subsidy_b_grid_smoothed = cashcrop_price_subsidy_b_grid_smoothed .- cashcrop_price_subsidy_b_grid_smoothed[no_taus]

# Relative land
relative_land_to_staples_subsidy_b_grid_smoothed = movmean((100*(relative_land_to_staples_subsidy_b_grid./relative_land_to_staples_subsidy_b_grid[no_taus]) .- 100),5)
relative_land_to_staples_subsidy_b_grid_smoothed = relative_land_to_staples_subsidy_b_grid_smoothed .- relative_land_to_staples_subsidy_b_grid_smoothed[no_taus]

# Smoothing for plotting tauSgrid values without balanced budgets
aggregate_consumption_subsidy_nb_grid_growth_smoothed = movmean((100*(aggregate_consumption_subsidy_nb_grid./aggregate_consumption_subsidy_nb_grid[no_taus]) .- 100),5)
aggregate_consumption_subsidy_nb_grid_growth_smoothed = aggregate_consumption_subsidy_nb_grid_growth_smoothed.-aggregate_consumption_subsidy_nb_grid_growth_smoothed[no_taus];
welfare_subsidy_nb_grid_real_smoothed = movmean(100*welfare_subsidy_nb_grid_real,5)
welfare_subsidy_nb_grid_real_smoothed = welfare_subsidy_nb_grid_real_smoothed .- welfare_subsidy_nb_grid_real_smoothed[no_taus];
welfare_subsidy_nb_grid_real_alt_smoothed = movmean(100*welfare_subsidy_nb_grid_real_alt,5)
welfare_subsidy_nb_grid_real_alt_smoothed = welfare_subsidy_nb_grid_real_alt_smoothed .- welfare_subsidy_nb_grid_real_alt_smoothed[no_taus];

avg_agri_prod_rural_subsidy_nb_grid_growth_smoothed = movmean((100*avg_agri_prod_rural_subsidy_nb_grid./(avg_agri_prod_rural_subsidy_nb_grid[no_taus]).-100),5)
avg_agri_prod_rural_subsidy_nb_grid_growth_smoothed = avg_agri_prod_rural_subsidy_nb_grid_growth_smoothed .- avg_agri_prod_rural_subsidy_nb_grid_growth_smoothed[no_taus]
avg_labor_prod_urban_subsidy_nb_grid_growth_smoothed = movmean((100*avg_labor_prod_urban_subsidy_nb_grid./(avg_labor_prod_urban_subsidy_nb_grid[no_taus]).-100),5)
avg_labor_prod_urban_subsidy_nb_grid_growth_smoothed = avg_labor_prod_urban_subsidy_nb_grid_growth_smoothed .- avg_labor_prod_urban_subsidy_nb_grid_growth_smoothed[no_taus]
current_worker_pop_subsidy_nb_grid_growth_smoothed = movmean((100*current_worker_pop_subsidy_nb_grid./(current_worker_pop_subsidy_nb_grid[no_taus]).-100),5)
current_worker_pop_subsidy_nb_grid_growth_smoothed = current_worker_pop_subsidy_nb_grid_growth_smoothed .- current_worker_pop_subsidy_nb_grid_growth_smoothed[no_taus]

staple_productivity_subsidy_nb_grid_smoothed = movmean((100*(staple_productivity_subsidy_nb_grid./staple_productivity_subsidy_nb_grid[no_taus]) .- 100),5)
staple_productivity_subsidy_nb_grid_smoothed = staple_productivity_subsidy_nb_grid_smoothed .- staple_productivity_subsidy_nb_grid_smoothed[no_taus]

cashcrop_productivity_subsidy_nb_grid_smoothed = movmean((100*(cashcrop_productivity_subsidy_nb_grid./cashcrop_productivity_subsidy_nb_grid[no_taus]) .- 100),5)
cashcrop_productivity_subsidy_nb_grid_smoothed = cashcrop_productivity_subsidy_nb_grid_smoothed .- cashcrop_productivity_subsidy_nb_grid_smoothed[no_taus]

prod_staple_subsidy_nb_grid_smoothed = movmean((100*(prod_staple_subsidy_nb_grid./prod_staple_subsidy_nb_grid[no_taus]) .- 100),5)
prod_staple_subsidy_nb_grid_smoothed = prod_staple_subsidy_nb_grid_smoothed .- prod_staple_subsidy_nb_grid_smoothed[no_taus]

prod_cashcrop_subsidy_nb_grid_smoothed = movmean((100*(prod_cashcrop_subsidy_nb_grid./prod_cashcrop_subsidy_nb_grid[no_taus]) .- 100),5)
prod_cashcrop_subsidy_nb_grid_smoothed = prod_cashcrop_subsidy_nb_grid_smoothed .- prod_cashcrop_subsidy_nb_grid_smoothed[no_taus]


savings_subsidy_nb_grid = zeros(no_taus);
for iterate = 1:no_taus
    savings_subsidy_nb_grid[iterate] = sum(stat_distr_subsidy_nb_grid[:, iterate] .* reshape(a_prime_fine_local_subsidy_nb_grid[:,:, iterate],Baseline_parameter.ns_fine*3))
end

savings_subsidy_nb_grid_smoothed = movmean((100*(savings_subsidy_nb_grid./savings_subsidy_nb_grid[no_taus]) .- 100),5)
savings_subsidy_nb_grid_smoothed = savings_subsidy_nb_grid_smoothed .- savings_subsidy_nb_grid_smoothed[no_taus]

var_MPX_subsidy_nb_grid_growth_smoothed = movmean((100*(var_MPX_subsidy_nb_grid./var_MPX_subsidy_nb_grid[no_taus]) .- 100),5)
var_MPX_subsidy_nb_grid_growth_smoothed = var_MPX_subsidy_nb_grid_growth_smoothed .- var_MPX_subsidy_nb_grid_growth_smoothed[no_taus]
#a_prime_difference_subsidy_nb_grid_growth_smoothed = movmean((100*a_prime_difference_subsidy_nb_grid ),5)


c_S_worker_sum_subsidy_nb_grid_growth_smoothed = movmean((100*(c_S_worker_sum_subsidy_nb_grid./c_S_worker_sum_subsidy_nb_grid[no_taus]) .- 100),5)
c_S_worker_sum_subsidy_nb_grid_growth_smoothed = c_S_worker_sum_subsidy_nb_grid_growth_smoothed .- c_S_worker_sum_subsidy_nb_grid_growth_smoothed[no_taus]
c_S_staple_sum_subsidy_nb_grid_growth_smoothed = movmean((100*(c_S_staple_sum_subsidy_nb_grid./c_S_staple_sum_subsidy_nb_grid[no_taus]) .- 100),5)
c_S_staple_sum_subsidy_nb_grid_growth_smoothed = c_S_staple_sum_subsidy_nb_grid_growth_smoothed .- c_S_staple_sum_subsidy_nb_grid_growth_smoothed[no_taus]
c_S_cashcrop_sum_subsidy_nb_grid_growth_smoothed = movmean((100*(c_S_cashcrop_sum_subsidy_nb_grid./c_S_cashcrop_sum_subsidy_nb_grid[no_taus]) .- 100),5)
c_S_cashcrop_sum_subsidy_nb_grid_growth_smoothed = c_S_cashcrop_sum_subsidy_nb_grid_growth_smoothed .- c_S_cashcrop_sum_subsidy_nb_grid_growth_smoothed[no_taus]
c_M_worker_sum_subsidy_nb_grid_growth_smoothed = movmean((100*(c_M_worker_sum_subsidy_nb_grid./c_M_worker_sum_subsidy_nb_grid[no_taus]) .- 100),5)
c_M_worker_sum_subsidy_nb_grid_growth_smoothed = c_M_worker_sum_subsidy_nb_grid_growth_smoothed .- c_M_worker_sum_subsidy_nb_grid_growth_smoothed[no_taus]
c_M_staple_sum_subsidy_nb_grid_growth_smoothed = movmean((100*(c_M_staple_sum_subsidy_nb_grid./c_M_staple_sum_subsidy_nb_grid[no_taus]) .- 100),5)
c_M_staple_sum_subsidy_nb_grid_growth_smoothed = c_M_staple_sum_subsidy_nb_grid_growth_smoothed .- c_M_staple_sum_subsidy_nb_grid_growth_smoothed[no_taus]
c_M_cashcrop_sum_subsidy_nb_grid_growth_smoothed = movmean((100*(c_M_cashcrop_sum_subsidy_nb_grid./c_M_cashcrop_sum_subsidy_nb_grid[no_taus]) .- 100),5)
c_M_cashcrop_sum_subsidy_nb_grid_growth_smoothed = c_M_cashcrop_sum_subsidy_nb_grid_growth_smoothed .- c_M_cashcrop_sum_subsidy_nb_grid_growth_smoothed[no_taus]
c_B_worker_sum_subsidy_nb_grid_growth_smoothed = movmean((100*(c_B_worker_sum_subsidy_nb_grid./c_B_worker_sum_subsidy_nb_grid[no_taus]) .- 100),5)
c_B_worker_sum_subsidy_nb_grid_growth_smoothed = c_B_worker_sum_subsidy_nb_grid_growth_smoothed .- c_B_worker_sum_subsidy_nb_grid_growth_smoothed[no_taus]
c_B_staple_sum_subsidy_nb_grid_growth_smoothed = movmean((100*(c_B_staple_sum_subsidy_nb_grid./c_B_staple_sum_subsidy_nb_grid[no_taus]) .- 100),5)
c_B_staple_sum_subsidy_nb_grid_growth_smoothed = c_B_staple_sum_subsidy_nb_grid_growth_smoothed .- c_B_staple_sum_subsidy_nb_grid_growth_smoothed[no_taus]
c_B_cashcrop_sum_subsidy_nb_grid_growth_smoothed = movmean((100*(c_B_cashcrop_sum_subsidy_nb_grid./c_B_cashcrop_sum_subsidy_nb_grid[no_taus]) .- 100),5)
c_B_cashcrop_sum_subsidy_nb_grid_growth_smoothed = c_B_cashcrop_sum_subsidy_nb_grid_growth_smoothed .- c_B_cashcrop_sum_subsidy_nb_grid_growth_smoothed[no_taus]

APG_subsidy_nb_grid_growth_smoothed = movmean((100*(APG_subsidy_nb_grid./APG_subsidy_nb_grid[no_taus]) .- 100),5)
APG_subsidy_nb_grid_growth_smoothed = APG_subsidy_nb_grid_growth_smoothed .- APG_subsidy_nb_grid_growth_smoothed[no_taus]

fertilizer_use_subsidy_nb_grid = Import_value_subsidy_nb_grid/Baseline_parameter.p_x;
fertilizer_use_subsidy_nb_grid_growth_smoothed = movmean((100*(fertilizer_use_subsidy_nb_grid./fertilizer_use_subsidy_nb_grid[no_taus]) .- 100),5)
fertilizer_use_subsidy_nb_grid_growth_smoothed = fertilizer_use_subsidy_nb_grid_growth_smoothed.-fertilizer_use_subsidy_nb_grid_growth_smoothed[no_taus];

cashcrop_productivity_value_subsidy_nb_grid_smoothed = movmean((100*(prices_subsidy_nb_grid[1,:].*cashcrop_productivity_subsidy_nb_grid./prices_subsidy_nb_grid[1,no_taus]./cashcrop_productivity_subsidy_nb_grid[no_taus]) .- 100),5)
cashcrop_productivity_value_subsidy_nb_grid_smoothed = cashcrop_productivity_value_subsidy_nb_grid_smoothed .- cashcrop_productivity_value_subsidy_nb_grid_smoothed[no_taus]

# Cashcrop price
cashcrop_price_subsidy_nb_grid_smoothed = movmean((100*(prices_subsidy_nb_grid[1,:]./prices_subsidy_nb_grid[1,no_taus]) .- 100),5)
cashcrop_price_subsidy_nb_grid_smoothed = cashcrop_price_subsidy_nb_grid_smoothed .- cashcrop_price_subsidy_nb_grid_smoothed[no_taus]

# Relative land
relative_land_to_staples_subsidy_nb_grid_smoothed = movmean((100*(relative_land_to_staples_subsidy_nb_grid./relative_land_to_staples_subsidy_nb_grid[no_taus]) .- 100),5)
relative_land_to_staples_subsidy_nb_grid_smoothed = relative_land_to_staples_subsidy_nb_grid_smoothed .- relative_land_to_staples_subsidy_nb_grid_smoothed[no_taus]


# Undernourishment requires additional calculations
# Staple consumption w a threshold value -
# assume the threshold is chosen such that the calibration hits the average between the 17.9% (Pauw, Beck, and Mussa (2014) )
# and 24.5 (Malawi NSO) extreme food poverty rate of HHs. for now, its 20% in 2010 in Malawi
cons_level_substinence = 0.109;
    # Distribution of past occupations
worker_past_dist_subsidy_b = stat_distr_subsidy_b[(Baseline_parameter.ns_fine *0 + 1):(Baseline_parameter.ns_fine *1)];
staple_past_dist_subsidy_b = stat_distr_subsidy_b[(Baseline_parameter.ns_fine *1 + 1):(Baseline_parameter.ns_fine *2)];
cash_crop_past_dist_subsidy_b = stat_distr_subsidy_b[(Baseline_parameter.ns_fine *2 + 1):(Baseline_parameter.ns_fine *3)];

stay_workers_subsidy_b = worker_past_dist_subsidy_b.*(future_occupation_fine_local_subsidy_b[:,1].==1);
exit_staple_to_work_subsidy_b = staple_past_dist_subsidy_b.*(future_occupation_fine_local_subsidy_b[:,2].==1);
exit_cashcrop_to_work_subsidy_b = cash_crop_past_dist_subsidy_b.*(future_occupation_fine_local_subsidy_b[:,3].==1);
current_workers_subsidy_b = stay_workers_subsidy_b + exit_staple_to_work_subsidy_b + exit_cashcrop_to_work_subsidy_b;

entrants_staple_from_workers_subsidy_b = worker_past_dist_subsidy_b.*(future_occupation_fine_local_subsidy_b[:,1].==2);
incumbents_staple_subsidy_b = staple_past_dist_subsidy_b.*(future_occupation_fine_local_subsidy_b[:,2].==2);
exit_cashcrop_to_staple_subsidy_b = cash_crop_past_dist_subsidy_b.*(future_occupation_fine_local_subsidy_b[:,3].==2);
current_staple_subsidy_b = entrants_staple_from_workers_subsidy_b + incumbents_staple_subsidy_b + exit_cashcrop_to_staple_subsidy_b;

entrants_cashcrop_from_workers_subsidy_b = worker_past_dist_subsidy_b.*(future_occupation_fine_local_subsidy_b[:,1].==3);
entrants_from_staple_to_cashcrop_subsidy_b= staple_past_dist_subsidy_b.*(future_occupation_fine_local_subsidy_b[:,2].==3);
incumbents_cashcrop_subsidy_b = cash_crop_past_dist_subsidy_b.*(future_occupation_fine_local_subsidy_b[:,3].==3);
current_cashcrop_subsidy_b = entrants_cashcrop_from_workers_subsidy_b + incumbents_cashcrop_subsidy_b + entrants_from_staple_to_cashcrop_subsidy_b;

undernutritioned_workers_subsidy_b = sum( (c_S_W_fine_subsidy_b[:,1].<cons_level_substinence).*stay_workers_subsidy_b
 + (c_S_W_fine_subsidy_b[:,2].<cons_level_substinence).*exit_staple_to_work_subsidy_b +
(c_S_W_fine_subsidy_b[:,3].<cons_level_substinence) .*exit_cashcrop_to_work_subsidy_b );
undernutritioned_staple_farmer_subsidy_b = sum( (c_S_S_fine_subsidy_b[:,1].<cons_level_substinence).*entrants_staple_from_workers_subsidy_b +
(c_S_S_fine_subsidy_b[:,2].<cons_level_substinence).*incumbents_staple_subsidy_b +
(c_S_S_fine_subsidy_b[:,3].<cons_level_substinence) .*exit_cashcrop_to_staple_subsidy_b );
undernutritioned_cashcrop_farmer_subsidy_b = sum( (c_S_B_fine_subsidy_b[:,1].<cons_level_substinence).*entrants_cashcrop_from_workers_subsidy_b
 + (c_S_B_fine_subsidy_b[:,2].<cons_level_substinence).*entrants_from_staple_to_cashcrop_subsidy_b +
(c_S_B_fine_subsidy_b[:,3].<cons_level_substinence) .*incumbents_cashcrop_subsidy_b );
undernourished_subsidy_b = undernutritioned_workers_subsidy_b + undernutritioned_staple_farmer_subsidy_b +undernutritioned_cashcrop_farmer_subsidy_b

# No subsidy undernourished
worker_past_dist_no_subsidy = stat_distr_no_subsidy[(Baseline_parameter.ns_fine *0 + 1):(Baseline_parameter.ns_fine *1)];
staple_past_dist_no_subsidy = stat_distr_no_subsidy[(Baseline_parameter.ns_fine *1 + 1):(Baseline_parameter.ns_fine *2)];
cash_crop_past_dist_no_subsidy = stat_distr_no_subsidy[(Baseline_parameter.ns_fine *2 + 1):(Baseline_parameter.ns_fine *3)];

stay_workers_no_subsidy = worker_past_dist_no_subsidy.*(future_occupation_fine_local_no_subsidy[:,1].==1);
exit_staple_to_work_no_subsidy = staple_past_dist_no_subsidy.*(future_occupation_fine_local_no_subsidy[:,2].==1);
exit_cashcrop_to_work_no_subsidy = cash_crop_past_dist_no_subsidy.*(future_occupation_fine_local_no_subsidy[:,3].==1);
current_workers_no_subsidy = stay_workers_no_subsidy + exit_staple_to_work_no_subsidy + exit_cashcrop_to_work_no_subsidy;

entrants_staple_from_workers_no_subsidy = worker_past_dist_no_subsidy.*(future_occupation_fine_local_no_subsidy[:,1].==2);
incumbents_staple_no_subsidy = staple_past_dist_no_subsidy.*(future_occupation_fine_local_no_subsidy[:,2].==2);
exit_cashcrop_to_staple_no_subsidy = cash_crop_past_dist_no_subsidy.*(future_occupation_fine_local_no_subsidy[:,3].==2);
current_staple_no_subsidy = entrants_staple_from_workers_no_subsidy + incumbents_staple_no_subsidy + exit_cashcrop_to_staple_no_subsidy;

entrants_cashcrop_from_workers_no_subsidy = worker_past_dist_no_subsidy.*(future_occupation_fine_local_no_subsidy[:,1].==3);
entrants_from_staple_to_cashcrop_no_subsidy= staple_past_dist_no_subsidy.*(future_occupation_fine_local_no_subsidy[:,2].==3);
incumbents_cashcrop_no_subsidy = cash_crop_past_dist_no_subsidy.*(future_occupation_fine_local_no_subsidy[:,3].==3);
current_cashcrop_no_subsidy = entrants_cashcrop_from_workers_no_subsidy + incumbents_cashcrop_no_subsidy + entrants_from_staple_to_cashcrop_no_subsidy;

undernutritioned_workers_no_subsidy = sum( (c_S_W_fine_no_subsidy[:,1].<cons_level_substinence).*stay_workers_no_subsidy
    + (c_S_W_fine_no_subsidy[:,2].<cons_level_substinence).*exit_staple_to_work_no_subsidy +
(c_S_W_fine_no_subsidy[:,3].<cons_level_substinence) .*exit_cashcrop_to_work_no_subsidy );
undernutritioned_staple_farmer_no_subsidy = sum( (c_S_S_fine_no_subsidy[:,1].<cons_level_substinence).*entrants_staple_from_workers_no_subsidy +
(c_S_S_fine_no_subsidy[:,2].<cons_level_substinence).*incumbents_staple_no_subsidy +
(c_S_S_fine_no_subsidy[:,3].<cons_level_substinence) .*exit_cashcrop_to_staple_no_subsidy );
undernutritioned_cashcrop_farmer_no_subsidy = sum( (c_S_B_fine_no_subsidy[:,1].<cons_level_substinence).*entrants_cashcrop_from_workers_no_subsidy
    + (c_S_B_fine_no_subsidy[:,2].<cons_level_substinence).*entrants_from_staple_to_cashcrop_no_subsidy +
(c_S_B_fine_no_subsidy[:,3].<cons_level_substinence) .*incumbents_cashcrop_no_subsidy );

undernourished_no_subsidy = undernutritioned_workers_no_subsidy + undernutritioned_staple_farmer_no_subsidy +undernutritioned_cashcrop_farmer_no_subsidy;

# PE undernourished
worker_past_dist_subsidy_partial = stat_distr_subsidy_partial[(Baseline_parameter.ns_fine *0 + 1):(Baseline_parameter.ns_fine *1)];
staple_past_dist_subsidy_partial = stat_distr_subsidy_partial[(Baseline_parameter.ns_fine *1 + 1):(Baseline_parameter.ns_fine *2)];
cash_crop_past_dist_subsidy_partial = stat_distr_subsidy_partial[(Baseline_parameter.ns_fine *2 + 1):(Baseline_parameter.ns_fine *3)];

stay_workers_subsidy_partial = worker_past_dist_subsidy_partial.*(future_occupation_fine_local_subsidy_partial[:,1].==1);
exit_staple_to_work_subsidy_partial = staple_past_dist_subsidy_partial.*(future_occupation_fine_local_subsidy_partial[:,2].==1);
exit_cashcrop_to_work_subsidy_partial = cash_crop_past_dist_subsidy_partial.*(future_occupation_fine_local_subsidy_partial[:,3].==1);
current_workers_subsidy_partial = stay_workers_subsidy_partial + exit_staple_to_work_subsidy_partial + exit_cashcrop_to_work_subsidy_partial;

entrants_staple_from_workers_subsidy_partial = worker_past_dist_subsidy_partial.*(future_occupation_fine_local_subsidy_partial[:,1].==2);
incumbents_staple_subsidy_partial = staple_past_dist_subsidy_partial.*(future_occupation_fine_local_subsidy_partial[:,2].==2);
exit_cashcrop_to_staple_subsidy_partial = cash_crop_past_dist_subsidy_partial.*(future_occupation_fine_local_subsidy_partial[:,3].==2);
current_staple_subsidy_partial = entrants_staple_from_workers_subsidy_partial + incumbents_staple_subsidy_partial + exit_cashcrop_to_staple_subsidy_partial;

entrants_cashcrop_from_workers_subsidy_partial = worker_past_dist_subsidy_partial.*(future_occupation_fine_local_subsidy_partial[:,1].==3);
entrants_from_staple_to_cashcrop_subsidy_partial= staple_past_dist_subsidy_partial.*(future_occupation_fine_local_subsidy_partial[:,2].==3);
incumbents_cashcrop_subsidy_partial = cash_crop_past_dist_subsidy_partial.*(future_occupation_fine_local_subsidy_partial[:,3].==3);
current_cashcrop_subsidy_partial = entrants_cashcrop_from_workers_subsidy_partial + incumbents_cashcrop_subsidy_partial + entrants_from_staple_to_cashcrop_subsidy_partial;

undernutritioned_workers_subsidy_partial = sum( (c_S_W_fine_subsidy_partial[:,1].<cons_level_substinence).*stay_workers_subsidy_partial
    + (c_S_W_fine_subsidy_partial[:,2].<cons_level_substinence).*exit_staple_to_work_subsidy_partial +
(c_S_W_fine_subsidy_partial[:,3].<cons_level_substinence) .*exit_cashcrop_to_work_subsidy_partial );
undernutritioned_staple_farmer_subsidy_partial = sum( (c_S_S_fine_subsidy_partial[:,1].<cons_level_substinence).*entrants_staple_from_workers_subsidy_partial +
(c_S_S_fine_subsidy_partial[:,2].<cons_level_substinence).*incumbents_staple_subsidy_partial +
(c_S_S_fine_subsidy_partial[:,3].<cons_level_substinence) .*exit_cashcrop_to_staple_subsidy_partial );
undernutritioned_cashcrop_farmer_subsidy_partial = sum( (c_S_B_fine_subsidy_partial[:,1].<cons_level_substinence).*entrants_cashcrop_from_workers_subsidy_partial
    + (c_S_B_fine_subsidy_partial[:,2].<cons_level_substinence).*entrants_from_staple_to_cashcrop_subsidy_partial +
(c_S_B_fine_subsidy_partial[:,3].<cons_level_substinence) .*incumbents_cashcrop_subsidy_partial );

undernourished_subsidy_partial = undernutritioned_workers_subsidy_partial + undernutritioned_staple_farmer_subsidy_partial +undernutritioned_cashcrop_farmer_subsidy_partial;

# aid-financed subsidy undernourished
worker_past_dist_subsidy_nb = stat_distr_subsidy_nb[(Baseline_parameter.ns_fine *0 + 1):(Baseline_parameter.ns_fine *1)];
staple_past_dist_subsidy_nb = stat_distr_subsidy_nb[(Baseline_parameter.ns_fine *1 + 1):(Baseline_parameter.ns_fine *2)];
cash_crop_past_dist_subsidy_nb = stat_distr_subsidy_nb[(Baseline_parameter.ns_fine *2 + 1):(Baseline_parameter.ns_fine *3)];

stay_workers_subsidy_nb = worker_past_dist_subsidy_nb.*(future_occupation_fine_local_subsidy_nb[:,1].==1);
exit_staple_to_work_subsidy_nb = staple_past_dist_subsidy_nb.*(future_occupation_fine_local_subsidy_nb[:,2].==1);
exit_cashcrop_to_work_subsidy_nb = cash_crop_past_dist_subsidy_nb.*(future_occupation_fine_local_subsidy_nb[:,3].==1);
current_workers_subsidy_nb = stay_workers_subsidy_nb + exit_staple_to_work_subsidy_nb + exit_cashcrop_to_work_subsidy_nb;

entrants_staple_from_workers_subsidy_nb = worker_past_dist_subsidy_nb.*(future_occupation_fine_local_subsidy_nb[:,1].==2);
incumbents_staple_subsidy_nb = staple_past_dist_subsidy_nb.*(future_occupation_fine_local_subsidy_nb[:,2].==2);
exit_cashcrop_to_staple_subsidy_nb = cash_crop_past_dist_subsidy_nb.*(future_occupation_fine_local_subsidy_nb[:,3].==2);
current_staple_subsidy_nb = entrants_staple_from_workers_subsidy_nb + incumbents_staple_subsidy_nb + exit_cashcrop_to_staple_subsidy_nb;

entrants_cashcrop_from_workers_subsidy_nb = worker_past_dist_subsidy_nb.*(future_occupation_fine_local_subsidy_nb[:,1].==3);
entrants_from_staple_to_cashcrop_subsidy_nb= staple_past_dist_subsidy_nb.*(future_occupation_fine_local_subsidy_nb[:,2].==3);
incumbents_cashcrop_subsidy_nb = cash_crop_past_dist_subsidy_nb.*(future_occupation_fine_local_subsidy_nb[:,3].==3);
current_cashcrop_subsidy_nb = entrants_cashcrop_from_workers_subsidy_nb + incumbents_cashcrop_subsidy_nb + entrants_from_staple_to_cashcrop_subsidy_nb;

undernutritioned_workers_subsidy_nb = sum( (c_S_W_fine_subsidy_nb[:,1].<cons_level_substinence).*stay_workers_subsidy_nb
    + (c_S_W_fine_subsidy_nb[:,2].<cons_level_substinence).*exit_staple_to_work_subsidy_nb +
(c_S_W_fine_subsidy_nb[:,3].<cons_level_substinence) .*exit_cashcrop_to_work_subsidy_nb );
undernutritioned_staple_farmer_subsidy_nb = sum( (c_S_S_fine_subsidy_nb[:,1].<cons_level_substinence).*entrants_staple_from_workers_subsidy_nb +
(c_S_S_fine_subsidy_nb[:,2].<cons_level_substinence).*incumbents_staple_subsidy_nb +
(c_S_S_fine_subsidy_nb[:,3].<cons_level_substinence) .*exit_cashcrop_to_staple_subsidy_nb );
undernutritioned_cashcrop_farmer_subsidy_nb = sum( (c_S_B_fine_subsidy_nb[:,1].<cons_level_substinence).*entrants_cashcrop_from_workers_subsidy_nb
    + (c_S_B_fine_subsidy_nb[:,2].<cons_level_substinence).*entrants_from_staple_to_cashcrop_subsidy_nb +
(c_S_B_fine_subsidy_nb[:,3].<cons_level_substinence) .*incumbents_cashcrop_subsidy_nb );

undernourished_subsidy_nb = undernutritioned_workers_subsidy_nb + undernutritioned_staple_farmer_subsidy_nb +undernutritioned_cashcrop_farmer_subsidy_nb;

# Infrastructure undernourished
worker_past_dist_inf = stat_distr_inf[(Baseline_parameter.ns_fine *0 + 1):(Baseline_parameter.ns_fine *1)];
staple_past_dist_inf = stat_distr_inf[(Baseline_parameter.ns_fine *1 + 1):(Baseline_parameter.ns_fine *2)];
cash_crop_past_dist_inf = stat_distr_inf[(Baseline_parameter.ns_fine *2 + 1):(Baseline_parameter.ns_fine *3)];

stay_workers_inf = worker_past_dist_inf.*(future_occupation_fine_local_inf[:,1].==1);
exit_staple_to_work_inf = staple_past_dist_inf.*(future_occupation_fine_local_inf[:,2].==1);
exit_cashcrop_to_work_inf = cash_crop_past_dist_inf.*(future_occupation_fine_local_inf[:,3].==1);
current_workers_inf = stay_workers_inf + exit_staple_to_work_inf + exit_cashcrop_to_work_inf;

entrants_staple_from_workers_inf = worker_past_dist_inf.*(future_occupation_fine_local_inf[:,1].==2);
incumbents_staple_inf = staple_past_dist_inf.*(future_occupation_fine_local_inf[:,2].==2);
exit_cashcrop_to_staple_inf = cash_crop_past_dist_inf.*(future_occupation_fine_local_inf[:,3].==2);
current_staple_inf = entrants_staple_from_workers_inf + incumbents_staple_inf + exit_cashcrop_to_staple_inf;

entrants_cashcrop_from_workers_inf = worker_past_dist_inf.*(future_occupation_fine_local_inf[:,1].==3);
entrants_from_staple_to_cashcrop_inf= staple_past_dist_inf.*(future_occupation_fine_local_inf[:,2].==3);
incumbents_cashcrop_inf = cash_crop_past_dist_inf.*(future_occupation_fine_local_inf[:,3].==3);
current_cashcrop_inf = entrants_cashcrop_from_workers_inf + incumbents_cashcrop_inf + entrants_from_staple_to_cashcrop_inf;

undernutritioned_workers_inf = sum( (c_S_W_fine_inf[:,1].<cons_level_substinence).*stay_workers_inf
    + (c_S_W_fine_inf[:,2].<cons_level_substinence).*exit_staple_to_work_inf +
(c_S_W_fine_inf[:,3].<cons_level_substinence) .*exit_cashcrop_to_work_inf );
undernutritioned_staple_farmer_inf = sum( (c_S_S_fine_inf[:,1].<cons_level_substinence).*entrants_staple_from_workers_inf +
(c_S_S_fine_inf[:,2].<cons_level_substinence).*incumbents_staple_inf +
(c_S_S_fine_inf[:,3].<cons_level_substinence) .*exit_cashcrop_to_staple_inf );
undernutritioned_cashcrop_farmer_inf = sum( (c_S_B_fine_inf[:,1].<cons_level_substinence).*entrants_cashcrop_from_workers_inf
    + (c_S_B_fine_inf[:,2].<cons_level_substinence).*entrants_from_staple_to_cashcrop_inf +
(c_S_B_fine_inf[:,3].<cons_level_substinence) .*incumbents_cashcrop_inf );

undernourished_inf = undernutritioned_workers_inf + undernutritioned_staple_farmer_inf +undernutritioned_cashcrop_farmer_inf;

# Infrastructure with spillovers undernourished
worker_past_dist_inf_sp = stat_distr_inf_sp[(Baseline_parameter.ns_fine *0 + 1):(Baseline_parameter.ns_fine *1)];
staple_past_dist_inf_sp = stat_distr_inf_sp[(Baseline_parameter.ns_fine *1 + 1):(Baseline_parameter.ns_fine *2)];
cash_crop_past_dist_inf_sp = stat_distr_inf_sp[(Baseline_parameter.ns_fine *2 + 1):(Baseline_parameter.ns_fine *3)];

stay_workers_inf_sp = worker_past_dist_inf_sp.*(future_occupation_fine_local_inf_sp[:,1].==1);
exit_staple_to_work_inf_sp = staple_past_dist_inf_sp.*(future_occupation_fine_local_inf_sp[:,2].==1);
exit_cashcrop_to_work_inf_sp = cash_crop_past_dist_inf_sp.*(future_occupation_fine_local_inf_sp[:,3].==1);
current_workers_inf_sp = stay_workers_inf_sp + exit_staple_to_work_inf_sp + exit_cashcrop_to_work_inf_sp;

entrants_staple_from_workers_inf_sp = worker_past_dist_inf_sp.*(future_occupation_fine_local_inf_sp[:,1].==2);
incumbents_staple_inf_sp = staple_past_dist_inf_sp.*(future_occupation_fine_local_inf_sp[:,2].==2);
exit_cashcrop_to_staple_inf_sp = cash_crop_past_dist_inf_sp.*(future_occupation_fine_local_inf_sp[:,3].==2);
current_staple_inf_sp = entrants_staple_from_workers_inf_sp + incumbents_staple_inf_sp + exit_cashcrop_to_staple_inf_sp;

entrants_cashcrop_from_workers_inf_sp = worker_past_dist_inf_sp.*(future_occupation_fine_local_inf_sp[:,1].==3);
entrants_from_staple_to_cashcrop_inf_sp= staple_past_dist_inf_sp.*(future_occupation_fine_local_inf_sp[:,2].==3);
incumbents_cashcrop_inf_sp = cash_crop_past_dist_inf_sp.*(future_occupation_fine_local_inf_sp[:,3].==3);
current_cashcrop_inf_sp = entrants_cashcrop_from_workers_inf_sp + incumbents_cashcrop_inf_sp + entrants_from_staple_to_cashcrop_inf_sp;

undernutritioned_workers_inf_sp = sum( (c_S_W_fine_inf_sp[:,1].<cons_level_substinence).*stay_workers_inf_sp
    + (c_S_W_fine_inf_sp[:,2].<cons_level_substinence).*exit_staple_to_work_inf_sp +
(c_S_W_fine_inf_sp[:,3].<cons_level_substinence) .*exit_cashcrop_to_work_inf_sp );
undernutritioned_staple_farmer_inf_sp = sum( (c_S_S_fine_inf_sp[:,1].<cons_level_substinence).*entrants_staple_from_workers_inf_sp +
(c_S_S_fine_inf_sp[:,2].<cons_level_substinence).*incumbents_staple_inf_sp +
(c_S_S_fine_inf_sp[:,3].<cons_level_substinence) .*exit_cashcrop_to_staple_inf_sp );
undernutritioned_cashcrop_farmer_inf_sp = sum( (c_S_B_fine_inf_sp[:,1].<cons_level_substinence).*entrants_cashcrop_from_workers_inf_sp
    + (c_S_B_fine_inf_sp[:,2].<cons_level_substinence).*entrants_from_staple_to_cashcrop_inf_sp +
(c_S_B_fine_inf_sp[:,3].<cons_level_substinence) .*incumbents_cashcrop_inf_sp );

undernourished_inf_sp =  undernutritioned_workers_inf_sp + undernutritioned_staple_farmer_inf_sp +undernutritioned_cashcrop_farmer_inf_sp;

#Undernourished optimal subsidy with balanced budget

undernutritioned_subsidy_b_grid = zeros(no_taus);
food_share_cons_subsidy_b_grid = zeros(no_taus);
food_share_worker_subsidy_b_grid = zeros(no_taus);
food_share_staple_subsidy_b_grid = zeros(no_taus);
food_share_cashcrop_subsidy_b_grid = zeros(no_taus);
food_share_worker_subsidy_b_grid_consbased = zeros(Baseline_parameter.ns_fine,no_taus);
food_share_staple_subsidy_b_grid_consbased = zeros(Baseline_parameter.ns_fine,no_taus);
food_share_cashcrop_subsidy_b_grid_consbased = zeros(Baseline_parameter.ns_fine,no_taus);

manuf_share_worker_subsidy_b_grid_consbased = zeros(Baseline_parameter.ns_fine,no_taus);
manuf_share_staple_subsidy_b_grid_consbased = zeros(Baseline_parameter.ns_fine,no_taus);
manuf_share_cashcrop_subsidy_b_grid_consbased = zeros(Baseline_parameter.ns_fine,no_taus);

worker_subsidy_b_grid_consbased_ordered = zeros(Baseline_parameter.ns_fine,no_taus);
staple_subsidy_b_grid_consbased_ordered = zeros(Baseline_parameter.ns_fine,no_taus);
cashcrop_subsidy_b_grid_consbased_ordered = zeros(Baseline_parameter.ns_fine,no_taus);
# Distribution of past occupations
for ii = 1:no_taus
    worker_past_dist_subsidy_b_grid_tmp= stat_distr_subsidy_b_grid[(Baseline_parameter.ns_fine *0 + 1):(Baseline_parameter.ns_fine *1),ii];
    staple_past_dist_subsidy_b_grid_tmp= stat_distr_subsidy_b_grid[(Baseline_parameter.ns_fine *1 + 1):(Baseline_parameter.ns_fine *2),ii];
    cash_crop_past_dist_subsidy_b_grid_tmp= stat_distr_subsidy_b_grid[(Baseline_parameter.ns_fine *2 + 1):(Baseline_parameter.ns_fine *3),ii];

    stay_workers_subsidy_b_grid_tmp= worker_past_dist_subsidy_b_grid_tmp.*(future_occupation_fine_local_subsidy_b_grid[:,1,ii].==1);
    exit_staple_to_work_subsidy_b_grid_tmp= staple_past_dist_subsidy_b_grid_tmp.*(future_occupation_fine_local_subsidy_b_grid[:,2,ii].==1);
    exit_cashcrop_to_work_subsidy_b_grid_tmp= cash_crop_past_dist_subsidy_b_grid_tmp.*(future_occupation_fine_local_subsidy_b_grid[:,3,ii].==1);
    current_workers_subsidy_b_grid_tmp= stay_workers_subsidy_b_grid_tmp+ exit_staple_to_work_subsidy_b_grid_tmp+ exit_cashcrop_to_work_subsidy_b_grid_tmp;

    entrants_staple_from_workers_subsidy_b_grid_tmp= worker_past_dist_subsidy_b_grid_tmp.*(future_occupation_fine_local_subsidy_b_grid[:,1,ii].==2);
    incumbents_staple_subsidy_b_grid_tmp= staple_past_dist_subsidy_b_grid_tmp.*(future_occupation_fine_local_subsidy_b_grid[:,2,ii].==2);
    exit_cashcrop_to_staple_subsidy_b_grid_tmp= cash_crop_past_dist_subsidy_b_grid_tmp.*(future_occupation_fine_local_subsidy_b_grid[:,3,ii].==2);
    current_staple_subsidy_b_grid_tmp= entrants_staple_from_workers_subsidy_b_grid_tmp+ incumbents_staple_subsidy_b_grid_tmp+ exit_cashcrop_to_staple_subsidy_b_grid_tmp;

    entrants_cashcrop_from_workers_subsidy_b_grid_tmp= worker_past_dist_subsidy_b_grid_tmp.*(future_occupation_fine_local_subsidy_b_grid[:,1,ii].==3);
    entrants_from_staple_to_cashcrop_subsidy_b_grid_tmp= staple_past_dist_subsidy_b_grid_tmp.*(future_occupation_fine_local_subsidy_b_grid[:,2,ii].==3);
    incumbents_cashcrop_subsidy_b_grid_tmp= cash_crop_past_dist_subsidy_b_grid_tmp.*(future_occupation_fine_local_subsidy_b_grid[:,3,ii].==3);
    current_cashcrop_subsidy_b_grid_tmp= entrants_cashcrop_from_workers_subsidy_b_grid_tmp+ incumbents_cashcrop_subsidy_b_grid_tmp+ entrants_from_staple_to_cashcrop_subsidy_b_grid_tmp;

    undernutritioned_workers_subsidy_b_grid_tmp= sum( (c_S_W_fine_subsidy_b_grid[:,1,ii].<cons_level_substinence).*stay_workers_subsidy_b_grid_tmp
        + (c_S_W_fine_subsidy_b_grid[:,2,ii].<cons_level_substinence).*exit_staple_to_work_subsidy_b_grid_tmp+
    (c_S_W_fine_subsidy_b_grid[:,3,ii].<cons_level_substinence) .*exit_cashcrop_to_work_subsidy_b_grid_tmp);
    undernutritioned_staple_farmer_subsidy_b_grid_tmp= sum( (c_S_S_fine_subsidy_b_grid[:,1,ii].<cons_level_substinence).*entrants_staple_from_workers_subsidy_b_grid_tmp+
    (c_S_S_fine_subsidy_b_grid[:,2,ii].<cons_level_substinence).*incumbents_staple_subsidy_b_grid_tmp+
    (c_S_S_fine_subsidy_b_grid[:,3,ii].<cons_level_substinence) .*exit_cashcrop_to_staple_subsidy_b_grid_tmp);
    undernutritioned_cashcrop_farmer_subsidy_b_grid_tmp= sum( (c_S_B_fine_subsidy_b_grid[:,1,ii].<cons_level_substinence).*entrants_cashcrop_from_workers_subsidy_b_grid_tmp
        + (c_S_B_fine_subsidy_b_grid[:,2,ii].<cons_level_substinence).*entrants_from_staple_to_cashcrop_subsidy_b_grid_tmp+
    (c_S_B_fine_subsidy_b_grid[:,3,ii].<cons_level_substinence) .*incumbents_cashcrop_subsidy_b_grid_tmp);
    
    food_share_cons_workers_subsidy_b_grid_tmp= sum( (c_S_W_fine_subsidy_b_grid[:,1,ii]./(cons_fine_local_subsidy_b_grid[:,1,ii] )).*stay_workers_subsidy_b_grid_tmp
        + (c_S_W_fine_subsidy_b_grid[:,2,ii]./(cons_fine_local_subsidy_b_grid[:,2,ii] )).*exit_staple_to_work_subsidy_b_grid_tmp+
        (c_S_W_fine_subsidy_b_grid[:,3,ii]./(cons_fine_local_subsidy_b_grid[:,3,ii] )) .*exit_cashcrop_to_work_subsidy_b_grid_tmp)./sum(current_workers_subsidy_b_grid_tmp);

    food_share_cons_staple_farmer_subsidy_b_grid_tmp= sum( (c_S_S_fine_subsidy_b_grid[:,1,ii]./(cons_fine_local_subsidy_b_grid[:,1,ii] )).*entrants_staple_from_workers_subsidy_b_grid_tmp
        + (c_S_S_fine_subsidy_b_grid[:,2,ii]./(cons_fine_local_subsidy_b_grid[:,2,ii] )).*incumbents_staple_subsidy_b_grid_tmp+
        (c_S_S_fine_subsidy_b_grid[:,3,ii]./(cons_fine_local_subsidy_b_grid[:,3,ii] )) .*exit_cashcrop_to_staple_subsidy_b_grid_tmp)./sum(current_staple_subsidy_b_grid_tmp);
    
    food_share_cons_cashcrop_farmer_subsidy_b_grid_tmp= sum( (c_S_B_fine_subsidy_b_grid[:,1,ii]./(cons_fine_local_subsidy_b_grid[:,1,ii] )).*entrants_cashcrop_from_workers_subsidy_b_grid_tmp
        + (c_S_B_fine_subsidy_b_grid[:,2,ii]./(cons_fine_local_subsidy_b_grid[:,2,ii] )).*entrants_from_staple_to_cashcrop_subsidy_b_grid_tmp+
        (c_S_B_fine_subsidy_b_grid[:,3,ii]./(cons_fine_local_subsidy_b_grid[:,3,ii] )) .*incumbents_cashcrop_subsidy_b_grid_tmp)./sum(current_cashcrop_subsidy_b_grid_tmp);
    
    
    consumption_sort_index_worker = sortperm(cons_fine_local_subsidy_b_grid[:,1,ii]);
    food_share_worker_subsidy_b_grid_consbased[:,ii] = c_S_W_fine_subsidy_b_grid[consumption_sort_index_worker,1,ii]./(cons_fine_local_subsidy_b_grid[consumption_sort_index_worker,1,ii] );
    manuf_share_worker_subsidy_b_grid_consbased[:,ii] = prices_subsidy_b_grid[2,ii]*c_M_W_fine_subsidy_b_grid[consumption_sort_index_worker,1,ii]./(cons_fine_local_subsidy_b_grid[consumption_sort_index_worker,1,ii] );
    consumption_sort_index_staple = sortperm(cons_fine_local_subsidy_b_grid[:,2,ii]);
    food_share_staple_subsidy_b_grid_consbased[:,ii] = c_S_S_fine_subsidy_b_grid[consumption_sort_index_staple,2,ii]./(cons_fine_local_subsidy_b_grid[consumption_sort_index_staple,2,ii] );
    manuf_share_staple_subsidy_b_grid_consbased[:,ii] = prices_subsidy_b_grid[2,ii]*c_M_S_fine_subsidy_b_grid[consumption_sort_index_staple,2,ii]./(cons_fine_local_subsidy_b_grid[consumption_sort_index_staple,2,ii] );
    consumption_sort_index_cashcrop = sortperm(cons_fine_local_subsidy_b_grid[:,3,ii]);
    food_share_cashcrop_subsidy_b_grid_consbased[:,ii] = c_S_B_fine_subsidy_b_grid[consumption_sort_index_cashcrop,3,ii]./(cons_fine_local_subsidy_b_grid[consumption_sort_index_cashcrop,3,ii] );
    manuf_share_cashcrop_subsidy_b_grid_consbased[:,ii] = prices_subsidy_b_grid[2,ii]*c_M_B_fine_subsidy_b_grid[consumption_sort_index_cashcrop,3,ii]./(cons_fine_local_subsidy_b_grid[consumption_sort_index_cashcrop,3,ii] );
    # plot((1+Q_S)*c_S_S_fine_subsidy_b_grid[:,1,ii]./(cons_fine_local_subsidy_b_grid[:,1,ii] )) at what price this would be fair, I have no idea
    # plot!(prices_subsidy_b_grid[1,ii]*c_B_S_fine_subsidy_b_grid[:,1,ii]./(cons_fine_local_subsidy_b_grid[:,1,ii] ))
    # plot!(prices_subsidy_b_grid[2,ii]*c_M_S_fine_subsidy_b_grid[:,1,ii]./(cons_fine_local_subsidy_b_grid[:,1,ii] ))
    undernutritioned_subsidy_b_grid[ii] = undernutritioned_workers_subsidy_b_grid_tmp+ undernutritioned_staple_farmer_subsidy_b_grid_tmp+undernutritioned_cashcrop_farmer_subsidy_b_grid_tmp
    food_share_cons_subsidy_b_grid[ii] = food_share_cons_workers_subsidy_b_grid_tmp*sum(current_workers_subsidy_b_grid_tmp)+ food_share_cons_staple_farmer_subsidy_b_grid_tmp*sum(current_staple_subsidy_b_grid_tmp) + food_share_cons_cashcrop_farmer_subsidy_b_grid_tmp*sum(current_cashcrop_subsidy_b_grid_tmp);
    food_share_worker_subsidy_b_grid[ii] = food_share_cons_workers_subsidy_b_grid_tmp;
    food_share_staple_subsidy_b_grid[ii] = food_share_cons_staple_farmer_subsidy_b_grid_tmp;
    food_share_cashcrop_subsidy_b_grid[ii] = food_share_cons_cashcrop_farmer_subsidy_b_grid_tmp;
end
undernutritioned_subsidy_b_grid_growth_smoothed = movmean((100*(undernutritioned_subsidy_b_grid./undernutritioned_subsidy_b_grid[no_taus]) .- 100),5)
undernutritioned_subsidy_b_grid_growth_smoothed = undernutritioned_subsidy_b_grid_growth_smoothed.-undernutritioned_subsidy_b_grid_growth_smoothed[no_taus];

food_share_cons_subsidy_b_grid_growth_smoothed = movmean((100*(food_share_cons_subsidy_b_grid./food_share_cons_subsidy_b_grid[no_taus]) .- 100),5)
food_share_cons_subsidy_b_grid_growth_smoothed = food_share_cons_subsidy_b_grid_growth_smoothed.-food_share_cons_subsidy_b_grid_growth_smoothed[no_taus];

food_share_worker_subsidy_b_grid_growth_smoothed = movmean((100*(food_share_worker_subsidy_b_grid./food_share_worker_subsidy_b_grid[no_taus]) .- 100),5)
food_share_worker_subsidy_b_grid_growth_smoothed = food_share_worker_subsidy_b_grid_growth_smoothed.-food_share_worker_subsidy_b_grid_growth_smoothed[no_taus];

food_share_staple_subsidy_b_grid_growth_smoothed = movmean((100*(food_share_staple_subsidy_b_grid./food_share_staple_subsidy_b_grid[no_taus]) .- 100),5)
food_share_staple_subsidy_b_grid_growth_smoothed = food_share_staple_subsidy_b_grid_growth_smoothed.-food_share_staple_subsidy_b_grid_growth_smoothed[no_taus];

food_share_cashcrop_subsidy_b_grid_growth_smoothed = movmean((100*(food_share_cashcrop_subsidy_b_grid./food_share_cashcrop_subsidy_b_grid[no_taus]) .- 100),5)
food_share_cashcrop_subsidy_b_grid_growth_smoothed = food_share_cashcrop_subsidy_b_grid_growth_smoothed.-food_share_cashcrop_subsidy_b_grid_growth_smoothed[no_taus];

#Undernourished optimal subsidy without balanced budget

undernutritioned_subsidy_nb_grid = zeros(no_taus);
food_share_cons_subsidy_nb_grid = zeros(no_taus);
food_share_worker_subsidy_nb_grid = zeros(no_taus);
food_share_staple_subsidy_nb_grid = zeros(no_taus);
food_share_cashcrop_subsidy_nb_grid = zeros(no_taus);
food_share_worker_subsidy_nb_grid_consbased = zeros(Baseline_parameter.ns_fine,no_taus);
food_share_staple_subsidy_nb_grid_consbased = zeros(Baseline_parameter.ns_fine,no_taus);
food_share_cashcrop_subsidy_nb_grid_consbased = zeros(Baseline_parameter.ns_fine,no_taus);

manuf_share_worker_subsidy_nb_grid_consbased = zeros(Baseline_parameter.ns_fine,no_taus);
manuf_share_staple_subsidy_nb_grid_consbased = zeros(Baseline_parameter.ns_fine,no_taus);
manuf_share_cashcrop_subsidy_nb_grid_consbased = zeros(Baseline_parameter.ns_fine,no_taus);

worker_subsidy_nb_grid_consbased_ordered = zeros(Baseline_parameter.ns_fine,no_taus);
staple_subsidy_nb_grid_consbased_ordered = zeros(Baseline_parameter.ns_fine,no_taus);
cashcrop_subsidy_nb_grid_consbased_ordered = zeros(Baseline_parameter.ns_fine,no_taus);
# Distribution of past occupations
for ii = 1:no_taus
    worker_past_dist_subsidy_nb_grid_tmp= stat_distr_subsidy_nb_grid[(Baseline_parameter.ns_fine *0 + 1):(Baseline_parameter.ns_fine *1),ii];
    staple_past_dist_subsidy_nb_grid_tmp= stat_distr_subsidy_nb_grid[(Baseline_parameter.ns_fine *1 + 1):(Baseline_parameter.ns_fine *2),ii];
    cash_crop_past_dist_subsidy_nb_grid_tmp= stat_distr_subsidy_nb_grid[(Baseline_parameter.ns_fine *2 + 1):(Baseline_parameter.ns_fine *3),ii];

    stay_workers_subsidy_nb_grid_tmp= worker_past_dist_subsidy_nb_grid_tmp.*(future_occupation_fine_local_subsidy_nb_grid[:,1,ii].==1);
    exit_staple_to_work_subsidy_nb_grid_tmp= staple_past_dist_subsidy_nb_grid_tmp.*(future_occupation_fine_local_subsidy_nb_grid[:,2,ii].==1);
    exit_cashcrop_to_work_subsidy_nb_grid_tmp= cash_crop_past_dist_subsidy_nb_grid_tmp.*(future_occupation_fine_local_subsidy_nb_grid[:,3,ii].==1);
    current_workers_subsidy_nb_grid_tmp= stay_workers_subsidy_nb_grid_tmp+ exit_staple_to_work_subsidy_nb_grid_tmp+ exit_cashcrop_to_work_subsidy_nb_grid_tmp;

    entrants_staple_from_workers_subsidy_nb_grid_tmp= worker_past_dist_subsidy_nb_grid_tmp.*(future_occupation_fine_local_subsidy_nb_grid[:,1,ii].==2);
    incumbents_staple_subsidy_nb_grid_tmp= staple_past_dist_subsidy_nb_grid_tmp.*(future_occupation_fine_local_subsidy_nb_grid[:,2,ii].==2);
    exit_cashcrop_to_staple_subsidy_nb_grid_tmp= cash_crop_past_dist_subsidy_nb_grid_tmp.*(future_occupation_fine_local_subsidy_nb_grid[:,3,ii].==2);
    current_staple_subsidy_nb_grid_tmp= entrants_staple_from_workers_subsidy_nb_grid_tmp+ incumbents_staple_subsidy_nb_grid_tmp+ exit_cashcrop_to_staple_subsidy_nb_grid_tmp;

    entrants_cashcrop_from_workers_subsidy_nb_grid_tmp= worker_past_dist_subsidy_nb_grid_tmp.*(future_occupation_fine_local_subsidy_nb_grid[:,1,ii].==3);
    entrants_from_staple_to_cashcrop_subsidy_nb_grid_tmp= staple_past_dist_subsidy_nb_grid_tmp.*(future_occupation_fine_local_subsidy_nb_grid[:,2,ii].==3);
    incumbents_cashcrop_subsidy_nb_grid_tmp= cash_crop_past_dist_subsidy_nb_grid_tmp.*(future_occupation_fine_local_subsidy_nb_grid[:,3,ii].==3);
    current_cashcrop_subsidy_nb_grid_tmp= entrants_cashcrop_from_workers_subsidy_nb_grid_tmp+ incumbents_cashcrop_subsidy_nb_grid_tmp+ entrants_from_staple_to_cashcrop_subsidy_nb_grid_tmp;

    undernutritioned_workers_subsidy_nb_grid_tmp= sum( (c_S_W_fine_subsidy_nb_grid[:,1,ii].<cons_level_substinence).*stay_workers_subsidy_nb_grid_tmp
        + (c_S_W_fine_subsidy_nb_grid[:,2,ii].<cons_level_substinence).*exit_staple_to_work_subsidy_nb_grid_tmp+
    (c_S_W_fine_subsidy_nb_grid[:,3,ii].<cons_level_substinence) .*exit_cashcrop_to_work_subsidy_nb_grid_tmp);
    undernutritioned_staple_farmer_subsidy_nb_grid_tmp= sum( (c_S_S_fine_subsidy_nb_grid[:,1,ii].<cons_level_substinence).*entrants_staple_from_workers_subsidy_nb_grid_tmp+
    (c_S_S_fine_subsidy_nb_grid[:,2,ii].<cons_level_substinence).*incumbents_staple_subsidy_nb_grid_tmp+
    (c_S_S_fine_subsidy_nb_grid[:,3,ii].<cons_level_substinence) .*exit_cashcrop_to_staple_subsidy_nb_grid_tmp);
    undernutritioned_cashcrop_farmer_subsidy_nb_grid_tmp= sum( (c_S_B_fine_subsidy_nb_grid[:,1,ii].<cons_level_substinence).*entrants_cashcrop_from_workers_subsidy_nb_grid_tmp
        + (c_S_B_fine_subsidy_nb_grid[:,2,ii].<cons_level_substinence).*entrants_from_staple_to_cashcrop_subsidy_nb_grid_tmp+
    (c_S_B_fine_subsidy_nb_grid[:,3,ii].<cons_level_substinence) .*incumbents_cashcrop_subsidy_nb_grid_tmp);
    
    food_share_cons_workers_subsidy_nb_grid_tmp= sum( (c_S_W_fine_subsidy_nb_grid[:,1,ii]./(cons_fine_local_subsidy_nb_grid[:,1,ii] )).*stay_workers_subsidy_nb_grid_tmp
        + (c_S_W_fine_subsidy_nb_grid[:,2,ii]./(cons_fine_local_subsidy_nb_grid[:,2,ii] )).*exit_staple_to_work_subsidy_nb_grid_tmp+
        (c_S_W_fine_subsidy_nb_grid[:,3,ii]./(cons_fine_local_subsidy_nb_grid[:,3,ii] )) .*exit_cashcrop_to_work_subsidy_nb_grid_tmp)./sum(current_workers_subsidy_nb_grid_tmp);

    food_share_cons_staple_farmer_subsidy_nb_grid_tmp= sum( (c_S_S_fine_subsidy_nb_grid[:,1,ii]./(cons_fine_local_subsidy_nb_grid[:,1,ii] )).*entrants_staple_from_workers_subsidy_nb_grid_tmp
        + (c_S_S_fine_subsidy_nb_grid[:,2,ii]./(cons_fine_local_subsidy_nb_grid[:,2,ii] )).*incumbents_staple_subsidy_nb_grid_tmp+
        (c_S_S_fine_subsidy_nb_grid[:,3,ii]./(cons_fine_local_subsidy_nb_grid[:,3,ii] )) .*exit_cashcrop_to_staple_subsidy_nb_grid_tmp)./sum(current_staple_subsidy_nb_grid_tmp);
    
    food_share_cons_cashcrop_farmer_subsidy_nb_grid_tmp= sum( (c_S_B_fine_subsidy_nb_grid[:,1,ii]./(cons_fine_local_subsidy_nb_grid[:,1,ii] )).*entrants_cashcrop_from_workers_subsidy_nb_grid_tmp
        + (c_S_B_fine_subsidy_nb_grid[:,2,ii]./(cons_fine_local_subsidy_nb_grid[:,2,ii] )).*entrants_from_staple_to_cashcrop_subsidy_nb_grid_tmp+
        (c_S_B_fine_subsidy_nb_grid[:,3,ii]./(cons_fine_local_subsidy_nb_grid[:,3,ii] )) .*incumbents_cashcrop_subsidy_nb_grid_tmp)./sum(current_cashcrop_subsidy_nb_grid_tmp);
    
    
    consumption_sort_index_worker = sortperm(cons_fine_local_subsidy_nb_grid[:,1,ii]);
    food_share_worker_subsidy_nb_grid_consbased[:,ii] = c_S_W_fine_subsidy_nb_grid[consumption_sort_index_worker,1,ii]./(cons_fine_local_subsidy_nb_grid[consumption_sort_index_worker,1,ii] );
    manuf_share_worker_subsidy_nb_grid_consbased[:,ii] = prices_subsidy_nb_grid[2,ii]*c_M_W_fine_subsidy_nb_grid[consumption_sort_index_worker,1,ii]./(cons_fine_local_subsidy_nb_grid[consumption_sort_index_worker,1,ii] );
    consumption_sort_index_staple = sortperm(cons_fine_local_subsidy_nb_grid[:,2,ii]);
    food_share_staple_subsidy_nb_grid_consbased[:,ii] = c_S_S_fine_subsidy_nb_grid[consumption_sort_index_staple,2,ii]./(cons_fine_local_subsidy_nb_grid[consumption_sort_index_staple,2,ii] );
    manuf_share_staple_subsidy_nb_grid_consbased[:,ii] = prices_subsidy_nb_grid[2,ii]*c_M_S_fine_subsidy_nb_grid[consumption_sort_index_staple,2,ii]./(cons_fine_local_subsidy_nb_grid[consumption_sort_index_staple,2,ii] );
    consumption_sort_index_cashcrop = sortperm(cons_fine_local_subsidy_nb_grid[:,3,ii]);
    food_share_cashcrop_subsidy_nb_grid_consbased[:,ii] = c_S_B_fine_subsidy_nb_grid[consumption_sort_index_cashcrop,3,ii]./(cons_fine_local_subsidy_nb_grid[consumption_sort_index_cashcrop,3,ii] );
    manuf_share_cashcrop_subsidy_nb_grid_consbased[:,ii] = prices_subsidy_nb_grid[2,ii]*c_M_B_fine_subsidy_nb_grid[consumption_sort_index_cashcrop,3,ii]./(cons_fine_local_subsidy_nb_grid[consumption_sort_index_cashcrop,3,ii] );
    # plot((1+Q_S)*c_S_S_fine_subsidy_nb_grid[:,1,ii]./(cons_fine_local_subsidy_nb_grid[:,1,ii] )) at what price this would be fair, I have no idea
    # plot!(prices_subsidy_nb_grid[1,ii]*c_B_S_fine_subsidy_nb_grid[:,1,ii]./(cons_fine_local_subsidy_nb_grid[:,1,ii] ))
    # plot!(prices_subsidy_nb_grid[2,ii]*c_M_S_fine_subsidy_nb_grid[:,1,ii]./(cons_fine_local_subsidy_nb_grid[:,1,ii] ))
    undernutritioned_subsidy_nb_grid[ii] = undernutritioned_workers_subsidy_nb_grid_tmp+ undernutritioned_staple_farmer_subsidy_nb_grid_tmp+undernutritioned_cashcrop_farmer_subsidy_nb_grid_tmp
    food_share_cons_subsidy_nb_grid[ii] = food_share_cons_workers_subsidy_nb_grid_tmp*sum(current_workers_subsidy_nb_grid_tmp)+ food_share_cons_staple_farmer_subsidy_nb_grid_tmp*sum(current_staple_subsidy_nb_grid_tmp) + food_share_cons_cashcrop_farmer_subsidy_nb_grid_tmp*sum(current_cashcrop_subsidy_nb_grid_tmp);
    food_share_worker_subsidy_nb_grid[ii] = food_share_cons_workers_subsidy_nb_grid_tmp;
    food_share_staple_subsidy_nb_grid[ii] = food_share_cons_staple_farmer_subsidy_nb_grid_tmp;
    food_share_cashcrop_subsidy_nb_grid[ii] = food_share_cons_cashcrop_farmer_subsidy_nb_grid_tmp;
end
undernutritioned_subsidy_nb_grid_growth_smoothed = movmean((100*(undernutritioned_subsidy_nb_grid./undernutritioned_subsidy_nb_grid[no_taus]) .- 100),5)
undernutritioned_subsidy_nb_grid_growth_smoothed = undernutritioned_subsidy_nb_grid_growth_smoothed.-undernutritioned_subsidy_nb_grid_growth_smoothed[no_taus];

food_share_cons_subsidy_nb_grid_growth_smoothed = movmean((100*(food_share_cons_subsidy_nb_grid./food_share_cons_subsidy_nb_grid[no_taus]) .- 100),5)
food_share_cons_subsidy_nb_grid_growth_smoothed = food_share_cons_subsidy_nb_grid_growth_smoothed.-food_share_cons_subsidy_nb_grid_growth_smoothed[no_taus];

food_share_worker_subsidy_nb_grid_growth_smoothed = movmean((100*(food_share_worker_subsidy_nb_grid./food_share_worker_subsidy_nb_grid[no_taus]) .- 100),5)
food_share_worker_subsidy_nb_grid_growth_smoothed = food_share_worker_subsidy_nb_grid_growth_smoothed.-food_share_worker_subsidy_nb_grid_growth_smoothed[no_taus];

food_share_staple_subsidy_nb_grid_growth_smoothed = movmean((100*(food_share_staple_subsidy_nb_grid./food_share_staple_subsidy_nb_grid[no_taus]) .- 100),5)
food_share_staple_subsidy_nb_grid_growth_smoothed = food_share_staple_subsidy_nb_grid_growth_smoothed.-food_share_staple_subsidy_nb_grid_growth_smoothed[no_taus];

food_share_cashcrop_subsidy_nb_grid_growth_smoothed = movmean((100*(food_share_cashcrop_subsidy_nb_grid./food_share_cashcrop_subsidy_nb_grid[no_taus]) .- 100),5)
food_share_cashcrop_subsidy_nb_grid_growth_smoothed = food_share_cashcrop_subsidy_nb_grid_growth_smoothed.-food_share_cashcrop_subsidy_nb_grid_growth_smoothed[no_taus];

# Storage for the optimal transition numbers
T = 50
price_trans_actual_optimal = zeros(3,T,no_taus)
residual_store_optimal = zeros(6,T,no_taus)
distr_store_optimal= zeros(3*Baseline_parameter.ns_fine,T+2,no_taus)
coeff_store_optimal= zeros(Baseline_parameter.ns,3,T+2,no_taus)
prod_staple_store_optimal = zeros(T,no_taus)
prod_cashcrop_store_optimal= zeros(T,no_taus)
prod_manuf_store_optimal= zeros(T,no_taus)
asset_supply_store_optimal= zeros(T,no_taus)
current_worker_pop_store_optimal= zeros(T,no_taus)
current_staple_pop_store_optimal= zeros(T,no_taus)
current_cashcrop_pop_store_optimal= zeros(T,no_taus)
current_account_residual_store_optimal= zeros(T,no_taus)
staple_productivity_store_optimal= zeros(T,no_taus)
cashcrop_productivity_store_optimal= zeros(T,no_taus)
manuf_productivity_store_optimal= zeros(T,no_taus)
aggregate_consumption_store_optimal= zeros(T,no_taus)
relative_land_to_cashcrop_store_optimal= zeros(T,no_taus)
mean_land_share_staples_store_optimal= zeros(T,no_taus)
undernourished_store_optimal= zeros(T,no_taus)
fertilizer_use_store_optimal= zeros(T,no_taus)
APG_store_optimal= zeros(T,no_taus)
var_APland_store_optimal= zeros(T,no_taus)
var_MPX_store_optimal= zeros(T,no_taus)
avg_labor_prod_rural_store_optimal= zeros(T,no_taus)
avg_labor_prod_urban_store_optimal= zeros(T,no_taus)
avg_agri_prod_rural_store_optimal= zeros(T,no_taus)
avg_agri_prod_urban_store_optimal= zeros(T,no_taus)
V_saved_store_optimal= zeros(3*Baseline_parameter.ns_fine,T+2,no_taus)
a_prime_fine_store_optimal= zeros(Baseline_parameter.ns_fine,3,T,no_taus)

welfare_trans_optimal= zeros(no_taus)
welfare_trans_optimal_real= zeros(no_taus)
welfare_trans_optimal_real_alt= zeros(no_taus)

# Storage for the epsilon transition numbers

price_trans_actual_epsilon_trans = zeros(4,T,no_epsilons)
residual_store_epsilon_trans = zeros(6,T,no_epsilons)
distr_store_epsilon_trans= zeros(3*Baseline_parameter.ns_fine,T+2,no_epsilons)
coeff_store_epsilon_trans= zeros(Baseline_parameter.ns,3,T+2,no_epsilons)
prod_staple_store_epsilon_trans = zeros(T,no_epsilons)
prod_cashcrop_store_epsilon_trans= zeros(T,no_epsilons)
prod_manuf_store_epsilon_trans= zeros(T,no_epsilons)
asset_supply_store_epsilon_trans= zeros(T,no_epsilons)
current_worker_pop_store_epsilon_trans= zeros(T,no_epsilons)
current_staple_pop_store_epsilon_trans= zeros(T,no_epsilons)
current_cashcrop_pop_store_epsilon_trans= zeros(T,no_epsilons)
current_account_residual_store_epsilon_trans= zeros(T,no_epsilons)
staple_productivity_store_epsilon_trans= zeros(T,no_epsilons)
cashcrop_productivity_store_epsilon_trans= zeros(T,no_epsilons)
manuf_productivity_store_epsilon_trans= zeros(T,no_epsilons)
aggregate_consumption_store_epsilon_trans= zeros(T,no_epsilons)
relative_land_to_cashcrop_store_epsilon_trans= zeros(T,no_epsilons)
mean_land_share_staples_store_epsilon_trans= zeros(T,no_epsilons)
undernourished_store_epsilon_trans= zeros(T,no_epsilons)
fertilizer_use_store_epsilon_trans= zeros(T,no_epsilons)
APG_store_epsilon_trans= zeros(T,no_epsilons)
var_APland_store_epsilon_trans= zeros(T,no_epsilons)
var_MPX_store_epsilon_trans= zeros(T,no_epsilons)
avg_labor_prod_rural_store_epsilon_trans= zeros(T,no_epsilons)
avg_labor_prod_urban_store_epsilon_trans= zeros(T,no_epsilons)
avg_agri_prod_rural_store_epsilon_trans= zeros(T,no_epsilons)
avg_agri_prod_urban_store_epsilon_trans= zeros(T,no_epsilons)
V_saved_store_epsilon_trans= zeros(3*Baseline_parameter.ns_fine,T+2,no_epsilons)
a_prime_fine_store_epsilon_trans= zeros(Baseline_parameter.ns_fine,3,T,no_epsilons)

welfare_trans_epsilon_trans= zeros(no_epsilons)
welfare_trans_epsilon_trans_real= zeros(no_epsilons)
welfare_trans_epsilon_trans_real_alt= zeros(no_epsilons)