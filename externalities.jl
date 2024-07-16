include("details_model_externality.jl")
include("externalities_transition.jl")
### Do it in a matrix so its easier to get a finer grid and run transition_functions
cons_level_substinence = 0.109;
No_subsidy_parameter = copy(Baseline_parameter);
No_subsidy_parameter.τ_S=-0;
balanced_share = 0.0
foreign_supply_capital_subsidy_b = 7.631597488627073;
prices_epsilon_grid = [1.188283177632008 0.18349666252335484 0.2832894357711001;
1.1906660926758061 0.18476248806422077 0.28395740826966315;
1.1925825472322047 0.18727339112307217 0.2848214354714461;
1.195476951761732 0.20085571586329348 0.28334980569254575;
1.2145275998907963 0.2147096472250037 0.2936424726929504;
1.2093001817161375 0.2130922258686436 0.2979426218123759]'
no_epsilons = size(prices_epsilon_grid)[2];
epsilon_grid = [0,0.1155, 0.25, 0.5, 1.0 ,1.5]

residual_goods_epsilon_grid = ones(6,no_epsilons);
stat_distr_epsilon_grid = ones(Baseline_parameter.ns_fine*3,no_epsilons);
cons_fine_local_epsilon_grid = ones(Baseline_parameter.ns_fine,3,no_epsilons);
a_prime_fine_local_epsilon_grid = ones(Baseline_parameter.ns_fine,3,no_epsilons);
future_occupation_fine_local_epsilon_grid = ones(Baseline_parameter.ns_fine,3,no_epsilons);
x_S_S_fine_epsilon_grid = ones(Baseline_parameter.ns_fine,3,no_epsilons);
x_SC_fine_epsilon_grid = ones(Baseline_parameter.ns_fine,3,no_epsilons);
x_BC_fine_epsilon_grid = ones(Baseline_parameter.ns_fine,3,no_epsilons);
coeff_epsilon_grid = ones(Baseline_parameter.ns,3,no_epsilons);
transaction_cost_loss_epsilon_grid = ones(no_epsilons);
nominal_GDP_epsilon_grid= ones(no_epsilons);
welfare_val_epsilon_grid =  ones(Baseline_parameter.ns_fine*3,no_epsilons);
Import_value_epsilon_grid = ones(no_epsilons);
Export_value_epsilon_grid = ones(no_epsilons);
current_worker_pop_epsilon_grid = ones(no_epsilons);
current_staple_pop_epsilon_grid = ones(no_epsilons);
current_cashcrop_pop_epsilon_grid = ones(no_epsilons);
marketable_agr_surplus_share_epsilon_grid = ones(no_epsilons);
exportshare_cashcrop_epsilon_grid = ones(no_epsilons);
fraction_model_epsilon_grid= ones(no_epsilons);
program_spending_epsilon_grid= ones(no_epsilons);
prod_value_improvement_epsilon_grid= ones(no_epsilons);
share_selling_increase_epsilon_grid= ones(no_epsilons);
exp_ratio_model_epsilon_grid= ones(no_epsilons);
mig_rate_model_epsilon_grid= ones(no_epsilons);
rural_pop_only_staples_model_epsilon_grid= ones(no_epsilons);
rural_pop_only_cashcrop_model_epsilon_grid= ones(no_epsilons);
mean_land_share_to_staples_among_cc_model_epsilon_grid= ones(no_epsilons);
urban_rural_inc_ratio_model_epsilon_grid= ones(no_epsilons);
urban_rural_wealth_ratio_model_epsilon_grid= ones(no_epsilons);
urban_rural_consumption_ratio_model_epsilon_grid= ones(no_epsilons);
p90_wealth_rural_epsilon_grid= ones(no_epsilons);
p90_wealth_urban_epsilon_grid= ones(no_epsilons);
p99_wealth_rural_epsilon_grid= ones(no_epsilons);
p99_wealth_urban_epsilon_grid= ones(no_epsilons);
p90_cons_tmp_epsilon_grid= ones(no_epsilons);
p90_income_tmp_epsilon_grid= ones(no_epsilons);
p99_cons_tmp_epsilon_grid= ones(no_epsilons);
p99_income_tmp_epsilon_grid= ones(no_epsilons);
staple_productivity_epsilon_grid= ones(no_epsilons);
cashcrop_productivity_epsilon_grid= ones(no_epsilons);
manuf_productivity_epsilon_grid= ones(no_epsilons);
relative_land_to_staples_epsilon_grid= ones(no_epsilons);
relative_land_to_cashcrop_epsilon_grid= ones(no_epsilons);
share_constrained_cashcrop_epsilon_grid= ones(no_epsilons);
var_MPX_staples_S_epsilon_grid= ones(no_epsilons);
var_MPX_cashcrop_B_epsilon_grid= ones(no_epsilons);
var_MPX_cashcrop_S_epsilon_grid= ones(no_epsilons);
share_constrained_staple_epsilon_grid= ones(no_epsilons);
APG_epsilon_grid= ones(no_epsilons);
urban_rural_consumption_ratio_model_real_epsilon_grid= ones(no_epsilons);
aggregate_consumption_epsilon_grid= ones(no_epsilons);
worker_pop_effective_epsilon_grid= ones(no_epsilons);
prod_manuf_epsilon_grid= ones(no_epsilons);
total_entry_cost_epsilon_grid= ones(no_epsilons);
prod_staple_epsilon_grid= ones(no_epsilons);
prod_cashcrop_epsilon_grid= ones(no_epsilons);
input_staple_epsilon_grid= ones(no_epsilons);
input_cashcrop_epsilon_grid= ones(no_epsilons);
total_maintenance_cost_epsilon_grid= ones(no_epsilons);
current_account_residual_epsilon_grid= ones(no_epsilons);
fraction_cashcrop_suboptimal_model_epsilon_grid = ones(no_epsilons);
welfare_epsilon_grid= ones(no_epsilons);
V_saved_epsilon_grid = ones(Baseline_parameter.ns_fine, 3, no_epsilons);
avg_labor_prod_rural_epsilon_grid= ones(no_epsilons);
avg_labor_prod_urban_epsilon_grid= ones(no_epsilons);
avg_agri_prod_rural_epsilon_grid= ones(no_epsilons);
avg_agri_prod_urban_epsilon_grid= ones(no_epsilons);
var_MPX_cashcrop_epsilon_grid= ones(no_epsilons);
var_MPX_epsilon_grid= ones(no_epsilons);
TFP_epsilon_grid= ones(no_epsilons);
YL_manuf_epsilon_grid= ones(no_epsilons);
YL_agr_epsilon_grid= ones(no_epsilons);
coeff_var_labor_prod_rural_epsilon_grid= ones(no_epsilons);
coeff_var_labor_prod_urban_epsilon_grid= ones(no_epsilons);
coeff_var_agri_prod_rural_epsilon_grid= ones(no_epsilons);
coeff_var_agri_prod_urban_epsilon_grid= ones(no_epsilons);
p90_wealth_epsilon_grid= ones(no_epsilons);
p99_wealth_epsilon_grid= ones(no_epsilons);
p90_cons_epsilon_grid_rural= ones(no_epsilons);
p99_cons_epsilon_grid_rural= ones(no_epsilons);
p90_cons_epsilon_grid_urban= ones(no_epsilons);
p99_cons_epsilon_grid_urban= ones(no_epsilons);
p90_income_epsilon_grid_rural= ones(no_epsilons);
p99_income_epsilon_grid_rural= ones(no_epsilons);
p90_income_epsilon_grid_urban= ones(no_epsilons);
p99_income_epsilon_grid_urban= ones(no_epsilons);
wealth_of_workers_epsilon_grid= ones(no_epsilons);
wealth_of_staples_epsilon_grid= ones(no_epsilons);
wealth_of_cashcrop_epsilon_grid= ones(no_epsilons);
c_B_worker_sum_epsilon_grid= ones(no_epsilons);
c_B_staple_sum_epsilon_grid= ones(no_epsilons);
c_B_cashcrop_sum_epsilon_grid= ones(no_epsilons);
c_S_worker_sum_epsilon_grid= ones(no_epsilons);
c_S_staple_sum_epsilon_grid = ones(no_epsilons);
c_S_cashcrop_sum_epsilon_grid= ones(no_epsilons);
transaction_cost_staple_sum_epsilon_grid= ones(no_epsilons);
transaction_cost_cashcrop_sum_epsilon_grid= ones(no_epsilons);
transaction_cost_worker_sum_epsilon_grid= ones(no_epsilons);
c_M_worker_sum_epsilon_grid= ones(no_epsilons);
c_M_staple_sum_epsilon_grid= ones(no_epsilons);
c_M_cashcrop_sum_epsilon_grid= ones(no_epsilons);
MPX_mean_log_epsilon_grid= ones(no_epsilons);
MPX_mean_staples_S_log_epsilon_grid= ones(no_epsilons);
MPX_mean_cashcrop_log_epsilon_grid= ones(no_epsilons);
APland_mean_log_epsilon_grid= ones(no_epsilons);
APland_mean_cashcrop_log_epsilon_grid= ones(no_epsilons);
APland_mean_staples_S_log_epsilon_grid= ones(no_epsilons);
var_APland_epsilon_grid= ones(no_epsilons);
var_APland_cashcrop_epsilon_grid= ones(no_epsilons);
var_APland_staples_S_epsilon_grid= ones(no_epsilons);
c_S_W_fine_epsilon_grid= ones(Baseline_parameter.ns_fine,3,no_epsilons);
c_B_W_fine_epsilon_grid= ones(Baseline_parameter.ns_fine,3,no_epsilons);
c_M_W_fine_epsilon_grid= ones(Baseline_parameter.ns_fine,3,no_epsilons);
c_S_S_fine_epsilon_grid= ones(Baseline_parameter.ns_fine,3,no_epsilons);
c_B_S_fine_epsilon_grid= ones(Baseline_parameter.ns_fine,3,no_epsilons);
c_M_S_fine_epsilon_grid= ones(Baseline_parameter.ns_fine,3,no_epsilons);
c_S_B_fine_epsilon_grid= ones(Baseline_parameter.ns_fine,3,no_epsilons);
c_B_B_fine_epsilon_grid= ones(Baseline_parameter.ns_fine,3,no_epsilons);
c_M_B_fine_epsilon_grid= ones(Baseline_parameter.ns_fine,3,no_epsilons);

for iterate_loc = 1:no_epsilons
    prices_tmp = prices_epsilon_grid[:,iterate_loc]
    epsilon_u_tmp = epsilon_grid[iterate_loc]
    epsilon_r_tmp = epsilon_grid[iterate_loc]
    println("Populate for " ,epsilon_grid[iterate_loc], " elasticity of productivity to undernourishment")
    (residual_goods_epsilon_grid[:,iterate_loc], stat_distr_epsilon_grid[:,iterate_loc], cons_fine_local_epsilon_grid[:,:,iterate_loc], a_prime_fine_local_epsilon_grid[:,:,iterate_loc],
    future_occupation_fine_local_epsilon_grid[:,:,iterate_loc],x_S_S_fine_epsilon_grid[:,:,iterate_loc],x_SC_fine_epsilon_grid[:,:,iterate_loc],x_BC_fine_epsilon_grid[:,:,iterate_loc], 
    coeff_epsilon_grid[:,:,iterate_loc],  transaction_cost_loss_epsilon_grid[iterate_loc],nominal_GDP_epsilon_grid[iterate_loc],welfare_val_epsilon_grid[:,iterate_loc],
    Import_value_epsilon_grid[iterate_loc],Export_value_epsilon_grid[iterate_loc],current_worker_pop_epsilon_grid[iterate_loc],current_staple_pop_epsilon_grid[iterate_loc],
    current_cashcrop_pop_epsilon_grid[iterate_loc], marketable_agr_surplus_share_epsilon_grid[iterate_loc],exportshare_cashcrop_epsilon_grid[iterate_loc],
    fraction_model_epsilon_grid[iterate_loc],program_spending_epsilon_grid[iterate_loc],prod_value_improvement_epsilon_grid[iterate_loc],share_selling_increase_epsilon_grid[iterate_loc],
    exp_ratio_model_epsilon_grid[iterate_loc],mig_rate_model_epsilon_grid[iterate_loc],rural_pop_only_staples_model_epsilon_grid[iterate_loc],
    rural_pop_only_cashcrop_model_epsilon_grid[iterate_loc],mean_land_share_to_staples_among_cc_model_epsilon_grid[iterate_loc],urban_rural_inc_ratio_model_epsilon_grid[iterate_loc],
    urban_rural_wealth_ratio_model_epsilon_grid[iterate_loc],urban_rural_consumption_ratio_model_epsilon_grid[iterate_loc],p90_wealth_rural_epsilon_grid[iterate_loc],
    p90_wealth_urban_epsilon_grid[iterate_loc],p99_wealth_rural_epsilon_grid[iterate_loc],p99_wealth_urban_epsilon_grid[iterate_loc],p90_cons_tmp_epsilon_grid[iterate_loc],
    p90_income_tmp_epsilon_grid[iterate_loc],p99_cons_tmp_epsilon_grid[iterate_loc],p99_income_tmp_epsilon_grid[iterate_loc],staple_productivity_epsilon_grid[iterate_loc],
    cashcrop_productivity_epsilon_grid[iterate_loc],manuf_productivity_epsilon_grid[iterate_loc],relative_land_to_staples_epsilon_grid[iterate_loc],
    relative_land_to_cashcrop_epsilon_grid[iterate_loc],share_constrained_cashcrop_epsilon_grid[iterate_loc],var_MPX_staples_S_epsilon_grid[iterate_loc],
    var_MPX_cashcrop_B_epsilon_grid[iterate_loc],var_MPX_cashcrop_S_epsilon_grid[iterate_loc],share_constrained_staple_epsilon_grid[iterate_loc],APG_epsilon_grid[iterate_loc],
    urban_rural_consumption_ratio_model_real_epsilon_grid[iterate_loc],aggregate_consumption_epsilon_grid[iterate_loc],
    worker_pop_effective_epsilon_grid[iterate_loc],prod_manuf_epsilon_grid[iterate_loc],total_entry_cost_epsilon_grid[iterate_loc],prod_staple_epsilon_grid[iterate_loc],prod_cashcrop_epsilon_grid[iterate_loc],input_staple_epsilon_grid[iterate_loc],input_cashcrop_epsilon_grid[iterate_loc],
    total_maintenance_cost_epsilon_grid[iterate_loc],current_account_residual_epsilon_grid[iterate_loc],fraction_cashcrop_suboptimal_model_epsilon_grid[iterate_loc],V_saved_epsilon_grid[:,:,iterate_loc],
    avg_labor_prod_rural_epsilon_grid[iterate_loc],avg_labor_prod_urban_epsilon_grid[iterate_loc],avg_agri_prod_rural_epsilon_grid[iterate_loc],avg_agri_prod_urban_epsilon_grid[iterate_loc],
    var_MPX_cashcrop_epsilon_grid[iterate_loc] ,var_MPX_epsilon_grid[iterate_loc] ,TFP_epsilon_grid[iterate_loc] ,YL_manuf_epsilon_grid[iterate_loc]  , 
    YL_agr_epsilon_grid[iterate_loc] ,coeff_var_labor_prod_rural_epsilon_grid[iterate_loc] ,coeff_var_labor_prod_urban_epsilon_grid[iterate_loc] ,
    coeff_var_agri_prod_rural_epsilon_grid[iterate_loc] , coeff_var_agri_prod_urban_epsilon_grid[iterate_loc] ,
    p90_wealth_epsilon_grid[iterate_loc] ,p99_wealth_epsilon_grid[iterate_loc] ,p90_cons_epsilon_grid_rural[iterate_loc] ,p99_cons_epsilon_grid_rural[iterate_loc] ,p90_cons_epsilon_grid_urban[iterate_loc],p99_cons_epsilon_grid_urban[iterate_loc],p90_income_epsilon_grid_rural[iterate_loc] ,
    p99_income_epsilon_grid_rural[iterate_loc],p90_income_epsilon_grid_urban[iterate_loc],p99_income_epsilon_grid_urban[iterate_loc],
    wealth_of_workers_epsilon_grid[iterate_loc] ,wealth_of_staples_epsilon_grid[iterate_loc] ,wealth_of_cashcrop_epsilon_grid[iterate_loc] ,
    c_B_worker_sum_epsilon_grid[iterate_loc],c_B_staple_sum_epsilon_grid[iterate_loc],c_B_cashcrop_sum_epsilon_grid[iterate_loc],c_S_worker_sum_epsilon_grid[iterate_loc],c_S_staple_sum_epsilon_grid[iterate_loc] ,c_S_cashcrop_sum_epsilon_grid[iterate_loc],
    transaction_cost_staple_sum_epsilon_grid[iterate_loc],transaction_cost_cashcrop_sum_epsilon_grid[iterate_loc],transaction_cost_worker_sum_epsilon_grid[iterate_loc],c_M_worker_sum_epsilon_grid[iterate_loc],c_M_staple_sum_epsilon_grid[iterate_loc],c_M_cashcrop_sum_epsilon_grid[iterate_loc]
    ,MPX_mean_log_epsilon_grid[iterate_loc], MPX_mean_staples_S_log_epsilon_grid[iterate_loc],MPX_mean_cashcrop_log_epsilon_grid[iterate_loc]
    , APland_mean_log_epsilon_grid[iterate_loc],APland_mean_cashcrop_log_epsilon_grid[iterate_loc], APland_mean_staples_S_log_epsilon_grid[iterate_loc],var_APland_epsilon_grid[iterate_loc],var_APland_cashcrop_epsilon_grid[iterate_loc],var_APland_staples_S_epsilon_grid[iterate_loc],
    c_S_W_fine_epsilon_grid[:,:,iterate_loc],c_B_W_fine_epsilon_grid[:,:,iterate_loc],c_M_W_fine_epsilon_grid[:,:,iterate_loc],c_S_S_fine_epsilon_grid[:,:,iterate_loc],c_B_S_fine_epsilon_grid[:,:,iterate_loc],c_M_S_fine_epsilon_grid[:,:,iterate_loc],c_S_B_fine_epsilon_grid[:,:,iterate_loc],c_B_B_fine_epsilon_grid[:,:,iterate_loc],c_M_B_fine_epsilon_grid[:,:,iterate_loc]) = details_model_externality(prices_tmp,No_subsidy_parameter,2,moments,balanced_share,foreign_supply_capital_subsidy_b,cons_level_substinence, epsilon_u_tmp, epsilon_r_tmp);
end

V_saved_epsilon_grid_reshaped=reshape(V_saved_epsilon_grid,Baseline_parameter.ns_fine*3,no_epsilons)
welfare_epsilon_grid_real = zeros(no_epsilons);
for iterate_loc = 1:no_epsilons
    welfare_epsilon_grid_real[iterate_loc] = sum(stat_distr_epsilon_grid[:, iterate_loc] .* (exp.((V_saved_subsidy_b_reshaped - V_saved_epsilon_grid_reshaped[:, iterate_loc]) * (1.0 - Baseline_parameter.β)))) - 1 # With balanced budget
end

welfare_epsilon_grid_real_alt = copy(welfare_epsilon_grid_real);
for iterate_loc = 1:no_epsilons
    welfare_epsilon_grid_real_alt[iterate_loc] = (1 - Baseline_parameter.β) * (sum(stat_distr_subsidy_b .* V_saved_subsidy_b_reshaped) - sum(stat_distr_epsilon_grid[:, iterate_loc] .* V_saved_epsilon_grid_reshaped[:, iterate_loc]) ) # With balanced budget
end

#Plotting of epsilonSgrid values
welfare_trans_epsilon_trans_real_corrected = copy(welfare_trans_epsilon_trans_real)
deleteat!(welfare_trans_epsilon_trans_real_corrected, 3)
epsilon_grid_corrected = copy(epsilon_grid)
deleteat!(epsilon_grid_corrected, 3)
welfare_trans_epsilon_trans_real_corrected[1] = welfare_trans_optimal_real[1];
welfare_epsilon_grid_real_corrected = copy(welfare_epsilon_grid_real);
deleteat!(welfare_epsilon_grid_real_corrected, 3)
# Unfortunately, epsilon[3] required to be replaced, the error is just too big 

epsilonSgrid = a_grid_fine_gen_midpoints(epsilon_grid_corrected,0.0,1.5,4,size(epsilon_grid_corrected))
welfare_real_grid = a_grid_fine_gen_midpoints(welfare_trans_epsilon_trans_real_corrected,welfare_trans_epsilon_trans_real_corrected[1],welfare_trans_epsilon_trans_real_corrected[end],4,size(epsilon_grid_corrected))
welfare_real_grid_lr = a_grid_fine_gen_midpoints(welfare_epsilon_grid_real_corrected,welfare_epsilon_grid_real_corrected[1],welfare_epsilon_grid_real_corrected[end],4,size(epsilon_grid_corrected))
#welfare_real_alt_grid = a_grid_fine_gen_midpoints(welfare_trans_epsilon_trans_real_alt,welfare_trans_epsilon_trans_real_alt[1],welfare_trans_epsilon_trans_real_alt[end],4,size(epsilon_grid))
# indifference: welfare_real_grid[13] - 0.5 
epsilonSgrid[9]


plot(epsilonSgrid,[100*welfare_real_grid,100*welfare_real_grid_lr],legends = true,
linewidth = 10,linestyle = [:solid :dot], xlabel = "",legend=(0.42,-0.15),
ylabel = "Welfare gains compared to no subsidy",marker = [:none :none],linecolor =[:blue :green], 
label=["Transition" "Long-run" ],
grid = true,xlims = [0.0,0.75],ylims = [-10.0,3.5],size = (800,800),
xticks = ([0.0,0.1155,0.3,0.5,0.7],["0.0","","0.3","0.5","0.7"]),
tickfontsize = 16,xguidefontsize=16,yguidefontsize=16,legendfontsize=16,fontfamily="times",
bottom_margin = 37*Plots.mm)

#plot!([0.1155,0.23],xticks = ([0.1155,0.23],["Strauss (1986)","γ*"]), xtickfontcolor=[:red,:black], seriestype="vline",linewidth = 6,linestyle = [:dot :dashdot],marker = [:none :circle],tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
annotate!(0.06, -10.43, text("Strauss (1986)", :left,"times",16,:red))
#annotate!(0.22, -1.06, text("γ*", :left,"times",14,:red))
annotate!(0.16, -11.00, text("Elasticity of productivity to undernourishment", :left,"times",16))
plot!([0.1155], seriestype="vline",color =[:red],line = :dash,label =false)
plot!([0.5], seriestype="vline",color =[:black],line = :dash,label =false)
savefig("Figure7b.svg")

using Colors
theme(:dark)
colors0 = [:black,:black,:red, :orange, :black, :black, :black]

function multicolor_xticks!(colors0)
    p = Plots.current()
    xticks, xlabels = Plots.xticks(p)[1]
    yl = Plots.ylims(p)
    y0 = @. zero(xticks) + yl[1] - 0.06*(yl[2] - yl[1])
    n = length(xticks)
    colors = colors0[1:n]
    xticks!(xticks, fill(" ",n))
    [annotate!(xi,yi,text(li,9,ci,:bottom)) for (xi,yi,li,ci) in zip(xticks, y0, xlabels, colors)]
    return Plots.current()
end

# useage
plot(1:20, randn(20), xlabel="x", ylabel="y")
multicolor_xticks!(colors0)


# indifference: welfare_real_alt_grid[7] - 0.18275
epsilonSgrid[7]


plot(epsilonSgrid,100*welfare_real_alt_grid,legends = false,
linewidth = 2,linestyle = [:solid :dot :dashdot], xlabel = "Elasticity of productivity",
ylabel = "Welfare gains compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),yticks = ([-10,-5,0,5,10,15,20,25 ],["-10","-5","0","5","10","15","20","25"]),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure7b.svg")

