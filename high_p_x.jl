include("analysis.jl")
##########################################################################################################################################################
##
##          ADDITIONAL: effects of doubling of fertilizer price
##
##########################################################################################################################################################

#    Subsidy equilibrium with balanced budget + double fertilizer price
prices_high_p_x_subsidy_b =  [ 1.5374346868502475,
0.22598290665699955,
0.2352290656580293]; #!
high_p_x_parameter_subsidy_b = copy(Baseline_parameter);
high_p_x_parameter_subsidy_b.p_x = 2*Baseline_parameter.p_x;
(residual_goods_high_p_x_subsidy_b, stat_distr_high_p_x_subsidy_b, cons_fine_local_high_p_x_subsidy_b, a_prime_fine_local_high_p_x_subsidy_b,future_occupation_fine_local_high_p_x_subsidy_b,x_S_S_fine_high_p_x_subsidy_b,x_SC_fine_high_p_x_subsidy_b,x_BC_fine_high_p_x_subsidy_b, coeff_high_p_x_subsidy_b,
transaction_cost_loss_high_p_x_subsidy_b,nominal_GDP_high_p_x_subsidy_b,welfare_val_high_p_x_subsidy_b,Import_value_high_p_x_subsidy_b,Export_value_high_p_x_subsidy_b,current_worker_pop_high_p_x_subsidy_b,current_staple_pop_high_p_x_subsidy_b,current_cashcrop_pop_high_p_x_subsidy_b,
marketable_agr_surplus_share_high_p_x_subsidy_b,exportshare_cashcrop_high_p_x_subsidy_b,
fraction_model_high_p_x_subsidy_b,program_spending_high_p_x_subsidy_b,prod_value_improvement_high_p_x_subsidy_b,share_selling_increase_high_p_x_subsidy_b,exp_ratio_model_high_p_x_subsidy_b,mig_rate_model_high_p_x_subsidy_b,rural_pop_only_staples_model_high_p_x_subsidy_b,rural_pop_only_cashcrop_model_high_p_x_subsidy_b,
mean_land_share_to_staples_among_cc_model_high_p_x_subsidy_b,urban_rural_inc_ratio_model_high_p_x_subsidy_b,urban_rural_wealth_ratio_model_high_p_x_subsidy_b,urban_rural_consumption_ratio_model_high_p_x_subsidy_b,
p90_wealth_rural_high_p_x_subsidy_b,p90_wealth_urban_high_p_x_subsidy_b,p99_wealth_rural_high_p_x_subsidy_b,p99_wealth_urban_high_p_x_subsidy_b,p90_cons_tmp_high_p_x_subsidy_b,p90_income_tmp_high_p_x_subsidy_b,p99_cons_tmp_high_p_x_subsidy_b,p99_income_tmp_high_p_x_subsidy_b,
staple_productivity_high_p_x_subsidy_b,cashcrop_productivity_high_p_x_subsidy_b,manuf_productivity_high_p_x_subsidy_b,relative_land_to_staples_high_p_x_subsidy_b,relative_land_to_cashcrop_high_p_x_subsidy_b,share_constrained_cashcrop_high_p_x_subsidy_b,var_MPX_staples_S_high_p_x_subsidy_b,var_MPX_cashcrop_B_high_p_x_subsidy_b,var_MPX_cashcrop_S_high_p_x_subsidy_b,
share_constrained_staple_high_p_x_subsidy_b,APG_high_p_x_subsidy_b,urban_rural_consumption_ratio_model_real_high_p_x_subsidy_b,aggregate_consumption_high_p_x_subsidy_b,
worker_pop_effective_high_p_x_subsidy_b,prod_manuf_high_p_x_subsidy_b,total_entry_cost_high_p_x_subsidy_b,prod_staple_high_p_x_subsidy_b,prod_cashcrop_high_p_x_subsidy_b,input_staple_high_p_x_subsidy_b,input_cashcrop_high_p_x_subsidy_b,
total_maintenance_cost_high_p_x_subsidy_b,current_account_residual_high_p_x_subsidy_b,fraction_cashcrop_suboptimal_model_high_p_x_subsidy_b,V_saved_high_p_x_subsidy_b
,avg_labor_prod_rural_high_p_x_subsidy_b,avg_labor_prod_urban_high_p_x_subsidy_b,avg_agri_prod_rural_high_p_x_subsidy_b,avg_agri_prod_urban_high_p_x_subsidy_b,
var_MPX_cashcrop_high_p_x_subsidy_b ,var_MPX_high_p_x_subsidy_b ,TFP_high_p_x_subsidy_b ,YL_manuf_high_p_x_subsidy_b  , 
YL_agr_high_p_x_subsidy_b ,coeff_var_labor_prod_rural_high_p_x_subsidy_b ,coeff_var_labor_prod_urban_high_p_x_subsidy_b ,
coeff_var_agri_prod_rural_high_p_x_subsidy_b , coeff_var_agri_prod_urban_high_p_x_subsidy_b ,
p90_wealth_high_p_x_subsidy_b,p99_wealth_high_p_x_subsidy_b,p90_cons_high_p_x_subsidy_b_rural,p99_cons_high_p_x_subsidy_b_rural,p90_cons_high_p_x_subsidy_b_urban,p99_cons_high_p_x_subsidy_b_urban,p90_income_high_p_x_subsidy_b_rural,
p99_income_high_p_x_subsidy_b_rural,p90_income_high_p_x_subsidy_b_urban,p99_income_high_p_x_subsidy_b_urban,
wealth_of_workers_high_p_x_subsidy_b,wealth_of_staples_high_p_x_subsidy_b,wealth_of_cashcrop_high_p_x_subsidy_b,
c_B_worker_sum_high_p_x_subsidy_b,c_B_staple_sum_high_p_x_subsidy_b,c_B_cashcrop_sum_high_p_x_subsidy_b,c_S_worker_sum_high_p_x_subsidy_b,c_S_staple_sum_high_p_x_subsidy_b ,c_S_cashcrop_sum_high_p_x_subsidy_b,
transaction_cost_staple_sum_high_p_x_subsidy_b,transaction_cost_cashcrop_sum_high_p_x_subsidy_b,transaction_cost_worker_sum_high_p_x_subsidy_b,c_M_worker_sum_high_p_x_subsidy_b,c_M_staple_sum_high_p_x_subsidy_b,c_M_cashcrop_sum_high_p_x_subsidy_b
,MPX_mean_log_high_p_x_subsidy_b, MPX_mean_staples_S_log_high_p_x_subsidy_b,MPX_mean_cashcrop_log_high_p_x_subsidy_b
, APland_mean_log_high_p_x_subsidy_b,APland_mean_cashcrop_log_high_p_x_subsidy_b, APland_mean_staples_S_log_high_p_x_subsidy_b,var_APland_high_p_x_subsidy_b,var_APland_cashcrop_high_p_x_subsidy_b,var_APland_staples_S_high_p_x_subsidy_b,
c_S_W_fine_high_p_x_subsidy_b,c_B_W_fine_high_p_x_subsidy_b,c_M_W_fine_high_p_x_subsidy_b,c_S_S_fine_high_p_x_subsidy_b,c_B_S_fine_high_p_x_subsidy_b,c_M_S_fine_high_p_x_subsidy_b,c_S_B_fine_high_p_x_subsidy_b,c_B_B_fine_high_p_x_subsidy_b,c_M_B_fine_high_p_x_subsidy_b) = details_model(prices_high_p_x_subsidy_b,high_p_x_parameter_subsidy_b,2,moments,1.0,foreign_supply_capital_subsidy_b);

#    Subsidy equilibrium with balanced budget + double fertilizer price and closed current account
prices_high_p_x_subsidy_b_CA0 =  [ 1.5374346868502475,
0.22598290665699955,
0.2352290656580293]; #!
high_p_x_parameter_subsidy_b_CA0 = copy(Baseline_parameter);
high_p_x_parameter_subsidy_b_CA0.p_x = 2*Baseline_parameter.p_x;
foreign_supply_capital_high_p_x_subsidy_b_CA0 = 8.053580609954484;
(residual_goods_high_p_x_subsidy_b_CA0, stat_distr_high_p_x_subsidy_b_CA0, cons_fine_local_high_p_x_subsidy_b_CA0, a_prime_fine_local_high_p_x_subsidy_b_CA0,future_occupation_fine_local_high_p_x_subsidy_b_CA0,x_S_S_fine_high_p_x_subsidy_b_CA0,x_SC_fine_high_p_x_subsidy_b_CA0,x_BC_fine_high_p_x_subsidy_b_CA0, coeff_high_p_x_subsidy_b_CA0,
transaction_cost_loss_high_p_x_subsidy_b_CA0,nominal_GDP_high_p_x_subsidy_b_CA0,welfare_val_high_p_x_subsidy_b_CA0,Import_value_high_p_x_subsidy_b_CA0,Export_value_high_p_x_subsidy_b_CA0,current_worker_pop_high_p_x_subsidy_b_CA0,current_staple_pop_high_p_x_subsidy_b_CA0,current_cashcrop_pop_high_p_x_subsidy_b_CA0,
marketable_agr_surplus_share_high_p_x_subsidy_b_CA0,exportshare_cashcrop_high_p_x_subsidy_b_CA0,
fraction_model_high_p_x_subsidy_b_CA0,program_spending_high_p_x_subsidy_b_CA0,prod_value_improvement_high_p_x_subsidy_b_CA0,share_selling_increase_high_p_x_subsidy_b_CA0,exp_ratio_model_high_p_x_subsidy_b_CA0,mig_rate_model_high_p_x_subsidy_b_CA0,rural_pop_only_staples_model_high_p_x_subsidy_b_CA0,rural_pop_only_cashcrop_model_high_p_x_subsidy_b_CA0,
mean_land_share_to_staples_among_cc_model_high_p_x_subsidy_b_CA0,urban_rural_inc_ratio_model_high_p_x_subsidy_b_CA0,urban_rural_wealth_ratio_model_high_p_x_subsidy_b_CA0,urban_rural_consumption_ratio_model_high_p_x_subsidy_b_CA0,
p90_wealth_rural_high_p_x_subsidy_b_CA0,p90_wealth_urban_high_p_x_subsidy_b_CA0,p99_wealth_rural_high_p_x_subsidy_b_CA0,p99_wealth_urban_high_p_x_subsidy_b_CA0,p90_cons_tmp_high_p_x_subsidy_b_CA0,p90_income_tmp_high_p_x_subsidy_b_CA0,p99_cons_tmp_high_p_x_subsidy_b_CA0,p99_income_tmp_high_p_x_subsidy_b_CA0,
staple_productivity_high_p_x_subsidy_b_CA0,cashcrop_productivity_high_p_x_subsidy_b_CA0,manuf_productivity_high_p_x_subsidy_b_CA0,relative_land_to_staples_high_p_x_subsidy_b_CA0,relative_land_to_cashcrop_high_p_x_subsidy_b_CA0,share_constrained_cashcrop_high_p_x_subsidy_b_CA0,var_MPX_staples_S_high_p_x_subsidy_b_CA0,var_MPX_cashcrop_B_high_p_x_subsidy_b_CA0,var_MPX_cashcrop_S_high_p_x_subsidy_b_CA0,
share_constrained_staple_high_p_x_subsidy_b_CA0,APG_high_p_x_subsidy_b_CA0,urban_rural_consumption_ratio_model_real_high_p_x_subsidy_b_CA0,aggregate_consumption_high_p_x_subsidy_b_CA0,
worker_pop_effective_high_p_x_subsidy_b_CA0,prod_manuf_high_p_x_subsidy_b_CA0,total_entry_cost_high_p_x_subsidy_b_CA0,prod_staple_high_p_x_subsidy_b_CA0,prod_cashcrop_high_p_x_subsidy_b_CA0,input_staple_high_p_x_subsidy_b_CA0,input_cashcrop_high_p_x_subsidy_b_CA0,
total_maintenance_cost_high_p_x_subsidy_b_CA0,current_account_residual_high_p_x_subsidy_b_CA0,fraction_cashcrop_suboptimal_model_high_p_x_subsidy_b_CA0,V_saved_high_p_x_subsidy_b_CA0
,avg_labor_prod_rural_high_p_x_subsidy_b_CA0,avg_labor_prod_urban_high_p_x_subsidy_b_CA0,avg_agri_prod_rural_high_p_x_subsidy_b_CA0,avg_agri_prod_urban_high_p_x_subsidy_b_CA0,
var_MPX_cashcrop_high_p_x_subsidy_b_CA0 ,var_MPX_high_p_x_subsidy_b_CA0 ,TFP_high_p_x_subsidy_b_CA0 ,YL_manuf_high_p_x_subsidy_b_CA0  , 
YL_agr_high_p_x_subsidy_b_CA0 ,coeff_var_labor_prod_rural_high_p_x_subsidy_b_CA0 ,coeff_var_labor_prod_urban_high_p_x_subsidy_b_CA0 ,
coeff_var_agri_prod_rural_high_p_x_subsidy_b_CA0 , coeff_var_agri_prod_urban_high_p_x_subsidy_b_CA0 ,
p90_wealth_high_p_x_subsidy_b_CA0,p99_wealth_high_p_x_subsidy_b_CA0,p90_cons_high_p_x_subsidy_b_CA0_rural,p99_cons_high_p_x_subsidy_b_CA0_rural,p90_cons_high_p_x_subsidy_b_CA0_urban,p99_cons_high_p_x_subsidy_b_CA0_urban,p90_income_high_p_x_subsidy_b_CA0_rural,
p99_income_high_p_x_subsidy_b_CA0_rural,p90_income_high_p_x_subsidy_b_CA0_urban,p99_income_high_p_x_subsidy_b_CA0_urban,
wealth_of_workers_high_p_x_subsidy_b_CA0,wealth_of_staples_high_p_x_subsidy_b_CA0,wealth_of_cashcrop_high_p_x_subsidy_b_CA0,
c_B_worker_sum_high_p_x_subsidy_b_CA0,c_B_staple_sum_high_p_x_subsidy_b_CA0,c_B_cashcrop_sum_high_p_x_subsidy_b_CA0,c_S_worker_sum_high_p_x_subsidy_b_CA0,c_S_staple_sum_high_p_x_subsidy_b_CA0 ,c_S_cashcrop_sum_high_p_x_subsidy_b_CA0,
transaction_cost_staple_sum_high_p_x_subsidy_b_CA0,transaction_cost_cashcrop_sum_high_p_x_subsidy_b_CA0,transaction_cost_worker_sum_high_p_x_subsidy_b_CA0,c_M_worker_sum_high_p_x_subsidy_b_CA0,c_M_staple_sum_high_p_x_subsidy_b_CA0,c_M_cashcrop_sum_high_p_x_subsidy_b_CA0
,MPX_mean_log_high_p_x_subsidy_b_CA0, MPX_mean_staples_S_log_high_p_x_subsidy_b_CA0,MPX_mean_cashcrop_log_high_p_x_subsidy_b_CA0
, APland_mean_log_high_p_x_subsidy_b_CA0,APland_mean_cashcrop_log_high_p_x_subsidy_b_CA0, APland_mean_staples_S_log_high_p_x_subsidy_b_CA0,var_APland_high_p_x_subsidy_b_CA0,var_APland_cashcrop_high_p_x_subsidy_b_CA0,var_APland_staples_S_high_p_x_subsidy_b_CA0,
c_S_W_fine_high_p_x_subsidy_b_CA0,c_B_W_fine_high_p_x_subsidy_b_CA0,c_M_W_fine_high_p_x_subsidy_b_CA0,c_S_S_fine_high_p_x_subsidy_b_CA0,c_B_S_fine_high_p_x_subsidy_b_CA0,c_M_S_fine_high_p_x_subsidy_b_CA0,c_S_B_fine_high_p_x_subsidy_b_CA0,c_B_B_fine_high_p_x_subsidy_b_CA0,c_M_B_fine_high_p_x_subsidy_b_CA0) = details_model(prices_high_p_x_subsidy_b_CA0,high_p_x_parameter_subsidy_b_CA0,2,moments,1.0,foreign_supply_capital_high_p_x_subsidy_b_CA0);


#   No subsidy equilibrium with double fertilizer price
prices_high_p_x_no_subsidy =  [ 1.1936833649221996,
0.17237129935387707]; #!
high_p_x_no_subsidy_parameter = copy(Baseline_parameter);
high_p_x_no_subsidy_parameter.p_x = 2*Baseline_parameter.p_x;
high_p_x_no_subsidy_parameter.τ_S= -0.0;
(residual_goods_high_p_x_no_subsidy, stat_distr_high_p_x_no_subsidy, cons_fine_local_high_p_x_no_subsidy, a_prime_fine_local_high_p_x_no_subsidy,future_occupation_fine_local_high_p_x_no_subsidy,x_S_S_fine_high_p_x_no_subsidy,x_SC_fine_high_p_x_no_subsidy,x_BC_fine_high_p_x_no_subsidy, coeff_high_p_x_no_subsidy,
transaction_cost_loss_high_p_x_no_subsidy,nominal_GDP_high_p_x_no_subsidy,welfare_val_high_p_x_no_subsidy,Import_value_high_p_x_no_subsidy,Export_value_high_p_x_no_subsidy,current_worker_pop_high_p_x_no_subsidy,current_staple_pop_high_p_x_no_subsidy,current_cashcrop_pop_high_p_x_no_subsidy,
marketable_agr_surplus_share_high_p_x_no_subsidy,exportshare_cashcrop_high_p_x_no_subsidy,
fraction_model_high_p_x_no_subsidy,program_spending_high_p_x_no_subsidy,prod_value_improvement_high_p_x_no_subsidy,share_selling_increase_high_p_x_no_subsidy,exp_ratio_model_high_p_x_no_subsidy,mig_rate_model_high_p_x_no_subsidy,rural_pop_only_staples_model_high_p_x_no_subsidy,rural_pop_only_cashcrop_model_high_p_x_no_subsidy,
mean_land_share_to_staples_among_cc_model_high_p_x_no_subsidy,urban_rural_inc_ratio_model_high_p_x_no_subsidy,urban_rural_wealth_ratio_model_high_p_x_no_subsidy,urban_rural_consumption_ratio_model_high_p_x_no_subsidy,
p90_wealth_rural_high_p_x_no_subsidy,p90_wealth_urban_high_p_x_no_subsidy,p99_wealth_rural_high_p_x_no_subsidy,p99_wealth_urban_high_p_x_no_subsidy,p90_cons_tmp_high_p_x_no_subsidy,p90_income_tmp_high_p_x_no_subsidy,p99_cons_tmp_high_p_x_no_subsidy,p99_income_tmp_high_p_x_no_subsidy,
staple_productivity_high_p_x_no_subsidy,cashcrop_productivity_high_p_x_no_subsidy,manuf_productivity_high_p_x_no_subsidy,relative_land_to_staples_high_p_x_no_subsidy,relative_land_to_cashcrop_high_p_x_no_subsidy,share_constrained_cashcrop_high_p_x_no_subsidy,var_MPX_staples_S_high_p_x_no_subsidy,var_MPX_cashcrop_B_high_p_x_no_subsidy,var_MPX_cashcrop_S_high_p_x_no_subsidy,
share_constrained_staple_high_p_x_no_subsidy,APG_high_p_x_no_subsidy,urban_rural_consumption_ratio_model_real_high_p_x_no_subsidy,aggregate_consumption_high_p_x_no_subsidy,
worker_pop_effective_high_p_x_no_subsidy,prod_manuf_high_p_x_no_subsidy,total_entry_cost_high_p_x_no_subsidy,prod_staple_high_p_x_no_subsidy,prod_cashcrop_high_p_x_no_subsidy,input_staple_high_p_x_no_subsidy,input_cashcrop_high_p_x_no_subsidy,
total_maintenance_cost_high_p_x_no_subsidy,current_account_residual_high_p_x_no_subsidy,fraction_cashcrop_suboptimal_model_high_p_x_no_subsidy,V_saved_high_p_x_no_subsidy
,avg_labor_prod_rural_high_p_x_no_subsidy,avg_labor_prod_urban_high_p_x_no_subsidy,avg_agri_prod_rural_high_p_x_no_subsidy,avg_agri_prod_urban_high_p_x_no_subsidy,
var_MPX_cashcrop_high_p_x_no_subsidy ,var_MPX_high_p_x_no_subsidy ,TFP_high_p_x_no_subsidy ,YL_manuf_high_p_x_no_subsidy  , 
YL_agr_high_p_x_no_subsidy ,coeff_var_labor_prod_rural_high_p_x_no_subsidy ,coeff_var_labor_prod_urban_high_p_x_no_subsidy ,
coeff_var_agri_prod_rural_high_p_x_no_subsidy , coeff_var_agri_prod_urban_high_p_x_no_subsidy ,
p90_wealth_high_p_x_no_subsidy,p99_wealth_high_p_x_no_subsidy,p90_cons_high_p_x_no_subsidy_rural,p99_cons_high_p_x_no_subsidy_rural,p90_cons_high_p_x_no_subsidy_urban,p99_cons_high_p_x_no_subsidy_urban,p90_income_high_p_x_no_subsidy_rural,
p99_income_high_p_x_no_subsidy_rural,p90_income_high_p_x_no_subsidy_urban,p99_income_high_p_x_no_subsidy_urban,
wealth_of_workers_high_p_x_no_subsidy,wealth_of_staples_high_p_x_no_subsidy,wealth_of_cashcrop_high_p_x_no_subsidy,
c_B_worker_sum_high_p_x_no_subsidy,c_B_staple_sum_high_p_x_no_subsidy,c_B_cashcrop_sum_high_p_x_no_subsidy,c_S_worker_sum_high_p_x_no_subsidy,c_S_staple_sum_high_p_x_no_subsidy ,c_S_cashcrop_sum_high_p_x_no_subsidy,
transaction_cost_staple_sum_high_p_x_no_subsidy,transaction_cost_cashcrop_sum_high_p_x_no_subsidy,transaction_cost_worker_sum_high_p_x_no_subsidy,c_M_worker_sum_high_p_x_no_subsidy,c_M_staple_sum_high_p_x_no_subsidy,c_M_cashcrop_sum_high_p_x_no_subsidy
,MPX_mean_log_high_p_x_no_subsidy, MPX_mean_staples_S_log_high_p_x_no_subsidy,MPX_mean_cashcrop_log_high_p_x_no_subsidy
, APland_mean_log_high_p_x_no_subsidy,APland_mean_cashcrop_log_high_p_x_no_subsidy, APland_mean_staples_S_log_high_p_x_no_subsidy,var_APland_high_p_x_no_subsidy,var_APland_cashcrop_high_p_x_no_subsidy,var_APland_staples_S_high_p_x_no_subsidy,
c_S_W_fine_high_p_x_no_subsidy,c_B_W_fine_high_p_x_no_subsidy,c_M_W_fine_high_p_x_no_subsidy,c_S_S_fine_high_p_x_no_subsidy,c_B_S_fine_high_p_x_no_subsidy,c_M_S_fine_high_p_x_no_subsidy,c_S_B_fine_high_p_x_no_subsidy,c_B_B_fine_high_p_x_no_subsidy,c_M_B_fine_high_p_x_no_subsidy) = details_model(prices_high_p_x_no_subsidy,high_p_x_no_subsidy_parameter,2,moments,0.0,foreign_supply_capital_subsidy_b);

#   No subsidy equilibrium with double fertilizer price and balanced current account
prices_high_p_x_no_subsidy_CA0 =  [  1.1987879728517163,
0.17139336979746067]; #!
high_p_x_no_subsidy_CA0_parameter = copy(Baseline_parameter);
high_p_x_no_subsidy_CA0_parameter.p_x = 2*Baseline_parameter.p_x;
high_p_x_no_subsidy_CA0_parameter.τ_S= -0.0;
foreign_supply_capital_high_p_x_no_subsidy_CA0 = 8.483430365528525;
(residual_goods_high_p_x_no_subsidy_CA0, stat_distr_high_p_x_no_subsidy_CA0, cons_fine_local_high_p_x_no_subsidy_CA0, a_prime_fine_local_high_p_x_no_subsidy_CA0,future_occupation_fine_local_high_p_x_no_subsidy_CA0,x_S_S_fine_high_p_x_no_subsidy_CA0,x_SC_fine_high_p_x_no_subsidy_CA0,x_BC_fine_high_p_x_no_subsidy_CA0, coeff_high_p_x_no_subsidy_CA0,
transaction_cost_loss_high_p_x_no_subsidy_CA0,nominal_GDP_high_p_x_no_subsidy_CA0,welfare_val_high_p_x_no_subsidy_CA0,Import_value_high_p_x_no_subsidy_CA0,Export_value_high_p_x_no_subsidy_CA0,current_worker_pop_high_p_x_no_subsidy_CA0,current_staple_pop_high_p_x_no_subsidy_CA0,current_cashcrop_pop_high_p_x_no_subsidy_CA0,
marketable_agr_surplus_share_high_p_x_no_subsidy_CA0,exportshare_cashcrop_high_p_x_no_subsidy_CA0,
fraction_model_high_p_x_no_subsidy_CA0,program_spending_high_p_x_no_subsidy_CA0,prod_value_improvement_high_p_x_no_subsidy_CA0,share_selling_increase_high_p_x_no_subsidy_CA0,exp_ratio_model_high_p_x_no_subsidy_CA0,mig_rate_model_high_p_x_no_subsidy_CA0,rural_pop_only_staples_model_high_p_x_no_subsidy_CA0,rural_pop_only_cashcrop_model_high_p_x_no_subsidy_CA0,
mean_land_share_to_staples_among_cc_model_high_p_x_no_subsidy_CA0,urban_rural_inc_ratio_model_high_p_x_no_subsidy_CA0,urban_rural_wealth_ratio_model_high_p_x_no_subsidy_CA0,urban_rural_consumption_ratio_model_high_p_x_no_subsidy_CA0,
p90_wealth_rural_high_p_x_no_subsidy_CA0,p90_wealth_urban_high_p_x_no_subsidy_CA0,p99_wealth_rural_high_p_x_no_subsidy_CA0,p99_wealth_urban_high_p_x_no_subsidy_CA0,p90_cons_tmp_high_p_x_no_subsidy_CA0,p90_income_tmp_high_p_x_no_subsidy_CA0,p99_cons_tmp_high_p_x_no_subsidy_CA0,p99_income_tmp_high_p_x_no_subsidy_CA0,
staple_productivity_high_p_x_no_subsidy_CA0,cashcrop_productivity_high_p_x_no_subsidy_CA0,manuf_productivity_high_p_x_no_subsidy_CA0,relative_land_to_staples_high_p_x_no_subsidy_CA0,relative_land_to_cashcrop_high_p_x_no_subsidy_CA0,share_constrained_cashcrop_high_p_x_no_subsidy_CA0,var_MPX_staples_S_high_p_x_no_subsidy_CA0,var_MPX_cashcrop_B_high_p_x_no_subsidy_CA0,var_MPX_cashcrop_S_high_p_x_no_subsidy_CA0,
share_constrained_staple_high_p_x_no_subsidy_CA0,APG_high_p_x_no_subsidy_CA0,urban_rural_consumption_ratio_model_real_high_p_x_no_subsidy_CA0,aggregate_consumption_high_p_x_no_subsidy_CA0,
worker_pop_effective_high_p_x_no_subsidy_CA0,prod_manuf_high_p_x_no_subsidy_CA0,total_entry_cost_high_p_x_no_subsidy_CA0,prod_staple_high_p_x_no_subsidy_CA0,prod_cashcrop_high_p_x_no_subsidy_CA0,input_staple_high_p_x_no_subsidy_CA0,input_cashcrop_high_p_x_no_subsidy_CA0,
total_maintenance_cost_high_p_x_no_subsidy_CA0,current_account_residual_high_p_x_no_subsidy_CA0,fraction_cashcrop_suboptimal_model_high_p_x_no_subsidy_CA0,V_saved_high_p_x_no_subsidy_CA0
,avg_labor_prod_rural_high_p_x_no_subsidy_CA0,avg_labor_prod_urban_high_p_x_no_subsidy_CA0,avg_agri_prod_rural_high_p_x_no_subsidy_CA0,avg_agri_prod_urban_high_p_x_no_subsidy_CA0,
var_MPX_cashcrop_high_p_x_no_subsidy_CA0 ,var_MPX_high_p_x_no_subsidy_CA0 ,TFP_high_p_x_no_subsidy_CA0 ,YL_manuf_high_p_x_no_subsidy_CA0  , 
YL_agr_high_p_x_no_subsidy_CA0 ,coeff_var_labor_prod_rural_high_p_x_no_subsidy_CA0 ,coeff_var_labor_prod_urban_high_p_x_no_subsidy_CA0 ,
coeff_var_agri_prod_rural_high_p_x_no_subsidy_CA0 , coeff_var_agri_prod_urban_high_p_x_no_subsidy_CA0 ,
p90_wealth_high_p_x_no_subsidy_CA0,p99_wealth_high_p_x_no_subsidy_CA0,p90_cons_high_p_x_no_subsidy_CA0_rural,p99_cons_high_p_x_no_subsidy_CA0_rural,p90_cons_high_p_x_no_subsidy_CA0_urban,p99_cons_high_p_x_no_subsidy_CA0_urban,p90_income_high_p_x_no_subsidy_CA0_rural,
p99_income_high_p_x_no_subsidy_CA0_rural,p90_income_high_p_x_no_subsidy_CA0_urban,p99_income_high_p_x_no_subsidy_CA0_urban,
wealth_of_workers_high_p_x_no_subsidy_CA0,wealth_of_staples_high_p_x_no_subsidy_CA0,wealth_of_cashcrop_high_p_x_no_subsidy_CA0,
c_B_worker_sum_high_p_x_no_subsidy_CA0,c_B_staple_sum_high_p_x_no_subsidy_CA0,c_B_cashcrop_sum_high_p_x_no_subsidy_CA0,c_S_worker_sum_high_p_x_no_subsidy_CA0,c_S_staple_sum_high_p_x_no_subsidy_CA0 ,c_S_cashcrop_sum_high_p_x_no_subsidy_CA0,
transaction_cost_staple_sum_high_p_x_no_subsidy_CA0,transaction_cost_cashcrop_sum_high_p_x_no_subsidy_CA0,transaction_cost_worker_sum_high_p_x_no_subsidy_CA0,c_M_worker_sum_high_p_x_no_subsidy_CA0,c_M_staple_sum_high_p_x_no_subsidy_CA0,c_M_cashcrop_sum_high_p_x_no_subsidy_CA0
,MPX_mean_log_high_p_x_no_subsidy_CA0, MPX_mean_staples_S_log_high_p_x_no_subsidy_CA0,MPX_mean_cashcrop_log_high_p_x_no_subsidy_CA0
, APland_mean_log_high_p_x_no_subsidy_CA0,APland_mean_cashcrop_log_high_p_x_no_subsidy_CA0, APland_mean_staples_S_log_high_p_x_no_subsidy_CA0,var_APland_high_p_x_no_subsidy_CA0,var_APland_cashcrop_high_p_x_no_subsidy_CA0,var_APland_staples_S_high_p_x_no_subsidy_CA0,
c_S_W_fine_high_p_x_no_subsidy_CA0,c_B_W_fine_high_p_x_no_subsidy_CA0,c_M_W_fine_high_p_x_no_subsidy_CA0,c_S_S_fine_high_p_x_no_subsidy_CA0,c_B_S_fine_high_p_x_no_subsidy_CA0,c_M_S_fine_high_p_x_no_subsidy_CA0,c_S_B_fine_high_p_x_no_subsidy_CA0,c_B_B_fine_high_p_x_no_subsidy_CA0,c_M_B_fine_high_p_x_no_subsidy_CA0) = details_model(prices_high_p_x_no_subsidy_CA0,high_p_x_no_subsidy_parameter,2,moments,0.0,foreign_supply_capital_high_p_x_no_subsidy_CA0);


##########################################################################################################################################################
##
##          ADDITIONAL analysis: effects of doubling of fertilizer price create additional variables
##
##########################################################################################################################################################



w_high_p_x_no_subsidy = prices_high_p_x_no_subsidy[2]*(1-Baseline_parameter.α)/(1+0)*(R/Baseline_parameter.α*1/prices_high_p_x_no_subsidy[2])^(Baseline_parameter.α/(Baseline_parameter.α-1));
w_high_p_x_subsidy_b = prices_high_p_x_subsidy_b[2]*(1-Baseline_parameter.α)/(1+prices_high_p_x_subsidy_b[3])*(R/Baseline_parameter.α*1/prices_high_p_x_subsidy_b[2])^(Baseline_parameter.α/(Baseline_parameter.α-1));

welfare_subsidy_b_high_p_x =sum(stat_distr_high_p_x_no_subsidy  .* (exp.((welfare_val_high_p_x_subsidy_b -welfare_val_high_p_x_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
welfare_subsidy_b_high_p_x_from_baseline =sum(stat_distr_subsidy_b  .* (exp.((welfare_val_high_p_x_subsidy_b -welfare_val_subsidy_b) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
welfare_no_subsidy_high_p_x_from_baseline =sum(stat_distr_no_subsidy  .* (exp.((welfare_val_high_p_x_no_subsidy -welfare_val_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget

welfare_subsidy_b_high_p_x_CA0 =sum(stat_distr_high_p_x_no_subsidy_CA0  .* (exp.((welfare_val_high_p_x_subsidy_b_CA0 -welfare_val_high_p_x_no_subsidy_CA0) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget

V_saved_high_p_x_subsidy_b_reshaped=reshape(V_saved_high_p_x_subsidy_b,Baseline_parameter.ns_fine*3)
V_saved_high_p_x_subsidy_b_CA0_reshaped=reshape(V_saved_high_p_x_subsidy_b_CA0,Baseline_parameter.ns_fine*3)
V_saved_high_p_x_no_subsidy_reshaped=reshape(V_saved_high_p_x_no_subsidy,Baseline_parameter.ns_fine*3)
V_saved_high_p_x_no_subsidy_CA0_reshaped=reshape(V_saved_high_p_x_no_subsidy_CA0,Baseline_parameter.ns_fine*3)

welfare_subsidy_b_high_p_x_real = sum(stat_distr_high_p_x_no_subsidy  .* (exp.((V_saved_high_p_x_subsidy_b_reshaped -V_saved_high_p_x_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
welfare_subsidy_b_high_p_x_CA0_real = sum(stat_distr_high_p_x_no_subsidy_CA0  .* (exp.((V_saved_high_p_x_subsidy_b_CA0_reshaped -V_saved_high_p_x_no_subsidy_CA0_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
welfare_high_p_x_CA0_real_from_subsidy = sum(stat_distr_subsidy_b  .* (exp.((V_saved_high_p_x_subsidy_b_CA0_reshaped -V_saved_subsidy_b_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
welfare_high_p_x_real_from_subsidy = sum(stat_distr_subsidy_b  .* (exp.((V_saved_high_p_x_subsidy_b_reshaped -V_saved_subsidy_b_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget

welfare_subsidy_b_high_p_x_real_alt =(1 - Baseline_parameter.β)* (sum(stat_distr_high_p_x_subsidy_b .* V_saved_high_p_x_subsidy_b_reshaped) - sum(stat_distr_high_p_x_no_subsidy.*V_saved_high_p_x_no_subsidy_reshaped)  ); # With balanced budget
welfare_high_p_x_CA0_real_alt_from_baseline =(1 - Baseline_parameter.β)* (sum(stat_distr_high_p_x_subsidy_b .* V_saved_high_p_x_subsidy_b_reshaped) - sum(stat_distr_subsidy_b.*V_saved_subsidy_b_reshaped)  ); # With balanced budget

welfare_high_p_x_real_from_baseline = sum(stat_distr_no_subsidy  .* (exp.((V_saved_high_p_x_no_subsidy_reshaped -V_saved_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
welfare_high_p_x_CA0_real_from_baseline = sum(stat_distr_no_subsidy  .* (exp.((V_saved_high_p_x_no_subsidy_CA0_reshaped -V_saved_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
welfare_high_p_x_real_to_subsidy_from_baseline = sum(stat_distr_no_subsidy  .* (exp.((V_saved_high_p_x_subsidy_b_reshaped-V_saved_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
welfare_high_p_x_CA0_real_to_subsidy_from_baseline = sum(stat_distr_no_subsidy  .* (exp.((V_saved_high_p_x_subsidy_b_CA0_reshaped -V_saved_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget

welfare_high_p_x_real_alt_from_baseline = (1 - Baseline_parameter.β)* (sum(stat_distr_high_p_x_no_subsidy.* V_saved_high_p_x_no_subsidy_reshaped) - sum(stat_distr_no_subsidy.*V_saved_no_subsidy_reshaped)  ); 
welfare_high_p_x_CA0_real_alt_from_baseline = (1 - Baseline_parameter.β)* (sum(stat_distr_high_p_x_no_subsidy_CA0.* V_saved_high_p_x_no_subsidy_CA0_reshaped) - sum(stat_distr_no_subsidy.*V_saved_no_subsidy_reshaped)  ); 
welfare_high_p_x_real_alt_to_subsidy_from_baseline = (1 - Baseline_parameter.β)* (sum(stat_distr_high_p_x_subsidy_b.* V_saved_high_p_x_subsidy_b_reshaped) - sum(stat_distr_no_subsidy.*V_saved_no_subsidy_reshaped)  ); 
welfare_high_p_x_CA0_real_alt_to_subsidy_from_baseline = (1 - Baseline_parameter.β)* (sum(stat_distr_high_p_x_subsidy_b_CA0.* V_saved_high_p_x_subsidy_b_CA0_reshaped) - sum(stat_distr_no_subsidy.*V_saved_no_subsidy_reshaped)  ); 

mean_land_share_staples_high_p_x_no_subsidy =  convert(Int64, round(100 * (mean_land_share_to_staples_among_cc_model_high_p_x_no_subsidy*current_cashcrop_pop_high_p_x_no_subsidy
 + current_staple_pop_high_p_x_no_subsidy)/(1 - current_worker_pop_high_p_x_no_subsidy ) ))
mean_land_share_staples_high_p_x_subsidy_b =  convert(Int64, round(100 * (mean_land_share_to_staples_among_cc_model_high_p_x_subsidy_b*current_cashcrop_pop_high_p_x_subsidy_b
 + current_staple_pop_high_p_x_subsidy_b)/(1 - current_worker_pop_high_p_x_subsidy_b ) ))
savings_no_subsidy = sum(stat_distr_no_subsidy .* reshape(a_prime_fine_local_no_subsidy,Baseline_parameter.ns_fine*3))
savings_subsidy_b = sum(stat_distr_subsidy_b .* reshape(a_prime_fine_local_subsidy_b,Baseline_parameter.ns_fine*3))
savings_high_p_x_no_subsidy = sum(stat_distr_high_p_x_no_subsidy .* reshape(a_prime_fine_local_high_p_x_no_subsidy,Baseline_parameter.ns_fine*3))
savings_high_p_x_subsidy_b = sum(stat_distr_high_p_x_subsidy_b .* reshape(a_prime_fine_local_high_p_x_subsidy_b,Baseline_parameter.ns_fine*3))

# Undernourishment requires additional calculations
# Staple consumption w a threshold value - assume the threshold is chosen such that the calibration hits the 47% of undernourished HHs
cons_level_substinence = 0.109;
    # Distribution of past occupations
worker_past_dist_high_p_x_subsidy_b = stat_distr_high_p_x_subsidy_b[(Baseline_parameter.ns_fine *0 + 1):(Baseline_parameter.ns_fine *1)];
staple_past_dist_high_p_x_subsidy_b = stat_distr_high_p_x_subsidy_b[(Baseline_parameter.ns_fine *1 + 1):(Baseline_parameter.ns_fine *2)];
cash_crop_past_dist_high_p_x_subsidy_b = stat_distr_high_p_x_subsidy_b[(Baseline_parameter.ns_fine *2 + 1):(Baseline_parameter.ns_fine *3)];

stay_workers_high_p_x_subsidy_b = worker_past_dist_high_p_x_subsidy_b.*(future_occupation_fine_local_high_p_x_subsidy_b[:,1].==1);
exit_staple_to_work_high_p_x_subsidy_b = staple_past_dist_high_p_x_subsidy_b.*(future_occupation_fine_local_high_p_x_subsidy_b[:,2].==1);
exit_cashcrop_to_work_high_p_x_subsidy_b = cash_crop_past_dist_high_p_x_subsidy_b.*(future_occupation_fine_local_high_p_x_subsidy_b[:,3].==1);
current_workers_high_p_x_subsidy_b = stay_workers_high_p_x_subsidy_b + exit_staple_to_work_high_p_x_subsidy_b + exit_cashcrop_to_work_high_p_x_subsidy_b;

entrants_staple_from_workers_high_p_x_subsidy_b = worker_past_dist_high_p_x_subsidy_b.*(future_occupation_fine_local_high_p_x_subsidy_b[:,1].==2);
incumbents_staple_high_p_x_subsidy_b = staple_past_dist_high_p_x_subsidy_b.*(future_occupation_fine_local_high_p_x_subsidy_b[:,2].==2);
exit_cashcrop_to_staple_high_p_x_subsidy_b = cash_crop_past_dist_high_p_x_subsidy_b.*(future_occupation_fine_local_high_p_x_subsidy_b[:,3].==2);
current_staple_high_p_x_subsidy_b = entrants_staple_from_workers_high_p_x_subsidy_b + incumbents_staple_high_p_x_subsidy_b + exit_cashcrop_to_staple_high_p_x_subsidy_b;

entrants_cashcrop_from_workers_high_p_x_subsidy_b = worker_past_dist_high_p_x_subsidy_b.*(future_occupation_fine_local_high_p_x_subsidy_b[:,1].==3);
entrants_from_staple_to_cashcrop_high_p_x_subsidy_b= staple_past_dist_high_p_x_subsidy_b.*(future_occupation_fine_local_high_p_x_subsidy_b[:,2].==3);
incumbents_cashcrop_high_p_x_subsidy_b = cash_crop_past_dist_high_p_x_subsidy_b.*(future_occupation_fine_local_high_p_x_subsidy_b[:,3].==3);
current_cashcrop_high_p_x_subsidy_b = entrants_cashcrop_from_workers_high_p_x_subsidy_b + incumbents_cashcrop_high_p_x_subsidy_b + entrants_from_staple_to_cashcrop_high_p_x_subsidy_b;

undernutritioned_workers_high_p_x_subsidy_b = sum( (c_S_W_fine_high_p_x_subsidy_b[:,1].<cons_level_substinence).*stay_workers_high_p_x_subsidy_b
 + (c_S_W_fine_high_p_x_subsidy_b[:,2].<cons_level_substinence).*exit_staple_to_work_high_p_x_subsidy_b +
(c_S_W_fine_high_p_x_subsidy_b[:,3].<cons_level_substinence) .*exit_cashcrop_to_work_high_p_x_subsidy_b );
undernutritioned_staple_farmer_high_p_x_subsidy_b = sum( (c_S_S_fine_high_p_x_subsidy_b[:,1].<cons_level_substinence).*entrants_staple_from_workers_high_p_x_subsidy_b +
(c_S_S_fine_high_p_x_subsidy_b[:,2].<cons_level_substinence).*incumbents_staple_high_p_x_subsidy_b +
(c_S_S_fine_high_p_x_subsidy_b[:,3].<cons_level_substinence) .*exit_cashcrop_to_staple_high_p_x_subsidy_b );
undernutritioned_cashcrop_farmer_high_p_x_subsidy_b = sum( (c_S_B_fine_high_p_x_subsidy_b[:,1].<cons_level_substinence).*entrants_cashcrop_from_workers_high_p_x_subsidy_b
 + (c_S_B_fine_high_p_x_subsidy_b[:,2].<cons_level_substinence).*entrants_from_staple_to_cashcrop_high_p_x_subsidy_b +
(c_S_B_fine_high_p_x_subsidy_b[:,3].<cons_level_substinence) .*incumbents_cashcrop_high_p_x_subsidy_b );
undernourished_high_p_x_subsidy_b = undernutritioned_workers_high_p_x_subsidy_b + undernutritioned_staple_farmer_high_p_x_subsidy_b +undernutritioned_cashcrop_farmer_high_p_x_subsidy_b

# No subsidy undernourished
worker_past_dist_high_p_x_no_subsidy = stat_distr_high_p_x_no_subsidy[(Baseline_parameter.ns_fine *0 + 1):(Baseline_parameter.ns_fine *1)];
staple_past_dist_high_p_x_no_subsidy = stat_distr_high_p_x_no_subsidy[(Baseline_parameter.ns_fine *1 + 1):(Baseline_parameter.ns_fine *2)];
cash_crop_past_dist_high_p_x_no_subsidy = stat_distr_high_p_x_no_subsidy[(Baseline_parameter.ns_fine *2 + 1):(Baseline_parameter.ns_fine *3)];

stay_workers_high_p_x_no_subsidy = worker_past_dist_high_p_x_no_subsidy.*(future_occupation_fine_local_high_p_x_no_subsidy[:,1].==1);
exit_staple_to_work_high_p_x_no_subsidy = staple_past_dist_high_p_x_no_subsidy.*(future_occupation_fine_local_high_p_x_no_subsidy[:,2].==1);
exit_cashcrop_to_work_high_p_x_no_subsidy = cash_crop_past_dist_high_p_x_no_subsidy.*(future_occupation_fine_local_high_p_x_no_subsidy[:,3].==1);
current_workers_high_p_x_no_subsidy = stay_workers_high_p_x_no_subsidy + exit_staple_to_work_high_p_x_no_subsidy + exit_cashcrop_to_work_high_p_x_no_subsidy;

entrants_staple_from_workers_high_p_x_no_subsidy = worker_past_dist_high_p_x_no_subsidy.*(future_occupation_fine_local_high_p_x_no_subsidy[:,1].==2);
incumbents_staple_high_p_x_no_subsidy = staple_past_dist_high_p_x_no_subsidy.*(future_occupation_fine_local_high_p_x_no_subsidy[:,2].==2);
exit_cashcrop_to_staple_high_p_x_no_subsidy = cash_crop_past_dist_high_p_x_no_subsidy.*(future_occupation_fine_local_high_p_x_no_subsidy[:,3].==2);
current_staple_high_p_x_no_subsidy = entrants_staple_from_workers_high_p_x_no_subsidy + incumbents_staple_high_p_x_no_subsidy + exit_cashcrop_to_staple_high_p_x_no_subsidy;

entrants_cashcrop_from_workers_high_p_x_no_subsidy = worker_past_dist_high_p_x_no_subsidy.*(future_occupation_fine_local_high_p_x_no_subsidy[:,1].==3);
entrants_from_staple_to_cashcrop_high_p_x_no_subsidy= staple_past_dist_high_p_x_no_subsidy.*(future_occupation_fine_local_high_p_x_no_subsidy[:,2].==3);
incumbents_cashcrop_high_p_x_no_subsidy = cash_crop_past_dist_high_p_x_no_subsidy.*(future_occupation_fine_local_high_p_x_no_subsidy[:,3].==3);
current_cashcrop_high_p_x_no_subsidy = entrants_cashcrop_from_workers_high_p_x_no_subsidy + incumbents_cashcrop_high_p_x_no_subsidy + entrants_from_staple_to_cashcrop_high_p_x_no_subsidy;

undernutritioned_workers_high_p_x_no_subsidy = sum( (c_S_W_fine_high_p_x_no_subsidy[:,1].<cons_level_substinence).*stay_workers_high_p_x_no_subsidy
    + (c_S_W_fine_high_p_x_no_subsidy[:,2].<cons_level_substinence).*exit_staple_to_work_high_p_x_no_subsidy +
(c_S_W_fine_high_p_x_no_subsidy[:,3].<cons_level_substinence) .*exit_cashcrop_to_work_high_p_x_no_subsidy );
undernutritioned_staple_farmer_high_p_x_no_subsidy = sum( (c_S_S_fine_high_p_x_no_subsidy[:,1].<cons_level_substinence).*entrants_staple_from_workers_high_p_x_no_subsidy +
(c_S_S_fine_high_p_x_no_subsidy[:,2].<cons_level_substinence).*incumbents_staple_high_p_x_no_subsidy +
(c_S_S_fine_high_p_x_no_subsidy[:,3].<cons_level_substinence) .*exit_cashcrop_to_staple_high_p_x_no_subsidy );
undernutritioned_cashcrop_farmer_high_p_x_no_subsidy = sum( (c_S_B_fine_high_p_x_no_subsidy[:,1].<cons_level_substinence).*entrants_cashcrop_from_workers_high_p_x_no_subsidy
    + (c_S_B_fine_high_p_x_no_subsidy[:,2].<cons_level_substinence).*entrants_from_staple_to_cashcrop_high_p_x_no_subsidy +
(c_S_B_fine_high_p_x_no_subsidy[:,3].<cons_level_substinence) .*incumbents_cashcrop_high_p_x_no_subsidy );

undernourished_high_p_x_no_subsidy = undernutritioned_workers_high_p_x_no_subsidy + undernutritioned_staple_farmer_high_p_x_no_subsidy +undernutritioned_cashcrop_farmer_high_p_x_no_subsidy;


# Main results table for the paper

println("\\hline & No FISP & FISP & No FISP & FISP \\\\ ")
println("& \$ p_x = ",  Baseline_parameter.p_x, "\$   &  \$ p_x = ",  Baseline_parameter.p_x, "\$ &  \$ p_x = ",  2*Baseline_parameter.p_x, "\$  & \$ p_x = ",  2*Baseline_parameter.p_x, "\$ \\\\ ")
println("\\hline 
\\\\[-1.8ex] \\textbf{Prices \\& Aggregates} & & &   \\\\[1ex]")
println("Cash crop, \$p_B\$ & ", round(prices_no_subsidy[1],digits = 1),"\\% & +"
,convert(Int64, round(100 * prices_subsidy_b[1]/prices_no_subsidy[1])) -100,"\\% &+"
,convert(Int64, round(100 * prices_high_p_x_no_subsidy[1]/prices_no_subsidy[1])) -100,"\\% &+"
,convert(Int64, round(100 * prices_high_p_x_subsidy_b[1]/prices_no_subsidy[1]))-100,
"\\%\\\\")
println("Manufacturing, \$p_M\$ &", round(prices_no_subsidy[2],digits = 1),
"\\%& +"
,convert(Int64, round(100 * prices_subsidy_b[2]/prices_no_subsidy[2])) -100,"\\% & +"
,convert(Int64, round(100 * prices_high_p_x_no_subsidy[2]/prices_no_subsidy[2]))-100,"\\%"
,"&+",convert(Int64, round(100 * prices_high_p_x_subsidy_b[2]/prices_no_subsidy[2]))-100,"\\%\\\\")
println("Wages, \$w\$ & ",round(w_no_subsidy,digits = 1)
,"& +",
convert(Int64, round(100 * w_subsidy_b/w_no_subsidy))-100,"\\% &+",
convert(Int64, round(100 * w_high_p_x_no_subsidy/w_no_subsidy))-100,"\\%  & +",
convert(Int64, round(100 * w_high_p_x_subsidy_b/w_no_subsidy))-100,"\\%\\\\")
println("Consumption & ",round(aggregate_consumption_no_subsidy,digits = 1), " & ",
convert(Int64, round(100 * aggregate_consumption_subsidy_b/aggregate_consumption_no_subsidy)) -100,"\\%  & +",
convert(Int64, round(100 * aggregate_consumption_high_p_x_no_subsidy/aggregate_consumption_no_subsidy)) -100,
"\\%  &", convert(Int64, round(100 * aggregate_consumption_high_p_x_subsidy_b/aggregate_consumption_no_subsidy)) -100,"\\% \\\\")
println("Nominal output & ",round(nominal_GDP_no_subsidy,digits = 1), " & ", convert(Int64, round(100 * nominal_GDP_subsidy_b/nominal_GDP_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * nominal_GDP_high_p_x_no_subsidy/nominal_GDP_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * nominal_GDP_high_p_x_subsidy_b/nominal_GDP_no_subsidy)) -100,"\\% \\\\")
println("Share of cash crop exported & ",convert(Int64,round(100 * exportshare_cashcrop_no_subsidy)), " & +", convert(Int64, round(100 * exportshare_cashcrop_subsidy_b/exportshare_cashcrop_no_subsidy)) -100,
"\\%  &",convert(Int64, round(100 * exportshare_cashcrop_high_p_x_no_subsidy/exportshare_cashcrop_no_subsidy)) -100,
"\\%  &-", convert(Int64, round(100 * exportshare_cashcrop_high_p_x_subsidy_b/exportshare_cashcrop_no_subsidy)) -100,"\\% \\\\")
println("Transaction cost & ",round(transaction_cost_loss_no_subsidy,digits = 1), " & ", convert(Int64, round(100 * transaction_cost_loss_subsidy_b/transaction_cost_loss_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * transaction_cost_loss_high_p_x_no_subsidy/transaction_cost_loss_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * transaction_cost_loss_high_p_x_subsidy_b/transaction_cost_loss_no_subsidy)) -100,"\\% \\\\")
println("Current account surplus \\% of GDP & ",convert(Int64,round(100*current_account_residual_no_subsidy,digits = 0)), "\\%&",convert(Int64,round(100*current_account_residual_subsidy_b,digits = 0)),
 "\\%&",convert(Int64,round(100*current_account_residual_high_p_x_no_subsidy,digits = 0)),
 "\\%&",convert(Int64,round(100*current_account_residual_high_p_x_subsidy_b,digits = 0)),"\\%  \\\\[1ex]")
println("\\textbf{Production}& & &   \\\\[1ex]")
#println("Share of staple-only farmers & ",convert(Int64,round(100 * rural_pop_only_staples_model_no_subsidy)), "\\% &+",convert(Int64,round(100 * rural_pop_only_staples_model_subsidy_b/rural_pop_only_staples_model_no_subsidy))-100,
#"\\% &+",convert(Int64,round(100 * rural_pop_only_staples_model_high_p_x_no_subsidy/rural_pop_only_staples_model_no_subsidy))-100,
#"\\% & +",convert(Int64,round(100 * rural_pop_only_staples_model_high_p_x_subsidy_b/rural_pop_only_staples_model_no_subsidy))-100,"\\% \\\\")
println("Staple production & ",round(prod_staple_no_subsidy,digits = 1)," & +", convert(Int64, round(100 * prod_staple_subsidy_b/prod_staple_no_subsidy)) -100 ,"\\%  &+",
convert(Int64, round(100 * prod_staple_high_p_x_no_subsidy/prod_staple_no_subsidy)) -100 ,"\\%  &+", convert(Int64, round(100 * prod_staple_high_p_x_subsidy_b/prod_staple_no_subsidy)) -100,"\\% \\\\")
println("Staple productivity &",round(staple_productivity_no_subsidy,digits = 1)," & +", convert(Int64, round(100 * staple_productivity_subsidy_b/staple_productivity_no_subsidy)) -100,"\\% &+",
convert(Int64, round(100 * staple_productivity_high_p_x_no_subsidy/staple_productivity_no_subsidy)) -100,"\\% &+", convert(Int64, round(100 * staple_productivity_high_p_x_subsidy_b/staple_productivity_no_subsidy)) -100,"\\% \\\\")
println("Staple input use & ",round(input_staple_no_subsidy,digits = 2)," & +", convert(Int64, round(100 * input_staple_subsidy_b/input_staple_no_subsidy)) -100 ,"\\%  &+",
convert(Int64, round(100 * input_staple_high_p_x_no_subsidy/input_staple_no_subsidy)) -100 ,"\\%  &+", convert(Int64, round(100 * input_staple_high_p_x_subsidy_b/input_staple_no_subsidy)) -100,"\\% \\\\")
println("Cash crop production & ",round(prod_cashcrop_no_subsidy,digits = 1)," & ", convert(Int64, round(100 * prod_cashcrop_subsidy_b/prod_cashcrop_no_subsidy)) -100,"\\%  & ",
convert(Int64, round(100 * prod_cashcrop_high_p_x_no_subsidy/prod_cashcrop_no_subsidy)) -100,"\\%  & +", convert(Int64, round(100 * prod_cashcrop_high_p_x_subsidy_b/prod_cashcrop_no_subsidy))-100 ,"\\% \\\\")
println("Cash crop productivity &",round(cashcrop_productivity_no_subsidy,digits = 1)," & +", convert(Int64, round(100 * cashcrop_productivity_subsidy_b/cashcrop_productivity_no_subsidy))-100 ,"\\% & +",
convert(Int64, round(100 * cashcrop_productivity_high_p_x_no_subsidy/cashcrop_productivity_no_subsidy))-100 ,"\\% & +", convert(Int64, round(100 * cashcrop_productivity_high_p_x_subsidy_b/cashcrop_productivity_no_subsidy))-100 ,"\\% \\\\")
println("Cash crop input use &",round(input_cashcrop_no_subsidy,digits = 2)," & +", convert(Int64, round(100 * input_cashcrop_subsidy_b/input_cashcrop_no_subsidy)) -100,"\\% &+",
convert(Int64, round(100 * input_cashcrop_high_p_x_no_subsidy/input_cashcrop_no_subsidy)) -100,"\\% &+", convert(Int64, round(100 * input_cashcrop_high_p_x_subsidy_b/input_cashcrop_no_subsidy)) -100,"\\% \\\\")
println("Savings &",round(savings_no_subsidy,digits = 2)," & +", convert(Int64, round(100 * savings_subsidy_b/savings_no_subsidy)) -100,"\\% &+",
convert(Int64, round(100 * savings_high_p_x_no_subsidy/savings_no_subsidy)) -100,"\\% &+", convert(Int64, round(100 * savings_high_p_x_subsidy_b/savings_no_subsidy)) -100,"\\% \\\\")
println("Share of land devoted to staples & ", mean_land_share_staples_no_subsidy, "\\% &+",convert(Int64,round(mean_land_share_staples_subsidy_b/mean_land_share_staples_no_subsidy*100,digits = 0) -100), "\\%  &+",
convert(Int64,round(mean_land_share_staples_high_p_x_no_subsidy/mean_land_share_staples_no_subsidy*100,digits = 0) -100), "\\%  &+",convert(Int64,round(mean_land_share_staples_high_p_x_subsidy_b/mean_land_share_staples_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Share of farmers without surplus & ", convert(Int64,round(100 * fraction_model_no_subsidy)), "\\% &",convert(Int64,round(fraction_model_subsidy_b/fraction_model_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(fraction_model_high_p_x_no_subsidy/fraction_model_no_subsidy*100,digits = 0) -100), "\\%  &",convert(Int64,round(fraction_model_high_p_x_subsidy_b/fraction_model_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")

println("Share of financially constrained farmers & ",convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_no_subsidy)), "\\% &+",convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_subsidy_b/fraction_cashcrop_suboptimal_model_no_subsidy))-100, "\\% &",
convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_high_p_x_no_subsidy/fraction_cashcrop_suboptimal_model_no_subsidy))-100, "\\% &",convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_high_p_x_subsidy_b/fraction_cashcrop_suboptimal_model_no_subsidy))-100,"\\% \\\\")
println("Manufacturing production & ",convert(Int64,round(prod_manuf_no_subsidy)), "   &",convert(Int64,round(100 * prod_manuf_subsidy_b/prod_manuf_no_subsidy)) -100, "\\% &+",
convert(Int64,round(100 * prod_manuf_high_p_x_no_subsidy/prod_manuf_no_subsidy)) -100, "\\% &+",convert(Int64,round(100 * prod_manuf_high_p_x_subsidy_b/prod_manuf_no_subsidy))-100,"\\% \\\\")
println("Urbanization rate & ",convert(Int64,round(100 * current_worker_pop_no_subsidy)), "\\% &",convert(Int64,round(100 * current_worker_pop_subsidy_b/current_worker_pop_no_subsidy)) -100, "\\% &",
convert(Int64,round(100 * current_worker_pop_high_p_x_no_subsidy/current_worker_pop_no_subsidy)) -100, "\\% &",convert(Int64,round(100 * current_worker_pop_high_p_x_subsidy_b/current_worker_pop_no_subsidy))-100,"\\% \\\\")
println("Agricultural productivity gap & ",round(APG_no_subsidy,digits = 1), "&+",convert(Int64,round(APG_subsidy_b/APG_no_subsidy*100,digits = 0) -100), "\\%  &+",
convert(Int64,round(APG_high_p_x_no_subsidy/APG_no_subsidy*100,digits = 0) -100), "\\%  &+",convert(Int64,round(APG_high_p_x_subsidy_b/APG_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Average agricultural ability & ",round(avg_agri_prod_rural_no_subsidy,digits = 1), "&",convert(Int64,round(avg_agri_prod_rural_subsidy_b/avg_agri_prod_rural_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(avg_agri_prod_rural_high_p_x_no_subsidy/avg_agri_prod_rural_no_subsidy*100,digits = 0) -100), "\\%  &",convert(Int64,round(avg_agri_prod_rural_high_p_x_subsidy_b/avg_agri_prod_rural_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Average worker ability & ",round(avg_labor_prod_urban_no_subsidy,digits = 1), "&+",convert(Int64,round(avg_labor_prod_urban_subsidy_b/avg_labor_prod_urban_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(avg_labor_prod_urban_high_p_x_no_subsidy/avg_labor_prod_urban_no_subsidy*100,digits = 0) -100), "\\%  &",convert(Int64,round(avg_labor_prod_urban_high_p_x_subsidy_b/avg_labor_prod_urban_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Dispersion in ARPX & ",round(var_MPX_no_subsidy,digits = 2), "&+",convert(Int64,round(var_MPX_subsidy_b/var_MPX_no_subsidy*100,digits = 0) -100), "\\%  &+",
convert(Int64,round(var_MPX_high_p_x_no_subsidy/var_MPX_no_subsidy*100 -100,digits = 0)), "\\%  &+",convert(Int64,round(var_MPX_high_p_x_subsidy_b/var_MPX_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Dispersion in ARPX for cash crop farmers & ",round(var_MPX_cashcrop_no_subsidy,digits = 2), "&",convert(Int64,round(var_MPX_cashcrop_subsidy_b/var_MPX_cashcrop_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(var_MPX_cashcrop_high_p_x_no_subsidy/var_MPX_cashcrop_no_subsidy*100,digits = 0) -100), "\\%  &",convert(Int64,round(var_MPX_cashcrop_high_p_x_subsidy_b/var_MPX_cashcrop_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Dispersion in ARPX for staple farmers & ",round(var_MPX_staples_S_no_subsidy,digits = 2), "&+",convert(Int64,round(var_MPX_staples_S_subsidy_b/var_MPX_staples_S_no_subsidy*100,digits = 0) -100), "\\%  &+",
convert(Int64,round(var_MPX_staples_S_high_p_x_no_subsidy/var_MPX_staples_S_no_subsidy*100,digits = 0) -100), "\\%  &+",convert(Int64,round(var_MPX_staples_S_high_p_x_subsidy_b/var_MPX_staples_S_no_subsidy*100,digits = 0) - 100),"\\% \\\\[1ex]")
println("\\textbf{Welfare and Inequality}& & & &   \\\\[1ex]")
println("Consumption equivalent welfare  & - &  +",round(100 * welfare_subsidy_b_real_alt,digits = 1), "\\% &+",
round(100 * welfare_high_p_x_real_alt_from_baseline,digits = 1), "\\% &", round(100 * welfare_high_p_x_real_alt_to_subsidy_from_baseline,digits = 1) ,"\\% \\\\ ")
println("Consumption equivalent welfare fixed distribution  & - &  +",round(100 * welfare_subsidy_b_real,digits = 1), "\\% &",
round(100 * welfare_high_p_x_real_from_baseline,digits = 1), "\\% &", round(100 * welfare_high_p_x_real_to_subsidy_from_baseline,digits = 1) ,"\\% \\\\ ")
println("Share of undernourished  & ", convert(Int64,round(100 * undernourished_no_subsidy)), "\\% &",convert(Int64,round(undernourished_subsidy_b/undernourished_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(undernourished_high_p_x_no_subsidy/undernourished_no_subsidy*100,digits = 0) -100), "\\%  &",convert(Int64,round(undernourished_high_p_x_subsidy_b/undernourished_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Avg urban-rural consumption ratio & ",round(urban_rural_consumption_ratio_model_no_subsidy,digits = 1), 
"&",convert(Int64,round(100* urban_rural_consumption_ratio_model_subsidy_b/urban_rural_consumption_ratio_model_no_subsidy-100)),
"\\%&+",convert(Int64,round(100* urban_rural_consumption_ratio_model_high_p_x_no_subsidy/urban_rural_consumption_ratio_model_no_subsidy-100)),
"\\%&+",convert(Int64,round(100* urban_rural_consumption_ratio_model_high_p_x_subsidy_b/urban_rural_consumption_ratio_model_no_subsidy-100)),"\\% \\\\")
println("Avg urban-rural income ratio & ",round(urban_rural_inc_ratio_model_no_subsidy,digits = 1), 
"&+",convert(Int64,round(100* urban_rural_inc_ratio_model_subsidy_b/urban_rural_inc_ratio_model_no_subsidy-100)),
"\\%&+",convert(Int64,round(100* urban_rural_inc_ratio_model_high_p_x_no_subsidy/urban_rural_inc_ratio_model_no_subsidy-100)),
"\\%&+",convert(Int64,round(100* urban_rural_inc_ratio_model_high_p_x_subsidy_b/urban_rural_inc_ratio_model_no_subsidy-100)),"\\% \\\\")
println("Avg urban-rural wealth ratio & ",round(urban_rural_wealth_ratio_model_no_subsidy,digits = 1), 
"&+",convert(Int64,round(100* urban_rural_wealth_ratio_model_subsidy_b/urban_rural_wealth_ratio_model_no_subsidy-100)),
"\\%&+",convert(Int64,round(100* urban_rural_wealth_ratio_model_high_p_x_no_subsidy/urban_rural_wealth_ratio_model_no_subsidy-100)),
"\\%&+",convert(Int64,round(100* urban_rural_wealth_ratio_model_high_p_x_subsidy_b/urban_rural_wealth_ratio_model_no_subsidy-100)),"\\% \\\\")
 println("Share of wealth of the top 10\$\\% \$   & ",convert(Int64, round(100 * p90_wealth_no_subsidy)), " & ", convert(Int64, round(100 * p90_wealth_subsidy_b/p90_wealth_no_subsidy)) -100,
 "\\%  &+", convert(Int64, round(100 * p90_wealth_high_p_x_no_subsidy/p90_wealth_no_subsidy)) -100,
 "\\%  &+", convert(Int64, round(100 * p90_wealth_high_p_x_subsidy_b/p90_wealth_no_subsidy)) -100,"\\% \\\\")
println("Share of income of the top 10\$\\% \$   & ",convert(Int64, round(100 * p90_income_tmp_no_subsidy)), " & +", convert(Int64, round(100 * p90_income_tmp_subsidy_b/p90_income_tmp_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * p90_income_tmp_high_p_x_no_subsidy/p90_income_tmp_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * p90_income_tmp_high_p_x_subsidy_b/p90_income_tmp_no_subsidy)) -100,"\\% \\\\")
println("Share of consumption of the top 10\$\\% \$   & ",convert(Int64, round(100 * p90_cons_tmp_no_subsidy)), " & ", convert(Int64, round(100 * p90_cons_tmp_subsidy_b/p90_cons_tmp_no_subsidy)) -100,
"\\%  &+",
convert(Int64, round(100 * p90_cons_tmp_high_p_x_no_subsidy/p90_cons_tmp_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * p90_cons_tmp_high_p_x_subsidy_b/p90_cons_tmp_no_subsidy)) -100,"\\% \\\\")



y_value = (exp.((V_saved_high_p_x_subsidy_b[:,1] -V_saved_high_p_x_no_subsidy[:,1]) * (1.0 - Baseline_parameter.β) ) ) .- 1
y_mat = 100* mat_creator(y_value);
plot([y_mat[1,1:5:end],y_mat[3,1:5:end],y_mat[13,1:5:end],y_mat[16,1:5:end]], label=["Low rural & urban \\theta" "High rural, low urban \\theta" "Low rural, high urban \\theta" "High rural & urban \\theta"]
,legend=:outerbottom ,
linewidth = 2,linestyle = [:solid :dash :dot :dashdot],ylims = [-20.0,35.0], xlabel = "Wealth",
ylabel = "Welfare change",xticks = ([5,14,18,21 ],["\$50","\$1000", "\$2500", "\$6000"]),yticks = ([-20,-10,0,10,20,30 ],["-20","-10","0","10","20","30"]),
grid = false,size = (800, 800),legendcolumns=2,marker = [:none :none :circle :none],
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure1a_high_p_x.svg")
#Farmers
y_value = (exp.((V_saved_high_p_x_subsidy_b[:,2] -V_saved_high_p_x_no_subsidy[:,2]) * (1.0 - Baseline_parameter.β) ) ) .- 1
y_mat = 100* mat_creator(y_value);
plot([y_mat[1,1:5:end],y_mat[3,1:5:end],y_mat[13,1:5:end],y_mat[16,1:5:end]], label=["Low rural & urban \\theta" "High rural, low urban \\theta" "Low rural, high urban \\theta" "High rural & urban \\theta"]
,legend=:outerbottom ,
linewidth = 2,linestyle = [:solid :dash :dot :dashdot],ylims = [-30.0,20.0], xlabel = "Wealth",
ylabel = "Welfare change",xticks = ([5,14,18,21 ],["\$50","\$1000", "\$2500", "\$6000"]),yticks = ([-30,-20,-10,0,10,20 ],["-30","-20","-10","0","10","20",]),
grid = false,size = (800, 800),legendcolumns=2,marker = [:none :none :circle :none],
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure1b_high_p_x.svg")