include("analysis.jl")
##########################################################################################################################################################
##
##         Plotting: 
##
##########################################################################################################################################################


#Figure 1 plot - lambda_S. requires to evaluate within the solve_model_calibration.jl with the following local definitions:
cons_span = 10:140
prices = copy(prices_subsidy_b);
parameters_tmp = copy(Baseline_parameter);
foreign_supply_capital = copy(foreign_supply_capital_subsidy_b);
out = 2; 

# Evaluate the inside of solve_model_calibration2.jl (apart from the first and last line):

#For lambda = 1, use q_S_c1
#For lambda = 1+Q_S, use q_S_c2 
# Construct q_S for the intermediate cases:

worker_shadow_price = (1 + parameters_tmp.Q_S) * ones(size(cons_span))
plot([coeff_λ_2_s_fine[cons_span,1],coeff_λ_2_cashcrop_residual_unconstrained_fine[cons_span,1],coeff_λ_2_cashcrop_residual_unconstrained_fine[cons_span,166],
    worker_shadow_price], dpi = 300, label=["Staple farmer" "Cash crop farmer" "Productive cash crop farmer" "Urban workers"],legend=:outerbottom,
     linewidth = 2,linestyle = [:solid :dash :dot :dashdot],linecolor = [:blue :red :purple :green],marker = [:none :none :circle :none],
     ylims = [1.0,3.0], xlims = [1.0,130.0],
     ylabel = "Shadow price of staples",automargins = true,xlabel = "Consumption bundle",
minorgrid = false,xticks =([30,49],["C1","C2"]),yticks =([1,3.0],[1,"1 + Q_S"]),xmirror = false,size = (800, 800),legendcolumns=2,
    tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure1a.svg")


# Create the staple production figure:

x_S_subsidized = ((1 + τ_S) * p_x./(coeff_λ_2_s_fine[cons_span,1] .* ζ .* ϕ_S .*θ_fine[1] )).^(1 / (ζ - 1));
q_S_subsidized = ϕ_S * x_S_subsidized.^ζ .*θ_fine[1];
x_S_notsubsidized = (p_x./(coeff_λ_2_s_fine[cons_span,1] .* ζ .* ϕ_S .*θ_fine[1] )).^(1 / (ζ - 1));
q_S_notsubsidized = ϕ_S * x_S_notsubsidized.^ζ .*θ_fine[1];

plot([q_S_subsidized,q_S_notsubsidized ], label=[ "With subsidy" "Without subsidy"],legend=:outerbottom,
     linewidth = 2,linecolor = [:blue],linestyle = [:solid :dash ],ylims = [0.02,cons_level_substinence],xlims = [1.0,130.0],
     ylabel = "Quantity of staples produced",automargins = true,xlabel = "Consumption bundle",
     minorgrid = false,xticks =([30,49],["C1","C2"]),yticks =([0.02,cons_level_substinence],["Subsistence","Undernourishement"]),size = (800, 800),
     tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")

savefig("Figure1b.svg")
# Calibration Tables - targeted moments

println("Labor augmenting TFP in manuf \$A_W\$ & ", Baseline_parameter.A_W, "& 20\\% urban population [LSMS2010] & 18\\% & ",convert(Int64, round(100 * current_worker_pop_subsidy_b)),"\$\\% \$ \\tabularnewline")
println("Cash crop export demand shifter \$a_D\$ & ",round(Baseline_parameter.a_D,digits =2), " & Share of 2010 exported cash crops in total prod. of cash crops [FAO 2010] & 73\$\\% \$ & ",convert(Int64, round(100 * exportshare_cashcrop_subsidy_b)),"\$\\% \$ \\tabularnewline")
println("Subsidy rate for staple inputs \$\\tau_S\$ & ",convert(Int64, round(-1 *100 * Baseline_parameter.τ_S)), "\\% & Aggregate cost of program (\\% GDP) \\cite{chirwadorwardbook}& 3\$\\% \$ & ",convert(Int64, round(100 * program_spending_subsidy_b)),"\$\\% \$ \\tabularnewline")
println("Urban entry cost \$F_{M}\$  & ",convert(Int64, round(Baseline_parameter.F_W)), " & Rural-urban migration rate \\cite{bicketalJME} & 1\$\\% \$ & ",convert(Int64, round(100 * mig_rate_model_subsidy_b)),"\$\\% \$ \\tabularnewline")
println("Cash crop maintenance cost \$FM_{B}\$  & ", round(Baseline_parameter.FM_B,digits = 1), "& Share of land devoted to staples [FAOStat] & 70 \\%  & ",convert(Int64, round(100 * (mean_land_share_to_staples_among_cc_model_subsidy_b*current_cashcrop_pop_subsidy_b
 + current_staple_pop_subsidy_b)/(1 - current_worker_pop_subsidy_b ) )),"\\%  \\tabularnewline" )
println("Working capital constraint \$\\kappa\$ & ",convert(Int64, round(100 * Baseline_parameter.κ)), "\\% & Share of cash crop farmers with suboptimal inputs \\cite{bruneetal2016} & 70\\% & ",convert(Int64, round(100 * fraction_cashcrop_suboptimal_model_subsidy_b)),"\$\\% \$ \\tabularnewline")
println("Returns to scale in land for farming \$\\rho\$  & ",round( Baseline_parameter.ρ,digits =2),"&Standard deviation of average product of farms \\cite{restuccia2017land} & 1.8& ",round((var_APland_subsidy_b).^(1/2), digits = 2)," \\tabularnewline")
println("Correlation of urban/rural shocks \$\\rho_{RU}\$  & ",round( Baseline_parameter.ρ_SW,digits =2), "& Agricultural productivity gap \\cite{gollin2014} & 6.5 & ", round(APG_subsidy_b, digits = 1)," \\tabularnewline")

# Change the returns to scale calibration line
# Calibration Tables -nontargeted moments

println("Agriculture output share in GDP [WB 2010] & 30\$\\% \$ & ",convert(Int64, round(100 * marketable_agr_surplus_share_subsidy_b)),"\$\\% \$ \\tabularnewline")
println("Production value improvement to cash grant \\cite{aggarwaletal} & 22\$\\% \$ & ",convert(Int64, round(100 * prod_value_improvement_subsidy_b)),"\$\\% \$ \\tabularnewline")
println("Share of land devoted to staples among cash crop farmers [LSMS2010] & 30\\%& ",convert(Int64, round(100 * mean_land_share_to_staples_among_cc_model_subsidy_b)),"\\% \\tabularnewline")
println("Standard deviation of average product of fertilizer [LSMS2010] & 0.76& ",round((var_MPX_subsidy_b).^(1/2), digits = 2)," \\tabularnewline")
println("Standard deviation of average product of fertilizer [LSMS2010] & 0.76& ",round((var_MPX_subsidy_b).^(1/2), digits = 2)," \\tabularnewline")
#println( "Aggregate cashcrop to staple fertilizer expenditure ratio [FAO 2010]   & 2.6 & ", round(exp_ratio_model_subsidy_b, digits = 1)," \\tabularnewline")

println("\\textbf{Inequality measures from \\cite{santaJDE}:}\\tabularnewline")
println("Urban-rural wealth ratio  & 3.0& ",round( urban_rural_wealth_ratio_model_subsidy_b, digits = 1)," \\tabularnewline")
println("Urban-rural income ratio & 2.4& ", round(urban_rural_inc_ratio_model_subsidy_b, digits = 1)," \\tabularnewline")
println("Urban-rural consumption ratio & 2.2& ",round( urban_rural_consumption_ratio_model_subsidy_b, digits = 1)," \\tabularnewline")
println("Share of wealth of the top 10\$\\% \$  & 58\\%& ",convert(Int64, round(100 * p90_wealth_subsidy_b)),"\\% \\tabularnewline")
println("Share of income of the top 10\$\\% \$  & 48\\%& ",convert(Int64, round(100 * p90_income_tmp_subsidy_b)),"\\% \\tabularnewline")
println("Share of consumption of the top 10\$\\% \$  & 34\\%& ",convert(Int64, round(100 * p90_cons_tmp_subsidy_b)),"\\% \\tabularnewline")
println("Share of wealth of the top 1\$\\% \$  & 25\\%& ",convert(Int64, round(100 * p99_wealth_subsidy_b)),"\\% \\tabularnewline")
println("Share of income of the top 1\$\\% \$  & 18\\%& ",convert(Int64, round(100 * p99_income_tmp_subsidy_b)),"\\% \\tabularnewline")
println("Share of consumption of the top 1\$\\% \$  & 8\\%& ",convert(Int64, round(100 * p99_cons_tmp_subsidy_b)),"\\% \\tabularnewline")

println("Rural inequality measures from De Magalhaes \\& Santaeulalia (2018)")
println("Share of wealth of the top 10\$\\% \$  [ De Magalhaes \\& Santaeulalia, 2018] & 58\\%& ",convert(Int64, round(100 * p90_wealth_rural_subsidy_b)),"\\% \\tabularnewline")
println("Share of income of the top 10\$\\% \$  [ De Magalhaes \\& Santaeulalia, 2018] & 48\\%& ",convert(Int64, round(100 * p90_income_subsidy_b_rural)),"\\% \\tabularnewline")
println("Share of consumption of the top 10\$\\% \$  [ De Magalhaes \\& Santaeulalia, 2018] & 34\\%& ",convert(Int64, round(100 * p90_cons_subsidy_b_rural)),"\\% \\tabularnewline")
println("Share of wealth of the top 1\$\\% \$  [ De Magalhaes \\& Santaeulalia, 2018] & 25\\%& ",convert(Int64, round(100 * p99_wealth_rural_subsidy_b)),"\\% \\tabularnewline")
println("Share of income of the top 1\$\\% \$  [ De Magalhaes \\& Santaeulalia, 2018] & 18\\%& ",convert(Int64, round(100 * p99_income_subsidy_b_rural)),"\\% \\tabularnewline")
println("Share of consumption of the top 1\$\\% \$  [ De Magalhaes \\& Santaeulalia, 2018] & 8\\%& ",convert(Int64, round(100 * p99_cons_subsidy_b_rural)),"\\% \\tabularnewline")

# Main results table for the presentations

println("\\hline & No FISP & FISP & FISP \\\\ ")
println("&   & aid-financed  & tax-financed \$ \\tau_W = ",  convert(Int64,round(100* prices_subsidy_b[3])), "\\%\$ \\\\ ")
println("\\hline \\\\[-1.8ex] ")
println("\\alert<1>{Consumption equivalent welfare}  & - &  ",round(100 * welfare_subsidy_nb,digits = 1), "\\% &", round(100 * welfare_subsidy_b,digits = 1) ,"\\% \\\\ ")
println("Prices: \$p_B/p_M/w\$ & ", round(prices_no_subsidy[1],digits = 1),"/",round(prices_no_subsidy[2],digits = 1),"/",round(w_no_subsidy,digits = 1)," & \$+",convert(Int64, round(100 * prices_subsidy_nb[1]/prices_no_subsidy[1])) -100,"\\%/"
,"+",convert(Int64, round(100 * prices_subsidy_nb[2]/prices_no_subsidy[2]))-100,"\\%/"
,"+",convert(Int64, round(100 * w_subsidy_nb/w_no_subsidy))-100,"\\%\$  & \$+",convert(Int64, round(100 * prices_subsidy_b[1]/prices_no_subsidy[1]))-100,"\\%/"
,"+",convert(Int64, round(100 * prices_subsidy_b[2]/prices_no_subsidy[2]))-100,"\\%/"
,"+",convert(Int64, round(100 * w_subsidy_b/w_no_subsidy))-100,"\\%\$\\\\")
#println("Prices: \$p_B/p_M/w\$ & ", round(prices_no_subsidy[1],digits = 1),"/",round(prices_no_subsidy[2],digits = 1),"&+",convert(Int64, round(100 * prices_subsidy_nb[1]/prices_no_subsidy[1])) -100,"\\%/"
#,"+",convert(Int64, round(100 * prices_subsidy_nb[2]/prices_no_subsidy[2]))-100,"\\% & +" ,convert(Int64, round(100 * prices_subsidy_b[1]/prices_no_subsidy[1]))-100,"\\%/"
#,"+",convert(Int64, round(100 * prices_subsidy_b[2]/prices_no_subsidy[2]))-100,"\\%\\\\")
println("Share of staple-only farmers & ",convert(Int64,round(100 * rural_pop_only_staples_model_no_subsidy)), "\\% &+",convert(Int64,round(100 * rural_pop_only_staples_model_subsidy_nb/rural_pop_only_staples_model_no_subsidy))-100, "\\% & +",convert(Int64,round(100 * rural_pop_only_staples_model_subsidy_b/rural_pop_only_staples_model_no_subsidy))-100,"\\% \\\\")
println("\\alert<2>{Staple production} & ",round(prod_staple_no_subsidy,digits = 1)," & +", convert(Int64, round(100 * prod_staple_subsidy_nb/prod_staple_no_subsidy)) -100 ,"\\%  &+", convert(Int64, round(100 * prod_staple_subsidy_b/prod_staple_no_subsidy)) -100,"\\% \\\\")
println("\\alert<2>{Staple productivity} &",round(staple_productivity_no_subsidy,digits = 1)," & +", convert(Int64, round(100 * staple_productivity_subsidy_nb/staple_productivity_no_subsidy)) -100,"\\% &+", convert(Int64, round(100 * staple_productivity_subsidy_b/staple_productivity_no_subsidy)) -100,"\\% \\\\")
println("\\alert<3>{Cash crop production} & ",round(prod_cashcrop_no_subsidy,digits = 1)," & ", convert(Int64, round(100 * prod_cashcrop_subsidy_nb/prod_cashcrop_no_subsidy)) -100,"\\%  & +", convert(Int64, round(100 * prod_cashcrop_subsidy_b/prod_cashcrop_no_subsidy))-100 ,"\\% \\\\")
println("\\alert<3>{Cash crop productivity} &",round(cashcrop_productivity_no_subsidy,digits = 1)," & +", convert(Int64, round(100 * cashcrop_productivity_subsidy_nb/cashcrop_productivity_no_subsidy))-100 ,"\\% & +", convert(Int64, round(100 * cashcrop_productivity_subsidy_b/cashcrop_productivity_no_subsidy))-100 ,"\\% \\\\")
#println("\\alert<4>{Share of financially constrained farmers} & ",convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_no_subsidy)), "\\% &",convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_subsidy_nb/fraction_cashcrop_suboptimal_model_no_subsidy))-100, "\\% &",convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_subsidy_b/fraction_cashcrop_suboptimal_model_no_subsidy))-100,"\\% \\\\")
println("\\alert<4>{Manufacturing production} & ",convert(Int64,round(prod_manuf_no_subsidy)), "\\% &+",convert(Int64,round(100 * prod_manuf_subsidy_nb/prod_manuf_no_subsidy)) -100, "\\% &+",convert(Int64,round(100 * prod_manuf_subsidy_b/prod_manuf_no_subsidy))-100,"\\% \\\\")
println("\\alert<4>{Urbanization rate} & ",convert(Int64,round(100 * current_worker_pop_no_subsidy)), "\\% &",convert(Int64,round(100 * current_worker_pop_subsidy_nb/current_worker_pop_no_subsidy)) -100, "\\% &",convert(Int64,round(100 * current_worker_pop_subsidy_b/current_worker_pop_no_subsidy))-100,"\\% \\\\")
println("\\alert<5>{Agricultural productivity gap} & ",round(APG_no_subsidy,digits = 1), "&+",convert(Int64,round(APG_subsidy_nb/APG_no_subsidy*100,digits = 0) -100), "\\%  &+",convert(Int64,round(APG_subsidy_b/APG_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("\\alert<5>{Avg urban-rural consumption rate} & ",round(urban_rural_consumption_ratio_model_no_subsidy,digits = 1), "&",round(urban_rural_consumption_ratio_model_subsidy_nb,digits = 1),
"&",round(urban_rural_consumption_ratio_model_subsidy_b,digits = 1)," \\\\")
println("\\alert<5>{Dispersion in ARPX} & ",round(var_MPX_no_subsidy,digits = 2), "&+",convert(Int64,round(var_MPX_subsidy_nb/var_MPX_no_subsidy*100,digits = 0) -100), "\\%  &+",convert(Int64,round(var_MPX_subsidy_b/var_MPX_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Consumption & ",round(aggregate_consumption_no_subsidy,digits = 1), " & +", convert(Int64, round(100 * aggregate_consumption_subsidy_nb/aggregate_consumption_no_subsidy)) -100,
"\\%  &", convert(Int64, round(100 * aggregate_consumption_subsidy_b/aggregate_consumption_no_subsidy)) -100,"\\% \\\\")
println("Transaction cost & ",round(transaction_cost_loss_no_subsidy,digits = 1), " & +", convert(Int64, round(100 * transaction_cost_loss_subsidy_nb/transaction_cost_loss_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * transaction_cost_loss_subsidy_b/transaction_cost_loss_no_subsidy)) -100,"\\% \\\\")
println("Current account surplus \\% of GDP & ",convert(Int64,round(100*current_account_residual_no_subsidy,digits = 0)), "\\%&",convert(Int64,round(100*current_account_residual_subsidy_nb,digits = 0)),
 "\\%&",convert(Int64,round(100*current_account_residual_subsidy_b,digits = 0)),"\\% \\\\")

 

 # Main results table for the paper

println("\\hline & No FISP & FISP & FISP & FISP \\\\ ")
println("&   & partial eqm. & aid-financed  & tax-financed \$ \\tau_W = ",  convert(Int64,round(100* prices_subsidy_b[3])), "\\%\$ \\\\ ")
println("\\hline 
\\\\[-1.8ex] \\textbf{Prices \\& Aggregates} & & &   \\\\[1ex]")
println("Cash crop, \$p_B\$ & ", round(prices_no_subsidy[1],digits = 1)," & 0\\% & +"
,convert(Int64, round(100 * prices_subsidy_nb[1]/prices_no_subsidy[1])) -100,"\\% &+"
,convert(Int64, round(100 * prices_subsidy_b[1]/prices_no_subsidy[1]))-100,
"\\%\\\\")
println("Manufacturing, \$p_M\$ &", round(prices_no_subsidy[2],digits = 1),
"& 0\\% & +",convert(Int64, round(100 * prices_subsidy_nb[2]/prices_no_subsidy[2]))-100,"\\%"
,"&+",convert(Int64, round(100 * prices_subsidy_b[2]/prices_no_subsidy[2]))-100,"\\%\\\\")
println("Wages, \$w\$ & ",round(w_no_subsidy,digits = 1)
,"& 0\\% &+",convert(Int64, round(100 * w_subsidy_nb/w_no_subsidy))-100,"\\%  & +",convert(Int64, round(100 * w_subsidy_b/w_no_subsidy))-100,"\\%\\\\")
println("Consumption & ",round(aggregate_consumption_no_subsidy,digits = 1), " & ",
convert(Int64, round(100 * aggregate_consumption_subsidy_partial/aggregate_consumption_no_subsidy)) -100,"\\%  & +",
convert(Int64, round(100 * aggregate_consumption_subsidy_nb/aggregate_consumption_no_subsidy)) -100,
"\\%  &", convert(Int64, round(100 * aggregate_consumption_subsidy_b/aggregate_consumption_no_subsidy)) -100,"\\% \\\\")
println("Nominal output & ",round(nominal_GDP_no_subsidy,digits = 1), " & ", convert(Int64, round(100 * nominal_GDP_subsidy_partial/nominal_GDP_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * nominal_GDP_subsidy_nb/nominal_GDP_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * nominal_GDP_subsidy_b/nominal_GDP_no_subsidy)) -100,"\\% \\\\")
println("Share of cash crop exported & ",convert(Int64,round(100 * exportshare_cashcrop_no_subsidy)), " & +", convert(Int64, round(100 * exportshare_cashcrop_subsidy_partial/exportshare_cashcrop_no_subsidy)) -100,
"\\%  &",convert(Int64, round(100 * exportshare_cashcrop_subsidy_nb/exportshare_cashcrop_no_subsidy)) -100,
"\\%  &-", convert(Int64, round(100 * exportshare_cashcrop_subsidy_b/exportshare_cashcrop_no_subsidy)) -100,"\\% \\\\")
println("Transaction cost & ",round(transaction_cost_loss_no_subsidy,digits = 1), " & ", convert(Int64, round(100 * transaction_cost_loss_subsidy_partial/transaction_cost_loss_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * transaction_cost_loss_subsidy_nb/transaction_cost_loss_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * transaction_cost_loss_subsidy_b/transaction_cost_loss_no_subsidy)) -100,"\\% \\\\")
println("Current account surplus \\% of GDP & ",convert(Int64,round(100*current_account_residual_no_subsidy,digits = 0)), "\\%&",convert(Int64,round(100*current_account_residual_subsidy_partial,digits = 0)),
 "\\%&",convert(Int64,round(100*current_account_residual_subsidy_nb,digits = 0)),
 "\\%&",convert(Int64,round(100*current_account_residual_subsidy_b,digits = 0)),"\\%  \\\\[1ex]")
println("\\textbf{Production}& & &   \\\\[1ex]")
#println("Share of staple-only farmers & ",convert(Int64,round(100 * rural_pop_only_staples_model_no_subsidy)), "\\% &+",convert(Int64,round(100 * rural_pop_only_staples_model_subsidy_partial/rural_pop_only_staples_model_no_subsidy))-100,
#"\\% &+",convert(Int64,round(100 * rural_pop_only_staples_model_subsidy_nb/rural_pop_only_staples_model_no_subsidy))-100,
#"\\% & +",convert(Int64,round(100 * rural_pop_only_staples_model_subsidy_b/rural_pop_only_staples_model_no_subsidy))-100,"\\% \\\\")
println("Staple production & ",round(prod_staple_no_subsidy,digits = 1)," & +", convert(Int64, round(100 * prod_staple_subsidy_partial/prod_staple_no_subsidy)) -100 ,"\\%  &+",
convert(Int64, round(100 * prod_staple_subsidy_nb/prod_staple_no_subsidy)) -100 ,"\\%  &+", convert(Int64, round(100 * prod_staple_subsidy_b/prod_staple_no_subsidy)) -100,"\\% \\\\")
println("Staple productivity &",round(staple_productivity_no_subsidy,digits = 1)," & +", convert(Int64, round(100 * staple_productivity_subsidy_partial/staple_productivity_no_subsidy)) -100,"\\% &+",
convert(Int64, round(100 * staple_productivity_subsidy_nb/staple_productivity_no_subsidy)) -100,"\\% &+", convert(Int64, round(100 * staple_productivity_subsidy_b/staple_productivity_no_subsidy)) -100,"\\% \\\\")
println("Cash crop production & ",round(prod_cashcrop_no_subsidy,digits = 1)," & ", convert(Int64, round(100 * prod_cashcrop_subsidy_partial/prod_cashcrop_no_subsidy)) -100,"\\%  & ",
convert(Int64, round(100 * prod_cashcrop_subsidy_nb/prod_cashcrop_no_subsidy)) -100,"\\%  & +", convert(Int64, round(100 * prod_cashcrop_subsidy_b/prod_cashcrop_no_subsidy))-100 ,"\\% \\\\")
println("Cash crop productivity &",round(cashcrop_productivity_no_subsidy,digits = 1)," & +", convert(Int64, round(100 * cashcrop_productivity_subsidy_partial/cashcrop_productivity_no_subsidy))-100 ,"\\% & +",
convert(Int64, round(100 * cashcrop_productivity_subsidy_nb/cashcrop_productivity_no_subsidy))-100 ,"\\% & +", convert(Int64, round(100 * cashcrop_productivity_subsidy_b/cashcrop_productivity_no_subsidy))-100 ,"\\% \\\\")
println("Share of land devoted to staples & ", mean_land_share_staples_no_subsidy, "\\% &+",convert(Int64,round(mean_land_share_staples_subsidy_partial/mean_land_share_staples_no_subsidy*100,digits = 0) -100), "\\%  &+",
convert(Int64,round(mean_land_share_staples_subsidy_nb/mean_land_share_staples_no_subsidy*100,digits = 0) -100), "\\%  &+",convert(Int64,round(mean_land_share_staples_subsidy_b/mean_land_share_staples_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Share of farmers without surplus & ", convert(Int64,round(100 * fraction_model_no_subsidy)), "\\% &",convert(Int64,round(fraction_model_subsidy_partial/fraction_model_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(fraction_model_subsidy_nb/fraction_model_no_subsidy*100,digits = 0) -100), "\\%  &",convert(Int64,round(fraction_model_subsidy_b/fraction_model_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")

println("Share of financially constrained farmers & ",convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_no_subsidy)), "\\% &+",convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_subsidy_partial/fraction_cashcrop_suboptimal_model_no_subsidy))-100, "\\% &",
convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_subsidy_nb/fraction_cashcrop_suboptimal_model_no_subsidy))-100, "\\% &",convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_subsidy_b/fraction_cashcrop_suboptimal_model_no_subsidy))-100,"\\% \\\\")
println("Manufacturing production & ",convert(Int64,round(prod_manuf_no_subsidy)), "   &",convert(Int64,round(100 * prod_manuf_subsidy_partial/prod_manuf_no_subsidy)) -100, "\\% &+",
convert(Int64,round(100 * prod_manuf_subsidy_nb/prod_manuf_no_subsidy)) -100, "\\% &+",convert(Int64,round(100 * prod_manuf_subsidy_b/prod_manuf_no_subsidy))-100,"\\% \\\\")
println("Urbanization rate & ",convert(Int64,round(100 * current_worker_pop_no_subsidy)), "\\% &",convert(Int64,round(100 * current_worker_pop_subsidy_partial/current_worker_pop_no_subsidy)) -100, "\\% &",
convert(Int64,round(100 * current_worker_pop_subsidy_nb/current_worker_pop_no_subsidy)) -100, "\\% &",convert(Int64,round(100 * current_worker_pop_subsidy_b/current_worker_pop_no_subsidy))-100,"\\% \\\\")
println("Agricultural productivity gap & ",round(APG_no_subsidy,digits = 1), "&+",convert(Int64,round(APG_subsidy_partial/APG_no_subsidy*100,digits = 0) -100), "\\%  &+",
convert(Int64,round(APG_subsidy_nb/APG_no_subsidy*100,digits = 0) -100), "\\%  &+",convert(Int64,round(APG_subsidy_b/APG_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Average agricultural ability & ",round(avg_agri_prod_rural_no_subsidy,digits = 1), "&",convert(Int64,round(avg_agri_prod_rural_subsidy_partial/avg_agri_prod_rural_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(avg_agri_prod_rural_subsidy_nb/avg_agri_prod_rural_no_subsidy*100,digits = 0) -100), "\\%  &",convert(Int64,round(avg_agri_prod_rural_subsidy_b/avg_agri_prod_rural_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Average worker ability & ",round(avg_labor_prod_urban_no_subsidy,digits = 1), "&+",convert(Int64,round(avg_labor_prod_urban_subsidy_partial/avg_labor_prod_urban_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(avg_labor_prod_urban_subsidy_nb/avg_labor_prod_urban_no_subsidy*100,digits = 0) -100), "\\%  &",convert(Int64,round(avg_labor_prod_urban_subsidy_b/avg_labor_prod_urban_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Dispersion in ARPX & ",round(var_MPX_no_subsidy,digits = 2), "&+",convert(Int64,round(var_MPX_subsidy_partial/var_MPX_no_subsidy*100,digits = 0) -100), "\\%  &+",
convert(Int64,round(var_MPX_subsidy_nb/var_MPX_no_subsidy*100 -100,digits = 0)), "\\%  &+",convert(Int64,round(var_MPX_subsidy_b/var_MPX_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Dispersion in ARPX for cash crop farmers & ",round(var_MPX_cashcrop_no_subsidy,digits = 2), "&",convert(Int64,round(var_MPX_cashcrop_subsidy_partial/var_MPX_cashcrop_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(var_MPX_cashcrop_subsidy_nb/var_MPX_cashcrop_no_subsidy*100,digits = 0) -100), "\\%  &",convert(Int64,round(var_MPX_cashcrop_subsidy_b/var_MPX_cashcrop_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Dispersion in ARPX for staple farmers & ",round(var_MPX_staples_S_no_subsidy,digits = 2), "&+",convert(Int64,round(var_MPX_staples_S_subsidy_partial/var_MPX_staples_S_no_subsidy*100,digits = 0) -100), "\\%  &+",
convert(Int64,round(var_MPX_staples_S_subsidy_nb/var_MPX_staples_S_no_subsidy*100,digits = 0) -100), "\\%  &+",convert(Int64,round(var_MPX_staples_S_subsidy_b/var_MPX_staples_S_no_subsidy*100,digits = 0) - 100),"\\% \\\\[1ex]")
println("\\textbf{Welfare and Inequality}& & & &   \\\\[1ex]")
println("Consumption equivalent welfare  & - &  +",round(100 * welfare_subsidy_partial_real,digits = 1), "\\% &+",
round(100 * welfare_transition_subsidy_nb,digits = 1), "\\% &", round(100 * welfare_transition_subsidy_b,digits = 1) ,"\\% \\\\ ")
println("Consumption equivalent welfare - urban  & - &  +",round(100 * welfare_subsidy_partial_real_workers_lr,digits = 1), "\\% &+",
round(100 * welfare_subsidy_nb_real_workers,digits = 1), "\\% &+", round(100 * welfare_subsidy_b_real_workers,digits = 1) ,"\\% \\\\ ")
println("Consumption equivalent welfare - rural  & - &  +",round(100 * welfare_subsidy_partial_real_farmers_lr,digits = 1), "\\% &+",
round(100 * welfare_subsidy_nb_real_farmers,digits = 1), "\\% &", round(100 * welfare_subsidy_b_real_farmers,digits = 1) ,"\\% \\\\ ")
println("Share of undernourished  & ", convert(Int64,round(100 * undernourished_no_subsidy)), "\\% &",convert(Int64,round(undernourished_subsidy_partial/undernourished_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(undernourished_subsidy_nb/undernourished_no_subsidy*100,digits = 0) -100), "\\%  &",convert(Int64,round(undernourished_subsidy_b/undernourished_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Avg urban-rural consumption ratio & ",round(urban_rural_consumption_ratio_model_no_subsidy,digits = 1), 
"&",convert(Int64,round(100* urban_rural_consumption_ratio_model_subsidy_partial/urban_rural_consumption_ratio_model_no_subsidy-100)),
"\\%&+",convert(Int64,round(100* urban_rural_consumption_ratio_model_subsidy_nb/urban_rural_consumption_ratio_model_no_subsidy-100)),
"\\%&+",convert(Int64,round(100* urban_rural_consumption_ratio_model_subsidy_b/urban_rural_consumption_ratio_model_no_subsidy-100)),"\\% \\\\")
println("Avg urban-rural income ratio & ",round(urban_rural_inc_ratio_model_no_subsidy,digits = 1), 
"&+",convert(Int64,round(100* urban_rural_inc_ratio_model_subsidy_partial/urban_rural_inc_ratio_model_no_subsidy-100)),
"\\%&+",convert(Int64,round(100* urban_rural_inc_ratio_model_subsidy_nb/urban_rural_inc_ratio_model_no_subsidy-100)),
"\\%&+",convert(Int64,round(100* urban_rural_inc_ratio_model_subsidy_b/urban_rural_inc_ratio_model_no_subsidy-100)),"\\% \\\\")
println("Avg urban-rural wealth ratio & ",round(urban_rural_wealth_ratio_model_no_subsidy,digits = 1), 
"&+",convert(Int64,round(100* urban_rural_wealth_ratio_model_subsidy_partial/urban_rural_wealth_ratio_model_no_subsidy-100)),
"\\%&+",convert(Int64,round(100* urban_rural_wealth_ratio_model_subsidy_nb/urban_rural_wealth_ratio_model_no_subsidy-100)),
"\\%&+",convert(Int64,round(100* urban_rural_wealth_ratio_model_subsidy_b/urban_rural_wealth_ratio_model_no_subsidy-100)),"\\% \\\\")
 println("Share of wealth of the top 10\$\\% \$   & ",convert(Int64, round(100 * p90_wealth_no_subsidy)), " & ", convert(Int64, round(100 * p90_wealth_subsidy_partial/p90_wealth_no_subsidy)) -100,
 "\\%  &+", convert(Int64, round(100 * p90_wealth_subsidy_nb/p90_wealth_no_subsidy)) -100,
 "\\%  &+", convert(Int64, round(100 * p90_wealth_subsidy_b/p90_wealth_no_subsidy)) -100,"\\% \\\\")
println("Share of income of the top 10\$\\% \$   & ",convert(Int64, round(100 * p90_income_tmp_no_subsidy)), " & +", convert(Int64, round(100 * p90_income_tmp_subsidy_partial/p90_income_tmp_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * p90_income_tmp_subsidy_nb/p90_income_tmp_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * p90_income_tmp_subsidy_b/p90_income_tmp_no_subsidy)) -100,"\\% \\\\")
println("Share of consumption of the top 10\$\\% \$   & ",convert(Int64, round(100 * p90_cons_tmp_no_subsidy)), " & ", convert(Int64, round(100 * p90_cons_tmp_subsidy_partial/p90_cons_tmp_no_subsidy)) -100,
"\\%  &+",
convert(Int64, round(100 * p90_cons_tmp_subsidy_nb/p90_cons_tmp_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * p90_cons_tmp_subsidy_b/p90_cons_tmp_no_subsidy)) -100,"\\% \\\\")


# println("Avg urban-rural wealth rate & ",round(urban_rural_wealth_ratio_model_no_subsidy,digits = 2), "&",round(urban_rural_wealth_ratio_model_subsidy_nb,digits = 2),
# "&",round(urban_rural_wealth_ratio_model_subsidy_b,digits = 2)," \\\\")

# Workers

y_value = (exp.((V_saved_subsidy_b[:,1] -V_saved_no_subsidy[:,1]) * (1.0 - Baseline_parameter.β) ) ) .- 1
y_mat = 100* mat_creator(y_value);
plot([y_mat[1,1:5:end],y_mat[3,1:5:end],y_mat[13,1:5:end],y_mat[16,1:5:end]], label=["Low rural & urban \\theta" "High rural, low urban \\theta" "Low rural, high urban \\theta" "High rural & urban \\theta"]
,legend=:outerbottom ,
linewidth = 2,linestyle = [:solid :dash :dot :dashdot],ylims = [-20.0,35.0], xlabel = "Wealth",
ylabel = "Welfare change",xticks = ([5,14,18,21 ],["\$50","\$1000", "\$2500", "\$6000"]),yticks = ([-20,-10,0,10,20,30 ],["-20","-10","0","10","20","30"]),
grid = false,size = (800, 800),legendcolumns=2,marker = [:none :none :circle :none],
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure2a.svg")
#Farmers
y_value = (exp.((V_saved_subsidy_b[:,2] -V_saved_no_subsidy[:,2]) * (1.0 - Baseline_parameter.β) ) ) .- 1
y_mat = 100* mat_creator(y_value);
plot([y_mat[1,1:5:end],y_mat[3,1:5:end],y_mat[13,1:5:end],y_mat[16,1:5:end]], label=["Low rural & urban \\theta" "High rural, low urban \\theta" "Low rural, high urban \\theta" "High rural & urban \\theta"]
,legend=:outerbottom ,
linewidth = 2,linestyle = [:solid :dash :dot :dashdot],ylims = [-30.0,10.0], xlabel = "Wealth",
ylabel = "Welfare change",xticks = ([5,14,18,21 ],["\$50","\$1000", "\$2500", "\$6000"]),yticks = ([-30,-20,-10,0,10,20 ],["-30","-20","-10","0","10","20",]),
grid = false,size = (800, 800),legendcolumns=2,marker = [:none :none :circle :none],
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure2b.svg")

#Plotting of tauSgrid values for balanced budget: 3b and 3d requires transition to be loaded
plot(-1*reverse(τ_grid),[reverse(aggregate_consumption_subsidy_b_grid_growth_smoothed)
 ,reverse(current_worker_pop_subsidy_b_grid_growth_smoothed), reverse(undernutritioned_subsidy_b_grid_growth_smoothed)]
, label=["Aggregate consumption"   "Urbanization rate"  "Undernourished households"],
legend=(0.42,-0.17) ,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-26.0,5.0], xlabel = "",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = true,size = (800,800),yticks = ([-25,-20,-15,-10,-5,0,5 ],["-25","-20","-15","-10","-5","0","5"]),
xticks = ([0.0,0.2,0.4,0.45,0.6,0.8],["0.0","0.2","0.4","","0.6","0.8"]),
tickfontsize = 18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,fontfamily="times",
bottom_margin = 37*Plots.mm)
annotate!(0.44, -27.10, text("τ*", :left,"times",18,:red))
annotate!(0.26, -28.7, text("Subsidy rate per unit of staple input", :left,"times",18))
savefig("Figure3a.svg")

plot(-1*reverse(τ_grid),[reverse(100*welfare_trans_optimal_real), reverse(avg_agri_prod_rural_subsidy_b_grid_growth_smoothed),
reverse(avg_labor_prod_urban_subsidy_b_grid_growth_smoothed)]
, label=["Welfare change" "Mean rural ability" "Mean urban ability" ],legend=(0.42,-0.17) ,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-5.2,2.0], xlabel = "",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = true,size = (800,800),xticks = ([0.0,0.2,0.4,0.45,0.6,0.8],["0.0","0.2","0.4","","0.6","0.8"]),
yticks = ([-5,-2.5,0,2.0 ],["-5","-2.5","0","2.0"]),
tickfontsize = 18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,fontfamily="times",
bottom_margin = 37*Plots.mm)
annotate!(0.44, -5.45, text("τ*", :left,"times",18,:red))
annotate!(0.26, -5.8, text("Subsidy rate per unit of staple input", :left,"times",18))

savefig("Figure3b.svg")

plot(-1*reverse(τ_grid),[reverse(staple_productivity_subsidy_b_grid_smoothed),
reverse(var_MPX_subsidy_b_grid_growth_smoothed), reverse(APG_subsidy_b_grid_growth_smoothed)]
, label=[ "Staple productivity" "Dispersion of returns on fertilizer" "APG"  ],
legend=(0.42,-0.17),
linewidth = 2,linestyle = [:solid :dash :dashdot :dot],ylims = [-2.0,100.0], xlabel = "",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],xticks = ([0.0,0.2,0.4,0.45,0.6,0.8],["0.0","0.2","0.4","","0.6","0.8"]),
grid = true,size = (800,800),
tickfontsize = 18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,fontfamily="times",
bottom_margin = 37*Plots.mm)
annotate!(0.44, -5.5, text("τ*", :left,"times",18,:red))
annotate!(0.26, -9.6, text("Subsidy rate per unit of staple input", :left,"times",18))
savefig("Figure3c.svg")

welfare_trans_optimal_real_workers_smoothed = movmean(100*welfare_trans_optimal_real_workers,5) 
welfare_trans_optimal_real_farmers_smoothed = movmean(100*welfare_trans_optimal_real_farmers,5)
welfare_trans_optimal_real_workers_smoothed = welfare_trans_optimal_real_workers_smoothed .- welfare_trans_optimal_real_workers_smoothed[21]
welfare_trans_optimal_real_farmers_smoothed = welfare_trans_optimal_real_farmers_smoothed .- welfare_trans_optimal_real_farmers_smoothed[21]


plot(-1*reverse(τ_grid),[reverse(100*welfare_trans_optimal_real),
 reverse(welfare_trans_optimal_real_workers_smoothed),
reverse(welfare_trans_optimal_real_farmers_smoothed)]
, label=["Welfare change" "Urban welfare change" "Rural welfare change" ],legend=(0.42,-0.13),
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-2.0,1.0], xlabel = "",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = true,size = (800,800),xticks = ([0.0,0.2,0.4,0.45,0.6,0.8],["0.0","0.2","0.4","","0.6","0.8"]),
yticks = ([-2,-1.0,0,1.0],["-2","-1","0","1"]),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times",
bottom_margin = 35*Plots.mm)
annotate!(0.44, -2.09, text("τ*", :left,"times",14,:red))
annotate!(0.26, -2.18, text("Subsidy rate per unit of staple input", :left,"times",14))
savefig("Figure3d.svg")


welfare_trans_optimal_real_workers_smoothed_lr = movmean(100*welfare_trans_optimal_real_workers_lr,5) 
welfare_trans_optimal_real_farmers_smoothed_lr = movmean(100*welfare_trans_optimal_real_farmers_lr,5)
welfare_subsidy_b_grid_real_smoothed_lr = movmean(100*welfare_subsidy_b_grid_real,5)
welfare_trans_optimal_real_workers_smoothed_lr = welfare_trans_optimal_real_workers_smoothed_lr .- welfare_trans_optimal_real_workers_smoothed_lr[21]
welfare_trans_optimal_real_farmers_smoothed_lr = welfare_trans_optimal_real_farmers_smoothed_lr .- welfare_trans_optimal_real_farmers_smoothed_lr[21]
welfare_subsidy_b_grid_real_smoothed_lr  = welfare_subsidy_b_grid_real_smoothed_lr .- welfare_subsidy_b_grid_real_smoothed_lr[21]

plot(-1*reverse(τ_grid),[reverse(welfare_subsidy_b_grid_real_smoothed_lr), reverse(welfare_trans_optimal_real_workers_smoothed_lr),
reverse(welfare_trans_optimal_real_farmers_smoothed_lr)]
, label=["Welfare change" "Urban welfare change" "Rural welfare change" ],legend=(0.42,-0.13),
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-16.0,8.0], xlabel = "",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = true,size = (800,800),xticks = ([0.0,0.2,0.4,0.45,0.6,0.8],["0.0","0.2","0.4","","0.6","0.8"]),
yticks = ([-15,-10,-5,0,5 ],["-15","-10","-5","0","5"]),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times",
bottom_margin = 35*Plots.mm)
annotate!(0.44, -16.7, text("τ*", :left,"times",14,:red))
annotate!(0.26, -17.7, text("Subsidy rate per unit of staple input", :left,"times",14))
savefig("Figure3e.svg")

plot(-1*reverse(τ_grid),[reverse(c_S_worker_sum_subsidy_b_grid_growth_smoothed),
reverse(c_S_staple_sum_subsidy_b_grid_growth_smoothed), reverse(c_S_cashcrop_sum_subsidy_b_grid_growth_smoothed)]
, label=[ "Workers" "Staple farmers" "Cash crop farmers"  ],legend=(0.42,-0.17),
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-20.0,45.0],
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),yticks = ([-20,-10,0,10,20,30,40 ],["-20","-10","0","10","20","30","40"]),
tickfontsize = 18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,fontfamily="times",
bottom_margin = 37*Plots.mm)
annotate!(0.44, -22.2, text("τ*", :left,"times",18,:red))
annotate!(0.26, -25.6, text("Subsidy rate per unit of staple input", :left,"times",18))
savefig("Figure3f.svg")

plot(-1*reverse(τ_grid),[reverse(c_B_worker_sum_subsidy_b_grid_growth_smoothed),
reverse(c_B_staple_sum_subsidy_b_grid_growth_smoothed), reverse(c_B_cashcrop_sum_subsidy_b_grid_growth_smoothed)]
, label=[ "Workers" "Staple farmers" "Cash crop farmers"  ],legend=(0.42,-0.17),
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-20.0,45.0],
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),yticks = ([-20,-10,0,10,20,30,40 ],["-20","-10","0","10","20","30","40"]),
tickfontsize = 18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,fontfamily="times",
bottom_margin = 37*Plots.mm)
annotate!(0.44, -22.2, text("τ*", :left,"times",18,:red))
annotate!(0.26, -25.6, text("Subsidy rate per unit of staple input", :left,"times",18))
savefig("Figure3g.svg")

plot(-1*reverse(τ_grid),[reverse(c_M_worker_sum_subsidy_b_grid_growth_smoothed),
reverse(c_M_staple_sum_subsidy_b_grid_growth_smoothed), reverse(c_M_cashcrop_sum_subsidy_b_grid_growth_smoothed)]
, label=[ "Workers" "Staple farmers" "Cash crop farmers"  ],legend=(0.42,-0.17),
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-20.0,45.0], 
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),yticks = ([-20,-10,0,10,20,30,40 ],["-20","-10","0","10","20","30","40"]),
tickfontsize = 18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,fontfamily="times",
bottom_margin = 37*Plots.mm)
annotate!(0.44, -22.2, text("τ*", :left,"times",18,:red))
annotate!(0.26, -25.6, text("Subsidy rate per unit of staple input", :left,"times",18))
savefig("Figure3h.svg")

plot(-1*reverse(τ_grid),[reverse(food_share_worker_subsidy_b_grid_growth_smoothed),
reverse(food_share_staple_subsidy_b_grid_growth_smoothed), reverse(food_share_cashcrop_subsidy_b_grid_growth_smoothed)]
, label=[ "Workers" "Staple farmers" "Cash crop farmers"  ],legend=:outerbottom,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [0.0,40.0], xlabel = "Subsidy rate per unit of staple input",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure3i.svg")



# Compared with the DiD data

println("& Input use & Staple yield & Cash crop yield & Relative price & \\% staple land & \\% undernourished & \\% rural pop.\\tabularnewline ")
println("\\hline Diff-in-diff & +37\\%{*}{*} & +20\\%{*}{*}{*} & -19\\% & +24\\%{*}{*} & -2\\% & -23\\%{*}{*}{*} & -1 \\% \\tabularnewline")
println("Model \$\\tau_{S}=0.41\$   &  +",convert(Int64, round(fertilizer_use_subsidy_b_grid_growth_smoothed[11])), "\\%&+",
 convert(Int64, round(staple_productivity_subsidy_b_grid_smoothed[11])) ,"\\%&+", 
 convert(Int64, round(cashcrop_productivity_value_subsidy_b_grid_smoothed[11])),"\\%&+", 
 convert(Int64, round(cashcrop_price_subsidy_b_grid_smoothed[11])),"\\%&+", 
 convert(Int64, round(relative_land_to_staples_subsidy_b_grid_smoothed[11])),"\\%&", 
 convert(Int64, round(undernutritioned_subsidy_b_grid_growth_smoothed[11])),"\\%&+", 
 convert(Int64, round(current_worker_pop_subsidy_b_grid_growth_smoothed[11])) ,"\\%\\\\ ")
println("Model \$\\tau_{S}=0.81\$   &  +",convert(Int64, round(fertilizer_use_subsidy_b_grid_growth_smoothed[1])), "\\%&+",
 convert(Int64, round(staple_productivity_subsidy_b_grid_smoothed[1])) ,"\\%&+", 
 convert(Int64, round(cashcrop_productivity_value_subsidy_b_grid_smoothed[1])),"\\%&+", 
 convert(Int64, round(cashcrop_price_subsidy_b_grid_smoothed[1])),"\\%&+", 
 convert(Int64, round(relative_land_to_staples_subsidy_b_grid_smoothed[1])),"\\%&", 
 convert(Int64, round(undernutritioned_subsidy_b_grid_growth_smoothed[1])),"\\%&", 
 convert(Int64, round(current_worker_pop_subsidy_b_grid_growth_smoothed[1])) ,"\\% \\\\ ")



# Infrastructure results table for the presentation

println("\\alert<1>{Consumption equivalent welfare}  &",round(100 * welfare_subsidy_nb,digits = 1), "&  ",round(100 * welfare_inf,digits = 1), "&", round(100 * welfare_inf_sp,digits = 1) ,"\\\\ ")
println("Prices: \$p_B, p_M\$ &",convert(Int64, round(100 * prices_subsidy_nb[1]/prices_no_subsidy[1]))," , "
,convert(Int64, round(100 * prices_subsidy_nb[2]/prices_no_subsidy[2]))," & ",convert(Int64, round(100 * prices_inf[1]/prices_no_subsidy[1]))," , "
,convert(Int64, round(100 * prices_inf[2]/prices_no_subsidy[2]))," & " ,convert(Int64, round(100 * prices_inf_sp[1]/prices_no_subsidy[1]))," , "
,convert(Int64, round(100 * prices_inf_sp[2]/prices_no_subsidy[2])),"\\\\")
println("\\alert<2>{Staple production} & ", convert(Int64, round(100 * prod_staple_subsidy_nb/prod_staple_no_subsidy)) ," & ", convert(Int64, round(100 * prod_staple_inf/prod_staple_no_subsidy)) ,"  &", convert(Int64, round(100 * prod_staple_inf_sp/prod_staple_no_subsidy)) ," \\\\")
println("\\alert<2>{Staple productivity} & ", convert(Int64, round(100 * staple_productivity_subsidy_nb/staple_productivity_no_subsidy)) ," & ", convert(Int64, round(100 * staple_productivity_inf/staple_productivity_no_subsidy)) ,"  &", convert(Int64, round(100 * staple_productivity_inf_sp/staple_productivity_no_subsidy)) ," \\\\")
println("\\alert<3>{Cash crop production} & ", convert(Int64, round(100 * prod_cashcrop_subsidy_nb/prod_cashcrop_no_subsidy)) ,"  & ", convert(Int64, round(100 * prod_cashcrop_inf/prod_cashcrop_no_subsidy)) ,"  &", convert(Int64, round(100 * prod_cashcrop_inf_sp/prod_cashcrop_no_subsidy)) ," \\\\")
println("\\alert<3>{Cash crop productivity} & ", convert(Int64, round(100 * cashcrop_productivity_subsidy_nb/cashcrop_productivity_no_subsidy)) ," & ", convert(Int64, round(100 * cashcrop_productivity_inf/cashcrop_productivity_no_subsidy)) ,"  &", convert(Int64, round(100 * cashcrop_productivity_inf/cashcrop_productivity_no_subsidy)) ," \\\\")
println("\\alert<4>{Share of financially constrained farmers} & ",convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_subsidy_nb)), "&",convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_inf)), "&",convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_inf_sp))," \\\\")
println("Share of staple-only farmers & ",convert(Int64,round(100 * rural_pop_only_staples_model_subsidy_nb)), "&",convert(Int64,round(100 * rural_pop_only_staples_model_inf)), "&",convert(Int64,round(100 * rural_pop_only_staples_model_inf_sp))," \\\\")
println("\\alert<5>{Urbanization rate} & ",convert(Int64,round(100 * current_worker_pop_subsidy_nb)), "&",convert(Int64,round(100 * current_worker_pop_inf)), "&",convert(Int64,round(100 * current_worker_pop_inf_sp))," \\\\")
println("APG & ",round(APG_subsidy_nb,digits = 2), "&",round(APG_inf,digits = 2), "&",round(APG_inf_sp,digits = 2)," \\\\")
println("Consumption & ", convert(Int64, round(100 * aggregate_consumption_subsidy_nb/aggregate_consumption_no_subsidy)) ,
" & ", convert(Int64, round(100 * aggregate_consumption_inf/aggregate_consumption_no_subsidy)) ,
"  &", convert(Int64, round(100 * aggregate_consumption_inf_sp/aggregate_consumption_no_subsidy)) ," \\\\")
println("Transaction cost & ", convert(Int64, round(100 * transaction_cost_loss_subsidy_nb/transaction_cost_loss_no_subsidy)) ,
" & ", convert(Int64, round(100 * transaction_cost_loss_inf/transaction_cost_loss_no_subsidy)) ,
"  &", convert(Int64, round(100 * transaction_cost_loss_inf_sp/transaction_cost_loss_no_subsidy)) ," \\\\")
println("Current account surplus \\% of GDP & ",round(100*current_account_residual_subsidy_nb,digits = 2), "&",round(100*current_account_residual_inf,digits = 2),
 "&",round(100*current_account_residual_inf_sp,digits = 2)," \\\\")

 # Infrastructure results table for the paper

 println("\\hline & FISP & Infrastructure & Infrastructure w/ spillovers\\\\ ")
 println("   & foreign aid  &  \$ Q_S=",  round(infra_parameter_nsp_nb.Q_S,digits = 1), " \$ &  \$ Q_S=",  round(infra_parameter_sp_nb.Q_S,digits = 1), " \$ \\& \$ F_W=",  convert(Int64,round(infra_parameter_sp_nb.F_W,digits = 0)), "\$  \\\\ ")
 println("\\hline 
 \\\\[-1.8ex] \\textbf{Prices and Aggregates} & & &   \\\\[1ex]")
 println("Cash crop price, \$p_B\$ & +"
 ,convert(Int64, round(100 * prices_subsidy_nb[1]/prices_no_subsidy[1])) -100,"\\% &"
,convert(Int64, round(100 * prices_inf[1]/prices_no_subsidy[1])) -100,"\\% &"
,convert(Int64, round(100 * prices_inf_sp[1]/prices_no_subsidy[1]))-100,
"\\% \\\\")
println("Manufacturing price, \$p_M\$ &+",convert(Int64, round(100 * prices_subsidy_nb[2]/prices_no_subsidy[2]))-100,"\\%"
," & ",convert(Int64, round(100 * prices_inf[2]/prices_no_subsidy[2]))-100,"\\%"
,"&",convert(Int64, round(100 * prices_inf_sp[2]/prices_no_subsidy[2]))-100,"\\%\\\\")
println("Wage rate, \$w\$ &+",convert(Int64, round(100 * w_subsidy_nb/w_no_subsidy))-100
,"\\%&",convert(Int64, round(100 * w_inf/w_no_subsidy))-100,"\\%  & \$",convert(Int64, round(100 * w_inf_sp/w_no_subsidy))-100,"\\%\$\\\\")
println("Consumption &+",
convert(Int64, round(100 * aggregate_consumption_subsidy_nb/aggregate_consumption_no_subsidy)) -100,"\\%  & +",
convert(Int64, round(100 * aggregate_consumption_inf/aggregate_consumption_no_subsidy)) -100,"\\%  & +",
convert(Int64, round(100 * aggregate_consumption_inf_sp/aggregate_consumption_no_subsidy)) -100,"\\% \\\\")
println("Nominal output & +", convert(Int64, round(100 * nominal_GDP_subsidy_nb/nominal_GDP_no_subsidy)) -100,"\\%  &", 
convert(Int64, round(100 * nominal_GDP_inf/nominal_GDP_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * nominal_GDP_inf_sp/nominal_GDP_no_subsidy)) -100,"\\% \\\\")
println("Share of cash crop exported & ",convert(Int64, round(100 * exportshare_cashcrop_subsidy_nb/exportshare_cashcrop_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * exportshare_cashcrop_inf/exportshare_cashcrop_no_subsidy)) -100,
"\\%  &", convert(Int64, round(100 * exportshare_cashcrop_inf_sp/exportshare_cashcrop_no_subsidy)) -100,"\\% \\\\")
println("Transaction cost & +", convert(Int64, round(100 * transaction_cost_loss_subsidy_nb/transaction_cost_loss_no_subsidy)) -100,
"\\%  &",convert(Int64, round(100 * transaction_cost_loss_inf/transaction_cost_loss_no_subsidy)) -100, 
"\\%  &", convert(Int64, round(100 * transaction_cost_loss_inf_sp/transaction_cost_loss_no_subsidy)) -100,"\\% \\\\")
println("Current account surplus \\% of GDP & ",convert(Int64,round(100*current_account_residual_subsidy_nb,digits = 0)),
 "\\%& +",convert(Int64,round(100*current_account_residual_inf,digits = 0)),
 "\\%& +",convert(Int64,round(100*current_account_residual_inf_sp,digits = 0)),"\\%  \\\\[1ex]")
println("\\textbf{Production}& & &   \\\\[1ex]")
#println("Share of staple-only farmers &+",convert(Int64,round(100 * rural_pop_only_staples_model_subsidy_nb/rural_pop_only_staples_model_no_subsidy))-100,
#"\\% &+",convert(Int64,round(100 * rural_pop_only_staples_model_inf/rural_pop_only_staples_model_no_subsidy))-100,
#"\\% & ",convert(Int64,round(100 * rural_pop_only_staples_model_inf_sp/rural_pop_only_staples_model_no_subsidy))-100,"\\% \\\\")
println("Staple production  & +", convert(Int64, round(100 * prod_staple_subsidy_nb/prod_staple_no_subsidy)) -100 ,"\\%  &",
convert(Int64, round(100 * prod_staple_inf/prod_staple_no_subsidy)) -100 ,"\\%  &+", 
convert(Int64, round(100 * prod_staple_inf_sp/prod_staple_no_subsidy)) -100,"\\% \\\\")
println("Staple productivity & +", convert(Int64, round(100 * staple_productivity_subsidy_nb/staple_productivity_no_subsidy)) -100,"\\% &",
convert(Int64, round(100 * staple_productivity_inf/staple_productivity_no_subsidy)) -100,"\\% &+",
convert(Int64, round(100 * staple_productivity_inf_sp/staple_productivity_no_subsidy)) -100,"\\% \\\\")
println("Cash crop production  & ", 
convert(Int64, round(100 * prod_cashcrop_subsidy_nb/prod_cashcrop_no_subsidy)) -100,"\\%  & ",
convert(Int64, round(100 * prod_cashcrop_inf/prod_cashcrop_no_subsidy)) -100,"\\%  & +",
convert(Int64, round(100 * prod_cashcrop_inf_sp/prod_cashcrop_no_subsidy))-100 ,"\\% \\\\")
println("Cash crop productivity & +", 
convert(Int64, round(100 * cashcrop_productivity_subsidy_nb/cashcrop_productivity_no_subsidy))-100 ,"\\% & +",
convert(Int64, round(100 * cashcrop_productivity_inf/cashcrop_productivity_no_subsidy))-100 ,"\\% & +",
convert(Int64, round(100 * cashcrop_productivity_inf_sp/cashcrop_productivity_no_subsidy))-100 ,"\\% \\\\")
println("Share of land devoted to staples & +",
convert(Int64,round(mean_land_share_staples_subsidy_nb/mean_land_share_staples_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(mean_land_share_staples_inf/mean_land_share_staples_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(mean_land_share_staples_inf_sp/mean_land_share_staples_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Share of farmers without surplus & ",
convert(Int64,round(fraction_model_subsidy_nb/fraction_model_no_subsidy*100,digits = 0) -100), "\\%  &+",
convert(Int64,round(fraction_model_inf/fraction_model_no_subsidy*100,digits = 0) -100), "\\%  &+",
convert(Int64,round(fraction_model_inf_sp/fraction_model_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")

println("Share of financially constrained farmers & ",
convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_subsidy_nb/fraction_cashcrop_suboptimal_model_no_subsidy))-100, "\\% &",
convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_inf/fraction_cashcrop_suboptimal_model_no_subsidy))-100, "\\% &",
convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_inf_sp/fraction_cashcrop_suboptimal_model_no_subsidy))-100,"\\% \\\\")
println("Manufacturing production & +",
convert(Int64,round(100 * prod_manuf_subsidy_nb/prod_manuf_no_subsidy)) -100, "\\% &+",
convert(Int64,round(100 * prod_manuf_inf/prod_manuf_no_subsidy)) -100,        "\\% &+",
convert(Int64,round(100 * prod_manuf_inf_sp/prod_manuf_no_subsidy))-100,"\\% \\\\")

println("Urbanization rate &",
convert(Int64,round(100 * current_worker_pop_subsidy_nb/current_worker_pop_no_subsidy)) -100, "\\% &+",
convert(Int64,round(100 * current_worker_pop_inf/current_worker_pop_no_subsidy)) -100, "\\% &+",
convert(Int64,round(100 * current_worker_pop_inf_sp/current_worker_pop_no_subsidy))-100,"\\% \\\\")

println("Agricultural productivity gap &+",
convert(Int64,round(APG_subsidy_nb/APG_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(APG_inf/APG_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(APG_inf_sp/APG_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Average agricultural ability & ",
convert(Int64,round(avg_agri_prod_rural_subsidy_nb/avg_agri_prod_rural_no_subsidy*100,digits = 0) -100), "\\%  &+",
convert(Int64,round(avg_agri_prod_rural_inf/avg_agri_prod_rural_no_subsidy*100,digits = 0) -100), "\\%  &+",
convert(Int64,round(avg_agri_prod_rural_inf_sp/avg_agri_prod_rural_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Average worker ability & ",
convert(Int64,round(avg_labor_prod_urban_subsidy_nb/avg_labor_prod_urban_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(avg_labor_prod_urban_inf/avg_labor_prod_urban_no_subsidy*100,digits = 0) -100), "\\%  & +",
convert(Int64,round(avg_labor_prod_urban_inf_sp/avg_labor_prod_urban_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Dispersion in ARPX &+",
convert(Int64,round(var_MPX_subsidy_nb/var_MPX_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(var_MPX_inf/var_MPX_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(var_MPX_inf_sp/var_MPX_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")

println("Dispersion in ARPX for cash crop farmers &",
convert(Int64,round(var_MPX_cashcrop_subsidy_nb/var_MPX_cashcrop_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(var_MPX_cashcrop_inf/var_MPX_cashcrop_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(var_MPX_cashcrop_inf_sp/var_MPX_cashcrop_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Dispersion in ARPX for staple farmers &+",
convert(Int64,round(var_MPX_staples_S_subsidy_nb/var_MPX_staples_S_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(var_MPX_staples_S_inf/var_MPX_staples_S_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(var_MPX_staples_S_inf_sp/var_MPX_staples_S_no_subsidy*100,digits = 0) - 100),"\\% \\\\[1ex]")
println("\\textbf{Welfare and Inequality}& & & &   \\\\[1ex]")
println("Consumption equivalent welfare   &  +",
    round(100 * welfare_transition_subsidy_nb, digits=1), "\\% &+",
    round(100 * welfare_transition_subsidy_infra, digits=1), "\\% &+",
    round(100 * welfare_transition_subsidy_infra_sp, digits=1), "\\% \\\\ ")
println("Consumption equivalent welfare - urban  & +",
round(100 * welfare_subsidy_nb_real_workers,digits = 1), "\\% &+",
round(100 * welfare_inf_real_workers,digits = 1), "\\% &+",
round(100 * welfare_inf_sp_real_workers,digits = 1) ,"\\% \\\\ ")
println("Consumption equivalent welfare - rural  & 
  ",round(100 * welfare_subsidy_nb_real_farmers_lr,digits = 1), "\\% &+",
round(100 * welfare_inf_real_farmers,digits = 1), "\\% &+",
round(100 * welfare_inf_sp_real_farmers,digits = 1) ,"\\% \\\\ ")
println("Share of undernourished  & ",
convert(Int64,round(undernourished_subsidy_nb/undernourished_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(undernourished_inf/undernourished_no_subsidy*100,digits = 0) -100), "\\%  &",
convert(Int64,round(undernourished_inf_sp/undernourished_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Avg urban-rural consumption ratio & +",
convert(Int64,round(100*urban_rural_consumption_ratio_model_subsidy_nb/urban_rural_consumption_ratio_model_no_subsidy - 100)),"\\%&",
convert(Int64,round(100*urban_rural_consumption_ratio_model_inf/urban_rural_consumption_ratio_model_no_subsidy - 100)),"\\%&",
convert(Int64,round(100*urban_rural_consumption_ratio_model_inf_sp/urban_rural_consumption_ratio_model_no_subsidy - 100)),"\\% \\\\")
println("Avg urban-rural income ratio &+",
convert(Int64,round(100*urban_rural_inc_ratio_model_subsidy_nb/urban_rural_inc_ratio_model_no_subsidy - 100)),"\\%&",
convert(Int64,round(100*urban_rural_inc_ratio_model_inf/urban_rural_inc_ratio_model_no_subsidy - 100)),"\\%&",
convert(Int64,round(100*urban_rural_inc_ratio_model_inf_sp/urban_rural_inc_ratio_model_no_subsidy - 100)),"\\% \\\\")
println("Avg urban-rural wealth ratio &+",
convert(Int64,round(100*urban_rural_wealth_ratio_model_subsidy_nb/urban_rural_wealth_ratio_model_no_subsidy - 100)),"\\%&",
convert(Int64,round(100*urban_rural_wealth_ratio_model_inf/urban_rural_wealth_ratio_model_no_subsidy - 100)),"\\%&",
convert(Int64,round(100*urban_rural_wealth_ratio_model_inf_sp/urban_rural_wealth_ratio_model_no_subsidy - 100)),"\\% \\\\")

 println("Share of wealth of the top 10\$\\% \$   & +", 
 convert(Int64, round(100 * p90_wealth_subsidy_nb/p90_wealth_no_subsidy)) -100, "\\%  &+",
 convert(Int64, round(100 * p90_wealth_inf/p90_wealth_no_subsidy)) -100,"\\%  &",
 convert(Int64, round(100 * p90_wealth_inf_sp/p90_wealth_no_subsidy)) -100,"\\% \\\\")
println("Share of income of the top 10\$\\% \$   & +", 
convert(Int64, round(100 * p90_income_tmp_subsidy_nb/p90_income_tmp_no_subsidy)) -100,"\\%  &", 
convert(Int64, round(100 * p90_income_tmp_inf/p90_income_tmp_no_subsidy)) -100,"\\%  &",
convert(Int64, round(100 * p90_income_tmp_inf_sp/p90_income_tmp_no_subsidy)) -100,"\\% \\\\")
println("Share of consumption of the top 10\$\\% \$  & +", 
convert(Int64, round(100 * p90_cons_tmp_subsidy_nb/p90_cons_tmp_no_subsidy)) -100,
"\\%  &",convert(Int64, round(100 * p90_cons_tmp_inf/p90_cons_tmp_no_subsidy)) -100,
"\\%  &",convert(Int64, round(100 * p90_cons_tmp_inf_sp/p90_cons_tmp_no_subsidy)) -100,"\\% \\\\")



# Calibration Tables - effects of shutting down frictions for the paper

println("Urban population & 18\\% & ",convert(Int64, round(100 * current_worker_pop_subsidy_b)),"\$\\% \$ &",convert(Int64, round(100 * current_worker_pop_no_cbar_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * current_worker_pop_no_QS_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * current_worker_pop_no_cbarQS_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * current_worker_pop_no_F_W_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * current_worker_pop_no_FM_B_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * current_worker_pop_no_κ_subsidy_b)),"\$\\% \$ & \\tabularnewline")
println("Share of 2010 exported cash crops in total prod. of cash crops & 73\$\\% \$ & ",convert(Int64, round(100 * exportshare_cashcrop_subsidy_b)),"\$\\% \$ &",convert(Int64, round(100 * exportshare_cashcrop_no_cbar_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * exportshare_cashcrop_no_QS_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * exportshare_cashcrop_no_cbarQS_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * exportshare_cashcrop_no_F_W_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * exportshare_cashcrop_no_FM_B_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * exportshare_cashcrop_no_κ_subsidy_b)),"\$\\% \$ & \\tabularnewline")
println("Aggregate cost of program (\\% GDP) & 3\$\\% \$ & ",convert(Int64, round(100 * program_spending_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * program_spending_no_cbar_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * program_spending_no_QS_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * program_spending_no_cbarQS_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * program_spending_no_F_W_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * program_spending_no_FM_B_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * program_spending_no_κ_subsidy_b)),"\$\\% \$ \\tabularnewline")
println("Rural-urban migration rate & 1\$\\% \$ & ",convert(Int64, round(100 * mig_rate_model_subsidy_b)),"\$\\% \$ &",convert(Int64, round(100 * mig_rate_model_no_cbar_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * mig_rate_model_no_QS_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * mig_rate_model_no_cbarQS_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * mig_rate_model_no_F_W_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * mig_rate_model_no_FM_B_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * mig_rate_model_no_κ_subsidy_b)),"\$\\% \$ & \\tabularnewline")
println("Share of land devoted to staples  & 70\\% & ",convert(Int64, round(100 * (mean_land_share_to_staples_among_cc_model_subsidy_b*current_cashcrop_pop_subsidy_b
+ current_staple_pop_subsidy_b)/(1 - current_worker_pop_subsidy_b ) )),"\$\\% \$ &",
convert(Int64, round(100 * (mean_land_share_to_staples_among_cc_model_no_cbar_subsidy_b*current_cashcrop_pop_subsidy_b
 + current_staple_pop_subsidy_b)/(1 - current_worker_pop_subsidy_b ) )),"\$\\% \$&",
 convert(Int64, round(100 * (mean_land_share_to_staples_among_cc_model_no_QS_subsidy_b*current_cashcrop_pop_subsidy_b
 + current_staple_pop_subsidy_b)/(1 - current_worker_pop_subsidy_b ) )),"\$\\% \$&",
 convert(Int64, round(100 * (mean_land_share_to_staples_among_cc_model_no_cbarQS_subsidy_b*current_cashcrop_pop_subsidy_b
 + current_staple_pop_subsidy_b)/(1 - current_worker_pop_subsidy_b ) )),"\$\\% \$&",
 convert(Int64, round(100 * (mean_land_share_to_staples_among_cc_model_no_F_W_subsidy_b*current_cashcrop_pop_subsidy_b
 + current_staple_pop_subsidy_b)/(1 - current_worker_pop_subsidy_b ) )),"\$\\% \$&",
 convert(Int64, round(100 * (mean_land_share_to_staples_among_cc_model_no_FM_B_subsidy_b*current_cashcrop_pop_subsidy_b
 + current_staple_pop_subsidy_b)/(1 - current_worker_pop_subsidy_b ) )),"\$\\% \$&",
 convert(Int64, round(100 * (mean_land_share_to_staples_among_cc_model_no_κ_subsidy_b*current_cashcrop_pop_subsidy_b
 + current_staple_pop_subsidy_b)/(1 - current_worker_pop_subsidy_b ) )),"\$\\% \$ & \\tabularnewline")
println("Share of cash crop farmers with suboptimal inputs & 70\\% & ",convert(Int64, round(100 * fraction_cashcrop_suboptimal_model_subsidy_b)),"\$\\% \$ &",convert(Int64, round(100 * fraction_cashcrop_suboptimal_model_no_cbar_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * fraction_cashcrop_suboptimal_model_no_QS_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * fraction_cashcrop_suboptimal_model_no_cbarQS_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * fraction_cashcrop_suboptimal_model_no_F_W_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * fraction_cashcrop_suboptimal_model_no_FM_B_subsidy_b)),"\$\\% \$&",convert(Int64, round(100 * fraction_cashcrop_suboptimal_model_no_κ_subsidy_b)),"\$\\% \$ & \\tabularnewline")
println("Standard deviation of average product of farms & 1.8& ",round((var_APland_subsidy_b).^(1/2), digits = 2)," &",round((var_APland_no_cbar_subsidy_b).^(1/2), digits = 2) ,"&",round((var_APland_no_QS_subsidy_b).^(1/2), digits = 2),"&",round((var_APland_no_cbarQS_subsidy_b).^(1/2), digits = 2),"&",round((var_APland_no_F_W_subsidy_b).^(1/2), digits = 2),"&",round((var_APland_no_FM_B_subsidy_b).^(1/2), digits = 2),"&",round((var_APland_no_κ_subsidy_b).^(1/2), digits = 2)," & \\tabularnewline")
println("Agricultural productivity gap  & 6.5 & ", round(APG_subsidy_b, digits = 1),"& ", round(APG_no_cbar_subsidy_b, digits = 1),"& ", round(APG_no_QS_subsidy_b, digits = 1),"& ", round(APG_no_cbarQS_subsidy_b, digits = 1),"& ", round(APG_no_F_W_subsidy_b, digits = 1),"& ", round(APG_no_FM_B_subsidy_b, digits = 1),"& ", round(APG_no_κ_subsidy_b, digits = 1)," \\tabularnewline")

#println("Aggregate cashcrop to staple fertilizer expenditure ratio s & 2.62 & ",round(exp_ratio_model_subsidy_b, digits = 2)," &",round(exp_ratio_model_no_cbar_subsidy_b, digits = 2),"&",round(exp_ratio_model_no_QS_subsidy_b, digits = 2),"&",round(exp_ratio_model_no_cbarQS_subsidy_b, digits = 2),"&",round(exp_ratio_model_no_F_W_subsidy_b, digits = 2),"&",round(exp_ratio_model_no_FM_B_subsidy_b, digits = 2),"&",round(exp_ratio_model_no_κ_subsidy_b, digits = 2)," & \\tabularnewline")

# Shutting down frictions table for the presentation

println("\\hline & No FISP & FISP & FISP \\\\ ")
println("&   & aid-financed  & gov-financed  \\\\ ")
println("\\hline \\\\[-1.8ex] ")
println("\\alert<1>{Consumption equivalent welfare}  & - &  ",round(100 * welfare_subsidy_nb,digits = 1), "&", round(100 * welfare_subsidy_b,digits = 1) ,"\\\\ ")
println("Prices: \$p_B, p_M\$ & 100, 100 & ",convert(Int64, round(100 * prices_subsidy_nb[1]/prices_no_subsidy[1]))," , "
,convert(Int64, round(100 * prices_subsidy_nb[2]/prices_no_subsidy[2]))," & " ,convert(Int64, round(100 * prices_subsidy_b[1]/prices_no_subsidy[1]))," , "
,convert(Int64, round(100 * prices_subsidy_b[2]/prices_no_subsidy[2])),"\\\\")
println("\\alert<2>{Staple production} & 100 & ", convert(Int64, round(100 * prod_staple_subsidy_nb/prod_staple_no_subsidy)) ,"  &", convert(Int64, round(100 * prod_staple_subsidy_b/prod_staple_no_subsidy)) ," \\\\")
println("\\alert<2>{Staple productivity} & 100 & ", convert(Int64, round(100 * staple_productivity_subsidy_nb/staple_productivity_no_subsidy)) ,"  &", convert(Int64, round(100 * staple_productivity_subsidy_b/staple_productivity_no_subsidy)) ," \\\\")
println("\\alert<3>{Cash crop production} & 100 & ", convert(Int64, round(100 * prod_cashcrop_subsidy_nb/prod_cashcrop_no_subsidy)) ,"  &", convert(Int64, round(100 * prod_cashcrop_subsidy_b/prod_cashcrop_no_subsidy)) ," \\\\")
println("\\alert<3>{Cash crop productivity} & 100 & ", convert(Int64, round(100 * cashcrop_productivity_subsidy_nb/cashcrop_productivity_no_subsidy)) ,"  &", convert(Int64, round(100 * cashcrop_productivity_subsidy_b/cashcrop_productivity_no_subsidy)) ," \\\\")
println("\\alert<4>{Share of financially constrained farmers} & ",convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_no_subsidy)), "&",convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_subsidy_nb)), "&",convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_subsidy_b))," \\\\")
println("Share of staple-only farmers & ",convert(Int64,round(100 * rural_pop_only_staples_model_no_subsidy)), "&",convert(Int64,round(100 * rural_pop_only_staples_model_subsidy_nb)), "&",convert(Int64,round(100 * rural_pop_only_staples_model_subsidy_b))," \\\\")
println("\\alert<5>{Urbanization rate} & ",convert(Int64,round(100 * current_worker_pop_no_subsidy)), "&",convert(Int64,round(100 * current_worker_pop_subsidy_nb)), "&",convert(Int64,round(100 * current_worker_pop_subsidy_b))," \\\\")
println("APG & ",round(APG_no_subsidy,digits = 2), "&",round(APG_subsidy_nb,digits = 2), "&",round(APG_subsidy_b,digits = 2)," \\\\")
println("Consumption & 100 & ", convert(Int64, round(100 * aggregate_consumption_subsidy_nb/aggregate_consumption_no_subsidy)) ,
"  &", convert(Int64, round(100 * aggregate_consumption_subsidy_b/aggregate_consumption_no_subsidy)) ," \\\\")
println("Transaction cost & 100 & ", convert(Int64, round(100 * transaction_cost_loss_subsidy_nb/transaction_cost_loss_no_subsidy)) ,
"  &", convert(Int64, round(100 * transaction_cost_loss_subsidy_b/transaction_cost_loss_no_subsidy)) ," \\\\")
println("Current account surplus \\% of GDP & ",round(100*current_account_residual_no_subsidy,digits = 2), "&",round(100*current_account_residual_subsidy_nb,digits = 2),
 "&",round(100*current_account_residual_subsidy_b,digits = 2)," \\\\")
 println("Avg urban-rural consumption rate & ",round(urban_rural_consumption_ratio_model_no_subsidy,digits = 2), "&",round(urban_rural_consumption_ratio_model_subsidy_nb,digits = 2),
 "&",round(urban_rural_consumption_ratio_model_subsidy_b,digits = 2)," \\\\")







### Stuff that is experimental and not included in the paper (1st period transition, additional figures, etc. ...) ###

#Plotting of tauSgrid values without balanced budget
plot(-1*reverse(τ_grid),[reverse(aggregate_consumption_subsidy_nb_grid_growth_smoothed)
 ,reverse(welfare_subsidy_nb_grid_real_alt_smoothed), reverse(undernutritioned_subsidy_nb_grid_growth_smoothed)]
, label=["Aggregate consumption"  "Welfare change" "Undernourished households"],legend=:outerbottom ,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-26.0,10.0], xlabel = "Subsidy rate per unit of staple input",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),yticks = ([-25,-20,-15,-10,-5,0,5,10 ],["-25","-20","-15","-10","-5","0","5","10"]),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure4a.svg")

plot(-1*reverse(τ_grid),[reverse(current_worker_pop_subsidy_nb_grid_growth_smoothed), reverse(avg_agri_prod_rural_subsidy_nb_grid_growth_smoothed),
reverse(avg_labor_prod_urban_subsidy_nb_grid_growth_smoothed)]
, label=[ "Urbanization rate" "Mean rural ability" "Mean urban ability" ],legend=:outerbottom ,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-15.0,10.0], xlabel = "Subsidy rate per unit of staple input",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),yticks = ([-15,-10,-5,0,5,10 ],["-15","-10","-5","0","5","10"]),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure4b.svg")


plot(-1*reverse(τ_grid),[reverse(staple_productivity_subsidy_nb_grid_smoothed),
reverse(var_MPX_subsidy_nb_grid_growth_smoothed), reverse(APG_subsidy_nb_grid_growth_smoothed)]
, label=[ "Staple productivity" "Dispersion of returns on fertilizer" "APG"  ],legend=:outerbottom,
linewidth = 2,linestyle = [:solid :dash :dashdot :dot],ylims = [-5.0,100.0], xlabel = "Subsidy rate per unit of staple input",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure4c.svg")

plot(-1*reverse(τ_grid),[reverse(c_S_worker_sum_subsidy_nb_grid_growth_smoothed),
reverse(c_S_staple_sum_subsidy_nb_grid_growth_smoothed), reverse(c_S_cashcrop_sum_subsidy_nb_grid_growth_smoothed)]
, label=[ "Workers" "Staple farmers" "Cash crop farmers"  ],legend=:outerbottom,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [0.0,70.0], xlabel = "Subsidy rate per unit of staple input",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),yticks = ([0,10,20,30,40 ,50,60],["0","10","20","30","40","50","60"]),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure4d.svg")

plot(-1*reverse(τ_grid),[reverse(c_B_worker_sum_subsidy_nb_grid_growth_smoothed),
reverse(c_B_staple_sum_subsidy_nb_grid_growth_smoothed), reverse(c_B_cashcrop_sum_subsidy_nb_grid_growth_smoothed)]
, label=[ "Workers" "Staple farmers" "Cash crop farmers"  ],legend=:outerbottom,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-20.0,45.0], xlabel = "Subsidy rate per unit of staple input",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),yticks = ([-20,-10,0,10,20,30,40 ],["-20","-10","0","10","20","30","40"]),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure4e.svg")

plot(-1*reverse(τ_grid),[reverse(c_M_worker_sum_subsidy_nb_grid_growth_smoothed),
reverse(c_M_staple_sum_subsidy_nb_grid_growth_smoothed), reverse(c_M_cashcrop_sum_subsidy_nb_grid_growth_smoothed)]
, label=[ "Workers" "Staple farmers" "Cash crop farmers"  ],legend=:outerbottom,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-20.0,45.0], xlabel = "Subsidy rate per unit of staple input",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),yticks = ([-20,-10,0,10,20,30,40 ],["-20","-10","0","10","20","30","40"]),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure4f.svg")

plot(-1*reverse(τ_grid),[reverse(food_share_worker_subsidy_nb_grid_growth_smoothed),
reverse(food_share_staple_subsidy_nb_grid_growth_smoothed), reverse(food_share_cashcrop_subsidy_nb_grid_growth_smoothed)]
, label=[ "Workers" "Staple farmers" "Cash crop farmers"  ],legend=:outerbottom,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [0.0,40.0], xlabel = "Subsidy rate per unit of staple input",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure4g.svg")

# Comparing migration rates under different tau_S
mig_rate_model_subsidy_b_grid
mig_rate_model_subsidy_nb_grid

mig_rate_model_subsidy_b_grid_smoothed = movmean((100*(mig_rate_model_subsidy_b_grid./mig_rate_model_subsidy_b_grid[no_taus]) .- 100),5)
mig_rate_model_subsidy_b_grid_smoothed = mig_rate_model_subsidy_b_grid_smoothed .- mig_rate_model_subsidy_b_grid_smoothed[no_taus]

mig_rate_model_subsidy_nb_grid_smoothed = movmean((100*(mig_rate_model_subsidy_nb_grid./mig_rate_model_subsidy_nb_grid[no_taus]) .- 100),5)
mig_rate_model_subsidy_nb_grid_smoothed = mig_rate_model_subsidy_nb_grid_smoothed .- mig_rate_model_subsidy_nb_grid_smoothed[no_taus]

plot(-1*reverse(τ_grid),[reverse(mig_rate_model_subsidy_b_grid_smoothed),
reverse(mig_rate_model_subsidy_nb_grid_smoothed)], label=[ "Balanced" "Aid" ],legend=:outerbottom,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-60.0,10.0], xlabel = "Subsidy rate per unit of staple input",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle],
grid = false,size = (800,800),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")

# Comparing under Q_S = 0 

# Lifecycle
matcheck= zeros(400,12)
matcheck[:,1:3]= future_occupation_fine_local
matcheck[:,4:5] = s_fine
matcheck[:,6:8] = cons_fine_local
matcheck[:,9:11] = solve_cash_crop_index_B_fine



plot!(τ_grid,100*avg_agri_prod_rural_subsidy_b_grid./(avg_agri_prod_rural_subsidy_b_grid[no_taus]).-100)
plot!(τ_grid,100*avg_labor_prod_urban_subsidy_b_grid./(avg_labor_prod_urban_subsidy_b_grid[no_taus]).-100)

# Staple farmers
y_value =  (cons_fine_local_subsidy_b[:,2]./cons_fine_local_no_subsidy[:,2].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,100)
cons_avg_staple = movmean(cons_avg[1,:],7);

# cashcrop farmers
y_value =  (cons_fine_local_subsidy_b[:,3]./cons_fine_local_no_subsidy[:,3].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,100)
cons_avg_cashcrop = movmean(cons_avg[1,:],7);

# Workers
y_value =  (cons_fine_local_subsidy_b[:,1]./cons_fine_local_no_subsidy[:,1].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,100)
cons_avg_worker = movmean(cons_avg[1,:],7);
plot([cons_avg_worker,cons_avg_cashcrop], label=["Workers" "Farmers"],legend=:topright ,
linewidth = 2,linestyle = [:solid :dash],ylims = [-20.0,25.0], xlabel = "Wealth",
ylabel = "Consumption change",
grid = false,
tickfontsize = 14,xguidefontsize=12,yguidefontsize=12,legendfontsize=12,fontfamily="Times_Roman")
savefig("Figure2a.pdf")

# Staple farmers
y_value =  (cons_fine_local_subsidy_b[:,2]./cons_fine_local_no_subsidy[:,2].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,100)
cons_avg_staple = movmean(cons_avg[1,:],7);

# cashcrop farmers
y_value =  (cons_fine_local_subsidy_b[:,3]./cons_fine_local_no_subsidy[:,3].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,100)
cons_avg_cashcrop = movmean(cons_avg[1,:],7);

# Workers
y_value =  (cons_fine_local_subsidy_b[:,1]./cons_fine_local_no_subsidy[:,1].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,100)
cons_avg_worker = movmean(cons_avg[1,:],7);
plot([cons_avg_worker,cons_avg_cashcrop], label=["Workers" "Farmers"],legend=:topright ,
linewidth = 2,linestyle = [:solid :dash],ylims = [-40.0,25.0], xlabel = "Wealth",
ylabel = "Consumption change",
grid = false,
tickfontsize = 14,xguidefontsize=12,yguidefontsize=12,legendfontsize=12,fontfamily="Times_Roman")
savefig("Figure2b.pdf")

# no_QS_no_subsidy_parameter

# Staple farmers
y_value = (exp.(min.(welfare_val_no_QS_subsidy_b - welfare_val_no_QS_no_subsidy,10000) * (1.0 - Baseline_parameter.β)) .- 1)
y_value =  (cons_fine_local_no_QS_subsidy_b[:,2]./cons_fine_local_no_QS_no_subsidy[:,2].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,100)
cons_avg_staple = movmean(cons_avg[1,:],7);

# cashcrop farmers
y_value =  (cons_fine_local_no_QS_subsidy_b[:,3]./cons_fine_local_no_QS_no_subsidy[:,3].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,100)
cons_avg_cashcrop = movmean(cons_avg[1,:],7);

# Workers
y_value =  (cons_fine_local_no_QS_subsidy_b[:,1]./cons_fine_local_no_QS_no_subsidy[:,1].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,100)
cons_avg_worker = movmean(cons_avg[1,:],7);
plot([cons_avg_worker,cons_avg_cashcrop], label=["Workers" "Farmers"],legend=:topright ,
linewidth = 2,linestyle = [:solid :dash],ylims = [-40.0,25.0], xlabel = "Wealth",
ylabel = "Consumption change",
grid = false,
tickfontsize = 14,xguidefontsize=12,yguidefontsize=12,legendfontsize=12,fontfamily="Times_Roman")


# Infrastructure project without spillovers

# Staple farmers
y_value =  (cons_fine_local_inf[:,2]./cons_fine_local_no_subsidy[:,2].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,100)
cons_avg_staple = movmean(cons_avg[1,:],7);

# cashcrop farmers
y_value =  (cons_fine_local_inf[:,3]./cons_fine_local_no_subsidy[:,3].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,100)
cons_avg_cashcrop = movmean(cons_avg[1,:],7);

# Workers
y_value =  (cons_fine_local_inf[:,1]./cons_fine_local_no_subsidy[:,1].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,100)
cons_avg_worker = movmean(cons_avg[1,:],7);
plot([cons_avg_worker,cons_avg_cashcrop], label=["Workers" "Farmers"],legend=:topright ,
linewidth = 2,linestyle = [:solid :dash],ylims = [-4.0,15.0], xlabel = "Wealth",
ylabel = "Consumption change",
grid = false,
tickfontsize = 14,xguidefontsize=12,yguidefontsize=12,legendfontsize=12,fontfamily="Times_Roman")
savefig("Figure3.pdf")

# Infrastructure project with spillovers

# Staple farmers
y_value =  (cons_fine_local_inf_sp[:,2]./cons_fine_local_no_subsidy[:,2].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,1000)
cons_avg_staple = movmean(cons_avg[1,:],7);

# cashcrop farmers
y_value =  (cons_fine_local_inf_sp[:,3]./cons_fine_local_no_subsidy[:,3].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,1000)
cons_avg_cashcrop = movmean(cons_avg[1,:],7);

# Workers
y_value =  (cons_fine_local_inf_sp[:,1]./cons_fine_local_no_subsidy[:,1].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,1000)
cons_avg_worker = movmean(cons_avg[1,:],7);
plot([cons_avg_worker,cons_avg_cashcrop], label=["Workers" "Farmers"],legend=:bottomright ,
linewidth = 2,linestyle = [:solid :dash],ylims = [-10.0,500.0], xlabel = "Wealth",
ylabel = "Consumption change",
grid = false,
tickfontsize = 14,xguidefontsize=12,yguidefontsize=12,legendfontsize=12,fontfamily="Times_Roman")
savefig("Figure4.pdf")













#Urban population
y_value = (exp.((V_saved_subsidy_b_reshaped -V_saved_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ) .- 1
y_value = (exp.((V_saved_no_cbar_subsidy_b_reshaped -V_saved_no_cbar_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ) .- 1
y_value = (exp.((V_saved_no_QS_subsidy_b_reshaped -V_saved_no_QS_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ) .- 1
y_value = cons_fine_local_subsidy_b[:,1]./cons_fine_local_no_subsidy[:,1] .-1

y_value = cons_fine_local_subsidy_b[:,2]./cons_fine_local_no_subsidy[:,2] .-1


y_value =  (cons_fine_local_no_subsidy[:,1]./cons_fine_local_no_subsidy[:,1].-1);
y_mat = mat_creator(y_value);
plot([future_occupation_fine_local_subsidy_b[:,2],y_value],legend = nothing, linewidth = 2,linestyle = [:solid ],grid = false,tickfontsize = 14)
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
#cons_avg = 100.0 * (sum(y_mat[1:convert(Int64,(Baseline_parameter.n_fine[2]/2)),:],dims=1) +
#    sum(y_mat[convert(Int64,(Baseline_parameter.n_fine[2]/2 +1)):Baseline_parameter.n_fine[2],:],dims=1))/size(y_mat)[1];
#y_mat_applied_employed = 100.0 * y_mat[1:convert(Int64,(Baseline_parameter.n_fine[2]/2)),:];
#contour(y_mat_applied_employed', fill=true,levels = 5, c =cgrad(:Blues_6, rev = true),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,tickfontsize = 18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
cons_avg = movmean(cons_avg[1,:],5);
y_value =  (a_prime_fine_local_no_subsidy[:,1]./a_prime_fine_local_initial[:,1].-1);
y_mat = mat_creator(y_value);
a_prime_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1];
a_prime_avg = movmean(a_prime_avg[1,:],5);
y_value = (exp.(min.(welfare_val_no_subsidy - welfare_val_initial,10) * (1.0 - Baseline_parameter.β)) .- 1)[1:Baseline_parameter.ns_fine];
y_mat = mat_creator(y_value);
welfare_avg = 100.0 * sum(y_mat,dims=1)/size(y_mat)[1];
welfare_avg = min.(welfare_avg,100)
welfare_avg = movmean(welfare_avg[1,:],5);
plot([cons_avg,a_prime_avg,welfare_avg],legend = nothing, linewidth = 2,linestyle = [:solid :dash :dot],grid = false,tickfontsize = 14)
#for i=1:ns_fine
#    if((sum(coeff_λ_2_cashcrop_residual_unconstrained_fine[:,i].<=coeff_λ_2_cashcrop_residual_unconstrained_fine[n_fine[1],i] - tol ) -sum(
#        coeff_λ_2_cashcrop_residual_unconstrained_fine[:,i].<=coeff_λ_2_cashcrop_residual_unconstrained_fine[1,i] + tol)) ==2)
#        println(i)
#    end
#    if((sum(coeff_λ_2_s_fine[:,i].<=coeff_λ_2_s_fine[n_fine[1],i] - tol ) -sum(coeff_λ_2_s_fine[:,i].<=coeff_λ_2_s_fine[1,i] + tol)) ==)
#        println(i)
#    end
#end
savefig("Figure2a.pdf")

# Staple farmers
y_value =  (cons_fine_local_subsidy_b[:,2]./cons_fine_local_no_subsidy[:,2].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,100)
cons_avg_staple = movmean(cons_avg[1,:],5);

# cashcrop farmers
y_value =  (cons_fine_local_subsidy_b[:,3]./cons_fine_local_no_subsidy[:,3].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,100)
cons_avg_cashcrop = movmean(cons_avg[1,:],5);

# Workers
y_value =  (cons_fine_local_subsidy_b[:,1]./cons_fine_local_no_subsidy[:,1].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,100)
cons_avg_worker = movmean(cons_avg[1,:],5);
plot([cons_avg_worker,cons_avg_cashcrop])



y_value =  (a_prime_fine_local_subsidy_b[:,2]./a_prime_fine_no_subsidy[:,2].-1);
y_mat = mat_creator(y_value);
a_prime_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1];
a_prime_avg = min.(a_prime_avg,100)
a_prime_avg = movmean(a_prime_avg[1,:],5);
y_value = (exp.(welfare_val_subsidy_nb- welfare_val_no_subsidy * (1.0 - Baseline_parameter.β)) .- 1)[(Baseline_parameter.ns_fine + 1):(2*Baseline_parameter.ns_fine)];
y_mat = mat_creator(y_value);
welfare_avg = 100.0 * sum(y_mat,dims=1)/size(y_mat)[1];
welfare_avg = min.(welfare_avg,100)
welfare_avg = movmean(welfare_avg[1,:],5);
plot([cons_avg,a_prime_avg,welfare_avg],legend = nothing, linewidth = 2,linestyle = [:solid :dash :dot],grid = false,tickfontsize = 14)
savefig("Figure2b.pdf")
#cashcrop producers
y_value =  (cons_fine_local_no_subsidy[:,3]./cons_fine_local_initial[:,3].-1);
y_mat = mat_creator(y_value);
cons_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1]
cons_avg = min.(cons_avg,100)
cons_avg = movmean(cons_avg[1,:],5);
y_value =  (a_prime_fine_local_no_subsidy[:,3]./a_prime_fine_local_initial[:,3].-1);
y_mat = mat_creator(y_value);
a_prime_avg = 100.0 *sum(y_mat,dims=1)/size(y_mat)[1];
a_prime_avg = min.(a_prime_avg,100)
a_prime_avg = movmean(a_prime_avg[1,:],5);
y_value = (exp.((welfare_val_no_subsidy - welfare_val_initial) * (1.0 - Baseline_parameter.β)) .- 1)[481:end];
y_mat = mat_creator(y_value);
welfare_avg = 100.0 * sum(y_mat,dims=1)/size(y_mat)[1];
welfare_avg = min.(welfare_avg,100)
welfare_avg = movmean(welfare_avg[1,:],5);
plot([cons_avg,a_prime_avg,welfare_avg],legend = nothing, linewidth = 2,linestyle = [:solid :dash :dot],grid = false,tickfontsize = 14)
savefig("Figure2c.pdf")
#contour(y_mat_applied', fill=true,levels = 5, c =cgrad(:Blues_6, rev = true),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,tickfontsize = 14,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
#;colorbar_cticks  = [0,20,40,60,80])


###SOME MORE PLOTS FROM KAROL:

plot(τ_grid, 100 * residual_goods_subsidy_b_grid[2,:])
plot!(τ_grid, 100 * residual_goods_subsidy_b_grid[1, :])
plot!(τ_grid, 100 * residual_goods_subsidy_b_grid[5, :])
RSS=sum(residual_goods_subsidy_b_grid.^2,dims=1)
plot!(τ_grid, 1000 * transpose(RSS))

###### PLOT corr staple-cash choice vs agri prod + (2) urban rural choice vs assets
#order of all 2560 fine states: (4x4=16 prod states) + 160 agrid_fine:
# first iterate all 40/160 asset states for a fixed prod (1 out of 16)
# repeat above for a higher level of Θ with fixed labor_prod, 
# then increase labor_prod by one, iterate first all assets, then Θ etc

sum_current_workers_by_assets=zeros(Float64,size(agrid_fine))
for i=1:size(agrid_fine)[1]
    sum_current_workers_by_assets[i] = sum(current_workers[i:size(agrid_fine)[1]:(no_prod_shock^2-1)*((size(agrid_fine))[1])+i])/current_worker_pop
end
plot(sum_current_workers_by_assets)

sum_current_staple_by_assets = zeros(Float64, size(agrid_fine))
for i = 1:size(agrid_fine)[1]
    sum_current_staple_by_assets[i] = sum(current_staple[i:size(agrid_fine)[1]:(no_prod_shock^2-1)*((size(agrid_fine))[1])+i]) / current_staple_pop
end
plot(sum_current_staple_by_assets)


sum_current_cashcrop_by_assets = zeros(Float64, size(agrid_fine))
for i = 1:size(agrid_fine)[1]
    sum_current_cashcrop_by_assets[i] = sum(current_cashcrop[i:size(agrid_fine)[1]:(no_prod_shock^2-1)*((size(agrid_fine))[1])+i]) / current_cashcrop_pop
end
plot(sum_current_cashcrop_by_assets)

sum_current_staple_by_prod = zeros(Float64, no_prod_shock)
for i = 1:no_prod_shock
    tmp_vec = sum(reshape(current_staple, size(agrid_fine)[1], no_prod_shock^2), dims=1)
    sum_current_staple_by_prod[i] = sum(tmp_vec[i:no_prod_shock:(no_prod_shock^2-no_prod_shock)+i]) / (1-current_worker_pop)
end
plot(sum_current_staple_by_prod)


sum_current_cashcrop_by_prod = zeros(Float64, no_prod_shock)
for i = 1:no_prod_shock
    tmp_vec = sum(reshape(current_cashcrop, size(agrid_fine)[1], no_prod_shock^2), dims=1)
    sum_current_cashcrop_by_prod[i] = sum(tmp_vec[i:no_prod_shock:(no_prod_shock^2-no_prod_shock)+i]) / (1-current_worker_pop)
end
plot!(sum_current_cashcrop_by_prod)

sum_current_worker_by_prod = zeros(Float64, no_prod_shock)
for i = 1:no_prod_shock
    tmp_vec = sum(reshape(current_workers, size(agrid_fine)[1], no_prod_shock^2), dims=1)
    sum_current_worker_by_prod[i] = sum(tmp_vec[i:no_prod_shock:(no_prod_shock^2-no_prod_shock)+i]) / (current_worker_pop)
end
plot!(sum_current_worker_by_prod)

plot(sum(reshape(current_cashcrop, 16, 160), dims=2))
plot!(sum(reshape(current_staple, 16, 160), dims=2))

plot(τ_grid, 10 * program_spending_subsidy_b_grid)

plot(τ_grid, 10 * mig_rate_model_subsidy_b_grid)

plot(τ_grid, prod_staple_subsidy_b_grid)
plot!(τ_grid, prod_cashcrop_subsidy_b_grid)

input_staple_subsidy_b_grid_pc = input_staple_subsidy_b_grid ./ (current_staple_pop_subsidy_b_grid + current_cashcrop_pop_subsidy_b_grid.*mean_land_share_to_staples_among_cc_model_subsidy_b_grid)
input_cashcrop_subsidy_b_grid_pc = input_cashcrop_subsidy_b_grid ./ (current_cashcrop_pop_subsidy_b_grid.*(1 .- mean_land_share_to_staples_among_cc_model_subsidy_b_grid))

effective_no_staple_prods=(current_staple_pop_subsidy_b_grid + current_cashcrop_pop_subsidy_b_grid .* mean_land_share_to_staples_among_cc_model_subsidy_b_grid)
effective_no_cashcrop_prods = (current_cashcrop_pop_subsidy_b_grid .*( 1 .- mean_land_share_to_staples_among_cc_model_subsidy_b_grid))


plot(τ_grid, prod_staple_subsidy_b_grid)
plot!(τ_grid, prod_cashcrop_subsidy_b_grid)

plot(τ_grid, prod_staple_subsidy_b_grid_pc)
plot!(τ_grid, prod_cashcrop_subsidy_b_grid_pc)

plot(τ_grid, 10*staple_productivity_subsidy_b_grid)
plot!(τ_grid, 10*cashcrop_productivity_subsidy_b_grid)
plot!(τ_grid, manuf_productivity_subsidy_b_grid)
plot!(τ_grid, 10 * APG_subsidy_b_grid)

plot(τ_grid, mean_land_share_to_staples_among_cc_model_subsidy_b_grid)

plot(τ_grid, input_staple_subsidy_b_grid)
plot!(τ_grid, input_cashcrop_subsidy_b_grid)

plot(τ_grid, input_staple_subsidy_b_grid_pc)
plot!(τ_grid, input_cashcrop_subsidy_b_grid_pc)

plot(τ_grid, effective_no_staple_prods)
plot!(τ_grid, effective_no_cashcrop_prods)

plot(τ_grid, rural_pop_only_staples_model_subsidy_b_grid)
plot(τ_grid, sum(stat_distr_subsidy_b_grid .* reshape(a_prime_fine_local_subsidy_b_grid, 7680, 21), dims=1))
plot(τ_grid, p_x .* input_staple_subsidy_b_grid)
plot(τ_grid, p_x .* input_cashcrop_subsidy_b_grid)


plot(τ_grid, 10 * current_worker_pop_subsidy_b_grid)
plot!(τ_grid, 10 * current_cashcrop_pop_subsidy_b_grid)
plot!(τ_grid, 10 * current_staple_pop_subsidy_b_grid)




plot(τ_grid, prices_subsidy_b_grid[1, :])
plot!(τ_grid, 10 * prices_subsidy_b_grid[2, :])
plot!(τ_grid, 10 * prices_subsidy_b_grid[3, :])

plot(τ_grid, 10 * aggregate_consumption_subsidy_b_grid)



prices_subsidy_b_grid_smooth = copy(prices_subsidy_b_grid);
prices_subsidy_b_grid_smooth[1,1:(end-1)] = movmean(prices_subsidy_b_grid[1,1:(end-1)],7);
prices_subsidy_b_grid_smooth[2,1:(end-1)] = movmean(prices_subsidy_b_grid[2,1:(end-1)],7);
prices_subsidy_b_grid_smooth[3,1:(end-1)] = movmean(prices_subsidy_b_grid[3,1:(end-1)],7);


# 1st period of transition 
Phi_fine_aug = kron(Matrix(1.0I, 3, 3),row_kron(funbase(fspace_a,s_fine[:,1]),Phi_z_fine));
welfare_val_b_constant_no_subsidy_price = (Phi_fine_aug * vcat(vcat(coeff[:,1],coeff[:,2]),coeff[:,3]));
welfare_subsidy_constant_no_subsidy_price1 = sum(stat_distr_no_subsidy  .* (exp.((welfare_val_b_constant_no_subsidy_price -welfare_val_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1; # Without balanced budget
welfare_subsidy_constant_no_subsidy_price2 = sum(stat_distr  .* (exp.((welfare_val_b_constant_no_subsidy_price -welfare_val_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1; # Without balanced budget
V_saved_constant_no_subsidy_price_reshaped=reshape(V_saved_local,Baseline_parameter.ns_fine*3)
welfare_subsidy_b_real_constant =sum(stat_distr_no_subsidy  .* (exp.((V_saved_constant_no_subsidy_price_reshaped -V_saved_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget

# 1st period of transition no QS
Phi_fine_aug = kron(Matrix(1.0I, 3, 3),row_kron(funbase(fspace_a,s_fine[:,1]),Phi_z_fine));
welfare_val_b_constant_no_QS_no_subsidy_price = (Phi_fine_aug * vcat(vcat(coeff[:,1],coeff[:,2]),coeff[:,3]));
welfare_subsidy_constant_no_QS_no_subsidy_price1 = sum(stat_distr_no_QS_no_subsidy  .* (exp.((welfare_val_b_constant_no_subsidy_price -welfare_val_no_QS_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1; # Without balanced budget
welfare_subsidy_constant_no_QS_no_subsidy_price2 = sum(stat_distr  .* (exp.((welfare_val_b_constant_no_subsidy_price -welfare_val_no_QS_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1; # Without balanced budget
V_saved_constant_no_QS_no_subsidy_price_reshaped=reshape(V_saved_local,Baseline_parameter.ns_fine*3)
welfare_subsidy_b_noQS_real_constant =sum(stat_distr_no_QS_no_subsidy  .* (exp.((V_saved_constant_no_QS_no_subsidy_price_reshaped -V_saved_no_QS_no_subsidy_reshaped) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget

# Get the consumption production graph for a given state

fixed_state_index = 366;
Y_W_fixedst = zeros(ns_fine,3,n_fine[1]);
Y_S_fixedst = zeros(ns_fine,3,n_fine[1]);
c_S_S_fixedst= zeros(ns_fine,3,n_fine[1]);
c_B_S_fixedst= zeros(ns_fine,3,n_fine[1]);
c_M_S_fixedst= zeros(ns_fine,3,n_fine[1]);
q_S_S_fixedst= zeros(ns_fine,3,n_fine[1]);
P_S_fixedst= zeros(ns_fine,3,n_fine[1]);
x_S_S_fixedst= zeros(ns_fine,3,n_fine[1]);
solve_staple_index_S_fixedst= zeros(Int64,ns_fine,3,n_fine[1]);
λ_2_S_fixedst= zeros(ns_fine,3,n_fine[1]);
future_asset_S_fixedst= zeros(ns_fine,3,n_fine[1]);
future_asset_C_fixedst= zeros(ns_fine,3,n_fine[1]);
c_S_B_fixedst= zeros(ns_fine,3,n_fine[1]);
c_B_B_fixedst= zeros(ns_fine,3,n_fine[1]);
c_M_B_fixedst= zeros(ns_fine,3,n_fine[1]);
x_SC_fixedst= zeros(ns_fine,3,n_fine[1]);
x_BC_fixedst= zeros(ns_fine,3,n_fine[1]);
land_C_fixedst= zeros(ns_fine,3,n_fine[1]);
λ_2_fixedst= zeros(ns_fine,3,n_fine[1]);
P_B_fixedst= zeros(ns_fine,3,n_fine[1]);
Y_B_fixedst= zeros(ns_fine,3,n_fine[1]);
q_S_B_fixedst= zeros(ns_fine,3,n_fine[1]);
q_B_B_fixedst= zeros(ns_fine,3,n_fine[1]);
solve_cash_crop_index_B_fixedst= zeros(Int64,ns_fine,3,n_fine[1]);
solve_staple_index_B_fixedst= zeros(Int64,ns_fine,3,n_fine[1]);
TC_fixedst= zeros(ns_fine,3,n_fine[1]);

for ii=1:n_fine[1]
s_fixedst = copy(s_fine);
s_fixedst[:,1] .= s_fine[fixed_state_index,1] 
s_fixedst[:,2] .= s_fine[fixed_state_index,2]
cons_fixedst = C_grid_fine[ii] * ones(n_fine[1] * n_fine[2])


for j = 1:3
    for jj =1:3
        if jj == 1
            Y_W_fixedst[:,j,ii] =  P_W *cons_fixedst  +Y_W_fine  .+ w*FM_W .+ w*F_W * (j != 1);
        end
        if jj == 2
        (future_asset_S_fixedst[:,j,ii],Y_S_fixedst[:,j,ii],c_S_S_fixedst[:,j,ii],c_B_S_fixedst[:,j,ii],c_M_S_fixedst[:,j,ii],q_S_S_fixedst[:,j,ii],
        P_S_fixedst[:,j,ii],x_S_S_fixedst[:,j,ii],solve_staple_index_S_fixedst[:,j,ii],λ_2_S_fixedst[:,j,ii]) = policy_function_creator(
            cons_fixedst,jj,j,P_W_fine,Y_W_fine,s_fixedst,r,ρ,w,coeff_λ_2_cashcrop_residual_unconstrained_fine,
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
                (future_asset_C_fixedst[:,j,ii],c_S_B_fixedst[:,j,ii],c_B_B_fixedst[:,j,ii],c_M_B_fixedst[:,j,ii],x_SC_fixedst[:,j,ii],x_BC_fixedst[:,j,ii],
                land_C_fixedst[:,j,ii],λ_2_fixedst[:,j,ii],P_B_fixedst[:,j,ii],Y_B_fixedst[:,j,ii],q_S_B_fixedst[:,j,ii],q_B_B_fixedst[:,j,ii],solve_cash_crop_index_B_fixedst[:,j,ii]
                ,solve_staple_index_B_fixedst[:,j,ii],TC_fixedst[:,j,ii]) = policy_function_creator(
                    cons_fixedst,jj,j,P_W_fine,Y_W_fine,s_fixedst,r,ρ,w,coeff_λ_2_cashcrop_residual_unconstrained_fine,
                    coeff_λ_2_cashcrop_residual_constrained_fine,θ_fine,fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,
                    p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,a_min,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
                    coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,
                    labor_allocated_interior_c3a_fine,λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,
                    q_S_c2_fine,q_B_c1_fine,q_S_c3a_fine,q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,
                    Y_B_c3b_fine,c_S_mat_fine,c_B_mat_fine,c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,q_B_mat_fine,land_B_mat_fine,λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,
                    feasibility_mat_fine,C_max_mat_fine,C_min_mat_fine,q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,x_S_staples_fine,
                    λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat_fine,C_max_staple_fine,C_min_staple_fine,C_max_staple_constrained_fine,
                    C_min_staple_constrained_fine,TC_S_c3_constrained_fine,x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine,
                    x_S_mat_3c_fine,x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine)
            end
        end
end
end



println("\\hline  & FISP & FISP \\\\ ")
println("&     tax-financed \$ \\tau_W = ",  convert(Int64,round(100* prices_subsidy_b[3])), "\\%\$  & \$ \\tau_W = ",  convert(Int64,round(100* prices_no_QS_subsidy_b[3])), "\\%\$\\\\ ")
println("\\hline \\\\[-1.8ex] ")
println("\\alert<1>{Consumption equivalent welfare}  &  ",round(100 * welfare_subsidy_b,digits = 1), "\\% &", round(100 * welfare_subsidy_b_noQS,digits = 1) ,"\\% \\\\ ")
println("Prices: \$p_B/p_M/w\$ & \$+",convert(Int64, round(100 * prices_subsidy_b[1]/prices_no_subsidy[1])) -100,"\\%/"
,"+",convert(Int64, round(100 * prices_subsidy_b[2]/prices_no_subsidy[2]))-100,"\\%/"
,"+",convert(Int64, round(100 * w_subsidy_b/w_no_subsidy))-100,"\\%\$  & \$+",convert(Int64, round(100 * prices_no_QS_subsidy_b[1]/prices_no_QS_no_subsidy[1]))-100,"\\%/"
,"+",convert(Int64, round(100 * prices_no_QS_subsidy_b[2]/prices_no_QS_no_subsidy[2]))-100,"\\%/"
,"+",convert(Int64, round(100 * w_no_QS_subsidy_b/w_no_QS_no_subsidy))-100,"\\%\$\\\\")
#println("Prices: \$p_B/p_M/w\$ & ", round(prices_no_subsidy[1],digits = 1),"/",round(prices_no_subsidy[2],digits = 1),"&+",convert(Int64, round(100 * prices_subsidy_b[1]/prices_no_subsidy[1])) -100,"\\%/"
#,"+",convert(Int64, round(100 * prices_subsidy_b[2]/prices_no_subsidy[2]))-100,"\\% & +" ,convert(Int64, round(100 * prices_no_QS_subsidy_b[1]/prices_no_subsidy[1]))-100,"\\%/"
#,"+",convert(Int64, round(100 * prices_no_QS_subsidy_b[2]/prices_no_subsidy[2]))-100,"\\%\\\\")
println("Share of staple-only farmers  &+",convert(Int64,round(100 * rural_pop_only_staples_model_subsidy_b/rural_pop_only_staples_model_no_subsidy))-100, "\\% & ",convert(Int64,round(100 * rural_pop_only_staples_model_no_QS_subsidy_b/rural_pop_only_staples_model_no_subsidy))-100,"\\% \\\\")
println("\\alert<2>{Staple production}  & +", convert(Int64, round(100 * prod_staple_subsidy_b/prod_staple_no_subsidy)) -100 ,"\\%  &+", convert(Int64, round(100 * prod_staple_no_QS_subsidy_b/prod_staple_no_QS_no_subsidy)) -100,"\\% \\\\")
println("\\alert<2>{Staple productivity} & +", convert(Int64, round(100 * staple_productivity_subsidy_b/staple_productivity_no_subsidy)) -100,"\\% &+", convert(Int64, round(100 * staple_productivity_no_QS_subsidy_b/staple_productivity_no_QS_no_subsidy)) -100,"\\% \\\\")
println("\\alert<3>{Cash crop production} & +", convert(Int64, round(100 * prod_cashcrop_subsidy_b/prod_cashcrop_no_subsidy)) -100,"\\%  & +", convert(Int64, round(100 * prod_cashcrop_no_QS_subsidy_b/prod_cashcrop_no_QS_no_subsidy))-100 ,"\\% \\\\")
println("\\alert<3>{Cash crop productivity} & +", convert(Int64, round(100 * cashcrop_productivity_subsidy_b/cashcrop_productivity_no_subsidy))-100 ,"\\% & +", convert(Int64, round(100 * cashcrop_productivity_no_QS_subsidy_b/cashcrop_productivity_no_QS_no_subsidy))-100 ,"\\% \\\\")
println("\\alert<4>{Share of financially constrained farmers} &",convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_subsidy_b/fraction_cashcrop_suboptimal_model_no_subsidy))-100, "\\% &+",convert(Int64,round(100 * fraction_cashcrop_suboptimal_model_no_QS_subsidy_b/fraction_cashcrop_suboptimal_model_no_QS_no_subsidy))-100,"\\% \\\\")
println("\\alert<5>{Urbanization rate} & ",convert(Int64,round(100 * current_worker_pop_subsidy_b/current_worker_pop_no_subsidy)) -100, "\\% &+",convert(Int64,round(100 * current_worker_pop_no_QS_subsidy_b/current_worker_pop_no_QS_no_subsidy))-100,"\\% \\\\")
println("APG & +",convert(Int64,round(APG_subsidy_b/APG_no_subsidy*100,digits = 0) -100), "\\%  &",convert(Int64,round(APG_no_QS_subsidy_b/APG_no_QS_no_subsidy*100,digits = 0) - 100),"\\%  \\\\")
println("Consumption & ", convert(Int64, round(100 * aggregate_consumption_subsidy_b/aggregate_consumption_no_subsidy)) -100,
"\\%  &+", convert(Int64, round(100 * aggregate_consumption_no_QS_subsidy_b/aggregate_consumption_no_QS_no_subsidy)) -100,"\\% \\\\")
println("Transaction cost & +", convert(Int64, round(100 * transaction_cost_loss_subsidy_b/transaction_cost_loss_no_subsidy)) -100,
"\\%  & - \\\\")
println("Current account surplus \\% of GDP & ",convert(Int64,round(100*current_account_residual_subsidy_b,digits = 0)) - convert(Int64,round(100*current_account_residual_no_subsidy,digits = 0)),
"\\%&",convert(Int64,round(100*current_account_residual_no_QS_subsidy_b,digits = 0))- convert(Int64,round(100*current_account_residual_no_QS_no_subsidy,digits = 0)),"\\% \\\\")
println("Fraction of farmers without surplus  & ",convert(Int64,round(100*fraction_model_subsidy_b/fraction_model_no_subsidy,digits = 0))-100,
"\\%&",convert(Int64,round(100*fraction_model_no_QS_subsidy_b/fraction_model_no_QS_no_subsidy,digits = 0))-100,"\\% \\\\")
println("Program cost \\% of GDP & ",convert(Int64,round(100*program_spending_subsidy_b,digits = 0)),
"\\%&",convert(Int64,round(100*program_spending_no_QS_subsidy_b,digits = 0)),"\\% \\\\")
println("Migration rate & ",convert(Int64,round(100*mig_rate_model_subsidy_b/mig_rate_model_no_subsidy,digits = 0)) -100,
"\\%&+",convert(Int64,round(100*mig_rate_model_no_QS_subsidy_b/mig_rate_model_no_QS_no_subsidy,digits = 0))-100,"\\% \\\\")


# With keeping prices constant its not much different:
[APG_no_subsidy,APG_nb_cp,APG_b_cp]
table_3_results[11,:] = [11,1,nominal_GDP_subsidy_nb/nominal_GDP_no_subsidy,nominal_GDP_subsidy_b/nominal_GDP_no_subsidy,
nominal_GDP_inf/nominal_GDP_no_subsidy,nominal_GDP_inf_sp/nominal_GDP_no_subsidy];
table_3_results[12,:] = [12,1,aggregate_consumption_subsidy_nb/aggregate_consumption_no_subsidy,aggregate_consumption_subsidy_b/aggregate_consumption_no_subsidy
,aggregate_consumption_inf/aggregate_consumption_no_subsidy,aggregate_consumption_inf_sp/aggregate_consumption_no_subsidy];
table_3_results[13,:] = [13,1,transaction_cost_loss_subsidy_nb/transaction_cost_loss_no_subsidy,transaction_cost_loss_subsidy_b/transaction_cost_loss_no_subsidy
,transaction_cost_loss_inf/transaction_cost_loss_no_subsidy,transaction_cost_loss_inf_sp/transaction_cost_loss_no_subsidy];
table_3_results[14,:] = [14, current_account_residual_no_subsidy/nominal_GDP_no_subsidy,current_account_residual_subsidy_nb/nominal_GDP_subsidy_nb,0
,current_account_residual_inf/nominal_GDP_subsidy_nb,current_account_residual_inf_sp/nominal_GDP_subsidy_nb]


plot(100*[movmean(food_share_staple_subsidy_b_grid_consbased[:,21],160),
movmean(food_share_staple_subsidy_b_grid_consbased[:,11],160),movmean(food_share_staple_subsidy_b_grid_consbased[:,1],160)]
, label=[ "\\tau_S = 81 %" "\\tau_S = 41 %" "\\tau_S = 0 %"  ],legend=:outerbottom,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [0.0,80.0], xlabel = "Consumption bundle",
ylabel = "Percent of staples to consumption",xticks = [],xlims = [0,2000],
grid = false,size = (800,800),marker = [:none :circle :none],
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure3h.svg")

plot(100*[food_share_worker_subsidy_b_grid_consbased[:,21],
food_share_worker_subsidy_b_grid_consbased[:,11],food_share_worker_subsidy_b_grid_consbased[:,1]]
, label=[ "\\tau_S = 81 %" "\\tau_S = 41 %" "\\tau_S = 0 %"  ],legend=:outerbottom,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [0.0,40.0], xlabel = "Consumption bundle",
ylabel = "Percent of staples to consumption",xticks = [],xlims = [0,2000],
grid = false,size = (800,800),marker = [:none :circle :none],
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure3i.svg")

# Plots.contour(100*food_share_worker_subsidy_b_grid_consbased', fill=true,levels = 5,
#  c =cgrad(:Blues_6, rev = true),xlabel = "Consumption bundle",
#  ylabel = "Percent of staples to consumption",axis=nothing,tickfontsize = 14, ylims = [0.0,20.0], 
#  xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="Times_Roman")
# plot(100*[manuf_share_worker_subsidy_b_grid_consbased[:,21],
# manuf_share_worker_subsidy_b_grid_consbased[:,11],manuf_share_worker_subsidy_b_grid_consbased[:,1]]
# , label=[ "High subsidy" "Medium subsidy" "No subsidy"  ],legend=:outerbottom,
# linewidth = 2,linestyle = [:solid :dash :dashdot],ylims = [40.0,80.0], xlabel = "Consumption bundle",
# ylabel = "Percent of manufacturing consumption expenditure to consumption",
# grid = false,size = (800,800),
# tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="Times_Roman")
# savefig("weird.pdf") doesnt look that interesting


# Welfare plots:


