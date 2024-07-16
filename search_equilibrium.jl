include("setup.jl")

##########################################################################################################################################################
##
##          FINDING EQUILIBRIUM PRICES:
##
##########################################################################################################################################################



# Search for solutions - skip this part for only the result
# Subsidy with balanced budget
balanced_share = 1.0
parameters_tmp = copy(Baseline_parameter);
prices = [ 1.518014329164464,
0.2726473607850526,
0.1886863641824446];#,2.505662787966534

println("Solving the model ...")
#out=1/2 - choose 1 for solver, 2 for analysis
function f_solve(prices)
    residual = solve_model_calibration1(prices,parameters_tmp,1,moments,balanced_share);
    return residual
end

function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;

stst_simplex = Optim.AffineSimplexer(0.025,0.05);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =30000.0);
prices = optim_res.minimizer;
prices_subsidy_b = copy(prices);
#residual = f_solve(prices);
(residual_subsidy_b, stat_distr_subsidy_b, cons_fine_local_subsidy_b, future_occupation_fine_local_subsidy_b,x_SC_fine_subsidy_b,x_BC_fine_subsidy_b, 
    coeff_subsidy_b,residual_goods_subsidy_b,model_moments_subsidy_b,foreign_supply_capital_subsidy_b)  = solve_model_calibration1(prices_subsidy_b,parameters_tmp,2,moments,balanced_share);

# Subsidy without balanced budget
balanced_share = 0.0
parameters_tmp = copy(Baseline_parameter);
prices = [ 1.5087253977851676,
0.25456229359350663];
#foreign_supply_capital_subsidy_b = 7.628070905421486; # through the calibration of subsidy_b
println("Solving the model ...")
#out=1/2 - choose 1 for solver, 2 for analysis
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end

function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;

stst_simplex = Optim.AffineSimplexer(0.025,0.05);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =10000.0);
prices = optim_res.minimizer;
prices_subsidy_nb = copy(prices);
#residual = f_solve(prices);
(residual_subsidy_nb, stat_distr_subsidy_nb, cons_fine_local_subsidy_nb, future_occupation_fine_local_subsidy_nb,x_SC_fine_subsidy_nb,x_BC_fine_subsidy_nb, 
    coeff_subsidy_nb,residual_goods_subsidy_nb,model_moments_subsidy_nb,foreign_supply_capital_subsidy_nb)  = solve_model_calibration2(prices_subsidy_nb,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);

#No subsidy equilibrium
No_subsidy_parameter = copy(Baseline_parameter);
No_subsidy_parameter.τ_S=-0.0;
parameters_tmp = copy(No_subsidy_parameter);
balanced_share = 0.0
#foreign_supply_capital_subsidy_b = 10.254571759208405;
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end

prices = [1.1882831767847994,
0.18349666056314162];
function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.5;
upper_guess_bound = 2.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.5,0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =12000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_subsidy = copy(prices);

(residual_no_subsidy, stat_distr_no_subsidy, cons_fine_local_no_subsidy, future_occupation_fine_local_no_subsidy,x_SC_fine_no_subsidy,x_BC_fine_no_subsidy, 
    coeff_no_subsidy,residual_goods_no_subsidy,model_moments_no_subsidy,foreign_supply_capital_no_subsidy)  = solve_model_calibration2(prices_no_subsidy,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);

########################
### NEED TO SOLVE HERE NO SUB + SUB EQM UNDER epsilon_u/r >0 (say 0.5 , 1.0, 1.5 for start?)
########################

# Undernourishment requires additional calculations
# Staple consumption w a threshold value -
# assume the threshold is chosen such that the calibration hits the average between the 17.9% (Pauw, Beck, and Mussa (2014) )
# and 24.5 (Malawi NSO) extreme food poverty rate of HHs. for now, its 20% in 2010 in Malawi
cons_level_substinence = 0.109;
foreign_supply_capital_subsidy_b = 7.631597488627073;
# Subsidy with balanced budget under all values of epsilon_r/u  == subsidy w balanced budget baseline.
#No subsidy equilibrium - epsilon_r/u 0.1155 - Strauss implied numbers
epsilon_r=0.1155
epsilon_u=0.1155
No_subsidy_parameter = copy(Baseline_parameter);
No_subsidy_parameter.τ_S = -0.0;
parameters_tmp = copy(No_subsidy_parameter);
balanced_share = 0.0
#foreign_supply_capital_subsidy_b = 10.254571759208405;
function f_solve(prices)
    residual = solve_model_calibration_externality(prices, parameters_tmp, 1, moments, balanced_share, foreign_supply_capital_subsidy_b, cons_level_substinence, epsilon_u, epsilon_r)
    return residual
end

prices = [ 1.1925825472322047,
0.18727339112307217,
0.2848214354714461];
function f_sum(prices)
    residual = f_solve(prices)
    return sum(residual .^ 2)
end
lower_guess_bound = 0.5;
upper_guess_bound = 2.0;
ls_res = LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(), show_trace=true, store_trace=true,
    x_tol=1e-9, f_tol=1e-5, iterations=20, lower=lower_guess_bound * prices,
    upper=upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.5, 0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method=NelderMead(initial_simplex=stst_simplex), store_trace=true, show_trace=true,
    extended_trace=true, time_limit=12000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_subsidyEps_Strauss = copy(prices);

(residual_no_subsidyEps_Strauss, stat_distr_no_subsidyEps_Strauss, cons_fine_local_no_subsidyEps_Strauss, future_occupation_fine_local_no_subsidyEps_Strauss, x_SC_fine_no_subsidyEps_Strauss, x_BC_fine_no_subsidyEps_Strauss,
    coeff_no_subsidyEps_Strauss, residual_goods_no_subsidyEps_Strauss, model_moments_no_subsidyEps_Strauss, foreign_supply_capital_no_subsidyEps_Strauss) = solve_model_calibration_externality(prices_no_subsidyEps_Strauss, parameters_tmp, 2, moments, balanced_share, foreign_supply_capital_subsidy_b, cons_level_substinence, epsilon_u, epsilon_r);


#No subsidy equilibrium - epsilon_r/u 0.25
epsilon_r=0.25
epsilon_u=0.25
No_subsidy_parameter = copy(Baseline_parameter);
No_subsidy_parameter.τ_S = -0.0;
parameters_tmp = copy(No_subsidy_parameter);
balanced_share = 0.0
#foreign_supply_capital_subsidy_b = 10.254571759208405;
function f_solve(prices)
    residual = solve_model_calibration_externality(prices, parameters_tmp, 1, moments, balanced_share, foreign_supply_capital_subsidy_b, cons_level_substinence, epsilon_u, epsilon_r)
    return residual
end

prices = [ 1.1925825472322047,
0.18727339112307217,
0.2848214354714461];
function f_sum(prices)
    residual = f_solve(prices)
    return sum(residual .^ 2)
end
lower_guess_bound = 0.5;
upper_guess_bound = 2.0;
ls_res = LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(), show_trace=true, store_trace=true,
    x_tol=1e-9, f_tol=1e-5, iterations=20, lower=lower_guess_bound * prices,
    upper=upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.5, 0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method=NelderMead(initial_simplex=stst_simplex), store_trace=true, show_trace=true,
    extended_trace=true, time_limit=12000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_subsidyEps025 = copy(prices);

(residual_no_subsidyEps025, stat_distr_no_subsidyEps025, cons_fine_local_no_subsidyEps025, future_occupation_fine_local_no_subsidyEps025, x_SC_fine_no_subsidyEps025, x_BC_fine_no_subsidyEps025,
    coeff_no_subsidyEps025, residual_goods_no_subsidyEps025, model_moments_no_subsidyEps025, foreign_supply_capital_no_subsidyEps025) = solve_model_calibration_externality(prices_no_subsidyEps025, parameters_tmp, 2, moments, balanced_share, foreign_supply_capital_subsidy_b, cons_level_substinence, epsilon_u, epsilon_r);



#No subsidy equilibrium - epsilon_r/u 0.5
epsilon_r=0.5
epsilon_u=0.5
No_subsidy_parameter = copy(Baseline_parameter);
No_subsidy_parameter.τ_S = -0.0;
parameters_tmp = copy(No_subsidy_parameter);
balanced_share = 0.0
#foreign_supply_capital_subsidy_b = 10.254571759208405;
function f_solve(prices)
    residual = solve_model_calibration_externality(prices, parameters_tmp, 1, moments, balanced_share, foreign_supply_capital_subsidy_b, cons_level_substinence, epsilon_u, epsilon_r)
    return residual
end

prices = [ 1.195476951761732,
0.20085571586329348,
0.28334980569254575];
function f_sum(prices)
    residual = f_solve(prices)
    return sum(residual .^ 2)
end
lower_guess_bound = 0.5;
upper_guess_bound = 2.0;
ls_res = LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(), show_trace=true, store_trace=true,
    x_tol=1e-9, f_tol=1e-5, iterations=20, lower=lower_guess_bound * prices,
    upper=upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.5, 0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method=NelderMead(initial_simplex=stst_simplex), store_trace=true, show_trace=true,
    extended_trace=true, time_limit=12000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_subsidyEps05 = copy(prices);

(residual_no_subsidyEps05, stat_distr_no_subsidyEps05, cons_fine_local_no_subsidyEps05, future_occupation_fine_local_no_subsidyEps05, x_SC_fine_no_subsidyEps05, x_BC_fine_no_subsidyEps05,
    coeff_no_subsidyEps05, residual_goods_no_subsidyEps05, model_moments_no_subsidyEps05, foreign_supply_capital_no_subsidyEps05) = solve_model_calibration_externality(prices_no_subsidyEps05, parameters_tmp, 2, moments, balanced_share, foreign_supply_capital_subsidy_b, cons_level_substinence, epsilon_u, epsilon_r);



#No subsidy equilibrium - epsilon_r/u 1.0
epsilon_r = 1.0
epsilon_u = 1.0
No_subsidy_parameter = copy(Baseline_parameter);
No_subsidy_parameter.τ_S = -0.0;
parameters_tmp = copy(No_subsidy_parameter);
balanced_share = 0.0
#foreign_supply_capital_subsidy_b = 10.254571759208405;
function f_solve(prices)
    residual = solve_model_calibration_externality(prices, parameters_tmp, 1, moments, balanced_share, foreign_supply_capital_subsidy_b, cons_level_substinence, epsilon_u, epsilon_r)
    return residual
end

prices = [ 1.2145275998907963,
0.2147096472250037,
0.2936424726929504];
function f_sum(prices)
    residual = f_solve(prices)
    return sum(residual .^ 2)
end
lower_guess_bound = 0.5;
upper_guess_bound = 2.0;
ls_res = LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(), show_trace=true, store_trace=true,
    x_tol=1e-9, f_tol=1e-5, iterations=20, lower=lower_guess_bound * prices,
    upper=upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.1, 0.05);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method=NelderMead(initial_simplex=stst_simplex), store_trace=true, show_trace=true,
    extended_trace=true, time_limit=12000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_subsidyEps10 = copy(prices);

(residual_no_subsidyEps10, stat_distr_no_subsidyEps10, cons_fine_local_no_subsidyEps10, future_occupation_fine_local_no_subsidyEps10, x_SC_fine_no_subsidyEps10, x_BC_fine_no_subsidyEps10,
    coeff_no_subsidyEps10, residual_goods_no_subsidyEps10, model_moments_no_subsidyEps10, foreign_supply_capital_no_subsidyEps10) = solve_model_calibration_externality(prices_no_subsidyEps10, parameters_tmp, 2, moments, balanced_share, foreign_supply_capital_subsidy_b, cons_level_substinence, epsilon_u, epsilon_r);



#No subsidy equilibrium - epsilon_r/u 1.5
epsilon_r = 1.5
epsilon_u = 1.5
No_subsidy_parameter = copy(Baseline_parameter);
No_subsidy_parameter.τ_S = -0.0;
parameters_tmp = copy(No_subsidy_parameter);
balanced_share = 0.0
#foreign_supply_capital_subsidy_b = 10.254571759208405;
function f_solve(prices)
    residual = solve_model_calibration_externality(prices, parameters_tmp, 1, moments, balanced_share, foreign_supply_capital_subsidy_b, cons_level_substinence, epsilon_u, epsilon_r)
    return residual
end

prices = [  1.2093001817161375,
0.2130922258686436,
0.2979426218123759];
function f_sum(prices)
    residual = f_solve(prices)
    return sum(residual .^ 2)
end
lower_guess_bound = 0.5;
upper_guess_bound = 2.0;
ls_res = LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(), show_trace=true, store_trace=true,
    x_tol=1e-9, f_tol=1e-5, iterations=20, lower=lower_guess_bound * prices,
    upper=upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.5, 0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method=NelderMead(initial_simplex=stst_simplex), store_trace=true, show_trace=true,
    extended_trace=true, time_limit=12000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_subsidyEps15 = copy(prices);

(residual_no_subsidyEps15, stat_distr_no_subsidyEps15, cons_fine_local_no_subsidyEps15, future_occupation_fine_local_no_subsidyEps15, x_SC_fine_no_subsidyEps15, x_BC_fine_no_subsidyEps15,
    coeff_no_subsidyEps15, residual_goods_no_subsidyEps15, model_moments_no_subsidyEps15, foreign_supply_capital_no_subsidyEps15) = solve_model_calibration_externality(prices_no_subsidyEps15, parameters_tmp, 2, moments, balanced_share, foreign_supply_capital_subsidy_b, cons_level_substinence, epsilon_u, epsilon_r);


### Robustness analysis - contribution of each friction

# κ = 10; (make sure that the fraction of constrained farms is very small) with subsidy
# Given the inaccuracies, we might want to increase the interest rate
no_κ_parameter = copy(Baseline_parameter);
no_κ_parameter.κ= 10.0;
parameters_tmp = copy(no_κ_parameter);
balanced_share = 1.0
#foreign_supply_capital_subsidy_b = 7.628070905421486;
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end
prices = [  1.5140725342775192,
0.280305851412856,
0.26174075283278103]; #!
# prices = 
#prices_bug = [  1.0012488427964574,1.304207153437122,0.06502460876457435,1.0524640422259952];
#prices = copy(prices_bug)
function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.05,0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =32000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_κ_subsidy_b = copy(prices);

(residual_no_κ_subsidy_b, stat_distr_no_κ_subsidy_b, cons_fine_local_no_κ_subsidy_b, future_occupation_fine_local_no_κ_subsidy_b,x_SC_fine_no_κ_subsidy_b,x_BC_fine_no_κ_subsidy_b, 
    coeff_no_κ_subsidy_b,residual_goods_no_κ_subsidy_b,model_moments_no_κ_subsidy_b,foreign_supply_capital_no_κ_subsidy_b)  = solve_model_calibration2(prices_no_κ_subsidy_b,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);


# κ = 10; (make sure that the fraction of constrained farms is very small) with no subsidy
no_κ_parameter_no_subsidy = copy(Baseline_parameter);
no_κ_parameter_no_subsidy.κ= 10.0;
no_κ_parameter_no_subsidy.τ_S=-0.0;
parameters_tmp = copy(no_κ_parameter_no_subsidy);
balanced_share = 0.0
#foreign_supply_capital_subsidy_b = 7.628070905421486;
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end
prices = [ 1.1696929055108973,
0.1890191391243493]#!
# prices = 
#prices_bug = [  1.0012488427964574,1.304207153437122,0.06502460876457435,1.0524640422259952];
#prices = copy(prices_bug)
function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.05,0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =12000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_κ_no_subsidy = copy(prices);

(residual_no_κ_no_subsidy, stat_distr_no_κ_no_subsidy, cons_fine_local_no_κ_no_subsidy, future_occupation_fine_local_no_κ_no_subsidy,x_SC_fine_no_κ_no_subsidy,x_BC_fine_no_κ_no_subsidy, 
    coeff_no_κ_no_subsidy,residual_goods_no_κ_no_subsidy,model_moments_no_κ_no_subsidy,foreign_supply_capital_no_κ_no_subsidy)  = solve_model_calibration2(prices_no_κ_no_subsidy,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);

# F_W = 0 subsidy equilibrium
no_F_W_parameter = copy(Baseline_parameter);
no_F_W_parameter.F_W= 0.0;
parameters_tmp = copy(no_F_W_parameter);
balanced_share = 1.0
#foreign_supply_capital_subsidy_b = 7.628070905421486;
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end
prices = [ 1.4796471371798785,
0.14138775226957234,
0.16052856909641777]
# prices = 
#prices_bug = [  1.0012488427964574,1.304207153437122,0.06502460876457435,1.0524640422259952];
#prices = copy(prices_bug)
function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.05,0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =30000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_F_W_subsidy_b = copy(prices);

(residual_no_F_W, stat_distr_no_F_W, cons_fine_local_no_F_W, future_occupation_fine_local_no_F_W,x_SC_fine_no_F_W,x_BC_fine_no_F_W, 
    coeff_no_F_W,residual_goods_no_F_W,model_moments_no_F_W,foreign_supply_capital_no_F_W)  = solve_model_calibration2(prices_no_F_W_subsidy_b,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);

# F_W = 0 no subsidy equilibrium
no_F_W_no_subsidy_parameter = copy(Baseline_parameter);
no_F_W_no_subsidy_parameter.F_W= 0.0;
no_F_W_no_subsidy_parameter.τ_S=-0.0;
parameters_tmp = copy(no_F_W_no_subsidy_parameter);
balanced_share = 0.0;

#foreign_supply_capital_subsidy_b = 7.628070905421486;
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end
prices = [   1.1678646848278007,
0.0909235713934842];#!
# prices = 
#prices_bug = [  1.0012488427964574,1.304207153437122,0.06502460876457435,1.0524640422259952];
#prices = copy(prices_bug)
function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.05,0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =30000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_F_W_no_subsidy = copy(prices);

(residual_no_F_W_no_subsidy, stat_distr_no_F_W_no_subsidy, cons_fine_local_no_F_W_no_subsidy, future_occupation_fine_local_no_F_W_no_subsidy,x_SC_fine_no_F_W_no_subsidy,x_BC_fine_no_F_W_no_subsidy, 
    coeff_no_F_W_no_subsidy,residual_goods_no_F_W_no_subsidy,model_moments_no_F_W_no_subsidy,foreign_supply_capital_no_F_W_no_subsidy)  = solve_model_calibration2(prices_no_F_W_no_subsidy,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);

# FM_B = 0; balanced subsidy
no_FM_B_parameter = copy(Baseline_parameter);
no_FM_B_parameter.FM_B= 0.0;
parameters_tmp = copy(no_FM_B_parameter);
balanced_share = 1.0
#foreign_supply_capital_subsidy_b = 7.628070905421486;
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end
prices = [  1.4805614993444252,
0.2696222374380153,
0.1668785715032078]
# prices = 
#prices_bug = [  1.0012488427964574,1.304207153437122,0.06502460876457435,1.0524640422259952];
#prices = copy(prices_bug)
function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.05,0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =30000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_FM_B_subsidy_b = copy(prices);

(residual_no_FM_B_subsidy_b, stat_distr_no_FM_B_subsidy_b, cons_fine_local_no_FM_B_subsidy_b, future_occupation_fine_local_no_FM_B_subsidy_b,x_SC_fine_no_FM_B_subsidy_b,x_BC_fine_no_FM_B_subsidy_b, 
    coeff_no_FM_B_subsidy_b,residual_goods_no_FM_B_subsidy_b,model_moments_no_FM_B_subsidy_b,foreign_supply_capital_no_FM_B_subsidy_b)  = solve_model_calibration2(prices_no_FM_B_subsidy_b,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);

# FM_B = 0; no subsidy
no_FM_B_no_subsidy_parameter = copy(Baseline_parameter);
no_FM_B_no_subsidy_parameter.FM_B= 0.0;
no_FM_B_no_subsidy_parameter.τ_S=-0.0;
parameters_tmp = copy(no_FM_B_no_subsidy_parameter);
balanced_share = 0.0
#foreign_supply_capital_subsidy_b = 7.628070905421486;
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end
prices = [   1.1685831465774807,
0.21571073567745985];
# prices = 
#prices_bug = [  1.0012488427964574,1.304207153437122,0.06502460876457435,1.0524640422259952];
#prices = copy(prices_bug)
function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.05,0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =30000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_FM_B_no_subsidy = copy(prices);

(residual_no_FM_B_no_subsidy, stat_distr_no_FM_B_no_subsidy, cons_fine_local_no_FM_B_no_subsidy, future_occupation_fine_local_no_FM_B_no_subsidy,x_SC_fine_no_FM_B_no_subsidy,x_BC_fine_no_FM_B_no_subsidy, 
    coeff_no_FM_B_no_subsidy,residual_goods_no_FM_B_no_subsidy,model_moments_no_FM_B_no_subsidy,foreign_supply_capital_no_FM_B_no_subsidy)  = solve_model_calibration2(prices_no_FM_B_no_subsidy,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);


# No Substinence param + subsidy
no_cbar_parameter = copy(Baseline_parameter);
no_cbar_parameter.c̄_S= 0.0;
parameters_tmp = copy(no_cbar_parameter);
balanced_share = 1.0
#foreign_supply_capital_subsidy_b = 7.628070905421486;
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end
prices = [  1.5435433787380388,
0.259856531535514,
0.19293976562040993]
# prices = 
#prices_bug = [  1.0012488427964574,1.304207153437122,0.06502460876457435,1.0524640422259952];
#prices = copy(prices_bug)
function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.05,0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =30000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_cbar_subsidy_b = copy(prices);

(residual_no_cbar_subsidy_b, stat_distr_no_cbar_subsidy_b, cons_fine_local_no_cbar_subsidy_b, future_occupation_fine_local_no_cbar_subsidy_b,x_SC_fine_no_cbar_subsidy_b,x_BC_fine_no_cbar_subsidy_b, 
    coeff_no_cbar_subsidy_b,residual_goods_no_cbar_subsidy_b,model_moments_no_cbar_subsidy_b,foreign_supply_capital_no_cbar_subsidy_b)  = solve_model_calibration2(prices_no_cbar_subsidy_b,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);

# No Substinence param - no subsidy
no_cbar_no_subsidy_parameter = copy(Baseline_parameter);
no_cbar_no_subsidy_parameter.c̄_S= 0.0;
no_cbar_no_subsidy_parameter.τ_S = 0.0
parameters_tmp = copy(no_cbar_no_subsidy_parameter);
balanced_share = 0.0
#foreign_supply_capital_subsidy_b = 7.628070905421486;
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end

prices = [   1.205937339361828,
0.18357263199475019];
# prices = 
#prices_bug = [  1.0012488427964574,1.304207153437122,0.06502460876457435,1.0524640422259952];
#prices = copy(prices_bug)
function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.05,0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =30000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_cbar_no_subsidy = copy(prices);

(residual_no_cbar_no_subsidy, stat_distr_no_cbar_no_subsidy, cons_fine_local_no_cbar_no_subsidy, future_occupation_fine_local_no_cbar_no_subsidy,x_SC_fine_no_cbar_no_subsidy,x_BC_fine_no_cbar_no_subsidy, 
    coeff_no_cbar_no_subsidy,residual_goods_no_cbar_no_subsidy,model_moments_no_cbar_no_subsidy,foreign_supply_capital_no_cbar_no_subsidy)  = solve_model_calibration2(prices_no_cbar_no_subsidy,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);


# No Substinence + Qs param - no subsidy
no_cbarQS_no_subsidy_parameter = copy(Baseline_parameter);
no_cbarQS_no_subsidy_parameter.c̄_S= 0.0;
no_cbarQS_no_subsidy_parameter.Q_S = 0.0
no_cbarQS_no_subsidy_parameter.τ_S = 0.0
parameters_tmp = copy(no_cbarQS_no_subsidy_parameter);
balanced_share = 0.0
#foreign_supply_capital_subsidy_b = 7.628070905421486;
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end
prices = [  1.1985922807360259,
0.16529854238410874];
# prices = 
#prices_bug = [  1.0012488427964574,1.304207153437122,0.06502460876457435,1.0524640422259952];
#prices = copy(prices_bug)
function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.05,0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =30000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_cbarQS_no_subsidy = copy(prices);

(residual_no_cbarQS_no_subsidy, stat_distr_no_cbarQS_no_subsidy, cons_fine_local_no_cbarQS_no_subsidy, future_occupation_fine_local_no_cbarQS_no_subsidy,x_SC_fine_no_cbarQS_no_subsidy,x_BC_fine_no_cbarQS_no_subsidy, 
    coeff_no_cbarQS_no_subsidy,residual_goods_no_cbarQS_no_subsidy,model_moments_no_cbarQS_no_subsidy,foreign_supply_capital_no_cbarQS_no_subsidy)  = solve_model_calibration2(prices_no_cbarQS_no_subsidy,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);


# No Substinence + Qs param -  subsidy
no_cbarQS_parameter = copy(Baseline_parameter);
no_cbarQS_parameter.c̄_S= 0.0;
no_cbarQS_parameter.Q_S = 0.0
parameters_tmp = copy(no_cbarQS_parameter);
balanced_share = 1.0
#foreign_supply_capital_subsidy_b = 7.628070905421486;
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end
prices = [  1.5168290383060747,
0.20967638385502368,
0.168568086189524]
# prices = 
#prices_bug = [  1.0012488427964574,1.304207153437122,0.06502460876457435,1.0524640422259952];
#prices = copy(prices_bug)
function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.05,0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =30000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_cbarQS_subsidy_b = copy(prices);

(residual_no_cbarQS_subsidy_b, stat_distr_no_cbarQS_subsidy_b, cons_fine_local_no_cbarQS_subsidy_b, future_occupation_fine_local_no_cbarQS_subsidy_b,x_SC_fine_no_cbarQS_subsidy_b,x_BC_fine_no_cbarQS_subsidy_b, 
    coeff_no_cbarQS_subsidy_b,residual_goods_no_cbarQS_subsidy_b,model_moments_no_cbarQS_subsidy_b,foreign_supply_capital_no_cbarQS_subsidy_b)  = solve_model_calibration2(prices_no_cbarQS_subsidy_b,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);



# No Qs param - NO subsidy
no_QS_no_subsidy_parameter = copy(Baseline_parameter);
no_QS_no_subsidy_parameter.Q_S = 0.0
no_QS_no_subsidy_parameter.τ_S=-0;
parameters_tmp = copy(no_QS_no_subsidy_parameter);
balanced_share = 0.0
#foreign_supply_capital_subsidy_b = 7.628070905421486;
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end
prices = [  1.1860599621803165,
0.16412849229735588]
# prices = 
#prices_bug = [  1.0012488427964574,1.304207153437122,0.06502460876457435,1.0524640422259952];
#prices = copy(prices_bug)
function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.05,0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =30000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_QS_no_subsidy = copy(prices);

(residual_no_QS_no_subsidy, stat_distr_no_QS_no_subsidy, cons_fine_local_no_QS_no_subsidy, future_occupation_fine_local_no_QS_no_subsidy,x_SC_fine_no_QS_no_subsidy,x_BC_fine_no_QS_no_subsidy, 
    coeff_no_QS_no_subsidy,residual_goods_no_QS_no_subsidy,model_moments_no_QS_no_subsidy,foreign_supply_capital_no_QS_no_subsidy)  = solve_model_calibration2(prices_no_QS_no_subsidy,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);

    # No Qs param - subsidy 
no_QS_parameter = copy(Baseline_parameter);
no_QS_parameter.Q_S= 0.0;
parameters_tmp = copy(no_QS_parameter);
balanced_share = 1.0
#foreign_supply_capital_subsidy_b = 7.628070905421486;
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end
prices = [   1.5082999164210797,
0.20593172383716568,
0.17768187159642707]
# prices = 
#prices_bug = [  1.0012488427964574,1.304207153437122,0.06502460876457435,1.0524640422259952];
#prices = copy(prices_bug)
function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.05,0.01);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =30000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_QS_subsidy_b = copy(prices);

(residual_no_QS_subsidy_b, stat_distr_no_QS_subsidy_b, cons_fine_local_no_QS_subsidy_b, future_occupation_fine_local_no_QS_subsidy_b,x_SC_fine_no_QS_subsidy_b,x_BC_fine_no_QS_subsidy_b, 
    coeff_no_QS_subsidy_b,residual_goods_no_QS_subsidy_b,model_moments_no_QS_subsidy_b,foreign_supply_capital_no_QS_subsidy_b)  = solve_model_calibration2(prices_no_QS_subsidy_b,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);



# No Qs param, higher migration cost - NO subsidy
no_QS_recab_no_subsidy_parameter = copy(Baseline_parameter);
no_QS_recab_no_subsidy_parameter.Q_S = 0.0
no_QS_recab_no_subsidy_parameter.τ_S=-0;
no_QS_recab_no_subsidy_parameter.F_W = 305.0
parameters_tmp = copy(no_QS_recab_no_subsidy_parameter);
balanced_share = 0.0
#foreign_supply_capital_subsidy_b = 7.628070905421486;
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end
prices = [ 1.190752722671664,
0.16670260063915554]
# prices = 
#prices_bug = [  1.0012488427964574,1.304207153437122,0.06502460876457435,1.0524640422259952];
#prices = copy(prices_bug)
function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.05,0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =30000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_QS_recab_no_subsidy = copy(prices);

(residual_no_QS_recab_no_subsidy, stat_distr_no_QS_recab_no_subsidy, cons_fine_local_no_QS_recab_no_subsidy, future_occupation_fine_local_no_QS_recab_no_subsidy,x_SC_fine_no_QS_recab_no_subsidy,x_BC_fine_no_QS_recab_no_subsidy, 
    coeff_no_QS_recab_no_subsidy,residual_goods_no_QS_recab_no_subsidy,model_moments_no_QS_recab_no_subsidy,foreign_supply_capital_no_QS_recab_no_subsidy)  = solve_model_calibration2(prices_no_QS_recab_no_subsidy,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);

    # No Qs param higher migration cost - subsidy 
no_QS_recab_parameter = copy(Baseline_parameter);
no_QS_recab_parameter.Q_S= 0.0;
no_QS_recab_parameter.F_W = 305.0
parameters_tmp = copy(no_QS_recab_parameter);
balanced_share = 1.0
#foreign_supply_capital_subsidy_b = 7.628070905421486;
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end
prices = [1.5163159754906252,
0.24936579281247012,
0.1755991002351221]
# prices = 
#prices_bug = [  1.0012488427964574,1.304207153437122,0.06502460876457435,1.0524640422259952];
#prices = copy(prices_bug)
function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.05,0.01);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =30000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_no_QS_recab_subsidy_b = copy(prices);

(residual_no_QS_recab_subsidy_b, stat_distr_no_QS_recab_subsidy_b, cons_fine_local_no_QS_recab_subsidy_b, future_occupation_fine_local_no_QS_recab_subsidy_b,x_SC_fine_no_QS_recab_subsidy_b,x_BC_fine_no_QS_recab_subsidy_b, 
    coeff_no_QS_recab_subsidy_b,residual_goods_no_QS_recab_subsidy_b,model_moments_no_QS_recab_subsidy_b,foreign_supply_capital_no_QS_recab_subsidy_b)  = solve_model_calibration2(prices_no_QS_recab_subsidy_b,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);


#### Infrastructure

#### Infrastructure without subsidy and balance budget no spillovers
balanced_share = 0.0
infra_parameter_nsp_nb = copy(Baseline_parameter);
infra_parameter_nsp_nb.τ_S=-0;
infra_parameter_nsp_nb.Q_S = 0.7
prices = [ 1.190206255203472,
0.1745527606893067];#,2.505662787966534
#foreign_supply_capital_subsidy_b = 7.628070905421486; # through the calibration of subsidy_b
#transaction_cost_loss_subsidy_b = 0.29321009828641925;
#transaction_cost_loss_no_subsidy = 0.28506602793158564;
println("Solving the model ...")
#out=1/2 - choose 1 for solver, 2 for analysis
function f_solve(prices)
    residual = solve_model_calibration2(prices,infra_parameter_nsp_nb,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end

function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;

stst_simplex = Optim.AffineSimplexer(0.05,0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =30000.0);
prices = optim_res.minimizer;
prices_inf = copy(prices);
#residual = f_solve(prices);
(residual_goods_inf, stat_distr_inf, cons_fine_local_inf, a_prime_fine_local_inf,future_occupation_fine_local_inf,x_S_S_fine_inf,x_SC_fine_inf,x_BC_fine_inf, coeff_inf,
        transaction_cost_loss_inf,nominal_GDP_inf,welfare_val_inf,Import_value_inf,Export_value_inf,current_worker_pop_inf,current_staple_pop_inf,current_cashcrop_pop_inf,
        marketable_agr_surplus_share_inf,exportshare_cashcrop_inf,
        fraction_model_inf,program_spending_inf,prod_value_improvement_inf,share_selling_increase_inf,exp_ratio_model_inf,mig_rate_model_inf,rural_pop_only_staples_model_inf,rural_pop_only_cashcrop_model_inf,
        mean_land_share_to_staples_among_cc_model_inf,urban_rural_inc_ratio_model_inf,urban_rural_wealth_ratio_model_inf,urban_rural_consumption_ratio_model_inf,
        p90_wealth_rural_inf,p90_wealth_urban_inf,p99_wealth_rural_inf,p99_wealth_urban_inf,p90_cons_tmp_inf,p90_income_tmp_inf,p99_cons_tmp_inf,p99_income_tmp_inf,
        staple_productivity_inf,cashcrop_productivity_inf,manuf_productivity_inf,relative_land_to_staples_inf,relative_land_to_cashcrop_inf,share_constrained_cashcrop_inf,var_MPX_staples_S_inf,var_MPX_cashcrop_B_inf,var_MPX_cashcrop_S_inf,
        share_constrained_staple_inf,APG_inf,urban_rural_consumption_ratio_model_real_inf,aggregate_consumption_inf,
        worker_pop_effective_inf,prod_manuf_inf,total_entry_cost_inf,prod_staple_inf,prod_cashcrop_inf,input_staple_inf,input_cashcrop_inf,
        total_maintenance_cost_inf,current_account_residual_inf,fraction_cashcrop_suboptimal_model_inf,
        c_B_worker_sum_no_subsidy,c_B_staple_sum_no_subsidy,c_B_cashcrop_sum_no_subsidy,c_S_worker_sum_no_subsidy,c_S_staple_sum_no_subsidy ,c_S_cashcrop_sum_no_subsidy,
        transaction_cost_staple_sum_no_subsidy,transaction_cost_cashcrop_sum_no_subsidy,transaction_cost_worker_sum_no_subsidy,c_M_worker_sum_no_subsidy,c_M_staple_sum_no_subsidy,c_M_cashcrop_sum_no_subsidy,
        MPX_mean_log_no_subsidy, MPX_mean_staples_S_log_no_subsidy,MPX_mean_cashcrop_log_no_subsidy,
c_S_W_fine_no_subsidy,c_B_W_fine_no_subsidy,c_M_W_fine_no_subsidy,c_S_S_fine_no_subsidy,c_B_S_fine_no_subsidy,c_M_S_fine_no_subsidy,c_S_B_fine_no_subsidy,c_B_B_fine_no_subsidy,c_M_B_fine_no_subsidy) = details_model(prices_inf,infra_parameter_nsp_nb,2,moments,0.0,foreign_supply_capital_subsidy_b);

transaction_cost_loss_inf*prod_staple_inf/ (infra_parameter_nsp_nb.Q_S)*(Baseline_parameter.Q_S-infra_parameter_nsp_nb.Q_S) / nominal_GDP_inf
transaction_cost_loss_inf*prod_staple_inf/ (infra_parameter_nsp_nb.Q_S)*(Baseline_parameter.Q_S-infra_parameter_nsp_nb.Q_S)
program_spending_subsidy_nb * nominal_GDP_subsidy_nb

# Infrastructure without balanced budget with spillovers resulting in the same reduction in F_W
balanced_share = 0.0
parameters_tmp = copy(Baseline_parameter);
parameters_tmp.Q_S = 0.7
parameters_tmp.F_W = Baseline_parameter.F_W * (parameters_tmp.Q_S / Baseline_parameter.Q_S)
parameters_tmp.τ_S=-0.0;
prices = [ 1.1723046609628551
0.12471486958333933];#,2.505662787966534
#foreign_supply_capital_subsidy_b = 7.628070905421486; # through the calibration of subsidy_b
#transaction_cost_loss_subsidy_b = 0.29321009828641925;
#transaction_cost_loss_no_subsidy = 0.28506602793158564;
println("Solving the model ...")
#out=1/2 - choose 1 for solver, 2 for analysis
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end

function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;

stst_simplex = Optim.AffineSimplexer(0.05,0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =30000.0);
prices = optim_res.minimizer;
prices_inf_sp = copy(prices);
#residual = f_solve(prices);
(residual_goods_inf_sp, stat_distr_inf_sp, cons_fine_local_inf_sp, a_prime_fine_local_inf_sp,future_occupation_fine_local_inf_sp,x_S_S_fine_inf_sp,x_SC_fine_inf_sp,x_BC_fine_inf_sp, coeff_inf_sp,
        transaction_cost_loss_inf_sp,nominal_GDP_inf_sp,welfare_val_inf_sp,Import_value_inf_sp,Export_value_inf_sp,current_worker_pop_inf_sp,current_staple_pop_inf_sp,current_cashcrop_pop_inf_sp,
        marketable_agr_surplus_share_inf_sp,exportshare_cashcrop_inf_sp,
        fraction_model_inf_sp,program_spending_inf_sp,prod_value_improvement_inf_sp,share_selling_increase_inf_sp,exp_ratio_model_inf_sp,mig_rate_model_inf_sp,rural_pop_only_staples_model_inf_sp,rural_pop_only_cashcrop_model_inf_sp,
        mean_land_share_to_staples_among_cc_model_inf_sp,urban_rural_inc_ratio_model_inf_sp,urban_rural_wealth_ratio_model_inf_sp,urban_rural_consumption_ratio_model_inf_sp,
        p90_wealth_rural_inf_sp,p90_wealth_urban_inf_sp,p99_wealth_rural_inf_sp,p99_wealth_urban_inf_sp,p90_cons_tmp_inf_sp,p90_income_tmp_inf_sp,p99_cons_tmp_inf_sp,p99_income_tmp_inf_sp,
        staple_productivity_inf_sp,cashcrop_productivity_inf_sp,manuf_productivity_inf_sp,relative_land_to_staples_inf_sp,relative_land_to_cashcrop_inf_sp,share_constrained_cashcrop_inf_sp,var_MPX_staples_S_inf_sp,var_MPX_cashcrop_B_inf_sp,var_MPX_cashcrop_S_inf_sp,
        share_constrained_staple_inf_sp,APG_inf_sp,urban_rural_consumption_ratio_model_real_inf_sp,aggregate_consumption_inf_sp,
        worker_pop_effective_inf_sp,prod_manuf_inf_sp,total_entry_cost_inf_sp,prod_staple_inf_sp,prod_cashcrop_inf_sp,input_staple_inf_sp,input_cashcrop_inf_sp,
        total_maintenance_cost_inf_sp,current_account_residual_inf_sp,fraction_cashcrop_suboptimal_model_inf_sp,
        c_B_worker_sum_no_subsidy,c_B_staple_sum_no_subsidy,c_B_cashcrop_sum_no_subsidy,c_S_worker_sum_no_subsidy,c_S_staple_sum_no_subsidy ,c_S_cashcrop_sum_no_subsidy,
        transaction_cost_staple_sum_no_subsidy,transaction_cost_cashcrop_sum_no_subsidy,transaction_cost_worker_sum_no_subsidy,c_M_worker_sum_no_subsidy,c_M_staple_sum_no_subsidy,c_M_cashcrop_sum_no_subsidy
        ,
        MPX_mean_log_no_subsidy, MPX_mean_staples_S_log_no_subsidy,MPX_mean_cashcrop_log_no_subsidy)  = details_model(prices_inf_sp,parameters_tmp,2,moments,0.0,foreign_supply_capital_subsidy_b);

#program_spending_subsidy_nb * nominal_GDP_subsidy_nb - transaction_cost_loss_inf_sp/ (parameters_tmp.Q_S)*(Baseline_parameter.Q_S-parameters_tmp.Q_S)
transaction_cost_loss_inf_sp*prod_staple_inf_sp* (parameters_tmp.Q_S)*(Baseline_parameter.Q_S-parameters_tmp.Q_S) / nominal_GDP_inf_sp
# Decrease Q_S until this is zero 
# Optimal subsidy graph

# Subsidy with balanced budget with different subsidies
no_taus= 21 # (the final tauS is always zero, and is solved elsewhere. Starting price guess is the previous price)
#foreign_supply_capital_subsidy_b = 7.628070905421486;
prices = [ 1.518014329164464,
0.2726473607850526,
0.1886863641824446];
prices_subsidy_b_grid =  ones(3,no_taus);
residual_goods_subsidy_b_grid = ones(6,no_taus);
iterate = 1;
for τ = range(Baseline_parameter.τ_S,-0.05,length = no_taus-1)
    balanced_share = 1.0
    parameters_tmp = copy(Baseline_parameter);
    parameters_tmp.τ_S = τ;
    println("Solving the model with " ,τ )
    #out=1/2 - choose 1 for solver, 2 for analysis
    function f_solve(prices)
        residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
        return residual
    end

    function f_sum(prices)
        residual = f_solve(prices);
        return sum(residual.^2)
    end
    lower_guess_bound = 0.1;
    upper_guess_bound = 5.0;
    ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
    x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
    upper = upper_guess_bound * prices);
    prices = ls_res.minimizer;

    stst_simplex = Optim.AffineSimplexer(0.05,0.1);
    optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
    extended_trace=true,time_limit =30000.0);
    prices = optim_res.minimizer;
    prices_subsidy_b = copy(prices);
    #residual = f_solve(prices);
    (residual_subsidy_b, stat_distr_subsidy_b, cons_fine_local_subsidy_b, future_occupation_fine_local_subsidy_b,x_SC_fine_subsidy_b,x_BC_fine_subsidy_b, 
        coeff_subsidy_b,residual_goods_subsidy_b,model_moments_subsidy_b,foreign_supply_capital_subsidy_b)  = solve_model_calibration1(prices_subsidy_b,parameters_tmp,2,moments,balanced_share);
    prices_subsidy_b_grid[:,iterate] = prices_subsidy_b
    residual_goods_subsidy_b_grid[:,iterate] = residual_goods_subsidy_b
    iterate = iterate + 1;
end
prices_subsidy_b_grid[1,no_taus] =  1.1882831767847994;
prices_subsidy_b_grid[2,no_taus] =   0.18349666056314162;
prices_subsidy_b_grid[3,no_taus] = 0.00;

# Improve on accuracy:
prices_subsidy_b_grid1 = [1.5182148598312672 1.4748615289082854 1.4389313172837825 1.410732954524286 1.390264709095618 1.3688114420290358 1.3436655727518898 1.3271570823425125 1.312778814101722 1.3003105888992221 1.284946210796012 1.2714091997394257 1.2599461971286845 1.2479933342455625 1.2388104401541415 1.2294244830748497 1.2217938803468038 1.2132098804982514 1.2099699888601618 1.1995655457878038 1.1882831767847994; 0.27262899815531444 0.26338424759605666 0.245822286882918 0.2442836837960706 0.21752367635207415 0.21670507564964162 0.2061107155376924 0.2036475392226677 0.20098618728681622 0.19862949627869536 0.19659831128032512 0.19476856375812163 0.19323999476869438 0.18933334499100343 0.19016660308519126 0.18807703128020226 0.1869710641961219 0.1857564570225632 0.18361768757377586 0.18414175648950032 0.18349666056314162; 0.18860594848809528 0.14681832007063747 0.12372332717865375 0.10211507272581569 0.09328546623799094 0.0770503430189069 0.06580378586580023 0.055857350780308554 0.04822701987252334 0.04189534576885402 0.03530977698358126 0.029721337552662572 0.0251658618411578 0.02075905911514478 0.016733799327583475 0.013639766449854727 0.010610113658055436 0.007772421378081955 0.005375731191350952 0.002777150926899589 0.0]
no_taus = 21;
τ_grid = range(Baseline_parameter.τ_S,-0.05,length = no_taus-1)
τ_grid1 = zeros(no_taus);
τ_grid1[1:(no_taus-1)] = τ_grid;
τ_grid = τ_grid1; # [-0.8099992851842415, -0.7726309017534919, -0.7352625183227424, -0.7255549201637702, -0.6978941348919928, -0.6605257514612433, -0.641110555143299, -0.6231573680304937, -0.5857889845997442, -0.5566661901228277, -0.5484206011689946, -0.511052217738245, -0.4736838343074955, -0.4722218251023564, -0.43631545087674595, -0.3989470674459964, -0.3877774600818851, -0.36157868401524684, -0.3242103005844973, -0.3033330950614138, -0.28684191715374774, -0.24947353372299821, -0.21888873004094256, -0.21210515029224866, -0.1747367668614991, -0.13736838343074956, -0.13444436502047127, -0.1, -0.05, 0.0]
foreign_supply_capital_subsidy_b = 7.631597383998027;
#residual_goods_subsidy_b_grid1= [-0.01956758283830878 -0.014970121497167551 -0.009634245958413988 -0.0050474924488701466 -0.008340837245202898 -0.006614690793662608 -0.008487685428446039 -3.963924061358579e-5 0.0017926087891720365 0.0027034153160704023 0.0022428494517553197 0.008498053424513092 0.0081455206302953 0.00477403738489211 0.008301229385871675 0.008074964133256075 0.008386173071680565 0.00878582463657438 0.009421960603803842 0.00520385127230347 0.011049880301909683 0.00949704140208864 -0.0016714837879940186 0.014509399539817272 0.007478009027943847 0.00849512064894457 0.007819783863096506 0.014947222211959432 0.015034256729904466 0.018364138555352894; -0.01671542556265158 -0.005092969879332773 -0.007182782153106135 -0.010429517322365876 -0.0034307651895043203 -0.0006085789970615635 0.0051576962369826165 0.006698569967892861 0.008766386935774995 0.007733037608850627 0.007911565264676814 0.00722800100242349 0.010308693988353747 0.008277626086222881 0.011741228993160486 0.013012726906530794 0.013131744360423516 0.014023815071873625 0.012359107063162188 0.018605276330044305 0.00989717442438452 0.012834066307019698 0.03489721739681127 0.015336218584723268 0.025990131993149063 0.02614573502016043 0.01789035892615273 0.0163332071680088 0.0184235330695609 0.014433668788282351; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; -0.0045164134782469945 -0.0032023125932608308 -0.0027020270845698117 -0.002036984644452973 -0.002212749801021731 -0.0019473040708817982 -0.0013385317320359785 0.0003190407204560143 0.0005637787919106929 0.0006543311567283322 0.0016835673750581013 0.001415896766939693 0.0021469721584954373 0.0005965438751970963 0.0021713594736638255 0.002412327838978579 0.0024810711427011634 0.00273763493443695 0.002507943392272483 0.0029755772946017515 0.0025242854831828507 0.003153057205792137 -0.00523415346013772 0.004243758019535307 0.0064507632295426965 0.0040125221775699094 0.005149198309804932 0.004865490098072679 0.004047477992794155 0.0056533092135844735; -0.001348120684910654 -0.00037321519163053494 -0.0004905637067318403 -0.0005996043792757603 -0.0003800913955564468 -0.0002622748940468766 8.960096479860029e-5 -0.00020739100956583792 3.23748839643273e-5 7.866789238453593e-5 0.002452617586017561 3.316695121996137e-5 4.8890586173804226e-5 2.585182787609997e-5 0.0001414965959856684 0.00011547355269980133 0.00015391070709058635 4.0854922485247e-5 2.3646204003389763e-5 8.319989757816823e-5 -3.092128624674716e-5 2.047559165204729e-5 -0.02976379804887917 0.00011404851677601807 0.0001248104149185 -0.00020123175633471427 4.4480776693604385e-5 -4.442345074475462e-6 0.00010453584154952253 NaN] ;
prices_subsidy_b_grid = copy(prices_subsidy_b_grid1);
#residual_goods_subsidy_b_grid = copy(residual_goods_subsidy_b_grid1);
residual_goods_subsidy_b_grid = ones(6,no_taus);
iterate = 1;
while iterate < no_taus
    balanced_share = 1.0
    parameters_tmp = copy(Baseline_parameter); 
    parameters_tmp.τ_S = τ_grid[iterate];
    println("Solving the model with " ,τ_grid[iterate], "which is the " , iterate ," iteration" )
    #out=1/2 - choose 1 for solver, 2 for analysis
    function f_solve(prices)
        residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
        return residual
    end

    function f_sum(prices)
        residual = f_solve(prices);
        return sum(residual.^2)
    end
    prices = prices_subsidy_b_grid[:,iterate];
    lower_guess_bound = 0.1;
    upper_guess_bound = 5.0;

    ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
    x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
    upper = upper_guess_bound * prices);
    prices = ls_res.minimizer;

    stst_simplex = Optim.AffineSimplexer(0.05,0.1);
    optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
    extended_trace=true,time_limit =30000.0);
    prices = optim_res.minimizer;
    prices_subsidy_b = copy(prices);
    #residual = f_solve(prices);
    (residual_subsidy_b, stat_distr_subsidy_b, cons_fine_local_subsidy_b, future_occupation_fine_local_subsidy_b,x_SC_fine_subsidy_b,x_BC_fine_subsidy_b, 
        coeff_subsidy_b,residual_goods_subsidy_b,model_moments_subsidy_b,foreign_supply_capital_subsidy_b)  = solve_model_calibration1(prices_subsidy_b,parameters_tmp,2,moments,balanced_share);
    prices_subsidy_b_grid[:,iterate] = prices_subsidy_b
    residual_goods_subsidy_b_grid[:,iterate] = residual_goods_subsidy_b
    iterate = iterate + 1;
end

# Optimal subsidy graph for Q_S = 0

# Subsidy with balanced budget with different subsidies
no_taus= 21 # (the final tauS is always zero, and is solved elsewhere. Starting price guess is the previous price)
#foreign_supply_capital_subsidy_b = 7.628070905421486;
prices = [ 1.5082999164210797,
0.20593172383716568,
0.17768187159642707];
prices_subsidy_b_grid =  ones(3,no_taus);
residual_goods_subsidy_b_grid = ones(6,no_taus);
iterate = 1;
for τ = range(Baseline_parameter.τ_S,-0.05,length = no_taus-1)
    balanced_share = 1.0
    parameters_tmp = copy(Baseline_parameter);
    parameters_tmp.τ_S = τ;
    parameters_tmp.Q_S = 0.0;
    println("Solving the model with " ,τ )
    #out=1/2 - choose 1 for solver, 2 for analysis
    function f_solve(prices)
        residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
        return residual
    end

    function f_sum(prices)
        residual = f_solve(prices);
        return sum(residual.^2)
    end
    lower_guess_bound = 0.1;
    upper_guess_bound = 5.0;
    ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
    x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
    upper = upper_guess_bound * prices);
    prices = ls_res.minimizer;

    stst_simplex = Optim.AffineSimplexer(0.05,0.1);
    optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
    extended_trace=true,time_limit =30000.0);
    prices = optim_res.minimizer;
    prices_subsidy_b = copy(prices);
    #residual = f_solve(prices);
    (residual_subsidy_b, stat_distr_subsidy_b, cons_fine_local_subsidy_b, future_occupation_fine_local_subsidy_b,x_SC_fine_subsidy_b,x_BC_fine_subsidy_b, 
        coeff_subsidy_b,residual_goods_subsidy_b,model_moments_subsidy_b,foreign_supply_capital_subsidy_b)  = solve_model_calibration2(prices_subsidy_b,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);
    prices_subsidy_b_grid[:,iterate] = prices_subsidy_b
    residual_goods_subsidy_b_grid[:,iterate] = residual_goods_subsidy_b
    iterate = iterate + 1;
end
prices_subsidy_b_grid[1,no_taus] =  1.1860599621803165;
prices_subsidy_b_grid[2,no_taus] =   0.16412849229735588;
prices_subsidy_b_grid[3,no_taus] = 0.00;

#prices_subsidy_b_grid= [1.4598772874688932 1.4385245014811625 1.4049741835258556 1.4077138037158219 1.3821166752679364 1.3651971203006512 1.3573364085808683 1.363686641369129 1.3414500347358624 1.3323226063843923 1.327875051226767 1.3206902627537573 1.31038149959868 1.299209888144449 1.295273367696625 1.286440471180219 1.280855187630534 1.2693788918619269 1.2505170635332994 1.3058767212231221 1.2528814375260455 1.2390136405180543 1.2366561653512491 1.2397711208427493 1.2369018943039145 1.2292651896985871 1.2234888879703227 1.223637192043632 1.2123437075134387 1.1761; 0.2519999236690381 0.23498396688130332 0.23332979499142487 0.2331762054176376 0.23381710497342015 0.23254851966653814 0.22970890394934568 0.21321983055377886 0.21472526959909455 0.21520463011066363 0.21512099960052214 0.20424937174792618 0.2030998391253458 0.21684048568198006 0.20338551740011912 0.2031107078028406 0.20457532622931957 0.20671971878014053 0.2078049839228001 0.20138832443176954 0.2059963050532268 0.20902560956191182 0.20847460704850557 0.19109415674986044 0.1903407252421793 0.19041346055148342 0.20864525308669363 0.18988145820182262 0.1898044796482799 0.192054; 0.2049632173803936 0.1411151956716539 0.1143093243650973 0.1027914440951735 0.09459381800661454 0.07766669772376758 0.0722107637510424 0.06298170571272088 0.0547728674793511 0.047763531492882316 0.04623130159584552 0.03929669248628405 0.033648141106931 0.03384699432359155 0.027747870442189175 0.023731464417786858 0.02272952750763354 0.0207361287258126 0.018317258948674173 0.01628515625104149 0.014790446150075959 0.01215828872230139 0.010102326008558985 0.00967354325689979 0.00756756417746829 0.005680704866949134 0.005586786834627149 0.003926480520014699 0.0018768700814345938 0.0];


# Higher prices of fertilizer

# Higher prices of fertilizer -  NO subsidy
high_p_x_no_subsidy_parameter = copy(Baseline_parameter);
high_p_x_no_subsidy_parameter.p_x = 2*Baseline_parameter.p_x;
high_p_x_no_subsidy_parameter.τ_S= -0.0;
#high_p_x_no_subsidy_parameter.c̄_S= 0.01975;
parameters_tmp = copy(high_p_x_no_subsidy_parameter);
balanced_share = 0.0
#foreign_supply_capital_subsidy_b = 7.631597488627073;
foreign_supply_capital_high_p_x_subsidy_b = 8.053580609954484;
function f_solve(prices)
    #residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    #residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_high_p_x_subsidy_b);
    residual = solve_model_calibration1(prices,parameters_tmp,1,moments,balanced_share);
    return residual
end
prices = [ 1.1936833649221996,
0.17237129935387707]
# prices = 
#prices_bug = [  1.0012488427964574,1.304207153437122,0.06502460876457435,1.0524640422259952];
#prices = copy(prices_bug)
function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.05,0.1);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =30000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_high_p_x_no_subsidy = copy(prices);

#(residual_high_p_x_no_subsidy, stat_distr_high_p_x_no_subsidy, cons_fine_local_high_p_x_no_subsidy, future_occupation_fine_local_high_p_x_no_subsidy,x_SC_fine_high_p_x_no_subsidy,x_BC_fine_high_p_x_no_subsidy, 
#    coeff_high_p_x_no_subsidy,residual_goods_high_p_x_no_subsidy,model_moments_high_p_x_no_subsidy,foreign_supply_capital_high_p_x_no_subsidy)  = solve_model_calibration2(prices_high_p_x_no_subsidy,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);

#(residual_high_p_x_no_subsidy, stat_distr_high_p_x_no_subsidy, cons_fine_local_high_p_x_no_subsidy, future_occupation_fine_local_high_p_x_no_subsidy,x_SC_fine_high_p_x_no_subsidy,x_BC_fine_high_p_x_no_subsidy, 
#    coeff_high_p_x_no_subsidy,residual_goods_high_p_x_no_subsidy,model_moments_high_p_x_no_subsidy,foreign_supply_capital_high_p_x_no_subsidy)  = solve_model_calibration2(prices_high_p_x_no_subsidy,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_high_p_x_subsidy_b);

(residual_high_p_x_no_subsidy, stat_distr_high_p_x_no_subsidy, cons_fine_local_high_p_x_no_subsidy, future_occupation_fine_local_high_p_x_no_subsidy,x_SC_fine_high_p_x_no_subsidy,x_BC_fine_high_p_x_no_subsidy, 
    coeff_high_p_x_no_subsidy,residual_goods_high_p_x_no_subsidy,model_moments_high_p_x_no_subsidy,foreign_supply_capital_high_p_x_no_subsidy)  = solve_model_calibration1(prices_high_p_x_no_subsidy,parameters_tmp,2,moments,balanced_share);

    #Higher prices of fertilizer  Subsidy without balanced budget
balanced_share = 0.0
high_p_x_no_subsidy_parameter = copy(Baseline_parameter);
high_p_x_no_subsidy_parameter.p_x = 2*Baseline_parameter.p_x;
high_p_x_no_subsidy_parameter.τ_S= -0.0;
parameters_tmp = copy(high_p_x_no_subsidy_parameter);
prices = [ 1.5087253977851676,
0.25456229359350663];
#foreign_supply_capital_subsidy_b = 7.628070905421486; # through the calibration of subsidy_b
println("Solving the model ...")
#out=1/2 - choose 1 for solver, 2 for analysis
function f_solve(prices)
    residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end

function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;

stst_simplex = Optim.AffineSimplexer(0.025,0.05);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =10000.0);
prices = optim_res.minimizer;
prices_high_p_x_subsidy_nb = copy(prices);
#residual = f_solve(prices);
(residual_high_p_x_subsidy_nb, stat_distr_high_p_x_subsidy_nb, cons_fine_local_high_p_x_subsidy_nb, future_occupation_fine_local_high_p_x_subsidy_nb,x_SC_fine_high_p_x_subsidy_nb,x_BC_fine_high_p_x_subsidy_nb, 
    coeff_high_p_x_subsidy_nb,residual_goods_high_p_x_subsidy_nb,model_moments_high_p_x_subsidy_nb,foreign_supply_capital_high_p_x_subsidy_nb)  = solve_model_calibration2(prices_high_p_x_subsidy_nb,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);


# Higher prices of fertilizer - subsidy 
high_p_x_parameter = copy(Baseline_parameter);
high_p_x_parameter.p_x = 2*Baseline_parameter.p_x;
parameters_tmp = copy(high_p_x_parameter);
balanced_share = 1.0
#foreign_supply_capital_subsidy_b = 7.631597488627073;
function f_solve(prices)
    residual = solve_model_calibration1(prices,parameters_tmp,1,moments,balanced_share);
    #residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
    return residual
end
#prices = [1.5374346868502475,
#0.22598290665699955,
#0.2352290656580293]

prices = [1.5379749200856943
0.22585638848772757
0.23665783406135696]
# prices = 
#prices_bug = [  1.0012488427964574,1.304207153437122,0.06502460876457435,1.0524640422259952];
#prices = copy(prices_bug)
function f_sum(prices)
    residual = f_solve(prices);
    return sum(residual.^2)
end
lower_guess_bound = 0.1;
upper_guess_bound = 5.0;
ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.05,0.01);
optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =30000.0);
prices = optim_res.minimizer;
residual = f_solve(prices);
prices_high_p_x_subsidy_b = copy(prices);

#(residual_high_p_x_subsidy_b, stat_distr_high_p_x_subsidy_b, cons_fine_local_high_p_x_subsidy_b, future_occupation_fine_local_high_p_x_subsidy_b,x_SC_fine_high_p_x_subsidy_b,x_BC_fine_high_p_x_subsidy_b, 
#    coeff_high_p_x_subsidy_b,residual_goods_high_p_x_subsidy_b,model_moments_high_p_x_subsidy_b,foreign_supply_capital_high_p_x_subsidy_b)  = solve_model_calibration2(prices_high_p_x_subsidy_b,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);

(residual_high_p_x_subsidy_b, stat_distr_high_p_x_subsidy_b, cons_fine_local_high_p_x_subsidy_b, future_occupation_fine_local_high_p_x_subsidy_b,x_SC_fine_high_p_x_subsidy_b,x_BC_fine_high_p_x_subsidy_b, 
    coeff_high_p_x_subsidy_b,residual_goods_high_p_x_subsidy_b,model_moments_high_p_x_subsidy_b,foreign_supply_capital_high_p_x_subsidy_b)  = solve_model_calibration1(prices_high_p_x_subsidy_b,parameters_tmp,2,moments,balanced_share);


# Optimal subsidy graph without balanced budget

# Subsidy without balanced budget
no_taus= 21 # (the final tauS is always zero, and is solved elsewhere. Starting price guess is the previous price)
foreign_supply_capital_subsidy_b = 7.631597488627073;
prices = [ 1.5087253977851676,
0.25456229359350663];
prices_subsidy_nb_grid =  ones(2,no_taus);
residual_goods_subsidy_nb_grid = ones(5,no_taus);
iterate = 1;
for τ = range(Baseline_parameter.τ_S,-0.05,length = no_taus-1)
    balanced_share = 0.0
    parameters_tmp = copy(Baseline_parameter);
    parameters_tmp.τ_S = τ;
    println("Solving the model with " ,τ )
    #out=1/2 - choose 1 for solver, 2 for analysis
    function f_solve(prices)
        residual = solve_model_calibration2(prices,parameters_tmp,1,moments,balanced_share,foreign_supply_capital_subsidy_b);
        return residual
    end

    function f_sum(prices)
        residual = f_solve(prices);
        return sum(residual.^2)
    end
    lower_guess_bound = 0.1;
    upper_guess_bound = 5.0;
    ls_res= LeastSquaresOptim.optimize(f_solve, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
    x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
    upper = upper_guess_bound * prices);
    prices = ls_res.minimizer;

    stst_simplex = Optim.AffineSimplexer(0.05,0.1);
    optim_res = LeastSquaresOptim.optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
    extended_trace=true,time_limit =30000.0);
    prices = optim_res.minimizer;
    prices_subsidy_b = copy(prices);
    #residual = f_solve(prices);
    (residual_subsidy_b, stat_distr_subsidy_b, cons_fine_local_subsidy_b, future_occupation_fine_local_subsidy_b,x_SC_fine_subsidy_b,x_BC_fine_subsidy_b, 
        coeff_subsidy_b,residual_goods_subsidy_b,model_moments_subsidy_b,foreign_supply_capital_subsidy_b)  = solve_model_calibration2(prices_subsidy_b,parameters_tmp,2,moments,balanced_share,foreign_supply_capital_subsidy_b);
    prices_subsidy_nb_grid[:,iterate] = prices_subsidy_b
    residual_goods_subsidy_nb_grid[:,iterate] = residual_goods_subsidy_b
    iterate = iterate + 1;
end
prices_subsidy_nb_grid[1,no_taus] =  1.1882831767847994;
prices_subsidy_nb_grid[2,no_taus] =   0.18349666056314162;