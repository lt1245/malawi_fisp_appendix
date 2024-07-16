#include("setup.jl")
#include("transition_functions.jl") to be included later - now all functions are here
# Transition paths
case_final =6; # Put to the main script
τ_index = 21 # for the optimal subsidy rate exercise , if case =! 5, no effect
epsilon_index = 6 # for the epsilon exercise, if case =! 6, no effect
load_solution = 1 ; # 1 if load saved solutions,
# 0 from scratch, -1 for case1 as initial guess, -2 is for the previous tau_rate/epsilon, -3 is for the next tau_rate/epsilon
generate_results = 1; # To use the loaded results for generating plots and results. 
#If 0, it starts searching for solution
save_solution = 1; # 1 to save solution at the end, only matters if generate_results==0
if case_final == 1
    #Effect of introducing the subsidy program with govt budget - no_subs to subsidy_b
    price_start = copy(prices_no_subsidy);
    coeff_begin = copy(coeff_no_subsidy);
    init_distr = copy(stat_distr_no_subsidy);
    init_capital = copy(foreign_supply_capital_subsidy_b);
    price_finish = copy(prices_subsidy_b);
    coeff_end = copy(coeff_subsidy_b);
    fini_distr = copy(stat_distr_subsidy_b );
    fini_capital = copy(foreign_supply_capital_subsidy_b);
    T = 50; #Initial guess for the transition length
    T_end = 15; # Initial guess - convergence after
    parameter_initial = copy(No_subsidy_parameter);
    parameter_end = copy(Baseline_parameter );
    ### Initialization based on the initial and final guess
    # Slow introduction ofthe subsidy rate
    T_tau_transition = 5; # 
    τ_trans = ones(T+2);
    τ_trans = parameter_end.τ_S * τ_trans[:];
    τ_trans[1:T_tau_transition] = range(parameter_initial.τ_S,parameter_end.τ_S,length = T_tau_transition);

    coeff_store = zeros(size(coeff_end,1),size(coeff_end,2),T+2);
    distr_store = zeros(size(init_distr,1),T+2);
    distr_store[:,1] = init_distr;
    distr_store[:,2] = init_distr;
    distr_store[:,T + 2] = fini_distr;

    coeff_store[:,:,1] = coeff_begin;
    coeff_store[:,:,T+2] = coeff_end;

    capital_trans = zeros(T+2);
    capital_trans[1] = init_capital;
    capital_trans[2] = init_capital;
    capital_trans[T+2] = fini_capital;
    capital_trans[3:(T+1)] .= fini_capital;
    V_saved_store = zeros(parameter_initial.ns_fine*3,T+2);
    V_saved_store[:,1] = V_saved_no_subsidy_reshaped;
    V_saved_store[:,T+2] = V_saved_subsidy_b_reshaped; 

    if load_solution == 1
        dataframe_prices = CSV.read("vec_price_intro_subsidy.csv", DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1;
    elseif load_solution == -1 
        # longer guesses
        T_guess = 30
        price_trans = zeros(3,T+2);
        price_trans[1:3,end] = price_finish;
        price_trans[1:2,1] = price_start;
        price_trans[:,2:(T_guess+1)] = reshape(vec_price_trans,3,T_guess);
        price_trans[1,(T_guess+2):end - 1] .= price_finish[1];
        price_trans[2,(T_guess+2):end - 1] .= price_finish[2]; 
        price_trans[3,(T_guess+2):end - 1] .= price_finish[3];
        price_trans_actual = price_trans[:,2:T+1]
        vec_price_trans = copy(reshape(price_trans_actual,(T)*3)); 
    else
        price_trans = zeros(3,T+2);
        price_trans[1:3,end] = price_finish;
        price_trans[1:2,1] = price_start;
        price_trans[1:2,2] = price_start;
        fraction_guess = range(parameter_end.τ_S,parameter_initial.τ_S,length = T_end)/parameter_end.τ_S
        fraction_guess = repeat(fraction_guess,1,3)'
        fraction_guess_tau = range(price_finish[3]-0.01,0.0,length = T_tau_transition)/price_finish[3]
        fraction_guess[3,1:T_tau_transition] =fraction_guess_tau
        fraction_guess[3,(T_tau_transition+1):end] .= 0.0;
        price_trans[1:3,2:(T_end+1)] = fraction_guess.*repeat(price_trans[:,2],1,T_end)  + (1 .- fraction_guess).*repeat(price_trans[:,end],1,T_end)
        price_trans[1,(T_end+2):end - 1] .= price_finish[1];
        price_trans[2,(T_end+2):end - 1] .= price_finish[2];
        price_trans[3,(T_tau_transition+2):end - 1] .= price_finish[3];
        price_trans_actual = price_trans[:,2:T+1]
        vec_price_trans = copy(reshape(price_trans_actual,(T)*3));
    end
elseif case_final == 2
    #Effect of introducing the subsidy program with aid finance
    price_start = copy(prices_no_subsidy);
    coeff_begin = copy(coeff_no_subsidy);
    init_distr = copy(stat_distr_no_subsidy);
    init_capital = copy(foreign_supply_capital_subsidy_b);
    price_finish = copy(prices_subsidy_nb);
    coeff_end = copy(coeff_subsidy_nb);
    fini_distr = copy(stat_distr_subsidy_nb );
    fini_capital = copy(foreign_supply_capital_subsidy_b);
    T = 50; #Initial guess for the transition length
    T_end = 15; # Initial guess - convergence after
    parameter_initial = copy(No_subsidy_parameter);
    parameter_end = copy(Baseline_parameter );
    ### Initialization based on the initial and final guess
    # Slow shutting down of subsidy rate
    T_tau_transition = 5; # 
    τ_trans = ones(T+2);
    τ_trans = parameter_end.τ_S * τ_trans[:];
    τ_trans[1:T_tau_transition] = range(parameter_initial.τ_S,parameter_end.τ_S,length = T_tau_transition);

    coeff_store = zeros(size(coeff_end,1),size(coeff_end,2),T+2);
    distr_store = zeros(size(init_distr,1),T+2);
    distr_store[:,1] = init_distr;
    distr_store[:,2] = init_distr;
    distr_store[:,T + 2] = fini_distr;

    coeff_store[:,:,1] = coeff_begin;
    coeff_store[:,:,T+2] = coeff_end;

    capital_trans = zeros(T+2);
    capital_trans[1] = init_capital;
    capital_trans[2] = init_capital;
    capital_trans[T+2] = fini_capital;
    capital_trans[3:(T+1)] .= fini_capital;
    V_saved_store = zeros(parameter_initial.ns_fine*3,T+2);
    V_saved_store[:,1] = V_saved_no_subsidy_reshaped;
    V_saved_store[:,T+2] = V_saved_subsidy_nb_reshaped; 
    if load_solution == 1
        dataframe_prices = CSV.read("vec_price_intro_subsidy_nb.csv", DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1;
    else
        price_trans = zeros(3,T+2);
        price_trans[1:2,end] = price_finish;
        price_trans[1:2,1] = price_start;
        price_trans[1:2,2] = price_start;
        fraction_guess = range(parameter_end.τ_S,parameter_initial.τ_S,length = T_end)/parameter_end.τ_S
        fraction_guess = repeat(fraction_guess,1,2)'
        price_trans[1:2,2:(T_end+1)] = fraction_guess.*repeat(price_trans[1:2,2],1,T_end)  + (1 .- fraction_guess).*repeat(price_trans[1:2,end],1,T_end)
        price_trans[1,(T_end+2):end - 1] .= price_finish[1];
        price_trans[2,(T_end+2):end - 1] .= price_finish[2];
        price_trans_actual = price_trans[:,2:T+1]
        vec_price_trans = copy(reshape(price_trans_actual,(T)*3));
    end
elseif case_final == 3
    #Effect of introducing the infrastructure program with aid finance
    price_start = copy(prices_no_subsidy);
    coeff_begin = copy(coeff_no_subsidy);
    init_distr = copy(stat_distr_no_subsidy);
    init_capital = copy(foreign_supply_capital_subsidy_b);
    price_finish = copy(prices_inf);
    coeff_end = copy(coeff_inf);
    fini_distr = copy(stat_distr_inf );
    fini_capital = copy(foreign_supply_capital_subsidy_b);
    T = 50; #Initial guess for the transition length
    T_end = 15; # Initial guess - convergence after
    parameter_initial = copy(No_subsidy_parameter);
    parameter_end = copy(infra_parameter_nsp_nb );
    ### Initialization based on the initial and final guess
    # Slow shutting down of subsidy rate
    T_tau_transition = 5; # 
    Q_S_trans = ones(T+2);
    Q_S_trans = parameter_end.Q_S * Q_S_trans[:];
    Q_S_trans[1:T_tau_transition] = range(parameter_initial.Q_S,parameter_end.Q_S,length = T_tau_transition);
    F_W_trans = parameter_end.F_W * ones(T+2);
    coeff_store = zeros(size(coeff_end,1),size(coeff_end,2),T+2);
    distr_store = zeros(size(init_distr,1),T+2);
    distr_store[:,1] = init_distr;
    distr_store[:,2] = init_distr;
    distr_store[:,T + 2] = fini_distr;

    coeff_store[:,:,1] = coeff_begin;
    coeff_store[:,:,T+2] = coeff_end;

    capital_trans = zeros(T+2);
    capital_trans[1] = init_capital;
    capital_trans[2] = init_capital;
    capital_trans[T+2] = fini_capital;
    capital_trans[3:(T+1)] .= fini_capital;
    V_saved_store = zeros(parameter_initial.ns_fine*3,T+2);
    V_saved_store[:,1] = V_saved_no_subsidy_reshaped;
    V_saved_store[:,T+2] = V_saved_inf_reshaped; 
    if load_solution == 1
        dataframe_prices = CSV.read("vec_price_intro_infra.csv", DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1;
    else
        price_trans = zeros(3,T+2);
        price_trans[1:2,end] = price_finish;
        price_trans[1:2,1] = price_start;
        price_trans[1:2,2] = price_start;
        fraction_guess = range(parameter_end.Q_S,parameter_initial.Q_S,length = T_end)/parameter_end.Q_S
        fraction_guess = repeat(fraction_guess,1,2)'
        price_trans[1:2,2:(T_end+1)] = fraction_guess.*repeat(price_trans[1:2,2],1,T_end)  + (1 .- fraction_guess).*repeat(price_trans[1:2,end],1,T_end)
        price_trans[1,(T_end+2):end - 1] .= price_finish[1];
        price_trans[2,(T_end+2):end - 1] .= price_finish[2];
        price_trans_actual = price_trans[:,2:T+1]
        vec_price_trans = copy(reshape(price_trans_actual,(T)*3));
    end
elseif case_final == 4
    #Effect of introducing the infrastructure program with aid finance
    price_start = copy(prices_no_subsidy);
    coeff_begin = copy(coeff_no_subsidy);
    init_distr = copy(stat_distr_no_subsidy);
    init_capital = copy(foreign_supply_capital_subsidy_b);
    price_finish = copy(prices_inf_sp);
    coeff_end = copy(coeff_inf_sp);
    fini_distr = copy(stat_distr_inf_sp );
    fini_capital = copy(foreign_supply_capital_subsidy_b);
    T = 50; #Initial guess for the transition length
    T_end = 15; # Initial guess - convergence after
    parameter_initial = copy(No_subsidy_parameter);
    parameter_end = copy(infra_parameter_sp_nb );
    ### Initialization based on the initial and final guess
    # Slow shutting down of subsidy rate
    T_tau_transition = 5; # 
    Q_S_trans = ones(T+2);
    Q_S_trans = parameter_end.Q_S * Q_S_trans[:];
    Q_S_trans[1:T_tau_transition] = range(parameter_initial.Q_S,parameter_end.Q_S,length = T_tau_transition);
    F_W_trans = ones(T+2);
    F_W_trans = parameter_end.F_W * F_W_trans[:];
    F_W_trans[1:T_tau_transition] = range(parameter_initial.F_W,parameter_end.F_W,length = T_tau_transition);
    coeff_store = zeros(size(coeff_end,1),size(coeff_end,2),T+2);
    distr_store = zeros(size(init_distr,1),T+2);
    distr_store[:,1] = init_distr;
    distr_store[:,2] = init_distr;
    distr_store[:,T + 2] = fini_distr;

    coeff_store[:,:,1] = coeff_begin;
    coeff_store[:,:,T+2] = coeff_end;

    capital_trans = zeros(T+2);
    capital_trans[1] = init_capital;
    capital_trans[2] = init_capital;
    capital_trans[T+2] = fini_capital;
    capital_trans[3:(T+1)] .= fini_capital;
    V_saved_store = zeros(parameter_initial.ns_fine*3,T+2);
    V_saved_store[:,1] = V_saved_no_subsidy_reshaped;
    V_saved_store[:,T+2] = V_saved_inf_sp_reshaped;
    if load_solution == 1
        dataframe_prices = CSV.read("vec_price_intro_infra_sp.csv", DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1;
    else
        price_trans = zeros(3,T+2);
        price_trans[1:2,end] = price_finish;
        price_trans[1:2,1] = price_start;
        price_trans[1:2,2] = price_start;
        fraction_guess = range(parameter_end.Q_S,parameter_initial.Q_S,length = T_end)/parameter_end.Q_S
        fraction_guess = repeat(fraction_guess,1,2)'
        price_trans[1:2,2:(T_end+1)] = fraction_guess.*repeat(price_trans[1:2,2],1,T_end)  + (1 .- fraction_guess).*repeat(price_trans[1:2,end],1,T_end)
        price_trans[1,(T_end+2):end - 1] .= price_finish[1];
        price_trans[2,(T_end+2):end - 1] .= price_finish[2];
        price_trans_actual = price_trans[:,2:T+1]
        vec_price_trans = copy(reshape(price_trans_actual,(T)*3));
    end
elseif case_final == 5

    #Effect of introducing the subsidy program with govt budget - no_subs to subsidy_b
    price_start = copy(prices_no_subsidy);
    coeff_begin = copy(coeff_no_subsidy);
    init_distr = copy(stat_distr_no_subsidy);
    init_capital = copy(foreign_supply_capital_subsidy_b);
    price_finish = copy(prices_subsidy_b_grid[:,τ_index]);
    coeff_end = copy(coeff_subsidy_b_grid[:,:,τ_index]);
    fini_distr = copy(stat_distr_subsidy_b_grid[:,τ_index] );
    fini_capital = copy(foreign_supply_capital_subsidy_b);
    T = 50; #Initial guess for the transition length
    T_end = 15; # Initial guess - convergence after
    parameter_initial = copy(No_subsidy_parameter);
    parameter_end = copy(Baseline_parameter );
    parameter_end.τ_S = τ_grid[τ_index]
    ### Initialization based on the initial and final guess
    # Slow introduction ofthe subsidy rate
    T_tau_transition = 5; # 
    τ_trans = ones(T+2);
    τ_trans = parameter_end.τ_S * τ_trans[:];
    τ_trans[1:T_tau_transition] = range(parameter_initial.τ_S,parameter_end.τ_S,length = T_tau_transition);

    coeff_store = zeros(size(coeff_end,1),size(coeff_end,2),T+2);
    distr_store = zeros(size(init_distr,1),T+2);
    distr_store[:,1] = init_distr;
    distr_store[:,2] = init_distr;
    distr_store[:,T + 2] = fini_distr;

    coeff_store[:,:,1] = coeff_begin;
    coeff_store[:,:,T+2] = coeff_end;

    capital_trans = zeros(T+2);
    capital_trans[1] = init_capital;
    capital_trans[2] = init_capital;
    capital_trans[T+2] = fini_capital;
    capital_trans[3:(T+1)] .= fini_capital;
    V_saved_store = zeros(parameter_initial.ns_fine*3,T+2);
    V_saved_store[:,1] = V_saved_no_subsidy_reshaped;
    V_saved_store[:,T+2] = V_saved_b_grid_reshaped[:,τ_index];
    if load_solution == 1
        name_file = string("vec_price_intro_subsidy", τ_index, ".csv")
        dataframe_prices = CSV.read(name_file, DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1;
    elseif load_solution == -1 
        # case1 guess
        dataframe_prices = CSV.read("vec_price_intro_subsidy.csv", DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1;
    elseif load_solution == -2 
        # previous guess
        name_prev_file = string("vec_price_intro_subsidy", τ_index-1, ".csv")
        dataframe_prices = CSV.read(name_prev_file, DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1; 
    elseif load_solution == -3 
        # next guess
        name_prev_file = string("vec_price_intro_subsidy", τ_index+1, ".csv")
        dataframe_prices = CSV.read(name_prev_file, DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1; 
    else
        price_trans = zeros(3,T+2);
        price_trans[1:3,end] = price_finish;
        price_trans[1:2,1] = price_start;
        price_trans[1:2,2] = price_start;
        fraction_guess = range(parameter_end.τ_S,parameter_initial.τ_S,length = T_end)/parameter_end.τ_S
        fraction_guess = repeat(fraction_guess,1,3)'
        fraction_guess_tau = range(price_finish[3]-0.01,0.0,length = T_tau_transition)/price_finish[3]
        fraction_guess[3,1:T_tau_transition] =fraction_guess_tau
        fraction_guess[3,(T_tau_transition+1):end] .= 0.0;
        price_trans[1:3,2:(T_end+1)] = fraction_guess.*repeat(price_trans[:,2],1,T_end)  + (1 .- fraction_guess).*repeat(price_trans[:,end],1,T_end)
        price_trans[1,(T_end+2):end - 1] .= price_finish[1];
        price_trans[2,(T_end+2):end - 1] .= price_finish[2];
        price_trans[3,(T_tau_transition+2):end - 1] .= price_finish[3];
        price_trans_actual = price_trans[:,2:T+1]
        vec_price_trans = copy(reshape(price_trans_actual,(T)*3));
    end
elseif case_final == 6
    # For this case, externalities.jl must be evaluated
    #no subsidy to subsidy with epsilon effects 
    price_start = copy(prices_epsilon_grid[:,epsilon_index]);
    coeff_begin = copy(coeff_epsilon_grid[:,:,epsilon_index]);
    init_distr = copy(stat_distr_epsilon_grid[:,epsilon_index]);
    init_capital = copy(foreign_supply_capital_subsidy_b);
    price_finish = copy(prices_subsidy_b);
    coeff_end = copy(coeff_subsidy_b);
    fini_distr = copy(stat_distr_subsidy_b);
    fini_capital = copy(foreign_supply_capital_subsidy_b);
    T = 50; #Initial guess for the transition length
    T_end = 15; # Initial guess - convergence after
    parameter_initial = copy(No_subsidy_parameter);
    parameter_end = copy(Baseline_parameter );
    epsilon_u = epsilon_grid[epsilon_index];
    epsilon_r = epsilon_grid[epsilon_index];
    ### Initialization based on the initial and final guess
    # Slow introduction ofthe subsidy rate
    T_tau_transition = 5; # 
    τ_trans = ones(T+2);
    τ_trans = parameter_end.τ_S * τ_trans[:];
    τ_trans[1:T_tau_transition] = range(parameter_initial.τ_S,parameter_end.τ_S,length = T_tau_transition);

    coeff_store = zeros(size(coeff_end,1),size(coeff_end,2),T+2);
    distr_store = zeros(size(init_distr,1),T+2);
    distr_store[:,1] = init_distr;
    distr_store[:,2] = init_distr;
    distr_store[:,T + 2] = fini_distr;

    coeff_store[:,:,1] = coeff_begin;
    coeff_store[:,:,T+2] = coeff_end;

    capital_trans = zeros(T+2);
    capital_trans[1] = init_capital;
    capital_trans[2] = init_capital;
    capital_trans[T+2] = fini_capital;
    capital_trans[3:(T+1)] .= fini_capital;
    V_saved_store = zeros(parameter_initial.ns_fine*3,T+2);
    V_saved_store[:,1] = V_saved_epsilon_grid_reshaped[:,epsilon_index];
    V_saved_store[:,T+2] = V_saved_subsidy_b_reshaped; 
    if load_solution == 1
        name_file = string("vec_price_subsidy_epsilon", epsilon_index, ".csv")
        dataframe_prices = CSV.read(name_file, DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1;
    elseif load_solution == -1 
        # case1 guess
        dataframe_prices = CSV.read("vec_price_intro_subsidy.csv", DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1;
        price_trans_actual = reshape(vec_price_trans,3,T);
        price_trans_actual1 = zeros(size(price_trans_actual)[1] + 1,size(price_trans_actual)[2])
        price_trans_actual1[1:3,:] = price_trans_actual
        price_trans_actual1[4,:] = range(price_start[3],undernourished_subsidy_b,50);
        price_trans_actual = price_trans_actual1;
        vec_price_trans = reshape(price_trans_actual,4*T);
    elseif load_solution == -2 
        # previous guess
        name_prev_file = string("vec_price_subsidy_epsilon", epsilon_index-1, ".csv")
        dataframe_prices = CSV.read(name_prev_file, DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1; 
    elseif load_solution == -3 
        # next guess
        name_prev_file = string("vec_price_subsidy_epsilon", epsilon_index+1, ".csv")
        dataframe_prices = CSV.read(name_prev_file, DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1; 
    else
        price_trans = zeros(3,T+2);
        price_trans[1:3,end] = price_finish;
        price_trans[1:2,1] = price_start;
        price_trans[1:2,2] = price_start;
        fraction_guess = range(parameter_end.τ_S,parameter_initial.τ_S,length = T_end)/parameter_end.τ_S
        fraction_guess = repeat(fraction_guess,1,3)'
        fraction_guess_tau = range(price_finish[3]-0.01,0.0,length = T_tau_transition)/price_finish[3]
        fraction_guess[3,1:T_tau_transition] =fraction_guess_tau
        fraction_guess[3,(T_tau_transition+1):end] .= 0.0;
        price_trans[1:3,2:(T_end+1)] = fraction_guess.*repeat(price_trans[:,2],1,T_end)  + (1 .- fraction_guess).*repeat(price_trans[:,end],1,T_end)
        price_trans[1,(T_end+2):end - 1] .= price_finish[1];
        price_trans[2,(T_end+2):end - 1] .= price_finish[2];
        price_trans[3,(T_tau_transition+2):end - 1] .= price_finish[3];
        price_trans_actual = price_trans[:,2:T+1]
        vec_price_trans = copy(reshape(price_trans_actual,(T)*3));
    end
end

# This might be different case by case!

# Constrain the per period price guesses
p_B_max = 2.0;
p_B_min = 0.1;
p_M_max = 0.4;
p_M_min = 0.1;
τ_W_min = 0.0;
τ_W_max = 0.5; 

# Functions

# Backward iteration step - solve for the policy functions given a price guess. Build in the spirit of solve_model_calibration

function local_var_creator_policy(ns::Int64,ns_fine::Int64,T::Int64,C_grid_fine_no::Int64)
    # Create storage of the NECESSARY policy functions for market clearing 
    
    a_prime_fine_store = zeros(ns_fine,3,T);
    future_occupation_fine_local_store = zeros(ns_fine,3,T);
    cons_fine_local_store = zeros(ns_fine,3,T);#Array{Float64}(undef, ns_tmp_fine,3,country_no,size_iid_cost_val,T);#
    # Policy functions:
    c_S_W_fine_store= zeros(ns_fine,3,T);
    c_B_W_fine_store= zeros(ns_fine,3,T);
    c_M_W_fine_store= zeros(ns_fine,3,T);
    Y_W_fine_policy_store = zeros(ns_fine,3,T);
    Y_S_fine_store = zeros(ns_fine,3,T);
    c_S_S_fine_store= zeros(ns_fine,3,T);
    c_B_S_fine_store= zeros(ns_fine,3,T);
    c_M_S_fine_store= zeros(ns_fine,3,T);
    q_S_S_fine_store= zeros(ns_fine,3,T);
    P_S_fine_store= zeros(ns_fine,3,T);
    x_S_S_fine_store= zeros(ns_fine,3,T);
    solve_staple_index_S_fine_store= zeros(Int64,ns_fine,3,T);
    λ_2_S_fine_store= zeros(ns_fine,3,T);
    future_asset_S_store= zeros(ns_fine,3,T);
    future_asset_C_store= zeros(ns_fine,3,T);
    c_S_B_fine_store= zeros(ns_fine,3,T);
    c_B_B_fine_store= zeros(ns_fine,3,T);
    c_M_B_fine_store= zeros(ns_fine,3,T);
    x_SC_fine_store= zeros(ns_fine,3,T);
    x_BC_fine_store= zeros(ns_fine,3,T);
    land_C_fine_store= zeros(ns_fine,3,T);
    λ_2_fine_store= zeros(ns_fine,3,T);
    P_B_fine_store= zeros(ns_fine,3,T);
    Y_B_fine_store= zeros(ns_fine,3,T);
    q_S_B_fine_store= zeros(ns_fine,3,T);
    q_B_B_fine_store= zeros(ns_fine,3,T);
    solve_cash_crop_index_B_fine_store= zeros(Int64,ns_fine,3,T);
    solve_staple_index_B_fine_store= zeros(Int64,ns_fine,3,T);
    TC_fine_store= zeros(ns_fine,3,T);
    residual_store = zeros(6,T);

    coeff_λ_2_cashcrop_residual_unconstrained_store = zeros(C_grid_fine_no,ns,T);
    coeff_λ_2_cashcrop_residual_constrained_store = zeros(C_grid_fine_no,ns,T);
    C_max_unconstrained_store = zeros(ns,T);
    C_max_constrained_store = zeros(ns,T);
    C_min_unconstrained_store = zeros(ns,T);
    C_min_constrained_store = zeros(ns,T);
    coeff_λ_2_s_store = zeros(C_grid_fine_no,ns,T);
    C_max_staple_store = zeros(ns,T);
    C_min_staple_store = zeros(ns,T);
    C_max_staple_constrained_store = zeros(ns,T);
    C_min_staple_constrained_store = zeros(ns,T);
    TC_S_c3_constrained_store = zeros(ns,T);
    x_S_c3_constrained_store = zeros(ns,T);
    q_S_c3_constrained_store = zeros(ns,T);
    c_S_c3_constrained_store = zeros(ns,T);


    return (   a_prime_fine_store,future_occupation_fine_local_store,
    cons_fine_local_store,c_S_W_fine_store,c_B_W_fine_store,c_M_W_fine_store,Y_W_fine_policy_store,Y_S_fine_store,
    c_S_S_fine_store,c_B_S_fine_store,c_M_S_fine_store,q_S_S_fine_store,P_S_fine_store,
    x_S_S_fine_store,solve_staple_index_S_fine_store,λ_2_S_fine_store,future_asset_S_store,
    future_asset_C_store,c_S_B_fine_store,c_B_B_fine_store,c_M_B_fine_store,x_SC_fine_store,x_BC_fine_store,
    land_C_fine_store,λ_2_fine_store,P_B_fine_store,Y_B_fine_store,q_S_B_fine_store,
    q_B_B_fine_store,solve_cash_crop_index_B_fine_store,solve_staple_index_B_fine_store,TC_fine_store,residual_store,
    coeff_λ_2_cashcrop_residual_unconstrained_store,coeff_λ_2_cashcrop_residual_constrained_store,C_max_unconstrained_store,
    C_max_constrained_store,C_min_unconstrained_store,C_min_constrained_store,coeff_λ_2_s_store,C_max_staple_store,C_min_staple_store,
    C_max_staple_constrained_store,C_min_staple_constrained_store,TC_S_c3_constrained_store,x_S_c3_constrained_store,q_S_c3_constrained_store,
    c_S_c3_constrained_store)
end

function local_var_creator_policy_details(ns::Int64,ns_fine::Int64,T::Int64,C_grid_fine_no::Int64)
    # Create storage of the NECESSARY policy functions for market clearing 
    
    prod_staple_store = zeros(T);
    prod_cashcrop_store = zeros(T);
    prod_manuf_store = zeros(T);
    asset_supply_store = zeros(T);
    current_worker_pop_store = zeros(T);
    current_staple_pop_store = zeros(T);
    current_cashcrop_pop_store = zeros(T);
    current_account_residual_store = zeros(T);
    staple_productivity_store = zeros(T);
    cashcrop_productivity_store = zeros(T);
    manuf_productivity_store = zeros(T);
    aggregate_consumption_store = zeros(T);
    relative_land_to_cashcrop_store = zeros(T);
    mean_land_share_staples_store = zeros(T);
    undernourished_store = zeros(T);

    fertilizer_use_store = zeros(T);
    APG_store = zeros(T);
    var_APland_store = zeros(T);
    var_MPX_store = zeros(T);
    avg_labor_prod_rural_store = zeros(T);
    avg_labor_prod_urban_store = zeros(T);
    avg_agri_prod_rural_store = zeros(T);
    avg_agri_prod_urban_store = zeros(T);
    return (prod_staple_store,prod_cashcrop_store,prod_manuf_store,asset_supply_store,current_worker_pop_store,
    current_staple_pop_store,current_cashcrop_pop_store,current_account_residual_store,staple_productivity_store,
    cashcrop_productivity_store,manuf_productivity_store,aggregate_consumption_store,relative_land_to_cashcrop_store,
    mean_land_share_staples_store,undernourished_store,fertilizer_use_store,APG_store,var_APland_store,var_MPX_store,
    avg_labor_prod_rural_store,avg_labor_prod_urban_store,avg_agri_prod_rural_store,avg_agri_prod_urban_store)
end


function Residual_transition_backward(s::Array{Float64,2},ns::Int64,
    z::Array{Float64,1},z_W::Array{Float64,1},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,
    p_B::Float64,p_M::Float64,ϕ_B::Float64,τ_B::Float64,ρ::Float64,w::Float64,r::Float64,
    c̄_S::Float64,a_min::Float64,a_max::Float64,γ::Float64,n::Array{Int64,1},κ::Float64,Q_S::Float64,
    ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,agrid_fine::Array{Float64,1},
    fspace_C_fine::Dict{Symbol,Any},agrid::Array{Float64,1},
    coefficients_next_tmp::Array{Float64,2},coeff_next::Array{Float64,2},C_grid_fine::Array{Float64,1},
    F_W::Float64,F_S::Float64,F_B::Float64,FM_W::Float64,FM_S::Float64,FM_B::Float64,P_kron::SparseMatrixCSC{Float64, Int64},Phi::SparseMatrixCSC{Float64, Int64},
    Phi_z::SparseMatrixCSC{Float64, Int64},β::Float64,
    fspace_a::Dict{Symbol,Any},σ::Float64,tol::Float64 = 1e-8)

    (θ,labor_prod,tol,P_W,Y_W,coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
    coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,λ_B_interior_c3a,x_SC_interior_c3a,x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,
    q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b, c_S_mat,c_B_mat,
    c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat, λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,q_S_staples,c_S_staples,c_B_staples,
    c_M_staples,P_S_staples,x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,C_max_unconstrained ,C_max_constrained,C_min_unconstrained,
    C_min_constrained,TC_mat,C_max_staple,C_min_staple,C_max_staple_constrained,C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,
    q_S_c3_constrained,c_S_c3_constrained,cbar_violated, x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c) = income_creator(s,ns,z,z_W,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,
        c̄_S,a_min,a_max,γ,n,κ,Q_S,ϵ,ψ_S,ψ_B,ψ_M,C_grid_fine,fspace_C_fine,agrid);
    #println("cbar: ", cbar_violated)


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

        x_tmp = zeros(ns,3);
        V_tmp = zeros(ns,3);
        conv = 10.0

    (coeff, conv, conv_ind)  = Bellman_iteration(coefficients_next_tmp,coeff_next,x_tmp,V_tmp,P_kron,Phi,min_C_applied,max_C_applied,Phi_z,β,
        fspace_a,σ,P_W,Y_W,s,r,ρ,w,coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,θ,
        fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns,κ,tol,a_min,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
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

function Residual_transition_forward(coeff::Array{Float64,2},distr_previous::Array{Float64,1},s_fine::Array{Float64,2},ns_fine::Int64,
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
    R::Float64,δ::Float64,α::Float64,tol::Float64 = 1e-8)
   
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

    # Aggregates
        (θ_fine,labor_prod_fine,tol,P_W_fine,Y_W_fine,coeff_λ_2_cashcrop_residual_unconstrained_fine,coeff_λ_2_cashcrop_residual_constrained_fine,
        x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
        coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,labor_allocated_interior_c3a_fine,
        λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,q_S_c2_fine,q_B_c1_fine,q_S_c3a_fine,
        q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,Y_B_c3b_fine, c_S_mat_fine,c_B_mat_fine,
        c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,q_B_mat_fine,land_B_mat_fine, λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,feasibility_mat_fine,C_max_mat_fine,
        C_min_mat_fine,q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,
        Y_S_potential_fine,TC_mat_fine,C_max_staple_fine,C_min_staple_fine,C_max_staple_constrained_fine,
        C_min_staple_constrained_fine,TC_S_c3_constrained_fine,x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine, x_S_mat_3c_fine,x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine) =  income_creator_no_approx(s_fine,ns_fine,
        z,z_W,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n_fine,κ,Q_S,ϵ,ψ_S,ψ_B,ψ_M,coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,
        C_max_unconstrained ,C_max_constrained,C_min_unconstrained,C_min_constrained, coeff_λ_2_s,C_grid_fine,fspace_C_fine,C_max_staple,C_min_staple,C_max_staple_constrained,
    C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained);


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
       Q_trans = spzeros(ns_fine*3,ns_fine*3);
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
        residual_goods[4] = 0;#(w - w_new)/(w+w_new);

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


function Residual_transition_forward_detailed(coeff::Array{Float64,2},distr_previous::Array{Float64,1},s_fine::Array{Float64,2},ns_fine::Int64,
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
    R::Float64,δ::Float64,α::Float64,tol::Float64 = 1e-8)
   
    #   coeff = coeff_store[:,:,t]; 
    #   distr_previous = distr_store[:,t]; 
    #   τ_S = τ_S_loc;
    #   coeff_λ_2_cashcrop_residual_unconstrained = coeff_λ_2_cashcrop_residual_unconstrained_store[:,:,t-1]
    #   coeff_λ_2_cashcrop_residual_constrained = coeff_λ_2_cashcrop_residual_constrained_store[:,:,t-1]
    #   C_max_unconstrained = C_max_unconstrained_store[:,t-1]
    #   C_max_constrained = C_max_constrained_store[:,t-1]
    #   C_min_unconstrained = C_min_unconstrained_store[:,t-1]
    #   C_min_constrained = C_min_constrained_store[:,t-1]
    #   coeff_λ_2_s = coeff_λ_2_s_store[:,:,t-1]
    #   C_max_staple = C_max_staple_store[:,t-1]
    #   C_min_staple = C_min_staple_store[:,t-1]
    #   C_max_staple_constrained = C_max_staple_constrained_store[:,t-1]
    #   C_min_staple_constrained = C_min_staple_constrained_store[:,t-1]
    #   TC_S_c3_constrained = TC_S_c3_constrained_store[:,t-1]
    #   x_S_c3_constrained = x_S_c3_constrained_store[:,t-1]
    #   q_S_c3_constrained = q_S_c3_constrained_store[:,t-1]
    #   c_S_c3_constrained = c_S_c3_constrained_store[:,t-1]
    #   foreign_supply_capital = capital_trans[t]

    # Aggregates
        (θ_fine,labor_prod_fine,tol,P_W_fine,Y_W_fine,coeff_λ_2_cashcrop_residual_unconstrained_fine,coeff_λ_2_cashcrop_residual_constrained_fine,
        x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
        coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,labor_allocated_interior_c3a_fine,
        λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,q_S_c2_fine,q_B_c1_fine,q_S_c3a_fine,
        q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,Y_B_c3b_fine, c_S_mat_fine,c_B_mat_fine,
        c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,q_B_mat_fine,land_B_mat_fine, λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,feasibility_mat_fine,C_max_mat_fine,
        C_min_mat_fine,q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,
        Y_S_potential_fine,TC_mat_fine,C_max_staple_fine,C_min_staple_fine,C_max_staple_constrained_fine,
        C_min_staple_constrained_fine,TC_S_c3_constrained_fine,x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine, x_S_mat_3c_fine,x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine) =  income_creator_no_approx(s_fine,ns_fine,
        z,z_W,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n_fine,κ,Q_S,ϵ,ψ_S,ψ_B,ψ_M,coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,
        C_max_unconstrained ,C_max_constrained,C_min_unconstrained,C_min_constrained, coeff_λ_2_s,C_grid_fine,fspace_C_fine,C_max_staple,C_min_staple,C_max_staple_constrained,
    C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained);


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
       Q_trans = spzeros(ns_fine*3,ns_fine*3);
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
        capital_used = copy(asset_supply) + foreign_supply_capital;
        Net_Foreign_Factor_Income = R*(copy(asset_supply) - capital_used)
        current_account_residual = Export_value - Import_value + Net_Foreign_Factor_Income;

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

        residual_goods[2] = (c_M_worker_sum + c_M_staple_sum + c_M_cashcrop_sum - prod_manuf)/(c_M_worker_sum + 
        c_M_staple_sum + c_M_cashcrop_sum + prod_manuf);

        #residual[2] = (manuf_cons + input_staple + input_cashcrop - prod_manuf)/(manuf_cons + input_staple + input_cashcrop + prod_manuf); # -- KAROL ADDED THIS
        residual_goods[3] = 0 #No longer needed for the calibration as the foreign_supply_capital offsets capital markets(r - r_new)/(r+r_new);
        residual_goods[4] = 0;#(w - w_new)/(w+w_new);

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
        cons_level_substinence = 0.109; # DONT FORGET to reset this when needed! 
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

function Residual_transition_sequential(price_trans_actual::Array{Float64,2},capital_trans::Array{Float64,1},
    distr_store::Array{Float64,2},T::Int64,parameter_end::Parameter_type,coeff_store::Array{Float64,3},τ_trans::Array{Float64,1},r::Float64)
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
        if (τ_S_loc==0 || prices_loc[3]==0.0)
            balanced_share = 0.0;
        else
            balanced_share = 1.0;
        end
        p_B,p_M,R,r,w,τ_W= price_reshaper_fixed_r(prices_loc,δ,ϵ,ψ_S,ψ_B,ψ_M,p_x,τ_S_loc,τ_B,α,balanced_share,r);
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
            c_S_c3_constrained_store[:,t-1] )= Residual_transition_backward(s,ns,
            z,z_W,ϕ_S,ζ,τ_S_loc,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n,κ,Q_S,ϵ,ψ_S,ψ_B,ψ_M,agrid_fine,
            fspace_C_fine,agrid,coefficients_next_tmp,coeff_next,C_grid_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,P_kron,Phi,Phi_z,β,fspace_a,σ)
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
            if (τ_S_loc==0 || prices_loc[3]==0.0)
                balanced_share = 0.0;
            else
                balanced_share = 1.0;
            end
            p_B,p_M,R,r,w,τ_W= price_reshaper_fixed_r(prices_loc,δ,ϵ,ψ_S,ψ_B,ψ_M,p_x,τ_S_loc,τ_B,α,balanced_share,r);
             (residual_store[:,t-1],distr_store[:,t]) = Residual_transition_forward(coeff_store[:,:,t],distr_store[:,t-1],s_fine,ns_fine,
                 z,z_W,ϕ_S,ζ,τ_S_loc,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n_fine,κ,Q_S,
                 ϵ,ψ_S,ψ_B,ψ_M,coeff_λ_2_cashcrop_residual_unconstrained_store[:,:,t-1],
                 coeff_λ_2_cashcrop_residual_constrained_store[:,:,t-1],
                 C_max_unconstrained_store[:,t-1] ,C_max_constrained_store[:,t-1],C_min_unconstrained_store[:,t-1],
                 C_min_constrained_store[:,t-1],coeff_λ_2_s_store[:,:,t-1],agrid_fine,fspace_C_fine,C_max_staple_store[:,t-1],
                 C_min_staple_store[:,t-1],C_max_staple_constrained_store[:,t-1],
                 C_min_staple_constrained_store[:,t-1],TC_S_c3_constrained_store[:,t-1],
                 x_S_c3_constrained_store[:,t-1],q_S_c3_constrained_store[:,t-1],c_S_c3_constrained_store[:,t-1],capital_trans[t],
                 balanced_share,τ_W,C_grid_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,P_kron_fine,Phi_z_fine,β,fspace_a,σ,P_kron1,fspace_a_fine,a_D,b_D,R,δ,α);
            print( t-1,",")
        end
    end
    return residual_store,distr_store,coeff_store
end


function Residual_transition_sequential_infra(price_trans_actual::Array{Float64,2},capital_trans::Array{Float64,1},
    distr_store::Array{Float64,2},T::Int64,parameter_end::Parameter_type,coeff_store::Array{Float64,3},Q_S_trans::Array{Float64,1},F_W_trans::Array{Float64,1},r::Float64)
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
        τ_S_loc = τ_S 
        Q_S_loc = Q_S_trans[t];
        F_W_loc = F_W_trans[t];
        prices_loc = price_trans_actual[:,t-1];
        if (τ_S_loc==0 || prices_loc[3]==0.0)
            balanced_share = 0.0;
        else
            balanced_share = 1.0;
        end
        p_B,p_M,R,r,w,τ_W= price_reshaper_fixed_r(prices_loc,δ,ϵ,ψ_S,ψ_B,ψ_M,p_x,τ_S_loc,τ_B,α,balanced_share,r);
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
            c_S_c3_constrained_store[:,t-1] )= Residual_transition_backward(s,ns,
            z,z_W,ϕ_S,ζ,τ_S_loc,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n,κ,Q_S_loc,ϵ,ψ_S,ψ_B,ψ_M,agrid_fine,
            fspace_C_fine,agrid,coefficients_next_tmp,coeff_next,C_grid_fine,F_W_loc,F_S,F_B,FM_W,FM_S,FM_B,P_kron,Phi,Phi_z,β,fspace_a,σ)
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
            τ_S_loc = τ_S 
            Q_S_loc = Q_S_trans[t];
            F_W_loc = F_W_trans[t];
            prices_loc = price_trans_actual[:,t-1];
            if (τ_S_loc==0 || prices_loc[3]==0.0)
                balanced_share = 0.0;
            else
                balanced_share = 1.0;
            end
            p_B,p_M,R,r,w,τ_W= price_reshaper_fixed_r(prices_loc,δ,ϵ,ψ_S,ψ_B,ψ_M,p_x,τ_S_loc,τ_B,α,balanced_share,r);
             (residual_store[:,t-1],distr_store[:,t]) = Residual_transition_forward(coeff_store[:,:,t],distr_store[:,t-1],s_fine,ns_fine,
                 z,z_W,ϕ_S,ζ,τ_S_loc,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n_fine,κ,Q_S_loc,
                 ϵ,ψ_S,ψ_B,ψ_M,coeff_λ_2_cashcrop_residual_unconstrained_store[:,:,t-1],
                 coeff_λ_2_cashcrop_residual_constrained_store[:,:,t-1],
                 C_max_unconstrained_store[:,t-1] ,C_max_constrained_store[:,t-1],C_min_unconstrained_store[:,t-1],
                 C_min_constrained_store[:,t-1],coeff_λ_2_s_store[:,:,t-1],agrid_fine,fspace_C_fine,C_max_staple_store[:,t-1],
                 C_min_staple_store[:,t-1],C_max_staple_constrained_store[:,t-1],
                 C_min_staple_constrained_store[:,t-1],TC_S_c3_constrained_store[:,t-1],
                 x_S_c3_constrained_store[:,t-1],q_S_c3_constrained_store[:,t-1],c_S_c3_constrained_store[:,t-1],capital_trans[t],
                 balanced_share,τ_W,C_grid_fine,F_W_loc,F_S,F_B,FM_W,FM_S,FM_B,P_kron_fine,Phi_z_fine,β,fspace_a,σ,P_kron1,fspace_a_fine,a_D,b_D,R,δ,α);
            print( t-1,",")
        end
    end
    return residual_store,distr_store,coeff_store
end

function Residual_transition_sequential_detailed(price_trans_actual::Array{Float64,2},capital_trans::Array{Float64,1},
    distr_store::Array{Float64,2},T::Int64,parameter_end::Parameter_type,coeff_store::Array{Float64,3},τ_trans::Array{Float64,1},r::Float64,V_saved_store::Array{Float64,2})
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
        if τ_S_loc==0
            balanced_share = 0.0;
        else
            balanced_share = 1.0;
        end
        p_B,p_M,R,r,w,τ_W= price_reshaper_fixed_r(prices_loc,δ,ϵ,ψ_S,ψ_B,ψ_M,p_x,τ_S_loc,τ_B,α,balanced_share,r);
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
            c_S_c3_constrained_store[:,t-1] )= Residual_transition_backward(s,ns,
            z,z_W,ϕ_S,ζ,τ_S_loc,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n,κ,Q_S,ϵ,ψ_S,ψ_B,ψ_M,agrid_fine,
            fspace_C_fine,agrid,coefficients_next_tmp,coeff_next,C_grid_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,P_kron,Phi,Phi_z,β,fspace_a,σ)
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
            if τ_S_loc==0
                balanced_share = 0.0;
            else
                balanced_share = 1.0;
            end
             p_B,p_M,R,r,w,τ_W= price_reshaper_fixed_r(prices_loc,δ,ϵ,ψ_S,ψ_B,ψ_M,p_x,τ_S_loc,τ_B,α,balanced_share,r);
            (residual_store[:,t-1],distr_store[:,t],prod_staple_store[t-1],prod_cashcrop_store[t-1],prod_manuf_store[t-1],asset_supply_store[t-1],
            current_worker_pop_store[t-1],current_staple_pop_store[t-1],current_cashcrop_pop_store[t-1],
            current_account_residual_store[t-1],staple_productivity_store[t-1],cashcrop_productivity_store[t-1],manuf_productivity_store[t-1],aggregate_consumption_store[t-1],relative_land_to_cashcrop_store[t-1],
            mean_land_share_staples_store[t-1],undernourished_store[t-1],
            fertilizer_use_store[t-1],APG_store[t-1],var_APland_store[t-1],var_MPX_store[t-1],avg_labor_prod_rural_store[t-1],avg_labor_prod_urban_store[t-1],avg_agri_prod_rural_store[t-1],avg_agri_prod_urban_store[t-1],
            V_saved_store[:,t],a_prime_fine_store[:,:,t-1]        ) = Residual_transition_forward_detailed(coeff_store[:,:,t],distr_store[:,t-1],s_fine,ns_fine,
            z,z_W,ϕ_S,ζ,τ_S_loc,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n_fine,κ,Q_S,
            ϵ,ψ_S,ψ_B,ψ_M,coeff_λ_2_cashcrop_residual_unconstrained_store[:,:,t-1],
            coeff_λ_2_cashcrop_residual_constrained_store[:,:,t-1],
            C_max_unconstrained_store[:,t-1] ,C_max_constrained_store[:,t-1],C_min_unconstrained_store[:,t-1],
            C_min_constrained_store[:,t-1],coeff_λ_2_s_store[:,:,t-1],agrid_fine,fspace_C_fine,C_max_staple_store[:,t-1],
            C_min_staple_store[:,t-1],C_max_staple_constrained_store[:,t-1],
            C_min_staple_constrained_store[:,t-1],TC_S_c3_constrained_store[:,t-1],
            x_S_c3_constrained_store[:,t-1],q_S_c3_constrained_store[:,t-1],c_S_c3_constrained_store[:,t-1],capital_trans[t],
            balanced_share,τ_W,C_grid_fine,F_W,F_S,F_B,FM_W,FM_S,FM_B,P_kron_fine,Phi_z_fine,β,fspace_a,σ,P_kron1,fspace_a_fine,a_D,b_D,R,δ,α);
            print( t-1,",")
        end
    end
    return (residual_store,distr_store,coeff_store,prod_staple_store,prod_cashcrop_store,prod_manuf_store,asset_supply_store,current_worker_pop_store,current_staple_pop_store,current_cashcrop_pop_store,
        current_account_residual_store,staple_productivity_store,cashcrop_productivity_store,manuf_productivity_store,aggregate_consumption_store,relative_land_to_cashcrop_store,
        mean_land_share_staples_store,undernourished_store,fertilizer_use_store,APG_store,var_APland_store,var_MPX_store,avg_labor_prod_rural_store,avg_labor_prod_urban_store,avg_agri_prod_rural_store,
        avg_agri_prod_urban_store,V_saved_store,a_prime_fine_store)
end
function Residual_transition_sequential_infra_detailed(price_trans_actual::Array{Float64,2},capital_trans::Array{Float64,1},
    distr_store::Array{Float64,2},T::Int64,parameter_end::Parameter_type,coeff_store::Array{Float64,3},Q_S_trans::Array{Float64,1},F_W_trans::Array{Float64,1},r::Float64,V_saved_store::Array{Float64,2})
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
       τ_S_loc = τ_S 
       Q_S_loc = Q_S_trans[t];
       F_W_loc = F_W_trans[t];
       prices_loc = price_trans_actual[:,t-1];
       if (τ_S_loc==0 || prices_loc[3]==0.0)
           balanced_share = 0.0;
       else
           balanced_share = 1.0;
       end
       p_B,p_M,R,r,w,τ_W= price_reshaper_fixed_r(prices_loc,δ,ϵ,ψ_S,ψ_B,ψ_M,p_x,τ_S_loc,τ_B,α,balanced_share,r);
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
           c_S_c3_constrained_store[:,t-1] )= Residual_transition_backward(s,ns,
           z,z_W,ϕ_S,ζ,τ_S_loc,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n,κ,Q_S_loc,ϵ,ψ_S,ψ_B,ψ_M,agrid_fine,
           fspace_C_fine,agrid,coefficients_next_tmp,coeff_next,C_grid_fine,F_W_loc,F_S,F_B,FM_W,FM_S,FM_B,P_kron,Phi,Phi_z,β,fspace_a,σ)
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
            τ_S_loc = τ_S 
            Q_S_loc = Q_S_trans[t];
            F_W_loc = F_W_trans[t];
            prices_loc = price_trans_actual[:,t-1];
            if (τ_S_loc==0 || prices_loc[3]==0.0)
                balanced_share = 0.0;
            else
                balanced_share = 1.0;
            end
            p_B,p_M,R,r,w,τ_W= price_reshaper_fixed_r(prices_loc,δ,ϵ,ψ_S,ψ_B,ψ_M,p_x,τ_S_loc,τ_B,α,balanced_share,r);
            (residual_store[:,t-1],distr_store[:,t],prod_staple_store[t-1],prod_cashcrop_store[t-1],prod_manuf_store[t-1],asset_supply_store[t-1],
            current_worker_pop_store[t-1],current_staple_pop_store[t-1],current_cashcrop_pop_store[t-1],
            current_account_residual_store[t-1],staple_productivity_store[t-1],cashcrop_productivity_store[t-1],manuf_productivity_store[t-1],aggregate_consumption_store[t-1],relative_land_to_cashcrop_store[t-1],
            mean_land_share_staples_store[t-1],undernourished_store[t-1],
            fertilizer_use_store[t-1],APG_store[t-1],var_APland_store[t-1],var_MPX_store[t-1],avg_labor_prod_rural_store[t-1],avg_labor_prod_urban_store[t-1],avg_agri_prod_rural_store[t-1],avg_agri_prod_urban_store[t-1],
            V_saved_store[:,t],a_prime_fine_store[:,:,t-1]        ) = Residual_transition_forward_detailed(coeff_store[:,:,t],distr_store[:,t-1],s_fine,ns_fine,
            z,z_W,ϕ_S,ζ,τ_S_loc,p_x,p_B,p_M,ϕ_B,τ_B,ρ,w,r,c̄_S,a_min,a_max,γ,n_fine,κ,Q_S_loc,
            ϵ,ψ_S,ψ_B,ψ_M,coeff_λ_2_cashcrop_residual_unconstrained_store[:,:,t-1],
            coeff_λ_2_cashcrop_residual_constrained_store[:,:,t-1],
            C_max_unconstrained_store[:,t-1] ,C_max_constrained_store[:,t-1],C_min_unconstrained_store[:,t-1],
            C_min_constrained_store[:,t-1],coeff_λ_2_s_store[:,:,t-1],agrid_fine,fspace_C_fine,C_max_staple_store[:,t-1],
            C_min_staple_store[:,t-1],C_max_staple_constrained_store[:,t-1],
            C_min_staple_constrained_store[:,t-1],TC_S_c3_constrained_store[:,t-1],
            x_S_c3_constrained_store[:,t-1],q_S_c3_constrained_store[:,t-1],c_S_c3_constrained_store[:,t-1],capital_trans[t],
            balanced_share,τ_W,C_grid_fine,F_W_loc,F_S,F_B,FM_W,FM_S,FM_B,P_kron_fine,Phi_z_fine,β,fspace_a,σ,P_kron1,fspace_a_fine,a_D,b_D,R,δ,α);
            print( t-1,",")
        end
    end
    return (residual_store,distr_store,coeff_store,prod_staple_store,prod_cashcrop_store,prod_manuf_store,asset_supply_store,current_worker_pop_store,current_staple_pop_store,current_cashcrop_pop_store,
    current_account_residual_store,staple_productivity_store,cashcrop_productivity_store,manuf_productivity_store,aggregate_consumption_store,relative_land_to_cashcrop_store,
    mean_land_share_staples_store,undernourished_store,fertilizer_use_store,APG_store,var_APland_store,var_MPX_store,avg_labor_prod_rural_store,avg_labor_prod_urban_store,avg_agri_prod_rural_store,
    avg_agri_prod_urban_store,V_saved_store,a_prime_fine_store)
end

function Residual_transition_iterative(vec_price_trans::Array{Float64,1},capital_trans::Array{Float64,1},distr_store::Array{Float64,2},
    T::Int64,parameter_end::Parameter_type,coeff_store::Array{Float64,3},τ_trans::Array{Float64,1})
    price_trans_actual = reshape(vec_price_trans,3,T);
    # price_trans_actual = copy(price_trans_next)
    residual_store,distr_store,coeff_store =  Residual_transition_sequential(price_trans_actual,capital_trans,distr_store,T,parameter_end,coeff_store,τ_trans,moments[1])
    maxval, maxval_ind = findmax(abs.(residual_store));
    println("Max residual: ",maxval," in time period: ",maxval_ind[2]," for market: ",maxval_ind[1])
    max_T = maxval_ind[2];
    price_trans_next = zeros(3,T);
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


    cross_terms = zeros(3,4); # This is where you can meddle with cross terms, for now, only contemporenous prices of the same markets matter
    cross_terms[1,1] =1;
    cross_terms[2,2] =1;
    cross_terms[3,3] =-1;
    cross_terms[1,4] =-1.0;
    cross_terms[2,4] =-1.0;
    cross_terms[3,4] = 1.0;
    # 4th column contains the staple market, which I add with a negative sign

    price_trans_next[1,:] = (ones(T)+cross_terms[1,1]*res_cc+cross_terms[1,2]*res_manu +
    cross_terms[1,3]*res_tax + cross_terms[1,4] * res_staple).*price_trans_actual[1,:];
    price_trans_next[2,:] =  (ones(T)+cross_terms[2,1]*res_cc+cross_terms[2,2]*res_manu +
    cross_terms[2,3]*res_tax+ cross_terms[2,4] * res_staple).*price_trans_actual[2,:];
    price_trans_next[3,:] =  (ones(T)+cross_terms[3,1]*res_cc+cross_terms[3,2]*res_manu +
    cross_terms[3,3]*res_tax+ cross_terms[3,4] * res_staple).*price_trans_actual[3,:];
    residual_adj = zeros(6,T);
    residual_adj[1,:] = res_cc;
    residual_adj[2,:] = res_manu;
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
    vec_price_trans_next_maxT = reshape(price_trans_next_maxT,3*T);
    vec_price_trans_next = reshape(price_trans_next,3*T);
    vec_price_trans_next_single_update_plus = reshape(price_trans_next_single_update_plus,3*T);
    return vec_price_trans_next,vec_price_trans_next_maxT,vec_price_trans_next_single_update_plus,residual_store,distr_store,coeff_store
end
function Residual_transition_iterative_infra(vec_price_trans::Array{Float64,1},capital_trans::Array{Float64,1},distr_store::Array{Float64,2},
    T::Int64,parameter_end::Parameter_type,coeff_store::Array{Float64,3},Q_S_trans::Array{Float64,1},F_W_trans::Array{Float64,1})
    price_trans_actual = reshape(vec_price_trans,3,T);
    # price_trans_actual = copy(price_trans_next)
    residual_store,distr_store,coeff_store =  Residual_transition_sequential_infra(price_trans_actual,capital_trans,distr_store,T,parameter_end,coeff_store,Q_S_trans,F_W_trans,moments[1])
    maxval, maxval_ind = findmax(abs.(residual_store));
    println("Max residual: ",maxval," in time period: ",maxval_ind[2]," for market: ",maxval_ind[1])
    max_T = maxval_ind[2];
    price_trans_next = zeros(3,T);
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


    cross_terms = zeros(3,4); # This is where you can meddle with cross terms, for now, only contemporenous prices of the same markets matter
    cross_terms[1,1] =1;
    cross_terms[2,2] =1;
    cross_terms[3,3] =-1;
    cross_terms[1,4] =-1.0;
    cross_terms[2,4] =-1.0;
    cross_terms[3,4] = 1.0;
    # 4th column contains the staple market, which I add with a negative sign

    price_trans_next[1,:] = (ones(T)+cross_terms[1,1]*res_cc+cross_terms[1,2]*res_manu +
    cross_terms[1,3]*res_tax + cross_terms[1,4] * res_staple).*price_trans_actual[1,:];
    price_trans_next[2,:] =  (ones(T)+cross_terms[2,1]*res_cc+cross_terms[2,2]*res_manu +
    cross_terms[2,3]*res_tax+ cross_terms[2,4] * res_staple).*price_trans_actual[2,:];
    price_trans_next[3,:] =  (ones(T)+cross_terms[3,1]*res_cc+cross_terms[3,2]*res_manu +
    cross_terms[3,3]*res_tax+ cross_terms[3,4] * res_staple).*price_trans_actual[3,:];
    residual_adj = zeros(6,T);
    residual_adj[1,:] = res_cc;
    residual_adj[2,:] = res_manu;
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
    vec_price_trans_next_maxT = reshape(price_trans_next_maxT,3*T);
    vec_price_trans_next = reshape(price_trans_next,3*T);
    vec_price_trans_next_single_update_plus = reshape(price_trans_next_single_update_plus,3*T);
    return vec_price_trans_next,vec_price_trans_next_maxT,vec_price_trans_next_single_update_plus,residual_store,distr_store,coeff_store
end
if generate_results == 0
    # The script to get a better price guess:
    dampen_start = 0.90;
    dampen_limit = 0.99;
    dampen = copy(dampen_start);
    #iter_end = 10;
    #dampen_grid = range(dampen,0.995,length = iter_end);
    #for iterate1 = iter_end:-1:50
    iterate1 = 1;
    conv_res_past = 2000;
    conv_res = 1;

    #vec_price_trans = copy(reshape(price_trans[:,2:(end-1)],(T)*3));
    vec_price_trans_smallest_res = copy(vec_price_trans);
    vec_price_trans_next_applied = copy(vec_price_trans);
    simple_search = 1;
    for i = 1:20
        if (case_final == 3 || case_final == 4 )
            (vec_price_trans_next,vec_price_trans_next_maxT,vec_price_trans_next_single_update_plus,residual_store,distr_store,coeff_store) = Residual_transition_iterative_infra(
                vec_price_trans,capital_trans,distr_store,T,parameter_end,coeff_store,Q_S_trans,F_W_trans);
        elseif (case_final == 6 )
            (vec_price_trans_next,vec_price_trans_next_maxT,vec_price_trans_next_single_update_plus,residual_store,distr_store,coeff_store) = Residual_transition_iterative_epsilon(
                vec_price_trans,capital_trans,distr_store,T,parameter_end,coeff_store,τ_trans,cons_level_substinence,epsilon_u,epsilon_r);
        else
            (vec_price_trans_next,vec_price_trans_next_maxT,vec_price_trans_next_single_update_plus,residual_store,distr_store,coeff_store) = Residual_transition_iterative(
                vec_price_trans,capital_trans,distr_store,T,parameter_end,coeff_store,τ_trans);
        end
        conv = maximum(abs.(vec_price_trans_next - vec_price_trans));
        conv_res = maximum(abs.(residual_store));
    #        conv = sum((vec_price_trans_next - vec_price_trans).^2);
    #        conv_res = sum((residual_store).^2);
        if simple_search == 1
            println("Simple search step")
            conv_res_past = copy(conv_res)
            vec_price_trans_smallest_res = copy(vec_price_trans);
            vec_price_trans_next_applied = copy(vec_price_trans_next);
        else
            if (conv_res<conv_res_past || dampen>dampen_limit)
                if conv_res<conv_res_past
                    println("Smaller residual found")
                    conv_res_past = copy(conv_res)
                    dampen = copy(dampen_start);
                    vec_price_trans_smallest_res = copy(vec_price_trans);
                    vec_price_trans_next_applied = copy(vec_price_trans_next);
        #             #vec_price_trans_next_applied = copy(vec_price_trans_next_maxT);
        #             #vec_price_trans_next_applied = copy(vec_price_trans_next_jac_based);
        #             #vec_price_trans_next_applied = copy(vec_price_trans_next_single_update_plus);
                elseif dampen>dampen_limit
                    println("Dampening too low")
                    vec_price_trans_next_single_update_plus_dampened =  (1 - dampen)*copy(vec_price_trans_next_single_update_plus)+dampen*vec_price_trans_smallest_res;
                    if (case_final == 3 || case_final == 4 )
                        (vec_price_trans_next1,residual_store1,vec_price_trans_next_maxT1, vec_price_trans_next_single_update_plus1)= Residual_transition_iterative_infra(vec_price_trans_next_single_update_plus_dampened,
                        capital_trans,distr_store,T,parameter_end,coeff_store,Q_S_trans,F_W_trans);
                    else
                        (vec_price_trans_next1,residual_store1,vec_price_trans_next_maxT1, vec_price_trans_next_single_update_plus1)= Residual_transition_iterative(vec_price_trans_next_single_update_plus_dampened,
                        capital_trans,distr_store,T,parameter_end,coeff_store,τ_trans);
                    end
                    conv_res_plus = maximum(abs.(residual_store1));

        # #                conv_res_minus = sum((residual_store2).^2);
                    if (conv_res_plus <conv_res)
                        println("Positive price updating")
                        conv_res_past = copy(conv_res_plus)
                        dampen = copy(dampen_start);
                        vec_price_trans_next_applied = copy(vec_price_trans_next_single_update_plus1);
                        vec_price_trans_smallest_res = copy(vec_price_trans_next_single_update_plus_dampened);
                    else
                        println("Move to a new guess")
                        conv_res_past = copy(conv_res)
                        dampen = copy(dampen_start);
                        vec_price_trans_next_applied = copy(vec_price_trans_next);
                        #vec_price_trans_next_applied = copy(vec_price_trans_next_single_update_plus_dampened);
                    end
                end
            else
                dampen = (1.0 + 5 * dampen)/6
                println("Current dampening: ",dampen)
            end
        end
        #vec_price_trans = (1 - dampen_grid[iterate1])*copy(vec_price_trans_next)+dampen_grid[iterate1]*vec_price_trans;
        vec_price_trans = (1 - dampen)*copy(vec_price_trans_next_applied)+dampen*vec_price_trans_smallest_res;
        #vec_price_trans = (1 - dampen)*copy(vec_price_trans_next)+dampen*vec_price_trans_smallest_res;
        println("Iteration:",iterate1," Prices:",conv," Residuals:",conv_res)
        iterate1 = iterate1 + 1;
    end
    if save_solution==1
        # This assumes that the main script has been evaluated
        if case_final == 1
            open("vec_price_intro_subsidy.csv", "a") do io
                writedlm(io, vec_price_trans_smallest_res,',')
            end
        elseif case_final == 2
            open("vec_price_intro_subsidy_nb.csv", "a") do io
                writedlm(io, vec_price_trans_smallest_res,',')
            end
        elseif case_final == 3
            open("vec_price_intro_infra.csv", "a") do io
                writedlm(io, vec_price_trans_smallest_res,',')
            end
        elseif case_final == 4
            open("vec_price_intro_infra_sp.csv", "a") do io
                writedlm(io, vec_price_trans_smallest_res,',')
            end
        elseif case_final == 5
            name_file = string("vec_price_intro_subsidy", τ_index, ".csv")
            open(name_file, "a") do io
                writedlm(io, vec_price_trans_smallest_res,',')
            end
        elseif case_final == 6
            name_file = string("vec_price_subsidy_epsilon", epsilon_index, ".csv")
            open(name_file, "a") do io
                writedlm(io, vec_price_trans_smallest_res,',')
            end
        end
    end
elseif generate_results == 1 
    if case_final == 1
        price_trans_actual_baseline = reshape(vec_price_trans,3,T);
         (residual_store_baseline,distr_store_baseline,coeff_store_baseline,prod_staple_store_baseline,prod_cashcrop_store_baseline,prod_manuf_store_baseline,asset_supply_store_baseline,
        current_worker_pop_store_baseline,current_staple_pop_store_baseline,current_cashcrop_pop_store_baseline,
        current_account_residual_store_baseline,staple_productivity_store_baseline,cashcrop_productivity_store_baseline,
        manuf_productivity_store_baseline, aggregate_consumption_store_baseline,relative_land_to_cashcrop_store_baseline,
         mean_land_share_staples_store_baseline,
        undernourished_store_baseline,fertilizer_use_store_baseline,APG_store_baseline,var_APland_store_baseline,
        var_MPX_store_baseline,avg_labor_prod_rural_store_baseline,avg_labor_prod_urban_store_baseline,
        avg_agri_prod_rural_store_baseline,avg_agri_prod_urban_store_baseline,V_saved_store_baseline,a_prime_fine_store_baseline) = Residual_transition_sequential_detailed(price_trans_actual_baseline,capital_trans,
    distr_store,T,parameter_end,coeff_store,τ_trans,moments[1],V_saved_store);

    Phi_fine_aug = kron(Matrix(1.0I, 3, 3),row_kron(funbase(Baseline_parameter.fspace_a,Baseline_parameter.s_fine[:,1]),Baseline_parameter.Phi_z_fine));
    welfare_val_tmp =  (Phi_fine_aug * vcat(vcat(coeff_store_baseline[:,1,2],coeff_store_baseline[:,2,2]),coeff_store_baseline[:,3,2]));
    
    #sum(distr_store_baseline[:,51]  .* (exp.((V_saved_store_baseline[:,51] -V_saved_store_baseline[:,1]) * (1.0 - Baseline_parameter.β) ) ))  - 1
    welfare_subsidy_b_trans = sum(distr_store_baseline[:,1] .* (exp.((welfare_val_tmp -welfare_val_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
    welfare_transition_subsidy_b = sum(distr_store_baseline[:,1]  .* (exp.((V_saved_store_baseline[:,2] -V_saved_store_baseline[:,1]) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget    
    welfare_transition_subsidy_b_alt = (1 - Baseline_parameter.β)* (sum(distr_store_baseline[:,1] .* V_saved_store_baseline[:,2]) - sum(distr_store_baseline[:,1].*V_saved_store_baseline[:,1])  ); # With balanced budget
    # Plot similar to optimal tau 

        staple_productivity_store_plotcorrection_baseline = staple_productivity_subsidy_b/ staple_productivity_store_baseline[end];
        staple_productivity_store_plotcorrection_baseline_range  = range(1,staple_productivity_store_plotcorrection_baseline,T)
        staple_productivity_store_smooth_baseline = 100 * (staple_productivity_store_plotcorrection_baseline_range .* staple_productivity_store_baseline[:]/staple_productivity_no_subsidy.-1);
        staple_productivity_store_smooth_baseline1 = 100* (staple_productivity_subsidy_b/staple_productivity_no_subsidy.-1) * ones(T+20)
        staple_productivity_store_smooth_baseline1[1:T] = staple_productivity_store_smooth_baseline;
        staple_productivity_store_smooth_baseline = copy(staple_productivity_store_smooth_baseline1);
        staple_productivity_store_smooth_baseline= movmean(staple_productivity_store_smooth_baseline,10);
        var_MPX_store_plotcorrection_baseline = var_MPX_subsidy_b/ var_MPX_store_baseline[end];
        var_MPX_store_plotcorrection_baseline_range  = range(1,var_MPX_store_plotcorrection_baseline,T)
        var_MPX_store_smooth_baseline = 100 * (var_MPX_store_plotcorrection_baseline_range .* var_MPX_store_baseline[:]/var_MPX_no_subsidy.-1);

        var_MPX_store_smooth_baseline1 = 100* (var_MPX_subsidy_b/var_MPX_no_subsidy.-1) * ones(T+20)
        var_MPX_store_smooth_baseline1[1:T] = var_MPX_store_smooth_baseline;
        var_MPX_store_smooth_baseline = copy(var_MPX_store_smooth_baseline1);
        var_MPX_store_smooth_baseline= movmean(var_MPX_store_smooth_baseline,10);

        APG_store_plotcorrection_baseline = APG_subsidy_b/ APG_store_baseline[end];
        APG_store_plotcorrection_baseline_range  = range(1,APG_store_plotcorrection_baseline,T)
        APG_store_smooth_baseline = 100 * (APG_store_plotcorrection_baseline_range .* APG_store_baseline[:]/APG_no_subsidy.-1);

        APG_store_smooth_baseline1 = 100* (APG_subsidy_b/APG_no_subsidy.-1) * ones(T+20)
        APG_store_smooth_baseline1[1:T] = APG_store_smooth_baseline;
        APG_store_smooth_baseline = copy(APG_store_smooth_baseline1);
        APG_store_smooth_baseline= movmean(APG_store_smooth_baseline,10);

        plot([staple_productivity_store_smooth_baseline
        ,var_MPX_store_smooth_baseline, APG_store_smooth_baseline]
        , label=[ "Staple productivity" "Dispersion of returns on fertilizer" "APG"  ],legend=:outerbottom,
        linewidth = 2,linestyle = [:solid :dash :dashdot :dot],ylims = [-2.0,150.0], xlabel = "Years",
        ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
        grid = false,size = (800,800),
        tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
        savefig("Figure_trans_5a.svg")


        current_worker_pop_store_plotcorrection_baseline = current_worker_pop_subsidy_b/ current_worker_pop_store_baseline[end];
        current_worker_pop_store_plotcorrection_baseline_range  = range(1,current_worker_pop_store_plotcorrection_baseline,T)
        current_worker_pop_store_smooth_baseline = 100 * (current_worker_pop_store_plotcorrection_baseline_range .* current_worker_pop_store_baseline[:]/current_worker_pop_no_subsidy.-1);        
        current_worker_pop_store_smooth_baseline1 = 100* (current_worker_pop_subsidy_b/current_worker_pop_no_subsidy.-1) * ones(T+20)
        current_worker_pop_store_smooth_baseline1[1:T] = current_worker_pop_store_smooth_baseline;
        current_worker_pop_store_smooth_baseline = copy(current_worker_pop_store_smooth_baseline1);
        current_worker_pop_store_smooth_baseline= movmean(current_worker_pop_store_smooth_baseline,10);


        avg_agri_prod_rural_store_plotcorrection_baseline = avg_agri_prod_rural_subsidy_b/ avg_agri_prod_rural_store_baseline[end];
        avg_agri_prod_rural_store_plotcorrection_baseline_range  = range(1,avg_agri_prod_rural_store_plotcorrection_baseline,T)
        avg_agri_prod_rural_store_smooth_baseline = 100 * (avg_agri_prod_rural_store_plotcorrection_baseline_range .* avg_agri_prod_rural_store_baseline[:]/avg_agri_prod_rural_no_subsidy.-1);
        avg_agri_prod_rural_store_smooth_baseline1 = 100* (avg_agri_prod_rural_subsidy_b/avg_agri_prod_rural_no_subsidy.-1) * ones(T+20)
        avg_agri_prod_rural_store_smooth_baseline1[1:T] = avg_agri_prod_rural_store_smooth_baseline;
        avg_agri_prod_rural_store_smooth_baseline = copy(avg_agri_prod_rural_store_smooth_baseline1);
        avg_agri_prod_rural_store_smooth_baseline= movmean(avg_agri_prod_rural_store_smooth_baseline,10);

        avg_labor_prod_urban_store_plotcorrection_baseline = avg_labor_prod_urban_subsidy_b/ avg_labor_prod_urban_store_baseline[end];
        avg_labor_prod_urban_store_plotcorrection_baseline_range  = range(1,avg_labor_prod_urban_store_plotcorrection_baseline,T)
        avg_labor_prod_urban_store_smooth_baseline = 100 * (avg_labor_prod_urban_store_plotcorrection_baseline_range .* avg_labor_prod_urban_store_baseline[:]/avg_labor_prod_urban_no_subsidy.-1);
        avg_labor_prod_urban_store_smooth_baseline1 = 100* (avg_labor_prod_urban_subsidy_b/avg_labor_prod_urban_no_subsidy.-1) * ones(T+20)
        avg_labor_prod_urban_store_smooth_baseline1[1:T] = avg_labor_prod_urban_store_smooth_baseline;
        avg_labor_prod_urban_store_smooth_baseline = copy(avg_labor_prod_urban_store_smooth_baseline1);
        avg_labor_prod_urban_store_smooth_baseline= movmean(avg_labor_prod_urban_store_smooth_baseline,10);

        plot([current_worker_pop_store_smooth_baseline
        ,avg_agri_prod_rural_store_smooth_baseline, avg_labor_prod_urban_store_smooth_baseline]
        , label=[ "Urbanization rate" "Mean rural ability" "Mean urban ability" ],legend=:outerbottom ,
        linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-20.0,5.0], xlabel = "Years",
        ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
        grid = false,size = (800,800),yticks = ([-20,-15,-10,-5,0,5 ],["-20","-15","-10","-5","0","5"]),
        tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
        savefig("Figure_trans_5b.svg")

        aggregate_consumption_store_plotcorrection_baseline = aggregate_consumption_subsidy_b/ aggregate_consumption_store_baseline[end];
        aggregate_consumption_store_plotcorrection_baseline_range  = range(1,aggregate_consumption_store_plotcorrection_baseline,T)
        aggregate_consumption_store_smooth_baseline = 100 * (aggregate_consumption_store_plotcorrection_baseline_range .* aggregate_consumption_store_baseline[:]/aggregate_consumption_no_subsidy.-1);
        aggregate_consumption_store_smooth_baseline1 = 100* (aggregate_consumption_subsidy_b/aggregate_consumption_no_subsidy.-1) * ones(T+20)
        aggregate_consumption_store_smooth_baseline1[1:T] = aggregate_consumption_store_smooth_baseline;
        aggregate_consumption_store_smooth_baseline = copy(aggregate_consumption_store_smooth_baseline1);
        aggregate_consumption_store_smooth_baseline= movmean(aggregate_consumption_store_smooth_baseline,10);

        undernourished_store_plotcorrection_baseline = undernourished_subsidy_b/ undernourished_store_baseline[end];
        undernourished_store_plotcorrection_baseline_range  = range(1,undernourished_store_plotcorrection_baseline,T)
        undernourished_store_smooth_baseline = 100 * (undernourished_store_plotcorrection_baseline_range .* undernourished_store_baseline[:]/undernourished_no_subsidy.-1);
        undernourished_store_smooth_baseline1 = 100* (undernourished_subsidy_b/undernourished_no_subsidy.-1) * ones(T+20)
        undernourished_store_smooth_baseline1[1:T] = undernourished_store_smooth_baseline;
        undernourished_store_smooth_baseline = copy(undernourished_store_smooth_baseline1);
        undernourished_store_smooth_baseline= movmean(undernourished_store_smooth_baseline,10);

        current_account_residual_store_plotcorrection_baseline = current_account_residual_subsidy_b/ current_account_residual_store_baseline[end];
        current_account_residual_store_plotcorrection_baseline_range  = range(1,current_account_residual_store_plotcorrection_baseline,T)
        current_account_residual_store_smooth_baseline = 100 * (current_account_residual_store_plotcorrection_baseline_range .* current_account_residual_store_baseline[:]/current_account_residual_no_subsidy.-1);
        current_account_residual_store_smooth_baseline1 = 100* (current_account_residual_subsidy_b/current_account_residual_no_subsidy.-1) * ones(T+20)
        current_account_residual_store_smooth_baseline1[1:T] = current_account_residual_store_smooth_baseline;
        current_account_residual_store_smooth_baseline = copy(current_account_residual_store_smooth_baseline1);
        current_account_residual_store_smooth_baseline= movmean(current_account_residual_store_smooth_baseline,10);

        plot([aggregate_consumption_store_smooth_baseline
        ,undernourished_store_smooth_baseline]
        , label=["Aggregate consumption" "Undernourished households"],legend=:outerbottom ,
        linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-34.0,5.0], xlabel = "Years",
        ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
        grid = false,size = (800,800),yticks = ([-30,-25,-20,-15,-10,-5,0,5 ],["-30","-25","-20","-15","-10","-5","0","5"]),
        tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
        savefig("Figure_trans_5c.svg")

        price_cc_plotcorrection_baseline = prices_subsidy_b[1]/ price_trans_actual_baseline[1,end];
        price_cc_plotcorrection_baseline_range  = range(1,price_cc_plotcorrection_baseline,T)
        price_cc_smooth_baseline = 100 * (price_cc_plotcorrection_baseline_range .* price_trans_actual_baseline[1,:]/prices_no_subsidy[1].-1);
        price_cc_smooth_baseline1 = 100* (prices_subsidy_b[1]/prices_no_subsidy[1].-1) * ones(T+20)
        price_cc_smooth_baseline1[1:T] = price_cc_smooth_baseline;
        price_cc_smooth_baseline = copy(price_cc_smooth_baseline1);
        price_cc_smooth_baseline= movmean(price_cc_smooth_baseline,10);

        price_manuf_plotcorrection_baseline = prices_subsidy_b[2]/ price_trans_actual_baseline[2,end];
        price_manuf_plotcorrection_baseline_range  = range(1,price_manuf_plotcorrection_baseline,T)
        price_manuf_smooth_baseline = 100 * (price_manuf_plotcorrection_baseline_range .* price_trans_actual_baseline[2,:]/prices_no_subsidy[2].-1);
        price_manuf_smooth_baseline1 = 100* (prices_subsidy_b[2]/prices_no_subsidy[2].-1) * ones(T+20)
        price_manuf_smooth_baseline1[1:T] = price_manuf_smooth_baseline;
        price_manuf_smooth_baseline = copy(price_manuf_smooth_baseline1);
        price_manuf_smooth_baseline= movmean(price_manuf_smooth_baseline,10);
        
        price_tau_W_plotcorrection_baseline = prices_subsidy_b[3]/ price_trans_actual_baseline[3,end];
        price_tau_W_plotcorrection_baseline_range  = range(1,price_tau_W_plotcorrection_baseline,T)
        price_tau_W_smooth_baseline = 100 * (price_tau_W_plotcorrection_baseline_range .* price_trans_actual_baseline[3,:]);
        price_tau_W_smooth_baseline1 = 100* prices_subsidy_b[3] * ones(T+20)
        price_tau_W_smooth_baseline1[1:T] = price_tau_W_smooth_baseline;
        price_tau_W_smooth_baseline = copy(price_tau_W_smooth_baseline1);
        price_tau_W_smooth_baseline= movmean(price_tau_W_smooth_baseline,10);

        plot([price_cc_smooth_baseline
        ,price_manuf_smooth_baseline, price_tau_W_smooth_baseline]
        , label=[ "Cashcrop price" "Manufacturing price" "Labor tax rate" ],legend=:outerbottom ,
        linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [0.0,55.0], xlabel = "Years",
        ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
        grid = false,size = (800,800),yticks = ([0,10,20,30,40,50 ],["0","10","20","30","40","50"]),
        tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
        savefig("Figure_trans_5d.svg")
        # Now get the DiD numbers
        # Fertilizer use
        fertilizer_use_subsidy_b = input_cashcrop_subsidy_b +   input_staple_subsidy_b ;
        fertilizer_use_no_subsidy = input_cashcrop_no_subsidy +   input_staple_no_subsidy;
        # 2010 level 34.13kg/ha
        fertilizer_kg_ha_2010 = 34.13;
        fertilizer_transformation_kg = fertilizer_kg_ha_2010/(fertilizer_use_subsidy_b/(mean_land_share_staples_subsidy_b/100));
        fertilizer_use_baseline_kg_ha_trans = fertilizer_transformation_kg*fertilizer_use_store_baseline./(mean_land_share_staples_store_baseline/100)
        fertilizer_use_baseline_kg_ha_subsidy_b = fertilizer_transformation_kg*fertilizer_use_subsidy_b./(mean_land_share_staples_subsidy_b/100)
        fertilizer_use_baseline_kg_ha_no_subsidy = fertilizer_transformation_kg*fertilizer_use_no_subsidy./(mean_land_share_staples_no_subsidy/100)
        fertilizer_use_kg_ha_smooth_baseline= movmean(fertilizer_use_baseline_kg_ha_trans,10);
        # Staple productivity 
        staple_productivity_store_plotcorrection_baseline = log.( staple_productivity_subsidy_b) - log.(staple_productivity_store_baseline[end]);
        staple_productivity_store_plotcorrection_baseline_range  = range(0,staple_productivity_store_plotcorrection_baseline,T)
        staple_productivity_store_smooth_baseline = (staple_productivity_store_plotcorrection_baseline_range .+  log.(staple_productivity_store_baseline[:]) .- log.(staple_productivity_no_subsidy));
        staple_productivity_store_smooth_baseline1 =  (log(staple_productivity_subsidy_b) - log(staple_productivity_no_subsidy)) * ones(T+20)
        staple_productivity_store_smooth_baseline1[1:T] = staple_productivity_store_smooth_baseline;
        staple_productivity_store_smooth_baseline = copy(staple_productivity_store_smooth_baseline1);
        staple_productivity_store_smooth_baseline= movmean(staple_productivity_store_smooth_baseline,10);

        cashcrop_productivity_store_plotcorrection_baseline = log.( cashcrop_productivity_subsidy_b) - log.(cashcrop_productivity_store_baseline[end]);
        cashcrop_productivity_store_plotcorrection_baseline_range  = range(0,cashcrop_productivity_store_plotcorrection_baseline,T)
        cashcrop_productivity_store_smooth_baseline = (cashcrop_productivity_store_plotcorrection_baseline_range .+  log.(cashcrop_productivity_store_baseline[:]) .- log.(cashcrop_productivity_no_subsidy));
        cashcrop_productivity_store_smooth_baseline1 =  (log(cashcrop_productivity_subsidy_b) - log(cashcrop_productivity_no_subsidy)) * ones(T+20)
        cashcrop_productivity_store_smooth_baseline1[1:T] = cashcrop_productivity_store_smooth_baseline;
        cashcrop_productivity_store_smooth_baseline = copy(cashcrop_productivity_store_smooth_baseline1);
        cashcrop_productivity_store_smooth_baseline= movmean(cashcrop_productivity_store_smooth_baseline,10);
        
        price_cc_plotcorrection_baseline = log(prices_subsidy_b[1])- log(price_trans_actual_baseline[1,end]);
        price_cc_plotcorrection_baseline_range  = range(0,price_cc_plotcorrection_baseline,T)
        price_cc_smooth_baseline =  (price_cc_plotcorrection_baseline_range .+ log.(price_trans_actual_baseline[1,:]) .- log.(prices_no_subsidy[1]));
        price_cc_smooth_baseline1 =  (log(prices_subsidy_b[1]) - log(prices_no_subsidy[1])) * ones(T+20)
        price_cc_smooth_baseline1[1:T] = price_cc_smooth_baseline;
        price_cc_smooth_baseline = copy(price_cc_smooth_baseline1);
        price_cc_smooth_baseline= movmean(price_cc_smooth_baseline,10);

        relative_land_to_cashcrop_store_plotcorrection_baseline = relative_land_to_cashcrop_subsidy_b - relative_land_to_cashcrop_store_baseline[end];
        relative_land_to_cashcrop_store_plotcorrection_baseline_range  = range(0,relative_land_to_cashcrop_store_plotcorrection_baseline,T)
        relative_land_to_cashcrop_store_smooth_baseline =  (relative_land_to_cashcrop_store_plotcorrection_baseline_range .+ relative_land_to_cashcrop_store_baseline .- relative_land_to_cashcrop_no_subsidy);
        relative_land_to_cashcrop_store_smooth_baseline1 = (relative_land_to_cashcrop_subsidy_b - relative_land_to_cashcrop_no_subsidy) * ones(T+20)
        relative_land_to_cashcrop_store_smooth_baseline1[1:T] = relative_land_to_cashcrop_store_smooth_baseline;
        relative_land_to_cashcrop_store_smooth_baseline = copy(relative_land_to_cashcrop_store_smooth_baseline1);
        relative_land_to_cashcrop_store_smooth_baseline= movmean(relative_land_to_cashcrop_store_smooth_baseline,10);

        undernourished_store_plotcorrection_baseline = undernourished_subsidy_b - undernourished_store_baseline[end];
        undernourished_store_plotcorrection_baseline_range  = range(0,undernourished_store_plotcorrection_baseline,T)
        undernourished_store_smooth_baseline =  (undernourished_store_plotcorrection_baseline_range .+ undernourished_store_baseline .- undernourished_no_subsidy);
        undernourished_store_smooth_baseline1 = (undernourished_subsidy_b - undernourished_no_subsidy) * ones(T+20)
        undernourished_store_smooth_baseline1[1:T] = undernourished_store_smooth_baseline;
        undernourished_store_smooth_baseline = copy(undernourished_store_smooth_baseline1);
        undernourished_store_smooth_baseline= movmean(undernourished_store_smooth_baseline,10);

        current_worker_pop_store_plotcorrection_baseline = current_worker_pop_subsidy_b - current_worker_pop_store_baseline[end];
        current_worker_pop_store_plotcorrection_baseline_range  = range(0,current_worker_pop_store_plotcorrection_baseline,T)
        current_worker_pop_store_smooth_baseline =  (current_worker_pop_store_plotcorrection_baseline_range .+ current_worker_pop_store_baseline .- current_worker_pop_no_subsidy);
        current_worker_pop_store_smooth_baseline1 = (current_worker_pop_subsidy_b - current_worker_pop_no_subsidy) * ones(T+20)
        current_worker_pop_store_smooth_baseline1[1:T] = current_worker_pop_store_smooth_baseline;
        current_worker_pop_store_smooth_baseline = copy(current_worker_pop_store_smooth_baseline1);
        current_worker_pop_store_smooth_baseline= movmean(current_worker_pop_store_smooth_baseline,10);
        
        mat_for_stata = zeros(11,7)
        mat_for_stata[:,1] = fertilizer_use_kg_ha_smooth_baseline[1:11].-fertilizer_use_baseline_kg_ha_no_subsidy
        mat_for_stata[:,2] = staple_productivity_store_smooth_baseline[1:11]
        mat_for_stata[:,3] = cashcrop_productivity_store_smooth_baseline[1:11]
        mat_for_stata[:,4] = price_cc_smooth_baseline[1:11]
        mat_for_stata[:,5] = relative_land_to_cashcrop_store_smooth_baseline[1:11]
        mat_for_stata[:,6] = undernourished_store_smooth_baseline[1:11]
        mat_for_stata[:,7] = -current_worker_pop_store_smooth_baseline[1:11]

        # Workers
        y_value = (exp.((V_saved_store_baseline[1:2560,2] -V_saved_store_baseline[1:2560,1]) * (1.0 - Baseline_parameter.β) ) ) .- 1
        y_mat = 100* mat_creator(y_value);
        plot([y_mat[1,1:5:end],y_mat[3,1:5:end],y_mat[13,1:5:end],y_mat[16,1:5:end]], label=["Low rural & urban \\theta" "High rural, low urban \\theta" "Low rural, high urban \\theta" "High rural & urban \\theta"]
        ,legend=:outerbottom ,
        linewidth = 2,linestyle = [:solid :dash :dot :dashdot],ylims = [-10.0,10.0], xlabel = "Wealth",
        ylabel = "Welfare change",xticks = ([5,14,18,21 ],["\$50","\$1000", "\$2500", "\$6000"]),yticks = ([-10,-5,0,5,10],["-10","-5","0","5","10"]),
        grid = false,size = (800, 800),legendcolumns=2,marker = [:none :none :circle :none],
        tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
        savefig("Figure2c.svg")
        #Farmers
        y_value = (exp.((V_saved_store_baseline[2561:(2*2560),2] -V_saved_store_baseline[2561:(2*2560),1]) * (1.0 - Baseline_parameter.β) ) ) .- 1
        y_mat = 100* mat_creator(y_value);
        plot([y_mat[1,1:5:end],y_mat[3,1:5:end],y_mat[13,1:5:end],y_mat[16,1:5:end]], label=["Low rural & urban \\theta" "High rural, low urban \\theta" "Low rural, high urban \\theta" "High rural & urban \\theta"]
        ,legend=:outerbottom ,
        linewidth = 2,linestyle = [:solid :dash :dot :dashdot],ylims = [-10.0,10.0], xlabel = "Wealth",
        ylabel = "Welfare change",xticks = ([5,14,18,21 ],["\$50","\$1000", "\$2500", "\$6000"]),yticks = ([-10,-5,0,5,10],["-10","-5","0","5","10"]),
        grid = false,size = (800, 800),legendcolumns=2,marker = [:none :none :circle :none],
        tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
        savefig("Figure2d.svg")

    elseif case_final == 2
        price_trans_actual_trans_nb = reshape(vec_price_trans,3,T);
         (residual_store_trans_nb,distr_store_trans_nb,coeff_store_trans_nb,prod_staple_store_trans_nb,prod_cashcrop_store_trans_nb,prod_manuf_store_trans_nb,asset_supply_store_trans_nb,
        current_worker_pop_store_trans_nb,current_staple_pop_store_trans_nb,current_cashcrop_pop_store_trans_nb,
        current_account_residual_store_trans_nb,staple_productivity_store_trans_nb,cashcrop_productivity_store_trans_nb,
        manuf_productivity_store_trans_nb, aggregate_consumption_store_trans_nb,relative_land_to_cashcrop_store_trans_nb,
         mean_land_share_staples_store_trans_nb,
        undernourished_store_trans_nb,fertilizer_use_store_trans_nb,APG_store_trans_nb,var_APland_store_trans_nb,
        var_MPX_store_trans_nb,avg_labor_prod_rural_store_trans_nb,avg_labor_prod_urban_store_trans_nb,
        avg_agri_prod_rural_store_trans_nb,avg_agri_prod_urban_store_trans_nb,V_saved_store_trans_nb,a_prime_fine_store_trans_nb) = Residual_transition_sequential_detailed(price_trans_actual_trans_nb,capital_trans,
    distr_store,T,parameter_end,coeff_store,τ_trans,moments[1],V_saved_store);
    Phi_fine_aug = kron(Matrix(1.0I, 3, 3),row_kron(funbase(Baseline_parameter.fspace_a,Baseline_parameter.s_fine[:,1]),Baseline_parameter.Phi_z_fine));
    welfare_val_tmp =  (Phi_fine_aug * vcat(vcat(coeff_store_trans_nb[:,1,2],coeff_store_trans_nb[:,2,2]),coeff_store_trans_nb[:,3,2]));
    
    #sum(distr_store[:,51]  .* (exp.((V_saved_store_trans_nb[:,51] -V_saved_store_trans_nb[:,1]) * (1.0 - Baseline_parameter.β) ) ))  - 1
    welfare_subsidy_nb_trans = sum(distr_store_trans_nb[:,1] .* (exp.((welfare_val_tmp -welfare_val_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
    welfare_transition_subsidy_nb = sum(distr_store_trans_nb[:,1]  .* (exp.((V_saved_store_trans_nb[:,2] -V_saved_store_trans_nb[:,1]) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget    
    welfare_transition_subsidy_nb_alt = (1 - Baseline_parameter.β)* (sum(distr_store_trans_nb[:,1] .* V_saved_store_trans_nb[:,2]) - sum(distr_store_trans_nb[:,1].*V_saved_store_trans_nb[:,1])  ); # With balanced budget
        # Plot similar to optimal tau 

        staple_productivity_store_plotcorrection_trans_nb = staple_productivity_subsidy_b/ staple_productivity_store_trans_nb[end];
        staple_productivity_store_plotcorrection_trans_nb_range  = range(1,staple_productivity_store_plotcorrection_trans_nb,T)
        staple_productivity_store_smooth_trans_nb = 100 * (staple_productivity_store_plotcorrection_trans_nb_range .* staple_productivity_store_trans_nb[:]/staple_productivity_no_subsidy.-1);
        staple_productivity_store_smooth_trans_nb1 = 100* (staple_productivity_subsidy_b/staple_productivity_no_subsidy.-1) * ones(T+20)
        staple_productivity_store_smooth_trans_nb1[1:T] = staple_productivity_store_smooth_trans_nb;
        staple_productivity_store_smooth_trans_nb = copy(staple_productivity_store_smooth_trans_nb1);
        staple_productivity_store_smooth_trans_nb= movmean(staple_productivity_store_smooth_trans_nb,10);
        var_MPX_store_plotcorrection_trans_nb = var_MPX_subsidy_b/ var_MPX_store_trans_nb[end];
        var_MPX_store_plotcorrection_trans_nb_range  = range(1,var_MPX_store_plotcorrection_trans_nb,T)
        var_MPX_store_smooth_trans_nb = 100 * (var_MPX_store_plotcorrection_trans_nb_range .* var_MPX_store_trans_nb[:]/var_MPX_no_subsidy.-1);

        var_MPX_store_smooth_trans_nb1 = 100* (var_MPX_subsidy_b/var_MPX_no_subsidy.-1) * ones(T+20)
        var_MPX_store_smooth_trans_nb1[1:T] = var_MPX_store_smooth_trans_nb;
        var_MPX_store_smooth_trans_nb = copy(var_MPX_store_smooth_trans_nb1);
        var_MPX_store_smooth_trans_nb= movmean(var_MPX_store_smooth_trans_nb,10);

        APG_store_plotcorrection_trans_nb = APG_subsidy_b/ APG_store_trans_nb[end];
        APG_store_plotcorrection_trans_nb_range  = range(1,APG_store_plotcorrection_trans_nb,T)
        APG_store_smooth_trans_nb = 100 * (APG_store_plotcorrection_trans_nb_range .* APG_store_trans_nb[:]/APG_no_subsidy.-1);

        APG_store_smooth_trans_nb1 = 100* (APG_subsidy_b/APG_no_subsidy.-1) * ones(T+20)
        APG_store_smooth_trans_nb1[1:T] = APG_store_smooth_trans_nb;
        APG_store_smooth_trans_nb = copy(APG_store_smooth_trans_nb1);
        APG_store_smooth_trans_nb= movmean(APG_store_smooth_trans_nb,10);

        plot([staple_productivity_store_smooth_trans_nb
        ,var_MPX_store_smooth_trans_nb, APG_store_smooth_trans_nb]
        , label=[ "Staple productivity" "Dispersion of returns on fertilizer" "APG"  ],legend=:outerbottom,
        linewidth = 2,linestyle = [:solid :dash :dashdot :dot],ylims = [-2.0,150.0], xlabel = "Years",
        ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
        grid = false,size = (800,800),
        tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
        savefig("Figure_trans_6a.svg")


        current_worker_pop_store_plotcorrection_trans_nb = current_worker_pop_subsidy_b/ current_worker_pop_store_trans_nb[end];
        current_worker_pop_store_plotcorrection_trans_nb_range  = range(1,current_worker_pop_store_plotcorrection_trans_nb,T)
        current_worker_pop_store_smooth_trans_nb = 100 * (current_worker_pop_store_plotcorrection_trans_nb_range .* current_worker_pop_store_trans_nb[:]/current_worker_pop_no_subsidy.-1);        
        current_worker_pop_store_smooth_trans_nb1 = 100* (current_worker_pop_subsidy_b/current_worker_pop_no_subsidy.-1) * ones(T+20)
        current_worker_pop_store_smooth_trans_nb1[1:T] = current_worker_pop_store_smooth_trans_nb;
        current_worker_pop_store_smooth_trans_nb = copy(current_worker_pop_store_smooth_trans_nb1);
        current_worker_pop_store_smooth_trans_nb= movmean(current_worker_pop_store_smooth_trans_nb,10);


        avg_agri_prod_rural_store_plotcorrection_trans_nb = avg_agri_prod_rural_subsidy_b/ avg_agri_prod_rural_store_trans_nb[end];
        avg_agri_prod_rural_store_plotcorrection_trans_nb_range  = range(1,avg_agri_prod_rural_store_plotcorrection_trans_nb,T)
        avg_agri_prod_rural_store_smooth_trans_nb = 100 * (avg_agri_prod_rural_store_plotcorrection_trans_nb_range .* avg_agri_prod_rural_store_trans_nb[:]/avg_agri_prod_rural_no_subsidy.-1);
        avg_agri_prod_rural_store_smooth_trans_nb1 = 100* (avg_agri_prod_rural_subsidy_b/avg_agri_prod_rural_no_subsidy.-1) * ones(T+20)
        avg_agri_prod_rural_store_smooth_trans_nb1[1:T] = avg_agri_prod_rural_store_smooth_trans_nb;
        avg_agri_prod_rural_store_smooth_trans_nb = copy(avg_agri_prod_rural_store_smooth_trans_nb1);
        avg_agri_prod_rural_store_smooth_trans_nb= movmean(avg_agri_prod_rural_store_smooth_trans_nb,10);

        avg_labor_prod_urban_store_plotcorrection_trans_nb = avg_labor_prod_urban_subsidy_b/ avg_labor_prod_urban_store_trans_nb[end];
        avg_labor_prod_urban_store_plotcorrection_trans_nb_range  = range(1,avg_labor_prod_urban_store_plotcorrection_trans_nb,T)
        avg_labor_prod_urban_store_smooth_trans_nb = 100 * (avg_labor_prod_urban_store_plotcorrection_trans_nb_range .* avg_labor_prod_urban_store_trans_nb[:]/avg_labor_prod_urban_no_subsidy.-1);
        avg_labor_prod_urban_store_smooth_trans_nb1 = 100* (avg_labor_prod_urban_subsidy_b/avg_labor_prod_urban_no_subsidy.-1) * ones(T+20)
        avg_labor_prod_urban_store_smooth_trans_nb1[1:T] = avg_labor_prod_urban_store_smooth_trans_nb;
        avg_labor_prod_urban_store_smooth_trans_nb = copy(avg_labor_prod_urban_store_smooth_trans_nb1);
        avg_labor_prod_urban_store_smooth_trans_nb= movmean(avg_labor_prod_urban_store_smooth_trans_nb,10);

        plot([current_worker_pop_store_smooth_trans_nb
        ,avg_agri_prod_rural_store_smooth_trans_nb, avg_labor_prod_urban_store_smooth_trans_nb]
        , label=[ "Urbanization rate" "Mean rural ability" "Mean urban ability" ],legend=:outerbottom ,
        linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-20.0,5.0], xlabel = "Years",
        ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
        grid = false,size = (800,800),yticks = ([-20,-15,-10,-5,0,5 ],["-20","-15","-10","-5","0","5"]),
        tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
        savefig("Figure_trans_6b.svg")

        aggregate_consumption_store_plotcorrection_trans_nb = aggregate_consumption_subsidy_b/ aggregate_consumption_store_trans_nb[end];
        aggregate_consumption_store_plotcorrection_trans_nb_range  = range(1,aggregate_consumption_store_plotcorrection_trans_nb,T)
        aggregate_consumption_store_smooth_trans_nb = 100 * (aggregate_consumption_store_plotcorrection_trans_nb_range .* aggregate_consumption_store_trans_nb[:]/aggregate_consumption_no_subsidy.-1);
        aggregate_consumption_store_smooth_trans_nb1 = 100* (aggregate_consumption_subsidy_b/aggregate_consumption_no_subsidy.-1) * ones(T+20)
        aggregate_consumption_store_smooth_trans_nb1[1:T] = aggregate_consumption_store_smooth_trans_nb;
        aggregate_consumption_store_smooth_trans_nb = copy(aggregate_consumption_store_smooth_trans_nb1);
        aggregate_consumption_store_smooth_trans_nb= movmean(aggregate_consumption_store_smooth_trans_nb,10);

        undernourished_store_plotcorrection_trans_nb = undernourished_subsidy_b/ undernourished_store_trans_nb[end];
        undernourished_store_plotcorrection_trans_nb_range  = range(1,undernourished_store_plotcorrection_trans_nb,T)
        undernourished_store_smooth_trans_nb = 100 * (undernourished_store_plotcorrection_trans_nb_range .* undernourished_store_trans_nb[:]/undernourished_no_subsidy.-1);
        undernourished_store_smooth_trans_nb1 = 100* (undernourished_subsidy_b/undernourished_no_subsidy.-1) * ones(T+20)
        undernourished_store_smooth_trans_nb1[1:T] = undernourished_store_smooth_trans_nb;
        undernourished_store_smooth_trans_nb = copy(undernourished_store_smooth_trans_nb1);
        undernourished_store_smooth_trans_nb= movmean(undernourished_store_smooth_trans_nb,10);

        current_account_residual_store_plotcorrection_trans_nb = current_account_residual_subsidy_b/ current_account_residual_store_trans_nb[end];
        current_account_residual_store_plotcorrection_trans_nb_range  = range(1,current_account_residual_store_plotcorrection_trans_nb,T)
        current_account_residual_store_smooth_trans_nb = 100 * (current_account_residual_store_plotcorrection_trans_nb_range .* current_account_residual_store_trans_nb[:]/current_account_residual_no_subsidy.-1);
        current_account_residual_store_smooth_trans_nb1 = 100* (current_account_residual_subsidy_b/current_account_residual_no_subsidy.-1) * ones(T+20)
        current_account_residual_store_smooth_trans_nb1[1:T] = current_account_residual_store_smooth_trans_nb;
        current_account_residual_store_smooth_trans_nb = copy(current_account_residual_store_smooth_trans_nb1);
        current_account_residual_store_smooth_trans_nb= movmean(current_account_residual_store_smooth_trans_nb,10);

        plot([aggregate_consumption_store_smooth_trans_nb
        ,undernourished_store_smooth_trans_nb]
        , label=["Aggregate consumption" "Undernourished households"],legend=:outerbottom ,
        linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-34.0,7.0], xlabel = "Years",
        ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
        grid = false,size = (800,800),yticks = ([-30,-25,-20,-15,-10,-5,0,5 ],["-30","-25","-20","-15","-10","-5","0","5"]),
        tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
        savefig("Figure_trans_6c.svg")

        price_cc_plotcorrection_trans_nb = prices_subsidy_b[1]/ price_trans_actual_trans_nb[1,end];
        price_cc_plotcorrection_trans_nb_range  = range(1,price_cc_plotcorrection_trans_nb,T)
        price_cc_smooth_trans_nb = 100 * (price_cc_plotcorrection_trans_nb_range .* price_trans_actual_trans_nb[1,:]/prices_no_subsidy[1].-1);
        price_cc_smooth_trans_nb1 = 100* (prices_subsidy_b[1]/prices_no_subsidy[1].-1) * ones(T+20)
        price_cc_smooth_trans_nb1[1:T] = price_cc_smooth_trans_nb;
        price_cc_smooth_trans_nb = copy(price_cc_smooth_trans_nb1);
        price_cc_smooth_trans_nb= movmean(price_cc_smooth_trans_nb,10);

        price_manuf_plotcorrection_trans_nb = prices_subsidy_b[2]/ price_trans_actual_trans_nb[2,end];
        price_manuf_plotcorrection_trans_nb_range  = range(1,price_manuf_plotcorrection_trans_nb,T)
        price_manuf_smooth_trans_nb = 100 * (price_manuf_plotcorrection_trans_nb_range .* price_trans_actual_trans_nb[2,:]/prices_no_subsidy[2].-1);
        price_manuf_smooth_trans_nb1 = 100* (prices_subsidy_b[2]/prices_no_subsidy[2].-1) * ones(T+20)
        price_manuf_smooth_trans_nb1[1:T] = price_manuf_smooth_trans_nb;
        price_manuf_smooth_trans_nb = copy(price_manuf_smooth_trans_nb1);
        price_manuf_smooth_trans_nb= movmean(price_manuf_smooth_trans_nb,10);
        
        plot([price_cc_smooth_trans_nb
        ,price_manuf_smooth_trans_nb]
        , label=[ "Cashcrop price" "Manufacturing price" ],legend=:outerbottom ,
        linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [0.0,55.0], xlabel = "Years",
        ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
        grid = false,size = (800,800),yticks = ([0,10,20,30,40,50 ],["0","10","20","30","40","50"]),
        tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
        savefig("Figure_trans_6d.svg")
    elseif case_final == 3
        price_trans_actual_trans_infra = reshape(vec_price_trans,3,T);
         (residual_store_trans_infra,distr_store_trans_infra,coeff_store_trans_infra,prod_staple_store_trans_infra,prod_cashcrop_store_trans_infra,prod_manuf_store_trans_infra,asset_supply_store_trans_infra,
        current_worker_pop_store_trans_infra,current_staple_pop_store_trans_infra,current_cashcrop_pop_store_trans_infra,
        current_account_residual_store_trans_infra,staple_productivity_store_trans_infra,cashcrop_productivity_store_trans_infra,
        manuf_productivity_store_trans_infra, aggregate_consumption_store_trans_infra,relative_land_to_cashcrop_store_trans_infra,
         mean_land_share_staples_store_trans_infra,
        undernourished_store_trans_infra,fertilizer_use_store_trans_infra,APG_store_trans_infra,var_APland_store_trans_infra,
        var_MPX_store_trans_infra,avg_labor_prod_rural_store_trans_infra,avg_labor_prod_urban_store_trans_infra,
        avg_agri_prod_rural_store_trans_infra,avg_agri_prod_urban_store_trans_infra,V_saved_store_trans_infra,a_prime_fine_store_trans_infra) =     Residual_transition_sequential_infra_detailed(price_trans_actual_trans_infra,capital_trans,
        distr_store,T,parameter_end,coeff_store,Q_S_trans,F_W_trans,moments[1],V_saved_store)


    Phi_fine_aug = kron(Matrix(1.0I, 3, 3),row_kron(funbase(Baseline_parameter.fspace_a,Baseline_parameter.s_fine[:,1]),Baseline_parameter.Phi_z_fine));
    welfare_val_tmp =  (Phi_fine_aug * vcat(vcat(coeff_store_trans_infra[:,1,2],coeff_store_trans_infra[:,2,2]),coeff_store_trans_infra[:,3,2]));
    
    #sum(distr_store[:,51]  .* (exp.((V_saved_store_trans_infra[:,51] -V_saved_store_trans_infra[:,1]) * (1.0 - Baseline_parameter.β) ) ))  - 1
    welfare_subsidy_infra_trans = sum(distr_store_trans_infra[:,1] .* (exp.((welfare_val_tmp -welfare_val_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
    welfare_transition_subsidy_infra = sum(distr_store_trans_infra[:,1]  .* (exp.((V_saved_store_trans_infra[:,2] -V_saved_store_trans_infra[:,1]) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget    
    welfare_transition_subsidy_infra_alt = (1 - Baseline_parameter.β)* (sum(distr_store_trans_infra[:,1] .* V_saved_store_trans_infra[:,2]) - sum(distr_store_trans_infra[:,1].*V_saved_store_trans_infra[:,1])  ); # With balanced budget
elseif case_final == 4
    price_trans_actual_trans_infra_sp = reshape(vec_price_trans,3,T);
     (residual_store_trans_infra_sp,distr_store_trans_infra_sp,coeff_store_trans_infra_sp,prod_staple_store_trans_infra_sp,prod_cashcrop_store_trans_infra_sp,prod_manuf_store_trans_infra_sp,asset_supply_store_trans_infra_sp,
    current_worker_pop_store_trans_infra_sp,current_staple_pop_store_trans_infra_sp,current_cashcrop_pop_store_trans_infra_sp,
    current_account_residual_store_trans_infra_sp,staple_productivity_store_trans_infra_sp,cashcrop_productivity_store_trans_infra_sp,
    manuf_productivity_store_trans_infra_sp, aggregate_consumption_store_trans_infra_sp,relative_land_to_cashcrop_store_trans_infra_sp,
     mean_land_share_staples_store_trans_infra_sp,
    undernourished_store_trans_infra_sp,fertilizer_use_store_trans_infra_sp,APG_store_trans_infra_sp,var_APland_store_trans_infra_sp,
    var_MPX_store_trans_infra_sp,avg_labor_prod_rural_store_trans_infra_sp,avg_labor_prod_urban_store_trans_infra_sp,
    avg_agri_prod_rural_store_trans_infra_sp,avg_agri_prod_urban_store_trans_infra_sp,V_saved_store_trans_infra_sp,a_prime_fine_store_trans_infra_sp) =     Residual_transition_sequential_infra_detailed(price_trans_actual_trans_infra_sp,capital_trans,
    distr_store,T,parameter_end,coeff_store,Q_S_trans,F_W_trans,moments[1],V_saved_store)


    Phi_fine_aug = kron(Matrix(1.0I, 3, 3),row_kron(funbase(Baseline_parameter.fspace_a,Baseline_parameter.s_fine[:,1]),Baseline_parameter.Phi_z_fine));
    welfare_val_tmp =  (Phi_fine_aug * vcat(vcat(coeff_store_trans_infra_sp[:,1,2],coeff_store_trans_infra_sp[:,2,2]),coeff_store_trans_infra_sp[:,3,2]));

    #sum(distr_store[:,51]  .* (exp.((V_saved_store_trans_infra[:,51] -V_saved_store_trans_infra[:,1]) * (1.0 - Baseline_parameter.β) ) ))  - 1
    welfare_subsidy_infra_sp_trans = sum(distr_store_trans_infra_sp[:,1] .* (exp.((welfare_val_tmp -welfare_val_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
    welfare_transition_subsidy_infra_sp = sum(distr_store_trans_infra_sp[:,1]  .* (exp.((V_saved_store_trans_infra_sp[:,2] -V_saved_store_trans_infra_sp[:,1]) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget    
    welfare_transition_subsidy_infra_sp_alt = (1 - Baseline_parameter.β)* (sum(distr_store_trans_infra_sp[:,1] .* V_saved_store_trans_infra_sp[:,2]) - sum(distr_store_trans_infra_sp[:,1].*V_saved_store_trans_infra_sp[:,1])  ); # With balanced budget

elseif case_final == 5

    price_trans_actual_optimal_tmp = reshape(vec_price_trans,3,T);
     (residual_store_optimal_tmp,distr_store_optimal_tmp,coeff_store_optimal_tmp,prod_staple_store_optimal_tmp,prod_cashcrop_store_optimal_tmp,prod_manuf_store_optimal_tmp,asset_supply_store_optimal_tmp,
    current_worker_pop_store_optimal_tmp,current_staple_pop_store_optimal_tmp,current_cashcrop_pop_store_optimal_tmp,
    current_account_residual_store_optimal_tmp,staple_productivity_store_optimal_tmp,cashcrop_productivity_store_optimal_tmp,
    manuf_productivity_store_optimal_tmp, aggregate_consumption_store_optimal_tmp,relative_land_to_cashcrop_store_optimal_tmp,
     mean_land_share_staples_store_optimal_tmp,
    undernourished_store_optimal_tmp,fertilizer_use_store_optimal_tmp,APG_store_optimal_tmp,var_APland_store_optimal_tmp,
    var_MPX_store_optimal_tmp,avg_labor_prod_rural_store_optimal_tmp,avg_labor_prod_urban_store_optimal_tmp,
    avg_agri_prod_rural_store_optimal_tmp,avg_agri_prod_urban_store_optimal_tmp,V_saved_store_optimal_tmp,a_prime_fine_store_optimal_tmp) = Residual_transition_sequential_detailed(price_trans_actual_optimal_tmp,capital_trans,
    distr_store,T,parameter_end,coeff_store,τ_trans,moments[1],V_saved_store);
    price_trans_actual_optimal[:,:,τ_index] = price_trans_actual_optimal_tmp;
    residual_store_optimal[:,:,τ_index] = residual_store_optimal_tmp
    distr_store_optimal[:,:,τ_index] = distr_store_optimal_tmp
    coeff_store_optimal[:,:,:,τ_index] = coeff_store_optimal_tmp;
    prod_staple_store_optimal[:,τ_index] = prod_staple_store_optimal_tmp
    prod_cashcrop_store_optimal[:,τ_index] = prod_cashcrop_store_optimal_tmp
    prod_manuf_store_optimal[:,τ_index] = prod_manuf_store_optimal_tmp
    asset_supply_store_optimal[:,τ_index] = asset_supply_store_optimal_tmp
    current_worker_pop_store_optimal[:,τ_index] = current_worker_pop_store_optimal_tmp
    current_staple_pop_store_optimal[:,τ_index] = current_staple_pop_store_optimal_tmp
    current_cashcrop_pop_store_optimal[:,τ_index] = current_cashcrop_pop_store_optimal_tmp
    current_account_residual_store_optimal[:,τ_index] = current_account_residual_store_optimal_tmp
    staple_productivity_store_optimal[:,τ_index] = staple_productivity_store_optimal_tmp
    cashcrop_productivity_store_optimal[:,τ_index] = cashcrop_productivity_store_optimal_tmp
    manuf_productivity_store_optimal[:,τ_index] = manuf_productivity_store_optimal_tmp
    aggregate_consumption_store_optimal[:,τ_index] = aggregate_consumption_store_optimal_tmp
    relative_land_to_cashcrop_store_optimal[:,τ_index] = relative_land_to_cashcrop_store_optimal_tmp
    mean_land_share_staples_store_optimal[:,τ_index] = mean_land_share_staples_store_optimal_tmp
    undernourished_store_optimal[:,τ_index] = undernourished_store_optimal_tmp
    fertilizer_use_store_optimal[:,τ_index] = fertilizer_use_store_optimal_tmp
    APG_store_optimal[:,τ_index] = APG_store_optimal_tmp
    var_APland_store_optimal[:,τ_index] = var_APland_store_optimal_tmp
    var_MPX_store_optimal[:,τ_index] = var_MPX_store_optimal_tmp
    avg_labor_prod_rural_store_optimal[:,τ_index] = avg_labor_prod_rural_store_optimal_tmp
    avg_labor_prod_urban_store_optimal[:,τ_index] = avg_labor_prod_urban_store_optimal_tmp
    avg_agri_prod_rural_store_optimal[:,τ_index] = avg_agri_prod_rural_store_optimal_tmp
    avg_agri_prod_urban_store_optimal[:,τ_index] = avg_agri_prod_urban_store_optimal_tmp
    V_saved_store_optimal[:,:,τ_index] = V_saved_store_optimal_tmp
    a_prime_fine_store_optimal[:,:,:,τ_index] = a_prime_fine_store_optimal_tmp

    Phi_fine_aug = kron(Matrix(1.0I, 3, 3),row_kron(funbase(Baseline_parameter.fspace_a,Baseline_parameter.s_fine[:,1]),Baseline_parameter.Phi_z_fine));
    welfare_val_tmp =  (Phi_fine_aug * vcat(vcat(coeff_store_optimal_tmp[:,1,2],coeff_store_optimal_tmp[:,2,2]),coeff_store_optimal_tmp[:,3,2]));

    welfare_trans_optimal_tmp = sum(distr_store_optimal_tmp[:,1] .* (exp.((welfare_val_tmp -welfare_val_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1;
    welfare_trans_optimal_real_tmp = sum(distr_store_optimal_tmp[:,1]  .* (exp.((V_saved_store_optimal_tmp[:,2] -V_saved_store_optimal_tmp[:,1]) * (1.0 - Baseline_parameter.β) ) ))  - 1;    
    welfare_trans_optimal_real_alt_tmp = (1 - Baseline_parameter.β)* (sum(distr_store_optimal_tmp[:,1] .* V_saved_store_optimal_tmp[:,2]) - sum(distr_store_optimal_tmp[:,1].*V_saved_store_optimal_tmp[:,1])  ); 

    welfare_trans_optimal[τ_index] = welfare_trans_optimal_tmp;
    welfare_trans_optimal_real[τ_index] = welfare_trans_optimal_real_tmp;
    welfare_trans_optimal_real_alt[τ_index] = welfare_trans_optimal_real_alt_tmp;
elseif case_final == 6

    price_trans_actual_epsilon_trans_tmp = reshape(vec_price_trans,4,T);
     (residual_store_epsilon_trans_tmp,distr_store_epsilon_trans_tmp,coeff_store_epsilon_trans_tmp,prod_staple_store_epsilon_trans_tmp,prod_cashcrop_store_epsilon_trans_tmp,prod_manuf_store_epsilon_trans_tmp,asset_supply_store_epsilon_trans_tmp,
    current_worker_pop_store_epsilon_trans_tmp,current_staple_pop_store_epsilon_trans_tmp,current_cashcrop_pop_store_epsilon_trans_tmp,
    current_account_residual_store_epsilon_trans_tmp,staple_productivity_store_epsilon_trans_tmp,cashcrop_productivity_store_epsilon_trans_tmp,
    manuf_productivity_store_epsilon_trans_tmp, aggregate_consumption_store_epsilon_trans_tmp,relative_land_to_cashcrop_store_epsilon_trans_tmp,
     mean_land_share_staples_store_epsilon_trans_tmp,
    undernourished_store_epsilon_trans_tmp,fertilizer_use_store_epsilon_trans_tmp,APG_store_epsilon_trans_tmp,var_APland_store_epsilon_trans_tmp,
    var_MPX_store_epsilon_trans_tmp,avg_labor_prod_rural_store_epsilon_trans_tmp,avg_labor_prod_urban_store_epsilon_trans_tmp,
    avg_agri_prod_rural_store_epsilon_trans_tmp,avg_agri_prod_urban_store_epsilon_trans_tmp,V_saved_store_epsilon_trans_tmp,a_prime_fine_store_epsilon_trans_tmp) = Residual_transition_sequential_epsilon_detailed(price_trans_actual_epsilon_trans_tmp,capital_trans,
    distr_store,T,parameter_end,coeff_store,τ_trans,moments[1],V_saved_store,cons_level_substinence,epsilon_u,epsilon_r);
    price_trans_actual_epsilon_trans[:,:,epsilon_index] = price_trans_actual_epsilon_trans_tmp;
    residual_store_epsilon_trans[:,:,epsilon_index] = residual_store_epsilon_trans_tmp
    distr_store_epsilon_trans[:,:,epsilon_index] = distr_store_epsilon_trans_tmp
    coeff_store_epsilon_trans[:,:,:,epsilon_index] = coeff_store_epsilon_trans_tmp;
    prod_staple_store_epsilon_trans[:,epsilon_index] = prod_staple_store_epsilon_trans_tmp
    prod_cashcrop_store_epsilon_trans[:,epsilon_index] = prod_cashcrop_store_epsilon_trans_tmp
    prod_manuf_store_epsilon_trans[:,epsilon_index] = prod_manuf_store_epsilon_trans_tmp
    asset_supply_store_epsilon_trans[:,epsilon_index] = asset_supply_store_epsilon_trans_tmp
    current_worker_pop_store_epsilon_trans[:,epsilon_index] = current_worker_pop_store_epsilon_trans_tmp
    current_staple_pop_store_epsilon_trans[:,epsilon_index] = current_staple_pop_store_epsilon_trans_tmp
    current_cashcrop_pop_store_epsilon_trans[:,epsilon_index] = current_cashcrop_pop_store_epsilon_trans_tmp
    current_account_residual_store_epsilon_trans[:,epsilon_index] = current_account_residual_store_epsilon_trans_tmp
    staple_productivity_store_epsilon_trans[:,epsilon_index] = staple_productivity_store_epsilon_trans_tmp
    cashcrop_productivity_store_epsilon_trans[:,epsilon_index] = cashcrop_productivity_store_epsilon_trans_tmp
    manuf_productivity_store_epsilon_trans[:,epsilon_index] = manuf_productivity_store_epsilon_trans_tmp
    aggregate_consumption_store_epsilon_trans[:,epsilon_index] = aggregate_consumption_store_epsilon_trans_tmp
    relative_land_to_cashcrop_store_epsilon_trans[:,epsilon_index] = relative_land_to_cashcrop_store_epsilon_trans_tmp
    mean_land_share_staples_store_epsilon_trans[:,epsilon_index] = mean_land_share_staples_store_epsilon_trans_tmp
    undernourished_store_epsilon_trans[:,epsilon_index] = undernourished_store_epsilon_trans_tmp
    fertilizer_use_store_epsilon_trans[:,epsilon_index] = fertilizer_use_store_epsilon_trans_tmp
    APG_store_epsilon_trans[:,epsilon_index] = APG_store_epsilon_trans_tmp
    var_APland_store_epsilon_trans[:,epsilon_index] = var_APland_store_epsilon_trans_tmp
    var_MPX_store_epsilon_trans[:,epsilon_index] = var_MPX_store_epsilon_trans_tmp
    avg_labor_prod_rural_store_epsilon_trans[:,epsilon_index] = avg_labor_prod_rural_store_epsilon_trans_tmp
    avg_labor_prod_urban_store_epsilon_trans[:,epsilon_index] = avg_labor_prod_urban_store_epsilon_trans_tmp
    avg_agri_prod_rural_store_epsilon_trans[:,epsilon_index] = avg_agri_prod_rural_store_epsilon_trans_tmp
    avg_agri_prod_urban_store_epsilon_trans[:,epsilon_index] = avg_agri_prod_urban_store_epsilon_trans_tmp
    V_saved_store_epsilon_trans[:,:,epsilon_index] = V_saved_store_epsilon_trans_tmp
    a_prime_fine_store_epsilon_trans[:,:,:,epsilon_index] = a_prime_fine_store_epsilon_trans_tmp

    Phi_fine_aug = kron(Matrix(1.0I, 3, 3),row_kron(funbase(Baseline_parameter.fspace_a,Baseline_parameter.s_fine[:,1]),Baseline_parameter.Phi_z_fine));
    welfare_val_tmp =  (Phi_fine_aug * vcat(vcat(coeff_store_epsilon_trans_tmp[:,1,2],coeff_store_epsilon_trans_tmp[:,2,2]),coeff_store_epsilon_trans_tmp[:,3,2]));

    welfare_trans_epsilon_trans_tmp = sum(distr_store_epsilon_trans_tmp[:,1] .* (exp.((welfare_val_tmp -welfare_val_no_subsidy) * (1.0 - Baseline_parameter.β) ) ))  - 1;
    welfare_trans_epsilon_trans_real_tmp = sum(distr_store_epsilon_trans_tmp[:,1]  .* (exp.((V_saved_store_epsilon_trans_tmp[:,2] -V_saved_store_epsilon_trans_tmp[:,1]) * (1.0 - Baseline_parameter.β) ) ))  - 1;    
    welfare_trans_epsilon_trans_real_alt_tmp = (1 - Baseline_parameter.β)* (sum(distr_store_epsilon_trans_tmp[:,1] .* V_saved_store_epsilon_trans_tmp[:,2]) - sum(distr_store_epsilon_trans_tmp[:,1].*V_saved_store_epsilon_trans_tmp[:,1])  ); 

    welfare_trans_epsilon_trans[epsilon_index] = welfare_trans_epsilon_trans_tmp;
    welfare_trans_epsilon_trans_real[epsilon_index] = welfare_trans_epsilon_trans_real_tmp;
    welfare_trans_epsilon_trans_real_alt[epsilon_index] = welfare_trans_epsilon_trans_real_alt_tmp;
    end
end