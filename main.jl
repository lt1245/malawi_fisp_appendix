

#cd("/Users/paroofa/Dropbox/Malawi_project/Model/Model_code/transaction_cost_v7 -Laszlo Copy/")
cd("/Users/bpu063010/Dropbox/Malawi_project/Model/Model_code/transaction_cost_v7 -Laszlo Copy/")
#cd("\\Users\\phbs\\Dropbox\\Malawi_project\\Model\\Model_code\\transaction_cost_v7 -Laszlo Copy\\")

#using MKL
using Roots
using MultistartOptimization
using DelimitedFiles
using Optim
using DataFrames
using CSV

using SharedArrays, Distributed
using CompEcon, QuantEcon
using LinearAlgebra, Statistics
using SparseArrays, StructArrays
using BasisMatrices, Arpack
#using BlackBoxOptim, NLopt
using BenchmarkTools
using LeastSquaresOptim
using GLM
using StatsPlots
using GeometryTypes
using StatsBase, GLM
#using LaTeXStrings

grid(ranges::NTuple{N, <: AbstractRange}) where N = Point.(Iterators.product(ranges...))

#choose if solve with balanced gov buget: balanced_share \in [0.0,1.0] governs how much of FISP expenditures is financed through labor tax.
balanced_share=1.0
balanced_share=min(max(balanced_share,0.0),1.0)
#τ_W_guess0=0.5 #initial tau_w guess for calibrations with balanced_share>0.0

mutable struct Parameter_type
    #Production parameters:
    K_a::Float64 # Constant in the demand for capital
    K_b::Float64 # Slope in the demand for capital
    δ::Float64 # Depreciation rate
    ζ::Float64 # CD parameter on intermediate inputs
    ρ::Float64 # CD parameter on additional labor inputs
    α::Float64 # Capital share in the manufacturing sector
    #Utility parameters
    σ::Float64 # Elasticity of substitution
    β::Float64 # Discount factor
    ϵ::Float64 # Intratemporal elasticity between staple, cash_crop and manufacturing goods
    ψ_S::Float64 # Utility weight on staple crop
    ψ_B::Float64 # Utility weight on cash crop
    ψ_M::Float64 # Utility weight on cash crop
    ϕ_S::Float64 # TFP of staple crop
    ϕ_B::Float64 # TFP of cash crop
    c̄_S::Float64 # Subsistence level of staple crop
    # Entry costs
    F_W::Float64 # Fixed cost to the manufacturing urban sector
    F_S::Float64 # Fixed cost to staple crop
    F_B::Float64 # Fixed cost to cash crop
    # Maintenance cost
    FM_W::Float64 # Fixed cost to the manufacturing urban sector
    FM_S::Float64 # Fixed cost to staple crop
    FM_B::Float64 # Fixed cost to cash crop
    Q_S::Float64 # Transaction cost for accessing food on top of home grown
    # Govt and external demand
    p_x::Float64 # price of the fertilizer, exogenously given
    τ_S::Float64 # Subsidy/tax on staple crop
    τ_B::Float64 # Subsidy/tax on cash crop
    a_D::Float64 # Constant in the demand for exports of cash crops
    b_D::Float64 # Constant in the demand for exports of cash crops
    #Financial markets
    γ::Float64 # Fraction of profits that can be borrowed
    # State space
    A_W::Float64 # Productivity shifter
    ρ_W::Float64 # autoregressive part of the transitory urban productivity shock
    ρ_S::Float64 # autoregressive part of the transitory rural productivity shock
    ρ_SW::Float64 # correlation between productivity shocks
    σ_S::Float64 # Standard deviation of the transitory rural productivity shock
    n::Array{Int64,1} # joint approximation number of gridpoints for 1) assets and 2) productivity
    n_fine::Array{Int64,1} # joint simulation number of gridpoints for 1) assets and 2) productivity
    agrid::Array{Float64,1} # grid for assets
    agrid_fine::Array{Float64,1} # finer grid for assets for the distribution
    a_min::Float64 # Lower bound for the asset grid/ borrowing constraint
    a_max::Float64 # Upper bound for the asset grid/ borrowing constraint
    spliorder::Int64 # Order of the spline approximation
    fspace_a::Dict{Symbol,Any} # function space for the approximation of the value function
    fspace_a_fine::Dict{Symbol,Any} # function space for the approximation of the distribution
    s::Array{Float64,2} #matrix of states, with the first column for assets, second for productivity/exogenous variable
    ns::Int64 #total number of states, taken as a cartesian product of "s"
    s_fine::Array{Float64,2} #matrix of states for the distribution, with the first column for assets, second for productivity
    ns_fine::Int64 #total number of states for the distribution, taken as a cartesian product of "s_fine"
    z::Array{Float64,1} # Productivity grid
    l_z::Array{Float64,1} # Labor productivity grid ---- TO BE REMOVED here!!!
    Phi_z::SparseMatrixCSC{Float64,Int64} #basis matrix for exogenous variables
    Phi_z_fine::SparseMatrixCSC{Float64,Int64}  #basis matrix for exogenous variables in the distribution
    Phi::SparseMatrixCSC{Float64,Int64}  #joint basis matrix
    Phi_aug::SparseMatrixCSC{Float64,Int64}  #joint basis matrix with occupation choice
    P_kron::SparseMatrixCSC{Float64,Int64}
    P_kron1::SparseMatrixCSC{Float64,Int64}
    P_kron_fine::SparseMatrixCSC{Float64,Int64}
    κ::Float64
    #l_z_low::Float64
    σ_W::Float64 # Standard deviation of the transitory urban productivity shock
    z_W::Array{Float64,1} # Urban Productivity grid
    #Phi_z_W::SparseMatrixCSC{Float64,Int64} #basis matrix for exogenous variables (urban)
    #Phi_z_W_fine::SparseMatrixCSC{Float64,Int64}  #basis matrix for exogenous variables in the distribution (urban)
    no_labor_shocks::Int64
    C_grid_fine_no::Int64
    C_grid_fine::Array{Float64,1}
    fspace_C_fine::Dict{Symbol,Any} 
end
#Redefine copy to allow for copyiing the parameter object created.
Base.copy(x::T) where T = T([deepcopy(getfield(x, k)) for k ∈ fieldnames(T)]...)

include("functions.jl")
#include("solve_model.jl")
#include("solve_model_exports.jl")
include("solve_model_calibration1.jl")
include("solve_model_calibration2.jl")
include("details_model.jl")
#include("solve_model_nomanu_price.jl")

# Initialize the parameters
Simple_function_space = fundefn(:spli, 10, 0.1, 1.0,1);
Baseline_parameter = Parameter_type(0.1,-1.0,0.06,1/3,1/3,1/3,2.0,0.86,2.0,1/3,1/3,1/3,1.0,1.5,0.1,
0.3,0.0,0.2,0.0,0.0,0.0,0.05,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.85,0.92,0.0,0.045,[50,36],[200,36],
zeros(2),zeros(2),0.0005,10.0,1,Simple_function_space,Simple_function_space,
zeros(2,2),2,zeros(2,2),2,zeros(2),zeros(2),spzeros(2,2),spzeros(2,2),spzeros(2,2),spzeros(2,2),
spzeros(2,2),spzeros(2,2),spzeros(2,2),0.5,0.4,zeros(2),2,2,zeros(2),Simple_function_space);
Baseline_parameter.n = [40,16]
Baseline_parameter.n_fine = [40,16]
Baseline_parameter.no_labor_shocks=sqrt(Baseline_parameter.n[2])
if mod(Baseline_parameter.n[2],Baseline_parameter.no_labor_shocks)!=0
    error("n & no_labor_shocks not compatible")
end

Baseline_parameter.γ = 0.0;
if Baseline_parameter.γ == 0.0
    Baseline_parameter.a_min = 0.01
    #Baseline_parameter.a_min = 0.0
    Baseline_parameter.a_max = 50.0
end

curve = 0.2;
curve_fine  = 0.2;
Baseline_parameter.agrid = range(Baseline_parameter.a_min^curve,Baseline_parameter.a_max^curve,length=
    Baseline_parameter.n[1]).^(1/curve);

##########################################################################################################################################################
##
##          SET MOMENTS:
##
##########################################################################################################################################################



#Targets:
r = 0.05;
L = 0.2; #[0.2]Malawi urban pop share
K_Y_ratio = 3.84 #K/Y ratio = 3.84 from [UN, 2014], if we dont use K/Yurban
Y_agr_Y_tot = 0.3 #share of agriculture in GDP Y_agr/Y_tot=30% (World Bank)
exported_cash_crop_ratio = 0.73 # [0.65,0.8]Share of 2010 exported cash crops in total prod. of cash crops = 73% (from FAOStat for Malawi)
fraction_staple_producers_without_surplus = 0.9 #Share of staple producers that do not market their staples = 90%, i.e. with c_s<=q_s
G_Y_ratio = 0.03 #[0.03]fiscal cost of FISP = 3% of GDP target with tau_S "Agri Input subsidies (book) Chirwa et al"
Urban_unemp_rate = 0.36 #Unemployment rate in urban [de Janvry et al. 2020]
RCT_moment1_value = 0.11 #RCT evidence in Daidone_AJAE_CashTransfersRCTEvaluation.pdf on impact of cash grant equiv. to 20-25% of their consumption to ultra poor on value of agricultural production. I think we need to define ultra-poor arbitrary (see the real definition in line below) - e.g. bottom 10% of asset or staple-consumption.
RCT_moment2_share = 0.072
#In Malawian case a household is classified as "ultrapoor" if any of the following conditions is verified:
# (a) the household has an average of only one meal a day,
# (b) the household survives from begging,
# (c) the household is undernourished,
# (d) the household does not possess any valuable assets, and
#(e) the household does not receive any monetary help, food, or gifts from others.
# Laszlo's assessment: a-c are related to consumption, only d is related to assets. My proposal: order people based on their consumption
exp_ratio = 2.7 #Ratio of total expenditures (inc subsidies) of cash crop producers to staple producers = 2.0
RU_migration_rate = 0.01 #[0.01,0.03]avg btw migration rate rural-to-urban=0.9% and urban-to-rural 2.8% in 2010-2013 panel ( = (migrants from rural going to urban) / (rural pop) ) 0.8% from [Lagakos et al. 2020]

#For moment12 I have two proposals (probably could use both of them in overidentifying calibr, or pick the one that works best for us):
fraction_only_cashcrops = 0.06 #6% of rural HHs producing only cash crops. From LSMS micro-data.
mean_land_share_to_staples_among_cc = 0.3 #[0.25,0.35]30% is the mean land share devoted to staples among cash crop producers. From LSMS micro-data.

rural_pop_only_staples = 0.41 #[0.4,0.5]From micro-data: 41% of rural pop cultivates only staples (i.e. 50% of 80% of country's pop living in rural cultivates only staples, which includes staple and cash crop farmers) ---- not sure if this will be easy to hit...
urban_rural_inc_ratio = 2.44 #Ratio of income in rural vs urban = 2.44 (=2,795/1,142)   [Santaeulalia&Magalhaes, 2018]
# Laszlo additional distribution momemts
urban_rural_wealth_ratio = 3.03; #(3976/1309)
urban_rural_consumption_ratio = 2.238998482549317;#(2951/1318)

# # Top 10% share - overall:
# Top10_share_consumption = 0.34
# Top10_share_income = 0.48
# Top10_share_wealth = 0.58
# # Top 1% share - overall:
# Top1_share_consumption = 0.08
# Top1_share_income = 0.18
# Top1_share_wealth = 0.25

# Top 10% share - rural:
Top10_share_consumption_rural = 0.30
Top10_share_income_rural = 0.44
Top10_share_wealth_rural = 0.49
# Top 1% share - rural:
Top1_share_consumption_rural = 0.06
Top1_share_income_rural = 0.15
Top1_share_wealth_rural = 0.17
# Top 10% share - urban:
Top10_share_consumption_urban = 0.34
Top10_share_income_urban = 0.57
Top10_share_wealth_urban = 0.73
# Top 1% share - urban:
Top1_share_consumption_urban = 0.06
Top1_share_income_urban = 0.19
Top1_share_wealth_urban = 0.35
# # Var-Logs
# var_log_cons_rural=0.41
# var_log_cons_urban=0.55
# var_log_inc_rural=0.98
# var_log_inc_urban=1.56
# var_log_wealth_rural=1.49
# var_log_wealth_urban=4.52
fraction_cashcrop_suboptimal = 0.6 #[0.6,0.85]from Brune et al EDCC 2016 Footnote 9.
APG_data= 4.0 #[4,7]

moments = [r, L, K_Y_ratio, Y_agr_Y_tot, exported_cash_crop_ratio, fraction_staple_producers_without_surplus, G_Y_ratio,
    RCT_moment1_value, RCT_moment2_share, exp_ratio, RU_migration_rate, fraction_only_cashcrops, mean_land_share_to_staples_among_cc,
    rural_pop_only_staples, urban_rural_inc_ratio, urban_rural_wealth_ratio, urban_rural_consumption_ratio,
    Top1_share_wealth_rural, Top1_share_income_rural, Top1_share_consumption_rural, Top10_share_wealth_rural, Top10_share_income_rural, Top10_share_consumption_rural,
    Top1_share_wealth_urban, Top1_share_income_urban, Top1_share_consumption_urban, Top10_share_wealth_urban, Top10_share_income_urban, Top10_share_consumption_urban,
    0.025, 0.025, 0.025,fraction_cashcrop_suboptimal,APG_data];


##########################################################################################################################################################
##
##          SET PARAMETERS:
##
##########################################################################################################################################################

#Set parameters:
Baseline_parameter.ψ_S = 0.124; # US share of exp on food.
Baseline_parameter.ψ_M = 0.8;
Baseline_parameter.ψ_B = 1.0 - Baseline_parameter.ψ_M - Baseline_parameter.ψ_S; # assumed low share of exp on alochol, coffee, ciggarates, some of clothes (cotton)
Baseline_parameter.ϵ=0.95 #from [Herrendorf et al. 2013] although transferability of this parameter may be low for us, as we model preferences over staples vs cash crops vs manuf, and they focus on agr vs services vs manuf.
Baseline_parameter.ζ=0.15 #cost share estimate from FAO data (see DB folder in Calibr/)
Baseline_parameter.Q_S=2.0#half of transaction costs documented in Santaeulalia Llopis & Maghaleas 2018
#Baseline_parameter.ρ_S = 0.7525; # Using the AR1 calculations based on the rural panel of LSMS 2010-2013: 3 year rho_S coeff=0.426 & sigma_S=0.472. Annualization implies rho^(1/3)
Baseline_parameter.ρ_S = 0.5745; # Using the AR1 calculations based on the rural panel of LSMS 2010-2013: 3 year rho_S coeff=0.19. Annualization implies rho^(1/3)
#Baseline_parameter.ρ_W=0.6028;# Using the AR1 calculations based on the urban panel of LSMS 2010-2013: 3 year rho_W coeff=0.219 & sigma_W=1.0293. Annualization implies rho^(1/3)
Baseline_parameter.ρ_W=0.4932;# Using the AR1 calculations based on the urban panel of LSMS 2010-2013: 3 year rho_W coeff=0.12 & sigma_W=1.0293. Annualization implies rho^(1/3)
#Baseline_parameter.σ_S = 0.3436; # Using the AR1 calculations based on the rural panel of LSMS 2010-2013: 3 year sigma_S=0.472. Annualization implies sigma/(1 + rho^2 + rho^4)^(1/2)
Baseline_parameter.σ_S = 0.9418; # Using the AR1 calculations based on the rural panel of LSMS 2010-2013: 3 year sigma_S=1.13. Annualization implies sigma/(1 + rho^2 + rho^4)^(1/2)
#Baseline_parameter.σ_W = 0.8417; # Using the AR1 calculations based on the urban panel of LSMS 2010-2013: 3 year sigma_W=1.0293. Annualization implies sigma/(1 + rho^2 + rho^4)^(1/2)
Baseline_parameter.σ_W = 1.1128; # Using the AR1 calculations based on the urban panel of LSMS 2010-2013: 3 year sigma_W=1.27. Annualization implies sigma/(1 + rho^2 + rho^4)^(1/2)
Baseline_parameter.σ=1.0 #log preferences
Baseline_parameter.α = 0.4 # share of capital in manufacturing output

Baseline_parameter.τ_B=0.0 #no subsidy for cash crops at baseline
Baseline_parameter.β = 0.85; # 0.85 BKS2021 ; 0.9 Lowest in Donovan, following th logic of BKS2021 with r>0
Baseline_parameter.b_D = -0.45;
Baseline_parameter.FM_W=0.0;
Baseline_parameter.F_B = 0.0;
Baseline_parameter.p_x = 1.26; # from multiple data sources (see DB folder Calibr/)
#Because i guess alternative way to approx p_x is from FOC: p_x=zeta*y/x. We found zeta=0.1 in the LSMS data already (btw this was only for fertilizer, maybe we want to add seed cost there too). Then, with maize price 0.23USD/kg and 35kg of fert applied per ha, empirically we’d have p_x=0.1*0.23*2200/35=1.44.
#If we also think about seeds, then we could use 0.15 estimated above (coz it seems to me like cost share estimation), and use p_x=0.15*0.23*2200/(35+25)=1.26.
Baseline_parameter.δ=0.05;
Baseline_parameter.ϕ_B = 1.0;
Baseline_parameter.ϕ_S = 1.0;

# Internally calibrated parameters - starting values
Baseline_parameter.κ = 0.1; # working capital constraint
Baseline_parameter.a_D = 0.8593341451145834;
Baseline_parameter.c̄_S = 0.01975466941261332; # Will use an algorithm to set the maximal value 
Baseline_parameter.F_W = 285;#305.0;
Baseline_parameter.FM_B = 2.2199277084690983;#0.85;
Baseline_parameter.ρ = 0.7412814048806258;
Baseline_parameter.ρ_SW = 0.23999102235984843;
Baseline_parameter.A_W = 2.85;
Baseline_parameter.τ_S=-0.8099992851842415;


Baseline_parameter.agrid_fine = a_grid_fine_gen_midpoints(Baseline_parameter.agrid,
Baseline_parameter.a_min,Baseline_parameter.a_max,4,Baseline_parameter.n);
Baseline_parameter.n_fine[1] = size(Baseline_parameter.agrid_fine)[1];

Baseline_parameter.C_grid_fine = a_grid_fine_gen_midpoints(Baseline_parameter.agrid,
Baseline_parameter.a_min,Baseline_parameter.a_max,4,Baseline_parameter.n);
Baseline_parameter.C_grid_fine_no = size(Baseline_parameter.C_grid_fine)[1];
Baseline_parameter.fspace_C_fine = fundef((:spli, Baseline_parameter.C_grid_fine, 0,1));

Baseline_parameter.fspace_a = fundef((:spli, Baseline_parameter.agrid, 0,Baseline_parameter.spliorder
        ));# define function space for the approximation of the solution of vfi.
Baseline_parameter.fspace_a_fine = fundef((:spli, Baseline_parameter.agrid_fine, 0,1
        ));# define function space for the approximation of the stationary distribution.
#Baseline_parameter.σ_W = 0.01

(Baseline_parameter.s,Baseline_parameter.ns,Baseline_parameter.s_fine,
Baseline_parameter.ns_fine,Baseline_parameter.Phi_z,Baseline_parameter.Phi_z_fine,Baseline_parameter.Phi,
Baseline_parameter.Phi_aug,Baseline_parameter.P_kron,Baseline_parameter.P_kron1,Baseline_parameter.P_kron_fine,Baseline_parameter.z,
Baseline_parameter.z_W) = setup_state_space(Baseline_parameter);
if Baseline_parameter.σ==1.0
    function curr_util(cons_guess::Array{Float64,1},σ::Float64)
        return @fastmath log.(cons_guess);
    end
else
    function curr_util(cons_guess::Array{Float64,1},σ::Float64)
        return @fastmath cons_guess.^(1 - σ)./(1 - σ);
    end
end


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
        coeff_subsidy_b,residual_goods_subsidy_b,model_moments_subsidy_b,foreign_supply_capital_subsidy_b)  = solve_model_calibration1(prices_subsidy_b,parameters_tmp,2,moments,balanced_share);
    prices_subsidy_b_grid[:,iterate] = prices_subsidy_b
    residual_goods_subsidy_b_grid[:,iterate] = residual_goods_subsidy_b
    iterate = iterate + 1;
end
prices_subsidy_b_grid[1,no_taus] =  1.1882831767847994;
prices_subsidy_b_grid[2,no_taus] =   0.18349666056314162;
prices_subsidy_b_grid[3,no_taus] = 0.00;

prices_subsidy_b_grid= [1.4598772874688932 1.4385245014811625 1.4049741835258556 1.4077138037158219 1.3821166752679364 1.3651971203006512 1.3573364085808683 1.363686641369129 1.3414500347358624 1.3323226063843923 1.327875051226767 1.3206902627537573 1.31038149959868 1.299209888144449 1.295273367696625 1.286440471180219 1.280855187630534 1.2693788918619269 1.2505170635332994 1.3058767212231221 1.2528814375260455 1.2390136405180543 1.2366561653512491 1.2397711208427493 1.2369018943039145 1.2292651896985871 1.2234888879703227 1.223637192043632 1.2123437075134387 1.1761; 0.2519999236690381 0.23498396688130332 0.23332979499142487 0.2331762054176376 0.23381710497342015 0.23254851966653814 0.22970890394934568 0.21321983055377886 0.21472526959909455 0.21520463011066363 0.21512099960052214 0.20424937174792618 0.2030998391253458 0.21684048568198006 0.20338551740011912 0.2031107078028406 0.20457532622931957 0.20671971878014053 0.2078049839228001 0.20138832443176954 0.2059963050532268 0.20902560956191182 0.20847460704850557 0.19109415674986044 0.1903407252421793 0.19041346055148342 0.20864525308669363 0.18988145820182262 0.1898044796482799 0.192054; 0.2049632173803936 0.1411151956716539 0.1143093243650973 0.1027914440951735 0.09459381800661454 0.07766669772376758 0.0722107637510424 0.06298170571272088 0.0547728674793511 0.047763531492882316 0.04623130159584552 0.03929669248628405 0.033648141106931 0.03384699432359155 0.027747870442189175 0.023731464417786858 0.02272952750763354 0.0207361287258126 0.018317258948674173 0.01628515625104149 0.014790446150075959 0.01215828872230139 0.010102326008558985 0.00967354325689979 0.00756756417746829 0.005680704866949134 0.005586786834627149 0.003926480520014699 0.0018768700814345938 0.0];
##########################################################################################################################################################
##
##          SAVING RESULTS FOR ANALYSIS:
## #! indicates unsolved/inaccurate results
##########################################################################################################################################################

# Details from the initial steady state - evaluate from here if want to obtain only the results
foreign_supply_capital_subsidy_b = 7.631597488627073;
#No subsidy equilibrium
prices_no_subsidy = [ 1.188283177632008,
0.18349666252335484];
No_subsidy_parameter = copy(Baseline_parameter);
No_subsidy_parameter.τ_S=-0;
(residual_goods_no_subsidy, stat_distr_no_subsidy, cons_fine_local_no_subsidy, a_prime_fine_local_no_subsidy,future_occupation_fine_local_no_subsidy,x_S_S_fine_no_subsidy,x_SC_fine_no_subsidy,x_BC_fine_no_subsidy, coeff_no_subsidy,
transaction_cost_loss_no_subsidy,nominal_GDP_no_subsidy,welfare_val_no_subsidy,Import_value_no_subsidy,Export_value_no_subsidy,current_worker_pop_no_subsidy,current_staple_pop_no_subsidy,current_cashcrop_pop_no_subsidy,
marketable_agr_surplus_share_no_subsidy,exportshare_cashcrop_no_subsidy,
fraction_model_no_subsidy,program_spending_no_subsidy,prod_value_improvement_no_subsidy,share_selling_increase_no_subsidy,exp_ratio_model_no_subsidy,mig_rate_model_no_subsidy,rural_pop_only_staples_model_no_subsidy,rural_pop_only_cashcrop_model_no_subsidy,
mean_land_share_to_staples_among_cc_model_no_subsidy,urban_rural_inc_ratio_model_no_subsidy,urban_rural_wealth_ratio_model_no_subsidy,urban_rural_consumption_ratio_model_no_subsidy,
p90_wealth_rural_no_subsidy,p90_wealth_urban_no_subsidy,p99_wealth_rural_no_subsidy,p99_wealth_urban_no_subsidy,p90_cons_tmp_no_subsidy,p90_income_tmp_no_subsidy,p99_cons_tmp_no_subsidy,p99_income_tmp_no_subsidy,
staple_productivity_no_subsidy,cashcrop_productivity_no_subsidy,manuf_productivity_no_subsidy,relative_land_to_staples_no_subsidy,relative_land_to_cashcrop_no_subsidy,share_constrained_cashcrop_no_subsidy,var_MPX_staples_S_no_subsidy,var_MPX_cashcrop_B_no_subsidy,var_MPX_cashcrop_S_no_subsidy,
share_constrained_staple_no_subsidy,APG_no_subsidy,urban_rural_consumption_ratio_model_real_no_subsidy,aggregate_consumption_no_subsidy,
worker_pop_effective_no_subsidy,prod_manuf_no_subsidy,total_entry_cost_no_subsidy,prod_staple_no_subsidy,prod_cashcrop_no_subsidy,input_staple_no_subsidy,input_cashcrop_no_subsidy,total_maintenance_cost_no_subsidy,
current_account_residual_no_subsidy,fraction_cashcrop_suboptimal_model_no_subsidy,V_saved_no_subsidy
,avg_labor_prod_rural_no_subsidy,avg_labor_prod_urban_no_subsidy,avg_agri_prod_rural_no_subsidy,avg_agri_prod_urban_no_subsidy
,var_MPX_cashcrop_no_subsidy,var_MPX_no_subsidy,TFP_no_subsidy,YL_manuf_no_subsidy , YL_agr_no_subsidy,coeff_var_labor_prod_rural_no_subsidy,coeff_var_labor_prod_urban_no_subsidy,coeff_var_agri_prod_rural_no_subsidy, coeff_var_agri_prod_urban_no_subsidy
,p90_wealth_no_subsidy,p99_wealth_no_subsidy,p90_cons_no_subsidy_rural,p99_cons_no_subsidy_rural,p90_cons_no_subsidy_urban,p99_cons_no_subsidy_urban,p90_income_no_subsidy_rural,
p99_income_no_subsidy_rural,p90_income_no_subsidy_urban,p99_income_no_subsidy_urban,wealth_of_workers_no_subsidy,wealth_of_staples_no_subsidy,wealth_of_cashcrop_no_subsidy
,c_B_worker_sum_no_subsidy,c_B_staple_sum_no_subsidy,c_B_cashcrop_sum_no_subsidy,c_S_worker_sum_no_subsidy,c_S_staple_sum_no_subsidy ,c_S_cashcrop_sum_no_subsidy,
transaction_cost_staple_sum_no_subsidy,transaction_cost_cashcrop_sum_no_subsidy,transaction_cost_worker_sum_no_subsidy,c_M_worker_sum_no_subsidy,c_M_staple_sum_no_subsidy,c_M_cashcrop_sum_no_subsidy,
MPX_mean_log_no_subsidy, MPX_mean_staples_S_log_no_subsidy,MPX_mean_cashcrop_log_no_subsidy
, APland_mean_log_no_subsidy,APland_mean_cashcrop_log_no_subsidy, APland_mean_staples_S_log_no_subsidy,var_APland_no_subsidy,var_APland_cashcrop_no_subsidy,var_APland_staples_S_no_subsidy,
c_S_W_fine_no_subsidy,c_B_W_fine_no_subsidy,c_M_W_fine_no_subsidy,c_S_S_fine_no_subsidy,c_B_S_fine_no_subsidy,c_M_S_fine_no_subsidy,c_S_B_fine_no_subsidy,c_B_B_fine_no_subsidy,c_M_B_fine_no_subsidy) = details_model(prices_no_subsidy,No_subsidy_parameter,2,moments,0.0,foreign_supply_capital_subsidy_b);

#No cbar QS + no subsidy equilibrium
prices_no_cbarQS_no_subsidy = [ 1.1965420095105022,
0.1658018507502125];
no_cbarQS_no_subsidy_parameter = copy(Baseline_parameter);
no_cbarQS_no_subsidy_parameter.τ_S=-0;
no_cbarQS_no_subsidy_parameter.Q_S=0;
no_cbarQS_no_subsidy_parameter.c̄_S=0.0;
(residual_goods_no_cbarQS_no_subsidy, stat_distr_no_cbarQS_no_subsidy, cons_fine_local_no_cbarQS_no_subsidy, a_prime_fine_local_no_cbarQS_no_subsidy,future_occupation_fine_local_no_cbarQS_no_subsidy,x_S_S_fine_no_cbarQS_no_subsidy,x_SC_fine_no_cbarQS_no_subsidy,x_BC_fine_no_cbarQS_no_subsidy, coeff_no_cbarQS_no_subsidy,
transaction_cost_loss_no_cbarQS_no_subsidy,nominal_GDP_no_cbarQS_no_subsidy,welfare_val_no_cbarQS_no_subsidy,Import_value_no_cbarQS_no_subsidy,Export_value_no_cbarQS_no_subsidy,current_worker_pop_no_cbarQS_no_subsidy,current_staple_pop_no_cbarQS_no_subsidy,current_cashcrop_pop_no_cbarQS_no_subsidy,
marketable_agr_surplus_share_no_cbarQS_no_subsidy,exportshare_cashcrop_no_cbarQS_no_subsidy,
fraction_model_no_cbarQS_no_subsidy,program_spending_no_cbarQS_no_subsidy,prod_value_improvement_no_cbarQS_no_subsidy,share_selling_increase_no_cbarQS_no_subsidy,exp_ratio_model_no_cbarQS_no_subsidy,mig_rate_model_no_cbarQS_no_subsidy,rural_pop_only_staples_model_no_cbarQS_no_subsidy,rural_pop_only_cashcrop_model_no_cbarQS_no_subsidy,
mean_land_share_to_staples_among_cc_model_no_cbarQS_no_subsidy,urban_rural_inc_ratio_model_no_cbarQS_no_subsidy,urban_rural_wealth_ratio_model_no_cbarQS_no_subsidy,urban_rural_consumption_ratio_model_no_cbarQS_no_subsidy,
p90_wealth_rural_no_cbarQS_no_subsidy,p90_wealth_urban_no_cbarQS_no_subsidy,p99_wealth_rural_no_cbarQS_no_subsidy,p99_wealth_urban_no_cbarQS_no_subsidy,p90_cons_tmp_no_cbarQS_no_subsidy,p90_income_tmp_no_cbarQS_no_subsidy,p99_cons_tmp_no_cbarQS_no_subsidy,p99_income_tmp_no_cbarQS_no_subsidy,
staple_productivity_no_cbarQS_no_subsidy,cashcrop_productivity_no_cbarQS_no_subsidy,manuf_productivity_no_cbarQS_no_subsidy,relative_land_to_staples_no_cbarQS_no_subsidy,relative_land_to_cashcrop_no_cbarQS_no_subsidy,share_constrained_cashcrop_no_cbarQS_no_subsidy,var_MPX_staples_S_no_cbarQS_no_subsidy,var_MPX_cashcrop_B_no_cbarQS_no_subsidy,var_MPX_cashcrop_S_no_cbarQS_no_subsidy,
share_constrained_staple_no_cbarQS_no_subsidy,APG_no_cbarQS_no_subsidy,urban_rural_consumption_ratio_model_real_no_cbarQS_no_subsidy,aggregate_consumption_no_cbarQS_no_subsidy,
worker_pop_effective_no_cbarQS_no_subsidy,prod_manuf_no_cbarQS_no_subsidy,total_entry_cost_no_cbarQS_no_subsidy,prod_staple_no_cbarQS_no_subsidy,prod_cashcrop_no_cbarQS_no_subsidy,input_staple_no_cbarQS_no_subsidy,input_cashcrop_no_cbarQS_no_subsidy,total_maintenance_cost_no_cbarQS_no_subsidy,
current_account_residual_no_cbarQS_no_subsidy,fraction_cashcrop_suboptimal_model_no_cbarQS_no_subsidy,V_saved_no_cbarQS_no_subsidy,avg_labor_prod_rural_no_cbarQS_no_subsidy,
avg_labor_prod_urban_no_cbarQS_no_subsidy,avg_agri_prod_rural_no_cbarQS_no_subsidy,avg_agri_prod_urban_no_cbarQS_no_subsidy,
var_MPX_cashcrop_no_cbarQS_no_subsidy,var_MPX_no_cbarQS_no_subsidy,TFP_no_cbarQS_no_subsidy,YL_manuf_no_cbarQS_no_subsidy , YL_agr_no_cbarQS_no_subsidy,coeff_var_labor_prod_rural_no_cbarQS_no_subsidy,coeff_var_labor_prod_urban_no_cbarQS_no_subsidy,coeff_var_agri_prod_rural_no_cbarQS_no_subsidy, coeff_var_agri_prod_urban_no_cbarQS_no_subsidy,
p90_wealth_no_cbarQS_no_subsidy,p99_wealth_no_cbarQS_no_subsidy,p90_cons_no_cbarQS_no_subsidy_rural,p99_cons_no_cbarQS_no_subsidy_rural,p90_cons_no_cbarQS_no_subsidy_urban,p99_cons_no_cbarQS_no_subsidy_urban,p90_income_no_cbarQS_no_subsidy_rural,
p99_income_no_cbarQS_no_subsidy_rural,p90_income_no_cbarQS_no_subsidy_urban,p99_income_no_cbarQS_no_subsidy_urban,wealth_of_workers_no_cbarQS_no_subsidy,wealth_of_staples_no_cbarQS_no_subsidy,wealth_of_cashcrop_no_cbarQS_no_subsidy,
c_B_worker_sum_no_cbarQS_no_subsidy,c_B_staple_sum_no_cbarQS_no_subsidy,c_B_cashcrop_sum_no_cbarQS_no_subsidy,c_S_worker_sum_no_cbarQS_no_subsidy,c_S_staple_sum_no_cbarQS_no_subsidy ,c_S_cashcrop_sum_no_cbarQS_no_subsidy,
transaction_cost_staple_sum_no_cbarQS_no_subsidy,transaction_cost_cashcrop_sum_no_cbarQS_no_subsidy,transaction_cost_worker_sum_no_cbarQS_no_subsidy,c_M_worker_sum_no_cbarQS_no_subsidy,c_M_staple_sum_no_cbarQS_no_subsidy,c_M_cashcrop_sum_no_cbarQS_no_subsidy,
MPX_mean_log_no_cbarQS_no_subsidy, MPX_mean_staples_S_log_no_cbarQS_no_subsidy,MPX_mean_cashcrop_log_no_cbarQS_no_subsidy
, APland_mean_log_no_cbarQS_no_subsidy,APland_mean_cashcrop_log_no_cbarQS_no_subsidy, APland_mean_staples_S_log_no_cbarQS_no_subsidy,var_APland_no_cbarQS_no_subsidy,var_APland_cashcrop_no_cbarQS_no_subsidy,var_APland_staples_S_no_cbarQS_no_subsidy,
c_S_W_fine_no_cbarQS_no_subsidy,c_B_W_fine_no_cbarQS_no_subsidy,c_M_W_fine_no_cbarQS_no_subsidy,c_S_S_fine_no_cbarQS_no_subsidy,c_B_S_fine_no_cbarQS_no_subsidy,c_M_S_fine_no_cbarQS_no_subsidy,c_S_B_fine_no_cbarQS_no_subsidy,c_B_B_fine_no_cbarQS_no_subsidy,c_M_B_fine_no_cbarQS_no_subsidy) = details_model(prices_no_cbarQS_no_subsidy,no_cbarQS_no_subsidy_parameter,2,moments,0.0,foreign_supply_capital_subsidy_b);



#No cbar + no subsidy equilibrium
prices_no_cbar_no_subsidy = [   1.205937339361828,
0.18357263199475019];
no_cbar_no_subsidy_parameter = copy(Baseline_parameter);
no_cbar_no_subsidy_parameter.τ_S=-0;
no_cbar_no_subsidy_parameter.c̄_S=0.0;

(residual_goods_no_cbar_no_subsidy, stat_distr_no_cbar_no_subsidy, cons_fine_local_no_cbar_no_subsidy, a_prime_fine_local_no_cbar_no_subsidy,future_occupation_fine_local_no_cbar_no_subsidy,x_S_S_fine_no_cbar_no_subsidy,x_SC_fine_no_cbar_no_subsidy,x_BC_fine_no_cbar_no_subsidy, coeff_no_cbar_no_subsidy,
transaction_cost_loss_no_cbar_no_subsidy,nominal_GDP_no_cbar_no_subsidy,welfare_val_no_cbar_no_subsidy,Import_value_no_cbar_no_subsidy,Export_value_no_cbar_no_subsidy,current_worker_pop_no_cbar_no_subsidy,current_staple_pop_no_cbar_no_subsidy,current_cashcrop_pop_no_cbar_no_subsidy,
marketable_agr_surplus_share_no_cbar_no_subsidy,exportshare_cashcrop_no_cbar_no_subsidy,
fraction_model_no_cbar_no_subsidy,program_spending_no_cbar_no_subsidy,prod_value_improvement_no_cbar_no_subsidy,share_selling_increase_no_cbar_no_subsidy,exp_ratio_model_no_cbar_no_subsidy,mig_rate_model_no_cbar_no_subsidy,rural_pop_only_staples_model_no_cbar_no_subsidy,rural_pop_only_cashcrop_model_no_cbar_no_subsidy,
mean_land_share_to_staples_among_cc_model_no_cbar_no_subsidy,urban_rural_inc_ratio_model_no_cbar_no_subsidy,urban_rural_wealth_ratio_model_no_cbar_no_subsidy,urban_rural_consumption_ratio_model_no_cbar_no_subsidy,
p90_wealth_rural_no_cbar_no_subsidy,p90_wealth_urban_no_cbar_no_subsidy,p99_wealth_rural_no_cbar_no_subsidy,p99_wealth_urban_no_cbar_no_subsidy,p90_cons_tmp_no_cbar_no_subsidy,p90_income_tmp_no_cbar_no_subsidy,p99_cons_tmp_no_cbar_no_subsidy,p99_income_tmp_no_cbar_no_subsidy,
staple_productivity_no_cbar_no_subsidy,cashcrop_productivity_no_cbar_no_subsidy,manuf_productivity_no_cbar_no_subsidy,relative_land_to_staples_no_cbar_no_subsidy,relative_land_to_cashcrop_no_cbar_no_subsidy,share_constrained_cashcrop_no_cbar_no_subsidy,var_MPX_staples_S_no_cbar_no_subsidy,var_MPX_cashcrop_B_no_cbar_no_subsidy,var_MPX_cashcrop_S_no_cbar_no_subsidy,
share_constrained_staple_no_cbar_no_subsidy,APG_no_cbar_no_subsidy,urban_rural_consumption_ratio_model_real_no_cbar_no_subsidy,aggregate_consumption_no_cbar_no_subsidy,
worker_pop_effective_no_cbar_no_subsidy,prod_manuf_no_cbar_no_subsidy,total_entry_cost_no_cbar_no_subsidy,prod_staple_no_cbar_no_subsidy,prod_cashcrop_no_cbar_no_subsidy,input_staple_no_cbar_no_subsidy,input_cashcrop_no_cbar_no_subsidy,total_maintenance_cost_no_cbar_no_subsidy,
current_account_residual_no_cbar_no_subsidy,fraction_cashcrop_suboptimal_model_no_cbar_no_subsidy,V_saved_no_cbar_no_subsidy,avg_labor_prod_rural_no_cbar_no_subsidy 
,avg_labor_prod_urban_no_cbar_no_subsidy ,avg_agri_prod_rural_no_cbar_no_subsidy ,avg_agri_prod_urban_no_cbar_no_subsidy
,var_MPX_cashcrop_no_cbar_no_subsidy,var_MPX_no_cbar_no_subsidy,TFP_no_cbar_no_subsidy,YL_manuf_no_cbar_no_subsidy , YL_agr_no_cbar_no_subsidy,coeff_var_labor_prod_rural_no_cbar_no_subsidy,coeff_var_labor_prod_urban_no_cbar_no_subsidy,coeff_var_agri_prod_rural_no_cbar_no_subsidy, coeff_var_agri_prod_urban_no_cbar_no_subsidy,
p90_wealth_no_cbar_no_subsidy,p99_wealth_no_cbar_no_subsidy,p90_cons_no_cbar_no_subsidy_rural,p99_cons_no_cbar_no_subsidy_rural,p90_cons_no_cbar_no_subsidy_urban,p99_cons_no_cbar_no_subsidy_urban,p90_income_no_cbar_no_subsidy_rural,
p99_income_no_cbar_no_subsidy_rural,p90_income_no_cbar_no_subsidy_urban,p99_income_no_cbar_no_subsidy_urban,wealth_of_workers_no_cbar_no_subsidy
,wealth_of_staples_no_cbar_no_subsidy,wealth_of_cashcrop_no_cbar_no_subsidy,
c_B_worker_sum_no_cbar_no_subsidy,c_B_staple_sum_no_cbar_no_subsidy,c_B_cashcrop_sum_no_cbar_no_subsidy,c_S_worker_sum_no_cbar_no_subsidy,c_S_staple_sum_no_cbar_no_subsidy ,c_S_cashcrop_sum_no_cbar_no_subsidy,
transaction_cost_staple_sum_no_cbar_no_subsidy,transaction_cost_cashcrop_sum_no_cbar_no_subsidy,transaction_cost_worker_sum_no_cbar_no_subsidy,c_M_worker_sum_no_cbar_no_subsidy,c_M_staple_sum_no_cbar_no_subsidy,c_M_cashcrop_sum_no_cbar_no_subsidy,
MPX_mean_log_no_cbar_no_subsidy, MPX_mean_staples_S_log_no_cbar_no_subsidy,MPX_mean_cashcrop_log_no_cbar_no_subsidy
, APland_mean_log_no_cbar_no_subsidy,APland_mean_cashcrop_log_no_cbar_no_subsidy, APland_mean_staples_S_log_no_cbar_no_subsidy,var_APland_no_cbar_no_subsidy,var_APland_cashcrop_no_cbar_no_subsidy,var_APland_staples_S_no_cbar_no_subsidy,
c_S_W_fine_no_cbar_no_subsidy,c_B_W_fine_no_cbar_no_subsidy,c_M_W_fine_no_cbar_no_subsidy,c_S_S_fine_no_cbar_no_subsidy,c_B_S_fine_no_cbar_no_subsidy,c_M_S_fine_no_cbar_no_subsidy,c_S_B_fine_no_cbar_no_subsidy,c_B_B_fine_no_cbar_no_subsidy,c_M_B_fine_no_cbar_no_subsidy) = details_model(prices_no_cbar_no_subsidy,no_cbar_no_subsidy_parameter,2,moments,0.0,foreign_supply_capital_subsidy_b);

#No QS + no subsidy equilibrium
prices_no_QS_no_subsidy = [  1.186049115293932,
0.1641263316802302]
no_QS_no_subsidy_parameter = copy(Baseline_parameter);
no_QS_no_subsidy_parameter.τ_S=-0;
no_QS_no_subsidy_parameter.Q_S=0;

(residual_goods_no_QS_no_subsidy, stat_distr_no_QS_no_subsidy, cons_fine_local_no_QS_no_subsidy, a_prime_fine_local_no_QS_no_subsidy,future_occupation_fine_local_no_QS_no_subsidy,x_S_S_fine_no_QS_no_subsidy,x_SC_fine_no_QS_no_subsidy,x_BC_fine_no_QS_no_subsidy, coeff_no_QS_no_subsidy,
transaction_cost_loss_no_QS_no_subsidy,nominal_GDP_no_QS_no_subsidy,welfare_val_no_QS_no_subsidy,Import_value_no_QS_no_subsidy,Export_value_no_QS_no_subsidy,current_worker_pop_no_QS_no_subsidy,current_staple_pop_no_QS_no_subsidy,current_cashcrop_pop_no_QS_no_subsidy,
marketable_agr_surplus_share_no_QS_no_subsidy,exportshare_cashcrop_no_QS_no_subsidy,
fraction_model_no_QS_no_subsidy,program_spending_no_QS_no_subsidy,prod_value_improvement_no_QS_no_subsidy,share_selling_increase_no_QS_no_subsidy,exp_ratio_model_no_QS_no_subsidy,mig_rate_model_no_QS_no_subsidy,rural_pop_only_staples_model_no_QS_no_subsidy,rural_pop_only_cashcrop_model_no_QS_no_subsidy,
mean_land_share_to_staples_among_cc_model_no_QS_no_subsidy,urban_rural_inc_ratio_model_no_QS_no_subsidy,urban_rural_wealth_ratio_model_no_QS_no_subsidy,urban_rural_consumption_ratio_model_no_QS_no_subsidy,
p90_wealth_rural_no_QS_no_subsidy,p90_wealth_urban_no_QS_no_subsidy,p99_wealth_rural_no_QS_no_subsidy,p99_wealth_urban_no_QS_no_subsidy,p90_cons_tmp_no_QS_no_subsidy,p90_income_tmp_no_QS_no_subsidy,p99_cons_tmp_no_QS_no_subsidy,p99_income_tmp_no_QS_no_subsidy,
staple_productivity_no_QS_no_subsidy,cashcrop_productivity_no_QS_no_subsidy,manuf_productivity_no_QS_no_subsidy,relative_land_to_staples_no_QS_no_subsidy,relative_land_to_cashcrop_no_QS_no_subsidy,share_constrained_cashcrop_no_QS_no_subsidy,var_MPX_staples_S_no_QS_no_subsidy,var_MPX_cashcrop_B_no_QS_no_subsidy,var_MPX_cashcrop_S_no_QS_no_subsidy,
share_constrained_staple_no_QS_no_subsidy,APG_no_QS_no_subsidy,urban_rural_consumption_ratio_model_real_no_QS_no_subsidy,aggregate_consumption_no_QS_no_subsidy,
worker_pop_effective_no_QS_no_subsidy,prod_manuf_no_QS_no_subsidy,total_entry_cost_no_QS_no_subsidy,prod_staple_no_QS_no_subsidy,prod_cashcrop_no_QS_no_subsidy,input_staple_no_QS_no_subsidy,input_cashcrop_no_QS_no_subsidy,total_maintenance_cost_no_QS_no_subsidy,
current_account_residual_no_QS_no_subsidy,fraction_cashcrop_suboptimal_model_no_QS_no_subsidy,V_saved_no_QS_no_subsidy
,avg_labor_prod_rural_no_QS_no_subsidy,avg_labor_prod_urban_no_QS_no_subsidy,avg_agri_prod_rural_no_QS_no_subsidy,avg_agri_prod_urban_no_QS_no_subsidy
,var_MPX_cashcrop_no_QS_no_subsidy,var_MPX_no_QS_no_subsidy,TFP_no_QS_no_subsidy,YL_manuf_no_QS_no_subsidy , YL_agr_no_QS_no_subsidy,coeff_var_labor_prod_rural_no_QS_no_subsidy,coeff_var_labor_prod_urban_no_QS_no_subsidy,coeff_var_agri_prod_rural_no_QS_no_subsidy, coeff_var_agri_prod_urban_no_QS_no_subsidy,
p90_wealth_no_QS_no_subsidy,p99_wealth_no_QS_no_subsidy,p90_cons_no_QS_no_subsidy_rural,p99_cons_no_QS_no_subsidy_rural,p90_cons_no_QS_no_subsidy_urban,p99_cons_no_QS_no_subsidy_urban,p90_income_no_QS_no_subsidy_rural,
p99_income_no_QS_no_subsidy_rural,p90_income_no_QS_no_subsidy_urban,p99_income_no_QS_no_subsidy_urban,
wealth_of_workers_no_QS_no_subsidy,wealth_of_staples_no_QS_no_subsidy,wealth_of_cashcrop_no_QS_no_subsidy,
c_B_worker_sum_no_QS_no_subsidy,c_B_staple_sum_no_QS_no_subsidy,c_B_cashcrop_sum_no_QS_no_subsidy,c_S_worker_sum_no_QS_no_subsidy,c_S_staple_sum_no_QS_no_subsidy ,c_S_cashcrop_sum_no_QS_no_subsidy,
transaction_cost_staple_sum_no_QS_no_subsidy,transaction_cost_cashcrop_sum_no_QS_no_subsidy,transaction_cost_worker_sum_no_QS_no_subsidy,c_M_worker_sum_no_QS_no_subsidy,c_M_staple_sum_no_QS_no_subsidy,c_M_cashcrop_sum_no_QS_no_subsidy,
MPX_mean_log_no_QS_no_subsidy, MPX_mean_staples_S_log_no_QS_no_subsidy,MPX_mean_cashcrop_log_no_QS_no_subsidy
, APland_mean_log_no_QS_no_subsidy,APland_mean_cashcrop_log_no_QS_no_subsidy, APland_mean_staples_S_log_no_QS_no_subsidy,var_APland_no_QS_no_subsidy,var_APland_cashcrop_no_QS_no_subsidy,var_APland_staples_S_no_QS_no_subsidy,
c_S_W_fine_no_QS_no_subsidy,c_B_W_fine_no_QS_no_subsidy,c_M_W_fine_no_QS_no_subsidy,c_S_S_fine_no_QS_no_subsidy,c_B_S_fine_no_QS_no_subsidy,c_M_S_fine_no_QS_no_subsidy,c_S_B_fine_no_QS_no_subsidy,c_B_B_fine_no_QS_no_subsidy,c_M_B_fine_no_QS_no_subsidy) = details_model(prices_no_QS_no_subsidy,no_QS_no_subsidy_parameter,2,moments,0.0,foreign_supply_capital_subsidy_b);



#    No subsidy equilibrium + No F_W
prices_no_F_W_no_subsidy = [ 1.1678794680674407,
0.09092347139175247];#!
no_F_W_no_subsidy_parameter = copy(Baseline_parameter);
no_F_W_no_subsidy_parameter.F_W= 0.0;
no_F_W_no_subsidy_parameter.τ_S=-0.0;
(residual_goods_no_F_W_no_subsidy, stat_distr_no_F_W_no_subsidy, cons_fine_local_no_F_W_no_subsidy, a_prime_fine_local_no_F_W_no_subsidy,future_occupation_fine_local_no_F_W_no_subsidy,x_S_S_fine_no_F_W_no_subsidy,x_SC_fine_no_F_W_no_subsidy,x_BC_fine_no_F_W_no_subsidy, coeff_no_F_W_no_subsidy,
transaction_cost_loss_no_F_W_no_subsidy,nominal_GDP_no_F_W_no_subsidy,welfare_val_no_F_W_no_subsidy,Import_value_no_F_W_no_subsidy,Export_value_no_F_W_no_subsidy,current_worker_pop_no_F_W_no_subsidy,current_staple_pop_no_F_W_no_subsidy,current_cashcrop_pop_no_F_W_no_subsidy,
marketable_agr_surplus_share_no_F_W_no_subsidy,exportshare_cashcrop_no_F_W_no_subsidy,
fraction_model_no_F_W_no_subsidy,program_spending_no_F_W_no_subsidy,prod_value_improvement_no_F_W_no_subsidy,share_selling_increase_no_F_W_no_subsidy,exp_ratio_model_no_F_W_no_subsidy,mig_rate_model_no_F_W_no_subsidy,rural_pop_only_staples_model_no_F_W_no_subsidy,rural_pop_only_cashcrop_model_no_F_W_no_subsidy,
mean_land_share_to_staples_among_cc_model_no_F_W_no_subsidy,urban_rural_inc_ratio_model_no_F_W_no_subsidy,urban_rural_wealth_ratio_model_no_F_W_no_subsidy,urban_rural_consumption_ratio_model_no_F_W_no_subsidy,
p90_wealth_rural_no_F_W_no_subsidy,p90_wealth_urban_no_F_W_no_subsidy,p99_wealth_rural_no_F_W_no_subsidy,p99_wealth_urban_no_F_W_no_subsidy,p90_cons_tmp_no_F_W_no_subsidy,p90_income_tmp_no_F_W_no_subsidy,p99_cons_tmp_no_F_W_no_subsidy,p99_income_tmp_no_F_W_no_subsidy,
staple_productivity_no_F_W_no_subsidy,cashcrop_productivity_no_F_W_no_subsidy,manuf_productivity_no_F_W_no_subsidy,relative_land_to_staples_no_F_W_no_subsidy,relative_land_to_cashcrop_no_F_W_no_subsidy,share_constrained_cashcrop_no_F_W_no_subsidy,var_MPX_staples_S_no_F_W_no_subsidy,var_MPX_cashcrop_B_no_F_W_no_subsidy,var_MPX_cashcrop_S_no_F_W_no_subsidy,
share_constrained_staple_no_F_W_no_subsidy,APG_no_F_W_no_subsidy,urban_rural_consumption_ratio_model_real_no_F_W_no_subsidy,aggregate_consumption_no_F_W_no_subsidy,
worker_pop_effective_no_F_W_no_subsidy,prod_manuf_no_F_W_no_subsidy,total_entry_cost_no_F_W_no_subsidy,prod_staple_no_F_W_no_subsidy,prod_cashcrop_no_F_W_no_subsidy,input_staple_no_F_W_no_subsidy,input_cashcrop_no_F_W_no_subsidy,
total_maintenance_cost_no_F_W_no_subsidy,current_account_residual_no_F_W_no_subsidy,fraction_cashcrop_suboptimal_model_no_F_W_no_subsidy,V_saved_no_F_W_no_subsidy
,avg_labor_prod_rural_no_F_W_no_subsidy,avg_labor_prod_urban_no_F_W_no_subsidy,avg_agri_prod_rural_no_F_W_no_subsidy,avg_agri_prod_urban_no_F_W_no_subsidy
,var_MPX_cashcrop_no_F_W_no_subsidy,var_MPX_no_F_W_no_subsidy,TFP_no_F_W_no_subsidy,YL_manuf_no_F_W_no_subsidy , YL_agr_no_F_W_no_subsidy,coeff_var_labor_prod_rural_no_F_W_no_subsidy,coeff_var_labor_prod_urban_no_F_W_no_subsidy,coeff_var_agri_prod_rural_no_F_W_no_subsidy, coeff_var_agri_prod_urban_no_F_W_no_subsidy,
p90_wealth_no_F_W_no_subsidy,p99_wealth_no_F_W_no_subsidy,p90_cons_no_F_W_no_subsidy_rural,p99_cons_no_F_W_no_subsidy_rural,p90_cons_no_F_W_no_subsidy_urban,p99_cons_no_F_W_no_subsidy_urban,p90_income_no_F_W_no_subsidy_rural,
p99_income_no_F_W_no_subsidy_rural,p90_income_no_F_W_no_subsidy_urban,p99_income_no_F_W_no_subsidy_urban,
wealth_of_workers_no_F_W_no_subsidy,wealth_of_staples_no_F_W_no_subsidy,wealth_of_cashcrop_no_F_W_no_subsidy,
c_B_worker_sum_no_F_W_no_subsidy,c_B_staple_sum_no_F_W_no_subsidy,c_B_cashcrop_sum_no_F_W_no_subsidy,c_S_worker_sum_no_F_W_no_subsidy,c_S_staple_sum_no_F_W_no_subsidy ,c_S_cashcrop_sum_no_F_W_no_subsidy,
transaction_cost_staple_sum_no_F_W_no_subsidy,transaction_cost_cashcrop_sum_no_F_W_no_subsidy,transaction_cost_worker_sum_no_F_W_no_subsidy,c_M_worker_sum_no_F_W_no_subsidy,c_M_staple_sum_no_F_W_no_subsidy,c_M_cashcrop_sum_no_F_W_no_subsidy,
MPX_mean_log_no_F_W_no_subsidy, MPX_mean_staples_S_log_no_F_W_no_subsidy,MPX_mean_cashcrop_log_no_F_W_no_subsidy
, APland_mean_log_no_F_W_no_subsidy,APland_mean_cashcrop_log_no_F_W_no_subsidy, APland_mean_staples_S_log_no_F_W_no_subsidy,var_APland_no_F_W_no_subsidy,var_APland_cashcrop_no_F_W_no_subsidy,var_APland_staples_S_no_F_W_no_subsidy,
c_S_W_fine_no_F_W_no_subsidy,c_B_W_fine_no_F_W_no_subsidy,c_M_W_fine_no_F_W_no_subsidy,c_S_S_fine_no_F_W_no_subsidy,c_B_S_fine_no_F_W_no_subsidy,c_M_S_fine_no_F_W_no_subsidy,c_S_B_fine_no_F_W_no_subsidy,c_B_B_fine_no_F_W_no_subsidy,c_M_B_fine_no_F_W_no_subsidy) = details_model(prices_no_F_W_no_subsidy,no_F_W_no_subsidy_parameter,2,moments,0.0,foreign_supply_capital_subsidy_b);

#    No subsidy equilibrium + No FM_B
prices_no_FM_B_no_subsidy = [ 1.168581278103919,
0.21559107857631152];
no_FM_B_no_subsidy_parameter = copy(Baseline_parameter);
no_FM_B_no_subsidy_parameter.FM_B= 0.0;
no_FM_B_no_subsidy_parameter.τ_S=-0.0;
(residual_goods_no_FM_B_no_subsidy, stat_distr_no_FM_B_no_subsidy, cons_fine_local_no_FM_B_no_subsidy, a_prime_fine_local_no_FM_B_no_subsidy,future_occupation_fine_local_no_FM_B_no_subsidy,x_S_S_fine_no_FM_B_no_subsidy,x_SC_fine_no_FM_B_no_subsidy,x_BC_fine_no_FM_B_no_subsidy, coeff_no_FM_B_no_subsidy,
transaction_cost_loss_no_FM_B_no_subsidy,nominal_GDP_no_FM_B_no_subsidy,welfare_val_no_FM_B_no_subsidy,Import_value_no_FM_B_no_subsidy,Export_value_no_FM_B_no_subsidy,current_worker_pop_no_FM_B_no_subsidy,current_staple_pop_no_FM_B_no_subsidy,current_cashcrop_pop_no_FM_B_no_subsidy,
marketable_agr_surplus_share_no_FM_B_no_subsidy,exportshare_cashcrop_no_FM_B_no_subsidy,
fraction_model_no_FM_B_no_subsidy,program_spending_no_FM_B_no_subsidy,prod_value_improvement_no_FM_B_no_subsidy,share_selling_increase_no_FM_B_no_subsidy,exp_ratio_model_no_FM_B_no_subsidy,mig_rate_model_no_FM_B_no_subsidy,rural_pop_only_staples_model_no_FM_B_no_subsidy,rural_pop_only_cashcrop_model_no_FM_B_no_subsidy,
mean_land_share_to_staples_among_cc_model_no_FM_B_no_subsidy,urban_rural_inc_ratio_model_no_FM_B_no_subsidy,urban_rural_wealth_ratio_model_no_FM_B_no_subsidy,urban_rural_consumption_ratio_model_no_FM_B_no_subsidy,
p90_wealth_rural_no_FM_B_no_subsidy,p90_wealth_urban_no_FM_B_no_subsidy,p99_wealth_rural_no_FM_B_no_subsidy,p99_wealth_urban_no_FM_B_no_subsidy,p90_cons_tmp_no_FM_B_no_subsidy,p90_income_tmp_no_FM_B_no_subsidy,p99_cons_tmp_no_FM_B_no_subsidy,p99_income_tmp_no_FM_B_no_subsidy,
staple_productivity_no_FM_B_no_subsidy,cashcrop_productivity_no_FM_B_no_subsidy,manuf_productivity_no_FM_B_no_subsidy,relative_land_to_staples_no_FM_B_no_subsidy,relative_land_to_cashcrop_no_FM_B_no_subsidy,share_constrained_cashcrop_no_FM_B_no_subsidy,var_MPX_staples_S_no_FM_B_no_subsidy,var_MPX_cashcrop_B_no_FM_B_no_subsidy,var_MPX_cashcrop_S_no_FM_B_no_subsidy,
share_constrained_staple_no_FM_B_no_subsidy,APG_no_FM_B_no_subsidy,urban_rural_consumption_ratio_model_real_no_FM_B_no_subsidy,aggregate_consumption_no_FM_B_no_subsidy,
worker_pop_effective_no_FM_B_no_subsidy,prod_manuf_no_FM_B_no_subsidy,total_entry_cost_no_FM_B_no_subsidy,prod_staple_no_FM_B_no_subsidy,prod_cashcrop_no_FM_B_no_subsidy,input_staple_no_FM_B_no_subsidy,input_cashcrop_no_FM_B_no_subsidy,
total_maintenance_cost_no_FM_B_no_subsidy,current_account_residual_no_FM_B_no_subsidy,fraction_cashcrop_suboptimal_model_no_FM_B_no_subsidy,V_saved_no_FM_B_no_subsidy,avg_labor_prod_rural_no_FM_B_no_subsidy
,avg_labor_prod_urban_no_FM_B_no_subsidy,avg_agri_prod_rural_no_FM_B_no_subsidy,avg_agri_prod_urban_no_FM_B_no_subsidy,
var_MPX_cashcrop_no_FM_B_no_subsidy,var_MPX_no_FM_B_no_subsidy,TFP_no_FM_B_no_subsidy,YL_manuf_no_FM_B_no_subsidy ,
 YL_agr_no_FM_B_no_subsidy,coeff_var_labor_prod_rural_no_FM_B_no_subsidy,coeff_var_labor_prod_urban_no_FM_B_no_subsidy,
coeff_var_agri_prod_rural_no_FM_B_no_subsidy, coeff_var_agri_prod_urban_no_FM_B_no_subsidy,
p90_wealth_no_FM_B_no_subsidy,p99_wealth_no_FM_B_no_subsidy,p90_cons_no_FM_B_no_subsidy_rural,p99_cons_no_FM_B_no_subsidy_rural,p90_cons_no_FM_B_no_subsidy_urban,p99_cons_no_FM_B_no_subsidy_urban,p90_income_no_FM_B_no_subsidy_rural,
p99_income_no_FM_B_no_subsidy_rural,p90_income_no_FM_B_no_subsidy_urban,p99_income_no_FM_B_no_subsidy_urban,wealth_of_workers_no_FM_B_no_subsidy ,
wealth_of_staples_no_FM_B_no_subsidy ,wealth_of_cashcrop_no_FM_B_no_subsidy ,
c_B_worker_sum_no_FM_B_no_subsidy,c_B_staple_sum_no_FM_B_no_subsidy,c_B_cashcrop_sum_no_FM_B_no_subsidy,c_S_worker_sum_no_FM_B_no_subsidy,c_S_staple_sum_no_FM_B_no_subsidy ,c_S_cashcrop_sum_no_FM_B_no_subsidy,
transaction_cost_staple_sum_no_FM_B_no_subsidy,transaction_cost_cashcrop_sum_no_FM_B_no_subsidy,transaction_cost_worker_sum_no_FM_B_no_subsidy,c_M_worker_sum_no_FM_B_no_subsidy,c_M_staple_sum_no_FM_B_no_subsidy,c_M_cashcrop_sum_no_FM_B_no_subsidy,
MPX_mean_log_no_FM_B_no_subsidy, MPX_mean_staples_S_log_no_FM_B_no_subsidy,MPX_mean_cashcrop_log_no_FM_B_no_subsidy
, APland_mean_log_no_FM_B_no_subsidy,APland_mean_cashcrop_log_no_FM_B_no_subsidy, APland_mean_staples_S_log_no_FM_B_no_subsidy,var_APland_no_FM_B_no_subsidy,var_APland_cashcrop_no_FM_B_no_subsidy,var_APland_staples_S_no_FM_B_no_subsidy,
c_S_W_fine_no_FM_B_no_subsidy,c_B_W_fine_no_FM_B_no_subsidy,c_M_W_fine_no_FM_B_no_subsidy,c_S_S_fine_no_FM_B_no_subsidy,c_B_S_fine_no_FM_B_no_subsidy,c_M_S_fine_no_FM_B_no_subsidy,c_S_B_fine_no_FM_B_no_subsidy,c_B_B_fine_no_FM_B_no_subsidy,c_M_B_fine_no_FM_B_no_subsidy) = details_model(prices_no_FM_B_no_subsidy,no_FM_B_no_subsidy_parameter,2,moments,0.0,foreign_supply_capital_subsidy_b);

#    No subsidy equilibrium + high kappa
prices_no_κ_no_subsidy = [ 1.1696929055108973,
0.1890191391243493]#!
no_κ_parameter_no_subsidy = copy(Baseline_parameter);
no_κ_parameter_no_subsidy.κ= 10.0;
no_κ_parameter_no_subsidy.τ_S=-0.0;
(residual_goods_no_κ_no_subsidy, stat_distr_no_κ_no_subsidy, cons_fine_local_no_κ_no_subsidy, a_prime_fine_local_no_κ_no_subsidy,future_occupation_fine_local_no_κ_no_subsidy,x_S_S_fine_no_κ_no_subsidy,x_SC_fine_no_κ_no_subsidy,x_BC_fine_no_κ_no_subsidy, coeff_no_κ_no_subsidy,
transaction_cost_loss_no_κ_no_subsidy,nominal_GDP_no_κ_no_subsidy,welfare_val_no_κ_no_subsidy,Import_value_no_κ_no_subsidy,Export_value_no_κ_no_subsidy,current_worker_pop_no_κ_no_subsidy,current_staple_pop_no_κ_no_subsidy,current_cashcrop_pop_no_κ_no_subsidy,
marketable_agr_surplus_share_no_κ_no_subsidy,exportshare_cashcrop_no_κ_no_subsidy,
fraction_model_no_κ_no_subsidy,program_spending_no_κ_no_subsidy,prod_value_improvement_no_κ_no_subsidy,share_selling_increase_no_κ_no_subsidy,exp_ratio_model_no_κ_no_subsidy,mig_rate_model_no_κ_no_subsidy,rural_pop_only_staples_model_no_κ_no_subsidy,rural_pop_only_cashcrop_model_no_κ_no_subsidy,
mean_land_share_to_staples_among_cc_model_no_κ_no_subsidy,urban_rural_inc_ratio_model_no_κ_no_subsidy,urban_rural_wealth_ratio_model_no_κ_no_subsidy,urban_rural_consumption_ratio_model_no_κ_no_subsidy,
p90_wealth_rural_no_κ_no_subsidy,p90_wealth_urban_no_κ_no_subsidy,p99_wealth_rural_no_κ_no_subsidy,p99_wealth_urban_no_κ_no_subsidy,p90_cons_tmp_no_κ_no_subsidy,p90_income_tmp_no_κ_no_subsidy,p99_cons_tmp_no_κ_no_subsidy,p99_income_tmp_no_κ_no_subsidy,
staple_productivity_no_κ_no_subsidy,cashcrop_productivity_no_κ_no_subsidy,manuf_productivity_no_κ_no_subsidy,relative_land_to_staples_no_κ_no_subsidy,relative_land_to_cashcrop_no_κ_no_subsidy,share_constrained_cashcrop_no_κ_no_subsidy,var_MPX_staples_S_no_κ_no_subsidy,var_MPX_cashcrop_B_no_κ_no_subsidy,var_MPX_cashcrop_S_no_κ_no_subsidy,
share_constrained_staple_no_κ_no_subsidy,APG_no_κ_no_subsidy,urban_rural_consumption_ratio_model_real_no_κ_no_subsidy,aggregate_consumption_no_κ_no_subsidy,
worker_pop_effective_no_κ_no_subsidy,prod_manuf_no_κ_no_subsidy,total_entry_cost_no_κ_no_subsidy,prod_staple_no_κ_no_subsidy,prod_cashcrop_no_κ_no_subsidy,input_staple_no_κ_no_subsidy,input_cashcrop_no_κ_no_subsidy,
total_maintenance_cost_no_κ_no_subsidy,current_account_residual_no_κ_no_subsidy,fraction_cashcrop_suboptimal_model_no_κ_no_subsidy,V_saved_no_κ_no_subsidy,avg_labor_prod_rural_no_κ_no_subsidy
,avg_labor_prod_urban_no_κ_no_subsidy,avg_agri_prod_rural_no_κ_no_subsidy,avg_agri_prod_urban_no_κ_no_subsidy,
var_MPX_cashcrop_no_κ_no_subsidy,var_MPX_no_κ_no_subsidy,TFP_no_κ_no_subsidy,YL_manuf_no_κ_no_subsidy , 
YL_agr_no_κ_no_subsidy,coeff_var_labor_prod_rural_no_κ_no_subsidy,coeff_var_labor_prod_urban_no_κ_no_subsidy,
coeff_var_agri_prod_rural_no_κ_no_subsidy, coeff_var_agri_prod_urban_no_κ_no_subsidy,
p90_wealth_no_κ_no_subsidy,p99_wealth_no_κ_no_subsidy,p90_cons_no_κ_no_subsidy_rural,p99_cons_no_κ_no_subsidy_rural,p90_cons_no_κ_no_subsidy_urban,p99_cons_no_κ_no_subsidy_urban,p90_income_no_κ_no_subsidy_rural,
p99_income_no_κ_no_subsidy_rural,p90_income_no_κ_no_subsidy_urban,p99_income_no_κ_no_subsidy_urban,wealth_of_workers_no_κ_no_subsidy,
wealth_of_staples_no_κ_no_subsidy,wealth_of_cashcrop_no_κ_no_subsidy,
c_B_worker_sum_no_κ_no_subsidy,c_B_staple_sum_no_κ_no_subsidy,c_B_cashcrop_sum_no_κ_no_subsidy,c_S_worker_sum_no_κ_no_subsidy,c_S_staple_sum_no_κ_no_subsidy ,c_S_cashcrop_sum_no_κ_no_subsidy,
transaction_cost_staple_sum_no_κ_no_subsidy,transaction_cost_cashcrop_sum_no_κ_no_subsidy,transaction_cost_worker_sum_no_κ_no_subsidy,c_M_worker_sum_no_κ_no_subsidy,c_M_staple_sum_no_κ_no_subsidy,c_M_cashcrop_sum_no_κ_no_subsidy,
MPX_mean_log_no_κ_no_subsidy, MPX_mean_staples_S_log_no_κ_no_subsidy,MPX_mean_cashcrop_log_no_κ_no_subsidy
, APland_mean_log_no_κ_no_subsidy,APland_mean_cashcrop_log_no_κ_no_subsidy, APland_mean_staples_S_log_no_κ_no_subsidy,var_APland_no_κ_no_subsidy,var_APland_cashcrop_no_κ_no_subsidy,var_APland_staples_S_no_κ_no_subsidy,
c_S_W_fine_no_κ_no_subsidy,c_B_W_fine_no_κ_no_subsidy,c_M_W_fine_no_κ_no_subsidy,c_S_S_fine_no_κ_no_subsidy,c_B_S_fine_no_κ_no_subsidy,c_M_S_fine_no_κ_no_subsidy,c_S_B_fine_no_κ_no_subsidy,c_B_B_fine_no_κ_no_subsidy,c_M_B_fine_no_κ_no_subsidy) = details_model(prices_no_κ_no_subsidy,no_FM_B_no_subsidy_parameter,2,moments,0.0,foreign_supply_capital_subsidy_b);


## Subsidy equilibria
#    Subsidy equilibrium without balanced budget
prices_subsidy_nb = [ 1.50872539725031,
0.25456229194429497];

(residual_goods_subsidy_nb, stat_distr_subsidy_nb, cons_fine_local_subsidy_nb, a_prime_fine_local_subsidy_nb,future_occupation_fine_local_subsidy_nb,x_S_S_fine_subsidy_nb,x_SC_fine_subsidy_nb,x_BC_fine_subsidy_nb, coeff_subsidy_nb,
transaction_cost_loss_subsidy_nb,nominal_GDP_subsidy_nb,welfare_val_subsidy_nb,Import_value_subsidy_nb,Export_value_subsidy_nb,current_worker_pop_subsidy_nb,current_staple_pop_subsidy_nb,current_cashcrop_pop_subsidy_nb,
marketable_agr_surplus_share_subsidy_nb,exportshare_cashcrop_subsidy_nb,
fraction_model_subsidy_nb,program_spending_subsidy_nb,prod_value_improvement_subsidy_nb,share_selling_increase_subsidy_nb,exp_ratio_model_subsidy_nb,mig_rate_model_subsidy_nb,rural_pop_only_staples_model_subsidy_nb,rural_pop_only_cashcrop_model_subsidy_nb,
mean_land_share_to_staples_among_cc_model_subsidy_nb,urban_rural_inc_ratio_model_subsidy_nb,urban_rural_wealth_ratio_model_subsidy_nb,urban_rural_consumption_ratio_model_subsidy_nb,
p90_wealth_rural_subsidy_nb,p90_wealth_urban_subsidy_nb,p99_wealth_rural_subsidy_nb,p99_wealth_urban_subsidy_nb,p90_cons_tmp_subsidy_nb,p90_income_tmp_subsidy_nb,p99_cons_tmp_subsidy_nb,p99_income_tmp_subsidy_nb,
staple_productivity_subsidy_nb,cashcrop_productivity_subsidy_nb,manuf_productivity_subsidy_nb,relative_land_to_staples_subsidy_nb,relative_land_to_cashcrop_subsidy_nb,share_constrained_cashcrop_subsidy_nb,var_MPX_staples_S_subsidy_nb,var_MPX_cashcrop_B_subsidy_nb,var_MPX_cashcrop_S_subsidy_nb,
share_constrained_staple_subsidy_nb,APG_subsidy_nb,urban_rural_consumption_ratio_model_real_subsidy_nb,aggregate_consumption_subsidy_nb,
worker_pop_effective_subsidy_nb,prod_manuf_subsidy_nb,total_entry_cost_subsidy_nb,prod_staple_subsidy_nb,prod_cashcrop_subsidy_nb,input_staple_subsidy_nb,input_cashcrop_subsidy_nb,
total_maintenance_cost_subsidy_nb,current_account_residual_subsidy_nb,fraction_cashcrop_suboptimal_model_subsidy_nb,V_saved_subsidy_nb
,avg_labor_prod_rural_subsidy_nb,avg_labor_prod_urban_subsidy_nb,avg_agri_prod_rural_subsidy_nb,avg_agri_prod_urban_subsidy_nb,
var_MPX_cashcrop_subsidy_nb,var_MPX_subsidy_nb,TFP_subsidy_nb,YL_manuf_subsidy_nb , 
YL_agr_subsidy_nb,coeff_var_labor_prod_rural_subsidy_nb,coeff_var_labor_prod_urban_subsidy_nb,
coeff_var_agri_prod_rural_subsidy_nb, coeff_var_agri_prod_urban_subsidy_nb,
p90_wealth_subsidy_nb,p99_wealth_subsidy_nb,p90_cons_subsidy_nb_rural,p99_cons_subsidy_nb_rural,p90_cons_subsidy_nb_urban,p99_cons_subsidy_nb_urban,p90_income_subsidy_nb_rural,
p99_income_subsidy_nb_rural,p90_income_subsidy_nb_urban,p99_income_subsidy_nb_urban,wealth_of_workers_subsidy_nb,wealth_of_staples_subsidy_nb,
wealth_of_cashcrop_subsidy_nb,
c_B_worker_sum_subsidy_nb,c_B_staple_sum_subsidy_nb,c_B_cashcrop_sum_subsidy_nb,c_S_worker_sum_subsidy_nb,c_S_staple_sum_subsidy_nb ,c_S_cashcrop_sum_subsidy_nb,
transaction_cost_staple_sum_subsidy_nb,transaction_cost_cashcrop_sum_subsidy_nb,transaction_cost_worker_sum_subsidy_nb,c_M_worker_sum_subsidy_nb,c_M_staple_sum_subsidy_nb,c_M_cashcrop_sum_subsidy_nb
,MPX_mean_log_subsidy_nb, MPX_mean_staples_S_log_subsidy_nb,MPX_mean_cashcrop_log_subsidy_nb
, APland_mean_log_subsidy_nb,APland_mean_cashcrop_log_subsidy_nb, APland_mean_staples_S_log_subsidy_nb,var_APland_subsidy_nb,var_APland_cashcrop_subsidy_nb,var_APland_staples_S_subsidy_nb,
c_S_W_fine_subsidy_nb,c_B_W_fine_subsidy_nb,c_M_W_fine_subsidy_nb,c_S_S_fine_subsidy_nb,c_B_S_fine_subsidy_nb,c_M_S_fine_subsidy_nb,c_S_B_fine_subsidy_nb,c_B_B_fine_subsidy_nb,c_M_B_fine_subsidy_nb) = details_model(prices_subsidy_nb,Baseline_parameter,2,moments,0.0,foreign_supply_capital_subsidy_b);



#    Subsidy equilibrium with balanced budget
prices_subsidy_b = [ 1.518014329164464,
0.2726473607850526,
0.1886863641824446];

(residual_goods_subsidy_b, stat_distr_subsidy_b, cons_fine_local_subsidy_b, a_prime_fine_local_subsidy_b,future_occupation_fine_local_subsidy_b,x_S_S_fine_subsidy_b,x_SC_fine_subsidy_b,x_BC_fine_subsidy_b, coeff_subsidy_b,
transaction_cost_loss_subsidy_b,nominal_GDP_subsidy_b,welfare_val_subsidy_b,Import_value_subsidy_b,Export_value_subsidy_b,current_worker_pop_subsidy_b,current_staple_pop_subsidy_b,current_cashcrop_pop_subsidy_b,
marketable_agr_surplus_share_subsidy_b,exportshare_cashcrop_subsidy_b,
fraction_model_subsidy_b,program_spending_subsidy_b,prod_value_improvement_subsidy_b,share_selling_increase_subsidy_b,exp_ratio_model_subsidy_b,mig_rate_model_subsidy_b,rural_pop_only_staples_model_subsidy_b,rural_pop_only_cashcrop_model_subsidy_b,
mean_land_share_to_staples_among_cc_model_subsidy_b,urban_rural_inc_ratio_model_subsidy_b,urban_rural_wealth_ratio_model_subsidy_b,urban_rural_consumption_ratio_model_subsidy_b,
p90_wealth_rural_subsidy_b,p90_wealth_urban_subsidy_b,p99_wealth_rural_subsidy_b,p99_wealth_urban_subsidy_b,p90_cons_tmp_subsidy_b,p90_income_tmp_subsidy_b,p99_cons_tmp_subsidy_b,p99_income_tmp_subsidy_b,
staple_productivity_subsidy_b,cashcrop_productivity_subsidy_b,manuf_productivity_subsidy_b,relative_land_to_staples_subsidy_b,relative_land_to_cashcrop_subsidy_b,share_constrained_cashcrop_subsidy_b,var_MPX_staples_S_subsidy_b,var_MPX_cashcrop_B_subsidy_b,var_MPX_cashcrop_S_subsidy_b,
share_constrained_staple_subsidy_b,APG_subsidy_b,urban_rural_consumption_ratio_model_real_subsidy_b,aggregate_consumption_subsidy_b,
worker_pop_effective_subsidy_b,prod_manuf_subsidy_b,total_entry_cost_subsidy_b,prod_staple_subsidy_b,prod_cashcrop_subsidy_b,input_staple_subsidy_b,input_cashcrop_subsidy_b,
total_maintenance_cost_subsidy_b,current_account_residual_subsidy_b,fraction_cashcrop_suboptimal_model_subsidy_b,V_saved_subsidy_b
,avg_labor_prod_rural_subsidy_b,avg_labor_prod_urban_subsidy_b,avg_agri_prod_rural_subsidy_b,avg_agri_prod_urban_subsidy_b,
var_MPX_cashcrop_subsidy_b,var_MPX_subsidy_b,TFP_subsidy_b,YL_manuf_subsidy_b , 
YL_agr_subsidy_b,coeff_var_labor_prod_rural_subsidy_b,coeff_var_labor_prod_urban_subsidy_b,
coeff_var_agri_prod_rural_subsidy_b, coeff_var_agri_prod_urban_subsidy_b,
p90_wealth_subsidy_b,p99_wealth_subsidy_b,p90_cons_subsidy_b_rural,p99_cons_subsidy_b_rural,p90_cons_subsidy_b_urban,p99_cons_subsidy_b_urban,p90_income_subsidy_b_rural,
p99_income_subsidy_b_rural,p90_income_subsidy_b_urban,p99_income_tmp_urban,wealth_of_workers_subsidy_b,wealth_of_staples_subsidy_b,wealth_of_cashcrop_subsidy_b,
c_B_worker_sum_subsidy_b,c_B_staple_sum_subsidy_b,c_B_cashcrop_sum_subsidy_b,c_S_worker_sum_subsidy_b,c_S_staple_sum_subsidy_b ,c_S_cashcrop_sum_subsidy_b,
transaction_cost_staple_sum_subsidy_b,transaction_cost_cashcrop_sum_subsidy_b,transaction_cost_worker_sum_subsidy_b,c_M_worker_sum_subsidy_b,c_M_staple_sum_subsidy_b,c_M_cashcrop_sum_subsidy_b
,MPX_mean_log_subsidy_b, MPX_mean_staples_S_log_subsidy_b,MPX_mean_cashcrop_log_subsidy_b
, APland_mean_log_subsidy_b,APland_mean_cashcrop_log_subsidy_b, APland_mean_staples_S_log_subsidy_b,var_APland_subsidy_b,var_APland_cashcrop_subsidy_b,var_APland_staples_S_subsidy_b,
c_S_W_fine_subsidy_b,c_B_W_fine_subsidy_b,c_M_W_fine_subsidy_b,c_S_S_fine_subsidy_b,c_B_S_fine_subsidy_b,c_M_S_fine_subsidy_b,c_S_B_fine_subsidy_b,c_B_B_fine_subsidy_b,c_M_B_fine_subsidy_b)   = details_model(prices_subsidy_b,Baseline_parameter,2,moments,1.0,foreign_supply_capital_subsidy_b);

#    Subsidy partial equilibrium without any price adjustment
prices_subsidy_partial = [  1.188283177632008,
0.18349666252335484,
0.0];

(residual_goods_subsidy_partial, stat_distr_subsidy_partial, cons_fine_local_subsidy_partial, a_prime_fine_local_subsidy_partial,future_occupation_fine_local_subsidy_partial,x_S_S_fine_subsidy_partial,x_SC_fine_subsidy_partial,x_BC_fine_subsidy_partial, coeff_subsidy_partial,
transaction_cost_loss_subsidy_partial,nominal_GDP_subsidy_partial,welfare_val_subsidy_partial,Import_value_subsidy_partial,Export_value_subsidy_partial,current_worker_pop_subsidy_partial,current_staple_pop_subsidy_partial,current_cashcrop_pop_subsidy_partial,
marketable_agr_surplus_share_subsidy_partial,exportshare_cashcrop_subsidy_partial,
fraction_model_subsidy_partial,program_spending_subsidy_partial,prod_value_improvement_subsidy_partial,share_selling_increase_subsidy_partial,exp_ratio_model_subsidy_partial,mig_rate_model_subsidy_partial,rural_pop_only_staples_model_subsidy_partial,rural_pop_only_cashcrop_model_subsidy_partial,
mean_land_share_to_staples_among_cc_model_subsidy_partial,urban_rural_inc_ratio_model_subsidy_partial,urban_rural_wealth_ratio_model_subsidy_partial,urban_rural_consumption_ratio_model_subsidy_partial,
p90_wealth_rural_subsidy_partial,p90_wealth_urban_subsidy_partial,p99_wealth_rural_subsidy_partial,p99_wealth_urban_subsidy_partial,p90_cons_tmp_subsidy_partial,p90_income_tmp_subsidy_partial,p99_cons_tmp_subsidy_partial,p99_income_tmp_subsidy_partial,
staple_productivity_subsidy_partial,cashcrop_productivity_subsidy_partial,manuf_productivity_subsidy_partial,relative_land_to_staples_subsidy_partial,relative_land_to_cashcrop_subsidy_partial,share_constrained_cashcrop_subsidy_partial,var_MPX_staples_S_subsidy_partial,var_MPX_cashcrop_B_subsidy_partial,var_MPX_cashcrop_S_subsidy_partial,
share_constrained_staple_subsidy_partial,APG_subsidy_partial,urban_rural_consumption_ratio_model_real_subsidy_partial,aggregate_consumption_subsidy_partial,
worker_pop_effective_subsidy_partial,prod_manuf_subsidy_partial,total_entry_cost_subsidy_partial,prod_staple_subsidy_partial,prod_cashcrop_subsidy_partial,input_staple_subsidy_partial,input_cashcrop_subsidy_partial,
total_maintenance_cost_subsidy_partial,current_account_residual_subsidy_partial,fraction_cashcrop_suboptimal_model_subsidy_partial,V_saved_subsidy_partial
,avg_labor_prod_rural_subsidy_partial,avg_labor_prod_urban_subsidy_partial,avg_agri_prod_rural_subsidy_partial,avg_agri_prod_urban_subsidy_partial,
var_MPX_cashcrop_subsidy_partial,var_MPX_subsidy_partial,TFP_subsidy_partial,YL_manuf_subsidy_partial , 
YL_agr_subsidy_partial,coeff_var_labor_prod_rural_subsidy_partial,coeff_var_labor_prod_urban_subsidy_partial,
coeff_var_agri_prod_rural_subsidy_partial, coeff_var_agri_prod_urban_subsidy_partial,
p90_wealth_subsidy_partial,p99_wealth_subsidy_partial,p90_cons_subsidy_partial_rural,p99_cons_subsidy_partial_rural,p90_cons_subsidy_partial_urban,p99_cons_subsidy_partial_urban,p90_income_subsidy_partial_rural,
p99_income_subsidy_partial_rural,p90_income_subsidy_partial_urban,p99_income_tmp_urban,wealth_of_workers_subsidy_partial,wealth_of_staples_subsidy_partial,wealth_of_cashcrop_subsidy_partial,
c_B_worker_sum_subsidy_partial,c_B_staple_sum_subsidy_partial,c_B_cashcrop_sum_subsidy_partial,c_S_worker_sum_subsidy_partial,c_S_staple_sum_subsidy_partial ,c_S_cashcrop_sum_subsidy_partial,
transaction_cost_staple_sum_subsidy_partial,transaction_cost_cashcrop_sum_subsidy_partial,transaction_cost_worker_sum_subsidy_partial,c_M_worker_sum_subsidy_partial,c_M_staple_sum_subsidy_partial,c_M_cashcrop_sum_subsidy_partial
,MPX_mean_log_subsidy_partial, MPX_mean_staples_S_log_subsidy_partial,MPX_mean_cashcrop_log_subsidy_partial
, APland_mean_log_subsidy_partial,APland_mean_cashcrop_log_subsidy_partial, APland_mean_staples_S_log_subsidy_partial,var_APland_subsidy_partial,var_APland_cashcrop_subsidy_partial,var_APland_staples_S_subsidy_partial,
c_S_W_fine_subsidy_partial,c_B_W_fine_subsidy_partial,c_M_W_fine_subsidy_partial,c_S_S_fine_subsidy_partial,c_B_S_fine_subsidy_partial,c_M_S_fine_subsidy_partial,c_S_B_fine_subsidy_partial,c_B_B_fine_subsidy_partial,c_M_B_fine_subsidy_partial)    = details_model(prices_subsidy_partial,Baseline_parameter,2,moments,1.0,foreign_supply_capital_subsidy_b);


#    Subsidy equilibrium with balanced budget + No QS
prices_no_QS_subsidy_b =  [ 1.5073211207397117,
0.206306989703047,
0.17586675729012607]
no_QS_parameter = copy(Baseline_parameter);
no_QS_parameter.Q_S=0.0;
(residual_goods_no_QS_subsidy_b, stat_distr_no_QS_subsidy_b, cons_fine_local_no_QS_subsidy_b, a_prime_fine_local_no_QS_subsidy_b,future_occupation_fine_local_no_QS_subsidy_b,x_S_S_fine_no_QS_subsidy_b,x_SC_fine_no_QS_subsidy_b,x_BC_fine_no_QS_subsidy_b, coeff_no_QS_subsidy_b,
transaction_cost_loss_no_QS_subsidy_b,nominal_GDP_no_QS_subsidy_b,welfare_val_no_QS_subsidy_b,Import_value_no_QS_subsidy_b,Export_value_no_QS_subsidy_b,current_worker_pop_no_QS_subsidy_b,current_staple_pop_no_QS_subsidy_b,current_cashcrop_pop_no_QS_subsidy_b,
marketable_agr_surplus_share_no_QS_subsidy_b,exportshare_cashcrop_no_QS_subsidy_b,
fraction_model_no_QS_subsidy_b,program_spending_no_QS_subsidy_b,prod_value_improvement_no_QS_subsidy_b,share_selling_increase_no_QS_subsidy_b,exp_ratio_model_no_QS_subsidy_b,mig_rate_model_no_QS_subsidy_b,rural_pop_only_staples_model_no_QS_subsidy_b,rural_pop_only_cashcrop_model_no_QS_subsidy_b,
mean_land_share_to_staples_among_cc_model_no_QS_subsidy_b,urban_rural_inc_ratio_model_no_QS_subsidy_b,urban_rural_wealth_ratio_model_no_QS_subsidy_b,urban_rural_consumption_ratio_model_no_QS_subsidy_b,
p90_wealth_rural_no_QS_subsidy_b,p90_wealth_urban_no_QS_subsidy_b,p99_wealth_rural_no_QS_subsidy_b,p99_wealth_urban_no_QS_subsidy_b,p90_cons_tmp_no_QS_subsidy_b,p90_income_tmp_no_QS_subsidy_b,p99_cons_tmp_no_QS_subsidy_b,p99_income_tmp_no_QS_subsidy_b,
staple_productivity_no_QS_subsidy_b,cashcrop_productivity_no_QS_subsidy_b,manuf_productivity_no_QS_subsidy_b,relative_land_to_staples_no_QS_subsidy_b,relative_land_to_cashcrop_no_QS_subsidy_b,share_constrained_cashcrop_no_QS_subsidy_b,var_MPX_staples_S_no_QS_subsidy_b,var_MPX_cashcrop_B_no_QS_subsidy_b,var_MPX_cashcrop_S_no_QS_subsidy_b,
share_constrained_staple_no_QS_subsidy_b,APG_no_QS_subsidy_b,urban_rural_consumption_ratio_model_real_no_QS_subsidy_b,aggregate_consumption_no_QS_subsidy_b,
worker_pop_effective_no_QS_subsidy_b,prod_manuf_no_QS_subsidy_b,total_entry_cost_no_QS_subsidy_b,prod_staple_no_QS_subsidy_b,prod_cashcrop_no_QS_subsidy_b,input_staple_no_QS_subsidy_b,input_cashcrop_no_QS_subsidy_b,
total_maintenance_cost_no_QS_subsidy_b,current_account_residual_no_QS_subsidy_b,fraction_cashcrop_suboptimal_model_no_QS_subsidy_b,V_saved_no_QS_subsidy_b
,avg_labor_prod_rural_no_QS_subsidy_b,avg_labor_prod_urban_no_QS_subsidy_b,avg_agri_prod_rural_no_QS_subsidy_b,avg_agri_prod_urban_no_QS_subsidy_b,
var_MPX_cashcrop_no_QS_subsidy_b ,var_MPX_no_QS_subsidy_b ,TFP_no_QS_subsidy_b ,YL_manuf_no_QS_subsidy_b  , 
YL_agr_no_QS_subsidy_b ,coeff_var_labor_prod_rural_no_QS_subsidy_b ,coeff_var_labor_prod_urban_no_QS_subsidy_b ,
coeff_var_agri_prod_rural_no_QS_subsidy_b , coeff_var_agri_prod_urban_no_QS_subsidy_b ,
p90_wealth_no_QS_subsidy_b,p99_wealth_no_QS_subsidy_b,p90_cons_no_QS_subsidy_b_rural,p99_cons_no_QS_subsidy_b_rural,p90_cons_no_QS_subsidy_b_urban,p99_cons_no_QS_subsidy_b_urban,p90_income_no_QS_subsidy_b_rural,
p99_income_no_QS_subsidy_b_rural,p90_income_no_QS_subsidy_b_urban,p99_income_no_QS_subsidy_b_urban,wealth_of_workers_no_QS_subsidy_b,
wealth_of_staples_no_QS_subsidy_b,wealth_of_cashcrop_no_QS_subsidy_b,
c_B_worker_sum_no_QS_subsidy_b,c_B_staple_sum_no_QS_subsidy_b,c_B_cashcrop_sum_no_QS_subsidy_b,c_S_worker_sum_no_QS_subsidy_b,c_S_staple_sum_no_QS_subsidy_b ,c_S_cashcrop_sum_no_QS_subsidy_b,
transaction_cost_staple_sum_no_QS_subsidy_b,transaction_cost_cashcrop_sum_no_QS_subsidy_b,transaction_cost_worker_sum_no_QS_subsidy_b,c_M_worker_sum_no_QS_subsidy_b,c_M_staple_sum_no_QS_subsidy_b,c_M_cashcrop_sum_no_QS_subsidy_b
,MPX_mean_log_no_QS_subsidy_b, MPX_mean_staples_S_log_no_QS_subsidy_b,MPX_mean_cashcrop_log_no_QS_subsidy_b
, APland_mean_log_no_QS_subsidy_b,APland_mean_cashcrop_log_no_QS_subsidy_b, APland_mean_staples_S_log_no_QS_subsidy_b,var_APland_no_QS_subsidy_b,var_APland_cashcrop_no_QS_subsidy_b,var_APland_staples_S_no_QS_subsidy_b,
c_S_W_fine_no_QS_subsidy_b,c_B_W_fine_no_QS_subsidy_b,c_M_W_fine_no_QS_subsidy_b,c_S_S_fine_no_QS_subsidy_b,c_B_S_fine_no_QS_subsidy_b,c_M_S_fine_no_QS_subsidy_b,c_S_B_fine_no_QS_subsidy_b,c_B_B_fine_no_QS_subsidy_b,c_M_B_fine_no_QS_subsidy_b) = details_model(prices_no_QS_subsidy_b,no_QS_parameter,2,moments,1.0,foreign_supply_capital_subsidy_b);

#    Subsidy equilibrium with balanced budget + No cbar
prices_no_cbar_subsidy_b = [  1.5435433787380388,
0.259856531535514,
0.19293976562040993]
no_cbar_parameter = copy(Baseline_parameter);
no_cbar_parameter.c̄_S = 0.0
(residual_goods_no_cbar_subsidy_b, stat_distr_no_cbar_subsidy_b, cons_fine_local_no_cbar_subsidy_b, a_prime_fine_local_no_cbar_subsidy_b,future_occupation_fine_local_no_cbar_subsidy_b,x_S_S_fine_no_cbar_subsidy_b,x_SC_fine_no_cbar_subsidy_b,x_BC_fine_no_cbar_subsidy_b, coeff_no_cbar_subsidy_b,
transaction_cost_loss_no_cbar_subsidy_b,nominal_GDP_no_cbar_subsidy_b,welfare_val_no_cbar_subsidy_b,Import_value_no_cbar_subsidy_b,Export_value_no_cbar_subsidy_b,current_worker_pop_no_cbar_subsidy_b,current_staple_pop_no_cbar_subsidy_b,current_cashcrop_pop_no_cbar_subsidy_b,
marketable_agr_surplus_share_no_cbar_subsidy_b,exportshare_cashcrop_no_cbar_subsidy_b,
fraction_model_no_cbar_subsidy_b,program_spending_no_cbar_subsidy_b,prod_value_improvement_no_cbar_subsidy_b,share_selling_increase_no_cbar_subsidy_b,exp_ratio_model_no_cbar_subsidy_b,mig_rate_model_no_cbar_subsidy_b,rural_pop_only_staples_model_no_cbar_subsidy_b,rural_pop_only_cashcrop_model_no_cbar_subsidy_b,
mean_land_share_to_staples_among_cc_model_no_cbar_subsidy_b,urban_rural_inc_ratio_model_no_cbar_subsidy_b,urban_rural_wealth_ratio_model_no_cbar_subsidy_b,urban_rural_consumption_ratio_model_no_cbar_subsidy_b,
p90_wealth_rural_no_cbar_subsidy_b,p90_wealth_urban_no_cbar_subsidy_b,p99_wealth_rural_no_cbar_subsidy_b,p99_wealth_urban_no_cbar_subsidy_b,p90_cons_tmp_no_cbar_subsidy_b,p90_income_tmp_no_cbar_subsidy_b,p99_cons_tmp_no_cbar_subsidy_b,p99_income_tmp_no_cbar_subsidy_b,
staple_productivity_no_cbar_subsidy_b,cashcrop_productivity_no_cbar_subsidy_b,manuf_productivity_no_cbar_subsidy_b,relative_land_to_staples_no_cbar_subsidy_b,relative_land_to_cashcrop_no_cbar_subsidy_b,share_constrained_cashcrop_no_cbar_subsidy_b,var_MPX_staples_S_no_cbar_subsidy_b,var_MPX_cashcrop_B_no_cbar_subsidy_b,var_MPX_cashcrop_S_no_cbar_subsidy_b,
share_constrained_staple_no_cbar_subsidy_b,APG_no_cbar_subsidy_b,urban_rural_consumption_ratio_model_real_no_cbar_subsidy_b,aggregate_consumption_no_cbar_subsidy_b,
worker_pop_effective_no_cbar_subsidy_b,prod_manuf_no_cbar_subsidy_b,total_entry_cost_no_cbar_subsidy_b,prod_staple_no_cbar_subsidy_b,prod_cashcrop_no_cbar_subsidy_b,input_staple_no_cbar_subsidy_b,input_cashcrop_no_cbar_subsidy_b,
total_maintenance_cost_no_cbar_subsidy_b,current_account_residual_no_cbar_subsidy_b,fraction_cashcrop_suboptimal_model_no_cbar_subsidy_b,V_saved_no_cbar_subsidy_b
,avg_labor_prod_rural_no_cbar_subsidy_b,avg_labor_prod_urban_no_cbar_subsidy_b,avg_agri_prod_rural_no_cbar_subsidy_b,avg_agri_prod_urban_no_cbar_subsidy_b,
var_MPX_cashcrop_no_cbar_subsidy_b ,var_MPX_no_cbar_subsidy_b ,TFP_no_cbar_subsidy_b ,YL_manuf_no_cbar_subsidy_b  , 
YL_agr_no_cbar_subsidy_b ,coeff_var_labor_prod_rural_no_cbar_subsidy_b ,coeff_var_labor_prod_urban_no_cbar_subsidy_b ,
coeff_var_agri_prod_rural_no_cbar_subsidy_b , coeff_var_agri_prod_urban_no_cbar_subsidy_b ,
p90_wealth_no_cbar_subsidy_b,p99_wealth_no_cbar_subsidy_b,p90_cons_no_cbar_subsidy_b_rural,p99_cons_no_cbar_subsidy_b_rural,p90_cons_no_cbar_subsidy_b_urban,p99_cons_no_cbar_subsidy_b_urban,p90_income_no_cbar_subsidy_b_rural,
p99_income_no_cbar_subsidy_b_rural,p90_income_no_cbar_subsidy_b_urban,p99_income_no_cbar_subsidy_b_urban,wealth_of_workers_no_cbar_subsidy_b,
wealth_of_staples_no_cbar_subsidy_b,wealth_of_cashcrop_no_cbar_subsidy_b,
c_B_worker_sum_no_cbar_subsidy_b,c_B_staple_sum_no_cbar_subsidy_b,c_B_cashcrop_sum_no_cbar_subsidy_b,c_S_worker_sum_no_cbar_subsidy_b,c_S_staple_sum_no_cbar_subsidy_b ,c_S_cashcrop_sum_no_cbar_subsidy_b,
transaction_cost_staple_sum_no_cbar_subsidy_b,transaction_cost_cashcrop_sum_no_cbar_subsidy_b,transaction_cost_worker_sum_no_cbar_subsidy_b,c_M_worker_sum_no_cbar_subsidy_b,c_M_staple_sum_no_cbar_subsidy_b,c_M_cashcrop_sum_no_cbar_subsidy_b
,MPX_mean_log_no_cbar_subsidy_b, MPX_mean_staples_S_log_no_cbar_subsidy_b,MPX_mean_cashcrop_log_no_cbar_subsidy_b
, APland_mean_log_no_cbar_subsidy_b,APland_mean_cashcrop_log_no_cbar_subsidy_b, APland_mean_staples_S_log_no_cbar_subsidy_b,var_APland_no_cbar_subsidy_b,var_APland_cashcrop_no_cbar_subsidy_b,var_APland_staples_S_no_cbar_subsidy_b,
c_S_W_fine_no_cbar_subsidy_b,c_B_W_fine_no_cbar_subsidy_b,c_M_W_fine_no_cbar_subsidy_b,c_S_S_fine_no_cbar_subsidy_b,c_B_S_fine_no_cbar_subsidy_b,c_M_S_fine_no_cbar_subsidy_b,c_S_B_fine_no_cbar_subsidy_b,c_B_B_fine_no_cbar_subsidy_b,c_M_B_fine_no_cbar_subsidy_b) = details_model(prices_no_cbar_subsidy_b,no_cbar_parameter,2,moments,1.0,foreign_supply_capital_subsidy_b);

#    Subsidy equilibrium with balanced budget + No cbar QS
prices_no_cbarQS_subsidy_b = [ 1.516811586169806,
0.20967990557054922,
0.16861791381781926]
no_cbarQS_parameter = copy(Baseline_parameter);
no_cbarQS_parameter.Q_S=0.0;
no_cbarQS_parameter.c̄_S = 0.0
(residual_goods_no_cbarQS_subsidy_b, stat_distr_no_cbarQS_subsidy_b, cons_fine_local_no_cbarQS_subsidy_b, a_prime_fine_local_no_cbarQS_subsidy_b,future_occupation_fine_local_no_cbarQS_subsidy_b,x_S_S_fine_no_cbarQS_subsidy_b,x_SC_fine_no_cbarQS_subsidy_b,x_BC_fine_no_cbarQS_subsidy_b, coeff_no_cbarQS_subsidy_b,
transaction_cost_loss_no_cbarQS_subsidy_b,nominal_GDP_no_cbarQS_subsidy_b,welfare_val_no_cbarQS_subsidy_b,Import_value_no_cbarQS_subsidy_b,Export_value_no_cbarQS_subsidy_b,current_worker_pop_no_cbarQS_subsidy_b,current_staple_pop_no_cbarQS_subsidy_b,current_cashcrop_pop_no_cbarQS_subsidy_b,
marketable_agr_surplus_share_no_cbarQS_subsidy_b,exportshare_cashcrop_no_cbarQS_subsidy_b,
fraction_model_no_cbarQS_subsidy_b,program_spending_no_cbarQS_subsidy_b,prod_value_improvement_no_cbarQS_subsidy_b,share_selling_increase_no_cbarQS_subsidy_b,exp_ratio_model_no_cbarQS_subsidy_b,mig_rate_model_no_cbarQS_subsidy_b,rural_pop_only_staples_model_no_cbarQS_subsidy_b,rural_pop_only_cashcrop_model_no_cbarQS_subsidy_b,
mean_land_share_to_staples_among_cc_model_no_cbarQS_subsidy_b,urban_rural_inc_ratio_model_no_cbarQS_subsidy_b,urban_rural_wealth_ratio_model_no_cbarQS_subsidy_b,urban_rural_consumption_ratio_model_no_cbarQS_subsidy_b,
p90_wealth_rural_no_cbarQS_subsidy_b,p90_wealth_urban_no_cbarQS_subsidy_b,p99_wealth_rural_no_cbarQS_subsidy_b,p99_wealth_urban_no_cbarQS_subsidy_b,p90_cons_tmp_no_cbarQS_subsidy_b,p90_income_tmp_no_cbarQS_subsidy_b,p99_cons_tmp_no_cbarQS_subsidy_b,p99_income_tmp_no_cbarQS_subsidy_b,
staple_productivity_no_cbarQS_subsidy_b,cashcrop_productivity_no_cbarQS_subsidy_b,manuf_productivity_no_cbarQS_subsidy_b,relative_land_to_staples_no_cbarQS_subsidy_b,relative_land_to_cashcrop_no_cbarQS_subsidy_b,share_constrained_cashcrop_no_cbarQS_subsidy_b,var_MPX_staples_S_no_cbarQS_subsidy_b,var_MPX_cashcrop_B_no_cbarQS_subsidy_b,var_MPX_cashcrop_S_no_cbarQS_subsidy_b,
share_constrained_staple_no_cbarQS_subsidy_b,APG_no_cbarQS_subsidy_b,urban_rural_consumption_ratio_model_real_no_cbarQS_subsidy_b,aggregate_consumption_no_cbarQS_subsidy_b,
worker_pop_effective_no_cbarQS_subsidy_b,prod_manuf_no_cbarQS_subsidy_b,total_entry_cost_no_cbarQS_subsidy_b,prod_staple_no_cbarQS_subsidy_b,prod_cashcrop_no_cbarQS_subsidy_b,input_staple_no_cbarQS_subsidy_b,input_cashcrop_no_cbarQS_subsidy_b,
total_maintenance_cost_no_cbarQS_subsidy_b,current_account_residual_no_cbarQS_subsidy_b,fraction_cashcrop_suboptimal_model_no_cbarQS_subsidy_b,V_saved_no_cbarQS_subsidy_b
,avg_labor_prod_rural_no_cbar_subsidy_b,avg_labor_prod_urban_no_cbar_subsidy_b,avg_agri_prod_rural_no_cbar_subsidy_b,avg_agri_prod_urban_no_cbar_subsidy_b,
var_MPX_cashcrop_no_cbarQS_subsidy_b ,var_MPX_no_cbarQS_subsidy_b ,TFP_no_cbarQS_subsidy_b ,YL_manuf_no_cbarQS_subsidy_b  , 
YL_agr_no_cbarQS_subsidy_b ,coeff_var_labor_prod_rural_no_cbarQS_subsidy_b ,coeff_var_labor_prod_urban_no_cbarQS_subsidy_b ,
coeff_var_agri_prod_rural_no_cbarQS_subsidy_b , coeff_var_agri_prod_urban_no_cbarQS_subsidy_b ,
p90_wealth_no_cbarQS_subsidy_b,p99_wealth_no_cbarQS_subsidy_b,p90_cons_no_cbarQS_subsidy_b_rural,p99_cons_no_cbarQS_subsidy_b_rural,p90_cons_no_cbarQS_subsidy_b_urban,p99_cons_no_cbarQS_subsidy_b_urban,p90_income_no_cbarQS_subsidy_b_rural,
p99_income_no_cbarQS_subsidy_b_rural,p90_income_no_cbarQS_subsidy_b_urban,p99_income_tmp_urban,
wealth_of_workers_no_cbarQS_subsidy_b,wealth_of_staples_no_cbarQS_subsidy_b,wealth_of_cashcrop_no_cbarQS_subsidy_b,
c_B_worker_sum_no_cbarQS_subsidy_b,c_B_staple_sum_no_cbarQS_subsidy_b,c_B_cashcrop_sum_no_cbarQS_subsidy_b,c_S_worker_sum_no_cbarQS_subsidy_b,c_S_staple_sum_no_cbarQS_subsidy_b ,c_S_cashcrop_sum_no_cbarQS_subsidy_b,
transaction_cost_staple_sum_no_cbarQS_subsidy_b,transaction_cost_cashcrop_sum_no_cbarQS_subsidy_b,transaction_cost_worker_sum_no_cbarQS_subsidy_b,c_M_worker_sum_no_cbarQS_subsidy_b,c_M_staple_sum_no_cbarQS_subsidy_b,c_M_cashcrop_sum_no_cbarQS_subsidy_b
,MPX_mean_log_no_cbarQS_subsidy_b, MPX_mean_staples_S_log_no_cbarQS_subsidy_b,MPX_mean_cashcrop_log_no_cbarQS_subsidy_b
, APland_mean_log_no_cbarQS_subsidy_b,APland_mean_cashcrop_log_no_cbarQS_subsidy_b, APland_mean_staples_S_log_no_cbarQS_subsidy_b,var_APland_no_cbarQS_subsidy_b,var_APland_cashcrop_no_cbarQS_subsidy_b,var_APland_staples_S_no_cbarQS_subsidy_b,
c_S_W_fine_no_cbarQS_subsidy_b,c_B_W_fine_no_cbarQS_subsidy_b,c_M_W_fine_no_cbarQS_subsidy_b,c_S_S_fine_no_cbarQS_subsidy_b,c_B_S_fine_no_cbarQS_subsidy_b,c_M_S_fine_no_cbarQS_subsidy_b,c_S_B_fine_no_cbarQS_subsidy_b,c_B_B_fine_no_cbarQS_subsidy_b,c_M_B_fine_no_cbarQS_subsidy_b) = details_model(prices_no_cbarQS_subsidy_b,no_cbarQS_parameter,2,moments,1.0,foreign_supply_capital_subsidy_b);


#    Subsidy equilibrium with balanced budget + No F_W
prices_no_F_W_subsidy_b = [ 1.4796471392365003,
0.14138775404436385,
0.1605285734494441]
no_F_W_parameter = copy(Baseline_parameter);
no_F_W_parameter.F_W= 0.0;
(residual_goods_no_F_W_subsidy_b, stat_distr_no_F_W_subsidy_b, cons_fine_local_no_F_W_subsidy_b, a_prime_fine_local_no_F_W_subsidy_b,future_occupation_fine_local_no_F_W_subsidy_b,x_S_S_fine_no_F_W_subsidy_b,x_SC_fine_no_F_W_subsidy_b,x_BC_fine_no_F_W_subsidy_b, coeff_no_F_W_subsidy_b,
transaction_cost_loss_no_F_W_subsidy_b,nominal_GDP_no_F_W_subsidy_b,welfare_val_no_F_W_subsidy_b,Import_value_no_F_W_subsidy_b,Export_value_no_F_W_subsidy_b,current_worker_pop_no_F_W_subsidy_b,current_staple_pop_no_F_W_subsidy_b,current_cashcrop_pop_no_F_W_subsidy_b,
marketable_agr_surplus_share_no_F_W_subsidy_b,exportshare_cashcrop_no_F_W_subsidy_b,
fraction_model_no_F_W_subsidy_b,program_spending_no_F_W_subsidy_b,prod_value_improvement_no_F_W_subsidy_b,share_selling_increase_no_F_W_subsidy_b,exp_ratio_model_no_F_W_subsidy_b,mig_rate_model_no_F_W_subsidy_b,rural_pop_only_staples_model_no_F_W_subsidy_b,rural_pop_only_cashcrop_model_no_F_W_subsidy_b,
mean_land_share_to_staples_among_cc_model_no_F_W_subsidy_b,urban_rural_inc_ratio_model_no_F_W_subsidy_b,urban_rural_wealth_ratio_model_no_F_W_subsidy_b,urban_rural_consumption_ratio_model_no_F_W_subsidy_b,
p90_wealth_rural_no_F_W_subsidy_b,p90_wealth_urban_no_F_W_subsidy_b,p99_wealth_rural_no_F_W_subsidy_b,p99_wealth_urban_no_F_W_subsidy_b,p90_cons_tmp_no_F_W_subsidy_b,p90_income_tmp_no_F_W_subsidy_b,p99_cons_tmp_no_F_W_subsidy_b,p99_income_tmp_no_F_W_subsidy_b,
staple_productivity_no_F_W_subsidy_b,cashcrop_productivity_no_F_W_subsidy_b,manuf_productivity_no_F_W_subsidy_b,relative_land_to_staples_no_F_W_subsidy_b,relative_land_to_cashcrop_no_F_W_subsidy_b,share_constrained_cashcrop_no_F_W_subsidy_b,var_MPX_staples_S_no_F_W_subsidy_b,var_MPX_cashcrop_B_no_F_W_subsidy_b,var_MPX_cashcrop_S_no_F_W_subsidy_b,
share_constrained_staple_no_F_W_subsidy_b,APG_no_F_W_subsidy_b,urban_rural_consumption_ratio_model_real_no_F_W_subsidy_b,aggregate_consumption_no_F_W_subsidy_b,
worker_pop_effective_no_F_W_subsidy_b,prod_manuf_no_F_W_subsidy_b,total_entry_cost_no_F_W_subsidy_b,prod_staple_no_F_W_subsidy_b,prod_cashcrop_no_F_W_subsidy_b,input_staple_no_F_W_subsidy_b,input_cashcrop_no_F_W_subsidy_b,
total_maintenance_cost_no_F_W_subsidy_b,current_account_residual_no_F_W_subsidy_b,fraction_cashcrop_suboptimal_model_no_F_W_subsidy_b,V_saved_no_F_W_subsidy_b
,avg_labor_prod_rural_no_F_W_subsidy_b,avg_labor_prod_urban_no_F_W_subsidy_b,avg_agri_prod_rural_no_F_W_subsidy_b,avg_agri_prod_urban_no_F_W_subsidy_b,
var_MPX_cashcrop_no_F_W_subsidy_b ,var_MPX_no_F_W_subsidy_b ,TFP_no_F_W_subsidy_b ,YL_manuf_no_F_W_subsidy_b  , 
YL_agr_no_F_W_subsidy_b ,coeff_var_labor_prod_rural_no_F_W_subsidy_b ,coeff_var_labor_prod_urban_no_F_W_subsidy_b ,
coeff_var_agri_prod_rural_no_F_W_subsidy_b , coeff_var_agri_prod_urban_no_F_W_subsidy_b, 
p90_wealth_no_F_W_subsidy_b,p99_wealth_no_F_W_subsidy_b,p90_cons_no_F_W_subsidy_b_rural,p99_cons_no_F_W_subsidy_b_rural,p90_cons_no_F_W_subsidy_b_urban,p99_cons_no_F_W_subsidy_b_urban,p90_income_no_F_W_subsidy_b_rural,
p99_income_no_F_W_subsidy_b_rural,p90_income_no_F_W_subsidy_b_urban,p99_income_no_F_W_subsidy_b_urban,wealth_of_workers_no_F_W_subsidy_b,
wealth_of_staples_no_F_W_subsidy_b,wealth_of_cashcrop_no_F_W_subsidy_b,
c_B_worker_sum_no_F_W_subsidy_b,c_B_staple_sum_no_F_W_subsidy_b,c_B_cashcrop_sum_no_F_W_subsidy_b,c_S_worker_sum_no_F_W_subsidy_b,c_S_staple_sum_no_F_W_subsidy_b ,c_S_cashcrop_sum_no_F_W_subsidy_b,
transaction_cost_staple_sum_no_F_W_subsidy_b,transaction_cost_cashcrop_sum_no_F_W_subsidy_b,transaction_cost_worker_sum_no_F_W_subsidy_b,c_M_worker_sum_no_F_W_subsidy_b,c_M_staple_sum_no_F_W_subsidy_b,c_M_cashcrop_sum_no_F_W_subsidy_b
,MPX_mean_log_no_F_W_subsidy_b, MPX_mean_staples_S_log_no_F_W_subsidy_b,MPX_mean_cashcrop_log_no_F_W_subsidy_b
, APland_mean_log_no_F_W_subsidy_b,APland_mean_cashcrop_log_no_F_W_subsidy_b, APland_mean_staples_S_log_no_F_W_subsidy_b,var_APland_no_F_W_subsidy_b,var_APland_cashcrop_no_F_W_subsidy_b,var_APland_staples_S_no_F_W_subsidy_b,
c_S_W_fine_no_F_W_subsidy_b,c_B_W_fine_no_F_W_subsidy_b,c_M_W_fine_no_F_W_subsidy_b,c_S_S_fine_no_F_W_subsidy_b,c_B_S_fine_no_F_W_subsidy_b,c_M_S_fine_no_F_W_subsidy_b,c_S_B_fine_no_F_W_subsidy_b,c_B_B_fine_no_F_W_subsidy_b,c_M_B_fine_no_F_W_subsidy_b) = details_model(prices_no_F_W_subsidy_b,no_F_W_parameter,2,moments,1.0,foreign_supply_capital_subsidy_b);


#    Subsidy equilibrium with balanced budget + No FM_B
prices_no_FM_B_subsidy_b = [   1.4805615707512334,
0.26962104767974854,
0.16687937137711306]
no_FM_B_parameter = copy(Baseline_parameter);
no_FM_B_parameter.FM_B= 0.0;
(residual_goods_no_FM_B_subsidy_b, stat_distr_no_FM_B_subsidy_b, cons_fine_local_no_FM_B_subsidy_b, a_prime_fine_local_no_FM_B_subsidy_b,future_occupation_fine_local_no_FM_B_subsidy_b,x_S_S_fine_no_FM_B_subsidy_b,x_SC_fine_no_FM_B_subsidy_b,x_BC_fine_no_FM_B_subsidy_b, coeff_no_FM_B_subsidy_b,
transaction_cost_loss_no_FM_B_subsidy_b,nominal_GDP_no_FM_B_subsidy_b,welfare_val_no_FM_B_subsidy_b,Import_value_no_FM_B_subsidy_b,Export_value_no_FM_B_subsidy_b,current_worker_pop_no_FM_B_subsidy_b,current_staple_pop_no_FM_B_subsidy_b,current_cashcrop_pop_no_FM_B_subsidy_b,
marketable_agr_surplus_share_no_FM_B_subsidy_b,exportshare_cashcrop_no_FM_B_subsidy_b,
fraction_model_no_FM_B_subsidy_b,program_spending_no_FM_B_subsidy_b,prod_value_improvement_no_FM_B_subsidy_b,share_selling_increase_no_FM_B_subsidy_b,exp_ratio_model_no_FM_B_subsidy_b,mig_rate_model_no_FM_B_subsidy_b,rural_pop_only_staples_model_no_FM_B_subsidy_b,rural_pop_only_cashcrop_model_no_FM_B_subsidy_b,
mean_land_share_to_staples_among_cc_model_no_FM_B_subsidy_b,urban_rural_inc_ratio_model_no_FM_B_subsidy_b,urban_rural_wealth_ratio_model_no_FM_B_subsidy_b,urban_rural_consumption_ratio_model_no_FM_B_subsidy_b,
p90_wealth_rural_no_FM_B_subsidy_b,p90_wealth_urban_no_FM_B_subsidy_b,p99_wealth_rural_no_FM_B_subsidy_b,p99_wealth_urban_no_FM_B_subsidy_b,p90_cons_tmp_no_FM_B_subsidy_b,p90_income_tmp_no_FM_B_subsidy_b,p99_cons_tmp_no_FM_B_subsidy_b,p99_income_tmp_no_FM_B_subsidy_b,
staple_productivity_no_FM_B_subsidy_b,cashcrop_productivity_no_FM_B_subsidy_b,manuf_productivity_no_FM_B_subsidy_b,relative_land_to_staples_no_FM_B_subsidy_b,relative_land_to_cashcrop_no_FM_B_subsidy_b,share_constrained_cashcrop_no_FM_B_subsidy_b,var_MPX_staples_S_no_FM_B_subsidy_b,var_MPX_cashcrop_B_no_FM_B_subsidy_b,var_MPX_cashcrop_S_no_FM_B_subsidy_b,
share_constrained_staple_no_FM_B_subsidy_b,APG_no_FM_B_subsidy_b,urban_rural_consumption_ratio_model_real_no_FM_B_subsidy_b,aggregate_consumption_no_FM_B_subsidy_b,
worker_pop_effective_no_FM_B_subsidy_b,prod_manuf_no_FM_B_subsidy_b,total_entry_cost_no_FM_B_subsidy_b,prod_staple_no_FM_B_subsidy_b,prod_cashcrop_no_FM_B_subsidy_b,input_staple_no_FM_B_subsidy_b,input_cashcrop_no_FM_B_subsidy_b,
total_maintenance_cost_no_FM_B_subsidy_b,current_account_residual_no_FM_B_subsidy_b,fraction_cashcrop_suboptimal_model_no_FM_B_subsidy_b,V_saved_no_FM_B_subsidy_b
,avg_labor_prod_rural_no_FM_B_subsidy_b,avg_labor_prod_urban_no_FM_B_subsidy_b,avg_agri_prod_rural_no_FM_B_subsidy_b,avg_agri_prod_urban_no_FM_B_subsidy_b,
var_MPX_cashcrop_no_FM_B_subsidy_b ,var_MPX_no_FM_B_subsidy_b ,TFP_no_FM_B_subsidy_b ,YL_manuf_no_FM_B_subsidy_b  , 
YL_agr_no_FM_B_subsidy_b ,coeff_var_labor_prod_rural_no_FM_B_subsidy_b ,coeff_var_labor_prod_urban_no_FM_B_subsidy_b ,
coeff_var_agri_prod_rural_no_FM_B_subsidy_b , coeff_var_agri_prod_urban_no_FM_B_subsidy_b, 
p90_wealth_no_FM_B_subsidy_b,p99_wealth_no_FM_B_subsidy_b,p90_cons_no_FM_B_subsidy_b_rural,p99_cons_no_FM_B_subsidy_b_rural,p90_cons_no_FM_B_subsidy_b_urban,p99_cons_no_FM_B_subsidy_b_urban,p90_income_no_FM_B_subsidy_b_rural,
p99_income_no_FM_B_subsidy_b_rural,p90_income_no_FM_B_subsidy_b_urban,p99_income_no_FM_B_subsidy_b_urban,
wealth_of_workers_no_FM_B_subsidy_b,wealth_of_staples_no_FM_B_subsidy_b,wealth_of_cashcrop_no_FM_B_subsidy_b,
c_B_worker_sum_no_FM_B_subsidy_b,c_B_staple_sum_no_FM_B_subsidy_b,c_B_cashcrop_sum_no_FM_B_subsidy_b,c_S_worker_sum_no_FM_B_subsidy_b,c_S_staple_sum_no_FM_B_subsidy_b ,c_S_cashcrop_sum_no_FM_B_subsidy_b,
transaction_cost_staple_sum_no_FM_B_subsidy_b,transaction_cost_cashcrop_sum_no_FM_B_subsidy_b,transaction_cost_worker_sum_no_FM_B_subsidy_b,c_M_worker_sum_no_FM_B_subsidy_b,c_M_staple_sum_no_FM_B_subsidy_b,c_M_cashcrop_sum_no_FM_B_subsidy_b
,MPX_mean_log_no_FM_B_subsidy_b, MPX_mean_staples_S_log_no_FM_B_subsidy_b,MPX_mean_cashcrop_log_no_FM_B_subsidy_b
, APland_mean_log_no_FM_B_subsidy_b,APland_mean_cashcrop_log_no_FM_B_subsidy_b, APland_mean_staples_S_log_no_FM_B_subsidy_b,var_APland_no_FM_B_subsidy_b,var_APland_cashcrop_no_FM_B_subsidy_b,var_APland_staples_S_no_FM_B_subsidy_b,
c_S_W_fine_no_FM_B_subsidy_b,c_B_W_fine_no_FM_B_subsidy_b,c_M_W_fine_no_FM_B_subsidy_b,c_S_S_fine_no_FM_B_subsidy_b,c_B_S_fine_no_FM_B_subsidy_b,c_M_S_fine_no_FM_B_subsidy_b,c_S_B_fine_no_FM_B_subsidy_b,c_B_B_fine_no_FM_B_subsidy_b,c_M_B_fine_no_FM_B_subsidy_b) = details_model(prices_no_FM_B_subsidy_b,no_FM_B_parameter,2,moments,1.0,foreign_supply_capital_subsidy_b);

#    Subsidy equilibrium with balanced budget + No κ
prices_no_κ_subsidy_b =  [   1.514071468437403,
0.2803054849913619,
0.2617382599984764]; #!
no_κ_parameter_subsidy_b = copy(Baseline_parameter);
no_κ_parameter_subsidy_b.κ= 10.0;
(residual_goods_no_κ_subsidy_b, stat_distr_no_κ_subsidy_b, cons_fine_local_no_κ_subsidy_b, a_prime_fine_local_no_κ_subsidy_b,future_occupation_fine_local_no_κ_subsidy_b,x_S_S_fine_no_κ_subsidy_b,x_SC_fine_no_κ_subsidy_b,x_BC_fine_no_κ_subsidy_b, coeff_no_κ_subsidy_b,
transaction_cost_loss_no_κ_subsidy_b,nominal_GDP_no_κ_subsidy_b,welfare_val_no_κ_subsidy_b,Import_value_no_κ_subsidy_b,Export_value_no_κ_subsidy_b,current_worker_pop_no_κ_subsidy_b,current_staple_pop_no_κ_subsidy_b,current_cashcrop_pop_no_κ_subsidy_b,
marketable_agr_surplus_share_no_κ_subsidy_b,exportshare_cashcrop_no_κ_subsidy_b,
fraction_model_no_κ_subsidy_b,program_spending_no_κ_subsidy_b,prod_value_improvement_no_κ_subsidy_b,share_selling_increase_no_κ_subsidy_b,exp_ratio_model_no_κ_subsidy_b,mig_rate_model_no_κ_subsidy_b,rural_pop_only_staples_model_no_κ_subsidy_b,rural_pop_only_cashcrop_model_no_κ_subsidy_b,
mean_land_share_to_staples_among_cc_model_no_κ_subsidy_b,urban_rural_inc_ratio_model_no_κ_subsidy_b,urban_rural_wealth_ratio_model_no_κ_subsidy_b,urban_rural_consumption_ratio_model_no_κ_subsidy_b,
p90_wealth_rural_no_κ_subsidy_b,p90_wealth_urban_no_κ_subsidy_b,p99_wealth_rural_no_κ_subsidy_b,p99_wealth_urban_no_κ_subsidy_b,p90_cons_tmp_no_κ_subsidy_b,p90_income_tmp_no_κ_subsidy_b,p99_cons_tmp_no_κ_subsidy_b,p99_income_tmp_no_κ_subsidy_b,
staple_productivity_no_κ_subsidy_b,cashcrop_productivity_no_κ_subsidy_b,manuf_productivity_no_κ_subsidy_b,relative_land_to_staples_no_κ_subsidy_b,relative_land_to_cashcrop_no_κ_subsidy_b,share_constrained_cashcrop_no_κ_subsidy_b,var_MPX_staples_S_no_κ_subsidy_b,var_MPX_cashcrop_B_no_κ_subsidy_b,var_MPX_cashcrop_S_no_κ_subsidy_b,
share_constrained_staple_no_κ_subsidy_b,APG_no_κ_subsidy_b,urban_rural_consumption_ratio_model_real_no_κ_subsidy_b,aggregate_consumption_no_κ_subsidy_b,
worker_pop_effective_no_κ_subsidy_b,prod_manuf_no_κ_subsidy_b,total_entry_cost_no_κ_subsidy_b,prod_staple_no_κ_subsidy_b,prod_cashcrop_no_κ_subsidy_b,input_staple_no_κ_subsidy_b,input_cashcrop_no_κ_subsidy_b,
total_maintenance_cost_no_κ_subsidy_b,current_account_residual_no_κ_subsidy_b,fraction_cashcrop_suboptimal_model_no_κ_subsidy_b,V_saved_no_κ_subsidy_b
,avg_labor_prod_rural_no_κ_subsidy_b,avg_labor_prod_urban_no_κ_subsidy_b,avg_agri_prod_rural_no_κ_subsidy_b,avg_agri_prod_urban_no_κ_subsidy_b,
var_MPX_cashcrop_no_κ_subsidy_b ,var_MPX_no_κ_subsidy_b ,TFP_no_κ_subsidy_b ,YL_manuf_no_κ_subsidy_b  , 
YL_agr_no_κ_subsidy_b ,coeff_var_labor_prod_rural_no_κ_subsidy_b ,coeff_var_labor_prod_urban_no_κ_subsidy_b ,
coeff_var_agri_prod_rural_no_κ_subsidy_b , coeff_var_agri_prod_urban_no_κ_subsidy_b ,
p90_wealth_no_κ_subsidy_b,p99_wealth_no_κ_subsidy_b,p90_cons_no_κ_subsidy_b_rural,p99_cons_no_κ_subsidy_b_rural,p90_cons_no_κ_subsidy_b_urban,p99_cons_no_κ_subsidy_b_urban,p90_income_no_κ_subsidy_b_rural,
p99_income_no_κ_subsidy_b_rural,p90_income_no_κ_subsidy_b_urban,p99_income_no_κ_subsidy_b_urban,
wealth_of_workers_no_κ_subsidy_b,wealth_of_staples_no_κ_subsidy_b,wealth_of_cashcrop_no_κ_subsidy_b,
c_B_worker_sum_no_κ_subsidy_b,c_B_staple_sum_no_κ_subsidy_b,c_B_cashcrop_sum_no_κ_subsidy_b,c_S_worker_sum_no_κ_subsidy_b,c_S_staple_sum_no_κ_subsidy_b ,c_S_cashcrop_sum_no_κ_subsidy_b,
transaction_cost_staple_sum_no_κ_subsidy_b,transaction_cost_cashcrop_sum_no_κ_subsidy_b,transaction_cost_worker_sum_no_κ_subsidy_b,c_M_worker_sum_no_κ_subsidy_b,c_M_staple_sum_no_κ_subsidy_b,c_M_cashcrop_sum_no_κ_subsidy_b
,MPX_mean_log_no_κ_subsidy_b, MPX_mean_staples_S_log_no_κ_subsidy_b,MPX_mean_cashcrop_log_no_κ_subsidy_b
, APland_mean_log_no_κ_subsidy_b,APland_mean_cashcrop_log_no_κ_subsidy_b, APland_mean_staples_S_log_no_κ_subsidy_b,var_APland_no_κ_subsidy_b,var_APland_cashcrop_no_κ_subsidy_b,var_APland_staples_S_no_κ_subsidy_b,
c_S_W_fine_no_κ_subsidy_b,c_B_W_fine_no_κ_subsidy_b,c_M_W_fine_no_κ_subsidy_b,c_S_S_fine_no_κ_subsidy_b,c_B_S_fine_no_κ_subsidy_b,c_M_S_fine_no_κ_subsidy_b,c_S_B_fine_no_κ_subsidy_b,c_B_B_fine_no_κ_subsidy_b,c_M_B_fine_no_κ_subsidy_b) = details_model(prices_no_κ_subsidy_b,no_κ_parameter_subsidy_b,2,moments,1.0,foreign_supply_capital_subsidy_b);



# Infrastructure 
infra_parameter_nsp_nb = copy(Baseline_parameter);
infra_parameter_nsp_nb.τ_S=-0;
infra_parameter_nsp_nb.Q_S = 0.7;
prices_inf =  [ 1.190206255203472,
0.1745527606893067];#,2.505662787966534
(residual_goods_inf, stat_distr_inf, cons_fine_local_inf, a_prime_fine_local_inf,future_occupation_fine_local_inf,x_S_S_fine_inf,x_SC_fine_inf,x_BC_fine_inf, coeff_inf,
transaction_cost_loss_inf,nominal_GDP_inf,welfare_val_inf,Import_value_inf,Export_value_inf,current_worker_pop_inf,current_staple_pop_inf,current_cashcrop_pop_inf,
marketable_agr_surplus_share_inf,exportshare_cashcrop_inf,
fraction_model_inf,program_spending_inf,prod_value_improvement_inf,share_selling_increase_inf,exp_ratio_model_inf,mig_rate_model_inf,rural_pop_only_staples_model_inf,rural_pop_only_cashcrop_model_inf,
mean_land_share_to_staples_among_cc_model_inf,urban_rural_inc_ratio_model_inf,urban_rural_wealth_ratio_model_inf,urban_rural_consumption_ratio_model_inf,
p90_wealth_rural_inf,p90_wealth_urban_inf,p99_wealth_rural_inf,p99_wealth_urban_inf,p90_cons_tmp_inf,p90_income_tmp_inf,p99_cons_tmp_inf,p99_income_tmp_inf,
staple_productivity_inf,cashcrop_productivity_inf,manuf_productivity_inf,relative_land_to_staples_inf,relative_land_to_cashcrop_inf,share_constrained_cashcrop_inf,var_MPX_staples_S_inf,var_MPX_cashcrop_B_inf,var_MPX_cashcrop_S_inf,
share_constrained_staple_inf,APG_inf,urban_rural_consumption_ratio_model_real_inf,aggregate_consumption_inf,
worker_pop_effective_inf,prod_manuf_inf,total_entry_cost_inf,prod_staple_inf,prod_cashcrop_inf,input_staple_inf,input_cashcrop_inf,
total_maintenance_cost_inf,current_account_residual_inf,fraction_cashcrop_suboptimal_model_inf,V_saved_inf,
avg_labor_prod_rural_inf,avg_labor_prod_urban_inf,avg_agri_prod_rural_inf,avg_agri_prod_urban_inf,
var_MPX_cashcrop_inf ,var_MPX_inf ,TFP_inf ,YL_manuf_inf  , 
YL_agr_inf ,coeff_var_labor_prod_rural_inf ,coeff_var_labor_prod_urban_inf ,
coeff_var_agri_prod_rural_inf , coeff_var_agri_prod_urban_inf ,
p90_wealth_inf,p99_wealth_inf,p90_cons_inf_rural,p99_cons_inf_rural,p90_cons_inf_urban,p99_cons_inf_urban,p90_income_inf_rural,
p99_income_inf_rural,p90_income_inf_urban,p99_income_inf_urban,
wealth_of_workers_inf,wealth_of_staples_inf,wealth_of_cashcrop_inf,
c_B_worker_sum_inf,c_B_staple_sum_inf,c_B_cashcrop_sum_inf,c_S_worker_sum_inf,c_S_staple_sum_inf ,c_S_cashcrop_sum_inf,
transaction_cost_staple_sum_inf,transaction_cost_cashcrop_sum_inf,transaction_cost_worker_sum_inf,c_M_worker_sum_inf,c_M_staple_sum_inf,c_M_cashcrop_sum_inf
,MPX_mean_log_inf, MPX_mean_staples_S_log_inf,MPX_mean_cashcrop_log_inf
, APland_mean_log_inf,APland_mean_cashcrop_log_inf, APland_mean_staples_S_log_inf,var_APland_inf,var_APland_cashcrop_inf,var_APland_staples_S_inf,
c_S_W_fine_inf,c_B_W_fine_inf,c_M_W_fine_inf,c_S_S_fine_inf,c_B_S_fine_inf,c_M_S_fine_inf,c_S_B_fine_inf,c_B_B_fine_inf,c_M_B_fine_inf) = details_model(prices_inf,infra_parameter_nsp_nb,2,moments,0.0,foreign_supply_capital_subsidy_b);

#program_spending_subsidy_nb * nominal_GDP_subsidy_nb - transaction_cost_loss_inf/ (infra_parameter_nsp_nb.Q_S)*(Baseline_parameter.Q_S-infra_parameter_nsp_nb.Q_S)
transaction_cost_loss_inf*prod_staple_inf/ (infra_parameter_nsp_nb.Q_S)*(Baseline_parameter.Q_S-infra_parameter_nsp_nb.Q_S) / nominal_GDP_inf

infra_parameter_sp_nb = copy(Baseline_parameter);
infra_parameter_sp_nb.τ_S=-0;
infra_parameter_sp_nb.Q_S = 0.7;
infra_parameter_sp_nb.F_W = Baseline_parameter.F_W * (infra_parameter_sp_nb.Q_S / Baseline_parameter.Q_S)
prices_inf_sp = [1.1723046609628551,
0.12471486958333933]#
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
        total_maintenance_cost_inf_sp,current_account_residual_inf_sp,fraction_cashcrop_suboptimal_model_inf_sp,V_saved_inf_sp,
        avg_labor_prod_rural_inf_sp,avg_labor_prod_urban_inf_sp,avg_agri_prod_rural_inf_sp,avg_agri_prod_urban_inf_sp,
        var_MPX_cashcrop_inf_sp ,var_MPX_inf_sp ,TFP_inf_sp ,YL_manuf_inf_sp  , 
        YL_agr_inf_sp ,coeff_var_labor_prod_rural_inf_sp ,coeff_var_labor_prod_urban_inf_sp ,
        coeff_var_agri_prod_rural_inf_sp , coeff_var_agri_prod_urban_inf_sp ,
        p90_wealth_inf_sp,p99_wealth_inf_sp,p90_cons_inf_sp_rural,p99_cons_inf_sp_rural,p90_cons_inf_sp_urban,p99_cons_inf_sp_urban,p90_income_inf_sp_rural,
        p99_income_inf_sp_rural,p90_income_inf_sp_urban,p99_income_inf_sp_urban,
        wealth_of_workers_inf_sp,wealth_of_staples_inf_sp,wealth_of_cashcrop_inf_sp,
        c_B_worker_sum_inf_sp,c_B_staple_sum_inf_sp,c_B_cashcrop_sum_inf_sp,c_S_worker_sum_inf_sp,c_S_staple_sum_inf_sp ,c_S_cashcrop_sum_inf_sp,
        transaction_cost_staple_sum_inf_sp,transaction_cost_cashcrop_sum_inf_sp,transaction_cost_worker_sum_inf_sp,c_M_worker_sum_inf_sp,c_M_staple_sum_inf_sp,c_M_cashcrop_sum_inf_sp
        ,MPX_mean_log_inf_sp, MPX_mean_staples_S_log_inf_sp,MPX_mean_cashcrop_log_inf_sp
        , APland_mean_log_inf_sp,APland_mean_cashcrop_log_inf_sp, APland_mean_staples_S_log_inf_sp,var_APland_inf_sp,var_APland_cashcrop_inf_sp,var_APland_staples_S_inf_sp,
c_S_W_fine_inf_sp,c_B_W_fine_inf_sp,c_M_W_fine_inf_sp,c_S_S_fine_inf_sp,c_B_S_fine_inf_sp,c_M_S_fine_inf_sp,c_S_B_fine_inf_sp,c_B_B_fine_inf_sp,c_M_B_fine_inf_sp) = details_model(prices_inf_sp,infra_parameter_sp_nb,2,moments,0.0,foreign_supply_capital_subsidy_b);

# Obtain the results for optimal subsidy:

# # This part is only needed for different tauS vectors
# no_taus1 = 11
# prices_subsidy_b_grid1 = ones(3,no_taus1);
# prices_subsidy_b_grid1[1,:] = [1.50193 ,  1.41901 ,   1.37031  ,  1.32695  ,  1.29374  ,  1.2656  ,   1.24367  ,  1.22488 ,  1.20621  ,  1.18598 ,   1.1761 ];
# prices_subsidy_b_grid1[2,:] = [0.252278 , 0.23388  ,  0.229202 ,  0.215522 ,  0.217217 ,  0.205816 ,  0.205173 ,  0.203484 , 0.207593 ,  0.192035 ,  0.192054];
# prices_subsidy_b_grid1[3,:] = [0.152596 , 0.0931243 , 0.0682528 , 0.0489459 , 0.0348802 , 0.0247181 , 0.0176974 , 0.011634 , 0.0063424 , 0.0021779 , 0.0 ];


# τ_grid1 = zeros(no_taus1);
# τ_grid1[1:(no_taus1 - 1) ] = range(Baseline_parameter.τ_S,-0.05,no_taus1-1)

# no_taus2= 21
# prices_subsidy_b_grid2 = ones(3,no_taus2);
# prices_subsidy_b_grid2[1,:] = [1.50193 ,  1.46548 ,  1.42793  , 1.40183  , 1.37866  ,  1.35987  ,  1.34164  ,  1.32354  ,  1.30643 ,   1.29362  ,  1.28082   , 1.26905  ,  1.25875  ,  1.2467   ,  1.23701  ,  1.22909  ,  1.21865 ,   1.21536  ,   1.20605  ,   1.19481   ,  1.1761  ];       
# prices_subsidy_b_grid2[2,:] = [0.252278 , 0.233978 , 0.233412 , 0.232739 , 0.232361 ,  0.213076  , 0.212905 ,  0.215287 ,  0.204994  , 0.204275 ,  0.204633 ,  0.205771 ,  0.207009  , 0.207411 ,  0.209514 ,  0.209253  , 0.192362 ,  0.189848  ,  0.189807  ,  0.191641  ,  0.192054   ];    
# prices_subsidy_b_grid2[3,:] = [0.152596 , 0.129457 , 0.099764 , 0.084152 , 0.0725594 , 0.0642034 , 0.0556398 , 0.0474308 , 0.0409965 , 0.0358801 , 0.0299909 , 0.0258702 , 0.0222282 , 0.0189747 , 0.0157095 , 0.0131466 , 0.0108565 , 0.00901757 , 0.00679134 , 0.00452941 , 0.0 ];


# τ_grid2 = zeros(no_taus2);
# τ_grid2[1:(no_taus2 - 1) ] = range(Baseline_parameter.τ_S,-0.1,no_taus2-1)

# no_taus = no_taus2 + no_taus1 - 2;
# τ_grid = zeros(no_taus);
# prices_subsidy_b_grid = ones(3,no_taus); 
# prices_subsidy_b_grid[:,1:(no_taus1 - 1)] = prices_subsidy_b_grid1[:,1:(no_taus1 - 1) ] ;
# prices_subsidy_b_grid[:,no_taus1:(end)] = prices_subsidy_b_grid2[:,2:(no_taus2) ] ;
# τ_grid[1:(no_taus1 - 1)] = τ_grid1[1:(no_taus1 - 1) ] ;
# τ_grid[no_taus1:(end - 1)] = τ_grid2[2:(no_taus2 - 1) ] ;

# ordered_tauS_index = sortperm(τ_grid);
# τ_grid_ordered = τ_grid[ordered_tauS_index]
# prices_subsidy_b_grid_ordered = prices_subsidy_b_grid[:,ordered_tauS_index]
# τ_grid = copy(τ_grid_ordered);
# prices_subsidy_b_grid = copy(prices_subsidy_b_grid_ordered);

### Now the two objects are saved:
prices_subsidy_b_grid = [1.5182148598312672 1.4748615289082854 1.4389313172837825 1.410732954524286 1.390264709095618 1.3688114420290358 1.3436655727518898 1.3271570823425125 1.312778814101722 1.3003105888992221 1.284946210796012 1.2714091997394257 1.2599461971286845 1.2479933342455625 1.2388104401541415 1.2294244830748497 1.2217938803468038 1.2132098804982514 1.2099699888601618 1.1995655457878038 1.1882831767847994; 0.27262899815531444 0.26338424759605666 0.245822286882918 0.2442836837960706 0.21752367635207415 0.21670507564964162 0.2061107155376924 0.2036475392226677 0.20098618728681622 0.19862949627869536 0.19659831128032512 0.19476856375812163 0.19323999476869438 0.18933334499100343 0.19016660308519126 0.18807703128020226 0.1869710641961219 0.1857564570225632 0.18361768757377586 0.18414175648950032 0.18349666056314162; 0.18860594848809528 0.14681832007063747 0.12372332717865375 0.10211507272581569 0.09328546623799094 0.0770503430189069 0.06580378586580023 0.055857350780308554 0.04822701987252334 0.04189534576885402 0.03530977698358126 0.029721337552662572 0.0251658618411578 0.02075905911514478 0.016733799327583475 0.013639766449854727 0.010610113658055436 0.007772421378081955 0.005375731191350952 0.002777150926899589 0.0]
no_taus = 21;
τ_grid = range(Baseline_parameter.τ_S,-0.05,length = no_taus-1)
τ_grid1 = zeros(no_taus);
τ_grid1[1:(no_taus-1)] = τ_grid;
τ_grid = τ_grid1;
residual_goods_subsidy_b_grid = ones(6,no_taus);
stat_distr_subsidy_b_grid = ones(Baseline_parameter.ns_fine*3,no_taus);
cons_fine_local_subsidy_b_grid = ones(Baseline_parameter.ns_fine,3,no_taus);
a_prime_fine_local_subsidy_b_grid = ones(Baseline_parameter.ns_fine,3,no_taus);
future_occupation_fine_local_subsidy_b_grid = ones(Baseline_parameter.ns_fine,3,no_taus);
x_S_S_fine_subsidy_b_grid = ones(Baseline_parameter.ns_fine,3,no_taus);
x_SC_fine_subsidy_b_grid = ones(Baseline_parameter.ns_fine,3,no_taus);
x_BC_fine_subsidy_b_grid = ones(Baseline_parameter.ns_fine,3,no_taus);
coeff_subsidy_b_grid = ones(Baseline_parameter.ns,3,no_taus);
transaction_cost_loss_subsidy_b_grid = ones(no_taus);
nominal_GDP_subsidy_b_grid= ones(no_taus);
welfare_val_subsidy_b_grid =  ones(Baseline_parameter.ns_fine*3,no_taus);
Import_value_subsidy_b_grid = ones(no_taus);
Export_value_subsidy_b_grid = ones(no_taus);
current_worker_pop_subsidy_b_grid = ones(no_taus);
current_staple_pop_subsidy_b_grid = ones(no_taus);
current_cashcrop_pop_subsidy_b_grid = ones(no_taus);
marketable_agr_surplus_share_subsidy_b_grid = ones(no_taus);
exportshare_cashcrop_subsidy_b_grid = ones(no_taus);
fraction_model_subsidy_b_grid= ones(no_taus);
program_spending_subsidy_b_grid= ones(no_taus);
prod_value_improvement_subsidy_b_grid= ones(no_taus);
share_selling_increase_subsidy_b_grid= ones(no_taus);
exp_ratio_model_subsidy_b_grid= ones(no_taus);
mig_rate_model_subsidy_b_grid= ones(no_taus);
rural_pop_only_staples_model_subsidy_b_grid= ones(no_taus);
rural_pop_only_cashcrop_model_subsidy_b_grid= ones(no_taus);
mean_land_share_to_staples_among_cc_model_subsidy_b_grid= ones(no_taus);
urban_rural_inc_ratio_model_subsidy_b_grid= ones(no_taus);
urban_rural_wealth_ratio_model_subsidy_b_grid= ones(no_taus);
urban_rural_consumption_ratio_model_subsidy_b_grid= ones(no_taus);
p90_wealth_rural_subsidy_b_grid= ones(no_taus);
p90_wealth_urban_subsidy_b_grid= ones(no_taus);
p99_wealth_rural_subsidy_b_grid= ones(no_taus);
p99_wealth_urban_subsidy_b_grid= ones(no_taus);
p90_cons_tmp_subsidy_b_grid= ones(no_taus);
p90_income_tmp_subsidy_b_grid= ones(no_taus);
p99_cons_tmp_subsidy_b_grid= ones(no_taus);
p99_income_tmp_subsidy_b_grid= ones(no_taus);
staple_productivity_subsidy_b_grid= ones(no_taus);
cashcrop_productivity_subsidy_b_grid= ones(no_taus);
manuf_productivity_subsidy_b_grid= ones(no_taus);
relative_land_to_staples_subsidy_b_grid= ones(no_taus);
relative_land_to_cashcrop_subsidy_b_grid= ones(no_taus);
share_constrained_cashcrop_subsidy_b_grid= ones(no_taus);
var_MPX_staples_S_subsidy_b_grid= ones(no_taus);
var_MPX_cashcrop_B_subsidy_b_grid= ones(no_taus);
var_MPX_cashcrop_S_subsidy_b_grid= ones(no_taus);
share_constrained_staple_subsidy_b_grid= ones(no_taus);
APG_subsidy_b_grid= ones(no_taus);
urban_rural_consumption_ratio_model_real_subsidy_b_grid= ones(no_taus);
aggregate_consumption_subsidy_b_grid= ones(no_taus);
worker_pop_effective_subsidy_b_grid= ones(no_taus);
prod_manuf_subsidy_b_grid= ones(no_taus);
total_entry_cost_subsidy_b_grid= ones(no_taus);
prod_staple_subsidy_b_grid= ones(no_taus);
prod_cashcrop_subsidy_b_grid= ones(no_taus);
input_staple_subsidy_b_grid= ones(no_taus);
input_cashcrop_subsidy_b_grid= ones(no_taus);
total_maintenance_cost_subsidy_b_grid= ones(no_taus);
current_account_residual_subsidy_b_grid= ones(no_taus);
fraction_cashcrop_suboptimal_model_subsidy_b_grid = ones(no_taus);
welfare_subsidy_b_grid= ones(no_taus);
V_saved_b_grid = ones(Baseline_parameter.ns_fine, 3, no_taus);
avg_labor_prod_rural_subsidy_b_grid= ones(no_taus);
avg_labor_prod_urban_subsidy_b_grid= ones(no_taus);
avg_agri_prod_rural_subsidy_b_grid= ones(no_taus);
avg_agri_prod_urban_subsidy_b_grid= ones(no_taus);
var_MPX_cashcrop_subsidy_b_grid= ones(no_taus);
var_MPX_subsidy_b_grid= ones(no_taus);
TFP_subsidy_b_grid= ones(no_taus);
YL_manuf_subsidy_b_grid= ones(no_taus);
YL_agr_subsidy_b_grid= ones(no_taus);
coeff_var_labor_prod_rural_subsidy_b_grid= ones(no_taus);
coeff_var_labor_prod_urban_subsidy_b_grid= ones(no_taus);
coeff_var_agri_prod_rural_subsidy_b_grid= ones(no_taus);
coeff_var_agri_prod_urban_subsidy_b_grid= ones(no_taus);
p90_wealth_subsidy_b_grid= ones(no_taus);
p99_wealth_subsidy_b_grid= ones(no_taus);
p90_cons_subsidy_b_grid_rural= ones(no_taus);
p99_cons_subsidy_b_grid_rural= ones(no_taus);
p90_cons_subsidy_b_grid_urban= ones(no_taus);
p99_cons_subsidy_b_grid_urban= ones(no_taus);
p90_income_subsidy_b_grid_rural= ones(no_taus);
p99_income_subsidy_b_grid_rural= ones(no_taus);
p90_income_subsidy_b_grid_urban= ones(no_taus);
p99_income_subsidy_b_grid_urban= ones(no_taus);
wealth_of_workers_subsidy_b_grid= ones(no_taus);
wealth_of_staples_subsidy_b_grid= ones(no_taus);
wealth_of_cashcrop_subsidy_b_grid= ones(no_taus);
c_B_worker_sum_subsidy_b_grid= ones(no_taus);
c_B_staple_sum_subsidy_b_grid= ones(no_taus);
c_B_cashcrop_sum_subsidy_b_grid= ones(no_taus);
c_S_worker_sum_subsidy_b_grid= ones(no_taus);
c_S_staple_sum_subsidy_b_grid = ones(no_taus);
c_S_cashcrop_sum_subsidy_b_grid= ones(no_taus);
transaction_cost_staple_sum_subsidy_b_grid= ones(no_taus);
transaction_cost_cashcrop_sum_subsidy_b_grid= ones(no_taus);
transaction_cost_worker_sum_subsidy_b_grid= ones(no_taus);
c_M_worker_sum_subsidy_b_grid= ones(no_taus);
c_M_staple_sum_subsidy_b_grid= ones(no_taus);
c_M_cashcrop_sum_subsidy_b_grid= ones(no_taus);
MPX_mean_log_subsidy_b_grid= ones(no_taus);
MPX_mean_staples_S_log_subsidy_b_grid= ones(no_taus);
MPX_mean_cashcrop_log_subsidy_b_grid= ones(no_taus);
APland_mean_log_subsidy_b_grid= ones(no_taus);
APland_mean_cashcrop_log_subsidy_b_grid= ones(no_taus);
APland_mean_staples_S_log_subsidy_b_grid= ones(no_taus);
var_APland_subsidy_b_grid= ones(no_taus);
var_APland_cashcrop_subsidy_b_grid= ones(no_taus);
var_APland_staples_S_subsidy_b_grid= ones(no_taus);
c_S_W_fine_subsidy_b_grid= ones(Baseline_parameter.ns_fine,3,no_taus);
c_B_W_fine_subsidy_b_grid= ones(Baseline_parameter.ns_fine,3,no_taus);
c_M_W_fine_subsidy_b_grid= ones(Baseline_parameter.ns_fine,3,no_taus);
c_S_S_fine_subsidy_b_grid= ones(Baseline_parameter.ns_fine,3,no_taus);
c_B_S_fine_subsidy_b_grid= ones(Baseline_parameter.ns_fine,3,no_taus);
c_M_S_fine_subsidy_b_grid= ones(Baseline_parameter.ns_fine,3,no_taus);
c_S_B_fine_subsidy_b_grid= ones(Baseline_parameter.ns_fine,3,no_taus);
c_B_B_fine_subsidy_b_grid= ones(Baseline_parameter.ns_fine,3,no_taus);
c_M_B_fine_subsidy_b_grid= ones(Baseline_parameter.ns_fine,3,no_taus);
iterate = 1
#τ_grid = zeros(no_taus);
#τ_grid[1:(no_taus - 1) ] = range(Baseline_parameter.τ_S,-0.05,no_taus-1)
for τ = τ_grid
    parameters_tmp = copy(Baseline_parameter);
    parameters_tmp.τ_S = τ;
    println("Populate for " ,τ, "subsidy rate")
    (residual_goods_subsidy_b_grid[:,iterate], stat_distr_subsidy_b_grid[:,iterate], cons_fine_local_subsidy_b_grid[:,:,iterate], a_prime_fine_local_subsidy_b_grid[:,:,iterate],
    future_occupation_fine_local_subsidy_b_grid[:,:,iterate],x_S_S_fine_subsidy_b_grid[:,:,iterate],x_SC_fine_subsidy_b_grid[:,:,iterate],x_BC_fine_subsidy_b_grid[:,:,iterate], 
    coeff_subsidy_b_grid[:,:,iterate],  transaction_cost_loss_subsidy_b_grid[iterate],nominal_GDP_subsidy_b_grid[iterate],welfare_val_subsidy_b_grid[:,iterate],
    Import_value_subsidy_b_grid[iterate],Export_value_subsidy_b_grid[iterate],current_worker_pop_subsidy_b_grid[iterate],current_staple_pop_subsidy_b_grid[iterate],
    current_cashcrop_pop_subsidy_b_grid[iterate], marketable_agr_surplus_share_subsidy_b_grid[iterate],exportshare_cashcrop_subsidy_b_grid[iterate],
    fraction_model_subsidy_b_grid[iterate],program_spending_subsidy_b_grid[iterate],prod_value_improvement_subsidy_b_grid[iterate],share_selling_increase_subsidy_b_grid[iterate],
    exp_ratio_model_subsidy_b_grid[iterate],mig_rate_model_subsidy_b_grid[iterate],rural_pop_only_staples_model_subsidy_b_grid[iterate],
    rural_pop_only_cashcrop_model_subsidy_b_grid[iterate],mean_land_share_to_staples_among_cc_model_subsidy_b_grid[iterate],urban_rural_inc_ratio_model_subsidy_b_grid[iterate],
    urban_rural_wealth_ratio_model_subsidy_b_grid[iterate],urban_rural_consumption_ratio_model_subsidy_b_grid[iterate],p90_wealth_rural_subsidy_b_grid[iterate],
    p90_wealth_urban_subsidy_b_grid[iterate],p99_wealth_rural_subsidy_b_grid[iterate],p99_wealth_urban_subsidy_b_grid[iterate],p90_cons_tmp_subsidy_b_grid[iterate],
    p90_income_tmp_subsidy_b_grid[iterate],p99_cons_tmp_subsidy_b_grid[iterate],p99_income_tmp_subsidy_b_grid[iterate],staple_productivity_subsidy_b_grid[iterate],
    cashcrop_productivity_subsidy_b_grid[iterate],manuf_productivity_subsidy_b_grid[iterate],relative_land_to_staples_subsidy_b_grid[iterate],
    relative_land_to_cashcrop_subsidy_b_grid[iterate],share_constrained_cashcrop_subsidy_b_grid[iterate],var_MPX_staples_S_subsidy_b_grid[iterate],
    var_MPX_cashcrop_B_subsidy_b_grid[iterate],var_MPX_cashcrop_S_subsidy_b_grid[iterate],share_constrained_staple_subsidy_b_grid[iterate],APG_subsidy_b_grid[iterate],
    urban_rural_consumption_ratio_model_real_subsidy_b_grid[iterate],aggregate_consumption_subsidy_b_grid[iterate],
    worker_pop_effective_subsidy_b_grid[iterate],prod_manuf_subsidy_b_grid[iterate],total_entry_cost_subsidy_b_grid[iterate],prod_staple_subsidy_b_grid[iterate],prod_cashcrop_subsidy_b_grid[iterate],input_staple_subsidy_b_grid[iterate],input_cashcrop_subsidy_b_grid[iterate],
    total_maintenance_cost_subsidy_b_grid[iterate],current_account_residual_subsidy_b_grid[iterate],fraction_cashcrop_suboptimal_model_subsidy_b_grid[iterate],V_saved_b_grid[:,:,iterate],
    avg_labor_prod_rural_subsidy_b_grid[iterate],avg_labor_prod_urban_subsidy_b_grid[iterate],avg_agri_prod_rural_subsidy_b_grid[iterate],avg_agri_prod_urban_subsidy_b_grid[iterate],
    var_MPX_cashcrop_subsidy_b_grid[iterate] ,var_MPX_subsidy_b_grid[iterate] ,TFP_subsidy_b_grid[iterate] ,YL_manuf_subsidy_b_grid[iterate]  , 
    YL_agr_subsidy_b_grid[iterate] ,coeff_var_labor_prod_rural_subsidy_b_grid[iterate] ,coeff_var_labor_prod_urban_subsidy_b_grid[iterate] ,
    coeff_var_agri_prod_rural_subsidy_b_grid[iterate] , coeff_var_agri_prod_urban_subsidy_b_grid[iterate] ,
    p90_wealth_subsidy_b_grid[iterate] ,p99_wealth_subsidy_b_grid[iterate] ,p90_cons_subsidy_b_grid_rural[iterate] ,p99_cons_subsidy_b_grid_rural[iterate] ,p90_cons_subsidy_b_grid_urban[iterate],p99_cons_subsidy_b_grid_urban[iterate],p90_income_subsidy_b_grid_rural[iterate] ,
    p99_income_subsidy_b_grid_rural[iterate],p90_income_subsidy_b_grid_urban[iterate],p99_income_subsidy_b_grid_urban[iterate],
    wealth_of_workers_subsidy_b_grid[iterate] ,wealth_of_staples_subsidy_b_grid[iterate] ,wealth_of_cashcrop_subsidy_b_grid[iterate] ,
    c_B_worker_sum_subsidy_b_grid[iterate],c_B_staple_sum_subsidy_b_grid[iterate],c_B_cashcrop_sum_subsidy_b_grid[iterate],c_S_worker_sum_subsidy_b_grid[iterate],c_S_staple_sum_subsidy_b_grid[iterate] ,c_S_cashcrop_sum_subsidy_b_grid[iterate],
    transaction_cost_staple_sum_subsidy_b_grid[iterate],transaction_cost_cashcrop_sum_subsidy_b_grid[iterate],transaction_cost_worker_sum_subsidy_b_grid[iterate],c_M_worker_sum_subsidy_b_grid[iterate],c_M_staple_sum_subsidy_b_grid[iterate],c_M_cashcrop_sum_subsidy_b_grid[iterate]
    ,MPX_mean_log_subsidy_b_grid[iterate], MPX_mean_staples_S_log_subsidy_b_grid[iterate],MPX_mean_cashcrop_log_subsidy_b_grid[iterate]
    , APland_mean_log_subsidy_b_grid[iterate],APland_mean_cashcrop_log_subsidy_b_grid[iterate], APland_mean_staples_S_log_subsidy_b_grid[iterate],var_APland_subsidy_b_grid[iterate],var_APland_cashcrop_subsidy_b_grid[iterate],var_APland_staples_S_subsidy_b_grid[iterate],
    c_S_W_fine_subsidy_b_grid[:,:,iterate],c_B_W_fine_subsidy_b_grid[:,:,iterate],c_M_W_fine_subsidy_b_grid[:,:,iterate],c_S_S_fine_subsidy_b_grid[:,:,iterate],c_B_S_fine_subsidy_b_grid[:,:,iterate],c_M_S_fine_subsidy_b_grid[:,:,iterate],c_S_B_fine_subsidy_b_grid[:,:,iterate],c_B_B_fine_subsidy_b_grid[:,:,iterate],c_M_B_fine_subsidy_b_grid[:,:,iterate]) = details_model(prices_subsidy_b_grid[:,iterate],parameters_tmp,2,moments,1.0,foreign_supply_capital_subsidy_b);
    iterate = iterate + 1;
end
for iterate = 1:no_taus
    welfare_subsidy_b_grid[iterate] = sum(stat_distr_subsidy_b_grid[:,no_taus] .* (exp.((welfare_val_subsidy_b_grid[:,iterate] - welfare_val_subsidy_b_grid[:,no_taus]) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
end

welfare_subsidy_b_grid_alt = copy(welfare_subsidy_b_grid);
for iterate = 1:no_taus
    welfare_subsidy_b_grid_alt[iterate] = (1 - Baseline_parameter.β)* (sum(stat_distr_subsidy_b_grid[:,iterate] .* welfare_val_subsidy_b_grid[:,iterate] ) - sum(stat_distr_subsidy_b_grid[:,no_taus].*welfare_val_subsidy_b_grid[:,no_taus])  ); # With balanced budget
end

welfare_subsidy_b_grid_real = copy(welfare_subsidy_b_grid);

V_saved_b_grid_reshaped=reshape(V_saved_b_grid,Baseline_parameter.ns_fine*3,no_taus)
for iterate = 1:no_taus-1
    welfare_subsidy_b_grid_real[iterate] = sum(stat_distr_subsidy_b_grid[:, no_taus] .* (exp.((V_saved_b_grid_reshaped[:, iterate] - V_saved_b_grid_reshaped[:, no_taus]) * (1.0 - Baseline_parameter.β)))) - 1 # With balanced budget
end

welfare_subsidy_b_grid_real_alt = copy(welfare_subsidy_b_grid_real);
for iterate = 1:no_taus-1
    welfare_subsidy_b_grid_real_alt[iterate] = (1 - Baseline_parameter.β) * (sum(stat_distr_subsidy_b_grid[:, iterate] .* V_saved_b_grid_reshaped[:, iterate]) - sum(stat_distr_subsidy_b_grid[:, no_taus] .* V_saved_b_grid_reshaped[:, no_taus]) ) # With balanced budget

end

welfare_subsidy_b_grid_alt_stat = copy(welfare_subsidy_b_grid);
for iterate = 1:no_taus
    welfare_subsidy_b_grid_alt_stat[iterate] = (1 - Baseline_parameter.β)* (sum(stat_distr_subsidy_b_grid[:,no_taus] .* welfare_val_subsidy_b_grid[:,iterate] ) - sum(stat_distr_subsidy_b_grid[:,no_taus].*welfare_val_subsidy_b_grid[:,no_taus])  ); # With balanced budget
end

welfare_subsidy_b_grid_real_alt_stat = copy(welfare_subsidy_b_grid_real);
for iterate = 1:no_taus-1
    welfare_subsidy_b_grid_real_alt_stat[iterate] = (1 - Baseline_parameter.β) * (sum(stat_distr_subsidy_b_grid[:, no_taus] .* V_saved_b_grid_reshaped[:, iterate]) - sum(stat_distr_subsidy_b_grid[:, no_taus] .* V_saved_b_grid_reshaped[:, no_taus]) ) # With balanced budget

end

welfare_subsidy_b_grid_chg_stat = copy(welfare_subsidy_b_grid_real);

for iterate = 1:no_taus
    welfare_subsidy_b_grid_chg_stat[iterate] = sum(stat_distr_subsidy_b_grid[:,iterate] .* (exp.((welfare_val_subsidy_b_grid[:,iterate] - welfare_val_subsidy_b_grid[:,no_taus]) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget
end
welfare_subsidy_b_grid_real_chg_stat = copy(welfare_subsidy_b_grid_real);
for iterate = 1:no_taus-1
    welfare_subsidy_b_grid_real_chg_stat[iterate] = sum(stat_distr_subsidy_b_grid[:, iterate] .* (exp.((V_saved_b_grid_reshaped[:, iterate] - V_saved_b_grid_reshaped[:, no_taus]) * (1.0 - Baseline_parameter.β)))) - 1 # With balanced budget
end

cons_difference_grid =  copy(welfare_subsidy_b_grid_real);
for iterate = 1:no_taus-1
    cons_difference_grid[iterate] =  sum(stat_distr_subsidy_b_grid[:, no_taus] .* reshape(cons_fine_local_subsidy_b_grid[:,:, iterate] ./ cons_fine_local_subsidy_b_grid[:, :,no_taus],Baseline_parameter.ns_fine*3 )) - 1 
end
cons_difference_grid_alt =  copy(welfare_subsidy_b_grid_real);
for iterate = 1:no_taus-1
    cons_difference_grid_alt[iterate] =   (sum(stat_distr_subsidy_b_grid[:,iterate] .*reshape(cons_fine_local_subsidy_b_grid[:,:, iterate] ,Baseline_parameter.ns_fine*3 ) ) -sum(stat_distr_subsidy_b_grid[:,no_taus].*reshape(cons_fine_local_subsidy_b_grid[:, :,no_taus],Baseline_parameter.ns_fine*3 )));
end


a_prime_difference_grid =  copy(welfare_subsidy_b_grid_real);
for iterate = 1:no_taus-1
    a_prime_difference_grid[iterate] =  sum(stat_distr_subsidy_b_grid[:, no_taus] .* reshape(a_prime_fine_local_subsidy_b_grid[:,:, iterate] ./ a_prime_fine_local_subsidy_b_grid[:, :,no_taus],Baseline_parameter.ns_fine*3 )) - 1 
end


# plot!(τ_grid,100*)


# plot!(τ_grid,100*welfare_subsidy_b_grid_alt)
# plot!(τ_grid,100*welfare_subsidy_b_grid_real_alt)

# plot(τ_grid,100*welfare_subsidy_b_grid_real)
# plot!(τ_grid,100*p90_cons_tmp_subsidy_b_grid)
# plot!(τ_grid,10*var_MPX_subsidy_b_grid)
# plot!(τ_grid,100*cons_difference_grid)
# plot!(τ_grid,100*avg_agri_prod_rural_subsidy_b_grid./(avg_agri_prod_rural_subsidy_b_grid[no_taus]).-100)
# plot!(τ_grid,100*avg_labor_prod_urban_subsidy_b_grid./(avg_labor_prod_urban_subsidy_b_grid[no_taus]).-100)
# plot!(τ_grid,APG_subsidy_b_grid)
# plot(τ_grid,100*wealth_of_workers_subsidy_b_grid)
# plot(τ_grid,100*wealth_of_staples_subsidy_b_grid)
# plot(τ_grid,100*fraction_model_subsidy_b_grid)
# plot(τ_grid,100*transaction_cost_loss_subsidy_b_grid)
# plot(τ_grid,100*a_prime_difference_grid)
# plot!(τ_grid,100*p90_wealth_subsidy_b_grid)
# plot(τ_grid,100*fraction_cashcrop_suboptimal_model_subsidy_b_grid)
# plot(τ_grid,[100*welfare_subsidy_b_grid,100*welfare_subsidy_b_grid_real,100*welfare_subsidy_b_grid_alt,100*welfare_subsidy_b_grid_real_alt]
# , label=["Laszlo expected" "Laszlo current" "Karol expected" "Karol current"],legend=:bottomright ,
# linewidth = 2,linestyle = [:solid :dash],ylims = [-10.0,5.0], xlabel = "τ_S",
# ylabel = "Welfare change",
# grid = false,
# tickfontsize = 14,xguidefontsize=12,yguidefontsize=12,legendfontsize=12,fontfamily="Times_Roman")
# savefig("Welfare_measures.pdf")

# plot(τ_grid,[100*movmean(welfare_subsidy_b_grid,4),100*movmean(welfare_subsidy_b_grid_real,4),100*movmean(welfare_subsidy_b_grid_alt,4),100*movmean(welfare_subsidy_b_grid_real_alt,4)]
# , label=["Laszlo expected" "Laszlo current" "Karol expected" "Karol current"],legend=:bottomright ,
# linewidth = 2,linestyle = [:solid :dash],ylims = [-10.0,5.0], xlabel = "τ_S",
# ylabel = "Welfare change",
# grid = false,
# tickfontsize = 14,xguidefontsize=12,yguidefontsize=12,legendfontsize=12,fontfamily="Times_Roman")
# savefig("Welfare_measures_smoothed.pdf")

# plot(τ_grid,[100*welfare_subsidy_b_grid,100*welfare_subsidy_b_grid_real,100*welfare_subsidy_b_grid_alt_stat,100*welfare_subsidy_b_grid_real_alt_stat]
# , label=["Laszlo expected" "Laszlo current" "Karol expected" "Karol current"],legend=:bottomright ,
# linewidth = 2,linestyle = [:solid :dash],ylims = [-10.0,5.0], xlabel = "τ_S",
# ylabel = "Welfare change",
# grid = false,
# tickfontsize = 14,xguidefontsize=12,yguidefontsize=12,legendfontsize=12,fontfamily="Times_Roman")
# savefig("Welfare_measures_stat.pdf")

# plot(τ_grid,[100*welfare_subsidy_b_grid_chg_stat,100*welfare_subsidy_b_grid_real_chg_stat,100*welfare_subsidy_b_grid_alt,100*welfare_subsidy_b_grid_real_alt]
# , label=["Laszlo expected" "Laszlo current" "Karol expected" "Karol current"],legend=:bottomright ,
# linewidth = 2,linestyle = [:solid :dash],ylims = [-10.0,5.0], xlabel = "τ_S",
# ylabel = "Welfare change",
# grid = false,
# tickfontsize = 14,xguidefontsize=12,yguidefontsize=12,legendfontsize=12,fontfamily="Times_Roman")
# savefig("Welfare_measures_chg_distr.pdf")

# plot(τ_grid,movmean(100*welfare_subsidy_b_grid_alt,4))
# plot!(τ_grid,100*welfare_subsidy_b_grid_alt)
# plot(τ_grid,100*cons_difference_grid)
# plot!(τ_grid,100*cons_difference_grid_alt)
# savefig("Figure_welfare.pdf")

# plot!(τ_grid,100* sum(abs.(residual_goods_subsidy_b_grid), dims = 1)')
# plot!(τ_grid,100* maximum(abs.(residual_goods_subsidy_b_grid), dims = 1)')
# plot(τ_grid,100*current_worker_pop_subsidy_b_grid)
# plot(τ_grid,10*prices_subsidy_b_grid[3,:])
# plot(τ_grid,100*prices_subsidy_b_grid[2,:])
# sum((stat_distr_sub_b_grid[:,iterate]*exp(welfare_val_sub_b_grid[:,iterate])) / sum(stat_distr_sub_b_grid[:,no_taus]*exp(welfare_val_sub_b_grid[:,no_taus]))) * (1-beta)
# plot(τ_grid,100* residual_goods_subsidy_b_grid[2,:])

# prices_subsidy_b_grid_smooth = copy(prices_subsidy_b_grid);
# prices_subsidy_b_grid_smooth[2,1:(end-1)] = movmean(prices_subsidy_b_grid[2,1:(end-1)],9);
# prices_subsidy_b_grid = copy(prices_subsidy_b_grid_smooth);
# # Consider people 
# sum(stat_distr_subsidy_b_grid[:,no_taus] .* (exp.((welfare_val_subsidy_b_grid[:,iterate] - welfare_val_subsidy_b_grid[:,no_taus]) * (1.0 - Baseline_parameter.β) ) ))  - 1; # With balanced budget

# 
#No QS + recalibrated F_W + no subsidy equilibrium
# prices_no_QS_recab_no_subsidy = [ 1.190752722671664,
# 0.16670260063915554]
# no_QS_recab_no_subsidy_parameter = copy(Baseline_parameter);
# no_QS_recab_no_subsidy_parameter.τ_S=-0;
# no_QS_recab_no_subsidy_parameter.Q_S=0;
# no_QS_recab_no_subsidy_parameter.F_W = 305.0

# (residual_goods_no_QS_recab_no_subsidy, stat_distr_no_QS_recab_no_subsidy, cons_fine_local_no_QS_recab_no_subsidy, a_prime_fine_local_no_QS_recab_no_subsidy,future_occupation_fine_local_no_QS_recab_no_subsidy,x_S_S_fine_no_QS_recab_no_subsidy,x_SC_fine_no_QS_recab_no_subsidy,x_BC_fine_no_QS_recab_no_subsidy, coeff_no_QS_recab_no_subsidy,
# transaction_cost_loss_no_QS_recab_no_subsidy,nominal_GDP_no_QS_recab_no_subsidy,welfare_val_no_QS_recab_no_subsidy,Import_value_no_QS_recab_no_subsidy,Export_value_no_QS_recab_no_subsidy,current_worker_pop_no_QS_recab_no_subsidy,current_staple_pop_no_QS_recab_no_subsidy,current_cashcrop_pop_no_QS_recab_no_subsidy,
# marketable_agr_surplus_share_no_QS_recab_no_subsidy,exportshare_cashcrop_no_QS_recab_no_subsidy,
# fraction_model_no_QS_recab_no_subsidy,program_spending_no_QS_recab_no_subsidy,prod_value_improvement_no_QS_recab_no_subsidy,share_selling_increase_no_QS_recab_no_subsidy,exp_ratio_model_no_QS_recab_no_subsidy,mig_rate_model_no_QS_recab_no_subsidy,rural_pop_only_staples_model_no_QS_recab_no_subsidy,rural_pop_only_cashcrop_model_no_QS_recab_no_subsidy,
# mean_land_share_to_staples_among_cc_model_no_QS_recab_no_subsidy,urban_rural_inc_ratio_model_no_QS_recab_no_subsidy,urban_rural_wealth_ratio_model_no_QS_recab_no_subsidy,urban_rural_consumption_ratio_model_no_QS_recab_no_subsidy,
# p90_wealth_rural_no_QS_recab_no_subsidy,p90_wealth_urban_no_QS_recab_no_subsidy,p99_wealth_rural_no_QS_recab_no_subsidy,p99_wealth_urban_no_QS_recab_no_subsidy,p90_cons_tmp_no_QS_recab_no_subsidy,p90_income_tmp_no_QS_recab_no_subsidy,p99_cons_tmp_no_QS_recab_no_subsidy,p99_income_tmp_no_QS_recab_no_subsidy,
# staple_productivity_no_QS_recab_no_subsidy,cashcrop_productivity_no_QS_recab_no_subsidy,manuf_productivity_no_QS_recab_no_subsidy,relative_land_to_staples_no_QS_recab_no_subsidy,relative_land_to_cashcrop_no_QS_recab_no_subsidy,share_constrained_cashcrop_no_QS_recab_no_subsidy,var_MPX_staples_S_no_QS_recab_no_subsidy,var_MPX_cashcrop_B_no_QS_recab_no_subsidy,var_MPX_cashcrop_S_no_QS_recab_no_subsidy,
# share_constrained_staple_no_QS_recab_no_subsidy,APG_no_QS_recab_no_subsidy,urban_rural_consumption_ratio_model_real_no_QS_recab_no_subsidy,aggregate_consumption_no_QS_recab_no_subsidy,
# worker_pop_effective_no_QS_recab_no_subsidy,prod_manuf_no_QS_recab_no_subsidy,total_entry_cost_no_QS_recab_no_subsidy,prod_staple_no_QS_recab_no_subsidy,prod_cashcrop_no_QS_recab_no_subsidy,input_staple_no_QS_recab_no_subsidy,input_cashcrop_no_QS_recab_no_subsidy,total_maintenance_cost_no_QS_recab_no_subsidy,
# current_account_residual_no_QS_recab_no_subsidy,fraction_cashcrop_suboptimal_model_no_QS_recab_no_subsidy,V_saved_no_QS_recab_no_subsidy
# ,avg_labor_prod_rural_no_QS_recab_no_subsidy,avg_labor_prod_urban_no_QS_recab_no_subsidy,avg_agri_prod_rural_no_QS_recab_no_subsidy,avg_agri_prod_urban_no_QS_recab_no_subsidy
# ,var_MPX_cashcrop_no_QS_recab_no_subsidy,var_MPX_no_QS_recab_no_subsidy,TFP_no_QS_recab_no_subsidy,YL_manuf_no_QS_recab_no_subsidy , YL_agr_no_QS_recab_no_subsidy,coeff_var_labor_prod_rural_no_QS_recab_no_subsidy,coeff_var_labor_prod_urban_no_QS_recab_no_subsidy,coeff_var_agri_prod_rural_no_QS_recab_no_subsidy, coeff_var_agri_prod_urban_no_QS_recab_no_subsidy
# ,
# p90_wealth_no_QS_recab_no_subsidy,p99_wealth_no_QS_recab_no_subsidy,p90_cons_no_QS_recab_no_subsidy_rural,p99_cons_no_QS_recab_no_subsidy_rural,p90_cons_no_QS_recab_no_subsidy_urban,p99_cons_no_QS_recab_no_subsidy_urban,p90_income_no_QS_recab_no_subsidy_rural,
# p99_income_no_QS_recab_no_subsidy_rural,p90_income_no_QS_recab_no_subsidy_urban,p99_income_no_QS_recab_no_subsidy_urban,
# wealth_of_workers_no_QS_recab_no_subsidy,wealth_of_staples_no_QS_recab_no_subsidy,wealth_of_cashcrop_no_QS_recab_no_subsidy,
# c_B_worker_sum_no_QS_recab_no_subsidy,c_B_staple_sum_no_QS_recab_no_subsidy,c_B_cashcrop_sum_no_QS_recab_no_subsidy,c_S_worker_sum_no_QS_recab_no_subsidy,c_S_staple_sum_no_QS_recab_no_subsidy ,c_S_cashcrop_sum_no_QS_recab_no_subsidy,
# transaction_cost_staple_sum_no_QS_recab_no_subsidy,transaction_cost_cashcrop_sum_no_QS_recab_no_subsidy,transaction_cost_worker_sum_no_QS_recab_no_subsidy,c_M_worker_sum_no_QS_recab_no_subsidy,c_M_staple_sum_no_QS_recab_no_subsidy,c_M_cashcrop_sum_no_QS_recab_no_subsidy
# )  = details_model(prices_no_QS_recab_no_subsidy,no_QS_recab_no_subsidy_parameter,2,moments,0.0,foreign_supply_capital_subsidy_b);


# #    Subsidy equilibrium with balanced budget + No QS + recalibrated F_W 
# prices_no_QS_recab_subsidy_b =  [1.5163159754906252,
# 0.24936579281247012,
# 0.1755991002351221]
# no_QS_recab_parameter = copy(Baseline_parameter);
# no_QS_recab_parameter.Q_S=0.0;
# no_QS_recab_parameter.F_W = 305.0
# (residual_goods_no_QS_recab_subsidy_b, stat_distr_no_QS_recab_subsidy_b, cons_fine_local_no_QS_recab_subsidy_b, a_prime_fine_local_no_QS_recab_subsidy_b,future_occupation_fine_local_no_QS_recab_subsidy_b,x_S_S_fine_no_QS_recab_subsidy_b,x_SC_fine_no_QS_recab_subsidy_b,x_BC_fine_no_QS_recab_subsidy_b, coeff_no_QS_recab_subsidy_b,
# transaction_cost_loss_no_QS_recab_subsidy_b,nominal_GDP_no_QS_recab_subsidy_b,welfare_val_no_QS_recab_subsidy_b,Import_value_no_QS_recab_subsidy_b,Export_value_no_QS_recab_subsidy_b,current_worker_pop_no_QS_recab_subsidy_b,current_staple_pop_no_QS_recab_subsidy_b,current_cashcrop_pop_no_QS_recab_subsidy_b,
# marketable_agr_surplus_share_no_QS_recab_subsidy_b,exportshare_cashcrop_no_QS_recab_subsidy_b,
# fraction_model_no_QS_recab_subsidy_b,program_spending_no_QS_recab_subsidy_b,prod_value_improvement_no_QS_recab_subsidy_b,share_selling_increase_no_QS_recab_subsidy_b,exp_ratio_model_no_QS_recab_subsidy_b,mig_rate_model_no_QS_recab_subsidy_b,rural_pop_only_staples_model_no_QS_recab_subsidy_b,rural_pop_only_cashcrop_model_no_QS_recab_subsidy_b,
# mean_land_share_to_staples_among_cc_model_no_QS_recab_subsidy_b,urban_rural_inc_ratio_model_no_QS_recab_subsidy_b,urban_rural_wealth_ratio_model_no_QS_recab_subsidy_b,urban_rural_consumption_ratio_model_no_QS_recab_subsidy_b,
# p90_wealth_rural_no_QS_recab_subsidy_b,p90_wealth_urban_no_QS_recab_subsidy_b,p99_wealth_rural_no_QS_recab_subsidy_b,p99_wealth_urban_no_QS_recab_subsidy_b,p90_cons_tmp_no_QS_recab_subsidy_b,p90_income_tmp_no_QS_recab_subsidy_b,p99_cons_tmp_no_QS_recab_subsidy_b,p99_income_tmp_no_QS_recab_subsidy_b,
# staple_productivity_no_QS_recab_subsidy_b,cashcrop_productivity_no_QS_recab_subsidy_b,manuf_productivity_no_QS_recab_subsidy_b,relative_land_to_staples_no_QS_recab_subsidy_b,relative_land_to_cashcrop_no_QS_recab_subsidy_b,share_constrained_cashcrop_no_QS_recab_subsidy_b,var_MPX_staples_S_no_QS_recab_subsidy_b,var_MPX_cashcrop_B_no_QS_recab_subsidy_b,var_MPX_cashcrop_S_no_QS_recab_subsidy_b,
# share_constrained_staple_no_QS_recab_subsidy_b,APG_no_QS_recab_subsidy_b,urban_rural_consumption_ratio_model_real_no_QS_recab_subsidy_b,aggregate_consumption_no_QS_recab_subsidy_b,
# worker_pop_effective_no_QS_recab_subsidy_b,prod_manuf_no_QS_recab_subsidy_b,total_entry_cost_no_QS_recab_subsidy_b,prod_staple_no_QS_recab_subsidy_b,prod_cashcrop_no_QS_recab_subsidy_b,input_staple_no_QS_recab_subsidy_b,input_cashcrop_no_QS_recab_subsidy_b,
# total_maintenance_cost_no_QS_recab_subsidy_b,current_account_residual_no_QS_recab_subsidy_b,fraction_cashcrop_suboptimal_model_no_QS_recab_subsidy_b,V_saved_no_QS_recab_subsidy_b
# ,avg_labor_prod_rural_no_QS_recab_subsidy_b,avg_labor_prod_urban_no_QS_recab_subsidy_b,avg_agri_prod_rural_no_QS_recab_subsidy_b,avg_agri_prod_urban_no_QS_recab_subsidy_b,
# var_MPX_cashcrop_no_QS_recab_subsidy_b ,var_MPX_no_QS_recab_subsidy_b ,TFP_no_QS_recab_subsidy_b ,YL_manuf_no_QS_recab_subsidy_b  , 
# YL_agr_no_QS_recab_subsidy_b ,coeff_var_labor_prod_rural_no_QS_recab_subsidy_b ,coeff_var_labor_prod_urban_no_QS_recab_subsidy_b ,
# coeff_var_agri_prod_rural_no_QS_recab_subsidy_b , coeff_var_agri_prod_urban_no_QS_recab_subsidy_b,
# p90_wealth_no_QS_recab_subsidy_b,p99_wealth_no_QS_recab_subsidy_b,p90_cons_no_QS_recab_subsidy_b_rural,p99_cons_no_QS_recab_subsidy_b_rural,p90_cons_no_QS_recab_subsidy_b_urban,p99_cons_no_QS_recab_subsidy_b_urban,p90_income_no_QS_recab_subsidy_b_rural,
# p99_income_no_QS_recab_subsidy_b_rural,p90_income_no_QS_recab_subsidy_b_urban,p99_income_no_QS_recab_subsidy_b_urban,
# wealth_of_workers_no_QS_recab_subsidy_b,wealth_of_staples_no_QS_recab_subsidy_b,wealth_of_cashcrop_no_QS_recab_subsidy_b,
# c_B_worker_sum_no_QS_recab_subsidy_b,c_B_staple_sum_no_QS_recab_subsidy_b,c_B_cashcrop_sum_no_QS_recab_subsidy_b,c_S_worker_sum_no_QS_recab_subsidy_b,c_S_staple_sum_no_QS_recab_subsidy_b ,c_S_cashcrop_sum_no_QS_recab_subsidy_b,
# transaction_cost_staple_sum_no_QS_recab_subsidy_b,transaction_cost_cashcrop_sum_no_QS_recab_subsidy_b,transaction_cost_worker_sum_no_QS_recab_subsidy_b,c_M_worker_sum_no_QS_recab_subsidy_b,c_M_staple_sum_no_QS_recab_subsidy_b,c_M_cashcrop_sum_no_QS_recab_subsidy_b
# )  = details_model(prices_no_QS_recab_subsidy_b,no_QS_recab_parameter,2,moments,1.0,foreign_supply_capital_subsidy_b);



##########################################################################################################################################################
##
##          ANALYSIS: additional data is constructed and load packages for the plotting
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

# Smoothing for plotting tauSgrid values
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
#a_prime_difference_grid_growth_smoothed = movmean((100*a_prime_difference_grid ),5)


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

APG_subsidy_b_grid
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

# Undernourishment requires additional calculations
# Staple consumption w a threshold value - assume the threshold is chosen such that the calibration hits the 47% of undernourished HHs
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


##########################################################################################################################################################
##
##         Plotting: 
##
##########################################################################################################################################################


#Figure 1 plot - lambda_S. requires to evaluate on within the solve_model_calibration.jl 
cons_span = 10:140
prices = copy(prices_subsidy_b);
parameters_tmp = copy(Baseline_parameter);
foreign_supply_capital = copy(foreign_supply_capital_subsidy_b);
out = 2; 

# To make this graph, evaluate the inside of solve_model_calibration2.jl:

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
println("Production value improvement to cash grant \\cite{daidone} & 11\$\\% \$ & ",convert(Int64, round(100 * prod_value_improvement_subsidy_b)),"\$\\% \$ \\tabularnewline")
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
println("Consumption equivalent welfare  & - &  +",round(100 * welfare_subsidy_partial_real_alt,digits = 1), "\\% &+",
round(100 * welfare_subsidy_nb_real_alt,digits = 1), "\\% &", round(100 * welfare_subsidy_b_real_alt,digits = 1) ,"\\% \\\\ ")
println("Consumption equivalent welfare fixed distribution  & - &  +",round(100 * welfare_subsidy_partial_real,digits = 1), "\\% &",
round(100 * welfare_subsidy_nb_real,digits = 1), "\\% &", round(100 * welfare_subsidy_b_real,digits = 1) ,"\\% \\\\ ")
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
ylabel = "Welfare change",xticks = ([21,66,86,109 ],["\$50","\$1000", "\$2500", "\$6000"]),yticks = ([-20,-10,0,10,20,30 ],["-20","-10","0","10","20","30"]),
grid = false,size = (800, 800),legendcolumns=2,marker = [:none :none :circle :none],
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure2a.svg")
#Farmers
y_value = (exp.((V_saved_subsidy_b[:,2] -V_saved_no_subsidy[:,2]) * (1.0 - Baseline_parameter.β) ) ) .- 1
y_mat = 100* mat_creator(y_value);
plot([y_mat[1,1:5:end],y_mat[3,1:5:end],y_mat[13,1:5:end],y_mat[16,1:5:end]], label=["Low rural & urban \\theta" "High rural, low urban \\theta" "Low rural, high urban \\theta" "High rural & urban \\theta"]
,legend=:outerbottom ,
linewidth = 2,linestyle = [:solid :dash :dot :dashdot],ylims = [-30.0,10.0], xlabel = "Wealth",
ylabel = "Welfare change",xticks = ([21,66,86,109 ],["\$50","\$1000", "\$2500", "\$6000"]),yticks = ([-30,-20,-10,0,10,20 ],["-30","-20","-10","0","10","20",]),
grid = false,size = (800, 800),legendcolumns=2,marker = [:none :none :circle :none],
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure2b.svg")


#Plotting of tauSgrid values
plot(-1*reverse(τ_grid),[reverse(aggregate_consumption_subsidy_b_grid_growth_smoothed)
 ,reverse(welfare_subsidy_b_grid_real_alt_smoothed), reverse(undernutritioned_subsidy_b_grid_growth_smoothed)]
, label=["Aggregate consumption"  "Welfare change" "Undernourished households"],legend=:outerbottom ,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-26.0,5.0], xlabel = "Subsidy rate per unit of staple input",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),yticks = ([-25,-20,-15,-10,-5,0,5 ],["-25","-20","-15","-10","-5","0","5"]),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure3a.svg")

plot(-1*reverse(τ_grid),[reverse(current_worker_pop_subsidy_b_grid_growth_smoothed), reverse(avg_agri_prod_rural_subsidy_b_grid_growth_smoothed),
reverse(avg_labor_prod_urban_subsidy_b_grid_growth_smoothed)]
, label=[ "Urbanization rate" "Mean rural ability" "Mean urban ability" ],legend=:outerbottom ,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-15.0,5.0], xlabel = "Subsidy rate per unit of staple input",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),yticks = ([-15,-10,-5,0,5 ],["-15","-10","-5","0","5"]),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure3b.svg")


plot(-1*reverse(τ_grid),[reverse(staple_productivity_subsidy_b_grid_smoothed),
reverse(var_MPX_subsidy_b_grid_growth_smoothed), reverse(APG_subsidy_b_grid_growth_smoothed)]
, label=[ "Staple productivity" "Dispersion of returns on fertilizer" "APG"  ],legend=:outerbottom,
linewidth = 2,linestyle = [:solid :dash :dashdot :dot],ylims = [-2.0,100.0], xlabel = "Subsidy rate per unit of staple input",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure3c.svg")

plot(-1*reverse(τ_grid),[reverse(c_S_worker_sum_subsidy_b_grid_growth_smoothed),
reverse(c_S_staple_sum_subsidy_b_grid_growth_smoothed), reverse(c_S_cashcrop_sum_subsidy_b_grid_growth_smoothed)]
, label=[ "Workers" "Staple farmers" "Cash crop farmers"  ],legend=:outerbottom,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-20.0,45.0], xlabel = "Subsidy rate per unit of staple input",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),yticks = ([-20,-10,0,10,20,30,40 ],["-20","-10","0","10","20","30","40"]),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure3d.svg")

plot(-1*reverse(τ_grid),[reverse(c_B_worker_sum_subsidy_b_grid_growth_smoothed),
reverse(c_B_staple_sum_subsidy_b_grid_growth_smoothed), reverse(c_B_cashcrop_sum_subsidy_b_grid_growth_smoothed)]
, label=[ "Workers" "Staple farmers" "Cash crop farmers"  ],legend=:outerbottom,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-20.0,45.0], xlabel = "Subsidy rate per unit of staple input",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),yticks = ([-20,-10,0,10,20,30,40 ],["-20","-10","0","10","20","30","40"]),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure3e.svg")

plot(-1*reverse(τ_grid),[reverse(c_M_worker_sum_subsidy_b_grid_growth_smoothed),
reverse(c_M_staple_sum_subsidy_b_grid_growth_smoothed), reverse(c_M_cashcrop_sum_subsidy_b_grid_growth_smoothed)]
, label=[ "Workers" "Staple farmers" "Cash crop farmers"  ],legend=:outerbottom,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [-20.0,45.0], xlabel = "Subsidy rate per unit of staple input",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),yticks = ([-20,-10,0,10,20,30,40 ],["-20","-10","0","10","20","30","40"]),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure3f.svg")

plot(-1*reverse(τ_grid),[reverse(food_share_worker_subsidy_b_grid_growth_smoothed),
reverse(food_share_staple_subsidy_b_grid_growth_smoothed), reverse(food_share_cashcrop_subsidy_b_grid_growth_smoothed)]
, label=[ "Workers" "Staple farmers" "Cash crop farmers"  ],legend=:outerbottom,
linewidth = 2,linestyle = [:solid :dot :dashdot],ylims = [0.0,40.0], xlabel = "Subsidy rate per unit of staple input",
ylabel = " Percent change compared to no subsidy",marker = [:none :circle :none],
grid = false,size = (800,800),
tickfontsize = 14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14,fontfamily="times")
savefig("Figure3g.svg")



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
println("Cons.-equiv. welfare   &  +",
    round(100 * welfare_subsidy_nb_real_alt, digits=1), "\\% &+",
    round(100 * welfare_inf_real_alt, digits=1), "\\% &+",
    round(100 * welfare_inf_sp_real_alt, digits=1), "\\% \\\\ ")

println("Cons.-equiv. welfare w/ fixed distribution  &  ",
    round(100 * welfare_subsidy_nb_real, digits=1), "\\% &+",
    round(100 * welfare_inf_real, digits=1), "\\% &+",
    round(100 * welfare_inf_sp_real, digits=1), "\\% \\\\ ")
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