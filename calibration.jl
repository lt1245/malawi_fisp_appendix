
#cd("/Users/paroofa/Dropbox/Malawi_project/Model/Model_code/transaction_cost_v6-Laszlo Copy/")
#cd("\\Users\\phbs\\Dropbox\\Malawi_project\\Model\\Model_code\\transaction_cost_v7\\")
cd("/Users/bpu063010/Dropbox/Malawi_project/Model/Model_code/transaction_cost_v7 -Laszlo Copy/")
#using MKL
using Roots
using MultistartOptimization, NLopt
using DelimitedFiles
using Optim
using DataFrames
using CSV

using SharedArrays, Distributed
using CompEcon, QuantEcon
using LinearAlgebra, Statistics
using SparseArrays, StructArrays
using BasisMatrices, Arpack
using NLsolve, BlackBoxOptim, LeastSquaresOptim
using BenchmarkTools
using LeastSquaresOptim

using GeometryTypes
using StatsBase, GLM


grid(ranges::NTuple{N, <: AbstractRange}) where N = Point.(Iterators.product(ranges...))

#choose if solve with balanced gov buget: balanced_share \in [0.0,1.0] governs how much of FISP expenditures is financed through labor tax.
balanced_share=0.0
balanced_share=min(max(balanced_share,0.0),1.0)
τ_W_guess0=0.5 #initial tau_w guess for calibrations with balanced_share>0.0

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
#include("details_model.jl")
#include("solve_model_nomanu_price.jl")

# Initialize the parameters
Simple_function_space = fundefn(:spli, 10, 0.1, 1.0,1);
Baseline_parameter = Parameter_type(0.1,-1.0,0.06,1/3,1/3,1/3,2.0,0.86,2.0,1/3,1/3,1/3,1.0,1.5,0.1,
0.3,0.0,0.2,0.0,0.0,0.0,0.05,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.85,0.92,0.0,0.045,[50,36],[200,36],
zeros(2),zeros(2),0.0005,10.0,1,Simple_function_space,Simple_function_space,
zeros(2,2),2,zeros(2,2),2,zeros(2),zeros(2),spzeros(2,2),spzeros(2,2),spzeros(2,2),spzeros(2,2),
spzeros(2,2),spzeros(2,2),spzeros(2,2),0.5,0.4,zeros(2),2,2,zeros(2),Simple_function_space);
Baseline_parameter.n = [25,16]
Baseline_parameter.n_fine = [25,16]
Baseline_parameter.no_labor_shocks=sqrt(Baseline_parameter.n[2])

if mod(Baseline_parameter.n[2],Baseline_parameter.no_labor_shocks)!=0
    error("n & no_labor_shocks not compatible")
end

Baseline_parameter.γ = 0.0;
if Baseline_parameter.γ == 0.0
    Baseline_parameter.a_min = 0.01
    #Baseline_parameter.a_min = 0.0
    Baseline_parameter.a_max = 30.0
end

curve = 0.2;
curve_fine  = 0.2;
Baseline_parameter.agrid = range(Baseline_parameter.a_min^curve,Baseline_parameter.a_max^curve,length=
    Baseline_parameter.n[1]).^(1/curve);


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
Baseline_parameter.b_D = -0.2;
Baseline_parameter.FM_W=0.0;
Baseline_parameter.F_B = 0.0;
Baseline_parameter.p_x = 1.26; # from multiple data sources (see DB folder Calibr/)
#Because i guess alternative way to approx p_x is from FOC: p_x=zeta*y/x. We found zeta=0.1 in the LSMS data already (btw this was only for fertilizer, maybe we want to add seed cost there too). Then, with maize price 0.23USD/kg and 35kg of fert applied per ha, empirically we’d have p_x=0.1*0.23*2200/35=1.44.
#If we also think about seeds, then we could use 0.15 estimated above (coz it seems to me like cost share estimation), and use p_x=0.15*0.23*2200/(35+25)=1.26.
Baseline_parameter.δ=0.05;
Baseline_parameter.ϕ_B = 1.0;
Baseline_parameter.ϕ_S = 1.0;



Baseline_parameter.agrid_fine = a_grid_fine_gen_midpoints(Baseline_parameter.agrid,
Baseline_parameter.a_min,Baseline_parameter.a_max,1,Baseline_parameter.n);
Baseline_parameter.n_fine[1] = size(Baseline_parameter.agrid_fine)[1];

Baseline_parameter.C_grid_fine = a_grid_fine_gen_midpoints(Baseline_parameter.agrid,
Baseline_parameter.a_min,Baseline_parameter.a_max,8,Baseline_parameter.n);
Baseline_parameter.C_grid_fine_no = size(Baseline_parameter.C_grid_fine)[1];
Baseline_parameter.fspace_C_fine = fundef((:spli, Baseline_parameter.C_grid_fine, 0,1));

Baseline_parameter.fspace_a = fundef((:spli, Baseline_parameter.agrid, 0,Baseline_parameter.spliorder
        ));# define function space for the approximation of the solution of vfi.
Baseline_parameter.fspace_a_fine = fundef((:spli, Baseline_parameter.agrid_fine, 0,1
        ));# define function space for the approximation of the stationary distribution.


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

# Internally calibrated parameters - starting values
Baseline_parameter.κ = 0.05; # working capital constraint
Baseline_parameter.a_D = 0.77;
Baseline_parameter.c̄_S = 0.0; # Will use an algorithm to set the maximal value 
Baseline_parameter.F_W = 85.0;
Baseline_parameter.FM_B = 1.72;#0.85;
Baseline_parameter.ρ = 0.780128279883382;
Baseline_parameter.ρ_SW = 0.18;
Baseline_parameter.A_W = 1.65;
Baseline_parameter.τ_S=-0.49 #effective subsidy rate from IHS2010
# Include tau_S 
param_start = [Baseline_parameter.a_D, Baseline_parameter.F_W, Baseline_parameter.FM_B, Baseline_parameter.κ, 
Baseline_parameter.ρ, Baseline_parameter.ρ_SW, Baseline_parameter.A_W,Baseline_parameter.τ_S];
param = copy(param_start);
param = [   1.319964393460525,
300.0,
  2.3199277084690983,
  0.05,
  0.7412814048806258,
  0.23999102235984843,
  2.85,
 -0.8099992851842415]

price_initial = [1.22, 0.215241619328014]

# Targets:
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
# Only use some of these moments for the calibration
parameters_tmp = copy(Baseline_parameter);
moments_lower_bound = [r, L, K_Y_ratio, Y_agr_Y_tot, 0.65, fraction_staple_producers_without_surplus, G_Y_ratio,
RCT_moment1_value, RCT_moment2_share, exp_ratio, 0.01, fraction_only_cashcrops, 0.25,
0.4, urban_rural_inc_ratio, urban_rural_wealth_ratio, urban_rural_consumption_ratio,
Top1_share_wealth_rural, Top1_share_income_rural, Top1_share_consumption_rural, Top10_share_wealth_rural, Top10_share_income_rural, Top10_share_consumption_rural,
Top1_share_wealth_urban, Top1_share_income_urban, Top1_share_consumption_urban, Top10_share_wealth_urban, Top10_share_income_urban, Top10_share_consumption_urban,
0.025, 0.025, 0.025,fraction_cashcrop_suboptimal,APG_data]
moments_upper_bound = [r, L, K_Y_ratio, Y_agr_Y_tot, 0.8, fraction_staple_producers_without_surplus, G_Y_ratio,
RCT_moment1_value, RCT_moment2_share, exp_ratio, 0.03, fraction_only_cashcrops, 0.35,
0.5, urban_rural_inc_ratio, urban_rural_wealth_ratio, urban_rural_consumption_ratio,
Top1_share_wealth_rural, Top1_share_income_rural, Top1_share_consumption_rural, Top10_share_wealth_rural, Top10_share_income_rural, Top10_share_consumption_rural,
Top1_share_wealth_urban, Top1_share_income_urban, Top1_share_consumption_urban, Top10_share_wealth_urban, Top10_share_income_urban, Top10_share_consumption_urban,
0.025, 0.025, 0.025,0.85,7.0]

residual_lower_tolerance = (moments_lower_bound - moments)./moments;
residual_upper_tolerance = (moments_upper_bound - moments)./moments;
residual_lower_tolerance = residual_lower_tolerance[2:34];
residual_upper_tolerance = residual_upper_tolerance[2:34];
function outer_f(param)
    println("Trying param : ", param)
    balanced_share= 1.0
    
    #param_start = [Baseline_parameter.a_D, Baseline_parameter.F_W, Baseline_parameter.FM_B, Baseline_parameter.κ, Baseline_parameter.ρ, Baseline_parameter.ρ_SW, Baseline_parameter.A_W]
    parameters_tmp.a_D = param[1];
    parameters_tmp.F_W = param[2];
    parameters_tmp.FM_B = param[3];
    parameters_tmp.κ= param[4];
    parameters_tmp.ρ = param[5];
    parameters_tmp.ρ_SW = param[6]
    parameters_tmp.A_W = param[7];
    parameters_tmp.τ_S= param[8];
    # Because of new shock parameters, we have to reconstruct the state space
    (parameters_tmp.s,parameters_tmp.ns,parameters_tmp.s_fine,
    parameters_tmp.ns_fine,parameters_tmp.Phi_z,parameters_tmp.Phi_z_fine,parameters_tmp.Phi,
    parameters_tmp.Phi_aug,parameters_tmp.P_kron,parameters_tmp.P_kron1,parameters_tmp.P_kron_fine,parameters_tmp.z,
    parameters_tmp.z_W) = setup_state_space(parameters_tmp);

    # Check that c̄_S is high enough for the no subsidy equilibrium
    parameters_tmp.c̄_S = cbar_creator(parameters_tmp);
    # resetting curr_util is not neessary
    function within_f(prices_goods)
        #println(" [+] Trying prices: ", prices_goods)
        residual = solve_model_calibration1(prices_goods,parameters_tmp,1,moments,balanced_share);
        return residual
    end
    #prices_goods = [1.2,1.48];
    function within_f_sq(prices_goods)
        residual = within_f(prices_goods);
        return sum(residual.^2)
    end
    #prices_goods = [  1.9296454285670754,
    #0.43389271866737805,
    #0.11553815945013067,
    #0.30175093866821556];
    ##prices_goods = [1.54,1.11,0.05,0.01];

    prices_goods = [1.5312302133415963,
     0.2461718505305291,
     0.21770529259768423];
    #prices = copy(prices_goods)
    #prices_goods = [3.9975550027818048,3.6679971269884234,0.05677843731745866]
    #cbar_violated=within_f_checkcbar(prices_goods) No longer required

    #if cbar_violated==0
    # Look for a reasonable initial guess in our database from our previous solutions
    if balanced_share>0.0
        dataframe_prices = CSV.read("database_prices_balanced.csv",DataFrame,delim = ',',header=true, missingstring="?",
                threaded=false, silencewarnings=true); # Sadly this line depends on the machine, if it produces an error,unhash the next line
    else
        dataframe_prices = CSV.read("database_prices.csv",DataFrame,delim = ',',header=true, missingstring="?",
                threaded=false, silencewarnings=true); # Sadly this line depends on the machine, if it produces an error,unhash the next line
    end
    #dataframe_prices = CSV.read("database_prices.csv",DataFrame,delim = 't'); #
    #println(typeof(dataframe_prices))
    #df1=convert(Matrix,dataframe_prices);
    df1=Matrix(dataframe_prices)
    memory_csv = min(10,size(df1,1)-1);
    price_initial_mat_full = df1[:,1:2];
    prices_goods_upper_bound = maximum(price_initial_mat_full.*(price_initial_mat_full.<9.95),dims= 1) .+ parameters_tmp.a_min;
    prices_goods_lower_bound = minimum(price_initial_mat_full,dims= 1).- parameters_tmp.a_min;
    price_initial_mat_grid = grid((
    range(prices_goods_lower_bound[1],prices_goods_upper_bound[1],length = 2),range(prices_goods_lower_bound[2],prices_goods_upper_bound[2],length = 2)));
    price_initial_mat = zeros(4 + memory_csv,2);
    for j = 1:4
        price_initial_mat[j,:]= convert(Array{Float64},price_initial_mat_grid[:,:,1][j])[1:2]
    end
    residual_initial_mat= zeros(4+ memory_csv,5);
    for i= 1:4
        residual_initial_mat[i,:] =  within_f(price_initial_mat[i,:]);
    end
    for i= 1:memory_csv
        price_initial_mat[i+4,:] .= df1[end-i,1:2];
        residual_initial_mat[i+4,:] =  within_f(price_initial_mat[i+4,:]);
    end

    sum_residual_initial = sum(residual_initial_mat.^2,dims  = 2);
    sum_residual_initial[ isnan.(sum_residual_initial)] .= Inf;
    minval, minindx = findmin(sum_residual_initial);
    prices_goods = price_initial_mat[minindx[1],:];
    println(prices_goods)
    price_guess_iterator= 0;
    while minval>2.0
        prices_goods_upper_bound = max.(prices_goods_lower_bound,prices_goods_upper_bound*0.95);
        prices_goods_lower_bound = min.(prices_goods_upper_bound,prices_goods_lower_bound*1.05);
        price_initial_mat_grid = grid((range(prices_goods_lower_bound[1],prices_goods_upper_bound[1],length = 2),range(prices_goods_lower_bound[2],prices_goods_upper_bound[2],length = 2)));
        price_initial_mat = zeros(4,2);
        for j = 1:4
            price_initial_mat[j,:]= convert(Array{Float64},price_initial_mat_grid[:,:,1][j])[1:2]
        end
        residual_initial_mat= zeros(4,5);
        for i= 1:4
            residual_initial_mat[i,:] =  within_f(price_initial_mat[i,:]);
        end

        sum_residual_initial = sum(residual_initial_mat.^2,dims  = 2);
        sum_residual_initial[ isnan.(sum_residual_initial)] .= Inf;
        minval, minindx = findmin(sum_residual_initial);
        price_guess_iterator = price_guess_iterator + 1;
        println(price_guess_iterator)
        if price_guess_iterator > 10
            minval = 1.9
        end
        prices_goods = price_initial_mat[minindx[1],:];
        println(prices_goods)
        #model_moments=2*9.99.*ones(25);
    end
    if minval < 1.9
        

        #prices_goods = copy(prices);
        lower_guess_bound = 0.1;
        upper_guess_bound = 5.0;
        #prices_goods_star = copy(prices_goods);
        #prices = copy(prices_goods);
        #prices_goods = copy(prices_goods_star);
        #prices_goods = [3.0,0.9];
        #prices_goods = [1.489,1.46];
        ls_res= LeastSquaresOptim.optimize(within_f, prices_goods, LevenbergMarquardt(),show_trace = true, store_trace = true,
        x_tol = 1e-9, f_tol= 1e-5,iterations=10,lower = lower_guess_bound * prices_goods,
        upper = upper_guess_bound * prices_goods);
        prices_goods = ls_res.minimizer;
        if ls_res.ssr>0.05
            stst_simplex = Optim.AffineSimplexer(0.025,0.01);
            optim_res = LeastSquaresOptim.optimize(within_f_sq, prices_goods, method =NelderMead(initial_simplex = stst_simplex),show_trace = true, store_trace= true,
            extended_trace=true,time_limit =10000.0);
            prices_goods = optim_res.minimizer;
        end
        ls_res= LeastSquaresOptim.optimize(within_f, prices_goods, LevenbergMarquardt(),show_trace = true, store_trace = true,
        x_tol = 1e-9, f_tol= 1e-5,iterations=10,lower = lower_guess_bound * prices_goods,
        upper = upper_guess_bound * prices_goods);
        prices_goods = ls_res.minimizer;
        prices_goods_star = copy(prices_goods);
        # Rootfinding
        #p_B_star = find_zero(within_f, (0.25,1.0), Bisection())
        #(residual, stat_distr, cons_fine_local, future_occupation_fine_local,x_SC_fine,
        #x_BC_fine, coeff,foreign_demand,residual_cashcrop) = solve_model_calibration(p_B_star,parameters_tmp,2,r,L);
        println("[ + ] Trying prices: ", prices_goods_star)
        (residual, stat_distr, cons_fine_local, future_occupation_fine_local,x_SC_fine,
        x_BC_fine, coeff,residual_goods,model_moments) = solve_model_calibration1(prices_goods_star,parameters_tmp,2,moments,balanced_share);
        #@btime residual = within_f(prices_goods);
    else
        residual = ones(33)*1002;
        model_moments=2*9.99.*ones(34);
        if balanced_share>0.0
            prices_goods_star=2*[9.99,9.99,9.99];
            residual_goods=2*[999.0,999.0,999.0,999.0];
        else
            prices_goods_star=2*[9.99,9.99];
            residual_goods=2*[999.0,999.0,999.0];
        end
    end
    #end

    # Save some stuff to compare it manually:
    save_parameters = transpose(param);
    save_small = zeros(10);
    save_moments = transpose(round.(model_moments; digits=2));
    residual1 = residual.* (1.0 .- ((residual_lower_tolerance.<residual) .* (residual_upper_tolerance.>residual)))
    # Residual correction block - only care for residuals outside the bounds, if any
    if balanced_share>0.0
        open("database_parameters_balanced.csv", "a") do io
                    writedlm(io, save_parameters,',')#writedlm(io, save_parameters)
        end
        save_small[1:3] = prices_goods_star;
        save_small[10] = sum(residual_goods.^2);
        save_small = transpose(save_small);
        open("database_prices_balanced.csv", "a") do io
                    writedlm(io, save_small,',')#writedlm(io, save_small)
        end

        open("database_moments_balanced.csv", "a") do io
                writedlm(io, save_moments,',')#writedlm(io, save_moments)
        end
    else
        open("database_parameters.csv", "a") do io
                    writedlm(io, save_parameters,',')#writedlm(io, save_parameters)
        end
        save_small[1:2] = prices_goods_star;
        save_small[10] = sum(residual_goods.^2);
        save_small = transpose(save_small);
        open("database_prices.csv", "a") do io
                    writedlm(io, save_small,',')#writedlm(io, save_small)
        end

        open("database_moments.csv", "a") do io
                writedlm(io, save_moments,',')#writedlm(io, save_moments)
        end
    end

    return residual1
end
function outer_f_sq(param)
    if abs(param[6])>0.47
        return 999.9
    elseif abs(param[5])>0.85
        return 999.9
    elseif param[1]<0.0
        return 999.9
    elseif param[2]<0.0
        return 999.9
    elseif param[3]<0.0
        return 999.9
    elseif param[5]<0.01
        return 999.9
    elseif param[4]<0.001
        return 999.9
    else
        return sum(outer_f(param).^2)
    end
end
# #residual = outer_f(param);

#check if ρ_SW feasible - keep urban first! ---> I found that rho_SW must be within [-0.2, 0.2], e.g. 0.21 doesnt work anymore:
#A_mat_tmp = [parameters_tmp.ρ_W 0.47; 0.47 parameters_tmp.ρ_S]^100

#param = [0.77,85.0,1.72,0.05,0.68,0.18,1.7]
lb_lim = [0.0, 0.0, 0.0, 0.05, 0.1, -0.47, -10.0,-1.0]
ub_lim = [100.0, 200.0, 10.0, 1.0, 0.89, 0.47, 10.0,0.0]


# # Normal calibration:
ls_res= LeastSquaresOptim.optimize(outer_f, param, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-4,lower = lb_lim,
upper =ub_lim);
param = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.5,0.3);
optim_res = LeastSquaresOptim.optimize(outer_f_sq, param, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit =50000.0);
param = optim_res.minimizer;

# # TIK TAK calibration:
# P = MinimizationProblem(outer_f_sq, lb_lim, ub_lim)
# local_method = NLoptLocalMethod(NLopt.LN_NELDERMEAD)
# multistart_method = TikTak(2)
# p = multistart_minimization(multistart_method, local_method, P)
# p.location, p.value

#Laszlo:
# parameters_tmp.ψ_M = 0.44
#Testing:
param[3] = 0.44
param[7] = 1.0
param[12] = 2.5
param[16] = 0.15;
param[1] = 0.88;
#wrong guess param:param = [0.9143978478773584,0.5059765625,0.45113281250000004,0.9345703125,0.8560546875,3.6888671875000005,0.7058203125,0.09296875,0.9172851562500002,2.2744140625,0.400390625,0.615234375,4.7373046875,0.7146875,0.6865625,0.13345703125,0.03623046875,0.06677734375]
outer_f(param)

outer_f_sq(param)
# even lower doesnt work
#using MultistartOptimization, NLopt
#       [1β,2ψ_S,3ψ_M,4a_D,5ϕ_B,6ϕ_S,7l_z_low,8c̄_S,9p_x,10κ,11F_B,12F_W,13FM_B,14ρ_W,15ρ_S,16σ_S,17ρ]
#lb_lim=[0.7,    0.02,0.02,0.1,0.1,0.1,0.01,0.0,0.05,1.0,0.0,0.0,0.0,0.5,0.5,0.01,0.01]
#ub_lim=[1/(1+r),0.8, 0.8, 2.5,3.0,3.0,0.95,0.4,1.4,2.5,5.0,5.0,1.5,0.98,0.98,0.5,0.8]
#P = MinimizationProblem(outer_f_sq, lb_lim, ub_lim)
#local_method = NLoptLocalMethod(NLopt.LN_NELDERMEAD)
#multistart_method = TikTak(300)
#p = multistart_minimization(multistart_method, local_method, P)
#p.location, p.value
r,L,K_Y_ratio,Y_agr_Y_tot,exported_cash_crop_ratio,fraction_staple_producers_without_surplus,G_Y_ratio,
        RCT_moment1_value,RCT_moment2_share,exp_ratio,RU_migration_rate,fraction_only_cashcrops,mean_land_share_to_staples_among_cc,
        rural_pop_only_staples,urban_rural_inc_ratio,urban_rural_wealth_ratio,urban_rural_consumption_ratio,
        Top10_share_consumption_rural,Top10_share_consumption_urban,Top10_share_income_rural,Top10_share_income_urban,Top10_share_wealth_rural,Top10_share_wealth_urban,
        Top1_share_consumption_rural,Top1_share_consumption_urban,Top1_share_income_rural,Top1_share_income_urban,Top1_share_wealth_rural,Top1_share_wealth_urban,
        Top1_share_income,Top1_share_consumption,Top10_share_income,Top10_share_consumption

0.0,0.05,0.3 ,0.6 ,0.9 ,0.03,0.11,0.072,2.0,0.008,0.41,0.06,0.3 ,2.44,3.03,2.24,0.49,0.73,0.17,0.35,0.34,0.48,0.06,0.18

0.0,0.19,0.14,0.53,0.03,0.08,0.0 ,0.3  ,0.66,0.1  ,0.08,0.0 ,0.62,0.62,1.39,0.8 ,0.25,0.28,0.05,0.04,0.18,0.2 ,0.03,0.03


#Get the actually converted parameters and prices:
#K_a*r^K_b
index_final_param = 550;
dataframe_prices = CSV.read("database_prices.csv",DataFrame,delim = ',');
dataframe_parameters = CSV.read("database_parameters.csv",DataFrame,delim = ',');
dataframe_moments = CSV.read("database_moments.csv",DataFrame,delim = ',');
#dataframe_prices = CSV.read("database_prices.csv",DataFrame,delim = 't'); #
#println(typeof(dataframe_prices))
#df1=convert(Matrix,dataframe_prices);
df_prices=Matrix(dataframe_prices);
df_parameters=Matrix(dataframe_parameters);
df_moments=Matrix(dataframe_moments);
prices_goods = df_prices[index_final_param,1:2];
param = df_parameters[index_final_param,:];
moments = df_moments[index_final_param,:];
# Give values, modify parameters
r = moments[1];
K_Y_ratio = moments[3];

parameters_tmp.β = param[1];
#parameters_tmp.σ = param[2];
parameters_tmp.ψ_S = param[3-1];
parameters_tmp.ψ_M = param[4-1];
parameters_tmp.a_D = param[5-1];
parameters_tmp.ϕ_B = param[6-1];
parameters_tmp.ϕ_S = param[7-1];
parameters_tmp.l_z_low = param[8-1];
parameters_tmp.c̄_S = param[9-1];
parameters_tmp.p_x = param[10-1];
parameters_tmp.κ = param[11-1];
parameters_tmp.F_B = param[12-1]; #with F_MB=0; or FM_B=param[10] instead, and then set F_B=0
parameters_tmp.F_W = param[13-1]; #(or FM_W)
parameters_tmp.FM_B = param[14-1]; #with F_MB=0; or FM_B=param[10] instead, and then set F_B=0
#parameters_tmp.FM_W = param[15]; #(or FM_W)
parameters_tmp.ρ_W = param[15-1];
parameters_tmp.ρ_S = param[16-1];
parameters_tmp.σ_S= param[17-1];
parameters_tmp.ρ = param[17];
foreign_supply_capital =param[18];
# Construct other parameters:
parameters_tmp.δ = parameters_tmp.α/K_Y_ratio - r# Using K/Y, this can be evaluated even before equilibrium
parameters_tmp.ψ_B = 1.0 - parameters_tmp.ψ_M - parameters_tmp.ψ_S;

# Because of new shock parameters, we have to reconstruct the state space
(parameters_tmp.s,parameters_tmp.ns,parameters_tmp.s_fine,
parameters_tmp.ns_fine,parameters_tmp.Phi_z,parameters_tmp.Phi_z_fine,parameters_tmp.Phi,
parameters_tmp.Phi_aug,parameters_tmp.P_kron,parameters_tmp.P_kron1,parameters_tmp.P_kron_fine,parameters_tmp.z,
parameters_tmp.l_z) = setup_state_space(parameters_tmp);



foreign_supply_capital = K_a*r^K_b;

####################################
### Post-calibration validation:
####################################

#example:
param=[0.738744,0.337637,0.398574,0.453906,1.0,1.0,0.381777,0.301172,1.27739,1.85693,0.5,0.5,0.695801,0.507031,0.879219,0.0927832,0.626416,0.1]
prices=[1.5, 0.9]

(residual, stat_distr, cons_fine_local, future_occupation_fine_local,x_SC_fine,x_BC_fine, coeff,residual_goods,model_moments,
land_C_fine,q_S_S_fine,q_S_B_fine,q_B_C_fine,x_S_S_fine,x_SC_fine,x_BC_fine) = solve_model_calibration1(prices,parameters_tmp,2,moments,foreign_supply_capital,balanced_share);

ns_fine= Baseline_parameter.n_fine[1]*Baseline_parameter.n_fine[2]
worker_past_dist = stat_distr[(ns_fine *0 + 1):(ns_fine *1)];
staple_past_dist = stat_distr[(ns_fine *1 + 1):(ns_fine *2)];
cash_crop_past_dist = stat_distr[(ns_fine *2 + 1):(ns_fine *3)];
stat_distr_mat=[worker_past_dist staple_past_dist cash_crop_past_dist]

const prob_dist = vec(stat_distr_mat .* (future_occupation_fine_local .!= 1)) ./ sum(stat_distr_mat .* (future_occupation_fine_local .!= 1))
weights = Weights(prob_dist)
ShareMaize = zeros(3*ns_fine);
ValueProd = zeros(3*ns_fine);
ShareMaize =  vec((1 .- land_C_fine) .* (future_occupation_fine_local .== 3) .+ 1 .* (future_occupation_fine_local .== 2))
ValueProd = vec((1 .* q_S_S_fine .* (future_occupation_fine_local .== 2)) .+ ((prices[1] .* q_B_C_fine .+ 1 .* q_S_B_fine) .* (future_occupation_fine_local .== 3)))
SSC_term1=(c_S_S_fine .* (future_occupation_fine_local .== 2)) ./ (1 .* q_S_S_fine .* (future_occupation_fine_local .== 2));
replace!(SSC_term1, NaN=>0)
SSC_term2=(prices[1] .* c_C_C_fine .* (future_occupation_fine_local .== 3) .+ c_S_C_fine .* (future_occupation_fine_local .== 3)) ./ ((prices[1] .* q_B_C_fine .+ 1 .* q_S_B_fine) .* (future_occupation_fine_local .== 3));
replace!(SSC_term2, NaN=>0)
ShareSelfConsumed = vec(SSC_term1 .+ SSC_term2); #maybe: max.(0,min.(1.0,vec(SSC_term1 .+ SSC_term2)));
InputsUsed = vec((x_S_S_fine .* (future_occupation_fine_local .== 2)) .+ (( x_BC_fine .+ x_SC_fine) .* (future_occupation_fine_local .== 3)))

sampled_data=zeros(1000000,5);
for ii=1:1000000
    sampled_data[ii,1]=sample(ShareMaize,weights)
    sampled_data[ii,2]=sample(ValueProd,weights)
    sampled_data[ii,3]=sample(ShareSelfConsumed,weights)
    sampled_data[ii,4]=sample(InputsUsed,weights)
    sampled_data[ii,5]=sample(repeat(kron(z,ones(n_fine[1])),size(l_z)[1]*3),weights) #rural productivy
end

df_validate = DataFrame(share_maize=sampled_data[:,1], prod = sampled_data[:,5], value_prod = log.(sampled_data[:,2]), share_selfconsumed = sampled_data[:,3], inputs_used = log.(sampled_data[:,4]))

Validate_reg1 = coef(lm(@formula(value_prod ~ share_maize), df_validate)) #0.08 in the data
Validate_reg2 = coef(lm(@formula(inputs_used ~ share_maize), df_validate)) #0.49 in the data
Validate_reg3 = coef(lm(@formula(share_selfconsumed ~ share_maize), df_validate)) #24 in the data  --- this one is way off! :/

Validate_reg1 = coef(lm(@formula(value_prod ~ share_maize + prod), df_validate)) #0.08 in the data
Validate_reg2 = coef(lm(@formula(inputs_used ~ share_maize + prod), df_validate)) #0.49 in the data
Validate_reg3 = coef(lm(@formula(share_selfconsumed ~ share_maize + prod), df_validate)) #24 in the data  --- this one is way off! :/
