#cd("/Users/laszl/Dropbox/Malawi_project/Model/Model_code/github_version/")
#cd("/Users/paroofa/Dropbox/Malawi_project/Model/Model_code/github_version/")
cd("/Users/bpu063010/Dropbox/Malawi_project/Model/Model_code/github_version/")
#cd("\\Users\\phbs\\Dropbox\\Malawi_project\\Model\\Model_code\\github_version\\")

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
using StatsBase
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
include("solve_model_calibration_externality.jl")
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
APG_data= 6.3 #[4,7]

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
