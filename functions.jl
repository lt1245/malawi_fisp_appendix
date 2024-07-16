# START AR(1) CODE FROM gospodinov & Lkhagvasuren (2014 JAE)

# (P_Rouw, z_Rouw) = setup_rouwen(parameters_tmp.ρ_S, 0,
#     parameters_tmp.σ_S / sqrt(1 - parameters_tmp.ρ_S^2), no_prod_shocks);
# exp_egrid = exp.(z_Rouw);

# #Productivity shock - urban
# no_urban_prod_shocks = parameters_tmp.no_labor_shocks;
# (P_urban_Rouw, z_urban_Rouw) = setup_rouwen(parameters_tmp.ρ_W, parameters_tmp.A_W,
#     parameters_tmp.σ_W / sqrt(1 - parameters_tmp.ρ_W^2), no_urban_prod_shocks);

# % This function constructs a finite state Markov chain approximation 
# % using the MM method in Gospodinov and Lkhagvasuren, 2013, for the two
# % dimensional VAR(1) process considered in the numerical experiment
# % of the  paper.
# %
# % The VAR(1) process: y'=Ay+epsilon,
# % where var(epsilon) is given by a daigonal matrix Omega.
# %
# %   INPUT: 
# %          A0x stands for the 2X2 coefficient matrix A.
# %          vex is the 2X2 diagonal matrix, Omega, i.e.  
# %               Omega(1,1)=omega_{1,1}^2 and Omega(2,2)=omega_{2,2}^2
# %          nbar  is the number of grid points for each i.
# %          ntune is the control variable, where 
# %                 setting ntune=0 performs the baseline method (MM0), while
# %                 setting ntune>1 performs the full version of the method,
# %                 MM. For the examples considered in the paper, ntune was 
# %                 set to 1000. While higher values of ntune gives a better 
# %                 approximation, the gain becomes negligible beyond 
# %                 the value 10000. 
# % 
# %          
# %
# %   OUPUT: 
# %           PN is the N*-by-N* transition matrix, where N* = nbar^2. The
# %                [row k, column j] element is the probability the system 
# %                switches from state j to state k. So, the elements of each 
# %                column add up to 1.
# %           YN is the N*-by-2 matrix of the discrete values of y1 and y2 
# %                for N* states. 
# %          
function var_Markov_MM(A0x::Array{Float64,2}, vex::Array{Float64,2}, nbar::Int64, ntune::Int64)

    if ntune<0
        println("ntune has to be a positive integer (including zero)");
    end

    if mod(ntune,1)!=0 
        println("ntune has to be a positive integer");
    end

    nx=ntune+1;
    n=nbar;
    n1=n;
    n2=n;

    probtemp,z = setup_rouwen(0,0,1,n);  
    y1=z;
    y2=y1;

    A0=A0x;

    #normalize the initial var so unconditional variances are 1.
    A0new, vynew, vyold, venew = var_norm(A0, vex);  

    vy=vyold;
    A=A0new;
    ve=venew;

    pmat=zeros(2,n,n,n);
    px=zeros(2,n,n);

    for i=1:n
        for j=1:n
            for k=1:2
                
                mu=A[k,1]*y1[i]+ A[k,2]*y2[j];
                
                vact=ve[k,k];
                
                r=sqrt(1-vact);
                
                prob1, z = setup_rouwen(r,0,1,n);
                
                v1, p, na,nb, dummy_exceed =cal_mu_fast(mu,vact,n,z);
                
                if nx<2
                    
                    if na==nb  #% if mu is outside of the grids
                        pmat[k,i,j,:]=prob1[:,na];
                    else #% more relevant case
                        pmat[k,i,j,:]=p*prob1[:,na]+(1-p)*prob1[:,nb];
                    end
                                    
                else
                    
                    if na==nb  #% if mu is outside of the grids
                        pmat[k,i,j,:]=prob1[:,na];
                    else # % begining of the more relevane
                        
                        B=999*ones(nx,6);
                        ixx=0;
                        for ix=1:nx
                            vactx=max(0.00000000000001, vact*(1.0-(ix-1)/(nx-1)));
                            v1x, px, nax,nbx, dummy_exceedx=cal_mu_fast(mu,vactx,n,z);
                            if abs(dummy_exceedx)<0.5
                                ixx=ixx+1;
                                B[ixx,:]=[v1x px nax nbx dummy_exceedx vactx];
                            end
                        end
                        
                        
                        if ixx<1
                            pmat[k,i,j,:]=p*prob1[:,na]+(1-p)*prob1[:,nb];
                        else
                            bvectemp=B[:,1] .- vact;
                            dif1=abs.(bvectemp);
                            difx, iz=findmin(dif1);
                            
                            pz=B[iz,2];
                            naz=Int64(B[iz,3]);
                            nbz=Int64(B[iz,4]);
                            vactz=B[iz,6];
                            
                            rz=sqrt(1-vactz);
                            probz, z = setup_rouwen(rz,0,1,n);
                            pmat[k,i,j,:]=pz*probz[:,naz]+(1-pz)*probz[:,nbz];
                            
                        end  # end of the more relevane
                    end
                end
                
            end
        end
    end


    #% convert the transition probabilities into a conventional form
    PN = bigPPP(pmat,n);

    ynum=n*n;
    ix=0;
    Y=zeros(n*n,2);
    for i=1:n
        for j=1:n
            ix=ix+1;
            Y[ix,:]=[y1[i] y2[j]];
        end
    end

    YN=[Y[:,1]*sqrt(vy[1,1]) Y[:,2]*sqrt(vy[2,2])];

    return PN, YN
end

# %var_norm
# %   This code is written to normalize unconditional variance of components
# %   of a VAR(1) process: y'=Ay+epsilon.
# %
# %   INPUT: ve - covariance matrix of the error term. This is a diagonal
# %               matrix where the i-th diagonal element is var(epsilon_i).
# %          A  - the coef. matrix.
# %
# %   OUPUT: vynew - cov. matrix of normalized y
# %          vyold - initial cov. matrix of y
# %          venew - cov. matrix of the error term of the normalized process
# %          Anew   - the new coef. matrix (i.e. the matrix of the normalized 
# %          process)
# %  Damba Lkhagvasuren, 2009

function var_norm(A::Array{Float64,2}, ve::Array{Float64,2})

    dif=100;
    temp=size(A);
    nx=temp[1];
    V0=zeros(nx,nx);

    while dif>0.00000000001
        V=A*V0*A'+ve;
        dif=maximum(maximum(V .- V0));
        V0=V;
    end
    
    vyold=V0;

    venew=zeros(nx,nx);
    Anew=zeros(nx,nx);
    for i=1:nx
        venew[i,i]=ve[i,i]/vyold[i,i];
        for j=1:nx
            Anew[i,j]=A[i,j]*sqrt(vyold[j,j])/sqrt(vyold[i,i]);
        end
    end

    vynew=zeros(nx,nx);

    for i=1:nx
        for j=1:nx
            vynew[i,j]=vyold[i,j]/(sqrt(vyold[i,i])*sqrt(vyold[j,j]) );
        end
    end
    return Anew, vynew, vyold, venew
end

# % cal_mu_fast
# % This function calculates the conditional variance of the mixture
# % distribution given the conditional mean mu and the conditional variance 
# % v0 of the mass distributions on the n grids given by z.
# % 
# % For details, see Nikolay Gospodinov and Damba Lkhagvasuren, 2013

function cal_mu_fast(mu::Float64,v0::Float64,n::Int64,z::Array{Float64,1})

    r=sqrt(1-v0);

    zm=z*r;

    if mu>=zm[n]
        dummy_exceed=1;
        na=n;
        nb=n;
        p=0;
        v1=v0;
    elseif mu<=zm[1]
        dummy_exceed=-1;
        na=1;
        nb=1;
        p=1;
        v1=v0;
    else 
        dummy_exceed=0;
        
        #if isnan(floor((mu-zm[1])/(zm[2]-zm[1])))==1
        #    na=1
        #else
            na=1+floor((mu-zm[1])/(zm[2]-zm[1]));
        #end
        
        nb=na+1; 
        zax = zm[Int64(na)];
        zbx = zm[Int64(nb)];

        p=(zbx-mu)/(zbx-zax);

        v1=v0+p*(1-p)*(zbx-zax)^2;
        
    end
    return v1, p, na,nb, dummy_exceed
end

# % This function is used by the main code which constructs the 
# % Finite-state Markov chain using the MM method.
# % For details, see Nikolay Gospodinov and Damba Lkhagvasuren, 2013 
function bigPPP(pmatxxx::Array{Float64,4}, n::Int64)
    PPP=zeros(n^2,n^2);

        ix2=0;
        for i1=1:n
            for i2=1:n
                ix2=ix2+1;
                for i3=1:n
                    for i4=1:n
                        ix1=(i3-1)*n+i4;
                        PPP[ix1,ix2]=pmatxxx[1,i1,i2,i3]*pmatxxx[2,i1,i2,i4];
                    end
                end
            end
        end
        
        for i = 1:n*n
            PPP[:,i] = PPP[:,i] / sum(PPP[:,i]);
        end

    return PPP
end

function setup_rouwen(rho_Rouw::Number, mu_uncond::Number, sig_uncond::Number, n_R::Int)
    step_R = sig_uncond*sqrt(n_R - 1)
    z_Rouw = Array(-1:2/(n_R-1):1);
    z_Rouw = mu_uncond*ones(size(z_Rouw))+step_R.*z_Rouw;
    p=(rho_Rouw + 1)/2;
    q=p;
    P_Rouw=[ p  (1-p);
            (1-q) q];
    for i_R=2:n_R-1
        a1R=[P_Rouw zeros(i_R, 1); zeros(1, i_R+1)];
        a2R=[zeros(i_R, 1) P_Rouw; zeros(1, i_R+1)];
        a3R=[zeros(1,i_R+1); P_Rouw zeros(i_R,1)];
        a4R=[zeros(1,i_R+1); zeros(i_R,1) P_Rouw];
        P_Rouw=p*a1R+(1-p)*a2R+(1-q)*a3R+q*a4R;
        P_Rouw[2:i_R, :] = P_Rouw[2:i_R, :]/2;
    end
    P_Rouw=P_Rouw';
    for i_R = 1:n_R
        P_Rouw[:,i_R] = P_Rouw[:,i_R]/sum(P_Rouw[:,i_R]);
    end
    return (P_Rouw, z_Rouw)
end

# END AR(1) CODE FROM gospodinov & Lkhagvasuren (2014 JAE)


function goldenx(f::Function,a::Array{Float64,1},b::Array{Float64,1},tol::Float64=1e-10,arg...)
    """Vectorized golden section search to maximize univariate functions simultaneously
    Returns the maximum and the maximal value of f. Closer to Matlab code
    Parameters
    ----------
    f : function
        Function that maps from R^n to R^n such that it is an augmented univariate function
    a : array_like
        The lower bound for the maximization, for each dimension
    b : array_like
        The upper bound for the maximization, for each dimension
    tol : float
        Specifies the default tolerance value (for the argument of f) - default is 10**(-10)
    Returns
    -------
    x1 : ndarray
       Returns the n dimensional solution of the maximization.
    f1 : ndarray
       Returns the n dimensional maximum values.
    Notes
    -----
    To test in Julia:
    f(x::Array{Float64,1},c::Array{Float64,1},d::Float64) = -d*(x-c).^2;
    c = Array(range(1,10,length=10));
    a = zeros(10);
    b = ones(10) * 20;
    tol = 1e-10;
    d= 2.0;
    x, val = goldenx(f,a,b,tol,c,d);
    """
    alpha1 = (3.0 - sqrt(5)) / 2.0;
    alpha2 = 1.0 - alpha1;
    d  = b - a;
    x1 = a + alpha1 * d;
    x2 = a + alpha2 * d;
    n_tmp = size(x1);
    ones_tmp = ones(n_tmp);
    sign_tmp  = copy(ones_tmp);
    f1 = f(x1,arg...);
    f2 = f(x2,arg...);
    d = alpha1 * alpha2 * d;
    conv = 2.0;
    while conv > tol
        i = f2.>f1;
        not_i = ones_tmp - i;
        x1[i] = x2[i];
        f1[i] = f2[i];
        d = alpha2 * d;
        x2 = x1 + sign_tmp.*(i - not_i).*d;
        sign_tmp = sign.(x2-x1);
        f2 = f(x2,arg...);
        conv = maximum(d);
    end
    return x1, f1
end


function goldenx_2d(f::Function,a::Array{Float64,2},b::Array{Float64,2},tol::Float64=1e-10,arg...)
    """Vectorized golden section search to maximize bivariate functions simultaneously
    Returns the maximum and the maximal value of f.
    Experimental, based on G Sandhya Rani et al 2019 IOP Conf. Ser.: Mater. Sci. Eng. 577 012175
    Parameters
    ----------
    f : function
        Function that maps from R^n to R^n such that it is an augmented univariate function
    a : array_like
        The lower bound for the maximization, for each dimension
    b : array_like
        The upper bound for the maximization, for each dimension
    tol : float
        Specifies the default tolerance value (for the argument of f) - default is 10**(-10)
    Returns
    -------
    x1 : ndarray
       Returns the n dimensional solution of the maximization.
    f1 : ndarray
       Returns the n dimensional maximum values.
    Notes
    -----
    To test in Julia:
    f(x::Array{Float64,2},c1::Array{Float64,1},d1::Float64) = -d1*(x[:,1]-c1).^2 + -(x[:,2]-c1./2).^2;
    c1 = Array(range(1,10,length=1000));
    a = zeros(1000,2);
    b = ones(1000,2) * 20;
    x = ones(1000,2);
    d1= 2.0;
    f(x,c1,d1)
    tol = 1e-10;
    goldenx_2d(f,a,b,tol,c1,d1)
    """

    alpha1 = (3.0 - sqrt(5)) / 2.0;
    alpha2 = 1.0 - alpha1;
    d  = b - a;
    x1 = a + alpha1 * d;
    x4 = a + alpha2 * d;
    x2 = copy(x1);
    x3 = copy(x4);
    x3[:,1] = x1[:,1];
    x2[:,1] = x4[:,1];
    n_tmp = size(x1);
    ones_tmp = ones(n_tmp);
    f_mat  = zeros(n_tmp[1],4);
    f_mat[:,1] = f(x1,arg...);
    f_mat[:,2] = f(x2,arg...);
    f_mat[:,3] = f(x3,arg...);
    f_mat[:,4] = f(x4,arg...);

    d = alpha1 * alpha2 * d;
    conv = 2.0;
    while conv > tol
        f_tmp, i = findmax(f_mat; dims=2);
        index_max = last.(Tuple.(i));
        drop_index = 5 .- index_max;
        drop_index1 = (drop_index.==1)
        drop_index2 = (drop_index.==2)
        drop_index3 = (drop_index.==3)
        drop_index4 = (drop_index.==4)
        a[:,1] = drop_index1.*x1[:,1] + drop_index3.*x3[:,1] + (1 .- drop_index1 .-  drop_index3) .* a[:,1];
        a[:,2] = drop_index1.*x1[:,2] + drop_index2.*x2[:,2] + (1 .- drop_index1 .-  drop_index2) .* a[:,2];
        b[:,1] = drop_index4.*x4[:,1] + drop_index2.*x2[:,1] + (1 .- drop_index2 .-  drop_index4) .* b[:,1];
        b[:,2] = drop_index4.*x4[:,2] + drop_index3.*x3[:,2] + (1 .- drop_index3 .-  drop_index4) .* b[:,2];
        d  = b - a;
        x1 = a + alpha1 * d;
        x4 = a + alpha2 * d;
        x2 = copy(x1);
        x3 = copy(x4);
        x3[:,1] = x1[:,1];
        x2[:,1] = x4[:,1];
        f_mat  = zeros(n_tmp[1],4);
        f_mat[:,1] = f(x1,arg...);
        f_mat[:,2] = f(x2,arg...);
        f_mat[:,3] = f(x3,arg...);
        f_mat[:,4] = f(x4,arg...);
        conv = maximum(d);
        #println(conv,sum(drop_index))
    end
    return x1, f_mat[:,1]
end
function a_grid_fine_gen(agrid_fine_tmp,n_tmp,n_fine_tmp,agrid,a_min,a_max,curve_fine)
    size_agrid_fine = size(agrid_fine_tmp)[1]
    ii = 1;
    while n_fine_tmp> size_agrid_fine
        agrid_fine_tmp1 = range(a_min^curve_fine,a_max^curve_fine,length=
            (n_fine_tmp-n_tmp+ii)).^(1.0/curve_fine);
        agrid_fine_tmp = sort(union(agrid,agrid_fine_tmp1));
        size_agrid_fine = size(agrid_fine_tmp)[1];
        ii = ii+1;
    end
    return agrid_fine_tmp
end
function a_grid_fine_gen_midpoints(agrid,a_min,a_max,multiplier,n)
    agrid_tmp = zeros(multiplier*n[1]);
    agrid_tmp[1] =a_min;
    tmp_col_it_lo = 2;
    tmp_col_it_up = multiplier+1;
    for ii = 1:(n[1]-2)
        agrid_tmp[tmp_col_it_lo:tmp_col_it_up] =  range(agrid[ii],agrid[ii+1],
        length = tmp_col_it_up - tmp_col_it_lo+2)[2:end];
        tmp_col_it_lo = tmp_col_it_lo + multiplier;
        tmp_col_it_up = tmp_col_it_up + multiplier;
    end
    top_end_node_no = multiplier*n[1];
    agrid_tmp[tmp_col_it_lo:top_end_node_no]= (range(agrid[end-1],agrid[end],
    length = top_end_node_no-tmp_col_it_lo+2)[2:end]);
    return agrid_tmp
end

function setup_state_space(parameters_tmp::Parameter_type)

    # Labor productivity shock
    #parameters_tmp = copy(Baseline_parameter);
    #no_labor_shocks =2;
    #l_z = [1.0+parameters_tmp.l_z_low,0.0];
    #l_z = [1.0,parameters_tmp.l_z_low];
    #u_rate = 0.25
    #P_l = [0.8,0.6,0.2,0.4];
    #P_l = [parameters_tmp.ρ_W,parameters_tmp.ρ_W,1-parameters_tmp.ρ_W,1-parameters_tmp.ρ_W];
    #P_l = reshape(P_l,2,2)
    #P_l^100

    # #Productivity shock - rural
    # no_prod_shocks = convert(Int64,parameters_tmp.n[2]/parameters_tmp.no_labor_shocks);
    # (P_Rouw, z_Rouw) = setup_rouwen(parameters_tmp.ρ_S,0,
    #     parameters_tmp.σ_S /sqrt(1-parameters_tmp.ρ_S^2),no_prod_shocks );
    # exp_egrid = exp.(z_Rouw);

    # #Productivity shock - urban
    no_urban_prod_shocks = parameters_tmp.no_labor_shocks;
    # (P_urban_Rouw, z_urban_Rouw) = setup_rouwen( parameters_tmp.ρ_W,parameters_tmp.A_W,
    #     parameters_tmp.σ_W /sqrt(1-parameters_tmp.ρ_W^2),no_urban_prod_shocks);
    # exp_urban_egrid = exp.(z_urban_Rouw);
    # P = kron(P_urban_Rouw,P_Rouw)';

    A0x = [parameters_tmp.ρ_W parameters_tmp.ρ_SW; parameters_tmp.ρ_SW parameters_tmp.ρ_S];
    vex = [parameters_tmp.σ_W^2 0.0; 0.0 parameters_tmp.σ_S^2];

    P, Y = var_Markov_MM(A0x, vex, Int64(parameters_tmp.n[2] / parameters_tmp.no_labor_shocks), 1000)
    exp_egrid = exp.(Y[1:no_urban_prod_shocks,2])
    exp_urban_egrid = exp.(Y[1:no_urban_prod_shocks:end,1]);
    exp_urban_egrid = exp_urban_egrid .* exp(parameters_tmp.A_W);
    # Most rapid shock is urban prod
    # Combine the two shocks
#    P = zeros(parameters_tmp.n[2],parameters_tmp.n[2]);
#    for i= 1:parameters_tmp.no_labor_shocks
#        for j= 1:parameters_tmp.no_labor_shocks
#            tmp_ind_sta1 = (i-1)*no_prod_shocks + 1;
#            tmp_ind_sta2 = (j-1)*no_prod_shocks + 1;
#            tmp_ind_end1 = i*no_prod_shocks;
#            tmp_ind_end2 = j*no_prod_shocks;
#            P[tmp_ind_sta1:tmp_ind_end1,tmp_ind_sta2:tmp_ind_end2] = P_urban_Rouw[i,j]* P_Rouw'
#        end
#    end

    check_nonzero = P.>1e-10;
    P = P.*check_nonzero;
    P = P'
    # Set up the function spaces
    fspace = fundef((:spli, parameters_tmp.agrid, 0,parameters_tmp.spliorder),(:spli, range(1,parameters_tmp.n[2],step=1),0,1));# joint
    fspace_fine = fundef((:spli, parameters_tmp.agrid_fine, 0,1),(:spli, range(1,parameters_tmp.n[2],step=1),0,1));# joint distribution
    fspace_x = fundef((:spli, range(1,parameters_tmp.n[2],step=1),0,1));# only productivity
    # Set up the grids
    grid_fun = funnode(fspace);
    s    = grid_fun[1];
    ns = length(s[:,1]);
    grid_fine = funnode(fspace_fine);
    s_fine    = grid_fine[1];
    ns_fine   = length(s_fine[:,1]);
    # Set up the basis matrices
    Phi_z = funbase(fspace_x, s[:,2]);
    Phi_z_fine = funbase(fspace_x, s_fine[:,2]);
    Phi = funbase(fspace, s);
    check_nonzero = Phi.>1e-10;
    Phi = Phi.*check_nonzero;
    Phi_aug =  kron(Matrix(1.0I, 3, 3) ,Phi);
    P_kron = kron(P,Matrix(1.0I, parameters_tmp.n[1], parameters_tmp.n[1]));
    P_kron1 =  kron(Matrix(1.0I, 3, 3),P_kron);
    P_kron_fine = kron(ones(parameters_tmp.n_fine[1],1),P);
    #Eliminate small entries
    check_nonzero = P_kron.>1e-10;
    P_kron = P_kron.*check_nonzero;
    check_nonzero = P_kron1.>1e-10;
    P_kron1 = P_kron1.*check_nonzero;
    check_nonzero = P_kron_fine.>1e-10;
    P_kron_fine = P_kron_fine.*check_nonzero;
    #Force the rest to be in sparse format
    P_kron = SparseArrays.sparse(P_kron);
    P_kron1 = SparseArrays.sparse(P_kron1);
    P_kron_fine = SparseArrays.sparse(P_kron_fine);
    return (s,ns,s_fine,ns_fine,Phi_z,Phi_z_fine,Phi,Phi_aug,P_kron,P_kron1,P_kron_fine,exp_egrid,exp_urban_egrid)
end

function local_parameters(parameters_tmp::Parameter_type)
    δ = parameters_tmp.δ;
    ζ = parameters_tmp.ζ;
    ρ = parameters_tmp.ρ;
    α = parameters_tmp.α;
    σ = parameters_tmp.σ;
    β = parameters_tmp.β;
    ϵ = parameters_tmp.ϵ;
    ψ_S = parameters_tmp.ψ_S;
    ψ_B = parameters_tmp.ψ_B;
    ψ_M = parameters_tmp.ψ_M;
    ϕ_S = parameters_tmp.ϕ_S;
    ϕ_B = parameters_tmp.ϕ_B;
    c̄_S = parameters_tmp.c̄_S;
    F_W = parameters_tmp.F_W;
    F_S = parameters_tmp.F_S;
    F_B = parameters_tmp.F_B;
    FM_W = parameters_tmp.FM_W;
    FM_S = parameters_tmp.FM_S;
    FM_B = parameters_tmp.FM_B;
    Q_S = parameters_tmp.Q_S;
    p_x = parameters_tmp.p_x;
    τ_S = parameters_tmp.τ_S;
    τ_B = parameters_tmp.τ_B;
    a_D = parameters_tmp.a_D;
    b_D = parameters_tmp.b_D;
    K_a = parameters_tmp.K_a;
    K_b = parameters_tmp.K_b;
    γ = parameters_tmp.γ;
    A_W = parameters_tmp.A_W;
    ρ_S = parameters_tmp.ρ_S;
    ρ_SW= parameters_tmp.ρ_SW;
    σ_S = parameters_tmp.σ_S;
    ρ_W = parameters_tmp.ρ_W;
    σ_W = parameters_tmp.σ_W;
    n = parameters_tmp.n;
    n_fine = parameters_tmp.n_fine;
    C_grid_fine_no = parameters_tmp.C_grid_fine_no;
    agrid = parameters_tmp.agrid;
    agrid_fine = parameters_tmp.agrid_fine;
    C_grid_fine = parameters_tmp.C_grid_fine;
    a_min = parameters_tmp.a_min;
    a_max = parameters_tmp.a_max;
    spliorder = parameters_tmp.spliorder;
    fspace_a = parameters_tmp.fspace_a;
    fspace_a_fine = parameters_tmp.fspace_a_fine;
    fspace_C_fine = parameters_tmp.fspace_C_fine;
    s = parameters_tmp.s;
    ns = parameters_tmp.ns;
    s_fine = parameters_tmp.s_fine;
    ns_fine = parameters_tmp.ns_fine;
    z = parameters_tmp.z;
    z_W = parameters_tmp.z_W;
    #l_z = parameters_tmp.l_z;
    Phi_z = parameters_tmp.Phi_z;
    Phi_z_fine = parameters_tmp.Phi_z_fine;
    Phi = parameters_tmp.Phi;
    Phi_aug = parameters_tmp.Phi_aug;
    P_kron = parameters_tmp.P_kron;
    P_kron1 = parameters_tmp.P_kron1;
    P_kron_fine = parameters_tmp.P_kron_fine;
    κ = parameters_tmp.κ;
    return (δ,ζ,ρ,α,σ,β,ϵ,ψ_S,ψ_B,ψ_M,ϕ_S,ϕ_B,c̄_S,F_W,F_S,F_B,FM_W,FM_S,FM_B,Q_S,p_x,τ_S,τ_B,a_D,b_D,K_a,K_b,γ,A_W,
    ρ_S,ρ_SW,σ_S,ρ_W,σ_W,n,n_fine,agrid,agrid_fine,a_min,a_max,spliorder,fspace_a,fspace_a_fine,fspace_C_fine,C_grid_fine_no,C_grid_fine,s,ns,
    s_fine,ns_fine,z,z_W,Phi_z,Phi_z_fine,Phi,Phi_aug,P_kron,P_kron1,P_kron_fine,κ)
end

function local_parameters_ext(parameters_tmp::Parameter_type, epsilon_r::Float64, cttilde::Float64, ctilde::Float64)
    δ = parameters_tmp.δ
    ζ = parameters_tmp.ζ
    ρ = parameters_tmp.ρ
    α = parameters_tmp.α
    σ = parameters_tmp.σ
    β = parameters_tmp.β
    ϵ = parameters_tmp.ϵ
    ψ_S = parameters_tmp.ψ_S
    ψ_B = parameters_tmp.ψ_B
    ψ_M = parameters_tmp.ψ_M
    ϕ_S = parameters_tmp.ϕ_S * exp(-epsilon_r * (cttilde - 0.2))
    ϕ_B = parameters_tmp.ϕ_B * exp(-epsilon_r * (cttilde - 0.2))
    c̄_S = parameters_tmp.c̄_S
    F_W = parameters_tmp.F_W
    F_S = parameters_tmp.F_S
    F_B = parameters_tmp.F_B
    FM_W = parameters_tmp.FM_W
    FM_S = parameters_tmp.FM_S
    FM_B = parameters_tmp.FM_B
    Q_S = parameters_tmp.Q_S
    p_x = parameters_tmp.p_x
    τ_S = parameters_tmp.τ_S
    τ_B = parameters_tmp.τ_B
    a_D = parameters_tmp.a_D
    b_D = parameters_tmp.b_D
    K_a = parameters_tmp.K_a
    K_b = parameters_tmp.K_b
    γ = parameters_tmp.γ
    A_W = parameters_tmp.A_W
    ρ_S = parameters_tmp.ρ_S
    ρ_SW = parameters_tmp.ρ_SW
    σ_S = parameters_tmp.σ_S
    ρ_W = parameters_tmp.ρ_W
    σ_W = parameters_tmp.σ_W
    n = parameters_tmp.n
    n_fine = parameters_tmp.n_fine
    C_grid_fine_no = parameters_tmp.C_grid_fine_no
    agrid = parameters_tmp.agrid
    agrid_fine = parameters_tmp.agrid_fine
    C_grid_fine = parameters_tmp.C_grid_fine
    a_min = parameters_tmp.a_min
    a_max = parameters_tmp.a_max
    spliorder = parameters_tmp.spliorder
    fspace_a = parameters_tmp.fspace_a
    fspace_a_fine = parameters_tmp.fspace_a_fine
    fspace_C_fine = parameters_tmp.fspace_C_fine
    s = parameters_tmp.s
    ns = parameters_tmp.ns
    s_fine = parameters_tmp.s_fine
    ns_fine = parameters_tmp.ns_fine
    z = parameters_tmp.z
    z_W = parameters_tmp.z_W
    #l_z = parameters_tmp.l_z;
    Phi_z = parameters_tmp.Phi_z
    Phi_z_fine = parameters_tmp.Phi_z_fine
    Phi = parameters_tmp.Phi
    Phi_aug = parameters_tmp.Phi_aug
    P_kron = parameters_tmp.P_kron
    P_kron1 = parameters_tmp.P_kron1
    P_kron_fine = parameters_tmp.P_kron_fine
    κ = parameters_tmp.κ
    return (δ, ζ, ρ, α, σ, β, ϵ, ψ_S, ψ_B, ψ_M, ϕ_S, ϕ_B, c̄_S, F_W, F_S, F_B, FM_W, FM_S, FM_B, Q_S, p_x, τ_S, τ_B, a_D, b_D, K_a, K_b, γ, A_W,
        ρ_S, ρ_SW, σ_S, ρ_W, σ_W, n, n_fine, agrid, agrid_fine, a_min, a_max, spliorder, fspace_a, fspace_a_fine, fspace_C_fine, C_grid_fine_no, C_grid_fine, s, ns,
        s_fine, ns_fine, z, z_W, Phi_z, Phi_z_fine, Phi, Phi_aug, P_kron, P_kron1, P_kron_fine, κ)
end

function price_reshaper(prices::Array{Float64,1},δ::Float64,ϵ::Float64,
    ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,p_x::Float64,τ_S::Float64,
    τ_B::Float64,α::Float64,balanced_share::Float64)
    p_B= prices[1];
    p_M = prices[2];
    r = prices[3];

    if balanced_share>0.0
        τ_W = prices[4];
    else
        τ_W = 0.0
    end

    #w = prices[4];
    R = r + δ;
    w = p_M*(1-α)/(1+τ_W)*(R/α*1/p_M)^(α/(α-1));
    #residual = zeros(5);
    #P = (ψ_S^ϵ + p_B^(1 - ϵ) * ψ_B^ϵ + p_M^(1 - ϵ) * ψ_M^ϵ);
    return  p_B,p_M,R,r,w,τ_W#,residual
end
function price_reshaper_fixed_r(prices::Array{Float64,1},δ::Float64,ϵ::Float64,
    ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,p_x::Float64,τ_S::Float64,
    τ_B::Float64,α::Float64,balanced_share::Float64,r::Float64)
    p_B= prices[1];
    p_M = prices[2];
    if balanced_share>0.0
        τ_W = prices[3];
    else
        τ_W = 0.0
    end
    #p_x = p_x*p_M;
    #w = prices[4];
    R = r + δ;
    w = p_M*(1-α)/(1+τ_W)*(R/α/p_M)^(α/(α-1));
    #residual = zeros(5);
    #P = (ψ_S^ϵ + p_B^(1 - ϵ) * ψ_B^ϵ + p_M^(1 - ϵ) * ψ_M^ϵ);
    return  p_B,p_M,R,r,w,τ_W#,p_x#,residual
end
function price_reshaper_no_manu_price(prices::Array{Float64,1},δ::Float64,ϵ::Float64,
    ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,p_x::Float64,τ_S::Float64,
    τ_B::Float64)
    p_B= prices[1];
    r = prices[2];
    w = prices[3];
    residual = zeros(5);
    R = r + δ;
    p_M = R/ α * (w/(1 - α))^(1 - α)
    return  p_B,p_M,R,r,w,residual
end
function calibr_reshaper1(δ::Float64,ϵ::Float64,
    ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,p_x::Float64,τ_S::Float64,
    τ_B::Float64,L::Float64,r::Float64,K_L_ratio::Float64)
    capital_used = K_L_ratio * L
    τ_W = 0.0;
    p_M = (r + δ)/α * K_L_ratio^(1 - α);
    w = (1 - α)/(1 + τ_W)*p_M* K_L_ratio^α;
    R = r + δ;
    return  p_M,R,r,w,capital_used,τ_W
end
function calibr_reshaper2(δ::Float64,ϵ::Float64,
    ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,p_x::Float64,τ_S::Float64,
    τ_B::Float64,L::Float64,r::Float64)
    w = 1.0;
    R = r + δ;
    K_L_ratio = w/R  * α/(1 - α)
    capital_used = K_L_ratio * L
    τ_W = 0.0;
    p_M = R /α * K_L_ratio^(1 - α);
    return  p_M,R,w,capital_used,τ_W
end
function calibr_reshaper3(prices::Array{Float64,1}, δ::Float64,ϵ::Float64,
    ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,p_x::Float64,τ_S::Float64,
    τ_B::Float64,L::Float64,r::Float64,α::Float64)
    p_B= prices[1];
    p_M= prices[2];
    R = r + δ;
    K_L_ratio = (p_M/R*α)^(1/(1 - α))
    capital_used = K_L_ratio * L
    τ_W = 0.0;
    w = (1 - α)/(1 + τ_W)*p_M* K_L_ratio^α;
    return  p_B,p_M,R,w,capital_used,τ_W
end
# New staple producer problem:
function λ_2_staple_residual(λ_2::Array{Float64,1},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,
    p_B::Float64,p_M::Float64,ϕ_B::Float64,c̄_S::Float64,Q_S::Float64,z::Array{Float64,1},
    ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,C::Array{Float64,1},z_index::Int64)
    x_S_tmp = ((1 + τ_S) * p_x./(λ_2 * ζ * ϕ_S *z[z_index] )).^(1 / (ζ - 1));
    q_s_tmp = ϕ_S * x_S_tmp.^ζ *z[z_index];
    P_S_tmp = (λ_2.^(1 - ϵ)*ψ_S^ϵ .+ p_B^(1 -ϵ)*ψ_B^ϵ .+ p_M^(1 -ϵ)*ψ_M^ϵ).^(1 / (1 - ϵ));
    c_s_tmp = c̄_S .+ λ_2.^( -ϵ).*ψ_S.^ϵ.*C.*P_S_tmp.^ϵ
    return -(q_s_tmp - c_s_tmp).^2
end
function coeff_λ_2_s_approx_staples(coeff_λ_2_s::Array{Float64,2},C::Array{Float64,1},fspace_C_fine::Dict{Symbol,Any})
    Phi_prime_a_fine = funbase(fspace_C_fine, C);
    return diag(Phi_prime_a_fine * coeff_λ_2_s);
end
function λ_2_staple_constrained(C::Array{Float64,1},θ::Array{Float64,1},
    fspace_C_fine::Dict{Symbol,Any},s::Array{Float64,2},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,
        p_B::Float64,p_M::Float64,ϕ_B::Float64,c̄_S::Float64,Q_S::Float64,
        ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,κ::Float64,ns::Int64,TC_S_c3_constrained::Array{Float64,1},
        c_S_c3_constrained::Array{Float64,1},a_min::Float64)
    condition = max.((((c_S_c3_constrained .- c̄_S)./ψ_S.^ϵ./C).^((1 - ϵ)/ϵ) .- ψ_S.^ϵ),a_min);
    λ_2_constrained_tmp = ((p_B^(1 - ϵ) *ψ_B^ϵ + p_M^(1 - ϵ) *ψ_M^ϵ ) ./condition).^(1 /(1 - ϵ));
    λ_S_constrained_tmp = (λ_2_constrained_tmp .* ζ .* ϕ_S .* θ)./ ((TC_S_c3_constrained ).^(1 - ζ) .* ((1 + τ_S) * p_x).^ζ) .- 1;
    #λ_S_constrained_tmp1 = (((c_S_c3_constrained .- c̄_S).^((1 - ϵ)/ϵ)*ψ_S.^(ϵ-1).*C.^(( ϵ - 1)/ϵ) .-  ψ_S.^ϵ).*c_S_c3_constrained.^((1- ζ)*(1 -ϵ )/ζ
    #)./((p_B^(1 - ϵ) *ψ_B^ϵ + p_M^(1 - ϵ) *ψ_M^ϵ ) .* ((1 + τ_S) * p_x).^( ϵ - 1))).^(1/( ϵ - 1)) .* ζ .* (ϕ_S .* θ).^(1/ζ) .- 1
    return λ_2_constrained_tmp,λ_S_constrained_tmp
end
function staples_c3_objects(C::Array{Float64,1},ϕ_S::Float64,τ_S::Float64,p_x::Float64,
    p_B::Float64,p_M::Float64,ϕ_B::Float64,τ_B::Float64,Q_S::Float64,ϵ::Float64,ψ_S::Float64,
    ψ_B::Float64,ψ_M::Float64,ζ::Float64,θ::Array{Float64,1},coeff_λ_2_s::Array{Float64,2},fspace_C_fine::Dict{Symbol,Any},
    s::Array{Float64,2},κ::Float64,c̄_S::Float64,C_max_staple::Array{Float64,1},C_min_staple::Array{Float64,1},ns::Int64,C_max_staple_constrained::Array{Float64,1},
    C_min_staple_constrained::Array{Float64,1},TC_S_c3_constrained::Array{Float64,1},
    x_S_c3_constrained::Array{Float64,1},q_S_c3_constrained::Array{Float64,1},c_S_c3_constrained::Array{Float64,1},a_min::Float64)
    TC_mat_S_c3 = zeros(ns,2);
    λ_2_mat_S_c3 = copy(TC_mat_S_c3);
    λ_S_mat_S_c3 = copy(TC_mat_S_c3);
    x_S_mat_S_c3  = copy(TC_mat_S_c3);
    q_S_mat_S_c3  = copy(TC_mat_S_c3);
    c_S_mat_S_c3  = copy(TC_mat_S_c3);
    feasibility_mat_S_c3 = copy(TC_mat_S_c3);
    λ_2_tmp = coeff_λ_2_s_approx_staples(coeff_λ_2_s,C,fspace_C_fine);
    λ_2_mat_S_c3[:,1] = λ_2_tmp;
    ( λ_2_mat_S_c3[:,2],λ_S_mat_S_c3[:,2]) = λ_2_staple_constrained(C,θ,fspace_C_fine,s,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,κ,ns,TC_S_c3_constrained,c_S_c3_constrained,a_min);
    TC_mat_S_c3[:,2] = TC_S_c3_constrained;
    x_S_mat_S_c3[:,2] = TC_S_c3_constrained./ ((1 + τ_S) * p_x)
    x_S_tmp = ((1 + τ_S) * p_x./(λ_2_tmp .* ζ .* ϕ_S .*θ )).^(1 / (ζ - 1));
    x_S_mat_S_c3[:,1] =x_S_tmp;
    TC_mat_S_c3[:,1]  = (1 + τ_S) * p_x * x_S_tmp;
    TC_true, constrained_staple_index = findmin(TC_mat_S_c3; dims=2);
    x_S_tmp = x_S_mat_S_c3[constrained_staple_index];
    q_S_tmp = ϕ_S * x_S_tmp.^ζ .*θ;
    c_s_tmp = q_S_tmp;
    λ_2_tmp = λ_2_mat_S_c3[constrained_staple_index];
    λ_S_tmp = λ_S_mat_S_c3[constrained_staple_index];
    shadow_price = min.(max.(1,λ_2_tmp),(1 + Q_S));
    P_S_c3_tmp = (shadow_price.^(1 - ϵ)*ψ_S^ϵ .+ p_B^(1 -ϵ)*ψ_B^ϵ .+ p_M^(1 -ϵ)*ψ_M^ϵ).^(1 / (1 - ϵ));
    c_B_tmp = p_B^( -ϵ)*ψ_B^ϵ * C .* P_S_c3_tmp.^ϵ;
    c_M_tmp = p_M^( -ϵ)*ψ_M^ϵ * C .* P_S_c3_tmp.^ϵ;
    Y_S_tmp = p_B * c_B_tmp + c_M_tmp * p_M + (1 + τ_S) * p_x * x_S_tmp;
    feasibility_mat_S_c3[:,1] =  (C_max_staple.<C) +  ((C_min_staple).>C);
    feasibility_mat_S_c3[:,2] = (C_max_staple_constrained.<C) +  ((C_min_staple_constrained).>C);
    feasibility = feasibility_mat_S_c3[constrained_staple_index] + (λ_2_tmp.>(1 + Q_S)) + (λ_2_tmp.<1) + (λ_S_tmp.<0.0 );
    return Y_S_tmp,P_S_c3_tmp,c_s_tmp,c_B_tmp,c_M_tmp,x_S_tmp,shadow_price,constrained_staple_index,feasibility
end
function staples_objects(C::Array{Float64,1},ϕ_S::Float64,τ_S::Float64,p_x::Float64,
    p_B::Float64,p_M::Float64,ϕ_B::Float64,τ_B::Float64,Q_S::Float64,
    ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,ζ::Float64,θ::Array{Float64,1},coeff_λ_2_s::Array{Float64,2},fspace_C_fine::Dict{Symbol,Any},
    P_S_c1::Float64,P_S_c2::Float64,Y_S_c1::Array{Float64,1},Y_S_c2::Array{Float64,1},
    x_S_c1::Array{Float64,1}, x_S_c2::Array{Float64,1},s::Array{Float64,2},q_S_c1::Array{Float64,1},q_S_c2::Array{Float64,1},q_S_staples::Array{Float64,2},
    c_S_staples::Array{Float64,2},c_B_staples::Array{Float64,2},c_M_staples::Array{Float64,2},P_S_staples::Array{Float64,2},x_S_staples::Array{Float64,2},
    λ_2_S_staples::Array{Float64,2},unfeasible_mat::Array{Float64,2},Y_S_potential::Array{Float64,2},κ::Float64,c̄_S::Float64,C_max_staple::Array{Float64,1},
    C_min_staple::Array{Float64,1},ns::Int64,C_max_staple_constrained::Array{Float64,1},
    C_min_staple_constrained::Array{Float64,1},TC_S_c3_constrained::Array{Float64,1},
    x_S_c3_constrained::Array{Float64,1},q_S_c3_constrained::Array{Float64,1},c_S_c3_constrained::Array{Float64,1},a_min::Float64)
    Y_S_c3_tmp,P_S_c3_tmp,c_s_c3_tmp,c_B_C3_tmp,c_M_c3_tmp,x_S_c3_tmp,shadow_price_c3_tmp,constrained_staple_index_c3_tmp,feasibility_c3_tmp = staples_c3_objects(C,
    ϕ_S,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ζ,θ,coeff_λ_2_s,fspace_C_fine,s,κ,c̄_S,C_max_staple,C_min_staple,ns,C_max_staple_constrained,
    C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained,a_min);
    # First check production of staples - some paths are not feasible:
    q_S_staples[:,1] = q_S_c1;
    q_S_staples[:,2] = q_S_c2;
    q_S_staples[:,3] = ϕ_S *θ .* x_S_c3_tmp.^ζ;
    c_S_staples[:,1] = c̄_S .+ ψ_S.^ϵ.*C.*P_S_c1.^ϵ;
    c_S_staples[:,2] = c̄_S .+ (1 + Q_S)^( -ϵ) * ψ_S.^ϵ.*C.*P_S_c2.^ϵ;
    c_S_staples[:,3] = c_s_c3_tmp;
    c_B_staples[:,1] = p_B^( -ϵ)*ψ_B^ϵ * C .* P_S_c1.^ϵ;
    c_B_staples[:,2] = p_B^( -ϵ)*ψ_B^ϵ * C .* P_S_c2.^ϵ;
    c_B_staples[:,3] = c_B_C3_tmp;
    c_M_staples[:,1] = p_M^( -ϵ)*ψ_M^ϵ * C .* P_S_c1.^ϵ;
    c_M_staples[:,2] = p_M^( -ϵ)*ψ_M^ϵ * C .* P_S_c2.^ϵ;
    c_M_staples[:,3] = c_M_c3_tmp;
    P_S_staples[:,1] .= P_S_c1;
    P_S_staples[:,2] .= P_S_c2;
    P_S_staples[:,3] = P_S_c3_tmp;
    x_S_staples[:,1] = x_S_c1;
    x_S_staples[:,2] = x_S_c2;
    x_S_staples[:,3] = x_S_c3_tmp;
    λ_2_S_staples[:,1] .= 1;
    λ_2_S_staples[:,2] .= 1 + Q_S;
    λ_2_S_staples[:,3] = shadow_price_c3_tmp;
    unfeasible_mat[:,1] =  q_S_staples[:,1].< c_S_staples[:,1];
    Y_S_potential[:,1] = Y_S_c1 +  c_S_staples[:,1]  +  p_B * c_B_staples[:,1] + p_M * c_M_staples[:,1] + unfeasible_mat[:,1]*1000; # was : Y_S_potential[:,1] = Y_S_c1 +  P_S_staples[:,1].*C + unfeasible_mat[:,1]*1000; #Penalty
    #Y_S_potential[:,2] = Y_S_c2 +  (1 + Q_S) .*(c_S_staples[:,2] - q_S_c2) + unfeasible_mat_cbar[:,2].*1000 + unfeasible_mat_cs_qs_case2.*1000 +  p_B * c_B_staples[:,2] + p_M * c_M_staples[:,2]; #always feasible
    Y_S_potential[:,2] = Y_S_c2 +  Q_S .*(c_S_staples[:,2] - q_S_c2).*(c_S_staples[:,2].>q_S_c2) + (c_S_staples[:,2] - q_S_c2) +  p_B * c_B_staples[:,2] + p_M * c_M_staples[:,2]; #always feasible
    # The price index is still valid, but this is easier.
    Y_S_potential[:,3] = Y_S_c3_tmp + feasibility_c3_tmp*1000; # Y_S_tmp already contains the expenses, see staples_c3_objects
    Y_S, solve_staple_index = findmin(Y_S_potential; dims=2);
    c_S = c_S_staples[solve_staple_index];
    c_B = c_B_staples[solve_staple_index];
    c_M = c_M_staples[solve_staple_index];
    q_S = q_S_staples[solve_staple_index];
    P_S = P_S_staples[solve_staple_index];
    x_S = x_S_staples[solve_staple_index];
    λ_2_S = λ_2_S_staples[solve_staple_index];
    #sum(κ.*s[:,1]./p_x-x_S.>=0) check if any of the fertilizer use violates the conditions
    #sum(getindex.(solve_staple_index,2).==3)
    return Y_S,c_S,c_B,c_M,q_S,P_S,x_S,solve_staple_index,λ_2_S
end

#New cashcrop objects
function λ_2_residual_unconstrained(λ_2::Array{Float64,1},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,
    p_B::Float64,p_M::Float64,ϕ_B::Float64,τ_B::Float64,c̄_S::Float64,Q_S::Float64,z::Array{Float64,1},
    ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,C::Array{Float64,1},z_index::Int64,ρ::Float64)
    labor_allocated_interior_bar_c3c_tmp = (p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ))./((p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ)) .+ ((ϕ_S * λ_2) * (1+τ_B)^ζ).^(1/(1-ρ-ζ)))
    x_SC_interior_c3c_tmp = (((ϕ_S .* λ_2) .* ζ .* (1 .- labor_allocated_interior_bar_c3c_tmp).^ρ)./((1 + τ_S) * p_x)).^(1/(1 - ζ)).*z[z_index].^(1/(1 - ζ));
    q_S_c3c_tmp = ϕ_S *z[z_index] .* x_SC_interior_c3c_tmp.^ζ .* (1.0 .- labor_allocated_interior_bar_c3c_tmp).^ρ;
    P_B_c3c_tmp =  (λ_2.^(1 -ϵ) *ψ_S^ϵ .+ p_B^(1 -ϵ)*ψ_B^ϵ .+ p_M^(1 -ϵ)*ψ_M^ϵ).^(1 / (1 - ϵ));
    c_S_c3c_tmp = c̄_S .+ λ_2.^(-ϵ) .* P_B_c3c_tmp.^ϵ.*ψ_S.^ϵ.*C;
    #LASZLO Note:  goldenx should work, because c_S_c3c_tmp is decreasing in λ_2, while q_S_c3c_tmp is now increasing in λ_2.
    #Therefore there can only be either one solution, or zero
    return -(c_S_c3c_tmp - q_S_c3c_tmp).^2
end
function λ_2_residual_constrained(λ_2::Array{Float64,1},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,
    p_B::Float64,p_M::Float64,ϕ_B::Float64,τ_B::Float64,c̄_S::Float64,Q_S::Float64,z::Array{Float64,1},
    ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,C::Array{Float64,1},z_index::Int64,a::Array{Float64,1},a_index::Int64,κ::Float64,ρ::Float64)
    λ_B_tmp = ζ * z[z_index] ./ (κ*a[a_index]).^(1-ζ) .* ( (p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ)) .+ ((ϕ_S * λ_2) * (1+τ_B)^ζ).^(1/(1-ρ-ζ)) ).^(1-ρ-ζ) ./ (p_x*(1+τ_S)*(1+τ_B)).^ζ .- 1;
    λ_B_tmp = max.(λ_B_tmp,0.0);
    Lagrange_factor = (1 .+ λ_B_tmp).^(1/(ζ - 1));
    labor_allocated_interior_bar_c3c_tmp = (p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ))./((p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ)) .+ ((ϕ_S * λ_2) * (1+τ_B)^ζ).^(1/(1-ρ-ζ)))
    x_SC_interior_c3c_tmp = Lagrange_factor .* (((ϕ_S .* λ_2) .* ζ .* (1 .- labor_allocated_interior_bar_c3c_tmp).^ρ)./((1 + τ_S) * p_x)).^(1/(1 - ζ)).*z[z_index].^(1/(1 - ζ));
    q_S_c3c_tmp = ϕ_S *z[z_index] .* x_SC_interior_c3c_tmp.^ζ .* (1.0 .- labor_allocated_interior_bar_c3c_tmp).^ρ;
    P_B_c3c_tmp =  (λ_2.^(1 -ϵ) *ψ_S^ϵ .+ p_B^(1 -ϵ)*ψ_B^ϵ .+ p_M^(1 -ϵ)*ψ_M^ϵ).^(1 / (1 - ϵ));
    c_S_c3c_tmp = c̄_S .+ λ_2.^(-ϵ) .* P_B_c3c_tmp.^ϵ.*ψ_S.^ϵ.*C;
    return -(c_S_c3c_tmp - q_S_c3c_tmp).^2
end
function λ_2_approx_cashcrop(coeff_λ_2_cashcrop::Array{Float64,2},C::Array{Float64,1},fspace_C_fine::Dict{Symbol,Any})
    Phi_prime_a_fine = funbase(fspace_C_fine, C);
    return diag(Phi_prime_a_fine * coeff_λ_2_cashcrop);
end

function cashcrop_constrained_c3c(coeff_λ_2_cashcrop_residual_unconstrained::Array{Float64,2},coeff_λ_2_cashcrop_residual_constrained::Array{Float64,2},
    C::Array{Float64,1},θ::Array{Float64,1},fspace_C_fine::Dict{Symbol,Any},s::Array{Float64,2},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,p_B::Float64,
        p_M::Float64,ϕ_B::Float64,τ_B::Float64,c̄_S::Float64,Q_S::Float64,ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,ns::Int64,κ::Float64,tol::Float64,
        x_S_mat_tmp::Array{Float64,2},x_B_mat_tmp::Array{Float64,2},land_B_mat_tmp::Array{Float64,2},λ_2_mat_tmp::Array{Float64,2},TC_mat_tmp::Array{Float64,2},C_max_mat::Array{Float64,2},
        C_min_mat::Array{Float64,2},ρ::Float64,a_min::Float64)

    λ_2_mat_tmp[:,1] = λ_2_approx_cashcrop(coeff_λ_2_cashcrop_residual_unconstrained,C,fspace_C_fine);
    land_B_mat_tmp[:,1] = (p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ))./((p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ)) .+ ((ϕ_S * λ_2_mat_tmp[:,1]) * (1+τ_B)^ζ).^(1/(1-ρ-ζ)));
    x_S_mat_tmp[:,1] = (((ϕ_S .* λ_2_mat_tmp[:,1]) .* ζ .* (1 .- land_B_mat_tmp[:,1]).^ρ)./((1 + τ_S) * p_x)).^(1/(1 - ζ)).*θ.^(1/(1 - ζ));
    x_B_mat_tmp[:,1] = ((p_B * ϕ_B * ζ * land_B_mat_tmp[:,1].^ρ)./((1 + τ_B) * p_x)).^(1/(1 - ζ)).*θ.^(1/(1 - ζ));
    TC_mat_tmp[:,1] =p_x *  ((1 + τ_S) * x_S_mat_tmp[:,1] +  (1 + τ_B) * x_B_mat_tmp[:,1]);

    λ_2_constrained_tmp = λ_2_approx_cashcrop(coeff_λ_2_cashcrop_residual_constrained,C,fspace_C_fine);
    λ_2_mat_tmp[:,2] = λ_2_constrained_tmp;
    λ_B_tmp = ζ * θ ./ (κ*s[:,1]).^(1-ζ) .* ( (p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ)) .+ ((ϕ_S * λ_2_constrained_tmp) * (1+τ_B)^ζ).^(1/(1-ρ-ζ)) ).^(1-ρ-ζ) ./ (p_x*(1+τ_S)*(1+τ_B)).^ζ .- 1;
    λ_B_tmp = max.(λ_B_tmp,0.0);
    Lagrange_factor = (1 .+ λ_B_tmp).^(1/(ζ - 1));
    land_B_mat_tmp[:,2] = (p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ))./((p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ)) .+ ((ϕ_S * λ_2_mat_tmp[:,2]) * (1+τ_B)^ζ).^(1/(1-ρ-ζ)));
    x_S_mat_tmp[:,2] = (((ϕ_S .* λ_2_mat_tmp[:,2]) .* ζ .* (1 .- land_B_mat_tmp[:,2]).^ρ)./((1 + τ_S) * p_x)).^(1/(1 - ζ)).*θ.^(1/(1 - ζ)) .* Lagrange_factor;
    x_B_mat_tmp[:,2] = ((p_B * ϕ_B * ζ * land_B_mat_tmp[:,2].^ρ)./((1 + τ_B) * p_x)).^(1/(1 - ζ)).*θ.^(1/(1 - ζ)) .* Lagrange_factor;
    TC_mat_tmp[:,2] = p_x *  ((1 + τ_S) * x_S_mat_tmp[:,2] +  (1 + τ_B) * x_B_mat_tmp[:,2]);

    TC_mat_tmp_true, constrained_status = findmin(TC_mat_tmp, dims=2);
    labor_allocated_interior_true = land_B_mat_tmp[constrained_status];
    λ_2_true = λ_2_mat_tmp[constrained_status];
    x_S_true = x_S_mat_tmp[constrained_status];
    x_B_true = x_B_mat_tmp[constrained_status];
    q_S_true = ϕ_S *θ .* x_S_true.^ζ .* (1.0 .- labor_allocated_interior_true).^ρ;
    q_B_true = ϕ_B *θ .* x_B_true.^ζ .* labor_allocated_interior_true.^ρ;
    c_S_true = copy(q_S_true);
    P_B_true = (λ_2_true.^(1 -ϵ)*ψ_S^ϵ .+ p_B^(1 -ϵ)*ψ_B^ϵ .+ p_M^(1 -ϵ)*ψ_M^ϵ).^(1 / (1 - ϵ));
    c_B_true = p_B^( -ϵ)*ψ_B^ϵ * C .* P_B_true.^ϵ;
    c_M_true = p_M^( -ϵ)*ψ_M^ϵ * C .* P_B_true.^ϵ;
    C_max_true = C_max_mat[constrained_status];
    C_min_true = C_min_mat[constrained_status];
    feasibility =  (C_max_true.<C) +  ((C_min_true .- a_min).>C) + (λ_2_true.>(1 + Q_S)) + (λ_2_true.<1);

    Y_B_true =(- p_B * (q_B_true - c_B_true) + c_M_true * p_M + ((1 + τ_S) * p_x) * x_S_true + ((1 + τ_B) * p_x) * x_B_true);

    return c_S_true,c_B_true,c_M_true,x_S_true,labor_allocated_interior_true,λ_2_true,x_B_true,P_B_true,Y_B_true,q_S_true,q_B_true,feasibility,TC_mat_tmp_true
end

function cashcrop_objects(C::Array{Float64,1},coeff_λ_2_cashcrop_residual_unconstrained::Array{Float64,2},coeff_λ_2_cashcrop_residual_constrained::Array{Float64,2},
    θ::Array{Float64,1},fspace_C_fine::Dict{Symbol,Any},s::Array{Float64,2},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,p_B::Float64,
        p_M::Float64,ϕ_B::Float64,τ_B::Float64,c̄_S::Float64,Q_S::Float64,ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,ns::Int64,κ::Float64,tol::Float64,
        a_min::Float64,x_B_c1::Array{Float64,1},π_B_only_B_c1::Array{Float64,1},λ_B_only_B_c1::Array{Float64,1},P_B_c1::Float64,Y_B_c1::Array{Float64,1},
        coeff_λ_2_s::Array{Float64,2},P_S_c1::Float64,P_S_c2::Float64,Y_S_c1::Array{Float64,1},Y_S_c2::Array{Float64,1},x_S_c1::Array{Float64,1}, x_S_c2::Array{Float64,1},
        labor_allocated_interior_c3a::Float64,λ_B_interior_c3a::Array{Float64,1},x_SC_interior_c3a::Array{Float64,1},
        x_BC_interior_c3a::Array{Float64,1},Y_B_c3a::Array{Float64,1},P_B_c3a::Float64,P_B_c3b::Float64,q_S_c1::Array{Float64,1},q_S_c2::Array{Float64,1},
        q_B_c1::Array{Float64,1},q_S_c3a::Array{Float64,1},q_B_c3a::Array{Float64,1},q_S_c3b::Array{Float64,1},q_B_c3b::Array{Float64,1},
        x_SC_interior_c3b::Array{Float64,1},x_BC_interior_c3b::Array{Float64,1},labor_allocated_interior_c3b::Float64,Y_B_c3b::Array{Float64,1},
        c_S_mat::Array{Float64,2},c_B_mat::Array{Float64,2},
        c_M_mat::Array{Float64,2},x_S_mat::Array{Float64,2},x_B_mat::Array{Float64,2},q_S_mat::Array{Float64,2},q_B_mat::Array{Float64,2},land_B_mat::Array{Float64,2},
        λ_2_mat::Array{Float64,2},P_B_mat::Array{Float64,2},Y_B_mat::Array{Float64,2},feasibility_mat::Array{Float64,2},
        C_max_mat::Array{Float64,2},C_min_mat::Array{Float64,2},TC_mat::Array{Float64,2},
        q_S_staples::Array{Float64,2},c_S_staples::Array{Float64,2},c_B_staples::Array{Float64,2},c_M_staples::Array{Float64,2},P_S_staples::Array{Float64,2},
        x_S_staples::Array{Float64,2},λ_2_S_staples::Array{Float64,2},unfeasible_mat::Array{Float64,2},Y_S_potential::Array{Float64,2},ρ::Float64,C_max_staple::Array{Float64,1},
        C_min_staple::Array{Float64,1},C_max_staple_constrained::Array{Float64,1},
        C_min_staple_constrained::Array{Float64,1},TC_S_c3_constrained::Array{Float64,1},
        x_S_c3_constrained::Array{Float64,1},q_S_c3_constrained::Array{Float64,1},c_S_c3_constrained::Array{Float64,1},
        x_S_mat_3c::Array{Float64,2},x_B_mat_3c::Array{Float64,2},land_B_mat_3c::Array{Float64,2},λ_2_mat_3c::Array{Float64,2},TC_mat_3c::Array{Float64,2})




    # First check production of staples - some paths are not feasible:
    #1) Specialize in cash crop:
    c_S_mat[:,1] = c̄_S .+ (1 + Q_S)^( -ϵ) * ψ_S.^ϵ.*C.*P_B_c1.^ϵ;
    c_B_mat[:,1] = p_B^( -ϵ)*ψ_B^ϵ * C .* P_B_c1.^ϵ;
    c_M_mat[:,1] = p_M^( -ϵ)*ψ_M^ϵ * C .* P_B_c1.^ϵ; #  was before: c_M_mat[:,1] = p_M^( -ϵ)*ψ_B^ϵ * C .* P_B_c1.^ϵ;
    x_S_mat[:,1] .= 0;
    x_B_mat[:,1] = x_B_c1;
    land_B_mat[:,1] .= 1;
    λ_2_mat[:,1] .= (1 + Q_S);
    P_B_mat[:,1] .= P_B_c1;
    feasibility_mat[:,1] .= 0; # This case is always feasible
    Y_B_mat[:,1] = Y_B_c1 +  (1 + Q_S) * c_S_mat[:,1]  +  p_B * c_B_mat[:,1] + p_M * c_M_mat[:,1];
    q_S_mat[:,1] .= 0;
    q_B_mat[:,1] = q_B_c1;
    TC_mat[:,1] = p_x *   (1 + τ_B) * x_B_mat[:,1];
    #2) Specialize in staple crops:
    Y_B_mat[:,2],c_S_mat[:,2],c_B_mat[:,2],c_M_mat[:,2],q_S_mat[:,2],P_B_mat[:,2],x_S_mat[:,2],solve_staple_index,λ_2_mat[:,2] = staples_objects(C,ϕ_S,τ_S,p_x,
        p_B,p_M,ϕ_B,τ_B,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ζ,θ,coeff_λ_2_s,fspace_C_fine,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,
        x_S_c1, x_S_c2,s,q_S_c1,q_S_c2,q_S_staples,c_S_staples,c_B_staples,c_M_staples,P_S_staples,x_S_staples,
        λ_2_S_staples,unfeasible_mat,Y_S_potential,κ,c̄_S,C_max_staple,
        C_min_staple,ns,C_max_staple_constrained,C_min_staple_constrained,TC_S_c3_constrained,
        x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained,a_min);
    x_B_mat[:,2] .= 0;
    land_B_mat[:,2] .= 0;
    feasibility_mat[:,2] .= 0; # This case is always feasible
    q_B_mat[:,2] .= 0;
    TC_mat[:,2] = p_x *  (1 + τ_S) * x_S_mat[:,2];
    #3) Internal solution - produces both staples and cash crop
    #3a) q_S>c_S, hence no transaction cost paid


    λ_2_mat[:,3] .= 1;
    c_S_mat[:,3] = c̄_S .+ ψ_S.^ϵ.*C.*P_B_c3a.^ϵ;
    c_B_mat[:,3] = p_B^( -ϵ)*ψ_B^ϵ * C .* P_B_c3a.^ϵ;
    c_M_mat[:,3] = p_M^( -ϵ)*ψ_M^ϵ * C .* P_B_c3a.^ϵ; #similar as above
    x_S_mat[:,3] = x_SC_interior_c3a;
    x_B_mat[:,3] = x_BC_interior_c3a;
    land_B_mat[:,3] .=labor_allocated_interior_c3a;
    P_B_mat[:,3] .=P_B_c3a;
    q_S_mat[:,3] = q_S_c3a;
    q_B_mat[:,3] = q_B_c3a;
    feasibility_mat[:,3] .= q_S_mat[:,3].< c_S_mat[:,3];
    TC_mat[:,3] = p_x *  (1 + τ_S) * x_S_mat[:,3] + p_x *   (1 + τ_B) * x_B_mat[:,3];
    Y_B_mat[:,3] = Y_B_c3a +  c_S_mat[:,3]  +  p_B * c_B_mat[:,3] + p_M * c_M_mat[:,3];
    #3b) q_S<c_S, hence transaction cost is always paid marginally

    λ_2_mat[:,4] .= 1 + Q_S;
    c_S_mat[:,4] = c̄_S .+ (1 + Q_S)^( -ϵ) * ψ_S.^ϵ.*C.*P_B_c3b.^ϵ;
    c_B_mat[:,4] = p_B^( -ϵ)*ψ_B^ϵ * C .* P_B_c3b.^ϵ;
    c_M_mat[:,4] = p_M^( -ϵ)*ψ_M^ϵ * C .* P_B_c3b.^ϵ; # similar as above
    x_S_mat[:,4] = x_SC_interior_c3b;
    x_B_mat[:,4] = x_BC_interior_c3b;
    land_B_mat[:,4] .=labor_allocated_interior_c3b;
    P_B_mat[:,4] .=P_B_c3b;
    q_S_mat[:,4] = q_S_c3b;
    q_B_mat[:,4] = q_B_c3b;
    Y_B_mat[:,4] = Y_B_c3b +  (1 + Q_S) .*(c_S_mat[:,4] - q_S_mat[:,4])  +  p_B * c_B_mat[:,4] + p_M * c_M_mat[:,4] ;
    feasibility_mat[:,4] .= q_S_mat[:,4].> c_S_mat[:,4];
    TC_mat[:,4] = p_x *  (1 + τ_S) * x_S_mat[:,4] + p_x *   (1 + τ_B) * x_B_mat[:,4];
    #3c) q_S=c_S, hence transaction cost only distorts production and is not paid
    (c_S_mat[:,5],c_B_mat[:,5],c_M_mat[:,5],x_S_mat[:,5],land_B_mat[:,5],λ_2_mat[:,5],x_B_mat[:,5],P_B_mat[:,5],Y_B_mat[:,5],q_S_mat[:,5],q_B_mat[:,5],
    feasibility_mat[:,5],TC_mat[:,5])  = cashcrop_constrained_c3c(coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,
        C,θ,fspace_C_fine,s,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns,κ,tol,
        x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c,C_max_mat,C_min_mat,ρ,a_min);
    
    # For this case, we have additional conditions to check:
    # consumption is not equal to quantity produced
    #feasibility_mat[:,5] = feasibility_mat[:,5] + ((q_S_mat[:,5]- c_S_mat[:,5]).^2 .>tol);
    # consumption is not equal to quantity produced
    #    feasibility_mat[:,5] = feasibility_mat[:,5]; #+ (λ_2_mat[:,5] .> (1+ Q_S)) + (λ_2_mat[:,5] .< 1.0);
    Y_B_mat = Y_B_mat + 1000*feasibility_mat;
    Y_B, solve_cash_crop_index = findmin(Y_B_mat; dims=2);
    #solve_cash_crop_index_v = getindex.(solve_cash_crop_index,2);
    c_S = c_S_mat[solve_cash_crop_index];
    c_B = c_B_mat[solve_cash_crop_index];
    c_M = c_M_mat[solve_cash_crop_index];
    q_S = q_S_mat[solve_cash_crop_index];
    q_B = q_B_mat[solve_cash_crop_index];
    P_B = P_B_mat[solve_cash_crop_index];
    x_S = x_S_mat[solve_cash_crop_index];
    x_B = x_B_mat[solve_cash_crop_index];
    land_C= land_B_mat[solve_cash_crop_index];
    λ_2= λ_2_mat[solve_cash_crop_index];
    TC = TC_mat[solve_cash_crop_index];
    return c_S,c_B,c_M,x_S,x_B,land_C,λ_2,P_B,Y_B,q_S,q_B,solve_cash_crop_index,solve_staple_index,TC
end
function myCondition(y::Float64)
    return 0 .< y
end
function income_creator(s::Array{Float64,2},ns::Int64,
    z::Array{Float64,1},z_W::Array{Float64,1},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,
    p_B::Float64,p_M::Float64,ϕ_B::Float64,τ_B::Float64,ρ::Float64,w::Float64,r::Float64,
    c̄_S::Float64,a_min::Float64,a_max::Float64,γ::Float64,n::Array{Int64,1},κ::Float64,Q_S::Float64,
    ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,agrid_fine::Array{Float64,1},fspace_C_fine::Dict{Symbol,Any},agrid::Array{Float64,1},tol::Float64 = 1e-8)
    # Create the productivity of labor and production
    no_labor_shock = size(z_W)[1];
    no_prod_shock = convert(Int64,n[2]/no_labor_shock);
    θ = z[((convert(Array{Int64,1},s[:,2]).-1) .% no_prod_shock .+1)];
    labor_prod = z_W[convert(Array{Int64,1},floor.((s[:,2]/ no_prod_shock .-0.01) .+ 1))];
    ones_tmp = ones(ns);
    # Workers:
    P_W = ((1 + Q_S)^(1 -ϵ)*ψ_S^ϵ + p_B^(1 -ϵ)*ψ_B^ϵ + p_M^(1 -ϵ)*ψ_M^ϵ)^(1 / (1 - ϵ));
    Y_W = (1 + Q_S) * c̄_S .- w*labor_prod;#was: Y_W = (1 + Q_S) * c̄_S .- w*labor_prod; # consumption is added later as P_W * C
    # Staple farmer [22/03/2022 UPDATED LASZLO:]
    q_S_staples = zeros(ns,3);
    c_S_staples = zeros(ns,3);
    c_B_staples = zeros(ns,3);
    c_M_staples = zeros(ns,3);
    P_S_staples = zeros(ns,3);
    x_S_staples = zeros(ns,3);
    λ_2_S_staples = zeros(ns,3);
    unfeasible_mat = zeros(ns,3);
    Y_S_potential = zeros(ns,3);
    #Case 1: more than enough production and hence no transaction cost:
    x_S_bar_c1 = ((ϕ_S * ζ)/((1 + τ_S) * p_x))^(1/(1 - ζ))*θ.^(1/(1 - ζ));
    x_S_c1 = min.(κ.*s[:,1]./((1 + τ_S) * p_x),x_S_bar_c1);
    π_S_c1 = ϕ_S * θ .*x_S_c1.^ζ - ((1 + τ_S) * p_x) * x_S_c1;
    λ_S_c1 =   ϕ_S .* θ .*ζ./(p_x.*(1 + τ_S)) .* x_S_c1.^(ζ -1) .- 1;
    P_S_c1 = (ψ_S^ϵ + p_B^(1 -ϵ)*ψ_B^ϵ + p_M^(1 -ϵ)*ψ_M^ϵ)^(1 / (1 - ϵ));
    q_S_c1 = ϕ_S *θ .* x_S_c1.^ζ;
    Y_S_c1 = - π_S_c1; #was: Y_S_c1 = c̄_S .- π_S_c1;
    #Case 2: Impossibly to produce enough food:
    x_S_bar_c2 = ((ϕ_S*(1 + Q_S) * ζ)/((1 + τ_S) * p_x))^(1/(1 - ζ))*θ.^(1/(1 - ζ));
    x_S_c2 = min.(κ.*s[:,1]./((1 + τ_S) * p_x),x_S_bar_c2);
    q_S_c2 = ϕ_S * θ .*x_S_c2.^ζ;
    fertilizer_exp = ((1 + τ_S) * p_x) * x_S_c2;
    # λ_S_c2 =   ϕ_S *(1 + Q_S) .* θ .*ζ./p_x .* x_S_c2.^(ζ -1) .- (1 + τ_S); --- KAROL: I THINK THIS IS SLIGHLY WRONG TOO, CORRECTION BELOW (but this is subject to comment on OL):
    λ_S_c2 =    ϕ_S*(1 + Q_S) .* θ .*ζ./(p_x.*(1 + τ_S)) .* x_S_c2.^(ζ -1) .- 1;

    P_S_c2 = ((1 + Q_S)^(1 -ϵ)*ψ_S^ϵ + p_B^(1 -ϵ)*ψ_B^ϵ + p_M^(1 -ϵ)*ψ_M^ϵ)^(1 / (1 - ϵ));
    Y_S_c2 =  fertilizer_exp; # Might have to pay transaction cost even on parts of c̄_S

    #if sum(any.(-a_min .+ q_S_c2 .- fertilizer_exp .- c̄_S .<0))>0 ---Updated cbar_creator
    #    cbar_violated=1
    #else
    #    cbar_violated=0
    #end
    cbar_violated=0
    # handle the c_bar error here. if q_S_c2 - fertilizer_exp  + a_min < bar_c , break.

    #Case 3: Possible to produce enough food, but distorts the production [NOT UPDATED!]

    #c_s_solve_mat = zeros(size(agrid_fine)[1],ns); # Stores the consumption solutions at the double grid.
    # Each column is for different θ, row is for different C. Solution must be for each C
    # Each column is for different θ, row is for different C. Solution must be for each C
    ones_agrid_fine = ones(size(agrid_fine));
    Phi_a_fine = funbase(fspace_C_fine, agrid_fine);
    coeff_λ_2_s = zeros(size(agrid_fine)[1],ns);
    V_temp_staple_residual = zeros(size(agrid_fine)[1],ns);
    C_grid_fine = copy(agrid_fine);
    exit_flag_mat = zeros(size(agrid_fine)[1],ns);
    TC_S_c3_constrained = zeros(ns);
    x_S_c3_constrained = zeros(ns);
    q_S_c3_constrained = zeros(ns);
    c_S_c3_constrained = zeros(ns);

    if cbar_violated==0
        for z_index= 1:(convert(Int64,n[2]/no_labor_shock))

            (coeff_λ_2_s_tmp,V_temp) = goldenx(λ_2_staple_residual,ones_agrid_fine,(1 + Q_S) *ones_agrid_fine,tol,ϕ_S,ζ,τ_S,p_x,
                p_B,p_M,ϕ_B,c̄_S,Q_S,z,ϵ,ψ_S,ψ_B,ψ_M,C_grid_fine,z_index);
            #if sum(V_temp)<tol
            #    println("Problem with the intermediate case") - deleted because we take care of these issues from lines 1225
            #end
            θ_index_start = (z_index-1)*n[1] + 1;
            θ_index_end = z_index*n[1];
            labor_shock_same_index= convert(Int64,ns/no_labor_shock);
            for l_index=0:(no_labor_shock-1)
                coeff_λ_2_s[:,(labor_shock_same_index*l_index+ θ_index_start):(labor_shock_same_index*l_index+θ_index_end)].=coeff_λ_2_s_tmp;
                V_temp_staple_residual[:,(labor_shock_same_index*l_index+ θ_index_start):(labor_shock_same_index*l_index+θ_index_end)].=V_temp;
            end
            for a_index= 1:n[1]
                # Constrained case: bounds
                TC_S_constrained_tmp = κ*agrid[a_index];
                x_S_constrained_tmp = TC_S_constrained_tmp/((1 + τ_S) * p_x)
                q_S_constrained_tmp = ϕ_S * x_S_constrained_tmp.^ζ .*z[z_index];
                condition=(((q_S_constrained_tmp .- c̄_S)./ψ_S.^ϵ./C_grid_fine).^((1 - ϵ)/ϵ) .- ψ_S.^ϵ);
                exitflag_tmp = 0* copy(ones_agrid_fine);
                exitflag_tmp[condition.<0].=-1;
                for l_index=0:(no_labor_shock-1)
                    exit_flag_mat[:,(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))].=exitflag_tmp;
                    TC_S_c3_constrained[(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))] = TC_S_constrained_tmp;
                    x_S_c3_constrained[(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))] = x_S_constrained_tmp;
                    q_S_c3_constrained[(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))] = q_S_constrained_tmp;
                    c_S_c3_constrained[(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))] = q_S_constrained_tmp;
                end
            end
            #coeff_λ_2_s[:,θ_index] = Phi_a_fine\c_s_solve_mat[:,θ_index];
            # for this simple case, this isnt needed as c_s_solve_mat = coeff_λ_2_s!!!
        end
    end #end if cbarviolated

    # Create functions, use staples_objects later.

    # Cash crop farmer
    #π_possible = zeros(ns,3);
    #x_SC_possible  = zeros(ns,3);
    #x_BC_possible = zeros(ns,3);
    #labor_allocated_possible = zeros(ns,3);
    #λ_possible = zeros(ns,3);

    #1) Specialize in cash crop: [24/03/2022 UPDATED LASZLO:]
    x_B_bar_c1 = ((p_B * ϕ_B * ζ)/((1 + τ_B) * p_x))^(1/(1 - ζ))*θ.^(1/(1 - ζ));
    x_B_c1 = min.(κ.*s[:,1]./p_x./(1 + τ_B),x_B_bar_c1);
    q_B_c1 = ϕ_B *θ .* x_B_c1.^ζ;
    π_B_only_B_c1=  p_B * q_B_c1 - ((1 + τ_B) * p_x) * x_B_c1;
    #λ_B_only_B_c1_tmp =   p_B * ϕ_B .* θ .*ζ./p_x .* x_B_c1.^(ζ -1) ./ (1 + τ_B) .- 1;
    λ_B_only_B_c1 =   max.(p_B * ϕ_B .* θ .*ζ.*(κ.*s[:,1]).^(ζ -1) ./ (p_x*(1 + τ_B)).^ζ .- 1,0.0);
    P_B_c1 = ((1 + Q_S)^(1 -ϵ) *ψ_S^ϵ + p_B^(1 -ϵ)*ψ_B^ϵ + p_M^(1 -ϵ)*ψ_M^ϵ)^(1 / (1 - ϵ));
    Y_B_c1 = - π_B_only_B_c1; #was: Y_B_c1 = (1 + Q_S) * c̄_S .- π_B_only_B_c1;

    #2) Specialize in staple crop:
    # Exactly the same problem as before, so just use function staples_objects

    #π_possible[:,1] = π_B_only_B;
    #x_SC_possible[:,1] .= 0;
    #x_BC_possible[:,1] = x_B;
    #labor_allocated_possible[:,1] .= 1.0;
    #λ_possible[:,1] = λ_B_only_B;

    #3) Internal solution - produces both staples and cash crop
    #3a) q_S>c_S, hence no transaction cost paid [24/03/2022 UPDATED LASZLO:]
    labor_allocated_interior_c3a = (p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ))/((p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ)) + (ϕ_S * (1+τ_B)^ζ)^(1/(1-ρ-ζ)))
    x_SC_interior_c3a = ((ϕ_S * ζ * (1 .- labor_allocated_interior_c3a).^ρ)/((1 + τ_S) * p_x))^(1/(1 - ζ))*θ.^(1/(1 - ζ));
    x_BC_interior_c3a = ((p_B * ϕ_B * ζ * (labor_allocated_interior_c3a).^ρ)/((1 + τ_B) * p_x))^(1/(1 - ζ))*θ.^(1/(1 - ζ));
    total_cost_c3a =p_x *  ((1 + τ_S) * x_SC_interior_c3a +  (1 + τ_B) * x_BC_interior_c3a);
    λ_B_interior_c3a = ζ * θ ./ (κ.*s[:,1]).^(1-ζ) .* ( (p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ)) + (ϕ_S * (1+τ_B)^ζ)^(1/(1-ρ-ζ)) ).^(1-ρ-ζ) ./ (p_x*(1+τ_S)*(1+τ_B)).^ζ .- 1;
    λ_B_interior_c3a = max.(λ_B_interior_c3a,0.0);
    Lagrange_factor = (1.0 .+ λ_B_interior_c3a).^(1/(ζ - 1));
    #x_tot_interior_c3a = (1+τ_S).*x_SC_interior_c3a .+ (1+τ_B).*x_BC_interior_c3a
    x_SC_interior_c3a = x_SC_interior_c3a .*Lagrange_factor;
    x_BC_interior_c3a = x_BC_interior_c3a .* Lagrange_factor;
    total_cost_c3a = total_cost_c3a.*Lagrange_factor;
    q_S_c3a = ϕ_S *θ .* x_SC_interior_c3a.^ζ .* (1.0 .- labor_allocated_interior_c3a).^ρ;
    q_B_c3a = ϕ_B *θ .* x_BC_interior_c3a.^ζ .* labor_allocated_interior_c3a.^ρ;
    π_B_interior_c3a = q_S_c3a + p_B * q_B_c3a - total_cost_c3a;
    P_B_c3a =  (ψ_S^ϵ + p_B^(1 -ϵ)*ψ_B^ϵ + p_M^(1 -ϵ)*ψ_M^ϵ)^(1 / (1 - ϵ));
    Y_B_c3a = -π_B_interior_c3a #was: Y_B_c3a = c̄_S.-π_B_interior_c3a


    #3b) q_S<c_S, hence transaction cost is always paid marginally [24/03/2022 UPDATED LASZLO:]
    labor_allocated_interior_c3b = (p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ))/((p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ)) + ((ϕ_S * (1 + Q_S)) * (1+τ_B)^ζ)^(1/(1-ρ-ζ)))
    x_SC_interior_c3b = (((ϕ_S * (1 + Q_S)) * ζ * (1 .- labor_allocated_interior_c3b).^ρ)/((1 + τ_S) * p_x))^(1/(1 - ζ))*θ.^(1/(1 - ζ));
    x_BC_interior_c3b = ((p_B * ϕ_B * ζ * (labor_allocated_interior_c3b).^ρ)/((1 + τ_B) * p_x))^(1/(1 - ζ))*θ.^(1/(1 - ζ));
    total_cost_c3b =p_x *  ((1 + τ_S) * x_SC_interior_c3b +  (1 + τ_B) * x_BC_interior_c3b);
    λ_B_interior_c3b = ζ * θ ./ (κ.*s[:,1]).^(1-ζ) .* ( (p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ)) + ((ϕ_S * (1 + Q_S)) * (1+τ_B)^ζ)^(1/(1-ρ-ζ)) ).^(1-ρ-ζ) ./ (p_x*(1+τ_S)*(1+τ_B)).^ζ .- 1;
    λ_B_interior_c3b = max.(λ_B_interior_c3b,0.0);
    Lagrange_factor = (1.0 .+ λ_B_interior_c3b).^(1/(ζ - 1));
    x_SC_interior_c3b = x_SC_interior_c3b .*Lagrange_factor;
    x_BC_interior_c3b = x_BC_interior_c3b .* Lagrange_factor;
    total_cost_c3b = total_cost_c3b.*Lagrange_factor;
    q_S_c3b = ϕ_S *θ .* x_SC_interior_c3b.^ζ .* (1.0 .- labor_allocated_interior_c3b).^ρ;
    q_B_c3b = ϕ_B *θ .* x_BC_interior_c3b.^ζ .* labor_allocated_interior_c3b.^ρ;
    π_B_interior_c3b = p_B * q_B_c3b - total_cost_c3b
    Y_B_c3b = -π_B_interior_c3b; # market income for the csah crop producer
    P_B_c3b =  ((1 + Q_S)^(1 -ϵ) *ψ_S^ϵ + p_B^(1 -ϵ)*ψ_B^ϵ + p_M^(1 -ϵ)*ψ_M^ϵ)^(1 / (1 - ϵ));

    #3c) q_S=c_S, hence transaction cost is not paid, but production is distorted, see functions above [08/03/2022 NOT FINISHED - KAROL:]
    coeff_λ_2_cashcrop_residual_unconstrained = zeros(size(agrid_fine)[1],ns);
    coeff_λ_2_cashcrop_residual_constrained = zeros(size(agrid_fine)[1],ns);
    V_temp_cashcrop_residual_unconstrained = zeros(size(agrid_fine)[1],ns);
    V_temp_cashcrop_residual_constrained = zeros(size(agrid_fine)[1],ns);
    #    z_index = 1
    #    c_S_cashcrop_residual_unconstrained((c̄_S+tol) * ones_agrid_fine,ϕ_S,ζ,
    #    τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,z,ϵ,ψ_S,ψ_B,ψ_M,C_grid_fine,z_index,ones_agrid_fine,agrid_fine,tol,1)
    for z_index= 1:(convert(Int64,n[2]/no_labor_shock))
         (coeff_λ_2_cashcrop_residual_unconstrained_tmp,V_temp) = goldenx(λ_2_residual_unconstrained,  ones_agrid_fine
         ,(1 + Q_S) *ones_agrid_fine ,tol,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,z,ϵ,ψ_S,ψ_B,ψ_M,C_grid_fine,z_index,ρ);

        θ_index_start = (z_index-1)*n[1] + 1;
        θ_index_end = z_index*n[1];
        labor_shock_same_index= convert(Int64,ns/no_labor_shock);
        for l_index=0:(no_labor_shock-1)
            coeff_λ_2_cashcrop_residual_unconstrained[:,(labor_shock_same_index*l_index+ θ_index_start):(labor_shock_same_index*l_index+θ_index_end)].=coeff_λ_2_cashcrop_residual_unconstrained_tmp;
            V_temp_cashcrop_residual_unconstrained[:,(labor_shock_same_index*l_index+ θ_index_start):(labor_shock_same_index*l_index+θ_index_end)].=V_temp;
        end
        #unconstrained ends, constrained begins:
        for a_index= 1:n[1]
            (coeff_λ_2_cashcrop_residual_constrained_tmp,V_temp) = goldenx(λ_2_residual_constrained,  ones_agrid_fine
            ,(1 + Q_S) *ones_agrid_fine ,tol,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,z,ϵ,ψ_S,ψ_B,ψ_M,C_grid_fine,z_index,agrid,a_index,κ,ρ);
            for l_index=0:(no_labor_shock-1)
                coeff_λ_2_cashcrop_residual_constrained[:,(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))].=coeff_λ_2_cashcrop_residual_constrained_tmp;
                V_temp_cashcrop_residual_constrained[:,(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))].=V_temp;
            end
        end
    end
    # Initialization needed for the functions - better do it once here:
    c_S_mat = zeros(ns,5);
    c_B_mat = zeros(ns,5);
    c_M_mat = zeros(ns,5);
    x_S_mat = zeros(ns,5);
    x_B_mat = zeros(ns,5);
    q_S_mat = zeros(ns,5);
    q_B_mat = zeros(ns,5);
    land_B_mat = zeros(ns,5);
    λ_2_mat = zeros(ns,5);
    P_B_mat = zeros(ns,5);
    feasibility_mat = zeros(ns,5);
    Y_B_mat = zeros(ns,5);
    TC_mat = zeros(ns,5);
    # Initialization needed for the constrained case of 3c functions - better do it once here:

    x_S_mat_3c = zeros(ns,2);
    x_B_mat_3c = zeros(ns,2);
    land_B_mat_3c = zeros(ns,2);
    λ_2_mat_3c = zeros(ns,2);
    TC_mat_3c = zeros(ns,2);
    # Get consumption bounds
    C_grid_fine_mat = repeat(C_grid_fine,1,ns);
    C_max_unconstrained = maximum(C_grid_fine_mat .* (V_temp_cashcrop_residual_unconstrained.> -tol),dims = 1)[:] .+ a_min;
    C_max_constrained =  maximum(C_grid_fine_mat .* (V_temp_cashcrop_residual_constrained.> -tol),dims = 1)[:].+ a_min;
    C_max_staple = maximum(C_grid_fine_mat .* (V_temp_staple_residual.> -tol),dims = 1)[:] .+ a_min;
    C_max_staple_constrained = maximum(C_grid_fine_mat .* (exit_flag_mat.== 0),dims = 1)[:];
    C_min_staple_constrained = a_min * ones_tmp;
    tmp=zeros(Int64,ns);
    tmp1=zeros(Int64,ns);
    tmp2=zeros(Int64,ns);
    for ii=1:ns

        if isnothing(findfirst(myCondition, (C_grid_fine_mat .* (V_temp_cashcrop_residual_unconstrained.> -tol))[:,ii]))
            tmp[ii]=1
        else
            tmp[ii]=findfirst(myCondition, (C_grid_fine_mat .* (V_temp_cashcrop_residual_unconstrained.> -tol))[:,ii])
        end

        if isnothing(findfirst(myCondition, (C_grid_fine_mat .* (V_temp_cashcrop_residual_constrained.> -tol))[:,ii]))
            tmp1[ii]=1
        else
            tmp1[ii]=findfirst(myCondition, (C_grid_fine_mat .* (V_temp_cashcrop_residual_constrained.> -tol))[:,ii])
        end
        if isnothing(findfirst(myCondition, (C_grid_fine_mat .* (V_temp_staple_residual.> -tol))[:,ii]))
            tmp2[ii]=1
        else
            tmp2[ii]=findfirst(myCondition, (C_grid_fine_mat .* (V_temp_staple_residual.> -tol))[:,ii])
        end
    end
    C_min_unconstrained = C_grid_fine[tmp];
    C_min_constrained = C_grid_fine[tmp1];
    C_min_staple = C_grid_fine[tmp2];
    C_max_mat = zeros(ns,2);
    C_min_mat = zeros(ns,2);
    C_max_mat[:,1] =C_max_unconstrained;
    C_max_mat[:,2] =C_max_constrained;
    C_min_mat[:,1] =C_min_unconstrained;
    C_min_mat[:,2] =C_min_constrained;

    return (θ,labor_prod,tol,P_W,Y_W,coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,
            x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
            coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,
            λ_B_interior_c3a,x_SC_interior_c3a,x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,
            x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b, c_S_mat,c_B_mat,
            c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat, λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,q_S_staples,c_S_staples,c_B_staples,
            c_M_staples,P_S_staples,x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,C_max_unconstrained ,C_max_constrained,C_min_unconstrained,C_min_constrained,TC_mat,
            C_max_staple,C_min_staple,C_max_staple_constrained,C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained,cbar_violated,
            x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c)
end
function income_creator_no_approx(s::Array{Float64,2},ns::Int64,
    z::Array{Float64,1},z_W::Array{Float64,1},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,
    p_B::Float64,p_M::Float64,ϕ_B::Float64,τ_B::Float64,ρ::Float64,w::Float64,r::Float64,
    c̄_S::Float64,a_min::Float64,a_max::Float64,γ::Float64,n::Array{Int64,1},κ::Float64,Q_S::Float64,
    ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,
    coeff_λ_2_cashcrop_residual_unconstrained::Array{Float64,2},
    coeff_λ_2_cashcrop_residual_constrained::Array{Float64,2},
    C_max_unconstrained::Array{Float64,1} ,C_max_constrained::Array{Float64,1},C_min_unconstrained::Array{Float64,1},
    C_min_constrained::Array{Float64,1},
    coeff_λ_2_s::Array{Float64,2},agrid_fine::Array{Float64,1},fspace_C_fine::Dict{Symbol,Any},C_max_staple::Array{Float64,1},
    C_min_staple::Array{Float64,1},C_max_staple_constrained::Array{Float64,1},
    C_min_staple_constrained::Array{Float64,1},TC_S_c3_constrained::Array{Float64,1},
    x_S_c3_constrained::Array{Float64,1},q_S_c3_constrained::Array{Float64,1},c_S_c3_constrained::Array{Float64,1},tol::Float64 = 1e-8)
    # Create the productivity of labor and production
    no_labor_shock = size(z_W)[1];
    no_prod_shock = convert(Int64,n[2]/no_labor_shock);
    θ = z[((convert(Array{Int64,1},s[:,2]).-1) .% no_prod_shock .+1)];
    labor_prod = z_W[convert(Array{Int64,1},floor.((s[:,2]/ no_prod_shock .-0.01) .+ 1))];
    ones_tmp = ones(ns);
    # Workers:
    P_W = ((1 + Q_S)^(1 -ϵ)*ψ_S^ϵ + p_B^(1 -ϵ)*ψ_B^ϵ + p_M^(1 -ϵ)*ψ_M^ϵ)^(1 / (1 - ϵ));
    Y_W = (1 + Q_S) * c̄_S .- w*labor_prod; #was: Y_W = (1 + Q_S) * c̄_S .- w*labor_prod; # consumption is added later as P_W * C
    # Staple farmer
    q_S_staples = zeros(ns,3);
    c_S_staples = zeros(ns,3);
    c_B_staples = zeros(ns,3);
    c_M_staples = zeros(ns,3);
    P_S_staples = zeros(ns,3);
    x_S_staples = zeros(ns,3);
    λ_2_S_staples = zeros(ns,3);
    unfeasible_mat = zeros(ns,3);
    Y_S_potential = zeros(ns,3);
    #Case 1: more than enough production and hence no transaction cost:
    x_S_bar_c1 = ((ϕ_S * ζ)/((1 + τ_S) * p_x))^(1/(1 - ζ))*θ.^(1/(1 - ζ));
    x_S_c1 = min.(κ.*s[:,1]./((1 + τ_S) * p_x),x_S_bar_c1);
    π_S_c1 = ϕ_S * θ .*x_S_c1.^ζ - ((1 + τ_S) * p_x) * x_S_c1;
    λ_S_c1 =   ϕ_S .* θ .*ζ./(p_x.*(1 + τ_S)) .* x_S_c1.^(ζ -1) .- 1;
    P_S_c1 = (ψ_S^ϵ + p_B^(1 -ϵ)*ψ_B^ϵ + p_M^(1 -ϵ)*ψ_M^ϵ)^(1 / (1 - ϵ));
    q_S_c1 = ϕ_S *θ .* x_S_c1.^ζ;
    Y_S_c1 = - π_S_c1; # was: Y_S_c1 = c̄_S .- π_S_c1;
    #Case 2: Impossibly to produce enough food:
    x_S_bar_c2 = ((ϕ_S*(1 + Q_S) * ζ)/((1 + τ_S) * p_x))^(1/(1 - ζ))*θ.^(1/(1 - ζ));
    x_S_c2 = min.(κ.*s[:,1]./((1 + τ_S) * p_x),x_S_bar_c2);
    q_S_c2 = ϕ_S * θ .*x_S_c2.^ζ;
    fertilizer_exp = ((1 + τ_S) * p_x) * x_S_c2;
    # λ_S_c2 =   ϕ_S *(1 + Q_S) .* θ .*ζ./p_x .* x_S_c2.^(ζ -1) .- (1 + τ_S); --- KAROL: I THINK THIS IS SLIGHLY WRONG TOO, CORRECTION BELOW (but this is subject to comment on OL):
    λ_S_c2 =    ϕ_S*(1 + Q_S) .* θ .*ζ./(p_x.*(1 + τ_S)) .* x_S_c2.^(ζ -1) .- 1;

    P_S_c2 = ((1 + Q_S)^(1 -ϵ)*ψ_S^ϵ + p_B^(1 -ϵ)*ψ_B^ϵ + p_M^(1 -ϵ)*ψ_M^ϵ)^(1 / (1 - ϵ));
    Y_S_c2 =  fertilizer_exp; # Might have to pay transaction cost even on parts of c̄_S
    #Case 3: Possible to produce enough food, but distorts the production

    #c_s_solve_mat = zeros(size(agrid_fine)[1],ns); # Stores the consumption solutions at the double grid.
    # Each column is for different θ, row is for different C. Solution must be for each C
    ones_agrid_fine = ones(size(agrid_fine));
    C_grid_fine = copy(agrid_fine);
    exit_flag_mat_fine = zeros(size(agrid_fine)[1],ns);
    TC_S_c3_constrained_fine = zeros(ns);
    x_S_c3_constrained_fine  = zeros(ns);
    q_S_c3_constrained_fine  = zeros(ns);
    c_S_c3_constrained_fine  = zeros(ns);
    for z_index= 1:(convert(Int64,n[2]/no_labor_shock))
        labor_shock_same_index= convert(Int64,ns/no_labor_shock);
        for a_index= 1:n[1]
            # Constrained case: bounds
            TC_S_constrained_tmp = κ*agrid_fine[a_index];
            x_S_constrained_tmp = TC_S_constrained_tmp/((1 + τ_S) * p_x)
            q_S_constrained_tmp = ϕ_S * x_S_constrained_tmp.^ζ .*z[z_index];
            condition=(((q_S_constrained_tmp .- c̄_S)./ψ_S.^ϵ./C_grid_fine).^((1 - ϵ)/ϵ) .- ψ_S.^ϵ);
            exitflag_tmp_fine = 0* copy(ones_agrid_fine);
            exitflag_tmp_fine[condition.<0].=-1;
            for l_index=0:(no_labor_shock-1)
                exit_flag_mat_fine[:,(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))].=exitflag_tmp_fine;
                TC_S_c3_constrained_fine[(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))] = TC_S_constrained_tmp;
                x_S_c3_constrained_fine[(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))] = x_S_constrained_tmp;
                q_S_c3_constrained_fine[(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))] = q_S_constrained_tmp;
                c_S_c3_constrained_fine[(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))] = q_S_constrained_tmp;
            end
        end
        #coeff_λ_2_s[:,θ_index] = Phi_a_fine\c_s_solve_mat[:,θ_index];
        # for this simple case, this isnt needed as c_s_solve_mat = coeff_λ_2_s!!!
    end
    C_grid_fine_mat = repeat(C_grid_fine,1,ns);
    C_max_staple_constrained_fine = maximum(C_grid_fine_mat .* (exit_flag_mat_fine.== 0),dims = 1)[:];
    C_min_staple_constrained_fine = a_min * ones_tmp;
    # Create functions, use staples_objects later.

    # Cash crop farmer
    #π_possible = zeros(ns,3);
    #x_SC_possible  = zeros(ns,3);
    #x_BC_possible = zeros(ns,3);
    #labor_allocated_possible = zeros(ns,3);
    #λ_possible = zeros(ns,3);

    #1) Specialize in cash crop: [24/03/2022 UPDATED LASZLO:]
    x_B_bar_c1 = ((p_B * ϕ_B * ζ)/((1 + τ_B) * p_x))^(1/(1 - ζ))*θ.^(1/(1 - ζ));
    x_B_c1 = min.(κ.*s[:,1]./p_x./(1 + τ_B),x_B_bar_c1);
    q_B_c1 = ϕ_B *θ .* x_B_c1.^ζ;
    π_B_only_B_c1=  p_B * q_B_c1 - ((1 + τ_B) * p_x) * x_B_c1;
    #λ_B_only_B_c1_tmp =   p_B * ϕ_B .* θ .*ζ./p_x .* x_B_c1.^(ζ -1) ./ (1 + τ_B) .- 1;
    λ_B_only_B_c1 =   max.(p_B * ϕ_B .* θ .*ζ.*(κ.*s[:,1]).^(ζ -1) ./ (p_x*(1 + τ_B)).^ζ .- 1,0.0);
    P_B_c1 = ((1 + Q_S)^(1 -ϵ) *ψ_S^ϵ + p_B^(1 -ϵ)*ψ_B^ϵ + p_M^(1 -ϵ)*ψ_M^ϵ)^(1 / (1 - ϵ));
    Y_B_c1 = - π_B_only_B_c1; #was: Y_B_c1 = (1 + Q_S) * c̄_S .- π_B_only_B_c1;

    #2) Specialize in staple crop:
    # Exactly the same problem as before, so just use function staples_objects

    #π_possible[:,1] = π_B_only_B;
    #x_SC_possible[:,1] .= 0;
    #x_BC_possible[:,1] = x_B;
    #labor_allocated_possible[:,1] .= 1.0;
    #λ_possible[:,1] = λ_B_only_B;

    #3) Internal solution - produces both staples and cash crop
    #3a) q_S>c_S, hence no transaction cost paid [24/03/2022 UPDATED LASZLO:]
    labor_allocated_interior_c3a = (p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ))/((p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ)) + (ϕ_S * (1+τ_B)^ζ)^(1/(1-ρ-ζ)))
    x_SC_interior_c3a = ((ϕ_S * ζ * (1 .- labor_allocated_interior_c3a).^ρ)/((1 + τ_S) * p_x))^(1/(1 - ζ))*θ.^(1/(1 - ζ));
    x_BC_interior_c3a = ((p_B * ϕ_B * ζ * (labor_allocated_interior_c3a).^ρ)/((1 + τ_B) * p_x))^(1/(1 - ζ))*θ.^(1/(1 - ζ));
    total_cost_c3a =p_x *  ((1 + τ_S) * x_SC_interior_c3a +  (1 + τ_B) * x_BC_interior_c3a);
    λ_B_interior_c3a = ζ * θ ./ (κ.*s[:,1]).^(1-ζ) .* ( (p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ)) + (ϕ_S * (1+τ_B)^ζ)^(1/(1-ρ-ζ)) ).^(1-ρ-ζ) ./ (p_x*(1+τ_S)*(1+τ_B)).^ζ .- 1;
    λ_B_interior_c3a = max.(λ_B_interior_c3a,0.0);
    Lagrange_factor = (1.0 .+ λ_B_interior_c3a).^(1/(ζ - 1));
    #x_tot_interior_c3a = (1+τ_S).*x_SC_interior_c3a .+ (1+τ_B).*x_BC_interior_c3a
    x_SC_interior_c3a = x_SC_interior_c3a .*Lagrange_factor;
    x_BC_interior_c3a = x_BC_interior_c3a .* Lagrange_factor;
    total_cost_c3a = total_cost_c3a.*Lagrange_factor;
    q_S_c3a = ϕ_S *θ .* x_SC_interior_c3a.^ζ .* (1.0 .- labor_allocated_interior_c3a).^ρ;
    q_B_c3a = ϕ_B *θ .* x_BC_interior_c3a.^ζ .* labor_allocated_interior_c3a.^ρ;
    π_B_interior_c3a = q_S_c3a + p_B * q_B_c3a - total_cost_c3a;
    P_B_c3a =  (ψ_S^ϵ + p_B^(1 -ϵ)*ψ_B^ϵ + p_M^(1 -ϵ)*ψ_M^ϵ)^(1 / (1 - ϵ));
    Y_B_c3a = -π_B_interior_c3a #was: Y_B_c3a = c̄_S.-π_B_interior_c3a


    #3b) q_S<c_S, hence transaction cost is always paid marginally [24/03/2022 UPDATED LASZLO:]
    labor_allocated_interior_c3b = (p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ))/((p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ)) + ((ϕ_S * (1 + Q_S)) * (1+τ_B)^ζ)^(1/(1-ρ-ζ)))
    x_SC_interior_c3b = (((ϕ_S * (1 + Q_S)) * ζ * (1 .- labor_allocated_interior_c3b).^ρ)/((1 + τ_S) * p_x))^(1/(1 - ζ))*θ.^(1/(1 - ζ));
    x_BC_interior_c3b = ((p_B * ϕ_B * ζ * (labor_allocated_interior_c3b).^ρ)/((1 + τ_B) * p_x))^(1/(1 - ζ))*θ.^(1/(1 - ζ));
    total_cost_c3b =p_x *  ((1 + τ_S) * x_SC_interior_c3b +  (1 + τ_B) * x_BC_interior_c3b);
    λ_B_interior_c3b = ζ * θ ./ (κ.*s[:,1]).^(1-ζ) .* ( (p_B * ϕ_B * (1+τ_S)^ζ)^(1/(1-ρ-ζ)) + ((ϕ_S * (1 + Q_S)) * (1+τ_B)^ζ)^(1/(1-ρ-ζ)) ).^(1-ρ-ζ) ./ (p_x*(1+τ_S)*(1+τ_B)).^ζ .- 1;
    λ_B_interior_c3b = max.(λ_B_interior_c3b,0.0);
    Lagrange_factor = (1.0 .+ λ_B_interior_c3b).^(1/(ζ - 1));
    x_SC_interior_c3b = x_SC_interior_c3b .*Lagrange_factor;
    x_BC_interior_c3b = x_BC_interior_c3b .* Lagrange_factor;
    total_cost_c3b = total_cost_c3b.*Lagrange_factor;
    q_S_c3b = ϕ_S *θ .* x_SC_interior_c3b.^ζ .* (1.0 .- labor_allocated_interior_c3b).^ρ;
    q_B_c3b = ϕ_B *θ .* x_BC_interior_c3b.^ζ .* labor_allocated_interior_c3b.^ρ;
    π_B_interior_c3b = p_B * q_B_c3b - total_cost_c3b;
    Y_B_c3b = -π_B_interior_c3b; # market income for the csah crop producer
    P_B_c3b =  ((1 + Q_S)^(1 -ϵ) *ψ_S^ϵ + p_B^(1 -ϵ)*ψ_B^ϵ + p_M^(1 -ϵ)*ψ_M^ϵ)^(1 / (1 - ϵ));
    #3c) q_S=c_S, hence transaction cost is not paid, but production is distorted, see functions above
    C_max_mat = zeros(ns,2);
    C_min_mat = zeros(ns,2);
    conversion_to_higher_dim = convert(Int64,ns/size(C_max_unconstrained)[1]);
    C_max_mat[:,1] =kron(C_max_unconstrained,ones(conversion_to_higher_dim));
    #C_max_mat[:,2] =kron(C_max_constrained,ones(conversion_to_higher_dim));
    C_min_mat[:,1] =kron(C_min_unconstrained,ones(conversion_to_higher_dim));
    C_max_staple_fine =kron(C_max_staple,ones(conversion_to_higher_dim));
    #C_max_mat[:,2] =kron(C_max_constrained,ones(conversion_to_higher_dim));
    C_min_staple_fine =kron(C_min_staple,ones(conversion_to_higher_dim));
    #C_min_mat[:,2] =kron(C_min_constrained,ones(conversion_to_higher_dim));
    # Initialization needed for the functions - better do it once here:
    c_S_mat = zeros(ns,5);
    c_B_mat = zeros(ns,5);
    c_M_mat = zeros(ns,5);
    x_S_mat = zeros(ns,5);
    x_B_mat = zeros(ns,5);
    q_S_mat = zeros(ns,5);
    q_B_mat = zeros(ns,5);
    land_B_mat = zeros(ns,5);
    λ_2_mat = zeros(ns,5);
    P_B_mat = zeros(ns,5);
    feasibility_mat = zeros(ns,5);
    Y_B_mat = zeros(ns,5);
    TC_mat = zeros(ns,5);
    # Initialization needed for the constrained case of 3c functions - better do it once here:

    x_S_mat_3c_fine = zeros(ns,2);
    x_B_mat_3c_fine = zeros(ns,2);
    land_B_mat_3c_fine = zeros(ns,2);
    λ_2_mat_3c_fine = zeros(ns,2);
    TC_mat_3c_fine = zeros(ns,2);
    #Calculate the coefficients - the constrained must be rerun as a_grid changes:
    coeff_λ_2_s = kron(coeff_λ_2_s,ones(1,conversion_to_higher_dim));
    coeff_λ_2_cashcrop_residual_unconstrained= kron(coeff_λ_2_cashcrop_residual_unconstrained,ones(1,conversion_to_higher_dim));
    #3c) q_S=c_S, hence transaction cost is not paid, but production is distorted, see functions above
    coeff_λ_2_cashcrop_residual_constrained = zeros(size(agrid_fine)[1],ns);
    V_temp_cashcrop_residual_constrained = zeros(size(agrid_fine)[1],ns);
    for z_index= 1:(convert(Int64,n[2]/no_labor_shock))
        for a_index= 1:n[1]
            (coeff_λ_2_cashcrop_residual_constrained_tmp,V_temp) = goldenx(λ_2_residual_constrained,  ones_agrid_fine
            ,(1 + Q_S) * ones_agrid_fine,tol,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,z,ϵ,ψ_S,ψ_B,ψ_M,C_grid_fine,z_index,agrid_fine,a_index,κ,ρ);
            #V_temp,coeff_x_SC_interior_c3a_constrained_tmp = c_S_cashcrop_residual_constrained(coeff_c_S_cashcrop_residual_constrained_tmp,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,z,ϵ,ψ_S,ψ_B,ψ_M,C_grid_fine,z_index,agrid_fine,a_index,κ,1);
            labor_shock_same_index= convert(Int64,ns/no_labor_shock);
            for l_index=0:(no_labor_shock-1)
                coeff_λ_2_cashcrop_residual_constrained[:,(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))].=coeff_λ_2_cashcrop_residual_constrained_tmp;
                #coeff_x_S_cashcrop_residual_constrained[:,(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))].=coeff_x_SC_interior_c3a_constrained_tmp;
                V_temp_cashcrop_residual_constrained[:,(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))].=V_temp;
            end
        end
    end
    C_max_constrained =  maximum(C_grid_fine_mat .* (V_temp_cashcrop_residual_constrained.> -tol),dims = 1)[:];
    #C_min_unconstrained = minimum(C_grid_fine_mat .* (V_temp_cashcrop_residual_unconstrained.> -tol)+ 1000*(V_temp_cashcrop_residual_unconstrained.< -tol) ,dims = 1  )[:];
    tmp=zeros(Int64,ns);
    for ii=1:ns
        if isnothing(findfirst(myCondition, (C_grid_fine_mat .* (V_temp_cashcrop_residual_constrained.> -tol))[:,ii]))
            tmp[ii]=1
        else
            tmp[ii]=findfirst(myCondition, (C_grid_fine_mat .* (V_temp_cashcrop_residual_constrained.> -tol))[:,ii])
        end
    end
    C_min_constrained = C_grid_fine[tmp];
    #C_min_constrained =  minimum(C_grid_fine_mat .* (V_temp_cashcrop_residual_constrained.> -tol)+ 1000*(V_temp_cashcrop_residual_constrained.< -tol),dims = 1)[:];
    C_max_mat[:,2] =C_max_constrained
    C_min_mat[:,2] =C_min_constrained

    #C = s_fine[:,1]/2 .+ a_min;
    #    @btime Y_S,c_S_S,c_B_S,c_M_S,q_S_S,P_S,x_S_S,solve_staple_index_S,λ_2_S = staples_objects(C,ϕ_S,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ζ,c_s_residual,
    #        c_s_constrained,c_s_approx_staples,θ_fine,coeff_λ_2_s_fine,fspace_a_fine,staples_c3_objects,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,s_fine,q_S_c1_fine,
    #        q_S_c2_fine,q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine);
    #C = s_fine[:,1]/2 .+ a_min;
    #Y_S,c_S_S,c_B_S,c_M_S,q_S_S,P_S,x_S_S,solve_staple_index_S,λ_2_S = staples_objects(C,ϕ_S,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ζ,c_s_residual,
    #        c_s_constrained,c_s_approx_staples,θ_fine,coeff_λ_2_s_fine,fspace_a_fine,staples_c3_objects,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine,
    #        x_S_c2_fine,s_fine,q_S_c1_fine,q_S_c2_fine,q_S_staples_fine,
    #        c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine);

    #@btime  (c_S_B,c_B_C,c_M_B,x_SC,x_BC,land_C,λ_2,P_B,Y_B,q_S_C,q_B_B,solve_cash_crop_index_B,solve_staple_index_B) = cashcrop_objects(C,coeff_c_S_cashcrop_residual_unconstrained_fine,
    #        coeff_x_S_cashcrop_residual_unconstrained_fine,coeff_c_S_cashcrop_residual_constrained_fine,θ_fine,
    #         fspace_a_fine,c_S_approx_cashcrop,x_S_approx_cashcrop,s_fine,ϕ_S,ζ,τ_S,p_x,p_B,
    #             p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,
    #             a_min,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
    #             c_s_residual,c_s_constrained,c_s_approx_staples,coeff_λ_2_s_fine,
    #             staples_c3_objects,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,
    #             y_c3a_fine,labor_allocated_interior_c3a_fine,λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,
    #             x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,q_S_c2_fine,
    #             q_B_c1_fine,q_S_c3a_fine,q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,
    #             x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,Y_B_c3b_fine,
    #             c_S_mat_fine,c_B_mat_fine,
    #             c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,q_B_mat_fine,land_B_mat_fine,
    #             λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,feasibility_mat_fine,C_max_mat_fine,C_min_mat_fine,
    #             y_bar_c3c_constraint_mat_fine,
    #             q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,
    #             x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine)

    return (θ,labor_prod,tol,P_W,Y_W,coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,
            x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
            coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,
            λ_B_interior_c3a,x_SC_interior_c3a,x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,
            x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b, c_S_mat,c_B_mat,
            c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat, λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,q_S_staples,c_S_staples,c_B_staples,
            c_M_staples,P_S_staples,x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,TC_mat,C_max_staple_fine,C_min_staple_fine,C_max_staple_constrained_fine,
            C_min_staple_constrained_fine,TC_S_c3_constrained_fine,x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine,
            x_S_mat_3c_fine,x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine)
end

function income_creator_ext(s::Array{Float64,2}, ns::Int64,
    z::Array{Float64,1}, z_W::Array{Float64,1}, ϕ_S::Float64, ζ::Float64, τ_S::Float64, p_x::Float64,
    p_B::Float64, p_M::Float64, ϕ_B::Float64, τ_B::Float64, ρ::Float64, w::Float64, r::Float64,
    c̄_S::Float64, a_min::Float64, a_max::Float64, γ::Float64, n::Array{Int64,1}, κ::Float64, Q_S::Float64,
    ϵ::Float64, ψ_S::Float64, ψ_B::Float64, ψ_M::Float64, agrid_fine::Array{Float64,1}, fspace_C_fine::Dict{Symbol,Any}, agrid::Array{Float64,1}, epsilon_u::Float64, cttilde::Float64, ctilde::Float64, tol::Float64=1e-8)
    # Create the productivity of labor and production
    no_labor_shock = size(z_W)[1]
    no_prod_shock = convert(Int64, n[2] / no_labor_shock)
    θ = z[((convert(Array{Int64,1}, s[:, 2]).-1).%no_prod_shock.+1)]
    labor_penalty_ext = exp(-epsilon_u * (cttilde - 0.2))
    labor_prod = z_W[convert(Array{Int64,1}, floor.((s[:, 2] / no_prod_shock .- 0.01) .+ 1))]*labor_penalty_ext
    ones_tmp = ones(ns)
    # Workers:
    P_W = ((1 + Q_S)^(1 - ϵ) * ψ_S^ϵ + p_B^(1 - ϵ) * ψ_B^ϵ + p_M^(1 - ϵ) * ψ_M^ϵ)^(1 / (1 - ϵ))
    Y_W = (1 + Q_S) * c̄_S .- w * labor_prod#was: Y_W = (1 + Q_S) * c̄_S .- w*labor_prod; # consumption is added later as P_W * C
    # Staple farmer [22/03/2022 UPDATED LASZLO:]
    q_S_staples = zeros(ns, 3)
    c_S_staples = zeros(ns, 3)
    c_B_staples = zeros(ns, 3)
    c_M_staples = zeros(ns, 3)
    P_S_staples = zeros(ns, 3)
    x_S_staples = zeros(ns, 3)
    λ_2_S_staples = zeros(ns, 3)
    unfeasible_mat = zeros(ns, 3)
    Y_S_potential = zeros(ns, 3)
    #Case 1: more than enough production and hence no transaction cost:
    x_S_bar_c1 = ((ϕ_S * ζ) / ((1 + τ_S) * p_x))^(1 / (1 - ζ)) * θ .^ (1 / (1 - ζ))
    x_S_c1 = min.(κ .* s[:, 1] ./ ((1 + τ_S) * p_x), x_S_bar_c1)
    π_S_c1 = ϕ_S * θ .* x_S_c1 .^ ζ - ((1 + τ_S) * p_x) * x_S_c1
    λ_S_c1 = ϕ_S .* θ .* ζ ./ (p_x .* (1 + τ_S)) .* x_S_c1 .^ (ζ - 1) .- 1
    P_S_c1 = (ψ_S^ϵ + p_B^(1 - ϵ) * ψ_B^ϵ + p_M^(1 - ϵ) * ψ_M^ϵ)^(1 / (1 - ϵ))
    q_S_c1 = ϕ_S * θ .* x_S_c1 .^ ζ
    Y_S_c1 = -π_S_c1 #was: Y_S_c1 = c̄_S .- π_S_c1;
    #Case 2: Impossibly to produce enough food:
    x_S_bar_c2 = ((ϕ_S * (1 + Q_S) * ζ) / ((1 + τ_S) * p_x))^(1 / (1 - ζ)) * θ .^ (1 / (1 - ζ))
    x_S_c2 = min.(κ .* s[:, 1] ./ ((1 + τ_S) * p_x), x_S_bar_c2)
    q_S_c2 = ϕ_S * θ .* x_S_c2 .^ ζ
    fertilizer_exp = ((1 + τ_S) * p_x) * x_S_c2
    # λ_S_c2 =   ϕ_S *(1 + Q_S) .* θ .*ζ./p_x .* x_S_c2.^(ζ -1) .- (1 + τ_S); --- KAROL: I THINK THIS IS SLIGHLY WRONG TOO, CORRECTION BELOW (but this is subject to comment on OL):
    λ_S_c2 = ϕ_S * (1 + Q_S) .* θ .* ζ ./ (p_x .* (1 + τ_S)) .* x_S_c2 .^ (ζ - 1) .- 1

    P_S_c2 = ((1 + Q_S)^(1 - ϵ) * ψ_S^ϵ + p_B^(1 - ϵ) * ψ_B^ϵ + p_M^(1 - ϵ) * ψ_M^ϵ)^(1 / (1 - ϵ))
    Y_S_c2 = fertilizer_exp # Might have to pay transaction cost even on parts of c̄_S

    #if sum(any.(-a_min .+ q_S_c2 .- fertilizer_exp .- c̄_S .<0))>0 ---Updated cbar_creator
    #    cbar_violated=1
    #else
    #    cbar_violated=0
    #end
    cbar_violated = 0
    # handle the c_bar error here. if q_S_c2 - fertilizer_exp  + a_min < bar_c , break.

    #Case 3: Possible to produce enough food, but distorts the production [NOT UPDATED!]

    #c_s_solve_mat = zeros(size(agrid_fine)[1],ns); # Stores the consumption solutions at the double grid.
    # Each column is for different θ, row is for different C. Solution must be for each C
    # Each column is for different θ, row is for different C. Solution must be for each C
    ones_agrid_fine = ones(size(agrid_fine))
    Phi_a_fine = funbase(fspace_C_fine, agrid_fine)
    coeff_λ_2_s = zeros(size(agrid_fine)[1], ns)
    V_temp_staple_residual = zeros(size(agrid_fine)[1], ns)
    C_grid_fine = copy(agrid_fine)
    exit_flag_mat = zeros(size(agrid_fine)[1], ns)
    TC_S_c3_constrained = zeros(ns)
    x_S_c3_constrained = zeros(ns)
    q_S_c3_constrained = zeros(ns)
    c_S_c3_constrained = zeros(ns)

    if cbar_violated == 0
        for z_index = 1:(convert(Int64, n[2] / no_labor_shock))

            (coeff_λ_2_s_tmp, V_temp) = goldenx(λ_2_staple_residual, ones_agrid_fine, (1 + Q_S) * ones_agrid_fine, tol, ϕ_S, ζ, τ_S, p_x,
                p_B, p_M, ϕ_B, c̄_S, Q_S, z, ϵ, ψ_S, ψ_B, ψ_M, C_grid_fine, z_index)
            #if sum(V_temp)<tol
            #    println("Problem with the intermediate case") - deleted because we take care of these issues from lines 1225
            #end
            θ_index_start = (z_index - 1) * n[1] + 1
            θ_index_end = z_index * n[1]
            labor_shock_same_index = convert(Int64, ns / no_labor_shock)
            for l_index = 0:(no_labor_shock-1)
                coeff_λ_2_s[:, (labor_shock_same_index*l_index+θ_index_start):(labor_shock_same_index*l_index+θ_index_end)] .= coeff_λ_2_s_tmp
                V_temp_staple_residual[:, (labor_shock_same_index*l_index+θ_index_start):(labor_shock_same_index*l_index+θ_index_end)] .= V_temp
            end
            for a_index = 1:n[1]
                # Constrained case: bounds
                TC_S_constrained_tmp = κ * agrid[a_index]
                x_S_constrained_tmp = TC_S_constrained_tmp / ((1 + τ_S) * p_x)
                q_S_constrained_tmp = ϕ_S * x_S_constrained_tmp .^ ζ .* z[z_index]
                condition = (((q_S_constrained_tmp .- c̄_S) ./ ψ_S .^ ϵ ./ C_grid_fine) .^ ((1 - ϵ) / ϵ) .- ψ_S .^ ϵ)
                exitflag_tmp = 0 * copy(ones_agrid_fine)
                exitflag_tmp[condition.<0] .= -1
                for l_index = 0:(no_labor_shock-1)
                    exit_flag_mat[:, (labor_shock_same_index*l_index+a_index+((z_index-1)*n[1]))] .= exitflag_tmp
                    TC_S_c3_constrained[(labor_shock_same_index*l_index+a_index+((z_index-1)*n[1]))] = TC_S_constrained_tmp
                    x_S_c3_constrained[(labor_shock_same_index*l_index+a_index+((z_index-1)*n[1]))] = x_S_constrained_tmp
                    q_S_c3_constrained[(labor_shock_same_index*l_index+a_index+((z_index-1)*n[1]))] = q_S_constrained_tmp
                    c_S_c3_constrained[(labor_shock_same_index*l_index+a_index+((z_index-1)*n[1]))] = q_S_constrained_tmp
                end
            end
            #coeff_λ_2_s[:,θ_index] = Phi_a_fine\c_s_solve_mat[:,θ_index];
            # for this simple case, this isnt needed as c_s_solve_mat = coeff_λ_2_s!!!
        end
    end #end if cbarviolated

    # Create functions, use staples_objects later.

    # Cash crop farmer
    #π_possible = zeros(ns,3);
    #x_SC_possible  = zeros(ns,3);
    #x_BC_possible = zeros(ns,3);
    #labor_allocated_possible = zeros(ns,3);
    #λ_possible = zeros(ns,3);

    #1) Specialize in cash crop: [24/03/2022 UPDATED LASZLO:]
    x_B_bar_c1 = ((p_B * ϕ_B * ζ) / ((1 + τ_B) * p_x))^(1 / (1 - ζ)) * θ .^ (1 / (1 - ζ))
    x_B_c1 = min.(κ .* s[:, 1] ./ p_x ./ (1 + τ_B), x_B_bar_c1)
    q_B_c1 = ϕ_B * θ .* x_B_c1 .^ ζ
    π_B_only_B_c1 = p_B * q_B_c1 - ((1 + τ_B) * p_x) * x_B_c1
    #λ_B_only_B_c1_tmp =   p_B * ϕ_B .* θ .*ζ./p_x .* x_B_c1.^(ζ -1) ./ (1 + τ_B) .- 1;
    λ_B_only_B_c1 = max.(p_B * ϕ_B .* θ .* ζ .* (κ .* s[:, 1]) .^ (ζ - 1) ./ (p_x * (1 + τ_B)) .^ ζ .- 1, 0.0)
    P_B_c1 = ((1 + Q_S)^(1 - ϵ) * ψ_S^ϵ + p_B^(1 - ϵ) * ψ_B^ϵ + p_M^(1 - ϵ) * ψ_M^ϵ)^(1 / (1 - ϵ))
    Y_B_c1 = -π_B_only_B_c1 #was: Y_B_c1 = (1 + Q_S) * c̄_S .- π_B_only_B_c1;

    #2) Specialize in staple crop:
    # Exactly the same problem as before, so just use function staples_objects

    #π_possible[:,1] = π_B_only_B;
    #x_SC_possible[:,1] .= 0;
    #x_BC_possible[:,1] = x_B;
    #labor_allocated_possible[:,1] .= 1.0;
    #λ_possible[:,1] = λ_B_only_B;

    #3) Internal solution - produces both staples and cash crop
    #3a) q_S>c_S, hence no transaction cost paid [24/03/2022 UPDATED LASZLO:]
    labor_allocated_interior_c3a = (p_B * ϕ_B * (1 + τ_S)^ζ)^(1 / (1 - ρ - ζ)) / ((p_B * ϕ_B * (1 + τ_S)^ζ)^(1 / (1 - ρ - ζ)) + (ϕ_S * (1 + τ_B)^ζ)^(1 / (1 - ρ - ζ)))
    x_SC_interior_c3a = ((ϕ_S * ζ * (1 .- labor_allocated_interior_c3a) .^ ρ) / ((1 + τ_S) * p_x))^(1 / (1 - ζ)) * θ .^ (1 / (1 - ζ))
    x_BC_interior_c3a = ((p_B * ϕ_B * ζ * (labor_allocated_interior_c3a) .^ ρ) / ((1 + τ_B) * p_x))^(1 / (1 - ζ)) * θ .^ (1 / (1 - ζ))
    total_cost_c3a = p_x * ((1 + τ_S) * x_SC_interior_c3a + (1 + τ_B) * x_BC_interior_c3a)
    λ_B_interior_c3a = ζ * θ ./ (κ .* s[:, 1]) .^ (1 - ζ) .* ((p_B * ϕ_B * (1 + τ_S)^ζ)^(1 / (1 - ρ - ζ)) + (ϕ_S * (1 + τ_B)^ζ)^(1 / (1 - ρ - ζ))) .^ (1 - ρ - ζ) ./ (p_x * (1 + τ_S) * (1 + τ_B)) .^ ζ .- 1
    λ_B_interior_c3a = max.(λ_B_interior_c3a, 0.0)
    Lagrange_factor = (1.0 .+ λ_B_interior_c3a) .^ (1 / (ζ - 1))
    #x_tot_interior_c3a = (1+τ_S).*x_SC_interior_c3a .+ (1+τ_B).*x_BC_interior_c3a
    x_SC_interior_c3a = x_SC_interior_c3a .* Lagrange_factor
    x_BC_interior_c3a = x_BC_interior_c3a .* Lagrange_factor
    total_cost_c3a = total_cost_c3a .* Lagrange_factor
    q_S_c3a = ϕ_S * θ .* x_SC_interior_c3a .^ ζ .* (1.0 .- labor_allocated_interior_c3a) .^ ρ
    q_B_c3a = ϕ_B * θ .* x_BC_interior_c3a .^ ζ .* labor_allocated_interior_c3a .^ ρ
    π_B_interior_c3a = q_S_c3a + p_B * q_B_c3a - total_cost_c3a
    P_B_c3a = (ψ_S^ϵ + p_B^(1 - ϵ) * ψ_B^ϵ + p_M^(1 - ϵ) * ψ_M^ϵ)^(1 / (1 - ϵ))
    Y_B_c3a = -π_B_interior_c3a #was: Y_B_c3a = c̄_S.-π_B_interior_c3a


    #3b) q_S<c_S, hence transaction cost is always paid marginally [24/03/2022 UPDATED LASZLO:]
    labor_allocated_interior_c3b = (p_B * ϕ_B * (1 + τ_S)^ζ)^(1 / (1 - ρ - ζ)) / ((p_B * ϕ_B * (1 + τ_S)^ζ)^(1 / (1 - ρ - ζ)) + ((ϕ_S * (1 + Q_S)) * (1 + τ_B)^ζ)^(1 / (1 - ρ - ζ)))
    x_SC_interior_c3b = (((ϕ_S * (1 + Q_S)) * ζ * (1 .- labor_allocated_interior_c3b) .^ ρ) / ((1 + τ_S) * p_x))^(1 / (1 - ζ)) * θ .^ (1 / (1 - ζ))
    x_BC_interior_c3b = ((p_B * ϕ_B * ζ * (labor_allocated_interior_c3b) .^ ρ) / ((1 + τ_B) * p_x))^(1 / (1 - ζ)) * θ .^ (1 / (1 - ζ))
    total_cost_c3b = p_x * ((1 + τ_S) * x_SC_interior_c3b + (1 + τ_B) * x_BC_interior_c3b)
    λ_B_interior_c3b = ζ * θ ./ (κ .* s[:, 1]) .^ (1 - ζ) .* ((p_B * ϕ_B * (1 + τ_S)^ζ)^(1 / (1 - ρ - ζ)) + ((ϕ_S * (1 + Q_S)) * (1 + τ_B)^ζ)^(1 / (1 - ρ - ζ))) .^ (1 - ρ - ζ) ./ (p_x * (1 + τ_S) * (1 + τ_B)) .^ ζ .- 1
    λ_B_interior_c3b = max.(λ_B_interior_c3b, 0.0)
    Lagrange_factor = (1.0 .+ λ_B_interior_c3b) .^ (1 / (ζ - 1))
    x_SC_interior_c3b = x_SC_interior_c3b .* Lagrange_factor
    x_BC_interior_c3b = x_BC_interior_c3b .* Lagrange_factor
    total_cost_c3b = total_cost_c3b .* Lagrange_factor
    q_S_c3b = ϕ_S * θ .* x_SC_interior_c3b .^ ζ .* (1.0 .- labor_allocated_interior_c3b) .^ ρ
    q_B_c3b = ϕ_B * θ .* x_BC_interior_c3b .^ ζ .* labor_allocated_interior_c3b .^ ρ
    π_B_interior_c3b = p_B * q_B_c3b - total_cost_c3b
    Y_B_c3b = -π_B_interior_c3b # market income for the csah crop producer
    P_B_c3b = ((1 + Q_S)^(1 - ϵ) * ψ_S^ϵ + p_B^(1 - ϵ) * ψ_B^ϵ + p_M^(1 - ϵ) * ψ_M^ϵ)^(1 / (1 - ϵ))

    #3c) q_S=c_S, hence transaction cost is not paid, but production is distorted, see functions above [08/03/2022 NOT FINISHED - KAROL:]
    coeff_λ_2_cashcrop_residual_unconstrained = zeros(size(agrid_fine)[1], ns)
    coeff_λ_2_cashcrop_residual_constrained = zeros(size(agrid_fine)[1], ns)
    V_temp_cashcrop_residual_unconstrained = zeros(size(agrid_fine)[1], ns)
    V_temp_cashcrop_residual_constrained = zeros(size(agrid_fine)[1], ns)
    #    z_index = 1
    #    c_S_cashcrop_residual_unconstrained((c̄_S+tol) * ones_agrid_fine,ϕ_S,ζ,
    #    τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,z,ϵ,ψ_S,ψ_B,ψ_M,C_grid_fine,z_index,ones_agrid_fine,agrid_fine,tol,1)
    for z_index = 1:(convert(Int64, n[2] / no_labor_shock))
        (coeff_λ_2_cashcrop_residual_unconstrained_tmp, V_temp) = goldenx(λ_2_residual_unconstrained, ones_agrid_fine, (1 + Q_S) * ones_agrid_fine, tol, ϕ_S, ζ, τ_S, p_x, p_B, p_M, ϕ_B, τ_B, c̄_S, Q_S, z, ϵ, ψ_S, ψ_B, ψ_M, C_grid_fine, z_index, ρ)

        θ_index_start = (z_index - 1) * n[1] + 1
        θ_index_end = z_index * n[1]
        labor_shock_same_index = convert(Int64, ns / no_labor_shock)
        for l_index = 0:(no_labor_shock-1)
            coeff_λ_2_cashcrop_residual_unconstrained[:, (labor_shock_same_index*l_index+θ_index_start):(labor_shock_same_index*l_index+θ_index_end)] .= coeff_λ_2_cashcrop_residual_unconstrained_tmp
            V_temp_cashcrop_residual_unconstrained[:, (labor_shock_same_index*l_index+θ_index_start):(labor_shock_same_index*l_index+θ_index_end)] .= V_temp
        end
        #unconstrained ends, constrained begins:
        for a_index = 1:n[1]
            (coeff_λ_2_cashcrop_residual_constrained_tmp, V_temp) = goldenx(λ_2_residual_constrained, ones_agrid_fine, (1 + Q_S) * ones_agrid_fine, tol, ϕ_S, ζ, τ_S, p_x, p_B, p_M, ϕ_B, τ_B, c̄_S, Q_S, z, ϵ, ψ_S, ψ_B, ψ_M, C_grid_fine, z_index, agrid, a_index, κ, ρ)
            for l_index = 0:(no_labor_shock-1)
                coeff_λ_2_cashcrop_residual_constrained[:, (labor_shock_same_index*l_index+a_index+((z_index-1)*n[1]))] .= coeff_λ_2_cashcrop_residual_constrained_tmp
                V_temp_cashcrop_residual_constrained[:, (labor_shock_same_index*l_index+a_index+((z_index-1)*n[1]))] .= V_temp
            end
        end
    end
    # Initialization needed for the functions - better do it once here:
    c_S_mat = zeros(ns, 5)
    c_B_mat = zeros(ns, 5)
    c_M_mat = zeros(ns, 5)
    x_S_mat = zeros(ns, 5)
    x_B_mat = zeros(ns, 5)
    q_S_mat = zeros(ns, 5)
    q_B_mat = zeros(ns, 5)
    land_B_mat = zeros(ns, 5)
    λ_2_mat = zeros(ns, 5)
    P_B_mat = zeros(ns, 5)
    feasibility_mat = zeros(ns, 5)
    Y_B_mat = zeros(ns, 5)
    TC_mat = zeros(ns, 5)
    # Initialization needed for the constrained case of 3c functions - better do it once here:

    x_S_mat_3c = zeros(ns, 2)
    x_B_mat_3c = zeros(ns, 2)
    land_B_mat_3c = zeros(ns, 2)
    λ_2_mat_3c = zeros(ns, 2)
    TC_mat_3c = zeros(ns, 2)
    # Get consumption bounds
    C_grid_fine_mat = repeat(C_grid_fine, 1, ns)
    C_max_unconstrained = maximum(C_grid_fine_mat .* (V_temp_cashcrop_residual_unconstrained .> -tol), dims=1)[:] .+ a_min
    C_max_constrained = maximum(C_grid_fine_mat .* (V_temp_cashcrop_residual_constrained .> -tol), dims=1)[:] .+ a_min
    C_max_staple = maximum(C_grid_fine_mat .* (V_temp_staple_residual .> -tol), dims=1)[:] .+ a_min
    C_max_staple_constrained = maximum(C_grid_fine_mat .* (exit_flag_mat .== 0), dims=1)[:]
    C_min_staple_constrained = a_min * ones_tmp
    tmp = zeros(Int64, ns)
    tmp1 = zeros(Int64, ns)
    tmp2 = zeros(Int64, ns)
    for ii = 1:ns

        if isnothing(findfirst(myCondition, (C_grid_fine_mat.*(V_temp_cashcrop_residual_unconstrained.>-tol))[:, ii]))
            tmp[ii] = 1
        else
            tmp[ii] = findfirst(myCondition, (C_grid_fine_mat.*(V_temp_cashcrop_residual_unconstrained.>-tol))[:, ii])
        end

        if isnothing(findfirst(myCondition, (C_grid_fine_mat.*(V_temp_cashcrop_residual_constrained.>-tol))[:, ii]))
            tmp1[ii] = 1
        else
            tmp1[ii] = findfirst(myCondition, (C_grid_fine_mat.*(V_temp_cashcrop_residual_constrained.>-tol))[:, ii])
        end
        if isnothing(findfirst(myCondition, (C_grid_fine_mat.*(V_temp_staple_residual.>-tol))[:, ii]))
            tmp2[ii] = 1
        else
            tmp2[ii] = findfirst(myCondition, (C_grid_fine_mat.*(V_temp_staple_residual.>-tol))[:, ii])
        end
    end
    C_min_unconstrained = C_grid_fine[tmp]
    C_min_constrained = C_grid_fine[tmp1]
    C_min_staple = C_grid_fine[tmp2]
    C_max_mat = zeros(ns, 2)
    C_min_mat = zeros(ns, 2)
    C_max_mat[:, 1] = C_max_unconstrained
    C_max_mat[:, 2] = C_max_constrained
    C_min_mat[:, 1] = C_min_unconstrained
    C_min_mat[:, 2] = C_min_constrained

    return (θ, labor_prod, tol, P_W, Y_W, coeff_λ_2_cashcrop_residual_unconstrained, coeff_λ_2_cashcrop_residual_constrained,
        x_B_c1, π_B_only_B_c1, λ_B_only_B_c1, P_B_c1, Y_B_c1,
        coeff_λ_2_s, P_S_c1, P_S_c2, Y_S_c1, Y_S_c2, x_S_c1, x_S_c2, labor_allocated_interior_c3a,
        λ_B_interior_c3a, x_SC_interior_c3a, x_BC_interior_c3a, Y_B_c3a, P_B_c3a, P_B_c3b, q_S_c1, q_S_c2, q_B_c1, q_S_c3a, q_B_c3a, q_S_c3b, q_B_c3b,
        x_SC_interior_c3b, x_BC_interior_c3b, labor_allocated_interior_c3b, Y_B_c3b, c_S_mat, c_B_mat,
        c_M_mat, x_S_mat, x_B_mat, q_S_mat, q_B_mat, land_B_mat, λ_2_mat, P_B_mat, Y_B_mat, feasibility_mat, C_max_mat, C_min_mat, q_S_staples, c_S_staples, c_B_staples,
        c_M_staples, P_S_staples, x_S_staples, λ_2_S_staples, unfeasible_mat, Y_S_potential, C_max_unconstrained, C_max_constrained, C_min_unconstrained, C_min_constrained, TC_mat,
        C_max_staple, C_min_staple, C_max_staple_constrained, C_min_staple_constrained, TC_S_c3_constrained, x_S_c3_constrained, q_S_c3_constrained, c_S_c3_constrained, cbar_violated,
        x_S_mat_3c, x_B_mat_3c, land_B_mat_3c, λ_2_mat_3c, TC_mat_3c)
end
function income_creator_no_approx_ext(s::Array{Float64,2}, ns::Int64,
    z::Array{Float64,1}, z_W::Array{Float64,1}, ϕ_S::Float64, ζ::Float64, τ_S::Float64, p_x::Float64,
    p_B::Float64, p_M::Float64, ϕ_B::Float64, τ_B::Float64, ρ::Float64, w::Float64, r::Float64,
    c̄_S::Float64, a_min::Float64, a_max::Float64, γ::Float64, n::Array{Int64,1}, κ::Float64, Q_S::Float64,
    ϵ::Float64, ψ_S::Float64, ψ_B::Float64, ψ_M::Float64,
    coeff_λ_2_cashcrop_residual_unconstrained::Array{Float64,2},
    coeff_λ_2_cashcrop_residual_constrained::Array{Float64,2},
    C_max_unconstrained::Array{Float64,1}, C_max_constrained::Array{Float64,1}, C_min_unconstrained::Array{Float64,1},
    C_min_constrained::Array{Float64,1},
    coeff_λ_2_s::Array{Float64,2}, agrid_fine::Array{Float64,1}, fspace_C_fine::Dict{Symbol,Any}, C_max_staple::Array{Float64,1},
    C_min_staple::Array{Float64,1}, C_max_staple_constrained::Array{Float64,1},
    C_min_staple_constrained::Array{Float64,1}, TC_S_c3_constrained::Array{Float64,1},
    x_S_c3_constrained::Array{Float64,1}, q_S_c3_constrained::Array{Float64,1}, c_S_c3_constrained::Array{Float64,1}, epsilon_u::Float64, cttilde::Float64, ctilde::Float64, tol::Float64=1e-8)
    # Create the productivity of labor and production
    no_labor_shock = size(z_W)[1]
    no_prod_shock = convert(Int64, n[2] / no_labor_shock)
    θ = z[((convert(Array{Int64,1}, s[:, 2]).-1).%no_prod_shock.+1)]
    labor_penalty_ext = exp(-epsilon_u * (cttilde - 0.2))
    labor_prod = z_W[convert(Array{Int64,1}, floor.((s[:, 2] / no_prod_shock .- 0.01) .+ 1))] * labor_penalty_ext
    ones_tmp = ones(ns)
    # Workers:
    P_W = ((1 + Q_S)^(1 - ϵ) * ψ_S^ϵ + p_B^(1 - ϵ) * ψ_B^ϵ + p_M^(1 - ϵ) * ψ_M^ϵ)^(1 / (1 - ϵ))
    Y_W = (1 + Q_S) * c̄_S .- w * labor_prod #was: Y_W = (1 + Q_S) * c̄_S .- w*labor_prod; # consumption is added later as P_W * C
    # Staple farmer
    q_S_staples = zeros(ns, 3)
    c_S_staples = zeros(ns, 3)
    c_B_staples = zeros(ns, 3)
    c_M_staples = zeros(ns, 3)
    P_S_staples = zeros(ns, 3)
    x_S_staples = zeros(ns, 3)
    λ_2_S_staples = zeros(ns, 3)
    unfeasible_mat = zeros(ns, 3)
    Y_S_potential = zeros(ns, 3)
    #Case 1: more than enough production and hence no transaction cost:
    x_S_bar_c1 = ((ϕ_S * ζ) / ((1 + τ_S) * p_x))^(1 / (1 - ζ)) * θ .^ (1 / (1 - ζ))
    x_S_c1 = min.(κ .* s[:, 1] ./ ((1 + τ_S) * p_x), x_S_bar_c1)
    π_S_c1 = ϕ_S * θ .* x_S_c1 .^ ζ - ((1 + τ_S) * p_x) * x_S_c1
    λ_S_c1 = ϕ_S .* θ .* ζ ./ (p_x .* (1 + τ_S)) .* x_S_c1 .^ (ζ - 1) .- 1
    P_S_c1 = (ψ_S^ϵ + p_B^(1 - ϵ) * ψ_B^ϵ + p_M^(1 - ϵ) * ψ_M^ϵ)^(1 / (1 - ϵ))
    q_S_c1 = ϕ_S * θ .* x_S_c1 .^ ζ
    Y_S_c1 = -π_S_c1 # was: Y_S_c1 = c̄_S .- π_S_c1;
    #Case 2: Impossibly to produce enough food:
    x_S_bar_c2 = ((ϕ_S * (1 + Q_S) * ζ) / ((1 + τ_S) * p_x))^(1 / (1 - ζ)) * θ .^ (1 / (1 - ζ))
    x_S_c2 = min.(κ .* s[:, 1] ./ ((1 + τ_S) * p_x), x_S_bar_c2)
    q_S_c2 = ϕ_S * θ .* x_S_c2 .^ ζ
    fertilizer_exp = ((1 + τ_S) * p_x) * x_S_c2
    # λ_S_c2 =   ϕ_S *(1 + Q_S) .* θ .*ζ./p_x .* x_S_c2.^(ζ -1) .- (1 + τ_S); --- KAROL: I THINK THIS IS SLIGHLY WRONG TOO, CORRECTION BELOW (but this is subject to comment on OL):
    λ_S_c2 = ϕ_S * (1 + Q_S) .* θ .* ζ ./ (p_x .* (1 + τ_S)) .* x_S_c2 .^ (ζ - 1) .- 1

    P_S_c2 = ((1 + Q_S)^(1 - ϵ) * ψ_S^ϵ + p_B^(1 - ϵ) * ψ_B^ϵ + p_M^(1 - ϵ) * ψ_M^ϵ)^(1 / (1 - ϵ))
    Y_S_c2 = fertilizer_exp # Might have to pay transaction cost even on parts of c̄_S
    #Case 3: Possible to produce enough food, but distorts the production

    #c_s_solve_mat = zeros(size(agrid_fine)[1],ns); # Stores the consumption solutions at the double grid.
    # Each column is for different θ, row is for different C. Solution must be for each C
    ones_agrid_fine = ones(size(agrid_fine))
    C_grid_fine = copy(agrid_fine)
    exit_flag_mat_fine = zeros(size(agrid_fine)[1], ns)
    TC_S_c3_constrained_fine = zeros(ns)
    x_S_c3_constrained_fine = zeros(ns)
    q_S_c3_constrained_fine = zeros(ns)
    c_S_c3_constrained_fine = zeros(ns)
    for z_index = 1:(convert(Int64, n[2] / no_labor_shock))
        labor_shock_same_index = convert(Int64, ns / no_labor_shock)
        for a_index = 1:n[1]
            # Constrained case: bounds
            TC_S_constrained_tmp = κ * agrid_fine[a_index]
            x_S_constrained_tmp = TC_S_constrained_tmp / ((1 + τ_S) * p_x)
            q_S_constrained_tmp = ϕ_S * x_S_constrained_tmp .^ ζ .* z[z_index]
            condition = (((q_S_constrained_tmp .- c̄_S) ./ ψ_S .^ ϵ ./ C_grid_fine) .^ ((1 - ϵ) / ϵ) .- ψ_S .^ ϵ)
            exitflag_tmp_fine = 0 * copy(ones_agrid_fine)
            exitflag_tmp_fine[condition.<0] .= -1
            for l_index = 0:(no_labor_shock-1)
                exit_flag_mat_fine[:, (labor_shock_same_index*l_index+a_index+((z_index-1)*n[1]))] .= exitflag_tmp_fine
                TC_S_c3_constrained_fine[(labor_shock_same_index*l_index+a_index+((z_index-1)*n[1]))] = TC_S_constrained_tmp
                x_S_c3_constrained_fine[(labor_shock_same_index*l_index+a_index+((z_index-1)*n[1]))] = x_S_constrained_tmp
                q_S_c3_constrained_fine[(labor_shock_same_index*l_index+a_index+((z_index-1)*n[1]))] = q_S_constrained_tmp
                c_S_c3_constrained_fine[(labor_shock_same_index*l_index+a_index+((z_index-1)*n[1]))] = q_S_constrained_tmp
            end
        end
        #coeff_λ_2_s[:,θ_index] = Phi_a_fine\c_s_solve_mat[:,θ_index];
        # for this simple case, this isnt needed as c_s_solve_mat = coeff_λ_2_s!!!
    end
    C_grid_fine_mat = repeat(C_grid_fine, 1, ns)
    C_max_staple_constrained_fine = maximum(C_grid_fine_mat .* (exit_flag_mat_fine .== 0), dims=1)[:]
    C_min_staple_constrained_fine = a_min * ones_tmp
    # Create functions, use staples_objects later.

    # Cash crop farmer
    #π_possible = zeros(ns,3);
    #x_SC_possible  = zeros(ns,3);
    #x_BC_possible = zeros(ns,3);
    #labor_allocated_possible = zeros(ns,3);
    #λ_possible = zeros(ns,3);

    #1) Specialize in cash crop: [24/03/2022 UPDATED LASZLO:]
    x_B_bar_c1 = ((p_B * ϕ_B * ζ) / ((1 + τ_B) * p_x))^(1 / (1 - ζ)) * θ .^ (1 / (1 - ζ))
    x_B_c1 = min.(κ .* s[:, 1] ./ p_x ./ (1 + τ_B), x_B_bar_c1)
    q_B_c1 = ϕ_B * θ .* x_B_c1 .^ ζ
    π_B_only_B_c1 = p_B * q_B_c1 - ((1 + τ_B) * p_x) * x_B_c1
    #λ_B_only_B_c1_tmp =   p_B * ϕ_B .* θ .*ζ./p_x .* x_B_c1.^(ζ -1) ./ (1 + τ_B) .- 1;
    λ_B_only_B_c1 = max.(p_B * ϕ_B .* θ .* ζ .* (κ .* s[:, 1]) .^ (ζ - 1) ./ (p_x * (1 + τ_B)) .^ ζ .- 1, 0.0)
    P_B_c1 = ((1 + Q_S)^(1 - ϵ) * ψ_S^ϵ + p_B^(1 - ϵ) * ψ_B^ϵ + p_M^(1 - ϵ) * ψ_M^ϵ)^(1 / (1 - ϵ))
    Y_B_c1 = -π_B_only_B_c1 #was: Y_B_c1 = (1 + Q_S) * c̄_S .- π_B_only_B_c1;

    #2) Specialize in staple crop:
    # Exactly the same problem as before, so just use function staples_objects

    #π_possible[:,1] = π_B_only_B;
    #x_SC_possible[:,1] .= 0;
    #x_BC_possible[:,1] = x_B;
    #labor_allocated_possible[:,1] .= 1.0;
    #λ_possible[:,1] = λ_B_only_B;

    #3) Internal solution - produces both staples and cash crop
    #3a) q_S>c_S, hence no transaction cost paid [24/03/2022 UPDATED LASZLO:]
    labor_allocated_interior_c3a = (p_B * ϕ_B * (1 + τ_S)^ζ)^(1 / (1 - ρ - ζ)) / ((p_B * ϕ_B * (1 + τ_S)^ζ)^(1 / (1 - ρ - ζ)) + (ϕ_S * (1 + τ_B)^ζ)^(1 / (1 - ρ - ζ)))
    x_SC_interior_c3a = ((ϕ_S * ζ * (1 .- labor_allocated_interior_c3a) .^ ρ) / ((1 + τ_S) * p_x))^(1 / (1 - ζ)) * θ .^ (1 / (1 - ζ))
    x_BC_interior_c3a = ((p_B * ϕ_B * ζ * (labor_allocated_interior_c3a) .^ ρ) / ((1 + τ_B) * p_x))^(1 / (1 - ζ)) * θ .^ (1 / (1 - ζ))
    total_cost_c3a = p_x * ((1 + τ_S) * x_SC_interior_c3a + (1 + τ_B) * x_BC_interior_c3a)
    λ_B_interior_c3a = ζ * θ ./ (κ .* s[:, 1]) .^ (1 - ζ) .* ((p_B * ϕ_B * (1 + τ_S)^ζ)^(1 / (1 - ρ - ζ)) + (ϕ_S * (1 + τ_B)^ζ)^(1 / (1 - ρ - ζ))) .^ (1 - ρ - ζ) ./ (p_x * (1 + τ_S) * (1 + τ_B)) .^ ζ .- 1
    λ_B_interior_c3a = max.(λ_B_interior_c3a, 0.0)
    Lagrange_factor = (1.0 .+ λ_B_interior_c3a) .^ (1 / (ζ - 1))
    #x_tot_interior_c3a = (1+τ_S).*x_SC_interior_c3a .+ (1+τ_B).*x_BC_interior_c3a
    x_SC_interior_c3a = x_SC_interior_c3a .* Lagrange_factor
    x_BC_interior_c3a = x_BC_interior_c3a .* Lagrange_factor
    total_cost_c3a = total_cost_c3a .* Lagrange_factor
    q_S_c3a = ϕ_S * θ .* x_SC_interior_c3a .^ ζ .* (1.0 .- labor_allocated_interior_c3a) .^ ρ
    q_B_c3a = ϕ_B * θ .* x_BC_interior_c3a .^ ζ .* labor_allocated_interior_c3a .^ ρ
    π_B_interior_c3a = q_S_c3a + p_B * q_B_c3a - total_cost_c3a
    P_B_c3a = (ψ_S^ϵ + p_B^(1 - ϵ) * ψ_B^ϵ + p_M^(1 - ϵ) * ψ_M^ϵ)^(1 / (1 - ϵ))
    Y_B_c3a = -π_B_interior_c3a #was: Y_B_c3a = c̄_S.-π_B_interior_c3a


    #3b) q_S<c_S, hence transaction cost is always paid marginally [24/03/2022 UPDATED LASZLO:]
    labor_allocated_interior_c3b = (p_B * ϕ_B * (1 + τ_S)^ζ)^(1 / (1 - ρ - ζ)) / ((p_B * ϕ_B * (1 + τ_S)^ζ)^(1 / (1 - ρ - ζ)) + ((ϕ_S * (1 + Q_S)) * (1 + τ_B)^ζ)^(1 / (1 - ρ - ζ)))
    x_SC_interior_c3b = (((ϕ_S * (1 + Q_S)) * ζ * (1 .- labor_allocated_interior_c3b) .^ ρ) / ((1 + τ_S) * p_x))^(1 / (1 - ζ)) * θ .^ (1 / (1 - ζ))
    x_BC_interior_c3b = ((p_B * ϕ_B * ζ * (labor_allocated_interior_c3b) .^ ρ) / ((1 + τ_B) * p_x))^(1 / (1 - ζ)) * θ .^ (1 / (1 - ζ))
    total_cost_c3b = p_x * ((1 + τ_S) * x_SC_interior_c3b + (1 + τ_B) * x_BC_interior_c3b)
    λ_B_interior_c3b = ζ * θ ./ (κ .* s[:, 1]) .^ (1 - ζ) .* ((p_B * ϕ_B * (1 + τ_S)^ζ)^(1 / (1 - ρ - ζ)) + ((ϕ_S * (1 + Q_S)) * (1 + τ_B)^ζ)^(1 / (1 - ρ - ζ))) .^ (1 - ρ - ζ) ./ (p_x * (1 + τ_S) * (1 + τ_B)) .^ ζ .- 1
    λ_B_interior_c3b = max.(λ_B_interior_c3b, 0.0)
    Lagrange_factor = (1.0 .+ λ_B_interior_c3b) .^ (1 / (ζ - 1))
    x_SC_interior_c3b = x_SC_interior_c3b .* Lagrange_factor
    x_BC_interior_c3b = x_BC_interior_c3b .* Lagrange_factor
    total_cost_c3b = total_cost_c3b .* Lagrange_factor
    q_S_c3b = ϕ_S * θ .* x_SC_interior_c3b .^ ζ .* (1.0 .- labor_allocated_interior_c3b) .^ ρ
    q_B_c3b = ϕ_B * θ .* x_BC_interior_c3b .^ ζ .* labor_allocated_interior_c3b .^ ρ
    π_B_interior_c3b = p_B * q_B_c3b - total_cost_c3b
    Y_B_c3b = -π_B_interior_c3b # market income for the csah crop producer
    P_B_c3b = ((1 + Q_S)^(1 - ϵ) * ψ_S^ϵ + p_B^(1 - ϵ) * ψ_B^ϵ + p_M^(1 - ϵ) * ψ_M^ϵ)^(1 / (1 - ϵ))
    #3c) q_S=c_S, hence transaction cost is not paid, but production is distorted, see functions above
    C_max_mat = zeros(ns, 2)
    C_min_mat = zeros(ns, 2)
    conversion_to_higher_dim = convert(Int64, ns / size(C_max_unconstrained)[1])
    C_max_mat[:, 1] = kron(C_max_unconstrained, ones(conversion_to_higher_dim))
    #C_max_mat[:,2] =kron(C_max_constrained,ones(conversion_to_higher_dim));
    C_min_mat[:, 1] = kron(C_min_unconstrained, ones(conversion_to_higher_dim))
    C_max_staple_fine = kron(C_max_staple, ones(conversion_to_higher_dim))
    #C_max_mat[:,2] =kron(C_max_constrained,ones(conversion_to_higher_dim));
    C_min_staple_fine = kron(C_min_staple, ones(conversion_to_higher_dim))
    #C_min_mat[:,2] =kron(C_min_constrained,ones(conversion_to_higher_dim));
    # Initialization needed for the functions - better do it once here:
    c_S_mat = zeros(ns, 5)
    c_B_mat = zeros(ns, 5)
    c_M_mat = zeros(ns, 5)
    x_S_mat = zeros(ns, 5)
    x_B_mat = zeros(ns, 5)
    q_S_mat = zeros(ns, 5)
    q_B_mat = zeros(ns, 5)
    land_B_mat = zeros(ns, 5)
    λ_2_mat = zeros(ns, 5)
    P_B_mat = zeros(ns, 5)
    feasibility_mat = zeros(ns, 5)
    Y_B_mat = zeros(ns, 5)
    TC_mat = zeros(ns, 5)
    # Initialization needed for the constrained case of 3c functions - better do it once here:

    x_S_mat_3c_fine = zeros(ns, 2)
    x_B_mat_3c_fine = zeros(ns, 2)
    land_B_mat_3c_fine = zeros(ns, 2)
    λ_2_mat_3c_fine = zeros(ns, 2)
    TC_mat_3c_fine = zeros(ns, 2)
    #Calculate the coefficients - the constrained must be rerun as a_grid changes:
    coeff_λ_2_s = kron(coeff_λ_2_s, ones(1, conversion_to_higher_dim))
    coeff_λ_2_cashcrop_residual_unconstrained = kron(coeff_λ_2_cashcrop_residual_unconstrained, ones(1, conversion_to_higher_dim))
    #3c) q_S=c_S, hence transaction cost is not paid, but production is distorted, see functions above
    coeff_λ_2_cashcrop_residual_constrained = zeros(size(agrid_fine)[1], ns)
    V_temp_cashcrop_residual_constrained = zeros(size(agrid_fine)[1], ns)
    for z_index = 1:(convert(Int64, n[2] / no_labor_shock))
        for a_index = 1:n[1]
            (coeff_λ_2_cashcrop_residual_constrained_tmp, V_temp) = goldenx(λ_2_residual_constrained, ones_agrid_fine, (1 + Q_S) * ones_agrid_fine, tol, ϕ_S, ζ, τ_S, p_x, p_B, p_M, ϕ_B, τ_B, c̄_S, Q_S, z, ϵ, ψ_S, ψ_B, ψ_M, C_grid_fine, z_index, agrid_fine, a_index, κ, ρ)
            #V_temp,coeff_x_SC_interior_c3a_constrained_tmp = c_S_cashcrop_residual_constrained(coeff_c_S_cashcrop_residual_constrained_tmp,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,z,ϵ,ψ_S,ψ_B,ψ_M,C_grid_fine,z_index,agrid_fine,a_index,κ,1);
            labor_shock_same_index = convert(Int64, ns / no_labor_shock)
            for l_index = 0:(no_labor_shock-1)
                coeff_λ_2_cashcrop_residual_constrained[:, (labor_shock_same_index*l_index+a_index+((z_index-1)*n[1]))] .= coeff_λ_2_cashcrop_residual_constrained_tmp
                #coeff_x_S_cashcrop_residual_constrained[:,(labor_shock_same_index*l_index+ a_index + ((z_index-1)*n[1]))].=coeff_x_SC_interior_c3a_constrained_tmp;
                V_temp_cashcrop_residual_constrained[:, (labor_shock_same_index*l_index+a_index+((z_index-1)*n[1]))] .= V_temp
            end
        end
    end
    C_max_constrained = maximum(C_grid_fine_mat .* (V_temp_cashcrop_residual_constrained .> -tol), dims=1)[:]
    #C_min_unconstrained = minimum(C_grid_fine_mat .* (V_temp_cashcrop_residual_unconstrained.> -tol)+ 1000*(V_temp_cashcrop_residual_unconstrained.< -tol) ,dims = 1  )[:];
    tmp = zeros(Int64, ns)
    for ii = 1:ns
        if isnothing(findfirst(myCondition, (C_grid_fine_mat.*(V_temp_cashcrop_residual_constrained.>-tol))[:, ii]))
            tmp[ii] = 1
        else
            tmp[ii] = findfirst(myCondition, (C_grid_fine_mat.*(V_temp_cashcrop_residual_constrained.>-tol))[:, ii])
        end
    end
    C_min_constrained = C_grid_fine[tmp]
    #C_min_constrained =  minimum(C_grid_fine_mat .* (V_temp_cashcrop_residual_constrained.> -tol)+ 1000*(V_temp_cashcrop_residual_constrained.< -tol),dims = 1)[:];
    C_max_mat[:, 2] = C_max_constrained
    C_min_mat[:, 2] = C_min_constrained

    #C = s_fine[:,1]/2 .+ a_min;
    #    @btime Y_S,c_S_S,c_B_S,c_M_S,q_S_S,P_S,x_S_S,solve_staple_index_S,λ_2_S = staples_objects(C,ϕ_S,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ζ,c_s_residual,
    #        c_s_constrained,c_s_approx_staples,θ_fine,coeff_λ_2_s_fine,fspace_a_fine,staples_c3_objects,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,s_fine,q_S_c1_fine,
    #        q_S_c2_fine,q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine);
    #C = s_fine[:,1]/2 .+ a_min;
    #Y_S,c_S_S,c_B_S,c_M_S,q_S_S,P_S,x_S_S,solve_staple_index_S,λ_2_S = staples_objects(C,ϕ_S,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ζ,c_s_residual,
    #        c_s_constrained,c_s_approx_staples,θ_fine,coeff_λ_2_s_fine,fspace_a_fine,staples_c3_objects,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine,
    #        x_S_c2_fine,s_fine,q_S_c1_fine,q_S_c2_fine,q_S_staples_fine,
    #        c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine);

    #@btime  (c_S_B,c_B_C,c_M_B,x_SC,x_BC,land_C,λ_2,P_B,Y_B,q_S_C,q_B_B,solve_cash_crop_index_B,solve_staple_index_B) = cashcrop_objects(C,coeff_c_S_cashcrop_residual_unconstrained_fine,
    #        coeff_x_S_cashcrop_residual_unconstrained_fine,coeff_c_S_cashcrop_residual_constrained_fine,θ_fine,
    #         fspace_a_fine,c_S_approx_cashcrop,x_S_approx_cashcrop,s_fine,ϕ_S,ζ,τ_S,p_x,p_B,
    #             p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,
    #             a_min,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
    #             c_s_residual,c_s_constrained,c_s_approx_staples,coeff_λ_2_s_fine,
    #             staples_c3_objects,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,
    #             y_c3a_fine,labor_allocated_interior_c3a_fine,λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,
    #             x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,q_S_c2_fine,
    #             q_B_c1_fine,q_S_c3a_fine,q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,
    #             x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,Y_B_c3b_fine,
    #             c_S_mat_fine,c_B_mat_fine,
    #             c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,q_B_mat_fine,land_B_mat_fine,
    #             λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,feasibility_mat_fine,C_max_mat_fine,C_min_mat_fine,
    #             y_bar_c3c_constraint_mat_fine,
    #             q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,
    #             x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine)

    return (θ, labor_prod, tol, P_W, Y_W, coeff_λ_2_cashcrop_residual_unconstrained, coeff_λ_2_cashcrop_residual_constrained,
        x_B_c1, π_B_only_B_c1, λ_B_only_B_c1, P_B_c1, Y_B_c1,
        coeff_λ_2_s, P_S_c1, P_S_c2, Y_S_c1, Y_S_c2, x_S_c1, x_S_c2, labor_allocated_interior_c3a,
        λ_B_interior_c3a, x_SC_interior_c3a, x_BC_interior_c3a, Y_B_c3a, P_B_c3a, P_B_c3b, q_S_c1, q_S_c2, q_B_c1, q_S_c3a, q_B_c3a, q_S_c3b, q_B_c3b,
        x_SC_interior_c3b, x_BC_interior_c3b, labor_allocated_interior_c3b, Y_B_c3b, c_S_mat, c_B_mat,
        c_M_mat, x_S_mat, x_B_mat, q_S_mat, q_B_mat, land_B_mat, λ_2_mat, P_B_mat, Y_B_mat, feasibility_mat, C_max_mat, C_min_mat, q_S_staples, c_S_staples, c_B_staples,
        c_M_staples, P_S_staples, x_S_staples, λ_2_S_staples, unfeasible_mat, Y_S_potential, TC_mat, C_max_staple_fine, C_min_staple_fine, C_max_staple_constrained_fine,
        C_min_staple_constrained_fine, TC_S_c3_constrained_fine, x_S_c3_constrained_fine, q_S_c3_constrained_fine, c_S_c3_constrained_fine,
        x_S_mat_3c_fine, x_B_mat_3c_fine, land_B_mat_3c_fine, λ_2_mat_3c_fine, TC_mat_3c_fine)
end

function future_asset_creator(C::Array{Float64,1},jj::Int64,j::Int64,P_W::Float64,Y_W::Array{Float64,1},s::Array{Float64,2},r::Float64,ρ::Float64,w::Float64,
    coeff_λ_2_cashcrop_residual_unconstrained::Array{Float64,2},coeff_λ_2_cashcrop_residual_constrained::Array{Float64,2},θ::Array{Float64,1},
    fspace_C_fine::Dict{Symbol,Any},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,p_B::Float64,
    p_M::Float64,ϕ_B::Float64,τ_B::Float64,c̄_S::Float64,Q_S::Float64,ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,ns::Int64,κ::Float64,tol::Float64,
    a_min::Float64,x_B_c1::Array{Float64,1},π_B_only_B_c1::Array{Float64,1},λ_B_only_B_c1::Array{Float64,1},P_B_c1::Float64,Y_B_c1::Array{Float64,1},
    coeff_λ_2_s::Array{Float64,2},P_S_c1::Float64,P_S_c2::Float64,Y_S_c1::Array{Float64,1},Y_S_c2::Array{Float64,1},x_S_c1::Array{Float64,1}, x_S_c2::Array{Float64,1},
    labor_allocated_interior_c3a::Float64,λ_B_interior_c3a::Array{Float64,1},x_SC_interior_c3a::Array{Float64,1},
    x_BC_interior_c3a::Array{Float64,1},Y_B_c3a::Array{Float64,1},P_B_c3a::Float64,P_B_c3b::Float64,q_S_c1::Array{Float64,1},q_S_c2::Array{Float64,1},
    q_B_c1::Array{Float64,1},q_S_c3a::Array{Float64,1},q_B_c3a::Array{Float64,1},q_S_c3b::Array{Float64,1},q_B_c3b::Array{Float64,1},
    x_SC_interior_c3b::Array{Float64,1},x_BC_interior_c3b::Array{Float64,1},labor_allocated_interior_c3b::Float64,Y_B_c3b::Array{Float64,1},
    c_S_mat::Array{Float64,2},c_B_mat::Array{Float64,2},c_M_mat::Array{Float64,2},x_S_mat::Array{Float64,2},x_B_mat::Array{Float64,2},q_S_mat::Array{Float64,2},
    q_B_mat::Array{Float64,2},land_B_mat::Array{Float64,2},λ_2_mat::Array{Float64,2},P_B_mat::Array{Float64,2},Y_B_mat::Array{Float64,2},feasibility_mat::Array{Float64,2},
    C_max_mat::Array{Float64,2},C_min_mat::Array{Float64,2},
    q_S_staples::Array{Float64,2},c_S_staples::Array{Float64,2},c_B_staples::Array{Float64,2},c_M_staples::Array{Float64,2},P_S_staples::Array{Float64,2},
    x_S_staples::Array{Float64,2},λ_2_S_staples::Array{Float64,2},unfeasible_mat::Array{Float64,2},Y_S_potential::Array{Float64,2},F_W::Float64,F_S::Float64,F_B::Float64,
    FM_W::Float64,FM_S::Float64,FM_B::Float64,TC_mat::Array{Float64,2},C_max_staple::Array{Float64,1},
    C_min_staple::Array{Float64,1},C_max_staple_constrained::Array{Float64,1},
    C_min_staple_constrained::Array{Float64,1},TC_S_c3_constrained::Array{Float64,1},
    x_S_c3_constrained::Array{Float64,1},q_S_c3_constrained::Array{Float64,1},c_S_c3_constrained::Array{Float64,1},
    x_S_mat_3c::Array{Float64,2},x_B_mat_3c::Array{Float64,2},land_B_mat_3c::Array{Float64,2},λ_2_mat_3c::Array{Float64,2},TC_mat_3c::Array{Float64,2})
    if jj == 1
        future_asset = (1 + r ) * s[:,1] -P_W *C  -Y_W  .- w*FM_W .- w*F_W * (j != 1);
    elseif jj == 2
        Y_S,c_S_S,c_B_S,c_M_S,q_S_S,P_S,x_S_S,solve_staple_index_S,λ_2_S = staples_objects(C,ϕ_S,τ_S,p_x,
            p_B,p_M,ϕ_B,τ_B,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ζ,θ,coeff_λ_2_s,fspace_C_fine,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,
            x_S_c1, x_S_c2,s,q_S_c1,q_S_c2,q_S_staples,c_S_staples,c_B_staples,c_M_staples,P_S_staples,x_S_staples,
            λ_2_S_staples,unfeasible_mat,Y_S_potential,κ,c̄_S,C_max_staple,
            C_min_staple,ns,C_max_staple_constrained,
            C_min_staple_constrained,TC_S_c3_constrained,
            x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained,a_min);
        future_asset = (1 + r ) * s[:,1] - Y_S[:]  .- w*FM_S .- w*F_S * (j == 1);
    elseif jj == 3
        (c_S_B,c_B_C,c_M_B,x_SC,x_BC,land_C,λ_2,P_B,Y_B,q_S_C,q_B_B,solve_cash_crop_index_B,solve_staple_index_B,TC_B) = cashcrop_objects(C,
        coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,
            θ,fspace_C_fine,s,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns,κ,tol,a_min,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
                coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,λ_B_interior_c3a,x_SC_interior_c3a,
                x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,
                x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b,c_S_mat,c_B_mat,c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat,
                λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,TC_mat,q_S_staples,c_S_staples,c_B_staples,c_M_staples,P_S_staples,
                x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,ρ,C_max_staple,
                C_min_staple,C_max_staple_constrained,C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained,
                x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c);
        future_asset = (1 + r ) * s[:,1] - Y_B[:]  .- w*FM_B .- w*F_B * (j != 3);
    end
    return future_asset
end
function future_asset_creator_inv(C::Array{Float64,1},jj::Int64,j::Int64,P_W::Float64,Y_W::Array{Float64,1},s::Array{Float64,2},r::Float64,ρ::Float64,w::Float64,
    coeff_λ_2_cashcrop_residual_unconstrained::Array{Float64,2},coeff_λ_2_cashcrop_residual_constrained::Array{Float64,2},θ::Array{Float64,1},
    fspace_C_fine::Dict{Symbol,Any},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,p_B::Float64,
    p_M::Float64,ϕ_B::Float64,τ_B::Float64,c̄_S::Float64,Q_S::Float64,ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,ns::Int64,κ::Float64,tol::Float64,
    a_min::Float64,x_B_c1::Array{Float64,1},π_B_only_B_c1::Array{Float64,1},λ_B_only_B_c1::Array{Float64,1},P_B_c1::Float64,Y_B_c1::Array{Float64,1},
    coeff_λ_2_s::Array{Float64,2},P_S_c1::Float64,P_S_c2::Float64,Y_S_c1::Array{Float64,1},Y_S_c2::Array{Float64,1},x_S_c1::Array{Float64,1}, x_S_c2::Array{Float64,1},
    labor_allocated_interior_c3a::Float64,λ_B_interior_c3a::Array{Float64,1},x_SC_interior_c3a::Array{Float64,1},
    x_BC_interior_c3a::Array{Float64,1},Y_B_c3a::Array{Float64,1},P_B_c3a::Float64,P_B_c3b::Float64,q_S_c1::Array{Float64,1},q_S_c2::Array{Float64,1},
    q_B_c1::Array{Float64,1},q_S_c3a::Array{Float64,1},q_B_c3a::Array{Float64,1},q_S_c3b::Array{Float64,1},q_B_c3b::Array{Float64,1},
    x_SC_interior_c3b::Array{Float64,1},x_BC_interior_c3b::Array{Float64,1},labor_allocated_interior_c3b::Float64,Y_B_c3b::Array{Float64,1},
    c_S_mat::Array{Float64,2},c_B_mat::Array{Float64,2},c_M_mat::Array{Float64,2},x_S_mat::Array{Float64,2},x_B_mat::Array{Float64,2},q_S_mat::Array{Float64,2},
    q_B_mat::Array{Float64,2},land_B_mat::Array{Float64,2},λ_2_mat::Array{Float64,2},P_B_mat::Array{Float64,2},Y_B_mat::Array{Float64,2},feasibility_mat::Array{Float64,2},
    C_max_mat::Array{Float64,2},C_min_mat::Array{Float64,2},
    q_S_staples::Array{Float64,2},c_S_staples::Array{Float64,2},c_B_staples::Array{Float64,2},c_M_staples::Array{Float64,2},P_S_staples::Array{Float64,2},
    x_S_staples::Array{Float64,2},λ_2_S_staples::Array{Float64,2},unfeasible_mat::Array{Float64,2},Y_S_potential::Array{Float64,2},F_W::Float64,F_S::Float64,F_B::Float64,
    FM_W::Float64,FM_S::Float64,FM_B::Float64,TC_mat::Array{Float64,2},C_max_staple::Array{Float64,1},
    C_min_staple::Array{Float64,1},C_max_staple_constrained::Array{Float64,1},
    C_min_staple_constrained::Array{Float64,1},TC_S_c3_constrained::Array{Float64,1},
    x_S_c3_constrained::Array{Float64,1},q_S_c3_constrained::Array{Float64,1},c_S_c3_constrained::Array{Float64,1},
    x_S_mat_3c::Array{Float64,2},x_B_mat_3c::Array{Float64,2},land_B_mat_3c::Array{Float64,2},λ_2_mat_3c::Array{Float64,2},TC_mat_3c::Array{Float64,2})
    return -(future_asset_creator(C,jj,j,P_W,Y_W,s,r,ρ,w,
        coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,
        θ,fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns,κ,tol,a_min,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
                coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,
                λ_B_interior_c3a,x_SC_interior_c3a,x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,
                x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b,c_S_mat,c_B_mat,
                c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat, λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,q_S_staples,
                c_S_staples,c_B_staples,c_M_staples,P_S_staples,x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,F_W,F_S,F_B,
                FM_W,FM_S,FM_B,TC_mat,C_max_staple,C_min_staple,C_max_staple_constrained,C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,
                c_S_c3_constrained, x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c) .- a_min).^2;
end


function bounds_consumption(P_W::Float64,Y_W::Array{Float64,1},s::Array{Float64,2},r::Float64,ρ::Float64,w::Float64,
    coeff_λ_2_cashcrop_residual_unconstrained::Array{Float64,2},coeff_λ_2_cashcrop_residual_constrained::Array{Float64,2},θ::Array{Float64,1},
    fspace_C_fine::Dict{Symbol,Any},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,p_B::Float64,
    p_M::Float64,ϕ_B::Float64,τ_B::Float64,c̄_S::Float64,Q_S::Float64,ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,ns::Int64,κ::Float64,tol::Float64,
    a_min::Float64,a_max::Float64,x_B_c1::Array{Float64,1},π_B_only_B_c1::Array{Float64,1},λ_B_only_B_c1::Array{Float64,1},P_B_c1::Float64,Y_B_c1::Array{Float64,1},
    coeff_λ_2_s::Array{Float64,2},P_S_c1::Float64,P_S_c2::Float64,Y_S_c1::Array{Float64,1},Y_S_c2::Array{Float64,1},x_S_c1::Array{Float64,1},
    x_S_c2::Array{Float64,1},labor_allocated_interior_c3a::Float64,λ_B_interior_c3a::Array{Float64,1},x_SC_interior_c3a::Array{Float64,1},
    x_BC_interior_c3a::Array{Float64,1},Y_B_c3a::Array{Float64,1},P_B_c3a::Float64,P_B_c3b::Float64,q_S_c1::Array{Float64,1},q_S_c2::Array{Float64,1},
    q_B_c1::Array{Float64,1},q_S_c3a::Array{Float64,1},q_B_c3a::Array{Float64,1},q_S_c3b::Array{Float64,1},q_B_c3b::Array{Float64,1},
    x_SC_interior_c3b::Array{Float64,1},x_BC_interior_c3b::Array{Float64,1},labor_allocated_interior_c3b::Float64,Y_B_c3b::Array{Float64,1},
    c_S_mat::Array{Float64,2},c_B_mat::Array{Float64,2},c_M_mat::Array{Float64,2},x_S_mat::Array{Float64,2},x_B_mat::Array{Float64,2},q_S_mat::Array{Float64,2},
    q_B_mat::Array{Float64,2},land_B_mat::Array{Float64,2},λ_2_mat::Array{Float64,2},P_B_mat::Array{Float64,2},Y_B_mat::Array{Float64,2},feasibility_mat::Array{Float64,2},
    C_max_mat::Array{Float64,2},C_min_mat::Array{Float64,2},q_S_staples::Array{Float64,2},c_S_staples::Array{Float64,2},c_B_staples::Array{Float64,2},
    c_M_staples::Array{Float64,2},P_S_staples::Array{Float64,2},x_S_staples::Array{Float64,2},λ_2_S_staples::Array{Float64,2},unfeasible_mat::Array{Float64,2},
    Y_S_potential::Array{Float64,2},F_W::Float64,F_S::Float64,F_B::Float64,FM_W::Float64,FM_S::Float64,FM_B::Float64,TC_mat::Array{Float64,2},C_max_staple::Array{Float64,1},
    C_min_staple::Array{Float64,1},C_max_staple_constrained::Array{Float64,1},
    C_min_staple_constrained::Array{Float64,1},TC_S_c3_constrained::Array{Float64,1},
    x_S_c3_constrained::Array{Float64,1},q_S_c3_constrained::Array{Float64,1},c_S_c3_constrained::Array{Float64,1},
    x_S_mat_3c::Array{Float64,2},x_B_mat_3c::Array{Float64,2},land_B_mat_3c::Array{Float64,2},λ_2_mat_3c::Array{Float64,2},TC_mat_3c::Array{Float64,2})
    C_min_abs = a_min * ones(ns);
    C_max_abs = a_max * ones(ns);
    min_C_applied = zeros(ns,3,3);
    max_C_applied = zeros(ns,3,3);
    for j = 1:3
        for jj =1:3
            min_C_applied[:,jj,j] .= a_min
            max_C_applied[:,jj,j],V_tmp = goldenx(future_asset_creator_inv,C_min_abs,C_max_abs ,tol,jj,j,P_W,Y_W,s,r,ρ,w,
                coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,
                θ,fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns,κ,tol,a_min,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
                        coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,
                        λ_B_interior_c3a,x_SC_interior_c3a,x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,
                        x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b,c_S_mat,c_B_mat,
                        c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat, λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,q_S_staples,
                        c_S_staples,c_B_staples,c_M_staples,P_S_staples,x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,F_W,F_S,F_B,
                        FM_W,FM_S,FM_B,TC_mat,C_max_staple,C_min_staple,C_max_staple_constrained,C_min_staple_constrained,TC_S_c3_constrained,
                        x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained, x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c);
            #println(jj,' ',j,' ',minimum(V_tmp))
        end
    end
    return min_C_applied,max_C_applied
end

function BellmanF_stable_transaction(cons_guess::Array{Float64,1},coeff::Array{Float64,2},Phi_z_co::SparseMatrixCSC{Float64,Int64},
        fspace_a_co::Dict{Symbol,Any},β_applied::Float64,σ::Float64,jj::Int64,j::Int64,P_W::Float64,Y_W::Array{Float64,1},s::Array{Float64,2},r::Float64,ρ::Float64,w::Float64,
        coeff_λ_2_cashcrop_residual_unconstrained::Array{Float64,2},coeff_λ_2_cashcrop_residual_constrained::Array{Float64,2},θ::Array{Float64,1},
        fspace_C_fine::Dict{Symbol,Any},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,p_B::Float64,
        p_M::Float64,ϕ_B::Float64,τ_B::Float64,c̄_S::Float64,Q_S::Float64,ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,ns::Int64,κ::Float64,tol::Float64,
        a_min::Float64,x_B_c1::Array{Float64,1},π_B_only_B_c1::Array{Float64,1},λ_B_only_B_c1::Array{Float64,1},P_B_c1::Float64,Y_B_c1::Array{Float64,1},
        coeff_λ_2_s::Array{Float64,2},P_S_c1::Float64,P_S_c2::Float64,Y_S_c1::Array{Float64,1},Y_S_c2::Array{Float64,1},x_S_c1::Array{Float64,1}, x_S_c2::Array{Float64,1},
        labor_allocated_interior_c3a::Float64,λ_B_interior_c3a::Array{Float64,1},x_SC_interior_c3a::Array{Float64,1},
        x_BC_interior_c3a::Array{Float64,1},Y_B_c3a::Array{Float64,1},P_B_c3a::Float64,P_B_c3b::Float64,q_S_c1::Array{Float64,1},q_S_c2::Array{Float64,1},
        q_B_c1::Array{Float64,1},q_S_c3a::Array{Float64,1},q_B_c3a::Array{Float64,1},q_S_c3b::Array{Float64,1},q_B_c3b::Array{Float64,1},
        x_SC_interior_c3b::Array{Float64,1},x_BC_interior_c3b::Array{Float64,1},labor_allocated_interior_c3b::Float64,Y_B_c3b::Array{Float64,1},
        c_S_mat::Array{Float64,2},c_B_mat::Array{Float64,2},c_M_mat::Array{Float64,2},x_S_mat::Array{Float64,2},x_B_mat::Array{Float64,2},q_S_mat::Array{Float64,2},
        q_B_mat::Array{Float64,2},land_B_mat::Array{Float64,2},λ_2_mat::Array{Float64,2},P_B_mat::Array{Float64,2},Y_B_mat::Array{Float64,2},feasibility_mat::Array{Float64,2},
        C_max_mat::Array{Float64,2},C_min_mat::Array{Float64,2},
        q_S_staples::Array{Float64,2},c_S_staples::Array{Float64,2},c_B_staples::Array{Float64,2},c_M_staples::Array{Float64,2},P_S_staples::Array{Float64,2},
        x_S_staples::Array{Float64,2},λ_2_S_staples::Array{Float64,2},unfeasible_mat::Array{Float64,2},Y_S_potential::Array{Float64,2},F_W::Float64,F_S::Float64,F_B::Float64,
        FM_W::Float64,FM_S::Float64,FM_B::Float64,TC_mat::Array{Float64,2},a_max::Float64,C_max_staple::Array{Float64,1},
        C_min_staple::Array{Float64,1},C_max_staple_constrained::Array{Float64,1},
        C_min_staple_constrained::Array{Float64,1},TC_S_c3_constrained::Array{Float64,1},
        x_S_c3_constrained::Array{Float64,1},q_S_c3_constrained::Array{Float64,1},c_S_c3_constrained::Array{Float64,1}, 
        x_S_mat_3c::Array{Float64,2},x_B_mat_3c::Array{Float64,2},land_B_mat_3c::Array{Float64,2},λ_2_mat_3c::Array{Float64,2},TC_mat_3c::Array{Float64,2},
        Return_Phi::Int64=0)
    #Current utility conditional on occupation
    x = future_asset_creator(cons_guess,jj,j,P_W,Y_W,s,r,ρ,w,
        coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,
        θ,fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns,κ,tol,a_min,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
                coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,
                λ_B_interior_c3a,x_SC_interior_c3a,x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,
                x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b,c_S_mat,c_B_mat,
                c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat, λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,q_S_staples,
                c_S_staples,c_B_staples,c_M_staples,P_S_staples,x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,F_W,F_S,F_B,
                FM_W,FM_S,FM_B,TC_mat,C_max_staple,C_min_staple,C_max_staple_constrained,C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained,
                x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c);
    #println(sum(x.<a_min))
    below_index = x.<a_min;
    above_index = x.>a_max;
    future_asset = below_index.*a_min + (1 .-below_index).*(1 .-above_index).*x + above_index.*a_max;#x;#
    util_penalty_below = below_index.*((x.- a_min).^2 .+ 1000)
    util_penalty_above = above_index.*((a_max .- x).^2 .+ 1000);
    #println(sum(x.>a_max))
    #future_asset = min.(future_asset,a_max);
    #future_asset = copy(x);
    Phi_prime_a = funbase(fspace_a_co, future_asset);
    Phi_prime = row_kron(Phi_z_co,Phi_prime_a); # Note - check if RowKron==dprod
    #Phi_prime = row_kron(Phi_prime_a,Phi_z_co);
    V_fut_mat = Phi_prime * coeff[:,jj];
    #util= curr_util(cons_guess,σ) + β_applied*V_fut_mat - ((x.<a_min) .*(a_min .- x ).^2 * 100000 + (x.>a_max).*(x .- a_max).^2 * 100000);
    util_occupation= curr_util(cons_guess,σ) + β_applied*V_fut_mat - util_penalty_below*1000; # + (x.>a_max)*1000;
    util = util_occupation - util_penalty_above*1000;
    if Return_Phi==0
        return util
    elseif Return_Phi==1
        return util_occupation , Phi_prime, future_asset
    end
end

function Bellman_iteration(coeff::Array{Float64,2},coeff_next::Array{Float64,2},x_tmp::Array{Float64,2},V_tmp::Array{Float64,2},
    P_kron::SparseMatrixCSC{Float64,Int64},Phi::SparseMatrixCSC{Float64,Int64},
    min_C_applied::Array{Float64,3},max_C_applied::Array{Float64,3},Phi_z::SparseMatrixCSC{Float64,Int64},β_loc::Float64,
    fspace_a::Dict{Symbol,Any},σ::Float64,P_W::Float64,Y_W::Array{Float64,1},s::Array{Float64,2},r::Float64,ρ::Float64,w::Float64,
    coeff_λ_2_cashcrop_residual_unconstrained::Array{Float64,2},coeff_λ_2_cashcrop_residual_constrained::Array{Float64,2},θ::Array{Float64,1},
    fspace_C_fine::Dict{Symbol,Any},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,p_B::Float64,
    p_M::Float64,ϕ_B::Float64,τ_B::Float64,c̄_S::Float64,Q_S::Float64,ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,ns::Int64,κ::Float64,tol::Float64,
    a_min::Float64,x_B_c1::Array{Float64,1},π_B_only_B_c1::Array{Float64,1},λ_B_only_B_c1::Array{Float64,1},P_B_c1::Float64,Y_B_c1::Array{Float64,1},
    coeff_λ_2_s::Array{Float64,2},P_S_c1::Float64,P_S_c2::Float64,Y_S_c1::Array{Float64,1},Y_S_c2::Array{Float64,1},x_S_c1::Array{Float64,1}, x_S_c2::Array{Float64,1},
    labor_allocated_interior_c3a::Float64,λ_B_interior_c3a::Array{Float64,1},x_SC_interior_c3a::Array{Float64,1},
    x_BC_interior_c3a::Array{Float64,1},Y_B_c3a::Array{Float64,1},P_B_c3a::Float64,P_B_c3b::Float64,q_S_c1::Array{Float64,1},q_S_c2::Array{Float64,1},
    q_B_c1::Array{Float64,1},q_S_c3a::Array{Float64,1},q_B_c3a::Array{Float64,1},q_S_c3b::Array{Float64,1},q_B_c3b::Array{Float64,1},
    x_SC_interior_c3b::Array{Float64,1},x_BC_interior_c3b::Array{Float64,1},labor_allocated_interior_c3b::Float64,Y_B_c3b::Array{Float64,1},
    c_S_mat::Array{Float64,2},c_B_mat::Array{Float64,2},c_M_mat::Array{Float64,2},x_S_mat::Array{Float64,2},x_B_mat::Array{Float64,2},q_S_mat::Array{Float64,2},
    q_B_mat::Array{Float64,2},land_B_mat::Array{Float64,2},λ_2_mat::Array{Float64,2},P_B_mat::Array{Float64,2},Y_B_mat::Array{Float64,2},feasibility_mat::Array{Float64,2},
    C_max_mat::Array{Float64,2},C_min_mat::Array{Float64,2},
    q_S_staples::Array{Float64,2},c_S_staples::Array{Float64,2},c_B_staples::Array{Float64,2},c_M_staples::Array{Float64,2},P_S_staples::Array{Float64,2},
    x_S_staples::Array{Float64,2},λ_2_S_staples::Array{Float64,2},unfeasible_mat::Array{Float64,2},Y_S_potential::Array{Float64,2},F_W::Float64,F_S::Float64,F_B::Float64,
    FM_W::Float64,FM_S::Float64,FM_B::Float64,TC_mat::Array{Float64,2},a_max::Float64,C_max_staple::Array{Float64,1},
    C_min_staple::Array{Float64,1},C_max_staple_constrained::Array{Float64,1},
    C_min_staple_constrained::Array{Float64,1},TC_S_c3_constrained::Array{Float64,1},
    x_S_c3_constrained::Array{Float64,1},q_S_c3_constrained::Array{Float64,1},c_S_c3_constrained::Array{Float64,1}, 
    x_S_mat_3c::Array{Float64,2},x_B_mat_3c::Array{Float64,2},land_B_mat_3c::Array{Float64,2},λ_2_mat_3c::Array{Float64,2},TC_mat_3c::Array{Float64,2})
    for j = 1:3
        for jj =1:3
            max_x_tmp = max_C_applied[:,jj,j];
            min_x_tmp = min_C_applied[:,jj,j];
            x_tmp[:,jj],V_tmp[:,jj] = goldenx(BellmanF_stable_transaction,min_x_tmp,
                max_x_tmp,tol,coeff,Phi_z,fspace_a,β_loc,σ,jj,j,P_W,Y_W,s,r,ρ,w,
                    coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,θ,
                    fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,
                    p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns,κ,tol,a_min,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
                    coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,
                    λ_B_interior_c3a,x_SC_interior_c3a,x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,
                    x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b, c_S_mat,c_B_mat,
                    c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat, λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,q_S_staples,c_S_staples,c_B_staples,
                    c_M_staples,P_S_staples,x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat,a_max,C_max_staple,C_min_staple,
                    C_max_staple_constrained,C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained,
                    x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c);
            (V_tmp[:,jj], Phi_prime_tmp_loc, x_tmp[:,jj])= BellmanF_stable_transaction(x_tmp[:,jj],coeff,Phi_z,
                    fspace_a,β_loc,σ,jj,j,P_W,Y_W,s,r,ρ,w,
                        coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,θ,
                        fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,
                        p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns,κ,tol,a_min,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
                        coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,
                        λ_B_interior_c3a,x_SC_interior_c3a,x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,
                        x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b, c_S_mat,c_B_mat,
                        c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat, λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,q_S_staples,c_S_staples,c_B_staples,
                        c_M_staples,P_S_staples,x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat,a_max,C_max_staple,
                        C_min_staple,C_max_staple_constrained,C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained,
                        x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c,1);
        end
        V_tmp[isnan.(V_tmp)] .= -10000;
        #V_tmp[isinf.(V_tmp)] .= -10000;
        #V_tmp[:,j] = V_tmp[:,j] .+ 0.01;
        V_next, future_occupation = findmax(V_tmp; dims=2);
        E_V_next = P_kron * V_next;
        coeff_next[:,j] = Phi\E_V_next;#+ coeff[:,j];
        
        #println(future_occupation_index)
        # future_occupation_index = getindex.(future_occupation,2);
        #  if j==1
        #      println("choosing urban from urban: ",sum(future_occupation_index.==1))
        #      println("choosing staplefarmer from urban: ",sum(future_occupation_index.==2))
        #      println("choosing cashcropfarmer from urban: ",sum(future_occupation_index.==3))
        #  elseif j==2
        #      println("choosing urban from staplefarmer: ",sum(future_occupation_index.==1))
        #      println("choosing staplefarmer from staplefarmer: ",sum(future_occupation_index.==2))
        #      println("choosing cashcropfarmer from staplefarmer: ",sum(future_occupation_index.==3))
        #  elseif j==3
        #      println("choosing urban from cashcropfarmer: ",sum(future_occupation_index.==1))
        #      println("choosing staplefarmer from cashcropfarmer: ",sum(future_occupation_index.==2))
        #      println("choosing cashcropfarmer from cashcropfarmer: ",sum(future_occupation_index.==3))
        #  end 
    end
    conv = maximum(abs.(coeff_next - coeff));
    conv_ind = findmax(abs.(coeff_next - coeff));
    #conv = mean(abs.(coeff_next - coeff));

    #println("Bellman iteration convergence: ", conv)
    return (coeff_next,conv,conv_ind)
end

function Newton_iteration(coeff::Array{Float64,2},x_tmp::Array{Float64,2},V_tmp::Array{Float64,2},
    P_kron::SparseMatrixCSC{Float64,Int64},Phi::SparseMatrixCSC{Float64,Int64},
    min_C_applied::Array{Float64,3},max_C_applied::Array{Float64,3},Phi_z::SparseMatrixCSC{Float64,Int64},β_loc::Float64,
    fspace_a::Dict{Symbol,Any},σ::Float64,P_W::Float64,Y_W::Array{Float64,1},s::Array{Float64,2},r::Float64,ρ::Float64,w::Float64,
    coeff_λ_2_cashcrop_residual_unconstrained::Array{Float64,2},coeff_λ_2_cashcrop_residual_constrained::Array{Float64,2},θ::Array{Float64,1},
    fspace_C_fine::Dict{Symbol,Any},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,p_B::Float64,
    p_M::Float64,ϕ_B::Float64,τ_B::Float64,c̄_S::Float64,Q_S::Float64,ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,ns::Int64,κ::Float64,tol::Float64,
    a_min::Float64,x_B_c1::Array{Float64,1},π_B_only_B_c1::Array{Float64,1},λ_B_only_B_c1::Array{Float64,1},P_B_c1::Float64,Y_B_c1::Array{Float64,1},
    coeff_λ_2_s::Array{Float64,2},P_S_c1::Float64,P_S_c2::Float64,Y_S_c1::Array{Float64,1},Y_S_c2::Array{Float64,1},x_S_c1::Array{Float64,1}, x_S_c2::Array{Float64,1},
    labor_allocated_interior_c3a::Float64,λ_B_interior_c3a::Array{Float64,1},x_SC_interior_c3a::Array{Float64,1},
    x_BC_interior_c3a::Array{Float64,1},Y_B_c3a::Array{Float64,1},P_B_c3a::Float64,P_B_c3b::Float64,q_S_c1::Array{Float64,1},q_S_c2::Array{Float64,1},
    q_B_c1::Array{Float64,1},q_S_c3a::Array{Float64,1},q_B_c3a::Array{Float64,1},q_S_c3b::Array{Float64,1},q_B_c3b::Array{Float64,1},
    x_SC_interior_c3b::Array{Float64,1},x_BC_interior_c3b::Array{Float64,1},labor_allocated_interior_c3b::Float64,Y_B_c3b::Array{Float64,1},
    c_S_mat::Array{Float64,2},c_B_mat::Array{Float64,2},c_M_mat::Array{Float64,2},x_S_mat::Array{Float64,2},x_B_mat::Array{Float64,2},q_S_mat::Array{Float64,2},
    q_B_mat::Array{Float64,2},land_B_mat::Array{Float64,2},λ_2_mat::Array{Float64,2},P_B_mat::Array{Float64,2},Y_B_mat::Array{Float64,2},feasibility_mat::Array{Float64,2},
    C_max_mat::Array{Float64,2},C_min_mat::Array{Float64,2},
    q_S_staples::Array{Float64,2},c_S_staples::Array{Float64,2},c_B_staples::Array{Float64,2},c_M_staples::Array{Float64,2},P_S_staples::Array{Float64,2},
    x_S_staples::Array{Float64,2},λ_2_S_staples::Array{Float64,2},unfeasible_mat::Array{Float64,2},Y_S_potential::Array{Float64,2},F_W::Float64,F_S::Float64,F_B::Float64,
    FM_W::Float64,FM_S::Float64,FM_B::Float64,TC_mat::Array{Float64,2},a_max::Float64,
    V_next_stacked::Array{Float64,2},iterate1::Int64,
    Phi_prime_tmp::Array{SparseMatrixCSC{Float64,Int64}},
    D_deriv_tmp_block::Array{SparseMatrixCSC{Float64,Int64}},
    Phi_aug::SparseMatrixCSC{Float64,Int64},
    P_kron1::SparseMatrixCSC{Float64,Int64},exitflag_tmp::Int64,C_max_staple::Array{Float64,1},
    C_min_staple::Array{Float64,1},C_max_staple_constrained::Array{Float64,1},
    C_min_staple_constrained::Array{Float64,1},TC_S_c3_constrained::Array{Float64,1},
    x_S_c3_constrained::Array{Float64,1},q_S_c3_constrained::Array{Float64,1},c_S_c3_constrained::Array{Float64,1},
    x_S_mat_3c::Array{Float64,2},x_B_mat_3c::Array{Float64,2},land_B_mat_3c::Array{Float64,2},λ_2_mat_3c::Array{Float64,2},TC_mat_3c::Array{Float64,2})
    D_deriv = spzeros(ns*3,ns*3);
    for j = 1:3
        #future_asset_tmp = copy(x_tmp)*0
        for jj =1:3
            max_x_tmp = max_C_applied[:,jj,j];
            min_x_tmp = min_C_applied[:,jj,j];
            
            x_tmp[:,jj],fval = goldenx(BellmanF_stable_transaction,min_x_tmp,
                max_x_tmp,tol,coeff,Phi_z,fspace_a,β_loc,σ,jj,j,P_W,Y_W,s,r,ρ,w,
                    coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,θ,
                    fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,
                    p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns,κ,tol,a_min,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
                    coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,
                    λ_B_interior_c3a,x_SC_interior_c3a,x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,
                    x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b, c_S_mat,c_B_mat,
                    c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat, λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,q_S_staples,c_S_staples,c_B_staples,
                    c_M_staples,P_S_staples,x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat,a_max,C_max_staple,C_min_staple,
                    C_max_staple_constrained,C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained,
                    x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c);
            #(V_tmp[:,jj], Phi_prime_tmp_loc, future_asset_tmp[:,jj])= BellmanF_stable_transaction(x_tmp[:,jj],coeff,Phi_z,
            (V_tmp[:,jj], Phi_prime_tmp_loc, future_asset)= BellmanF_stable_transaction(x_tmp[:,jj],coeff,Phi_z,
                    fspace_a,β_loc,σ,jj,j,P_W,Y_W,s,r,ρ,w,
                        coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,θ,
                        fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,
                        p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns,κ,tol,a_min,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
                        coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,
                        λ_B_interior_c3a,x_SC_interior_c3a,x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,
                        x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b, c_S_mat,c_B_mat,
                        c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat, λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,q_S_staples,c_S_staples,c_B_staples,
                        c_M_staples,P_S_staples,x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,F_W,F_S,F_B,FM_W,FM_S,FM_B,TC_mat,a_max,C_max_staple,
                        C_min_staple,C_max_staple_constrained,C_min_staple_constrained,TC_S_c3_constrained,x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained,
                        x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c,1);
            check_nonzero = Phi_prime_tmp_loc.>tol;
            Phi_prime_tmp[jj] = Phi_prime_tmp_loc.*check_nonzero;
        end
        V_tmp[isnan.(V_tmp)] .= -10000;
        #V_tmp[isinf.(V_tmp)] .= -10000;
        V_next, future_occupation = findmax(V_tmp; dims=2);
        future_occupation_index = getindex.(future_occupation,2);
        for jj=1:3
            D_deriv_tmp_block[jj,1] = P_kron * row_kron(sparse(hcat(float(
                        future_occupation_index[:,1].==jj))),Phi_prime_tmp[jj,1]);
        end
        D_deriv[((j-1)*ns + 1):(j*ns),1:ns] = D_deriv_tmp_block[1,1];
        D_deriv[((j-1)*ns + 1):(j*ns),(ns+1):(2*ns)]= D_deriv_tmp_block[2,1];
        D_deriv[((j-1)*ns + 1):(j*ns),(2*ns+1):(3*ns)] = D_deriv_tmp_block[3,1];
        V_next_stacked[((j-1)*ns + 1):(j*ns),1] = V_next;
    end
    D = Phi_aug - β_loc * D_deriv;
    E_V_next = P_kron1 * V_next_stacked;
    g = Phi_aug *[coeff[:,1];coeff[:,2];coeff[:,3]] - E_V_next;
    improvement = D\g;
    improvement_mat = zeros(size(coeff));
    improvement_mat[:,1] = improvement[((1-1)*ns + 1):(1*ns)];
    improvement_mat[:,2] = improvement[((2-1)*ns + 1):(2*ns)];
    improvement_mat[:,3] = improvement[((3-1)*ns + 1):(3*ns)];
    coeff_next = coeff - improvement_mat;
    conv = maximum(abs.(coeff_next - coeff));
    conv_ind = findmax(abs.(coeff_next - coeff));
    #conv = mean(abs.(coeff_next - coeff));
    #println("Newton iteration convergence: ", conv)
    iterate1 = iterate1 + 1;
    if iterate1 ==20
        #println("Converge not achieved after ", iterate1, "Newton iterations")
        if conv>10^(-4)
            # Only stop if the imprecision is very high
            exitflag_tmp = 4;
        end
        #conv =0;
    end
    return (coeff_next,conv,iterate1,exitflag_tmp,conv_ind)
end

function Q_transition(coeff::Array{Float64,2},
    ns_fine::Int64,P_kron_fine::SparseMatrixCSC{Float64,Int64},
    min_C_applied_fine::Array{Float64,3},max_C_applied_fine::Array{Float64,3},Phi_z_fine::SparseMatrixCSC{Float64,Int64},β_loc::Float64,
    fspace_a::Dict{Symbol,Any},σ::Float64,P_W_fine::Float64,Y_W_fine::Array{Float64,1},s_fine::Array{Float64,2},r::Float64,ρ::Float64,w::Float64,
    coeff_λ_2_cashcrop_residual_unconstrained_fine::Array{Float64,2},coeff_λ_2_cashcrop_residual_constrained_fine::Array{Float64,2},θ_fine::Array{Float64,1},
    fspace_C_fine::Dict{Symbol,Any},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,p_B::Float64,p_M::Float64,ϕ_B::Float64,τ_B::Float64,c̄_S::Float64,Q_S::Float64,
    ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,κ::Float64,tol::Float64,a_min::Float64,x_B_c1_fine::Array{Float64,1},π_B_only_B_c1_fine::Array{Float64,1},
    λ_B_only_B_c1_fine::Array{Float64,1},P_B_c1_fine::Float64,Y_B_c1_fine::Array{Float64,1},coeff_λ_2_s_fine::Array{Float64,2},P_S_c1_fine::Float64,P_S_c2_fine::Float64,
    Y_S_c1_fine::Array{Float64,1},Y_S_c2_fine::Array{Float64,1},x_S_c1_fine::Array{Float64,1},x_S_c2_fine::Array{Float64,1},labor_allocated_interior_c3a_fine::Float64,
    λ_B_interior_c3a_fine::Array{Float64,1},x_SC_interior_c3a_fine::Array{Float64,1},x_BC_interior_c3a_fine::Array{Float64,1},Y_B_c3a_fine::Array{Float64,1},
    P_B_c3a_fine::Float64,P_B_c3b_fine::Float64,q_S_c1_fine::Array{Float64,1},q_S_c2_fine::Array{Float64,1},q_B_c1_fine::Array{Float64,1},q_S_c3a_fine::Array{Float64,1},
    q_B_c3a_fine::Array{Float64,1},q_S_c3b_fine::Array{Float64,1},q_B_c3b_fine::Array{Float64,1},x_SC_interior_c3b_fine::Array{Float64,1},x_BC_interior_c3b_fine::Array{Float64,1},
    labor_allocated_interior_c3b_fine::Float64,Y_B_c3b_fine::Array{Float64,1},c_S_mat_fine::Array{Float64,2},c_B_mat_fine::Array{Float64,2},
    c_M_mat_fine::Array{Float64,2},x_S_mat_fine::Array{Float64,2},x_B_mat_fine::Array{Float64,2},q_S_mat_fine::Array{Float64,2},q_B_mat_fine::Array{Float64,2},
    land_B_mat_fine::Array{Float64,2},λ_2_mat_fine::Array{Float64,2},P_B_mat_fine::Array{Float64,2},Y_B_mat_fine::Array{Float64,2},feasibility_mat_fine::Array{Float64,2},
    C_max_mat_fine::Array{Float64,2},C_min_mat_fine::Array{Float64,2},q_S_staples_fine::Array{Float64,2},c_S_staples_fine::Array{Float64,2},
    c_B_staples_fine::Array{Float64,2},c_M_staples_fine::Array{Float64,2},P_S_staples_fine::Array{Float64,2},x_S_staples_fine::Array{Float64,2},
    λ_2_S_staples_fine::Array{Float64,2},unfeasible_mat_fine::Array{Float64,2},Y_S_potential_fine::Array{Float64,2},F_W::Float64,F_S::Float64,F_B::Float64,
    FM_W::Float64,FM_S::Float64,FM_B::Float64,TC_mat_fine::Array{Float64,2},a_max::Float64,
    P_kron1::SparseMatrixCSC{Float64,Int64},Q_trans::SparseMatrixCSC{Float64,Int64},C_max_staple_fine::Array{Float64,1},C_min_staple_fine::Array{Float64,1},C_max_staple_constrained_fine::Array{Float64,1},
    C_min_staple_constrained_fine::Array{Float64,1},TC_S_c3_constrained_fine::Array{Float64,1},
    x_S_c3_constrained_fine::Array{Float64,1},q_S_c3_constrained_fine::Array{Float64,1},c_S_c3_constrained_fine::Array{Float64,1},fspace_a_fine::Dict{Symbol,Any},
    x_S_mat_3c_fine::Array{Float64,2},x_B_mat_3c_fine::Array{Float64,2},land_B_mat_3c_fine::Array{Float64,2},λ_2_mat_3c_fine::Array{Float64,2},TC_mat_3c_fine::Array{Float64,2})
    #tol= 1e-8
    x_tmp = zeros(ns_fine,3);
    V_tmp = zeros(ns_fine,3);
    V_saved = zeros(ns_fine,3);
    cons_tmp = zeros(ns_fine,3);
    Phi_prime_fine = Array{SparseMatrixCSC{Float64,Int64}}(undef, 3);
    Q_tmp_block = Array{SparseMatrixCSC{Float64,Int64}}(undef, 3,1);
    cons_fine_local = zeros(ns_fine,3);
    a_prime_fine_local = zeros(ns_fine,3);
    future_occupation_fine_local = zeros(ns_fine,3);
    for j = 1:3
        for jj =1:3
            max_x_tmp = max_C_applied_fine[:,jj,j];
            min_x_tmp = min_C_applied_fine[:,jj,j];
            xtmp_loc,fval = goldenx(BellmanF_stable_transaction,min_x_tmp,
                max_x_tmp,tol,coeff,Phi_z_fine,fspace_a,β_loc,σ,jj,j,P_W_fine,Y_W_fine,s_fine,r,ρ,w,
                coeff_λ_2_cashcrop_residual_unconstrained_fine,coeff_λ_2_cashcrop_residual_constrained_fine,θ_fine,
                fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,
                a_min,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
                coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,
                labor_allocated_interior_c3a_fine,λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,
                x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,q_S_c2_fine,
                q_B_c1_fine,q_S_c3a_fine,q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,
                x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,Y_B_c3b_fine,
                c_S_mat_fine,c_B_mat_fine,c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,q_B_mat_fine,land_B_mat_fine,
                λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,feasibility_mat_fine,C_max_mat_fine,C_min_mat_fine,
                q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,
                x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine,F_W,F_S,F_B,
                FM_W,FM_S,FM_B,TC_mat_fine,a_max,C_max_staple_fine,C_min_staple_fine,
                C_max_staple_constrained_fine,C_min_staple_constrained_fine,TC_S_c3_constrained_fine,x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine,
                x_S_mat_3c_fine,x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine);
            V_tmp[:,jj],Phi_prime,x_tmp[:,jj] = BellmanF_stable_transaction(xtmp_loc,coeff,Phi_z_fine,fspace_a,β_loc,σ,jj,j,P_W_fine,Y_W_fine,s_fine,r,ρ,w,
            coeff_λ_2_cashcrop_residual_unconstrained_fine,coeff_λ_2_cashcrop_residual_constrained_fine,θ_fine,
            fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,
            a_min,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
            coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,
            labor_allocated_interior_c3a_fine,λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,
            x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,q_S_c2_fine,
            q_B_c1_fine,q_S_c3a_fine,q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,
            x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,Y_B_c3b_fine,
            c_S_mat_fine,c_B_mat_fine,c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,q_B_mat_fine,land_B_mat_fine,
            λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,feasibility_mat_fine,C_max_mat_fine,C_min_mat_fine,
            q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,
            x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine,F_W,F_S,F_B,
            FM_W,FM_S,FM_B,TC_mat_fine,a_max,C_max_staple_fine,C_min_staple_fine,
            C_max_staple_constrained_fine,C_min_staple_constrained_fine,TC_S_c3_constrained_fine,
            x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine,x_S_mat_3c_fine,x_B_mat_3c_fine,
            land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine,1);
            #(V_tmp[:,jj,j_cost])= BellmanF_stable(x_tmp[:,jj,j_cost],coeff,
            #    Phi_z_fine,fspace_a,income_mat_fine[:,jj,j,j_cost],β_loc,jj,price_final_current);
            Phi_prime_fine_a_tmp = funbase(fspace_a_fine, x_tmp[:,jj]);
            #println(typeof(P_kron_fine),typeof(Phi_prime_fine_a_tmp))
            #check_nonzero = Phi_prime_fine_a_tmp.>tol;
            #Phi_prime_fine_a_tmp = check_nonzero.*Phi_prime_fine_a_tmp;
            #println(nnz(Phi_prime_fine_a_tmp))
            #Phi_prime_fine_local =  row_kron(P_kron_fine,Phi_prime_fine_a_tmp);
            #check_nonzero = Phi_prime_fine_local.>tol;
            #Phi_prime_fine[jj,j_cost] = Phi_prime_fine_local.*check_nonzero;
            Phi_prime_fine[jj] =  row_kron(P_kron_fine,Phi_prime_fine_a_tmp);
            cons_tmp[:,jj] = xtmp_loc;
        end
        V_tmp[isnan.(V_tmp)] .= -10000;
        #V_tmp[isinf.(V_tmp)] .= -10000;
        V_next, future_occupation = findmax(V_tmp; dims=2);
        future_occupation_index = getindex.(future_occupation,2);
        cons = zeros(ns_fine,1);
        a_prime_local = zeros(ns_fine,1);
        V_saved[:,j]=V_next
        for jj=1:3
            future_occupation_tmp_dummy = sparse(hcat(float(
                        future_occupation_index[:,1].==jj)));
            Q_tmp_block[jj,1] =  row_kron(future_occupation_tmp_dummy,Phi_prime_fine[jj,1]);
            cons = cons + future_occupation_tmp_dummy.*cons_tmp[:,jj,1];
            a_prime_local = a_prime_local + future_occupation_tmp_dummy.*x_tmp[:,jj,1];
        end
        cons_fine_local[:,j] = cons;
        a_prime_fine_local[:,j] = a_prime_local;
        future_occupation_fine_local[:,j] = future_occupation_index;
        j_prev_ns_fine = (j-1)*ns_fine+1;
        j_ns_fine = j*ns_fine;
        Q_trans[j_prev_ns_fine:j_ns_fine,1:ns_fine] = Q_tmp_block[1,1];
        Q_trans[j_prev_ns_fine:j_ns_fine,(ns_fine+1):(2*ns_fine)]= Q_tmp_block[2,1];
        Q_trans[j_prev_ns_fine:j_ns_fine,(2*ns_fine+1):(3*ns_fine)] = Q_tmp_block[3,1];
    end
    check_nonzero = Q_trans.>tol;
    Q_trans = Q_trans.*check_nonzero;
    I1,J1,V1 = findnz(Q_trans);
    I2 = [I1;3 * ns_fine];
    J2 = [J1;3 * ns_fine];
    V2 = [V1;0.0];
    # I2=I1;
    # J2=J1;
    # V2=V1;
    Q_trans_prime = sparse(J2,I2,V2);
    return Q_trans_prime,cons_fine_local,a_prime_fine_local,future_occupation_fine_local,V_saved
end

function predict_irreducibility(future_occupation_fine_local::Array{Float64,2},
    exitflag_tmp::Int64)
    from_worker_trans_up = maximum(future_occupation_fine_local[:,1,1]);
    from_staple_trans_up = maximum(future_occupation_fine_local[:,2,1]);
    from_staple_trans_down = minimum(future_occupation_fine_local[:,2,1]);
    from_cash_trans_down = minimum(future_occupation_fine_local[:,3,1]);
    if (from_worker_trans_up<3 && from_staple_trans_up<3)
#        println("No cash crop producer in the economy")
        exitflag_tmp = exitflag_tmp + 1;
    elseif from_worker_trans_up<2
#        println("No rural population")
        exitflag_tmp = exitflag_tmp + 2;
    elseif from_staple_trans_down>1 && from_cash_trans_down>1
#        println("No urban population")
            # Only stop if the imprecision is very high
        exitflag_tmp = exitflag_tmp + 3;
    else
        # Final check: unstructured irredcucible
        #        Q_trans = transpose(Q_trans_prime); # tradeoff between additional
        #        Q_trans1 = Q_trans./sum(Q_trans,dims =2);
        #        mc = MarkovChain(Q_trans1)
        #        if QuantEcon.is_irreducible(mc)==true
        exitflag_tmp = 0;
        #        else
        #            exitflag_tmp = exitflag_tmp + 4;
        #        end
    end
    return exitflag_tmp
end

function stationary_distribution(Q_trans_prime::SparseMatrixCSC{Float64,Int64},
    ns_fine::Int64,exitflag_tmp::Int64,tol::Float64 = 1e-9)
    stat_distr=zeros(Float64,ns_fine)
    try
        lambda,V,eigflag = eigs(Q_trans_prime,nev=60, v0 = ones(3*ns_fine)./(3*ns_fine), which=:LR,maxiter = 1000);
        if (real(lambda[1] - lambda[2]) < tol)
            println("Unstable solution: second dominant eigenvalue is close to 1")
        end
        stat_distr = real(V[:,1])/sum(real(V[:,1]));
        stat_distr = stat_distr .* (stat_distr.>1e-12);

        return stat_distr,exitflag_tmp
    catch y
        #warn("Exception: ", y) # What to do on error.
        #if String(string(y))==".ARPACKException(1)"

        #KAROL EDIT:
        trans=Q_trans_prime
           probst = ones(3*ns_fine)./(3*ns_fine)
           test=1
           iter=0
           while test > 1e-9 && iter<1000
               iter=iter+1
              probst1 = trans*probst
              test = maximum(abs.(probst1-probst))
              probst=probst1
           end
           # if iter==1000
           #     println("St st dist didnt converge")
           # end
        stat_distr=probst
        exitflag_tmp=0

        # if isa(y, Arpack.ARPACKException) && iter==1000
        if iter==1000
#            println("Other issues with the stationary distribution")
            exitflag_tmp = 5;
            return ones(3*ns_fine)./(3*ns_fine),exitflag_tmp
        end
    end
    return stat_distr,exitflag_tmp
end

function policy_function_creator(C::Array{Float64,1},jj::Int64,j::Int64,P_W::Float64,Y_W::Array{Float64,1},s::Array{Float64,2},r::Float64,ρ::Float64,w::Float64,
    coeff_λ_2_cashcrop_residual_unconstrained::Array{Float64,2},coeff_λ_2_cashcrop_residual_constrained::Array{Float64,2},θ::Array{Float64,1},fspace_C_fine::Dict{Symbol,Any},
    ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,p_B::Float64,p_M::Float64,ϕ_B::Float64,τ_B::Float64,c̄_S::Float64,Q_S::Float64,ϵ::Float64,ψ_S::Float64,ψ_B::Float64,
    ψ_M::Float64,ns::Int64,κ::Float64,tol::Float64,a_min::Float64,x_B_c1::Array{Float64,1},π_B_only_B_c1::Array{Float64,1},λ_B_only_B_c1::Array{Float64,1},P_B_c1::Float64,
    Y_B_c1::Array{Float64,1},coeff_λ_2_s::Array{Float64,2},P_S_c1::Float64,P_S_c2::Float64,Y_S_c1::Array{Float64,1},Y_S_c2::Array{Float64,1},x_S_c1::Array{Float64,1},
    x_S_c2::Array{Float64,1},labor_allocated_interior_c3a::Float64,λ_B_interior_c3a::Array{Float64,1},x_SC_interior_c3a::Array{Float64,1},
    x_BC_interior_c3a::Array{Float64,1},Y_B_c3a::Array{Float64,1},P_B_c3a::Float64,P_B_c3b::Float64,q_S_c1::Array{Float64,1},q_S_c2::Array{Float64,1},
    q_B_c1::Array{Float64,1},q_S_c3a::Array{Float64,1},q_B_c3a::Array{Float64,1},q_S_c3b::Array{Float64,1},q_B_c3b::Array{Float64,1},x_SC_interior_c3b::Array{Float64,1},
    x_BC_interior_c3b::Array{Float64,1},labor_allocated_interior_c3b::Float64,Y_B_c3b::Array{Float64,1},c_S_mat::Array{Float64,2},c_B_mat::Array{Float64,2},
    c_M_mat::Array{Float64,2},x_S_mat::Array{Float64,2},x_B_mat::Array{Float64,2},q_S_mat::Array{Float64,2},q_B_mat::Array{Float64,2},land_B_mat::Array{Float64,2},
    λ_2_mat::Array{Float64,2},P_B_mat::Array{Float64,2},Y_B_mat::Array{Float64,2},feasibility_mat::Array{Float64,2},C_max_mat::Array{Float64,2},C_min_mat::Array{Float64,2},
    q_S_staples::Array{Float64,2},c_S_staples::Array{Float64,2},c_B_staples::Array{Float64,2},c_M_staples::Array{Float64,2},P_S_staples::Array{Float64,2},
    x_S_staples::Array{Float64,2},λ_2_S_staples::Array{Float64,2},unfeasible_mat::Array{Float64,2},Y_S_potential::Array{Float64,2},F_W::Float64,F_S::Float64,F_B::Float64,
    FM_W::Float64,FM_S::Float64,FM_B::Float64,TC_mat::Array{Float64,2},C_max_staple::Array{Float64,1},
    C_min_staple::Array{Float64,1},C_max_staple_constrained::Array{Float64,1},
    C_min_staple_constrained::Array{Float64,1},TC_S_c3_constrained::Array{Float64,1},
    x_S_c3_constrained::Array{Float64,1},q_S_c3_constrained::Array{Float64,1},c_S_c3_constrained::Array{Float64,1},
    x_S_mat_3c::Array{Float64,2},x_B_mat_3c::Array{Float64,2},land_B_mat_3c::Array{Float64,2},λ_2_mat_3c::Array{Float64,2},TC_mat_3c::Array{Float64,2})
    if jj == 1
        future_asset = (1 + r ) * s[:,1] -P_W *C  -Y_W  .- w*FM_W .- w*F_W * (j != 1);
        return future_asset
    elseif jj == 2
        Y_S,c_S_S,c_B_S,c_M_S,q_S_S,P_S,x_S_S,solve_staple_index_S,λ_2_S = staples_objects(C,ϕ_S,τ_S,p_x,
            p_B,p_M,ϕ_B,τ_B,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ζ,θ,coeff_λ_2_s,fspace_C_fine,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,
            x_S_c1, x_S_c2,s,q_S_c1,q_S_c2,q_S_staples,c_S_staples,c_B_staples,c_M_staples,P_S_staples,x_S_staples,
            λ_2_S_staples,unfeasible_mat,Y_S_potential,κ,c̄_S,C_max_staple,
            C_min_staple,ns,C_max_staple_constrained,
            C_min_staple_constrained,TC_S_c3_constrained,
            x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained,a_min);
        future_asset = (1 + r ) * s[:,1] - Y_S[:]  .- w*FM_S .- w*F_S * (j == 1);
        return future_asset,Y_S[:],c_S_S[:],c_B_S[:],c_M_S[:],q_S_S[:],P_S[:],x_S_S[:],getindex.(solve_staple_index_S,2)[:],λ_2_S[:]
    elseif jj == 3
        (c_S_B,c_B_C,c_M_B,x_SC,x_BC,land_C,λ_2,P_B,Y_B,q_S_C,q_B_B,solve_cash_crop_index_B,solve_staple_index_B,TC_B) = cashcrop_objects(C,
                coeff_λ_2_cashcrop_residual_unconstrained,coeff_λ_2_cashcrop_residual_constrained,
                    θ,fspace_C_fine,s,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns,κ,tol,a_min,x_B_c1,π_B_only_B_c1,λ_B_only_B_c1,P_B_c1,Y_B_c1,
                        coeff_λ_2_s,P_S_c1,P_S_c2,Y_S_c1,Y_S_c2,x_S_c1, x_S_c2,labor_allocated_interior_c3a,λ_B_interior_c3a,x_SC_interior_c3a,
                        x_BC_interior_c3a,Y_B_c3a,P_B_c3a,P_B_c3b,q_S_c1,q_S_c2,q_B_c1,q_S_c3a,q_B_c3a,q_S_c3b,q_B_c3b,
                        x_SC_interior_c3b,x_BC_interior_c3b,labor_allocated_interior_c3b,Y_B_c3b,c_S_mat,c_B_mat,c_M_mat,x_S_mat,x_B_mat,q_S_mat,q_B_mat,land_B_mat,
                        λ_2_mat,P_B_mat,Y_B_mat,feasibility_mat,C_max_mat,C_min_mat,TC_mat,q_S_staples,c_S_staples,c_B_staples,c_M_staples,P_S_staples,
                        x_S_staples,λ_2_S_staples,unfeasible_mat,Y_S_potential,ρ,C_max_staple,
                        C_min_staple,C_max_staple_constrained,
                        C_min_staple_constrained,TC_S_c3_constrained,
                        x_S_c3_constrained,q_S_c3_constrained,c_S_c3_constrained, x_S_mat_3c,x_B_mat_3c,land_B_mat_3c,λ_2_mat_3c,TC_mat_3c);
        future_asset = (1 + r ) * s[:,1] - Y_B[:]  .- w*FM_B .- w*F_B * (j != 3);
        return future_asset,c_S_B[:],c_B_C[:],c_M_B[:],x_SC[:],x_BC[:],land_C[:],λ_2[:],P_B[:],Y_B[:],q_S_C[:],q_B_B[:],getindex.(solve_cash_crop_index_B,2)[:],getindex.(solve_staple_index_B,2)[:],TC_B
    end
end

function RCT_opt_cons(coeff::Array{Float64,2},
    ns_fine::Int64,P_kron_fine::SparseMatrixCSC{Float64,Int64},
    min_C_applied_fine::Array{Float64,3},max_C_applied_fine::Array{Float64,3},Phi_z_fine::SparseMatrixCSC{Float64,Int64},β_loc::Float64,
    fspace_a::Dict{Symbol,Any},σ::Float64,P_W_fine::Float64,Y_W_fine::Array{Float64,1},s_fine::Array{Float64,2},r::Float64,ρ::Float64,w::Float64,
    coeff_λ_2_cashcrop_residual_unconstrained_fine::Array{Float64,2},coeff_λ_2_cashcrop_residual_constrained_fine::Array{Float64,2},θ_fine::Array{Float64,1},
    fspace_C_fine::Dict{Symbol,Any},ϕ_S::Float64,ζ::Float64,τ_S::Float64,p_x::Float64,p_B::Float64,p_M::Float64,ϕ_B::Float64,τ_B::Float64,c̄_S::Float64,Q_S::Float64,
    ϵ::Float64,ψ_S::Float64,ψ_B::Float64,ψ_M::Float64,κ::Float64,tol::Float64,a_min::Float64,x_B_c1_fine::Array{Float64,1},π_B_only_B_c1_fine::Array{Float64,1},
    λ_B_only_B_c1_fine::Array{Float64,1},P_B_c1_fine::Float64,Y_B_c1_fine::Array{Float64,1},coeff_λ_2_s_fine::Array{Float64,2},P_S_c1_fine::Float64,P_S_c2_fine::Float64,
    Y_S_c1_fine::Array{Float64,1},Y_S_c2_fine::Array{Float64,1},x_S_c1_fine::Array{Float64,1},x_S_c2_fine::Array{Float64,1},labor_allocated_interior_c3a_fine::Float64,
    λ_B_interior_c3a_fine::Array{Float64,1},x_SC_interior_c3a_fine::Array{Float64,1},x_BC_interior_c3a_fine::Array{Float64,1},Y_B_c3a_fine::Array{Float64,1},
    P_B_c3a_fine::Float64,P_B_c3b_fine::Float64,q_S_c1_fine::Array{Float64,1},q_S_c2_fine::Array{Float64,1},q_B_c1_fine::Array{Float64,1},q_S_c3a_fine::Array{Float64,1},
    q_B_c3a_fine::Array{Float64,1},q_S_c3b_fine::Array{Float64,1},q_B_c3b_fine::Array{Float64,1},x_SC_interior_c3b_fine::Array{Float64,1},x_BC_interior_c3b_fine::Array{Float64,1},
    labor_allocated_interior_c3b_fine::Float64,Y_B_c3b_fine::Array{Float64,1},c_S_mat_fine::Array{Float64,2},c_B_mat_fine::Array{Float64,2},
    c_M_mat_fine::Array{Float64,2},x_S_mat_fine::Array{Float64,2},x_B_mat_fine::Array{Float64,2},q_S_mat_fine::Array{Float64,2},q_B_mat_fine::Array{Float64,2},
    land_B_mat_fine::Array{Float64,2},λ_2_mat_fine::Array{Float64,2},P_B_mat_fine::Array{Float64,2},Y_B_mat_fine::Array{Float64,2},feasibility_mat_fine::Array{Float64,2},
    C_max_mat_fine::Array{Float64,2},C_min_mat_fine::Array{Float64,2},q_S_staples_fine::Array{Float64,2},c_S_staples_fine::Array{Float64,2},
    c_B_staples_fine::Array{Float64,2},c_M_staples_fine::Array{Float64,2},P_S_staples_fine::Array{Float64,2},x_S_staples_fine::Array{Float64,2},
    λ_2_S_staples_fine::Array{Float64,2},unfeasible_mat_fine::Array{Float64,2},Y_S_potential_fine::Array{Float64,2},F_W::Float64,F_S::Float64,F_B::Float64,
    FM_W::Float64,FM_S::Float64,FM_B::Float64,TC_mat_fine::Array{Float64,2},a_max::Float64,
    P_kron1::SparseMatrixCSC{Float64,Int64},Q_trans::SparseMatrixCSC{Float64,Int64},C_max_staple_fine::Array{Float64,1},C_min_staple_fine::Array{Float64,1},C_max_staple_constrained_fine::Array{Float64,1},
    C_min_staple_constrained_fine::Array{Float64,1},TC_S_c3_constrained_fine::Array{Float64,1},
    x_S_c3_constrained_fine::Array{Float64,1},q_S_c3_constrained_fine::Array{Float64,1},c_S_c3_constrained_fine::Array{Float64,1},
    x_S_mat_3c_fine::Array{Float64,2},x_B_mat_3c_fine::Array{Float64,2},land_B_mat_3c_fine::Array{Float64,2},λ_2_mat_3c_fine::Array{Float64,2},TC_mat_3c_fine::Array{Float64,2})
    #tol= 1e-8
    x_tmp = zeros(ns_fine,3);
    V_tmp = zeros(ns_fine,3);
    cons_tmp = zeros(ns_fine,3);
    Phi_prime_fine = Array{SparseMatrixCSC{Float64,Int64}}(undef, 3);
    Q_tmp_block = Array{SparseMatrixCSC{Float64,Int64}}(undef, 3,1);
    cons_fine_local = zeros(ns_fine,3);
    a_prime_fine_local = zeros(ns_fine,3);
    future_occupation_fine_local = zeros(ns_fine,3);
    for j = 1:3
        for jj =1:3
            max_x_tmp = max_C_applied_fine[:,jj,j];
            min_x_tmp = min_C_applied_fine[:,jj,j];
            xtmp_loc,fval = goldenx(BellmanF_stable_transaction,min_x_tmp,
                max_x_tmp,tol,coeff,Phi_z_fine,fspace_a,β_loc,σ,jj,j,P_W_fine,Y_W_fine,s_fine,r,ρ,w,
                coeff_λ_2_cashcrop_residual_unconstrained_fine,coeff_λ_2_cashcrop_residual_constrained_fine,θ_fine,
                fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,
                a_min,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
                coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,
                labor_allocated_interior_c3a_fine,λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,
                x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,q_S_c2_fine,
                q_B_c1_fine,q_S_c3a_fine,q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,
                x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,Y_B_c3b_fine,
                c_S_mat_fine,c_B_mat_fine,c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,q_B_mat_fine,land_B_mat_fine,
                λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,feasibility_mat_fine,C_max_mat_fine,C_min_mat_fine,
                q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,
                x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine,F_W,F_S,F_B,
                FM_W,FM_S,FM_B,TC_mat_fine,a_max,C_max_staple_fine,C_min_staple_fine,
                C_max_staple_constrained_fine,C_min_staple_constrained_fine,TC_S_c3_constrained_fine,x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine,
                x_S_mat_3c_fine,x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine);
            V_tmp[:,jj],Phi_prime,x_tmp[:,jj] = BellmanF_stable_transaction(xtmp_loc,coeff,Phi_z_fine,fspace_a,β_loc,σ,jj,j,P_W_fine,Y_W_fine,s_fine,r,ρ,w,
            coeff_λ_2_cashcrop_residual_unconstrained_fine,coeff_λ_2_cashcrop_residual_constrained_fine,θ_fine,
            fspace_C_fine,ϕ_S,ζ,τ_S,p_x,p_B,p_M,ϕ_B,τ_B,c̄_S,Q_S,ϵ,ψ_S,ψ_B,ψ_M,ns_fine,κ,tol,
            a_min,x_B_c1_fine,π_B_only_B_c1_fine,λ_B_only_B_c1_fine,P_B_c1_fine,Y_B_c1_fine,
            coeff_λ_2_s_fine,P_S_c1_fine,P_S_c2_fine,Y_S_c1_fine,Y_S_c2_fine,x_S_c1_fine, x_S_c2_fine,
            labor_allocated_interior_c3a_fine,λ_B_interior_c3a_fine,x_SC_interior_c3a_fine,
            x_BC_interior_c3a_fine,Y_B_c3a_fine,P_B_c3a_fine,P_B_c3b_fine,q_S_c1_fine,q_S_c2_fine,
            q_B_c1_fine,q_S_c3a_fine,q_B_c3a_fine,q_S_c3b_fine,q_B_c3b_fine,
            x_SC_interior_c3b_fine,x_BC_interior_c3b_fine,labor_allocated_interior_c3b_fine,Y_B_c3b_fine,
            c_S_mat_fine,c_B_mat_fine,c_M_mat_fine,x_S_mat_fine,x_B_mat_fine,q_S_mat_fine,q_B_mat_fine,land_B_mat_fine,
            λ_2_mat_fine,P_B_mat_fine,Y_B_mat_fine,feasibility_mat_fine,C_max_mat_fine,C_min_mat_fine,
            q_S_staples_fine,c_S_staples_fine,c_B_staples_fine,c_M_staples_fine,P_S_staples_fine,
            x_S_staples_fine,λ_2_S_staples_fine,unfeasible_mat_fine,Y_S_potential_fine,F_W,F_S,F_B,
            FM_W,FM_S,FM_B,TC_mat_fine,a_max,C_max_staple_fine,C_min_staple_fine,
            C_max_staple_constrained_fine,C_min_staple_constrained_fine,TC_S_c3_constrained_fine,
            x_S_c3_constrained_fine,q_S_c3_constrained_fine,c_S_c3_constrained_fine,x_S_mat_3c_fine,
            x_B_mat_3c_fine,land_B_mat_3c_fine,λ_2_mat_3c_fine,TC_mat_3c_fine,1);
            #(V_tmp[:,jj,j_cost])= BellmanF_stable(x_tmp[:,jj,j_cost],coeff,
            #    Phi_z_fine,fspace_a,income_mat_fine[:,jj,j,j_cost],β_loc,jj,price_final_current);
            #Phi_prime_fine_a_tmp = funbase(fspace_a_fine, x_tmp[:,jj]);
            #println(typeof(P_kron_fine),typeof(Phi_prime_fine_a_tmp))
            #check_nonzero = Phi_prime_fine_a_tmp.>tol;
            #Phi_prime_fine_a_tmp = check_nonzero.*Phi_prime_fine_a_tmp;
            #println(nnz(Phi_prime_fine_a_tmp))
            #Phi_prime_fine_local =  row_kron(P_kron_fine,Phi_prime_fine_a_tmp);
            #check_nonzero = Phi_prime_fine_local.>tol;
            #Phi_prime_fine[jj,j_cost] = Phi_prime_fine_local.*check_nonzero;
            #Phi_prime_fine[jj] =  row_kron(P_kron_fine,Phi_prime_fine_a_tmp);
            cons_tmp[:,jj] = xtmp_loc;
        end
        V_tmp[isnan.(V_tmp)] .= -10000;
        #V_tmp[isinf.(V_tmp)] .= -10000;
        V_next, future_occupation = findmax(V_tmp; dims=2);
        future_occupation_index = getindex.(future_occupation,2);
        cons = zeros(ns_fine,1);
        a_prime_local = zeros(ns_fine,1);
        for jj=1:3
            future_occupation_tmp_dummy = sparse(hcat(float(
                        future_occupation_index[:,1].==jj)));
            #Q_tmp_block[jj,1] =  row_kron(future_occupation_tmp_dummy,Phi_prime_fine[jj,1]);
            cons = cons + future_occupation_tmp_dummy.*cons_tmp[:,jj,1];
            a_prime_local = a_prime_local + future_occupation_tmp_dummy.*x_tmp[:,jj,1];
        end
        cons_fine_local[:,j] = cons;
        a_prime_fine_local[:,j] = a_prime_local;
        future_occupation_fine_local[:,j] = future_occupation_index;
        #j_prev_ns_fine = (j-1)*ns_fine+1;
        #j_ns_fine = j*ns_fine;
        #Q_trans[j_prev_ns_fine:j_ns_fine,1:ns_fine] = Q_tmp_block[1,1];
        #Q_trans[j_prev_ns_fine:j_ns_fine,(ns_fine+1):(2*ns_fine)]= Q_tmp_block[2,1];
        #Q_trans[j_prev_ns_fine:j_ns_fine,(2*ns_fine+1):(3*ns_fine)] = Q_tmp_block[3,1];
    end
    #check_nonzero = Q_trans.>tol;
    #Q_trans = Q_trans.*check_nonzero;
    #I1,J1,V1 = findnz(Q_trans);
    #I2 = [I1;3 * ns_fine];
    #J2 = [J1;3 * ns_fine];
    #V2 = [V1;0.0];
    # I2=I1;
    # J2=J1;
    # V2=V1;
    #Q_trans_prime = sparse(J2,I2,V2);
    return cons_fine_local,a_prime_fine_local,future_occupation_fine_local# Q_trans_prime,cons_fine_local,a_prime_fine_local,future_occupation_fine_local
end

function cbar_creator(parameters_tmp::Parameter_type)
    x_S_minimum = parameters_tmp.κ*parameters_tmp.a_min/(parameters_tmp.p_x);  # Need a box for the actual manufacturing price! Assume that working capital binds for the poorest! 
    q_S_minimum = parameters_tmp.ϕ_S * parameters_tmp.z[1] .*x_S_minimum.^parameters_tmp.ζ;
    fertilizer_exp_minimum =  parameters_tmp.p_x * x_S_minimum;
    return -parameters_tmp.a_min .+ q_S_minimum .- fertilizer_exp_minimum 
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