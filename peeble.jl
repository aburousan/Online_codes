using DifferentialEquations
using Plots, LaTeXStrings
using Interpolations, Roots

zeta_3 = 1.20205690316
#All in eV
EI = 13.6 #Hydrogen Binding,
Me = 5.10999e5 #Mass of electron
zeq = 3395 # Redshift for matter-radiation equality
H0 = 1.5e-33
αₛ = 1/137 #Fine structure Constant
Ωm = 0.3089 #Matter critical density",
Γ2s = 8.227*6.58e-16 # Two photon decay rate",
z0 = 2.725/(1.160451812e4)
η = 6e-10 #baryon to photon ratio
nb = 0.76*(2*η*zeta_3)/π^2#baryon no.density

T(z) = BigFloat(z0*(1 + z))#Conversion of z to T in eV

#For Equilibrium Part "Saha Equation"
pre_fac_saha = nb*(2*π/Me)^1.5
function Saha_Xe(z)
    f = pre_fac_saha*T(z)^(3/2)*BigFloat(exp(EI/T(z)))
    return (-1 + sqrt(1 + 4*f))/(2*f)
end
#---------------------------------
#Finding Where Xe = 0.5
soln_saha(z) = Saha_Xe(z) - 0.5
sol1 = find_zero(soln_saha, (1300,1400))#Z value where Xe = 0.5 for Saha
T(sol1)
#---------------------------------

# For Non-Equilibrium Part
function ab_beta(z)#Defining α, b and β
    t = T(z); x = EI/t
    alpha = BigFloat(9.8*(αₛ/Me)^2*sqrt(x)*log(x))
    b = alpha*(Me*t/(2*π))^(1.5)*BigFloat(exp(-x))
    beta = b*exp(0.75*x)
    return alpha, b, beta
end

H(z) = H0*sqrt(Ωm)*(1 + z)^(1.5)*sqrt(1 + (1+z)/(1 + zeq))

# For Boltzmann without peeble factor
function boltz_diff(Xe, p, z)
    prefac = H(z)*(1 + z)
    alpha, b, beta = ab_beta(z)
    fact_main = (alpha*nb*Xe^2 *T(z)^3) - b*(1-Xe)
    return fact_main/prefac
end
Xe0 = BigFloat(1)#Initial Value of Xe
zspan = (BigFloat(1800), BigFloat(0))#Range of soln
prob_boltz = ODEProblem(boltz_diff, Xe0, zspan)#Define the problem
sol_boltz = solve(prob_boltz, alg_hints = [:stiff])#Solution
z_bol = sort(sol_boltz.t)#solution Z values
Xe_bol = sort(sol_boltz.u)#solution Xe values
#-------------------------------------------------------------------------------------------------------------
inter_boltz = linear_interpolation(z_bol,Xe_bol)#Linear Interpolation
bolt_05(z) = inter_boltz(z) - 0.5 #For our simplification
sol1 = find_zero(bolt_05, (1100,1500))#Z value at which Xe 0.5 for Boltzmann
#-------------------------------------------------------------------------------------------------------------
T(1368)
#Defining peeble factor
function Cr_p(z,u)
    hh = H(z); temp = T(z)
    alpha, b, beta = ab_beta(z)
    Γ2p_nu = BigFloat(27*hh*EI^3)
    Γ2p_di_with_out_xe = BigFloat(64*π^2*nb*temp^3)*(1-u)
    Γ2p = BigFloat((Γ2p_nu/Γ2p_di_with_out_xe))
    Cr = BigFloat((Γ2s + Γ2p)/(Γ2s + Γ2p + beta))
    return Cr
end

#DifferentialEquation with peeble factor
function peeble_diff(Xe, p, z)
    prefac = H(z)*(1 + z)
    alpha, b, beta = ab_beta(z)
    fact_main = (alpha*nb*Xe^2 *T(z)^3) - b*(1-Xe)
    return Cr_p(z,Xe)*fact_main/prefac
end
Xe0 = 0.9999# Should be 1 but due to discontinuity 0.9999 is used
prob_peeble = ODEProblem(peeble_diff, Xe0, zspan)
sol_peeble = solve(prob_peeble, alg_hints = [:stiff])#solution peeble Eqn
z_peb = sort(sol_peeble.t)#z values for peeble soln
Xe_peb = sort(sol_peeble.u)#Xe values for peeble soln
cr_array = broadcast(Cr_p,z_peb,Xe_pebb)#Cr values
# Finding z for which Xe = 0.5 for peeble eqn
inter_peeb = linear_interpolation(z_peb,Xe_pebb)
pebb_05(z) = inter_peeb(z) - 0.5
sol3 = find_zero(pebb_05, (1100,1300))#z for which Xe = 0.5
T(1270)
#---------------------------------------------------------------------------------

#Plotting Xe for all cases
zs = range(0, 1800, length=10_000)
saha = Saha_Xe.(zs)
plot(yscale=:log10, minorgrid=true)
xlims!(0, 1800)
ylims!(1e-4, 1.5e+0)
xlabel!(L"Redshift ($z$)")
ylabel!(L"Free Electron Fraction ($X_e$)")
plot!([0,1800],[0.5,0.5],lc=:green,lw=2.5,label="",ls=:dash)#0.5 line
annotate!(700, 0.3, L"$X_e = 0.5$")
#Saha eqn plot
plot!(zs, saha, label="",lc=:black, lw=2)
annotate!(1100, 2*1e-4, "Saha", :black)
#Boltzmann Eqn plot
plot!(z_bol, Xe_bol, label="",lc=:red, lw=2)
annotate!(300, 2*1e-4, "Boltzmann", :red)
#Peeble eqn plot
plot!(z_peb, Xe_peb, label="",lc=:blue, lw=2)
annotate!(250, 9*1e-4, "Peebles", :blue)
savefig("saha_boltz_peeble.png")
# savefig("saha_boltzmann.png")

# Plotting Cr
plot(yscale=:log10,minorgrid=true)
xlims!(0, 1800)
ylims!(1e-4, 1.5e+0)
xlabel!(L"Redshift ($z$)")
ylabel!(L"Peeble Constant ($C_r$)")
plot!(z_peb,cr_array, minorgrid=true, label="", lw=2.2)
#Plotting the line for min Cr
c_min = minimum(cr_array)
plot!([0,1800],[c_min, c_min],label="",lw=2, ls=:dash)
annotate!(1000, c_min - c_min/3, L"$C_r = 1.827\times 10^{-3}$", :black)
savefig("peeble_fac.png")

