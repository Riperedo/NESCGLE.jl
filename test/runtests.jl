using Test
using NESCGLE

# k-vector
Nk = 200
kmax = 15*π; dk = kmax/Nk
k = dk*(collect(1:Nk) .- 0.5)
ϕ = 0.5
T = 0.1
ν = 6.0
sm = SM_HS(ϕ, k)
#sm = SM_WCA(ϕ, T, k)
#sm = update_SM(sm, [0.6])
#S = structure_factor(sm)
#iterations, gammas, system = Asymptotic(ϕ, k, S, flag = true)
#pp = StaticProcess(sm.params)
#pp = StaticProcess([ϕ, T, ν])
pp = InstantaneousProcess([ϕ, 1.0], [0.6, 1.0])
#pp = InstantaneousProcess([ϕ, T, ν], [0.7, T, ν])
#pp = FiniteRate([ϕ, T, ν], [0.65, T, ν], 100.0)
#pp = Hysteresis([ϕ, T, ν], [0.65, T, ν], 1.0)
#save_data("test.dat", [k S], header="k\tS")
#τ, Fs, F, Δζ, ΔG, Δr² = SCGLE(ϕ, k, S)#; dp = dyn_params(), k_max = 7.2)
#save_data("test.dat", [τ Fs F Δζ ΔG Δr²])
#save_data("test.dat", [pp.time [pp.process[n][2] for n in 1:length(pp.time)]])
NESCGLEsolver(sm, pp, so = saving_options(false, false, true))
@test true