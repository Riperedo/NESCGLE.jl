module NESCGLE

export InstantaneousProcess, cooling_rate, hysteresis, ChebyshevGrids, Cheb_plus_U, grid_plus_grid, Input_SCGLE, wave_vector, structure_factor, structure_factor_function, parameters, volume_fraction, system, Input_HS, Input_WCA, Input_SW,  Input_Yukawa, Input_SALR, Input_AO, SCGLE, Asymptotic, save_data, num2text, g_HS

include("NonInstantaneousProcess/CoolingRate/API.jl")
include("NonInstantaneousProcess/Hysteresis/API.jl")

end
