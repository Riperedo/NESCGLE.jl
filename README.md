![Logo](docs/build/assets/Logo.png)

# NESCGLE.jl Project

Welcome to the NESCGLE.jl project! This repository, written in [Julia](https://julialang.org/), contains a set of routines designed to help users, whether new to liquid theory or experienced, evaluate the equations of the [NE-SCGLE](https://doi.org/10.1103/PhysRevE.82.061503) theoretical framework.

The code is designed with simplicity and ease of use in mind, operating through intuitive APIs.

## Installation

The easiest way to install this package is through the Julia package manager. To do this, press `]` in the Julia REPL to enter the package manager environment, and then run the following command:

```julia
add "https://github.com/Riperedo/NESCGLE.jl.git"
```

## Usage
Here's a basic example to get you started:

### Static Structure Factor

```julia
using NESCGLE

Nk = 200; kmax = 15*π; dk = kmax/Nk
k = dk*(collect(1:Nk) .- 0.5)
ϕ = 0.5
sm = SM_HS(ϕ, k)
S = structure_factor(sm)
save_data("sdk.dat", [k S], header="k\tS")
```

### Asymptotics

```julia
using NESCGLE

Nk = 200; kmax = 15*π; dk = kmax/Nk
k = dk*(collect(1:Nk) .- 0.5)
ϕ = 0.5
sm = SM_HS(ϕ, k)
S = structure_factor(sm)
iterations, gammas, system = Asymptotic(ϕ, k, S, flag = true)
```

### Dynamics

```julia
using NESCGLE

Nk = 200; kmax = 15*π; dk = kmax/Nk
k = dk*(collect(1:Nk) .- 0.5)
ϕ = 0.5
sm = SM_HS(ϕ, k)
S = structure_factor(sm)
τ, Fs, F, Δζ, ΔG, Δr² = SCGLE(ϕ, k, S)#; dp = dyn_params(), k_max = 7.2
save_data("test.dat", [τ Fs F Δζ ΔG Δr²], header = "τ\tFs\tF\tΔζ\tΔG\tΔr²")
```

## Preparation protocols

We had included a set of different  preparation protocols:
* `StaticProcess` for an equivalent equilibrium version of the NESCGLE theory.
* `InstantaneousProcess` for an instantaneous quench from an initial thermodynamic state to a target thermodynamic state.
* `FiniteRate` for a finite rate of change in the thermodynamic parameters.
* `Hysteresis` for a thermal hysteresis process.

```julia
using NESCGLE

Nk = 200; kmax = 15*π; dk = kmax/Nk
k = dk*(collect(1:Nk) .- 0.5)
ϕi = 0.5
ϕf = 0.6
sm = SM_HS(ϕi, k)
ℇ = SM_HS(ϕf, k)
pp = InstantaneousProcess(sm.params, ℇ.params)
sol = NESCGLEsolver(sm, pp)
save_files(sol, sm, pp)
```

## Features

* Stability Matrix Computation: Compute the thermodynamic stability matrix for a given colloidal system.
* Preparation Protocol: Define the preparation protocol for different physical scenarios.
* Dynamics Evaluation: Solve the NESCGLE equations to obtain transport properties.
* Flexible Saving Options: Choose what information to save, including structure factors, transport properties, and waiting times.

## Contributing

Contributions are welcome! If you have any suggestions, bug reports, or pull requests, feel free to submit them on our GitHub repository.

## Licence

This project is licensed under the MIT License. See the LICENSE file for details.

## Acknowledgments

We would like to thank the developers of Julia and the scientific community for providing the tools and inspiration to make this project possible.
