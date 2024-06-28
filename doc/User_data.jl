using NESCGLE
using JSON
using DelimitedFiles

function interpolation(X::Array{Float64}, Y::Array{Float64}, x)
    if x < X[1]
        return X[1]
    elseif x >= X[end]
        return 1.0
    else
        index = sum(X .< x)
        m = (Y[index] - Y[index+1])/(X[index] - X[index+1])
        b = Y[index] - m*X[index]
        return m*x + b
    end
end

function main(args...)
    # user inputs
    sk_path, ϕ_str, T_str =  args

    println(sk_path)
    println(ϕ_str)

    if isfile(sk_path)
        # reading input data
        data = readdlm(sk_path)
        k_input = data[:,1]
        S_input = data[:,2]
        # volume fraction input
        ϕ = parse(Float64, ϕ_str)
        # temperature input
        T = parse(Float64, T_str)
        # defining structure factor funtion
        S_function(q::Float64) = interpolation(k_input, S_input, q)
        # defining wave vector
        k = collect(0.0:0.1:15*π)
        I = Input_SCGLE([ϕ, T], S_function, k, "HS")

        # preparing saving folder
        save_path = NESCGLE.make_directory("userID")
        save_path = NESCGLE.make_directory(save_path*"T"*num2text(T))
        save_path = NESCGLE.make_directory(save_path*"phi"*num2text(ϕ))
        filename = save_path*"output.json"
        
        # Computer Static Sructure Factor
        S = structure_factor(I)
        # computing dynamics
        τ, Fs, F, Δζ, Δη, D, W = SCGLE(I)
        # parsing to json file
        structural_data = Dict("k"=>k, "S"=>S)
        #dynamics_data = Dict("tau"=>τ, "sISF"=>Fs, "ISF"=>F, "Dzeta"=>Δζ, "Deta"=>Δη, "D"=>D, "MSD"=>W)
        dynamics_data = Dict("tau"=>τ, "Deta"=>Δη, "MSD"=>W)
        data = Dict("Statics"=>structural_data, "Dynamics"=>dynamics_data)
        # saving data
        open(filename, "w") do file
            JSON.print(file, data)
        end
    end
    println("Calculation complete.")
end

@time main(ARGS...)
# julia script.jl ruta_de_sk.dat phi T
# julia .\user_data_Sk.jl .\path_to_data\Sk.dat 0.5 300
