include("gdr_HS.jl")

using DelimitedFiles


"""
    save_data(nombre, formato; header = "", flag = true)
Simple function to save data, this must be in the form `[x y z ...]` if you include a coma a transpose is saved.

# Example
```
x = collect(1:10)
y = x.^2
save_data("example.dat", [x y])
```
"""
function save_data(nombre, formato; header = "", flag = true)
	@assert typeof(nombre) == typeof("hola") "El primer argumento debe ser texto"
	# el formato debe estar entre parentesis cuadrados [x, y, z, ...]
	# readdlm(source, delim::AbstractChar, T::Type, eol::AbstractChar; header=false, skipstart=0, skipblanks=true, use_mmap, quotes=true, dims, comments=false, comment_char='#')
	open(nombre, "w") do io
		if header != ""
			write(io, "# "*header*"\n")
		end
		writedlm(io, formato)
	end
	if flag	println("Data saved as ", nombre) end
end

const exception_dir = ["np", "nu", "VW"]

"""
    make_directory(path::String) -> String

Create a directory if it does not already exist, and return the path with a trailing slash.

# Arguments
- `path::String`: The path of the directory to create.

# Returns
- `String`: The path of the directory, with a trailing slash.

# Example
```julia
dir_path = make_directory("data/experiment1")
println("Directory path: ", dir_path)
```
"""
function make_directory(path :: String)
	if !isdir(path) #if there are not directory "path"
		mkdir(path) #then make it
	end
	path *= "/"
	#println(path)
	return path
end
