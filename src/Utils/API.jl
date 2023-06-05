include("gdr.jl")

using DelimitedFiles

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

function make_directory(path :: String)
	if !isdir(path) #if there are not directory "path"
		mkdir(path) #then make it
	end
	path *= "/"
	#println(path)
	return path
end

function num2text(x :: Float64)
	txt = ""
	texto = string(x)
	if length(texto) > 8
		txt = texto[1:8]
	else
		txt = texto
	end
	if '.' in txt
		txt = replace(txt, "." => "p")
	end
	return txt
end

#index = findfirst(item-> item !=0, params_i-params_f)
