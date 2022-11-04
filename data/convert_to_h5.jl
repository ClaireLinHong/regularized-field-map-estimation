# convert_to_h5.jl
# convert data from proprietary .mat format to open .h5 format

using MAT
using HDF5

# version for file with multiple variables, splitting each into its own file
mfile = "input_object_40sl_3d_epi_snr40.mat"
alldata = matread(mfile)

function saveit(hfile::String, key::String, data)
    h5open(hfile, "w") do file
        write(file, key, data, compress=3, chunk=size(data))
    end
end

# single precision suffices for floats:
converter(data::Array{<:AbstractFloat}) = Float32.(data)
converter(data::Array{<:Complex}) = ComplexF32.(data)
converter(data::Array{Bool}) = data

function dict2files(data)
    for key in keys(data)
        @show key
        hfile = key * ".h5"
        if isfile(hfile)
            @info "$hfile exists already"
            test = h5read(hfile, key)
            conv = converter(data[key])
            @assert isequal(test, conv) # allows NaN to match
        else
            saveit(hfile, key, converter(data[key]))
        end
    end
    nothing
end

if true
    for allkey in keys(alldata)
        @show allkey
        data = alldata[allkey]
        if data isa Dict
            dict2files(data)
        else
            dict2files(Dict(allkey => data))
        end
    end
end


#=
# version for file with a single variable
mfile = "Ankle_sensemap.mat"
hfile = replace(mfile, ".mat" => ".h5")

data = matread(mfile)

if isfile(hfile)
    @info "$hfile exists already"
else
    h5open(hfile, "w") do file
        for key in keys(data)
            tmp = data[key]
            tmp = ComplexF32.(tmp) # single precision suffices
            write(file, key, tmp, compress=3, chunk=size(tmp))
        end
    end
end

if true # test
    key = first(keys(data))
    test = h5read(hfile, "$key")
    @assert test == ComplexF32.(data[key])
end
=#
