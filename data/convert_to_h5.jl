# convert_to_h5.jl
# convert data from proprietary .mat format to open .h5 format

using MAT
using HDF5

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
