using MIRTjim: jim
using Plots: savefig

using HDF5: h5read

jim(:prompt, true)

files = readdir()
ish5 = f -> f[end-2:end] == ".h5"
files = filter(ish5, files)
for file in files
    pngfile = replace(file, ".h5" => ".png")
    @show file, pngfile
    key = replace(file, ".h5" => "")
    data = h5read(file, key)
    @show typeof(data)
    jim(data, title = "$key $(typeof(data))")
    if isfile(pngfile)
        @info("$pngfile exists")
    else
        savefig(pngfile)
    end
end
