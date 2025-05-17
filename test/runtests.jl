using Bensolve
using Test

@testset "Bensolve.jl" begin
    #for ex in filter(f -> endswith(f, ".vlp"), readdir("examples/", join=true))
    for i in 1:6
        0 == solve("examples/ex0$(i).vlp")
    end
end
