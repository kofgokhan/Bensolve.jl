using Bensolve
using Test

function _txt_file_to_matrix_problem_data(filename)
    
end

@testset "Bensolve.jl using *.vlp files" begin
    #for ex in filter(f -> endswith(f, ".vlp"), readdir("examples/", join=true))
    for i in 1:6
        0 == solve("examples/ex0$(i).vlp")
    end
end

@testset "Bensolve.jl using matrices and vector for problem data" begin
    #for ex in filter(f -> endswith(f, ".vlp"), readdir("examples/", join=true))
    for i in 1:6
        data = _txt_file_to_matrix_problem_data("examples/ex0$(i).vlp")
        0 == solve(data...)
    end
end