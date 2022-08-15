@testset "AggregationStructTest" begin
    n = 8
    G = Aggregations(n)

    push!(G, 4, 8)
    push!(G, 0, 7)
    push!(G, 2, 6)
    push!(G, 1, 5)
    push!(G, 2, 4)
    push!(G, 2, 3)
    push!(G, 1, 2)
    push!(G, 1, 1)

    G_gt = Dict(4 => [8], 0 => [7], 2 => [6, 4, 3], 1 => [5, 2, 1], 3 => [])#ground truth
    for i = 0:4
        Gᵢ = Int[]
        for j in Aggregation(G, i)
            push!(Gᵢ, j)
        end
        @test eltype(Gᵢ) == Int
        @test isempty(symdiff(G_gt[i], Gᵢ))
    end
end

@testset "PairWiseAggregationTest" begin
    a = [-1 1 0 0 0; -1 2 -1 0 0; 0 1 -2 1 0; 0 0 -1 2 -1; 0 0 0 0 1.0]
    b = [0.25, 0, 0, 0, 1]

    G_gt = Dict{Int,Vector{Int}}(
        0 => [5],
        1 => [1, 2],
        2 => [3, 4]
    ) #ground truth

    G = pairwise_aggregation(a, 0.25, true)

    for i = 0:2
        Gᵢ = Int[]
        for j in Aggregation(G, i)
            push!(Gᵢ, j)
        end
        @test isempty(symdiff(G_gt[i], Gᵢ))
    end
end