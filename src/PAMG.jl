module PAMG
using SparseArrays
using LinearAlgebra

# Write your package code here.
include("Aggregation.jl")
export Aggregation, Aggregations,push!,iterate,length, pairwise_aggregation, prolongation_matrix, double_pairwise_aggregation
end
