module PAMG
using SparseArrays
using LinearAlgebra

# Write your package code here.
include("Aggregation.jl")
export Aggregation,Aggregations,push!,iterate, pairwise_aggregation
end
