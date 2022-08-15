# PAMG.jl

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Quesys-tech.github.io/PAMG.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Quesys-tech.github.io/PAMG.jl/dev/)
-->
[![Build Status](https://github.com/Quesys-tech/PAMG.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Quesys-tech/PAMG.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Quesys-tech/PAMG.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Quesys-tech/PAMG.jl)


**PAMG** (Plain Aggregation algebraic MultiGrid solver) is a implementation of 
PA-AMG presented by Notay[[1]](#references).

## About PA-AMG
Multigrid methods are algorithms to solve large systems of linear equations faster than conventional iterative algorithms (Gauss-Seidel, Conjugate Gradient, etc.) utilizing coarse and fine grids, which can be classified as geometric multigrid and algebraic multigrid (AMG). Geometric multigrid uses geometric information from discretized problems to set up coarse grids. Geometric multigrid is available only for structured gird. 
However, An unstructured grid, in which we cannot utilize geometric multigrid methods, is widely used for many numerical methods to express complicated geometries. Algebraic multigrid methods, which set up coarse grids from the information of coefficient matrices, are suitable for not only unstructured grids but also structured grids.
PA-AMG uses pairwise matching to set up the aggregation of unknowns to construct coarse grids. PA-AMG ensures relatively lower setup computational cost and memory usage.

## References
1. Notay, Yvan. "An aggregation-based algebraic multigrid method." Electronic transactions on numerical analysis 37.6 (2010): 123-146. [pdf](http://etna.mcs.kent.edu/vol.37.2010/pp123-146.dir/pp123-146.pdf)