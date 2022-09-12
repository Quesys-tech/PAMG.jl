struct Aggregations
    head::Vector{Int}
    list::Vector{Int}
    tail::Vector{Int}
    len::Vector{Int} #長さ
    function Aggregations(n)
        head = Vector{Int}(undef, 0)
        list = Vector{Int}(undef, n)
        tail = Vector{Int}(undef, 0)
        len = Vector{Int}(undef, 0)
        list .= 0
        new(head, list, tail, len)
    end
end

struct Aggregation
    G::Aggregations
    i::Int #注意! 0 based index
    function Aggregation(G::Aggregations, i::Int)
        @assert length(G.head) == length(G.tail)
        @assert 0 <= i
        if i + 1 > length(G.head) #新たなaggretation
            l_o = length(G.head)#追加前の長さ
            resize!(G.head, i + 1)
            resize!(G.tail, i + 1)
            resize!(G.len, i + 1)
            G.head[l_o+1:i+1] .= 0 #0埋め
            G.len[l_o+1:i+1] .= 0 #0埋め
        end
        new(G, i)
    end
end

function Base.iterate(Gᵢ::Aggregation)
    item = Gᵢ.G.head[Gᵢ.i+1]
    if 1 <= item <= length(Gᵢ.G.list)
        return (item, item)
    else
        return nothing
    end
end

function Base.iterate(Gᵢ::Aggregation, state::Int)
    item = Gᵢ.G.list[state]
    if 1 <= item <= length(Gᵢ.G.list)
        return (item, item)
    else
        return nothing
    end
end

Base.eltype(::Type{Aggregation}) = Int

function Base.push!(G::Aggregation, i::Int)
    @assert 1 <= i <= length(G.G.list) #iは有効な範囲
    if 1 <= G.G.head[G.i+1] <= length(G.G.list) #Gᵢ は空集合でない   
        G.G.list[G.G.tail[G.i+1]] = i #末尾に追加 
        G.G.tail[G.i+1] = i
    else
        G.G.head[G.i+1] = i
        G.G.list[i] = 0
        G.G.tail[G.i+1] = i
    end
    G.G.len[G.i+1] += 1
end

function Base.push!(G::Aggregations, i::Int, j::Int)
    push!(Aggregation(G, i), j)
end

function Base.length(G::Aggregation)
    return G.G.len[G.i+1]
end


"""
    pairwise_aggregation(a::AbstractMatrix{T}, β::T, finest::Bool=false) where {T}

execute pair wise aggregation
"""
function pairwise_aggregation(a::AbstractMatrix{T}, β::T, finest::Bool=false) where {T}
    n = size(a)[1]
    @assert size(a) == (n, n)
    @assert 0 <= β <= 1

    G = Aggregations(n)
    if finest == true
        for i = 1:n
            #|aᵢᵢ| > 5∑_{i≂̸j}|aᵢⱼ| then i∈G₀
            sum_aᵢⱼ = sum(abs, @view a[i, :]) - abs(a[i, i])
            abs(a[i, i]) > 5sum_aᵢⱼ && push!(G, 0, i)
        end
    end

    #U ← {1...n} \ G₀
    U = trues(n)
    for i in Aggregation(G, 0)
        U[i] = false
    end
    n_c = 0
    num_S = Vector{Int}(undef, n) #number of elements in Sᵢ

    # calculate coupling threshold 
    # aᵢⱼ < β min_{aᵢₖ < 0} aᵢₖ (aᵢᵢ > 0)
    # aᵢⱼ > β max_{aᵢₖ > 0} aᵢₖ (aᵢᵢ < 0)
    S_thresh = zeros(n)
    for k = 1:n, i = 1:n
        if a[i, i] > 0
            if a[i, k] < 0
                S_thresh[i] = min(a[i, k], S_thresh[i])
            end
        else
            if a[i, k] > 0
                S_thresh[i] = max(a[i, k], S_thresh[i])
            end
        end
    end
    S_thresh .*= β

    while maximum(U) #while U≂̸∅
        n_c += 1

        num_S .= n + 1
        for i in 1:n
            !U[i] && continue #ensure i ∈ U
            num_S[i] = 0
            for j = 1:n
                (!U[j] || j == i) && continue #ensure j∈ U\{i}
                if (a[i, i] > 0 && a[i, j] < S_thresh[i]) || (a[i, i] < 0 && a[i, j] > S_thresh[i])
                    num_S[i] += 1
                end
            end
        end
        i = argmin(num_S)
        push!(G, n_c, i)

        if num_S[i] != 0
            Sᵢ = Dict{Int,T}()
            for j = 1:n
                (!U[j] || j == i) && continue #ensure j∈ U\{i}
                if (a[i, i] > 0 && a[i, j] < S_thresh[i]) || (a[i, i] < 0 && a[i, j] > S_thresh[i])
                    Sᵢ[j] = a[i, j]
                end
            end
            if a[i, i] > 0
                j = argmin(Sᵢ)
            else
                j = argmax(Sᵢ)
            end
            push!(G, n_c, j)
        end

        for j in Aggregation(G, n_c)# U ← U \ G_n_c
            U[j] = false
        end
    end
    return G, n_c
end

"""
    prolongation_matrix(G::Aggregations, T=Float64)::SparseMatrixCSC{T,Int}

generate prolongation matrix P
"""
function prolongation_matrix(G::Aggregations, T=Float64)::SparseMatrixCSC{T,Int}
    n = length(G.list)
    nₘ = length(G.head) - 1
    @assert 0 < nₘ <= n

    is = Int[]
    js = Int[]
    vs = T[]
    sizehint!(is, n)
    sizehint!(js, n)
    sizehint!(vs, n)
    for j = 1:nₘ
        for i in Aggregation(G, j)
            push!(is, i)
            push!(js, j)
            push!(vs, one(T))
        end
    end
    sparse(is, js, vs, n, nₘ)
end

"""
    double_pairwise_aggregation(a::AbstractMatrix{T}, β::T, finest::Bool=false) where {T}

double pairwise aggregation 
"""
function double_pairwise_aggregation(a::AbstractMatrix{T}, β::T, finest::Bool=false) where {T}
    n = size(a)[1]
    @assert size(a) == (n, n)
    @assert 0 <= β <= 1
    G_ast, n_m = pairwise_aggregation(a, β, finest)
    P_ast = prolongation_matrix(G_ast, T)
    a_m = P_ast' * a * P_ast
    G_ast2, n_c = pairwise_aggregation(a_m, β)

    G_vec = Vector{Int}(undef, n)
    for i = 1:n_c
        for j in Aggregation(G_ast2, i)
            for k in Aggregation(G_ast, j)
                G_vec[k] = i
            end
        end
    end
    for i in Aggregation(G_ast, 0)
        G_vec[i] = 0
    end

    G = Aggregations(n)
    for j = 1:n
        push!(G, G_vec[j], j)
    end
    return G, n_c
end