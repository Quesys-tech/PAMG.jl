struct Aggregations
    head::Vector{Int}
    list::Vector{Int}
    tail::Vector{Int}
    function Aggregations(n)
        head = Vector{Int}(undef, 0)
        list = zeros(Int, n)
        tail = Vector{Int}(undef, 0)
        new(head, list, tail)
    end
end

struct Aggregation
    G::Aggregations
    i::Int #注意! 0 based index
    function Aggregation(G, i)
        @assert length(G.head) == length(G.tail)
        @assert 0 <= i
        if i + 1 > length(G.head) #新たなaggretation
            l_o = length(G.head)#追加前の長さ
            resize!(G.head, i + 1)
            resize!(G.tail, i + 1)
            G.head[l_o+1:i+1] .= 0 #0埋め
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
end

function Base.push!(G::Aggregations, i::Int, j::Int)
    push!(Aggregation(G, i), j)
end

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
                (!U[i] || j == i) && continue #ensure j∈ U\{i}
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
                (!U[i] || j == i) && continue #ensure j∈ U\{i}
                if (a[i, i] > 0 && a[i, j] < S_thresh[i]) || (a[i, i] < 0 && a[i, j] > S_thresh[i])
                    Sᵢ[j] = a[i, j]
                end
            end
            j = argmin(Sᵢ)
            push!(G, n_c, j)
        end

        for j in Aggregation(G, n_c)# U ← U \ G_n_c
            U[j] = false
        end
    end
    return G
end