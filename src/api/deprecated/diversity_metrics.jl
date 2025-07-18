using NearestNeighbors
using Distances
"""
    get_max_pairwise_distance(population::AbstractMatrix{Float64}, method=Euclidean())

Computes the maximum pairwise distance as a diversity metric for a population of individuals. Non-allocating.
"""
function get_max_pairwise_distance(population::AbstractMatrix{Float64}, method=Euclidean())
    # Delegate to the unified folding implementation defined below.  This keeps
    # the original API intact while removing code duplication.
    _fold_pairwise(population, 0.0, (acc, d) -> max(acc, d); method = method)
end

function get_max_pairwise_distance(population::Vector{Vector{Float64}}, method=Euclidean())
    isempty(population) && return NaN
    pop_matrix = stack(population)
    get_max_pairwise_distance(pop_matrix, method)
end

function get_max_pairwise_distance(df::AbstractDataFrame, exclude_cols::AbstractVector{Symbol}, method=Euclidean())
    dfmat = df_to_matrix(df, exclude_cols)
    return get_max_pairwise_distance(dfmat, method)
end

"""
    get_avg_pairwise_distance(population::AbstractMatrix{Float64}, method=Euclidean())

Computes the average pairwise distance as a diversity metric for a population of individuals.
Each column represents an individual. Non-allocating.
"""
function get_avg_pairwise_distance(population::AbstractMatrix{Float64}, method=Euclidean())
    n = size(population, 2)
    n < 2 && return 0.0
    
    total_dist = 0.0
    for i in 1:n
        for j in i+1:n
            # Calculate distance between the i-th and j-th individual using views
            dist = evaluate(method, view(population, :, i), view(population, :, j))
            total_dist += dist
        end
    end

    # Compute the average distance. The number of pairs is n*(n-1)/2
    return total_dist / (n * (n - 1) / 2)
end

# Vector method that converts to matrix format
function get_avg_pairwise_distance(population::Vector{Vector{Float64}}, method=Euclidean())
    isempty(population) && return 0.0
    pop_matrix = stack(population)
    return get_avg_pairwise_distance(pop_matrix, method)
end

# DataFrame method
function get_avg_pairwise_distance(df::AbstractDataFrame, exclude_cols::AbstractVector{Symbol}, method=Euclidean())
    dfmat = df_to_matrix(df, exclude_cols)
    return get_avg_pairwise_distance(dfmat, method)
end

"""
    get_nearest_neighbor_distance(S::AbstractMatrix{Float64})

Uses `KDTree` to get distance to nearest neighbor for each individual in a population of individuals.
"""
function get_nearest_neighbor_distance(S::AbstractMatrix)
    #- Create a KDTree with Euclidean metric
    kdtree = KDTree(S)

    k = 2 #* Find the first nearest neighbor for each point

    #- Find the nearest neighbor for each point
    _, dists = knn(kdtree, S, k, true)

    #- Retreive minimum distance for each individual
    getindex.(dists, 2)
end


function get_spread(S::AbstractMatrix) 
    n = size(S,2) #* number of individuals
    n == 1 && return NaN #* if there is only one individual, return NaN

    #- Compute minimum distance for each individual
    Δₖ = get_nearest_neighbor_distance(S)

    #- Compute mean minimum distance
    Δ = mean(Δₖ)

    #- Compute spread metric by summing the absolute difference between each minimum distance and the mean, scaled by the mean
    sum(abs.(Δₖ.-Δ))/n*Δ
end

"""
    get_spread(df::DataFrame, exclude_cols::Vector{Symbol})

Computes the spread metric for a population of individuals, which is the average deviation of all minimum distances to the nearest neighbor from the mean distance, scaled by the mean.
"""
function get_spread(df::AbstractDataFrame, exclude_cols::AbstractVector{Symbol})
    dfmat = df_to_matrix(df, exclude_cols)
    return get_spread(dfmat)
end


using Clustering

"""
    calculate_shannon_index(clustering_result::ClusteringResult)

Calculate the Shannon diversity index for a given clustering result.

# Arguments
- `clustering_result::ClusteringResult`: The result of a clustering operation.

# Returns
- `Float64`: The Shannon diversity index.
"""
function get_shannon_index(cr::ClusteringResult)
    # n = nclusters(cr) # get the number of clusters
    weighted_cluster_sizes = wcounts(cr) # get the weighted cluster sizes
    p = weighted_cluster_sizes / sum(weighted_cluster_sizes) # get the proportions of each cluster
    H = -sum(p .* log.(p)) # compute the Shannon Index
    return H |> abs
end

function get_shannon_index(df::AbstractDataFrame, exclude_cols::AbstractVector{Symbol})
    best_k = get_optimal_clusters(df, 20, exclude_cols)
    cr = get_kmeans(df, best_k, exclude_cols)
    return shannon_index(cr)
end


"""
    get_simpson_index(clustering_result::ClusteringResult)

Calculate the Simpson diversity index for a given clustering result.

# Arguments
- `cr::ClusteringResult`: The result of a clustering operation.

# Returns
- `Float64`: The Simpson diversity index.
"""
function get_simpson_index(cr::ClusteringResult)
    cluster_sizes = wcounts(cr) # get the cluster sizes
    total_count = sum(cluster_sizes) # total number of elements
    proportions = cluster_sizes / total_count # compute the proportions of each cluster
    simpson_index = 1.0 - sum(proportions.^2) # compute the Simpson Index
    return simpson_index
end


#< GPU implementations >#
using CUDA

"""
GPU-accelerated implementation of max and avg pairwise distances.
Returns (max_dist, avg_dist) tuple.
"""
function compute_pairwise_distances_gpu(X::AbstractMatrix{Float32}; chunk_size=1000)
    n = size(X, 2)
    n_params = size(X, 1)

    max_chunk_size = min(chunk_size, 
                        floor(Int, 3.5e9 / (4 * n_params)))

    max_dist = 0.0f0
    total_dist = 0.0f0
    n_pairs = 0
    
    for i in 1:max_chunk_size:n
        chunk_end_i = min(i + max_chunk_size - 1, n)
        chunk1 = cu(X[:, i:chunk_end_i])
        
        norms = sum(chunk1 .* chunk1, dims=1)
        dots = CUDA.CUBLAS.gemm('T', 'N', chunk1, chunk1)
        dists = sqrt.(max.(reshape(norms, 1, :) .+ reshape(norms, :, 1) .- 2.0f0 .* dots, 0.0f0))
        
        current_max = maximum(dists)
        current_sum = sum(triu(dists, 1))
        current_pairs = (size(chunk1, 2) * (size(chunk1, 2) - 1)) ÷ 2
        
        max_dist = max(max_dist, Float32(current_max))
        total_dist += Float32(current_sum)
        n_pairs += current_pairs
        
        CUDA.reclaim()
    end
    
    return (Float64(max_dist), Float64(total_dist / n_pairs))
end

# GPU methods for Float32 inputs only
function get_max_pairwise_distance(population::AbstractMatrix{Float32}, method=Euclidean())
    compute_pairwise_distances_gpu(population)[1]
end

function get_avg_pairwise_distance(population::AbstractMatrix{Float32}, method=Euclidean())
    compute_pairwise_distances_gpu(population)[2]
end

######################################################################
# Pairwise distance utilities (CPU path)
######################################################################

# Internal helper: generic fold over all unique pairwise distances.
# `acc` – initial accumulator,
# `update` – function (acc, d) -> newacc that consumes each computed distance.
# Returns the final accumulator.  The Euclidean() metric is assumed unless
# another is passed.
function _fold_pairwise(pop::AbstractMatrix{<:Real}, acc, update!; method = Euclidean())
    n = size(pop, 2)
    n < 2 && return acc
    for i in 1:n - 1
        vi = view(pop, :, i)
        for j in i + 1:n
            d = evaluate(method, vi, view(pop, :, j))
            acc = update!(acc, d)
        end
    end
    return acc
end

# Convenience wrappers -------------------------------------------------------

"""
    get_max_pairwise_distance(pop::AbstractMatrix; method = Euclidean())

Compute the maximum pairwise Euclidean distance (CPU).  Falls back to GPU
method if the input is a `Float32` matrix *and* CUDA is available.
"""
function get_max_pairwise_distance(pop::AbstractMatrix{Float64}; method = Euclidean())
    _fold_pairwise(pop, 0.0, (acc, d) -> max(acc, d); method = method)
end

# Keep the existing CUDA overrides for Float32 matrices untouched (defined
# further below), so we *rename* the old method to avoid ambiguity.

"""
    _get_max_pairwise_distance_old

Legacy implementation retained for benchmarking/back-compat purposes.
"""
const _get_max_pairwise_distance_old = get_max_pairwise_distance

# Median pairwise distance ---------------------------------------------------

"""
    get_median_pairwise_distance(pop::AbstractMatrix; method = Euclidean())

Compute the exact median of all pairwise Euclidean distances.  This is an
O(n²) operation and allocates a `Vector` of length n*(n-1)/2, but guarantees
deterministic results regardless of population size.
"""
function get_median_pairwise_distance(pop::AbstractMatrix; method = Euclidean())
    n = size(pop, 2)
    n < 2 && return 0.0

    n_pairs = (n * (n - 1)) >>> 1  # every unique (i,j) with j>i
    dvec = Vector{eltype(pop)}(undef, n_pairs)
    idx_ref = Ref(1)

    # deterministically collect all pairwise distances
    _fold_pairwise(pop, nothing, (nothing, d) -> begin dvec[idx_ref[]] = d; idx_ref[] += 1; nothing; end; method = method)

    return Statistics.median(dvec)
end

# Wrapper overloads -----------------------------------------------------------

function get_median_pairwise_distance(pop::Vector{Vector{Float64}})
    isempty(pop) && return 0.0
    pop_matrix = stack(pop)
    return get_median_pairwise_distance(pop_matrix)
end

function get_median_pairwise_distance(df::AbstractDataFrame, exclude_cols::AbstractVector{Symbol})
    dfmat = df_to_matrix(df, exclude_cols)
    return get_median_pairwise_distance(dfmat)
end

##############################################################################
