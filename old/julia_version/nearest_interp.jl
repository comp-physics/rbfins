using NearestNeighbors

"""
    nearest_interp(target::AbstractMatrix, source::AbstractMatrix, k::Int)

Find k nearest neighbors for each target point in source points.

INPUTS:
  target - Target node coordinates [N1 x 2] where N1 is number of target points
  source - Source node coordinates [N2 x 2] where N2 is number of source points  
  k      - Number of nearest neighbors to find

OUTPUTS:
  Nearest_Idx - Index matrix [N1 x k] containing indices into source for k nearest neighbors
"""
function nearest_interp(target::AbstractMatrix, source::AbstractMatrix, k::Int)
    @assert size(target, 2) == 2 "Target points must be N×2 matrix"
    @assert size(source, 2) == 2 "Source points must be N×2 matrix"
    @assert k > 0 "Number of neighbors k must be positive"
    
    # Build KD-tree from source points (transpose for NearestNeighbors.jl convention)
    tree = KDTree(permutedims(source))  # 2×N2
    
    # Find k nearest neighbors for each target point
    idxs, _ = knn(tree, permutedims(target), k)
    
    # Convert to N1×k matrix format
    N1 = size(target, 1)
    Nearest_Idx = Matrix{Int}(undef, N1, k)
    
    @inbounds for i in 1:N1
        @inbounds for j in 1:k
            if j <= length(idxs[i])
                Nearest_Idx[i, j] = idxs[i][j]
            else
                # Handle case where fewer than k neighbors exist
                Nearest_Idx[i, j] = idxs[i][end]
            end
        end
    end
    
    return Nearest_Idx
end
