using LinearAlgebra
using SparseArrays

"""
    rbf_phs_fd_all(xy1, xy_s, Nearest_Idx, k, m, d)

Generate global RBF-FD differentiation matrices using PHS + polynomials.

This function creates global differentiation matrices using Radial Basis Function
Finite Difference (RBF-FD) method with Polyharmonic Spline (PHS) basis functions
augmented with polynomial terms for enhanced accuracy.

INPUTS:
  xy1         - Target node coordinates [N1 x 2] where derivatives are evaluated
  xy_s        - Source node coordinates [N2 x 2] where function values are known
  Nearest_Idx - Index matrix [N1 x k] of k nearest neighbors for each target node
  k           - Stencil size (number of nearest neighbors used in local stencils)
  m           - Order of PHS-RBF (r^m, typically m=3 for 2D problems)
  d           - Degree of polynomial augmentation (-1=none, 0=constant, 1=linear, 2=quadratic)

OUTPUTS:
  Tuple of sparse differentiation matrices: (Dx, Dy, L, Dxx, Dyy, Dxy) where:
  Dx, Dy  = first derivatives (∂/∂x, ∂/∂y)
  L       = Laplacian (∂²/∂x² + ∂²/∂y²)
  Dxx,Dyy = second derivatives (∂²/∂x², ∂²/∂y²)
  Dxy     = mixed derivative (∂²/∂x∂y)
"""
function rbf_phs_fd_all(xy1::AbstractMatrix, xy_s::AbstractMatrix,
                        Nearest_Idx::AbstractMatrix, k::Int, m::Int, d::Int)
    
    X1 = view(xy1, :, 1); Y1 = view(xy1, :, 2)
    X2 = view(xy_s, :, 1); Y2 = view(xy_s, :, 2)
    N1 = length(X1); N2 = length(X2)

    # Preallocate weight matrices for each derivative operator
    weight_x  = zeros(N1, k)
    weight_y  = zeros(N1, k)
    weight_L  = zeros(N1, k)
    weight_xx = zeros(N1, k)
    weight_yy = zeros(N1, k)
    weight_xy = zeros(N1, k)

    # Main loop: Compute RBF-FD weights for each target node
    for i in 1:N1
        if i % 100 == 0
            println("Computing RBF-FD weights: node $i of $N1")
        end
        
        # Extract local stencil and evaluation point
        idx = @view Nearest_Idx[i, 1:k]
        xk = X2[idx]; yk = Y2[idx]
        xe = X1[i];   ye = Y1[i]

        # Coordinate transformation: translate and scale for numerical stability
        x = xk .- xe
        y = yk .- ye

        # Scale coordinates using distance to farthest stencil node
        scale = hypot(x[end], y[end])
        if scale < eps(Float64)
            scale = 1.0  # Avoid division by zero
        end
        scale_x = scale; scale_y = scale
        x ./= scale_x
        y ./= scale_y

        # Build PHS interpolation matrix A (k×k)
        # A[p,q] = ||(x_p,y_p)-(x_q,y_q)||^m
        dx = x .- x'
        dy = y .- y'
        A = (dx.^2 .+ dy.^2).^(0.5*m)  # r^m

        # Build right-hand side matrix for derivative operators
        r = sqrt.(x.^2 .+ y.^2)
        r_m2 = r .^ (m - 2)

        # RHS L0: columns [dx, dy, Lap, dxx, dyy, dxy]
        L0 = zeros(k, 6)
        
        # First derivatives: d/dx(r^m) = m*r^(m-2)*x, d/dy(r^m) = m*r^(m-2)*y
        L0[:, 1] .= -m .* r_m2 .* x
        L0[:, 2] .= -m .* r_m2 .* y
        
        # Laplacian: Δ(r^m) = m*(m+1)*r^(m-2) in 2D for PHS r^m
        L0[:, 3] .= m * (m + 1) .* r_m2
        
        # Second derivatives (handle r=0 safely)
        r2 = r.^2
        mask = r .> eps(Float64)
        L0[mask, 4] .= m .* r_m2[mask] .* (1 .+ (m-2) .* (x[mask].^2) ./ r2[mask])
        L0[mask, 5] .= m .* r_m2[mask] .* (1 .+ (m-2) .* (y[mask].^2) ./ r2[mask])
        L0[mask, 6] .= m .* (m-2) .* r_m2[mask] .* (x[mask] .* y[mask]) ./ r2[mask]

        # Polynomial augmentation for enhanced accuracy
        if d == -1  # Pure RBF case (no polynomial augmentation)
            A_aug = A
            L = L0
        else  # Augment RBF with polynomial terms up to degree d
            # Build polynomial basis: 1, x, y, x^2, xy, y^2, x^3, x^2*y, xy^2, y^3, ...
            np = (d+1)*(d+2) ÷ 2  # Total number of polynomial terms for degree d
            XY = zeros(k, np)
            col = 1
            
            # Assemble polynomial terms: x^i * y^j for i+j <= d
            for j in 0:d
                for p in 0:j
                    XY[:, col] .= (x.^(j-p)) .* (y.^p)
                    col += 1
                end
            end
            
            # Create polynomial matching conditions (derivatives of polynomial terms)
            L1 = zeros(np, 6)
            if d >= 1  # Linear terms: derivatives of x and y
                L1[2, 1] = 1  # d/dx(x) = 1
                L1[3, 2] = 1  # d/dy(y) = 1
            end
            if d >= 2  # Quadratic terms: derivatives of x^2, xy, y^2
                L1[4, 3] = 2; L1[4, 4] = 2  # Laplacian of x^2: d²/dx²(x²) = 2, d²/dy²(x²) = 0 -> sum = 2
                L1[6, 3] = 2; L1[6, 5] = 2  # Laplacian of y^2: d²/dx²(y²) = 0, d²/dy²(y²) = 2 -> sum = 2
                L1[5, 6] = 1               # d²/dxdy(xy) = 1
            end
            
            # Assemble augmented system: [RBF + polynomials; polynomial constraints]
            A_aug = [A XY; XY' zeros(col-1, col-1)]
            L = [L0; L1]
        end

        # Solve linear system for RBF-FD weights
        # Ensure matrix is symmetric (numerical symmetry)
        A_aug = (A_aug + A_aug') * 0.5
        
        # Solve augmented system: A_aug * weights = L
        # Add small regularization to handle near-singular cases
        regularization = 1e-12 * maximum(diag(A_aug)) * I
        A_aug_reg = A_aug + regularization
        
        W = try
            A_aug_reg \ L
        catch e
            if isa(e, SingularException)
                # If still singular, use pseudo-inverse
                pinv(A_aug) * L
            else
                rethrow(e)
            end
        end
        
        # Extract RBF weights (ignore polynomial constraint multipliers)
        w = @view W[1:k, :]
        
        # Scale weights back to original coordinate system
        weight_x[i, 1:k]  .= w[:, 1] ./ scale_x
        weight_y[i, 1:k]  .= w[:, 2] ./ scale_y
        weight_L[i, 1:k]  .= w[:, 3] ./ (scale^2)
        weight_xx[i, 1:k] .= w[:, 4] ./ (scale_x^2)
        weight_yy[i, 1:k] .= w[:, 5] ./ (scale_y^2)
        weight_xy[i, 1:k] .= w[:, 6] ./ (scale_x * scale_y)
    end

    # Assemble global sparse differentiation matrices
    # Create index arrays for sparse matrix construction
    rows = repeat(collect(1:N1), inner=k)
    cols = reshape(permutedims(Nearest_Idx[:, 1:k]), :)

    # Flatten weight matrices into vectors and construct sparse matrices
    Dx  = sparse(rows, cols, reshape(permutedims(weight_x), :),  N1, N2)
    Dy  = sparse(rows, cols, reshape(permutedims(weight_y), :),  N1, N2)
    L   = sparse(rows, cols, reshape(permutedims(weight_L), :),  N1, N2)
    Dxx = sparse(rows, cols, reshape(permutedims(weight_xx), :), N1, N2)
    Dyy = sparse(rows, cols, reshape(permutedims(weight_yy), :), N1, N2)
    Dxy = sparse(rows, cols, reshape(permutedims(weight_xy), :), N1, N2)
    return Dx, Dy, L, Dxx, Dyy, Dxy
end
