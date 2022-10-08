module ParaPDEs

export CBEM, CBIM, CBCNM, NLBIM, NLBCNM

using LinearAlgebra

"""
    CBEM(x, t, c1, c2, f)
Constant Boundary Explicit Method for Parabolic PDEs.\\
Computes the numerical solution over the domain ranges of `x` and `t`, with 
`c1` and `c2` as the left and right constant boundary values. `f` is the function
defining the initial state at t=0.
Returns a length(x) x length(t) Matrix containing the computed values
"""
function CBEM( x, t, b1, b2, f, k)
    nx = length(x)
    nt = length(t)
    r = k * Float64(t.step) / Float64(x.step)^2

    M = SymTridiagonal(fill(1-2r, nx-2),
                       fill(r, nx-3))

    U = [reshape(b1.(t),1,nt);
         [f.(x[2:end-1]) zeros(nx-2, nt-1)];
         reshape(b2.(t),1,nt)]

    for step ∈ 2:length(t)
        b = [r*U[1,step-1]; zeros(nx-4); r*U[end, step-1]]
        U[2:end-1, step] = M * U[2:end-1, step-1] + b
    end

    return U
end

"""
    CBIM(x, t, c1, c2, f)
Constant Boundary Implicit Method for Parabolic PDEs.\\
Computes the numerical solution over the domain ranges of `x` and `t`, with 
`c1` and `c2` as the left and right constant boundary values. `f` is the function
defining the initial state at t=0.
Returns a length(x) x length(t) Matrix containing the computed values
"""
function CBIM( x, t, b1, b2, f, k)
    nx = length(x)
    nt = length(t)
    r = k * Float64(t.step) / Float64(x.step)^2

    M = SymTridiagonal(fill(1+2r, nx-2),
                       fill(-r, nx-3))

    U = [reshape(b1.(t),1,nt);
         [f.(x[2:end-1]) zeros(nx-2, nt-1)];
         reshape(b2.(t),1,nt);]

    for step ∈ 2:length(t)
        c = [r*U[1,step]; zeros(nx-4); r*U[end, step]]
        U[2:end-1, step] = M \ (U[2:end-1, step-1] + c)
    end

    return U
end

function NLBIM(x, t, f1, f2, i, k)
    nx = length(x)
    nt = length(t)
    δx = Float64(x.step)
    δt = Float64(t.step)
    r = k * δt / δx^2

    M = Tridiagonal([fill(-r, nx-2); -2r],
                    [1+2r+2r*δx; fill(1+2r, nx-2); 1+2r+2r*δx],
                    [-2r; fill(-r, nx-2)])

    U = [i.(x) zeros(nx, nt-1)]

    for step ∈ 2:nt
        b = [2r*δx*f1.(t[step]); zeros(nx-2); -2r*δx*f2.(t[step])]
        U[:,step] = M \ (U[:,step-1] - b)
    end

    return U
end

"""
    CBCNM(x, t, c1, c2, f)
Constant Boundary Crank-Nicholson Method for Parabolic PDEs.\\
Computes the numerical solution over the domain ranges of `x` and `t`, with 
`c1` and `c2` as the left and right constant boundary values. `f` is the function
defining the initial state at t=0.
Returns a length(x) x length(t) Matrix containing the computed values
"""
function CBCNM(x, t, b1, b2, f, k)
    nx = length(x)
    nt = length(t)
    r = k * Float64(t.step) / Float64(x.step)^2
    
    # Explicit Array
    Mx = SymTridiagonal(fill(1-r, nx-2),
                        fill(r/2, nx-3))

    #Implicit Array
    Mi = SymTridiagonal(fill(1+r, nx-2),
                        fill(-r/2, nx-3))

    U = [reshape(b1.(t),1,nt);
         [f.(x[2:end-1]) zeros(nx-2, nt-1)];
         reshape(b2.(t),1,nt);]

    for step ∈ 2:length(t)
        b = [(r/2)*U[1,step-1]; zeros(nx-4); (r/2)*U[end, step-1]]
        c = [(r/2)*U[1,step]; zeros(nx-4); (r/2)*U[end, step]]
        U[2:end-1, step] = Mi \ ((Mx * U[2:end-1, step-1] + b) + c)
    end

    return U
end

function NLBCNM(x, t, f1, f2, i, k)
    nx = length(x)
    nt = length(t)
    δx = Float64(x.step)
    δt = Float64(t.step)
    r = k * δt / δx^2

    Mx = Tridiagonal([fill(r, nx-2); 2r],
                     [2-2r-2r*δx; fill(2-2r, nx-2); 2-2r-2r*δx],
                     [2r; fill(r,nx-2)])

    Mi = Tridiagonal([fill(-r, nx-2); -2r],
                     [2+2r+2r*δx; fill(2+2r, nx-2); 2+2r+2r*δx],
                     [-2r; fill(-r, nx-2)])

    U = [i.(x) zeros(nx, nt-1)]

    for step ∈ 2:nt
        b = [-2r*δx*f1.(t[step-1]); zeros(nx-2); 2r*δx*f2.(t[step-1])]
        c = [2r*δx*f1.([step]); zeros(nx-2); -2r*δx*f2.(t[step])]
        U[:,step] = Mi \ ((Mx * U[:,step-1] + b) - c)
    end

    return U
end

end