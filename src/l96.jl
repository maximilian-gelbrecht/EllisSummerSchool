"""
$SIGNATURES

Basic one layer Lorenz96 model.

```math
\\frac{du_{j}}{dt} = (u_{j+1} - u_{j-2}) * u_{j-1} - u_{j} + F
```
"""
struct OneLayerLorenz96 
end 

function (model::OneLayerLorenz96)(du, u, p, t)
    F = p[1]  # forcing parameter
        
    # Vectorized Lorenz96 equation: du/dt = (u[j+1] - u[j-2]) * u[j-1] - u[j] + F
    du .= (circshift(u, -1) .- circshift(u, 2)) .* circshift(u, 1) .- u .+ F
end

"""
$SIGNATURES

Basic two layer Lorenz96 model with fast-slow coupling, K slow variables and J fast variables per slow variable.

```math
\\frac{dX_{i}}{dt} = (X_{i+1} - X_{i-2}) * X_{i-1} - X_{i} + F - \\frac{hc}{b} \\sum_{j=1}^J Y_{i,j}
\\frac{dY_{i,j}}{dt} = cb\\cdot (Y_{i,j+1} - Y_{i,j-2}) * Y_{i,j+1} - cY_{i,j} + (hc/b) X_{i}
```
"""
@kwdef struct TwoLayerLorenz96 
    K::Integer # 1st layer: number of grid points of the slow variable
    J::Integer # 2nd layer: number of grid points per fast variable
    N_J::Integer=K*J # 2nd layer: number of total fast variables = K*J
    N::Integer=K+K*J # total number of variables
end 

function (model::TwoLayerLorenz96)(du, u, p, t)
    (; K, J, N) = model
    (h, c, b, F) = p 
    hcb = h*c/b

    X = view(u, 1:K)
    Y = view(u, (K+1):N)
    
    # first layer (slow) 
    du[1:K] .= (circshift(X, -1) .- circshift(X, 2)) .* circshift(X, 1) .- X .+ F .- (hcb .* vec(sum(reshape(Y, J, K), dims=1)))

    # second layer (fast)
    du[K+1:end] .= (c*b) * (circshift(Y, 1) .- circshift(Y, -2)) .* circshift(Y, -1) .- c*Y .+ (hcb .* vec(repeat(reshape(X, 1, :), J)))

    return nothing
end

default_parameters() = (h=1., c=10., b=10., F=12.)
default_initial_condition(model::TwoLayerLorenz96) = [0.5 * sin.(2π * 3*(0:model.K-1) / model.K); 0.01 * randn(model.N_J)]

"""
Subgrid saving function for use with `SavingCallback`.
"""
function subgrid_save_func(u, t, integrator) 
    (h, c, b, F) = integrator.p
    (; K, J, N) = SciMLBase.unwrapped_f(integrator.f.f)
    hcb = h*c/b
    Y = view(u, (K+1):N)
    return (hcb .* vec(sum(reshape(Y, J, K), dims=1)))
end 
