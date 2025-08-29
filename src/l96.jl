"""
$SIGNATURES
Vectorized Lorenz96 equation without forcing term for one layer.
```math
\\frac{du}{dt} = (c*b) (u[j+1] - u[j-2]) * u[j-1] - c*u[j]
```
"""
lorenz69_layer(u; c::Number=1, b::Number=1) = (c*b) .* (circshift(u, -1) .- circshift(u, 2)) .* circshift(u, 1) .- c.* u

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
    (h, c, b, F) = p # forcing parameter
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

# some utilities for plotting 
layer1_θ(model::TwoLayerLorenz96) = range(start=0, step=2π/model.K, length=model.K)
layer2_θ(model::TwoLayerLorenz96) = range(start=0, step=2π/model.N_J, length=model.N_J)

"""
$SIGNATURES

Create an animated polar plot for the two-layer Lorenz96 system.

# Arguments
- `model::TwoLayerLorenz96`: The model structure
- `sol`: ODE solution object
- `filename::String="lorenz96_animation.mp4"`: Output filename for animation
- `framerate::Int=15`: Animation framerate
- `rlimits::Tuple=(-10,10)`: Radial limits for polar plot

# Returns
- `fig`: Makie figure object
"""
function animate_polar_plot(model::TwoLayerLorenz96, sol; 
                           filename::String="lorenz96_animation.mp4", 
                           framerate::Int=15, 
                           rlimits::Tuple=(-10,10))
    
    sol_data = Array(sol)
    
    # Create animated polar plot for both layers
    fig = Figure(size = (800, 800))
    ax = PolarAxis(fig[1,1], 
        rlimits=rlimits,
        radius_at_origin=rlimits[1],
        title = "Two-Layer Lorenz96 Dynamics"
    )

    # Time observable for animation
    itime = Observable(1)

    # Angular coordinates for both layers
    θ_slow = layer1_θ(model)
    θ_fast = layer2_θ(model)

    # Observables for data at current time
    slow_data = @lift sol_data[1:model.K, $itime]
    fast_data = @lift sol_data[(model.K+1):end, $itime]

    # Plot both layers with different colors and styles
    lines!(ax, θ_slow, slow_data, 
        color = :blue,
        linewidth = 3,
        label = "Slow Variables (X)")

    lines!(ax, θ_fast, fast_data, 
        color = :red, 
        linewidth = 1.5, 
        alpha = 0.7,
        label = "Fast Variables (Y)")

    # Add legend
    axislegend(ax, position = :rt)

    # Create animation
    record(fig, filename, 1:length(sol.t), framerate=framerate) do t
        itime[] = t
    end

    return fig
end

