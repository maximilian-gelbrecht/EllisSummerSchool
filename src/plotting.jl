
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
    slow_data = @lift [sol_data[1:model.K, $itime]; sol_data[1, $itime]] # add the first point again to close the circle 
    fast_data = @lift [sol_data[(model.K+1):end, $itime]; sol_data[(model.K+1), $itime]]

    # Plot both layers with different colors and styles
    lines!(ax, [θ_slow; θ_slow[1]], slow_data, 
        color = :blue,
        linewidth = 3,
        label = "Slow Variables (X)")

    lines!(ax, [θ_fast; θ_fast[1]], fast_data, 
        color = :red, 
        linewidth = 1.5, 
        alpha = 0.7,
        label = "Fast Variables (Y)")

    # Create animation
    record(fig, filename, 1:length(sol.t), framerate=framerate) do t
        itime[] = t
    end

    return fig
end

