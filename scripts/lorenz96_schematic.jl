using CairoMakie, Colors

"""
Create a schematic diagram of the two-layer Lorenz96 system structure.
Shows K slow variables in outer ring, each coupled to J fast variables in inner segments.
"""
function create_lorenz96_schematic(K=8, J=5)
    fig = Figure(size = (800, 800))
    ax = Axis(fig[1,1], 
        aspect = DataAspect(),
        title = "Two-Layer Lorenz96 System Structure\nK=$K slow variables, J=$J fast variables per slow variable",
        titlesize = 16
    )
    
    # Hide axes for clean schematic
    hidedecorations!(ax)
    hidespines!(ax)
    
    # Outer ring radius for slow variables
    R_slow = 1.0
    # Inner ring radius for fast variables  
    R_fast = 0.6
    
    # Colors
    slow_color = :blue
    fast_color = :red
    coupling_color = :gray
    
    # Draw slow variables (X_i) in outer ring
    for i in 1:K
        θ_slow = 2π * (i-1) / K
        x_slow = R_slow * cos(θ_slow)
        y_slow = R_slow * sin(θ_slow)
        
        # Draw slow variable node
        scatter!(ax, [x_slow], [y_slow], 
            color = slow_color, 
            markersize = 20,
            strokewidth = 2,
            strokecolor = :black)
        
        # Label slow variable
        text!(ax, x_slow * 1.15, y_slow * 1.15, 
            text = "X$i", 
            fontsize = 12,
            align = (:center, :center))
        
        # Draw J fast variables for this slow variable
        for j in 1:J
            # Angle for fast variable j associated with slow variable i
            θ_base = θ_slow - 2π/(K*J)
            θ_offset = 2π * (j-1) / (K*J)  # Spread fast variables around their slow variable
            θ_fast = θ_base + θ_offset - π/(K*J)
            
            x_fast = R_fast * cos(θ_fast)
            y_fast = R_fast * sin(θ_fast)
            
            # Draw fast variable node
            scatter!(ax, [x_fast], [y_fast], 
                color = fast_color, 
                markersize = 8,
                strokewidth = 1,
                strokecolor = :black)
            
            # Draw coupling line between X_i and Y_{i,j}
            lines!(ax, [x_slow, x_fast], [y_slow, y_fast], 
                color = coupling_color, 
                linewidth = 1,
                alpha = 0.6)
        end
    end
    
    # Draw fast variable interactions as a single circle (inner circle dynamics)
    circle_points = [Point2f(R_fast * 0.9 * cos(θ), R_fast * 0.9 * sin(θ)) for θ in 0:0.01:2π]
    lines!(ax, circle_points,
        color = fast_color,
        linewidth = 2,
        alpha = 0.6)
    
    # Draw connections between adjacent slow variables (L96 dynamics)
    for i in 1:(K+1)
        θ1 = 2π * (i-1) / K
        θ2 = 2π * (i % K) / K  # Next variable (with wraparound)
        
        x1, y1 = R_slow * cos(θ1), R_slow * sin(θ1)
        x2, y2 = R_slow * cos(θ2), R_slow * sin(θ2)
        
        # Curved arrow between adjacent slow variables
        arc!(ax, Point2f(0, 0), R_slow * 0.95, θ1, θ2, 
            color = slow_color, 
            linewidth = 2,
            alpha = 0.7)
    end
    
    # Add legend
    legend_elements = [
        MarkerElement(color = slow_color, marker = :circle, markersize = 15, 
                     strokewidth = 2, strokecolor = :black),
        MarkerElement(color = fast_color, marker = :circle, markersize = 8, 
                     strokewidth = 1, strokecolor = :black),
        LineElement(color = coupling_color, linewidth = 2, alpha = 0.6),
        LineElement(color = slow_color, linewidth = 2, alpha = 0.7),
        LineElement(color = fast_color, linewidth = 1.5, alpha = 0.5)
    ]
    
    legend_labels = [
        "Slow Variables (Xᵢ)",
        "Fast Variables (Yᵢⱼ)", 
        "Fast-Slow Coupling",
        "Slow Variable Interactions",
        "Fast Variable Interactions"
    ]
    
    Legend(fig[1,2], legend_elements, legend_labels, 
           framevisible = true,
           titletext = "System Components")
    
    # Add equations as text
    eq_text = """
    Slow: dXᵢ/dt = (Xᵢ₊₁ - Xᵢ₋₂)Xᵢ₋₁ - Xᵢ + F - (hc/b)∑ⱼYᵢⱼ
    Fast: dYᵢⱼ/dt = cb(Yᵢⱼ₊₁ - Yᵢⱼ₋₂)Yᵢⱼ₋₁ - cYᵢⱼ + (hc/b)Xᵢ
    """
    
    Label(fig[2,1:2], eq_text, 
          fontsize = 10,
          tellwidth = false,
          justification = :left)
    
    # Set axis limits
    xlims!(ax, -1.4, 1.4)
    ylims!(ax, -1.4, 1.4)
    
    return fig
end

# Create the schematic with default parameters
fig = create_lorenz96_schematic(8, 6)
save("lorenz96_schematic.png", fig)
fig
