module EllisSummerSchool

using DocStringExtensions, Makie 

export OneLayerLorenz96, TwoLayerLorenz96, lorenz69_layer, default_parameters
export subgrid_save_func, layer1_θ, layer2_θ, animate_polar_plot

include("l96.jl")

end
