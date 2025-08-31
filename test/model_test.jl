# test if the model initialize and run without blowing up 

@testset "TwoLayerLorenz96" begin
    model = TwoLayerLorenz96(K=36, J=10)
    p = default_parameters()
    u0 = default_initial_condition(model)
    tspan = (0., 5.)
    prob = ODEProblem(model, u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=0.005)
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test sol.u[end] != sol.u[1]         # something happended 
    @test all(isfinite.(sol.u[end]))   # no blow up for default config 
    @test (sum(abs.(sol.u[end])) .> 1) # no fixed point for the default config 
end

@testset "OneLayerLorenz96" begin
    model = OneLayerLorenz96(K=36)
    p = default_parameters()
    u0 = default_initial_condition(model)
    tspan = (0., 5.)
    prob = ODEProblem(model, u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=0.005)
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test sol.u[end] != sol.u[1]         # something happended 
    @test all(isfinite.(sol.u[end]))   # no blow up for default config 
    @test (sum(abs.(sol.u[end])) .> 1) # no fixed point for the default config 
end