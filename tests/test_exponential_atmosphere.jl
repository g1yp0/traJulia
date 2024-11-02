using Test
using Unitful
using traJulia.ExponentialAtmosphere  # Adjusted to match your module name

@testset "Exponential Atmosphere Tests" begin
    # Test case for sea level (altitude = 0 m)
    altitude = 0.0u"m"
    density = atmosphere(altitude)
    
    @test density ≈ 1.225u"kg/m^3" atol=0.001u"kg/m^3"  # Use tolerance for floating-point comparison
    
    # Test case for an altitude of 6700 m (the scale height)
    altitude = 6700.0u"m"
    density = atmosphere(altitude)
    
    @test density ≈ 1.225u"kg/m^3" * exp(-1.0) atol=0.001u"kg/m^3"  # Expected density at scale height
end
