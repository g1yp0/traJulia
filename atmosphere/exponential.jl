module AtmosphereExponential

"""
This module provides a low-fidelity model of Earth's atmopsheric density. 

It utilises a simple exponential atmosphere model from:
    Regan et al., "Dynamics of Atmospheric Re-Entry"

# Arguments
- `altitude`: Altitude above sea level in metres (m).

# Returns
- Atmospheric density in kilograms per cubic metre (kg/m³).
"""
function atmosphere(altitude)
    ρ0 = 1.225          # Sea level standard density (kg/m³)
    h_scale = 8434      # Scale height (m) from Regan
    ρ = ρ0 * exp(-altitude / h_scale)
    return ρ
end

export atmosphere

end
