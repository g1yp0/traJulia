module AtmosphereExponential

using Unitful
using Unitful.DefaultSymbols  # For unit symbols like m, kg, etc.


"""
This module provides a low-fidelity model of Earth's atmopsheric density. 

It utilises a simple exponential atmosphere model from:
- Regan, F. J. (1993). 'Dynamics of Atmospheric Re-Entry'.

# Arguments
- `altitude`: Altitude above sea level in metres (m), with units (e.g., `1000u"m"`).

Note, this altitude is geometric, not geopotential.

# Returns
- Atmospheric density (ρ), with units (e.g., 1.225u"kg/m^3)
"""

function atmosphere(altitude::Quantity{<:Real, Unitful.𝐋})
    ρ0 = 1.225u"kg/m^3"          # Sea level standard density (kg/m³)
    h_scale = 6700u"m"           # Scale height from Regan
    ρ = ρ0 * exp(-altitude / h_scale)
    return ρ
end

export atmosphere

end
