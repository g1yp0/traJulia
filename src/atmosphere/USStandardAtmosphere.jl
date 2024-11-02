module AtmosphereReganStandard

"""
This module provides atmospheric properties based on the 1976 US Standard Atmosphere (US76)
up to 86 km and the 1962 US Standard Atmosphere (US62) above 86 km.

The model is described in:
Regan, F. J., & Anandakrishnan, S. M. (1993). "Dynamics of Atmospheric Re-Entry".

### Arguments:

- `altitude`: Altitude at which to calculate atmospheric properties, with units (e.g., `1000u"m"`).
- `T₀`: Ground level temperature, with units (default: `288.15u"K"`).
- `P₀`: Ground level pressure, with units (default: `101325.0u"Pa"`).

### Returns:

A named tuple containing:

- `pressure`: Atmospheric pressure at the specified altitude.
- `density`: Atmospheric density at the specified altitude.
- `temperature`: Kinetic temperature at the specified altitude.
- `speed_of_sound`: Speed of sound at the specified altitude.
- `mean_free_path`: Mean free path length at the specified altitude.

### Notes:

- If the specified altitude is outside the range of the model (up to 700 km), the function
  returns `missing` for all properties.
- The model assumes that within each atmospheric layer, the temperature varies linearly
  with altitude or remains constant (isothermal layer).
- The use of Unitful.jl means that units are specified in variable definition

"""

using Unitful
using Unitful.DefaultSymbols  # For unit symbols like m, K, kg, etc.

# Constants
const R_UNIVERSAL = 8_314.4621u"J/(kmol*K)" # Universal gas constant
const G0 = 9.80665u"m/s^2"                  # Acceleration due to gravity
const NA = 6.0221e26u"1/kmol"               # Avogadro's number
const SIGMA = 3.65e-10u"m"                  # Effective diameter of atmospheric gas molecule
const M0 = 28.9644u"kg/kmol"                # Molecular weight of air at sea level
const RP = 6.3781e6u"m"                     # Planetary radius
const R_SPECIFIC = R_UNIVERSAL / M0         # Specific gas constant for air
const GAMMA = 1.4u"1"                       # Ratio of specific heats for air (unitless)

# Atmospheric layer data: (Altitude, Temperature, Molecular Weight)
const DATA = [
    (11_000u"m",     216.65u"K",   28.9644u"kg/kmol"),
    (20_000u"m",     216.65u"K",   28.9644u"kg/kmol"),
    (32_000u"m",     228.65u"K",   28.9644u"kg/kmol"),
    (47_000u"m",     270.65u"K",   28.9644u"kg/kmol"),
    (51_000u"m",     270.65u"K",   28.9644u"kg/kmol"),
    (71_000u"m",     214.65u"K",   28.9644u"kg/kmol"),
    (86_000u"m",     186.946u"K",  28.952u"kg/kmol"),
    (100_000u"m",    210.65u"K",   28.88u"kg/kmol"),
    (110_000u"m",    260.65u"K",   28.56u"kg/kmol"),
    (120_000u"m",    360.65u"K",   28.07u"kg/kmol"),
    (150_000u"m",    960.65u"K",   26.92u"kg/kmol"),
    (160_000u"m",   1110.65u"K",   26.66u"kg/kmol"),
    (170_000u"m",   1210.65u"K",   26.5u"kg/kmol"),
    (190_000u"m",   1350.65u"K",   25.85u"kg/kmol"),
    (230_000u"m",   1550.65u"K",   24.69u"kg/kmol"),
    (300_000u"m",   1830.65u"K",   22.66u"kg/kmol"),
    (400_000u"m",   2160.65u"K",   19.94u"kg/kmol"),
    (500_000u"m",   2420.65u"K",   17.94u"kg/kmol"),
    (600_000u"m",   2590.65u"K",   16.84u"kg/kmol"),
    (700_000u"m",   2700.00u"K",   16.17u"kg/kmol")
]

# Function to calculate atmospheric properties
function atmosphere(altitude::Quantity{<:Real, Unitful.𝐋}; T₀=288.15u"K", P₀=101325.0u"Pa")
    # Ensure altitude has length units
    altitude = Unitful.upreferred(altitude)
    @assert u"m" == unit(altitude) "Altitude must have length units (e.g., meters)."

    # Check if altitude is within model range
    if altitude < 0u"m" || altitude > 700_000u"m"
        return (
            pressure = missing,
            density = missing,
            temperature = missing,
            speed_of_sound = missing,
            mean_free_path = missing
        )
    end

    # Initialise arrays with units
    num_layers = length(DATA)
    z_breakpoints = [0.0u"m"; [d[1] for d in DATA]]
    T_breakpoints = [T₀; [d[2] for d in DATA]]
    M_breakpoints = [M0; [d[3] for d in DATA]]

    # Compute lapse rates
    lapse_rates = diff(T_breakpoints) ./ diff(z_breakpoints)

    # Initialise pressure and density arrays
    P_breakpoints = zeros(num_layers + 1) .* u"Pa"
    rho_breakpoints = zeros(num_layers + 1) .* u"kg/m^3"
    P_breakpoints[1] = P₀
    rho_breakpoints[1] = P₀ / (R_SPECIFIC * T₀)

    # Calculate pressure and density at breakpoints
    for i in 1:num_layers
        if lapse_rates[i] != 0u"K/m"
            # Non-isothermal layer
            exponent = -G0 / (R_SPECIFIC * lapse_rates[i])
            temperature_ratio = T_breakpoints[i + 1] / T_breakpoints[i]
            P_breakpoints[i + 1] = P_breakpoints[i] * temperature_ratio^exponent
            rho_breakpoints[i + 1] = rho_breakpoints[i] * temperature_ratio^(exponent - 1)
        else
            # Isothermal layer
            exponent = -G0 * (z_breakpoints[i + 1] - z_breakpoints[i]) / (R_SPECIFIC * T_breakpoints[i])
            P_breakpoints[i + 1] = P_breakpoints[i] * exp(exponent)
            rho_breakpoints[i + 1] = rho_breakpoints[i] * exp(exponent)
        end
    end

    # Find the atmospheric layer for the given altitude
    i = findfirst(x -> altitude < x, z_breakpoints[2:end])
    if i === nothing
        return (
            pressure = missing,
            density = missing,
            temperature = missing,
            speed_of_sound = missing,
            mean_free_path = missing
        )
    end

    # Interpolate properties at the desired altitude
    if lapse_rates[i] != 0u"K/m"
        # Non-isothermal layer
        temperature = T_breakpoints[i] + lapse_rates[i] * (altitude - z_breakpoints[i])
        exponent = -G0 / (R_SPECIFIC * lapse_rates[i])
        temperature_ratio = temperature / T_breakpoints[i]
        pressure = P_breakpoints[i] * temperature_ratio^exponent
        density = rho_breakpoints[i] * temperature_ratio^(exponent - 1)
    else
        # Isothermal layer
        temperature = T_breakpoints[i]
        exponent = -G0 * (altitude - z_breakpoints[i]) / (R_SPECIFIC * T_breakpoints[i])
        pressure = P_breakpoints[i] * exp(exponent)
        density = rho_breakpoints[i] * exp(exponent)
    end

    # Interpolate molecular weight
    M_start = M_breakpoints[i]
    M_end = M_breakpoints[i + 1]
    z_start = z_breakpoints[i]
    z_end = z_breakpoints[i + 1]
    molecular_weight = M_start + (M_end - M_start) * (altitude - z_start) / (z_end - z_start)

    # Adjusted specific gas constant
    R_specific_adjusted = R_UNIVERSAL / molecular_weight

    # Speed of sound
    speed_of_sound = sqrt(GAMMA * R_specific_adjusted * temperature)

    # Mean free path
    average_molecular_speed = sqrt(8 * R_SPECIFIC * temperature / π)
    collision_frequency = sqrt(2) * π * NA * SIGMA^2 * average_molecular_speed * density / molecular_weight
    mean_free_path = average_molecular_speed / collision_frequency

    # Return results with units
    return (
        pressure = pressure,
        density = density,
        temperature = temperature,
        speed_of_sound = speed_of_sound,
        mean_free_path = mean_free_path
    )
end

export atmosphere

end  # module AtmosphereReganStandard
