module AtmosphereReganStandard

"""
This module provides atmospheric properties.

The 1976 US Standard Atmosphere (US76) is used up to 86 km; for altitudes above 86 km
the model is based upon the 1962 US Standard Atmosphere (US62).

This approach is described in:
    Regan, F. J., & Anandakrishnan, S. M. (1993). "Dynamics of Atmospheric Re-Entry".

Regan notes:
    "The model may easily be changed to meet another standard or to accommodate different
    views on the value that should be set for the exospheric temperature. The model
    specifies that thermal layers be identified and that within each layer the temperature
    varies at most linearly with altitude."

### Arguments:

- `altitude`: Altitude in metres at which to calculate atmospheric properties.
- `T₀`: Ground level temperature in Kelvin. Optional.
- `P₀` Ground level pressure in N/m². Optional.

### Returns:

A tuple containing:

- `pressure`:Atmospheric pressure at the specified altitude (N/m²).
- `density`:Atmospheric density at the specified altitude (kg/m³).
- `temperature`:Kinetic temperature at the specified altitude (K).
- `speed_of_sound`: Speed of sound at the specified altitude (m/s).
- `mean_free_path`: Mean free path length at the specified altitude (m).

### Notes:

- If the specified altitude is outside the range of the model (up to 700 km), the function
  returns `NaN` for all properties.
- The model assumes that within each atmospheric layer, the temperature varies linearly
  with altitude or remains constant (isothermal layer).

"""

# Function to calculate atmospheric properties
function atmosphere(altitude, T₀=288.15, P₀=101325.0)
    # Constants
    R_universal = 8313.432        # Universal gas constant (J·kg⁻¹·K⁻¹)
    g₀ = 9.7803                   # Acceleration due to gravity (m/s²)
    Nₐ = 6.0221e26                # Avogadro's number (1/kmol)
    σ = 3.65e-10                  # Effective diameter of atmospheric gas molecule (m)
    M₀ = 28.9644                  # Molecular weight of air at sea level (kg/kmol)
    Rₚ = 6.3781e6                 # Planetary radius (m)
    β = 2.0 / Rₚ                  # Reciprocal of planetary radius
    R_specific = R_universal / M₀ # Specific gas constant for air (J·kg⁻¹·K⁻¹)

    # Initialize arrays
    z_breakpoints = zeros(21)     # Altitude breakpoints (m)
    T_breakpoints = zeros(21)     # Temperature at breakpoints (K)
    ρ_breakpoints = zeros(21)     # Density at breakpoints (kg/m³)
    P_breakpoints = zeros(21)     # Pressure at breakpoints (N/m²)
    M_breakpoints = zeros(21)     # Molecular weight at breakpoints (kg/kmol)
    lapse_rates = zeros(20)       # Lapse rates (K/m)

    # Set initial conditions at sea level
    z_breakpoints[1] = 0.0
    T_breakpoints[1] = T₀
    M_breakpoints[1] = M₀
    P_breakpoints[1] = P₀
    ρ_breakpoints[1] = P_breakpoints[1] / (R_specific * T_breakpoints[1])

    # Data for atmospheric layers (converted from km to m)
    data = [
        (11.0191e3, 216.65, 28.964),
        (20.0631e3, 216.65, 28.964),
        (32.1619e3, 228.65, 28.964),
        (47.3501e3, 270.65, 28.964),
        (51.4125e3, 270.65, 28.964),
        (71.8020e3, 214.65, 28.964),
        (86.00e3,   186.946,28.962),
        (100.00e3,  210.65, 28.880),
        (110.00e3,  260.65, 28.560),
        (120.00e3,  360.65, 28.070),
        (150.00e3,  960.65, 26.920),
        (160.00e3, 1110.60, 26.660),
        (170.00e3, 1210.65, 26.500),
        (190.00e3, 1350.65, 25.850),
        (230.00e3, 1550.65, 24.690),
        (300.00e3, 1830.65, 22.660),
        (400.00e3, 2160.65, 19.940),
        (500.00e3, 2420.65, 17.940),
        (600.00e3, 2590.65, 16.840),
        (700.00e3, 2700.00, 16.170)
    ]

    # Populate arrays with data
    for i in 2:21
        z_breakpoints[i], T_breakpoints[i], M_breakpoints[i] = data[i - 1]
    end

    # Calculate lapse rates for each atmospheric layer
    for i in 1:20
        lapse_rates[i] = (T_breakpoints[i + 1] - T_breakpoints[i]) / (z_breakpoints[i + 1] - z_breakpoints[i])
    end

    # Calculate pressure and density at each breakpoint
    for i in 1:20
        if abs(lapse_rates[i]) > 0.001
            # Non-isothermal layer
            A = 1 + β * ((T_breakpoints[i] / lapse_rates[i]) - z_breakpoints[i])
            exponent = (A * g₀) / (R_specific * lapse_rates[i])
            temperature_ratio = T_breakpoints[i + 1] / T_breakpoints[i]
            pressure_ratio = temperature_ratio^(-exponent)
            exp_factor = exp((g₀ * β * (z_breakpoints[i + 1] - z_breakpoints[i])) / (R_specific * lapse_rates[i]))
            P_breakpoints[i + 1] = P_breakpoints[i] * pressure_ratio * exp_factor
            density_exponent = exponent + 1
            ρ_breakpoints[i + 1] = ρ_breakpoints[i] * exp_factor * temperature_ratio^(-density_exponent)
        else
            # Isothermal layer
            exponent = -g₀ * (z_breakpoints[i + 1] - z_breakpoints[i]) * (1 - (β / 2) * (z_breakpoints[i + 1] + z_breakpoints[i])) / (R_specific * T_breakpoints[i])
            P_breakpoints[i + 1] = P_breakpoints[i] * exp(exponent)
            ρ_breakpoints[i + 1] = ρ_breakpoints[i] * exp(exponent)
        end
    end

    # Initialise output variables
    pressure = NaN
    density = NaN
    temperature = NaN
    speed_of_sound = NaN
    mean_free_path = NaN

    # Interpolate properties at the desired altitude
    for i in 1:20
        if altitude < z_breakpoints[i + 1]
            if abs(lapse_rates[i]) >= 0.001
                # Non-isothermal layer
                temperature = T_breakpoints[i] + lapse_rates[i] * (altitude - z_breakpoints[i])
                A = 1 + β * ((T_breakpoints[i] / lapse_rates[i]) - z_breakpoints[i])
                exponent = (A * g₀) / (R_specific * lapse_rates[i])
                temperature_ratio = temperature / T_breakpoints[i]
                pressure_ratio = temperature_ratio^(-exponent)
                exp_factor = exp((β * g₀ * (altitude - z_breakpoints[i])) / (R_specific * lapse_rates[i]))
                pressure = P_breakpoints[i] * pressure_ratio * exp_factor
                density_exponent = exponent + 1.0
                density = ρ_breakpoints[i] * (temperature_ratio^(-density_exponent)) * exp_factor
            else
                # Isothermal layer
                temperature = T_breakpoints[i]
                exponent = -g₀ * (altitude - z_breakpoints[i]) * (1.0 - (β / 2) * (altitude + z_breakpoints[i])) / (R_specific * T_breakpoints[i])
                pressure = P_breakpoints[i] * exp(exponent)
                density = ρ_breakpoints[i] * exp(exponent)
            end

            # Speed of sound calculation
            speed_of_sound = sqrt(1.4 * R_specific * temperature)

            # Molecular weight interpolation
            molecular_weight = M_breakpoints[i] + (M_breakpoints[i + 1] - M_breakpoints[i]) / (z_breakpoints[i + 1] - z_breakpoints[i]) * (altitude - z_breakpoints[i])

            # Adjusted temperature
            temperature = molecular_weight * temperature / M₀

            # Dynamic viscosity (Sutherland's formula)
            μ = 1.458e-6 * temperature^(1.5) / (temperature + 110.4)

            # Kinematic viscosity
            ν = μ / density

            # Thermal conductivity
            κ = 2.64638e-3 * temperature^(1.5) / (temperature + 245.4 * (10^(-12 / temperature)))

            # Average molecular speed
            average_molecular_speed = sqrt(8 * R_specific * temperature / π)

            # Collision frequency
            collision_frequency = sqrt(2) * π * Nₐ * σ^2 * average_molecular_speed

            # Mean free path
            mean_free_path = average_molecular_speed / collision_frequency
            break
        end
    end

    return pressure, density, temperature, speed_of_sound, mean_free_path
end

export atmosphere

end 
