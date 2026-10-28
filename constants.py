CONV = {
    # Pressure [Pa]
    "Pa": 1.0,
    "kPa": 1e3,
    "MPa": 1e6,
    "bar": 1e5,
    "psi": 6894.757293,

    # Length [m]
    "m": 1.0,
    "cm": 1e-2,
    "mm": 1e-3,
    "in": 0.0254,
    "ft": 0.3048,

    # Area [m^2]
    "m^2": 1.0,
    "cm^2": 1e-4,
    "mm^2": 1e-6,
    "in^2": 0.00064516,

    # Volume [m^3]
    "m^3": 1.0,
    "L": 1e-3,
    "mL": 1e-6,
    "in^3": 1.6387064e-5,
    "ft^3": 0.0283168466,
    "gal": 0.003785411784,

    # Mass [kg]
    "kg": 1.0,
    "g": 1e-3,
    "lbm": 0.45359237,

    # Mass Flow [kg/s]
    "kg/s": 1.0,
    "g/s": 1e-3,
    "lbm/s": 0.45359237,

    # Force [N]
    "N": 1.0,
    "lbf": 4.448221615,

    # Energy [J]
    "J": 1.0,
    "kJ": 1e3,

    # Specific Energy [J/kg]
    "J/kg": 1.0,
    "kJ/kg": 1e3,

    # Dynamic Viscosity [Pa*s]
    "Pa*s": 1.0,
    "cP": 1e-3,

    # Density [kg/m^3]
    "kg/m^3": 1.0,
    "g/cm^3": 1000.0,
}

CONST = {
    "g": 9.80665, # The gravitational constant. [m/s^2]
    "rho_STD": 999.0, # The density of H_(2)O at 4 + 273.15 K. [kg/m^3]
    "R": 8.314462618 # The ideal gas constant. [J/(mol*K)]
}