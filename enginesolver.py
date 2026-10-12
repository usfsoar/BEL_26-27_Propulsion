import numpy as np
from scipy.optimize import root

G0 = 9.81  # m/s^2

# Line parameters (same values as the MATLAB script)
OX_LINE = dict(f=0.04, D=0.0535088, L=0.762, k=10, Cv=9,     A_inj=2.20893e-5, dz=0.9144)
F_LINE  = dict(f=0.04, D=0.01270,   L=0.762, k=10, Cv=5.228, A_inj=1.8816e-3,  dz=0.9144)


def line_loss(m_dot, rho, f, D, L, k, Cv, A_inj, dz):
    """Total pressure drop (Pa) across a feed line."""
    rho_safe = max(rho, 1e-3)

    # velocity and dynamic pressure
    A_pipe = (np.pi / 4) * D**2
    v = m_dot / (rho_safe * A_pipe)
    q = 0.5 * rho_safe * v**2

    dp_pipe = f * (L / D) * q   # pipe friction loss
    dp_k = k * q                # minor losses

    # valve loss: convert m_dot to GPM for the Cv equation
    SG = rho_safe / 1000.0
    Q_gpm = (m_dot / rho_safe) * 15850.32   # m^3/s -> GPM
    dp_v = SG * (Q_gpm / Cv) ** 2 * 6894.76  # psi -> Pa

    dp_inj = 0.5 * (m_dot / A_inj) ** 2 / rho_safe  # injector loss
    dp_g = rho_safe * G0 * dz                       # hydrostatic loss

    return dp_pipe + dp_k + dp_v + dp_inj + dp_g


def engine_residuals(x, P_ox, P_f, rho_ox, rho_f, A_t, c_star):
    m_dot_ox = abs(x[0])
    m_dot_f = abs(x[1])
    Pc = abs(x[2])

    dp_ox = line_loss(m_dot_ox, rho_ox, **OX_LINE)
    dp_f = line_loss(m_dot_f, rho_f, **F_LINE)

    R = np.array([
        P_ox - dp_ox - Pc,
        P_f - dp_f - Pc,
        (m_dot_ox + m_dot_f) - (Pc * A_t / c_star),
    ])

    # keep NaN/inf away from the solver
    if not np.all(np.isfinite(R)):
        R = np.array([1e6, 1e6, 1e6])
    return R


def engine_solver(m_ox, m_f, P_ox, P_f, rho_ox, rho_f):
    """Returns (m_dot_ox, m_dot_f, Pc, F, O_F, Isp). m_ox and m_f are unused, as in the original."""
    # constants
    c_star = 2639.4          # m/s
    A_t = 7.85398e-5         # m^2
    A_e = 1.8816e-3          # m^2
    C_f = A_e / A_t          # kept identical to the MATLAB script

    # initial guesses [m_dot_ox (kg/s), m_dot_f (kg/s), Pc (Pa)]
    x0 = np.array([1.66337, 0.478249, 435.113 * 6894.76])

    with np.errstate(all="ignore"):
        sol = root(engine_residuals, x0,
                   args=(P_ox, P_f, rho_ox, rho_f, A_t, c_star),
                   method="lm")

    x_sol = sol.x if sol.success else x0   # fall back to initial guess on failure

    m_dot_ox = abs(x_sol[0])
    m_dot_f = abs(x_sol[1])
    Pc = abs(x_sol[2])

    # performance outputs
    m_dot_tot = m_dot_ox + m_dot_f
    F = C_f * Pc * A_t
    O_F = m_dot_ox / max(m_dot_f, 1e-6)
    Isp = F / (max(m_dot_tot, 1e-6) * G0)

    return float(m_dot_ox), float(m_dot_f), float(Pc), float(F), float(O_F), float(Isp)


if __name__ == "__main__":
    # example call (replace with your own values)
    out = engine_solver(0, 0, 600 * 6894.76, 600 * 6894.76, 1141.0, 800.0)
    print(dict(zip(["m_dot_ox", "m_dot_f", "Pc", "F", "O_F", "Isp"], out)))