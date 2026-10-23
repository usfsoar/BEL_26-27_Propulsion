import math
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

import constants as cn
import functions as fn

OX = "NitrousOxide"
FUEL = "Ethanol"
PSI = cn.CONV["psi"]
G0 = cn.CONST["g"]

# --------- SETTINGS ---------
dt = 0.05          # [s]
t_max = 25.0       # [s] (stops when tank runs out)

USE_CONSTANT_CSTAR = False
CSTAR_CONST = 2639.4   # [m/s]

C_f = 1.0          # replace with a real value

# --------- DESIGN (from notes) ---------
T_o = 15.0 + 273.15                      # [K] both tanks
V_tot = 3265 * cn.CONV["in^3"]           # [m^3] total tank volume (piston-separated)
V_OX_o = 0.0413082                       # [m^3] initial oxidizer volume
dP_piston = 15 * PSI                     # [Pa] fuel pressure = N2O pressure - 15 psi

A_t = 1.8816e-3                          # [m^2] throat
C_d = 0.7                                # injector discharge coefficient
L_pipe = 0.762                           # [m] line length
dz = 0.9144                              # [m] vertical rise (3 ft)
K_fit = 1.0                              # fitting loss coefficient

# n = number of parallel lines (oxidizer flow splits over 2), f = Darcy friction factor, Cv in gpm/sqrt(psi)
OX_LINE = dict(A_inj=7.85398e-5, D=0.5 * cn.CONV["in"], n=2, f=0.0133, Cv=9.0)
F_LINE = dict(A_inj=2.20893e-5, D=0.375 * cn.CONV["in"], n=1, f=0.01975, Cv=5.28)

# --------- LINE & INJECTOR MODEL ---------
def line_dP(m, rho, ln):
    """Flow-dependent pressure loss [Pa] (pipe friction + fitting + valve + injector) at total flow m [kg/s]."""
    Q = (m / ln["n"]) / rho                                   # volumetric flow per line [m^3/s]
    q = 0.5 * rho * (Q / (math.pi * ln["D"] ** 2 / 4)) ** 2   # dynamic pressure [Pa]
    dP_pipe = ln["f"] * (L_pipe / ln["D"]) * q
    dP_K = K_fit * q
    dP_valve = fn.SG(rho) * (Q * 15850.32 / ln["Cv"]) ** 2 * PSI   # SG*(Q[gpm]/Cv)^2 [psi]
    dP_inj = 0.5 / rho * (m / (C_d * ln["A_inj"])) ** 2
    return dP_pipe + dP_K + dP_valve + dP_inj


def mdot_at(P_tank, P_c, rho, ln):
    """Mass flow [kg/s] that the tank pressure can push against chamber pressure P_c."""
    P_avail = P_tank - rho * G0 * dz - P_c     # after hydrostatic head
    if P_avail <= 0:
        return 0.0
    m_hi = C_d * ln["A_inj"] * math.sqrt(2 * rho * P_avail)   # injector alone -> upper bound
    return root_scalar(lambda m: P_avail - line_dP(m, rho, ln), bracket=[0.0, m_hi], method="brentq").root


def cstar_val(P_c, OF):
    return CSTAR_CONST if USE_CONSTANT_CSTAR else fn.cstar(P_c, OF)


def tank_state(rho, u):
    err = None
    for r in (rho, rho * (1 + 1e-7), rho * (1 - 1e-7)):
        try:
            return (CP.PropsSI("P", "Dmass", r, "Umass", u, OX),
                    CP.PropsSI("T", "Dmass", r, "Umass", u, OX))
        except Exception as e:
            err = e
    raise err

# --------- INITIAL STATES ---------
# N2O: saturated liquid at T_o  
rho_l0 = CP.PropsSI("D", "T", T_o, "Q", 0, OX)
rho_v0 = CP.PropsSI("D", "T", T_o, "Q", 1, OX)
u_l0 = CP.PropsSI("U", "T", T_o, "Q", 0, OX)
u_v0 = CP.PropsSI("U", "T", T_o, "Q", 1, OX)
x0 = 1e-5
rho_ox0 = 1.0 / ((1 - x0) / rho_l0 + x0 / rho_v0)
u_ox0 = (1 - x0) * u_l0 + x0 * u_v0
P_OX_o = CP.PropsSI("P", "T", T_o, "Q", 0, OX)
m_ox = rho_ox0 * V_OX_o
U_ox = m_ox * u_ox0
T_crit_ox = CP.PropsSI("Tcrit", OX)

# for the ethanol
V_F_o = V_tot - V_OX_o
rho_F = CP.PropsSI("D", "T", T_o, "P", P_OX_o - dP_piston, FUEL)
m_f = rho_F * V_F_o

print(f"Initial N2O: P = {P_OX_o/PSI:.1f} psia, rho = {rho_ox0:.1f} kg/m3, m = {m_ox:.2f} kg   (notes: 666 psia, 820.5, 33.89)")
print(f"Initial fuel: P = {(P_OX_o-dP_piston)/PSI:.1f} psia, rho = {rho_F:.1f} kg/m3, m = {m_f:.2f} kg   (notes: 651 psia, 794, 9.68)")
print(f"c* source: {'constant ' + str(CSTAR_CONST) + ' m/s' if USE_CONSTANT_CSTAR else 'CEA table in functions.py'}")

# --------- TIME LOOP ---------
time = 0.0
history = []
stop_reason = "reached t_max"

while time <= t_max:
    if m_f <= 0.05:
        stop_reason = "fuel depleted"
        break

    # --- oxidizer tank state (piston moves as fuel leaves) ---
    V_ox = V_tot - m_f / rho_F
    rho_ox = m_ox / V_ox
    u_ox = U_ox / m_ox
    try:
        P_ox, T_ox = tank_state(rho_ox, u_ox)
        rho_inj_ox = rho_ox
        h_out_ox = CP.PropsSI("H", "Dmass", rho_ox, "Umass", u_ox, OX)
        if T_ox < T_crit_ox:
            rho_l = CP.PropsSI("D", "T", T_ox, "Q", 0, OX)
            rho_v = CP.PropsSI("D", "T", T_ox, "Q", 1, OX)
            if rho_ox <= rho_v:                      # all vapor: no liquid left to flow
                stop_reason = "N2O liquid depleted"
                break
            if rho_ox < rho_l:                       # two-phase: liquid leaves
                rho_inj_ox = rho_l
                h_out_ox = CP.PropsSI("H", "T", T_ox, "Q", 0, OX)
                x = (1 / rho_ox - 1 / rho_l) / (1 / rho_v - 1 / rho_l)
                if x > 0.99:
                    stop_reason = "N2O liquid depleted"
                    break
    except Exception as e:
        stop_reason = f"N2O property error: {e}"
        break
    P_f = P_ox - dP_piston

    # --- solve chamber pressure (mass balance), flows from tank -> lines -> injector ---
    def balance(P_c):
        mo = mdot_at(P_ox, P_c, rho_inj_ox, OX_LINE)
        mf = mdot_at(P_f, P_c, rho_F, F_LINE)
        OF = mo / mf if mf > 0 else 4.0
        return mo + mf - P_c * A_t / cstar_val(P_c, OF)

    P_lo = 0.1e6
    P_hi = min(P_ox - rho_inj_ox * G0 * dz, P_f - rho_F * G0 * dz) - 1e3
    try:
        P_c = root_scalar(balance, bracket=[P_lo, P_hi], method="brentq").root
    except Exception as e:
        stop_reason = f"chamber pressure solve failed: {e}"
        break

    mdot_ox = mdot_at(P_ox, P_c, rho_inj_ox, OX_LINE)
    mdot_f = mdot_at(P_f, P_c, rho_F, F_LINE)
    mdot_tot = mdot_ox + mdot_f
    F = C_f * P_c * A_t
    history.append(dict(time=time, P_c=P_c, P_ox=P_ox, P_f=P_f, mdot_OX=mdot_ox, mdot_F=mdot_f,
                        OF=mdot_ox / mdot_f, F=F, Isp=F / (mdot_tot * G0), T_ox=T_ox,
                        m_ox=m_ox, m_f=m_f))

    # --- advance: mass, and N2O energy (outflow enthalpy + work on the piston) ---
    dV_ox = mdot_f / rho_F * dt                      # N2O volume growth [m^3]
    U_ox -= mdot_ox * dt * h_out_ox + P_ox * dV_ox
    m_ox -= mdot_ox * dt
    m_f -= mdot_f * dt
    time += dt

print(f"\nSimulation stopped: {stop_reason}. Timesteps recorded: {len(history)}")
if not history:
    raise SystemExit("No timesteps recorded.")

# --------- SUMMARY ---------
h0, hN = history[0], history[-1]
print("\n--- Summary ---")
print(f"Burn time:  {hN['time'] + dt:.2f} s")
print(f"P_c:        {h0['P_c']/1e6:.2f} -> {hN['P_c']/1e6:.2f} MPa")
print(f"N2O tank:   {h0['P_ox']/1e6:.2f} -> {hN['P_ox']/1e6:.2f} MPa")
print(f"mdot_ox:    {h0['mdot_OX']:.3f} -> {hN['mdot_OX']:.3f} kg/s")
print(f"mdot_F:     {h0['mdot_F']:.3f} -> {hN['mdot_F']:.3f} kg/s")
print(f"O/F:        {h0['OF']:.2f} -> {hN['OF']:.2f}")
print(f"Thrust:     {h0['F']:.0f} -> {hN['F']:.0f} N (C_f = {C_f})")
print(f"Isp:        {h0['Isp']:.1f} -> {hN['Isp']:.1f} s")
print(f"Propellant used: N2O {history[0]['m_ox'] - m_ox:.2f} kg, ethanol {history[0]['m_f'] - m_f:.2f} kg")

# --------- PLOTS ---------
t = [h["time"] for h in history]
fig, ax = plt.subplots(2, 2, figsize=(12, 8))

ax[0, 0].plot(t, [h["P_c"] / 1e6 for h in history], label="P_c", color="crimson", lw=2)
ax[0, 0].plot(t, [h["P_ox"] / 1e6 for h in history], label="N2O tank", color="dodgerblue", ls="--")
ax[0, 0].plot(t, [h["P_f"] / 1e6 for h in history], label="Fuel tank", color="forestgreen", ls="--")
ax[0, 0].set_ylabel("Pressure [MPa]"); ax[0, 0].set_title("Pressures"); ax[0, 0].legend(); ax[0, 0].grid(True)

ax[0, 1].plot(t, [h["mdot_OX"] for h in history], label="Oxidizer", color="dodgerblue", lw=2)
ax[0, 1].plot(t, [h["mdot_F"] for h in history], label="Fuel", color="forestgreen", lw=2)
ax[0, 1].set_ylabel("Mass flow [kg/s]"); ax[0, 1].set_title("Mass flows"); ax[0, 1].legend(); ax[0, 1].grid(True)

ax[1, 0].plot(t, [h["OF"] for h in history], color="purple", lw=2)
ax[1, 0].set_ylabel("O/F [-]"); ax[1, 0].set_xlabel("Time [s]"); ax[1, 0].set_title("Mixture ratio"); ax[1, 0].grid(True)

ax[1, 1].plot(t, [h["F"] for h in history], color="darkorange", lw=2)
ax[1, 1].set_ylabel("Thrust [N]", color="darkorange"); ax[1, 1].set_xlabel("Time [s]"); ax[1, 1].grid(True)
ax2 = ax[1, 1].twinx()
ax2.plot(t, [h["Isp"] for h in history], color="black", lw=2, ls=":")
ax2.set_ylabel("Isp [s]")
ax[1, 1].set_title(f"Thrust and Isp (C_f = {C_f})")

plt.tight_layout()
plt.show()

