import math, fluids
import CoolProp as CP

import functions as fn
import constants as cn

# State Classes
class FluidState:
    def __init__(self, P: float, T: float, rho: float, u: float, h: float, mu: float, x: float, phase: str):
        self.P = P
        self.T = T
        self.rho = rho
        self.u = u
        self.h = h
        self.mu = mu
        self.x = x
        self.phase = phase

class TankState:
    def __init__(self, m: float, V: float, U: float, fluid_state: FluidState):
        self.m = m
        self.V = V
        self.U = U
        self.fluid_state = fluid_state

# Fluid Class
class Fluid:
    def __init__(self, name: str, REF: str):
        self.name = name
        self.REF = REF
        self.RP = CP.AbstractState("HEOS", self.REF)

    def state_Ph(self, P: float, h: float) -> FluidState:
        self.RP.update(CP.HmassP_INPUTS, h, P)

        return FluidState(
            P = self.RP.p(),
            T = self.RP.T(),
            rho = self.RP.rhomass(),
            u = self.RP.umass(),
            h = self.RP.hmass(),
            mu = self.RP.viscosity(),
            x = self.RP.Q(),
            phase = str(self.RP.phase())
        )

    def state_PT(self, P: float, T: float) -> FluidState:
        self.RP.update(CP.PT_INPUTS, P, T)

        try:
            mu_val = self.RP.viscosity()
        except ValueError:
            mu_val = 0.0001  # Fallback dynamic viscosity [Pa*s] if transport model is missing

        return FluidState(
            P = self.RP.p(),
            T = self.RP.T(),
            rho = self.RP.rhomass(),
            u = self.RP.umass(),
            h = self.RP.hmass(),
            mu = mu_val,
            x = self.RP.Q(),
            phase = str(self.RP.phase())
        )

    def state_rhou(self, rho: float, u: float) -> FluidState:
      self.RP.update(CP.DmassUmass_INPUTS, rho, u)

      try:
        mu_val = self.RP.viscosity()
      except ValueError:
        mu_val = (
            0.0001  # Fallback dynamic viscosity [Pa*s] if transport model is missing
        )

      return FluidState(
          P=self.RP.p(),
          T=self.RP.T(),
          rho=self.RP.rhomass(),
          u=self.RP.umass(),
          h=self.RP.hmass(),
          mu=mu_val,
          x=self.RP.Q(),
          phase=str(self.RP.phase()),
      )

# Tank Class
class Tank:
    def __init__(self, name: str, fluid: Fluid, tank_state: TankState):
        self.name = name
        self.fluid = fluid
        self.tank_state = tank_state

# Piping Component Class - L0
class PipingComponent:
    def solve(self, mDot: float, fluid: Fluid, fluid_state_i: FluidState, RAE: float, i_max: int) -> FluidState:
        raise NotImplementedError

# Pipe Class - L1
class Pipe(PipingComponent):
    def __init__(self, name: str, L: float, dZ: float, ID: float, rgh: float):
        self.name = name
        self.L = L
        self.ID = ID
        self.rgh = rgh
        self.dZ = dZ

        self.A = math.pi * self.ID**2 / 4

    def dP_F(self, mDot: float, fluid_state_i: FluidState, fluid_state_f: FluidState) -> float:
        rho = (fluid_state_i.rho + fluid_state_f.rho) / 2
        mu = (fluid_state_i.mu + fluid_state_f.mu) / 2

        v = fn.v(mDot, rho, self.A)
        Re = fluids.Reynolds(V = v, D = self.ID, rho = rho, mu = mu)
        DFF = fluids.friction_factor(Re = Re, eD = self.rgh / self.ID)
        DP = fn.DP(rho, v)

        return DFF * self.L * DP / self.ID

    def dP_A(self, mDot: float, fluid_state_i: FluidState, fluid_state_f: FluidState) -> float:
        v_i = fn.v(mDot, fluid_state_i.rho, self.A)
        v_f = fn.v(mDot, fluid_state_f.rho, self.A)

        DP_i = fn.DP(fluid_state_i.rho, v_i)
        DP_f = fn.DP(fluid_state_f.rho, v_f)

        return DP_f - DP_i

    def dP_Z(self, fluid_state_i: FluidState, fluid_state_f: FluidState) -> float:
        rho = (fluid_state_i.rho + fluid_state_f.rho) / 2

        return rho * cn.CONST["g"] * self.dZ

    def dP(self, mDot: float, fluid_state_i: FluidState, fluid_state_f: FluidState) -> float:
        return (
            self.dP_F(mDot, fluid_state_i, fluid_state_f)
            + self.dP_A(mDot, fluid_state_i, fluid_state_f)
            + self.dP_Z(fluid_state_i, fluid_state_f)
        )

    def solve(self, mDot: float, fluid: Fluid, fluid_state_i: FluidState, RAE: float, i_max: int) -> FluidState:
        v_i = fn.v(mDot, fluid_state_i.rho, self.A)

        # Initial Guess
        v_f = v_i
        h_f = fluid_state_i.h - cn.CONST["g"] * self.dZ

        Re = fluids.Reynolds(
            V = v_i,
            D = self.ID,
            rho = fluid_state_i.rho,
            mu = fluid_state_i.mu
        )

        DFF = fluids.friction_factor(
            Re = Re,
            eD = self.rgh / self.ID
        )

        DP = fn.DP(fluid_state_i.rho, v_i)

        dP_F = DFF * self.L * DP / self.ID
        dP_Z = fluid_state_i.rho * cn.CONST["g"] * self.dZ

        P_f = fluid_state_i.P - dP_F - dP_Z

        fluid_state_f = fluid.state_Ph(P_f, h_f)

        error = float("inf")
        i = 0

        while error > RAE and i < i_max:
            P_prev = P_f
            h_prev = h_f

            v_f = fn.v(mDot, fluid_state_f.rho, self.A)

            h_f = (
                fluid_state_i.h
                + (v_i**2 - v_f**2) / 2
                - cn.CONST["g"] * self.dZ
            )

            dP = self.dP(mDot, fluid_state_i, fluid_state_f)
            P_f = fluid_state_i.P - dP

            fluid_state_f = fluid.state_Ph(P_f, h_f)

            RAE_P = abs((P_f - P_prev) / P_f)
            RAE_h = abs((h_f - h_prev) / h_f)

            error = max(RAE_P, RAE_h)

            i += 1

        return fluid_state_f

# Valve Class - L1
class Valve(PipingComponent):
    def __init__(self, name: str, C_v: float):
        self.name = name
        self.C_v = C_v

        self.C_v_SI = self.C_v * (cn.CONV["gal"] / 60) / math.sqrt(cn.CONV["psi"])

    def dP(self, mDot: float, fluid_state: FluidState) -> float:
        Q = fn.Q(mDot, fluid_state.rho)
        SG = fn.SG(fluid_state.rho)

        return SG * (Q / self.C_v_SI)**2

    def solve(self, mDot: float, fluid: Fluid, fluid_state_i: FluidState, RAE: float, i_max: int) -> FluidState:
        P_f = fluid_state_i.P - self.dP(mDot, fluid_state_i)
        h_f = fluid_state_i.h

        return fluid.state_Ph(P_f, h_f)

# Fitting Class - L1
class Fitting(PipingComponent):
    def __init__(self, name: str, ID_i: float, ID_f: float, dZ: float):
        self.name = name
        self.ID_i = ID_i
        self.ID_f = ID_f
        self.dZ = dZ

        self.A_i = math.pi * self.ID_i**2 / 4
        self.A_f = math.pi * self.ID_f**2 / 4

    def dP_K(self, mDot: float, fluid_state_i: FluidState, fluid_state_f: FluidState) -> float:
        raise NotImplementedError

    def dP_A(self, mDot: float, fluid_state_i: FluidState, fluid_state_f: FluidState) -> float:
        v_i = fn.v(mDot, fluid_state_i.rho, self.A_i)
        v_f = fn.v(mDot, fluid_state_f.rho, self.A_f)

        DP_i = fn.DP(fluid_state_i.rho, v_i)
        DP_f = fn.DP(fluid_state_f.rho, v_f)

        return DP_f - DP_i

    def dP_Z(self, fluid_state_i: FluidState, fluid_state_f: FluidState) -> float:
        rho = (fluid_state_i.rho + fluid_state_f.rho) / 2

        return rho * cn.CONST["g"] * self.dZ

    def dP(self, mDot: float, fluid_state_i: FluidState, fluid_state_f: FluidState) -> float:
        return (
            self.dP_K(mDot, fluid_state_i, fluid_state_f)
            + self.dP_A(mDot, fluid_state_i, fluid_state_f)
            + self.dP_Z(fluid_state_i, fluid_state_f)
        )

    def solve(self, mDot: float, fluid: Fluid, fluid_state_i: FluidState, RAE: float, i_max: int) -> FluidState:
        v_i = fn.v(mDot, fluid_state_i.rho, self.A_i)

        # Initial Guess
        h_f = fluid_state_i.h - cn.CONST["g"] * self.dZ
        P_f = fluid_state_i.P - fluid_state_i.rho * cn.CONST["g"] * self.dZ

        fluid_state_f = fluid.state_Ph(P_f, h_f)

        error = float("inf")
        i = 0

        while error > RAE and i < i_max:
            P_prev = P_f
            h_prev = h_f

            v_f = fn.v(mDot, fluid_state_f.rho, self.A_f)

            h_f = (
                fluid_state_i.h
                + (v_i**2 - v_f**2) / 2
                - cn.CONST["g"] * self.dZ
            )

            dP = self.dP(mDot, fluid_state_i, fluid_state_f)
            P_f = fluid_state_i.P - dP

            fluid_state_f = fluid.state_Ph(P_f, h_f)

            RAE_P = abs((P_f - P_prev) / P_f)
            RAE_h = abs((h_f - h_prev) / h_f)

            error = max(RAE_P, RAE_h)

            i += 1

        return fluid_state_f

# Rounded Entrance Class - L2
class RoundedEntrance(Fitting):
    def __init__(self, name: str, ID_i: float, ID_f: float, R: float, dZ: float, method: str = "Rennels"):
        super().__init__(name, ID_i, ID_f, dZ)

        self.R = R
        self.method = method

    def dP_K(self, mDot: float, fluid_state_i: FluidState, fluid_state_f: FluidState) -> float:
        v = fn.v(mDot, fluid_state_f.rho, self.A_f)

        K = fluids.fittings.entrance_rounded(
            Di = self.ID_f,
            rc = self.R,
            method = self.method
        )

        return K * fn.DP(fluid_state_f.rho, v)

# Sharp Entrance Class - L2
class SharpEntrance(Fitting):
    def __init__(self, name: str, ID_i: float, ID_f: float, dZ: float, method: str = "Rennels"):
        super().__init__(name, ID_i, ID_f, dZ)

        self.method = method

    def dP_K(self, mDot: float, fluid_state_i: FluidState, fluid_state_f: FluidState) -> float:
        v = fn.v(mDot, fluid_state_f.rho, self.A_f)

        K = fluids.fittings.entrance_sharp(
            method = self.method
        )

        return K * fn.DP(fluid_state_f.rho, v)

# Beveled Entrance Class - L2
class BeveledEntrance(Fitting):
    def __init__(self, name: str, ID_i: float, ID_f: float, L: float, angle: float, dZ: float, method: str = "Rennels"):
        super().__init__(name, ID_i, ID_f, dZ)

        self.L = L
        self.angle = angle
        self.method = method

    def dP_K(self, mDot: float, fluid_state_i: FluidState, fluid_state_f: FluidState) -> float:
        v = fn.v(mDot, fluid_state_f.rho, self.A_f)

        K = fluids.fittings.entrance_beveled(
            Di = self.ID_f,
            l = self.L,
            angle = self.angle,
            method = self.method
        )

        return K * fn.DP(fluid_state_f.rho, v)

# Rounded Bend Class - L2
class Bend(Fitting):
    def __init__(self, name: str, ID_i: float, ID_f: float, ID_K: float, R: float, angle: float, rgh: float, L_unimpeded: float, dZ: float, method: str = "Rennels"):
        super().__init__(name, ID_i, ID_f, dZ)

        self.ID_K = ID_K
        self.R = R
        self.angle = angle
        self.rgh = rgh
        self.L_unimpeded = L_unimpeded
        self.method = method

        self.A_K = math.pi * self.ID_K**2 / 4

    def dP_K(self, mDot: float, fluid_state_i: FluidState, fluid_state_f: FluidState) -> float:
        rho = (fluid_state_i.rho + fluid_state_f.rho) / 2
        mu = (fluid_state_i.mu + fluid_state_f.mu) / 2

        v = fn.v(mDot, rho, self.A_K)
        Re = fluids.Reynolds(V = v, D = self.ID_K, rho = rho, mu = mu)

        K = fluids.fittings.bend_rounded(
            Di = self.ID_K,
            angle = self.angle,
            rc = self.R,
            Re = Re,
            roughness = self.rgh,
            L_unimpeded = self.L_unimpeded,
            method = self.method
        )

        return K * fn.DP(rho, v)

# Miter Bend Class - L2
class MiterBend(Fitting):
    def __init__(self, name: str, ID_i: float, ID_f: float, ID_K: float, angle: float, rgh: float, L_unimpeded: float, dZ: float, method: str = "Rennels"):
        super().__init__(name, ID_i, ID_f, dZ)

        self.ID_K = ID_K
        self.angle = angle
        self.rgh = rgh
        self.L_unimpeded = L_unimpeded
        self.method = method

        self.A_K = math.pi * self.ID_K**2 / 4

    def dP_K(self, mDot: float, fluid_state_i: FluidState, fluid_state_f: FluidState) -> float:
        rho = (fluid_state_i.rho + fluid_state_f.rho) / 2
        mu = (fluid_state_i.mu + fluid_state_f.mu) / 2

        v = fn.v(mDot, rho, self.A_K)
        Re = fluids.Reynolds(V = v, D = self.ID_K, rho = rho, mu = mu)

        K = fluids.fittings.bend_miter(
            angle = self.angle,
            Di = self.ID_K,
            Re = Re,
            roughness = self.rgh,
            L_unimpeded = self.L_unimpeded,
            method = self.method
        )

        return K * fn.DP(rho, v)

# Contraction Class - L2
class Contraction(Fitting):
    def __init__(self, name: str, ID_i: float, ID_f: float, L: float, angle: float, rgh: float, dZ: float, method: str = "Rennels"):
        super().__init__(name, ID_i, ID_f, dZ)

        self.L = L
        self.angle = angle
        self.rgh = rgh
        self.method = method

    def dP_K(self, mDot: float, fluid_state_i: FluidState, fluid_state_f: FluidState) -> float:
        v = fn.v(mDot, fluid_state_f.rho, self.A_f)

        Re = fluids.Reynolds(
            V = v,
            D = self.ID_f,
            rho = fluid_state_f.rho,
            mu = fluid_state_f.mu
        )

        K = fluids.fittings.contraction_conical(
            Di1 = self.ID_i,
            Di2 = self.ID_f,
            l = self.L,
            angle = self.angle,
            Re = Re,
            roughness = self.rgh,
            method = self.method
        )

        return K * fn.DP(fluid_state_f.rho, v)

# Expansion Class - L2
class Expansion(Fitting):
    def __init__(self, name: str, ID_i: float, ID_f: float, L: float, angle: float, rgh: float, dZ: float, method: str = "Rennels"):
        super().__init__(name, ID_i, ID_f, dZ)

        self.L = L
        self.angle = angle
        self.rgh = rgh
        self.method = method

    def dP_K(self, mDot: float, fluid_state_i: FluidState, fluid_state_f: FluidState) -> float:
        v = fn.v(mDot, fluid_state_i.rho, self.A_i)

        Re = fluids.Reynolds(
            V = v,
            D = self.ID_i,
            rho = fluid_state_i.rho,
            mu = fluid_state_i.mu
        )

        K = fluids.fittings.diffuser_conical(
            Di1 = self.ID_i,
            Di2 = self.ID_f,
            l = self.L,
            angle = self.angle,
            Re = Re,
            roughness = self.rgh,
            method = self.method
        )

        return K * fn.DP(fluid_state_i.rho, v)

# Exit Class - L2
class Fitting_Exit(Fitting):
    def __init__(self, name: str, ID_i: float, ID_f: float, dZ: float):
        super().__init__(name, ID_i, ID_f, dZ)

    def dP_K(self, mDot: float, fluid_state_i: FluidState, fluid_state_f: FluidState) -> float:
        v = fn.v(mDot, fluid_state_i.rho, self.A_i)

        K = fluids.fittings.exit_normal()

        return K * fn.DP(fluid_state_i.rho, v)

# Injector Class
class Injector:

  def __init__(self, C_d: float, A: float):
    self.C_d = C_d
    self.A = A

  def solve(self, *args, **kwargs):
    # Extract and return just the fluid state object expected by Line.solve()
    state_i = args[1] if len(args) > 1 else (args[0] if len(args) > 0 else None)
    return state_i

  def dP(self, mDot: float, fluid_state) -> float:
    return (mDot / (self.C_d * self.A)) ** 2 / (2 * fluid_state.rho)

# Line Class
class Line:
    def __init__(self, name: str, source: Tank):
        self.name = name
        self.source = source
        self.fluid = source.fluid
        self.components = []

    def add(self, component):
        if isinstance(component, list):
            self.components.extend(component)
        else:
            self.components.append(component)

    def solve(self, mDot: float, fluid_state_i: FluidState, RAE: float, i_max: int):
        fluid_state = fluid_state_i
        fluid_state_history = [fluid_state_i]

        for component in self.components:
            fluid_state = component.solve(
                mDot,
                self.fluid,
                fluid_state,
                RAE,
                i_max
            )

            fluid_state_history.append(fluid_state)

        return fluid_state, fluid_state_history

# Engine Class
class Engine:
    def __init__(self, injector_OX: Injector, injector_F: Injector, A_t: float):
        self.injector_OX = injector_OX
        self.injector_F = injector_F
        self.A_t = A_t

