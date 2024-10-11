"""Module containing PerkinsChokeModel() to implement the Perkins choke model for single phase flow from the 1993 paper:

https://onepetro.org/DC/article/8/04/271/70440?casa_token=S0CnmdUh9ksAAAAA:tEaYx5JMdo4ZKz9bTK1Vnr7Jzu3-koOyLBNf51Nq558mUix4FifueBQ-4Ju7HzBO_az5KyM


TODO: This implementation needs to be properly tested against known mass flow rates (kg/h) and also the existing estimation methodology 
for different PFVs before using it in a production environment.
"""

import math

from src.root_finding import (
    bisection_method,
    secant,
)
from src.units import (
    convert_barg_to_psia,
    convert_psia_to_barg,
    convert_kg_h_to_lb_s,
    convert_lb_s_to_kg_h,
    convert_deg_c_to_deg_f,
    convert_deg_f_to_deg_rankine,
    convert_m_to_ft,
)


# PLEASE BE CAREFUL WITH UNITS IN THIS MODULE - CALCULATIONS ARE DONE IN AMERICAN "FREEDOM" UNITS BUT
# CONSTRUCTOR IS DONE TO CONSUME UNITS USED BY THE OTHER LIBRARY COMPONENTS.
GC = 32.2  # acceleration due to gravity.
RN = 10.736 * 144  # gas constant in freedom units.
RHO_AIR = 0.0764  # Lb/cu.ft - reference air density at standard conditions.
MW_AIR = 28.96  # air molecular weight.
CRITICAL_CONVERGENCE_FACTOR = 1.5  # convergence criteria.


class PerkinsChokeModel:
    """
    Implements the Perkins choke model for single-phase flow.

    This class provides calculations based on the Perkins choke model, which is commonly
    used in the oil and gas industry for estimating flow through chokes. The model is
    implemented using oilfield units to align with industry standards.
    The entered units are barg, deg C, metres, mm, specific gravity, J/mol/K, and dimensionless.

    Defining "critical flow" as no physical way for more fluid to pass through the orifice given the upstream.
    i.e. the fluid velocity is travelling at the local speed of sound, which is the speed limit for a fluid.
    """

    # root-finding method params.
    tol = 1e-3
    max_iter = 400

    def __init__(
        self,
        inlet_pressure,
        inlet_temperature,
        upstream_diameter,
        orifice_diameter,
        gas_gravity,
        gas_cv,
        gas_z_factor,
    ):
        """Initialisation."""
        # ** gas params **
        self.p1 = convert_barg_to_psia(inlet_pressure)
        self.t1 = convert_deg_c_to_deg_f(inlet_temperature)
        self.gas_gravity = gas_gravity
        self.mw = gas_gravity * MW_AIR  # molecular weight.
        self.cv = (
            1e3 * (gas_cv / self.mw)
        ) * 0.1858625352  # cv from  J/mol/K to lb.ft/lb/R.

        self.z = gas_z_factor
        self.gas_density = (
            62.43
            * 1.2238e-3
            * gas_gravity
            / (
                2.8301e-2
                * self.z
                * convert_deg_f_to_deg_rankine(self.t1)
                / inlet_pressure  # pressure in barg.
            )
        )
        self.v1 = 1.0 / self.gas_density
        self.RO_std = self.gas_gravity * RHO_AIR
        self.n = ((1.3 - 0.31 * (self.gas_gravity - 0.55)) * self.cv) / self.cv
        # weight fractions of gas, oil and water.
        # for this use case gas is 1 and gas + oil + water = 1
        # f_g + f_o + f_w = 1
        # therefore, oil = water = 0.
        self.f_g, self.f_w, self.f_o = 1.0, 0.0, 0.0
        self.alpha = 0.0  # 1/v * (fo/ρo + fw/ρw) = 0.
        self.choke_cd = 1.00
        # λ = f_g + [(f_g.C_vg + f_o.C_vo + f_w.C_vw)]M/zR
        self.lam = self.f_g + ((self.f_g * gas_cv * self.mw) / (self.z * RN))

        # ** pipe dims params **
        # all pipe dims are in ft.
        # get diameter in ft from m.
        self.orifice_diameter = convert_m_to_ft(orifice_diameter)
        self.upstream_diameter = convert_m_to_ft(upstream_diameter)
        self.aperture_area = self.calculate_pipe_area(diameter=self.orifice_diameter)
        self.upstream_area = self.calculate_pipe_area(diameter=self.upstream_diameter)
        # edge cases: (protect for weird cases - in secant and bisection).
        self.a1 = self.upstream_area
        self.a2 = self.aperture_area
        if self.a1 >= 3 * self.a2:
            self.a1 = 3 * self.a2
        if self.a1 <= self.a2:
            self.a2 = self.a1
        self.areas = self.a2 / self.a1

        # critical flow calculated params.
        self.p_r_critical = None
        self.m_critical = None

    @staticmethod
    def calculate_pipe_area(diameter):
        """Calculate pipe area using a diameter in any unit."""
        return math.pi * (diameter / 2) ** 2

    def update_orifice_size(self, x):
        """update class attr self.areas (orifice size) using an updated orifice diameter, x (m).

        This updates the areas ratio: A2 / A1.
        """
        # get diameter in ft (from m) and get the new area.
        updated_orifice_diameter = convert_m_to_ft(x)
        self.aperture_area = self.calculate_pipe_area(diameter=updated_orifice_diameter)
        a1 = self.upstream_area
        a2 = self.aperture_area
        if a1 >= 3 * a2:
            a1 = 3 * a2
        if a1 <= a2:
            a2 = a1
        self.areas = a2 / a1

    def perkins_pressure_ratio(self, p_r):
        """Pressure ratio, p_r = p2/p1 or p3/p1.

        This is an implementation of equation A-30 from Perkins 1993.
        For passing flare valve use case, p2 > p3 meaning that the flow is critical and p_r = p2/p1.
        For readability, the equation has been broken down into terms. The RHS terms are subtracted from the LHS
        terms and the entire equation is set equal to 0 so that it can be solved using a root finder algorithm.
        """
        # set out some common terms.
        # f_g.p_r^(-1/n) + alpha
        common_term1 = (self.f_g * pow(p_r, -1 / self.n)) + self.alpha
        # p_r^(-(1+n)/n)
        common_term2 = pow(p_r, (-(self.n + 1) / self.n))

        # ** build equation A-30 **
        # 2λ[1 - p_r^(n-1)/n]
        term1 = 2 * self.lam * (1 - pow(p_r, (self.n - 1) / self.n))
        # 2α_1(1 - p_r)
        term2 = 2 * self.alpha * (1 - p_r)
        lhs_bracket1 = term1 + term2

        # [1 - (A2/A1)^2 . ((f_g + α_1) / (f_g.p_r^(-1/n) + α_1))^2] . [f_g / n . p_r^(-(1+n)/n)]
        term3 = (
            1 - self.areas**2 * pow(((self.f_g + self.alpha) / common_term1), 2)
        ) * ((self.f_g / self.n) * common_term2)

        # (A2/A1)^2 . ((f_g + α_1) / (f_g.p_r^(-1/n) + α_1))^2
        term4 = (
            self.areas**2
            * (self.f_g / self.n)
            * ((pow(self.f_g + self.alpha, 2) * common_term2) / common_term1**2)
        )
        lhs_bracket2 = term3 + term4

        # (1 - (A2/A1)[(f_g + α_1) / (f_g.p_r^(-1/n) + α_1))]^2 . f_g.p_r^(-1/n) + α_1)
        term5 = 1 - (self.areas * pow((self.f_g + self.alpha) / common_term1, 2))
        rhs_bracket1 = term5 * common_term1

        # λ(n-1/n)(p_r^(-1/n) + α_1)
        rhs_bracket2 = (
            self.lam * ((self.n - 1) / self.n) * (pow(p_r, -1 / self.n) + self.alpha)
        )

        # collect terms LHS and set = 0.
        return (lhs_bracket1 * lhs_bracket2) - (rhs_bracket1 * rhs_bracket2)  # = 0.

    def get_estimate_critical_pressure_ratio(self):
        """Get an estimate of the critical pressure ratio, p_r (critical) by
        finding the root of the equation A-30 in Perkins 1993. using the Secant
        root finding method.

        The secant method finds the pressure ratio, p_r that satisfies equation A-30.
        """
        # initial x0 and x1 to find root.
        # mostly, liquid will be vastly different to all gas.
        # a small lower bound accounts for this but realistically,
        # for this use case we expect ~ 0.5.
        prf, p_st = 0.005, 0.0045

        return secant(
            func=self.perkins_pressure_ratio,
            x0=prf,
            x1=p_st,
            tol=self.tol,
            max_iterations=self.max_iter,
        )

    def downstream_velocity(self, p_r):
        """Calculate the downstream velocity (V2 in Perkins 1993.) from the Perkins paper.

        This implements equation A-25 from the paper.
        """
        numerator = (
            self.lam * self.p1 * self.v1 * (1 - pow(p_r, (self.n - 1) / self.n))
        )  # + 0 (f_o and f_w = 0).
        denominator = 1 - self.areas**2 * pow(
            (self.f_g + self.alpha) / ((self.f_g * pow(p_r, -1 / self.n)) + self.alpha),
            2,
        )
        denominator = max(
            denominator, 1e-6
        )  # 1e-6 is a sensible lower bound (ft/s) for speed.
        # adding in the LB constraint prevents small or large values causing issues.

        return (288 * GC * (numerator / denominator)) ** 0.5

    def mass_flow_rate(self, p_r):
        """Calculate the mass flow rate, w_c in lb/s from the Perkins 1993. paper.

        The method implements equation A-27.

        The velocity (V2) is first found using A-26. ρ_2 is then found by rearranging
        A-27 for ρ_2. w_i is then found using A2 * V2 * ρ_2. The predicted mass flow rate
        is then found using w_c = K * w_i where K = 0.826 (K taken directly from the paper).
        """
        v2 = self.downstream_velocity(p_r)  # eq. A-26.
        if v2 == 0.0 and p_r < 1.0:
            v2 = p_r * self.v1

        # rearranged A-27 to calculate ρ2.
        rho2 = 1.0 / (self.f_g * self.v1 * pow(p_r, -1 / self.n))

        # A-27.
        w_i = self.a2 * v2 * rho2

        # return the suggested predicted flow by multiplying w_i by K.
        K = 0.826
        w_c = K * w_i

        return w_c

    def calc_perkins_dp(self, m):
        """calculate the outlet pressure, p2 (barg) for a given gas mass rate (kg/h) by solving
        A-30 for p_r and then substituting into A-27 for w_i.

        Steps:
        1. convert the observed mass flow rate, w_o into lbs/s from kg/h.
        2. estimate the critical pressure ratio, p_r using the secant method.
        3. calculate the mass flow rate, w_c.
        4. calculate the outlet pressure, p2 using the mass flow rates, and return it.
        """

        # observed mass flow rate in lb/s.
        w_o = convert_kg_h_to_lb_s(m)

        def _mass_flow_rate_delta(x):
            """wrapper to calculate mass flux, w_c and return the delta with it and observed.
            x is p2 and it is used to find p_r by doing, p_r = p2 / p1.
            """
            p_r = x / self.p1
            w_c = self.mass_flow_rate(p_r=p_r)

            return w_c - w_o

        # solve Eq. A-30 for p_r.
        self.p_r_critical = self.get_estimate_critical_pressure_ratio()
        # calculated mass flow rate in lb/s.
        w_c = self.mass_flow_rate(p_r=self.p_r_critical)
        # update calculate mass flow rate in kg/h.
        self.m_critical = convert_lb_s_to_kg_h(w_c)

        if w_c > 0 and w_o > w_c:
            p2 = 0.0 - (((w_o - w_c) / w_c) * CRITICAL_CONVERGENCE_FACTOR)
        else:
            # update x0 and x1 and find a new mass flow rate using secant method.
            # trying to get observed and calculated the same.
            # assuming p2 < p1 **always**, we want to make sure that we have p1 as the upper bound.
            x0 = 0.95 * self.p1
            x1 = 0.85 * x0
            p2 = secant(
                func=_mass_flow_rate_delta,
                x0=x0,
                x1=x1,
                tol=self.tol,
                max_iterations=self.max_iter,
            )  # psia
            p2 = convert_psia_to_barg(p2)

        return p2


# TODO: Need to fix this. Currently it is not converging.
def calculate_perkins_equivalent_orifice(
    inlet_pressure,
    inlet_temperature,
    mass_flow_rate,
    upstream_diameter,
    gas_gravity,
    gas_cv,
    gas_z_factor,
):
    """
    Calculate the equivalent orifice diameter (mm) using the Perkins choke model.

    This function attempts to find the orifice diameter that produces a critical
    flow rate matching the target mass flow rate. It uses the bisection method
    to iteratively adjust the orifice size.

    Args:
        inlet_pressure (float): Inlet pressure of the gas.
        inlet_temperature (float): Inlet temperature of the gas.
        mass_flow_rate (float): Target mass flow rate in kg/h.
        upstream_diameter (float): Upstream pipe diameter in meters.
        gas_gravity (float): Specific gravity of the gas.
        gas_cv (float): Specific heat capacity at constant volume of the gas.
        gas_z_factor (float): Compressibility factor of the gas.

    Returns:
        float: Calculated equivalent orifice diameter in mm.
    """
    w_o_target = convert_kg_h_to_lb_s(mass_flow_rate)
    orifice_diameter_start = 0.999 * (upstream_diameter * 1000)  # m to mm.

    choke = PerkinsChokeModel(
        inlet_pressure=inlet_pressure,
        inlet_temperature=inlet_temperature,
        upstream_diameter=upstream_diameter,
        orifice_diameter=orifice_diameter_start,
        gas_gravity=gas_gravity,
        gas_cv=gas_cv,
        gas_z_factor=gas_z_factor,
    )

    def _critical_flow_delta(x):
        """
        Calculate the difference between the critical flow rate and the target flow rate.

        This nested function is used by the bisection method to find the orifice diameter
        that produces the target flow rate under critical flow conditions.

        Args:
            x (float): Current orifice diameter in mm.

        Returns:
            float: Difference between calculated critical flow rate and target flow rate.
        """
        choke.update_orifice_size(x)
        p_r_critical = choke.get_estimate_critical_pressure_ratio()
        w_c_critical = choke.mass_flow_rate(p_r=p_r_critical)
        delta_w = w_c_critical - w_o_target

        return delta_w

    # mm this is a mad number, the correlation really doesn't work beyond 1/64".
    min_orifice = 0.1
    return bisection_method(
        func=_critical_flow_delta, a=orifice_diameter_start, b=min_orifice
    )
