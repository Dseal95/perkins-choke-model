"""Example use of the PerkinsChokeModel()."""

from src.orifice import PerkinsChokeModel

if __name__ == "__main__":
    # inputs:
    m = 54.102  # kg/h
    p1 = 70.39  # barg -- fluid pressure (inlet conditions)
    t1 = 25.7  # deg C -- fluid temperature (inlet conditions)
    pipe_diameter = 0.21584  # m
    choke_diameter = 0.003  # m
    sg = 0.6
    cv_gas = 39.76  # J/mol/K
    z_gas = 0.8743

    # ** perkins choke model **
    choke = PerkinsChokeModel(
        inlet_pressure=p1,
        inlet_temperature=t1,
        upstream_diameter=pipe_diameter,
        orifice_diameter=choke_diameter,
        gas_gravity=sg,
        gas_cv=cv_gas,
        gas_z_factor=z_gas,
    )
    # calculate p2, outlet pressure.
    p2 = choke.calc_perkins_dp(m=m)

    # debug info.
    print("\n\n** Perkins Choke Model Finished **")
    print(f"input params: t1={t1} degC, p1={p1} barg, m_observed = {m} kg/h")
    print(f"p_r (critical pressure ratio) = {choke.p_r_critical}")
    print(f"calculated critical mass flow rate, m_critical = {choke.m_critical} kg/h")
    print(f"calculated outlet pressure, p2 = {p2} barg")
    print(f"m_critical > m_observed (kg/h) = {choke.m_critical > m}\n\n")
