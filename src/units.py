"""units.py module contains functionality to convert units."""


def convert_barg_to_bar(p: float):
    """Convert input pressure from barg to bar.

    Args:
        p (float): Input pressure in units barg.

    Returns:
        (float)): Pressure in bar.
    """
    # barg to bar
    p_bar = p + 1.01325

    return p_bar


def convert_bar_to_kpa(p: float):
    """Convert input pressure from bar to kPa.

    Args:
        p (float): Input pressure in units bar.

    Returns:
        (float)): Pressure in kPa.
    """
    # bar to Kpa
    p_kpa = p * 100

    return p_kpa


def convert_kpa_to_bar(p: float):
    """Convert input pressure from kPa to bar.

    Args:
        p (float): Input pressure in units Kpa.

    Returns:
        (float)): Pressure in bar.
    """
    # Kpa to bar
    p_bar = p / 100

    return p_bar


def convert_bar_to_barg(p: float):
    """Convert input pressure from bar to barg.

    Args:
        p (float): Input pressure in units bar.

    Returns:
        (float)): Pressure in barg.
    """
    # bar to barg
    p_barg = p - 1.01325

    return p_barg


def convert_degc_to_k(t: float):
    """Convert input temperature from degC to K.

    Args:
        t (float): Input temperature in degC.

    Returns:
        float: Temperature in K.
    """
    return t + 273.15


def convert_k_to_degc(t: float):
    """Convert input temperature from K to degC.

    Args:
        t (float): Input temperature in K.

    Returns:
        float: Temperature in degC.
    """
    return t - 273.15


def convert_kg_h_to_kg_s(M: float):
    """Convert mass flow from kg/h to kg/s.

    Args:
        M (float): Input mass flow in kg/h.

    Returns:
        float: Converted mass flow in kg/s.
    """
    return M / 3600


def convert_barg_to_psia(p: float):
    """Convert pressure from barg to psia.

    Args:
        p (float): Input pressure in barg.

    Returns:
        float: Converted pressure in psia.
    """
    return (p * 14.5038) + 14.696


def convert_psia_to_barg(p: float):
    """Convert pressure from psia to barg.

    Args:
        p (float): Input pressure in psia.

    Returns:
        float: Converted pressure in barg.
    """
    return (p - 14.696) / 14.5038


def convert_kg_h_to_lb_s(m: float):
    """Convert mass flow from kg/h to lbs/s.

    Args:
        m (float): Input mass flow in kg/h.

    Returns:
        float: Converted mass flow in lbs/s.
    """
    return (2.20462 * m) / 3600


def convert_lb_s_to_kg_h(m: float):
    """Convert mass flow from lbs/s to kg/h.

    Args:
        m (float): Input mass flow in lbs/s.

    Returns:
        float: Converted mass flow in kg/h.
    """
    return (m * 3600) / 2.20462


def convert_deg_c_to_deg_f(t: float) -> float:
    """Convert temperature from Celsius to Fahrenheit.

    Args:
        t (float): Temperature in Celsius.

    Returns:
        float: Temperature in Fahrenheit.
    """
    return ((t * 9.0) / 5.0) + 32.0


def convert_deg_f_to_deg_c(t: float) -> float:
    """Convert temperature from Fahrenheit to Celsius.

    Args:
        t (float): Temperature in Fahrenheit.

    Returns:
        float: Temperature in Celsius.
    """
    return ((t - 32.0) * 5.0) / 9.0


def convert_deg_f_to_deg_rankine(t: float) -> float:
    """Convert temperature from Fahrenheit to Rankine.

    Args:
        t (float): Temperature in Fahrenheit.

    Returns:
        float: Temperature in Rankine.
    """
    return t + 459.67


def convert_deg_rankine_to_deg_f(t: float) -> float:
    """Convert temperature from Rankine to Fahrenheit.

    Args:
        t (float): Temperature in Rankine.

    Returns:
        float: Temperature in Fahrenheit.
    """
    return t - 459.67


def convert_m_to_ft(d: float) -> float:
    """Convert length from meters to feet.

    Args:
        d (float): Length in meters.

    Returns:
        float: Length in feet.
    """
    return d * 3.28084


def convert_ft_to_m(d: float) -> float:
    """Convert length from feet to meters.

    Args:
        d (float): Length in feet.

    Returns:
        float: Length in meters.
    """
    return d / 3.28084


def convert_knots_to_m_s(v):
    """
    Convert speed from knots to meters per second.

    This function takes a speed value in knots and returns the equivalent speed in meters per second.

    Args:
        v (float): The speed in knots to be converted.

    Returns:
        float: The equivalent speed in meters per second.
    """

    return v / 1.94
