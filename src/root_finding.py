"""Module containing root finding algorithm implementations."""


# bisection for choke size calculation because it works - no need to do initial calcs...
def bisection_method(func, a, b, tol=1e-7, max_iterations=100):
    """
    Find the root of the function `func` in the interval [a, b] using the bisection method.

    Args:
        func (callable): The function for which to find the root.
        a (float): The starting point of the interval.
        b (float): The ending point of the interval.
        tol (float, optional): The tolerance for the stopping criterion. Defaults to 1e-7.
        max_iter (int, optional): The maximum number of iterations. Defaults to 100.

    Returns:
        float: The approximate root of the function within the given tolerance.

    Raises:
        ValueError: If the function values at a and b have the same sign.
    """
    if func(a) * func(b) >= 0:
        raise ValueError(
            "The function must have different signs at the endpoints a and b."
        )

    iter_count = 0
    while (b - a) / 2 > tol and iter_count < max_iterations:
        c = (a + b) / 2  # midpoint.
        if func(c) == 0:
            return c  # exact root found.
        elif func(a) * func(c) < 0:
            b = c  # root is in the left subinterval.
        else:
            a = c  # root is in the right subinterval.
        iter_count += 1

    return (a + b) / 2


def secant(func, x0, x1, tol, max_iterations=100):
    """
    Implements the secant method to find the root of a function.

    The secant method is a root-finding algorithm that uses a succession of roots
    of secant lines to better approximate a root of a function f.

    Args:
        f (callable): The function for which to find the root.
        x0 (float): First initial guess.
        x1 (float): Second initial guess.
        tol (float): Tolerance for convergence. The function will return when |f(x)| < tol.
        max_iterations (int, optional): Maximum number of iterations. Defaults to 100.

    Returns:
        float: Approximation of the root of the function.

    Raises:
        ValueError: If the method fails to converge within the maximum number of iterations.
    """
    # keep initial values for error reporting
    init_x0 = x0
    init_x1 = x1

    # store y values instead of recomputing them.
    fx0 = func(x0)
    fx1 = func(x1)

    # iterate up to maximum number of times
    for _ in range(max_iterations):
        if abs(fx1) < tol:
            # return if converged.
            return x1

        # do calculation
        x2 = (x0 * fx1 - x1 * fx0) / (fx1 - fx0)

        # shift variables (prepare for next loop)
        x0, x1 = x1, x2
        fx0, fx1 = fx1, func(x2)

    # failed to converge.
    raise ValueError(
        f"call to secant(func={func.__name__}, x0={init_x0}, x1={init_x1}, tol={tol},"
        f" max_iterations={max_iterations}) has failed to converge"
    )
