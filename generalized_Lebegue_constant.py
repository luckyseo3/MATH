import sympy as sp


def generalized_Lebegue_function(basis):
    """Return the generalized Lebesgue function as a SymPy expression in x."""
    return sp.Add(*(sp.Abs(f) for f in basis))
