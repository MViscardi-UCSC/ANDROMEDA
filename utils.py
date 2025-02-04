"""
utils.py
Marcus Viscardi,    February 04, 2025

Random helper functions that don't discreetly fit into any other category.
"""

def clamp(n: float | int,
          smallest: float | int,
          largest: float | int) -> float | int:
    """
    Clamp a number between two values.
    """
    return max(smallest, min(n, largest))
