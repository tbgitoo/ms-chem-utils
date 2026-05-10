from __future__ import annotations
from collections import Counter
import re

_FORM_RE = re.compile(r"([A-Z][a-z]?)(\d*)")

def formula_to_counter(formula: str) -> Counter:
    """Parse formula like 'C2H6O' into Counter({'C':2,'H':6,'O':1}). No parentheses support."""
    formula = (formula or "").strip()
    if not formula:
        return Counter()
    out = Counter()
    for elem, num in _FORM_RE.findall(formula):
        out[elem] += int(num) if num else 1
    return out

def counter_to_formula(c: Counter) -> str:
    """Format Counter back to string, with C then H then alphabetical."""
    parts = []
    if c.get("C", 0) > 0:
        parts.append(("C", c["C"]))
    if c.get("H", 0) > 0:
        parts.append(("H", c["H"]))
    for elem in sorted(k for k in c.keys() if k not in ("C", "H")):
        if c.get(elem, 0) > 0:
            parts.append((elem, c[elem]))
    return "".join([f"{e}{n if n != 1 else ''}" for e, n in parts])