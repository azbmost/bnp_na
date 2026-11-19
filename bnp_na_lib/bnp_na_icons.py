"""bnp_na-only additions to the shared AZBMOST icon set.

azbmost_ui.py and azbmost_icons.py are identical copies shared by every
AZBMOST module, so tools that exist only in bnp_na register their glyph and
one-line purpose here at import time instead of editing those files.
"""
from __future__ import annotations

import math

from . import azbmost_icons, azbmost_ui


def g_bp_builder(p, k) -> None:
    """B-DNA (wide, below) joined on one axis to P-DNA above: a narrow
    phosphate core with its bases sticking out."""
    lo, hi, mid = 800, 236, 512

    def amplitude(y: float) -> float:
        # Wide B-DNA strands below, narrowing into the P-DNA core above.
        if y >= 560:
            return 132.0
        if y <= 470:
            return 40.0
        return 40.0 + 92.0 * (y - 470) / 90.0

    def phase(y: float) -> float:
        # Two B turns below; tighter turns above.
        if y >= 515:
            return 2 * math.pi * (lo - y) / (lo - 520)
        return 2 * math.pi * (lo - 515) / (lo - 520) + 2 * math.pi * (515 - y) / 120

    ys = [lo + (hi - lo) * i / 120 for i in range(121)]
    for sign, ink, width in ((-1, k.b, 34), (1, k.a, 46)):
        p.line([(mid + sign * amplitude(y) * math.sin(phase(y)), y) for y in ys], ink, width)
    for y, side in ((430, 1), (372, -1), (314, 1), (256, -1)):
        p.line([(mid + side * 62, y), (mid + side * 170, y - 14)], k.a, 30)


def register() -> None:
    azbmost_icons.TOOLS.setdefault("bp_builder", (g_bp_builder, "orange"))
    azbmost_ui.TOOL_INFO.setdefault(
        "bp_builder",
        ("bnp_na", "Join alternating B- and P-DNA PDBs; the P side's first two levels become the transition."),
    )


register()
