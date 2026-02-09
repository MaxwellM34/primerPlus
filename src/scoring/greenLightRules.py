"""Rule-based pass/fail checks for TaqMan-style qPCR trios."""

from __future__ import annotations

from typing import Any


def isGreenLightTaqman(trio: dict[str, Any]) -> bool:
    """Return True when a trio meets the conservative GreenLight criteria."""
    left = trio.get('left', {})
    right = trio.get('right', {})
    internal = trio.get('internal', {})
    pair = trio.get('pair', {})

    return bool(
        (70 <= pair.get('PRODUCT_SIZE', -1) <= 150)
        and (80.0 <= pair.get('PRODUCT_TM', -1.0) <= 90.0)
        and (57.0 <= left.get('TM', -1.0) <= 60.5)
        and (57.0 <= right.get('TM', -1.0) <= 60.5)
        and (abs(left.get('TM', 0.0) - right.get('TM', 0.0)) <= 1.5)
        and (35.0 <= left.get('GC_PERCENT', -1.0) <= 60.0)
        and (35.0 <= right.get('GC_PERCENT', -1.0) <= 60.0)
        and (25.0 <= left.get('BOUND', -1.0) <= 40.0)
        and (25.0 <= right.get('BOUND', -1.0) <= 40.0)
        and (2.0 <= left.get('END_STABILITY', -1.0) <= 6.0)
        and (2.0 <= right.get('END_STABILITY', -1.0) <= 6.0)
        and (left.get('SELF_ANY_TH', 999.0) <= 35.0)
        and (left.get('SELF_END_TH', 999.0) <= 35.0)
        and (left.get('HAIRPIN_TH', 999.0) <= 35.0)
        and (right.get('SELF_ANY_TH', 999.0) <= 35.0)
        and (right.get('SELF_END_TH', 999.0) <= 35.0)
        and (right.get('HAIRPIN_TH', 999.0) <= 35.0)
        and (pair.get('COMPL_ANY_TH', 999.0) <= 35.0)
        and (pair.get('COMPL_END_TH', 999.0) <= 35.0)
        and (65.0 <= internal.get('TM', -1.0) <= 75.0)
        and (internal.get('TM', -999.0) >= max(left.get('TM', 0.0), right.get('TM', 0.0)) + 5.0)
        and (internal.get('TM', -999.0) <= max(left.get('TM', 0.0), right.get('TM', 0.0)) + 15.0)
        and (40.0 <= internal.get('GC_PERCENT', -1.0) <= 65.0)
        and (60.0 <= internal.get('BOUND', -1.0) <= 85.0)
        and (internal.get('SELF_ANY_TH', 999.0) <= 35.0)
        and (internal.get('SELF_END_TH', 999.0) <= 35.0)
        and (internal.get('HAIRPIN_TH', 999.0) <= 35.0)
    )
