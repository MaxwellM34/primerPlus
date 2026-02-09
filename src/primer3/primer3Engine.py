"""Adapter around primer3-py bindings.

This module centralizes Primer3 request assembly so all pipelines use identical
argument construction and validation.
"""

from __future__ import annotations

from typing import Any

import primer3


def _resolveTunableUpdates(tunableKey: str, spec: dict[str, Any], levelIndex: int) -> dict[str, Any]:
    levels = spec['levels']
    if levelIndex < 0 or levelIndex >= len(levels):
        raise ValueError(f'Invalid level index for {tunableKey}: {levelIndex}')

    level = levels[levelIndex]
    blockType = spec['type']

    if blockType == 'tm_block':
        return {
            'PRIMER_MIN_TM': level['MIN'],
            'PRIMER_OPT_TM': level['OPT'],
            'PRIMER_MAX_TM': level['MAX'],
            'PRIMER_MAX_TM_DIFF': level['DIFF'],
        }
    if blockType == 'size_block':
        return {
            'PRIMER_MIN_SIZE': level['MIN'],
            'PRIMER_OPT_SIZE': level['OPT'],
            'PRIMER_MAX_SIZE': level['MAX'],
        }
    if blockType == 'product_size_block':
        return {'PRIMER_PRODUCT_SIZE_RANGE': level}
    if blockType == 'gc_block':
        return {
            'PRIMER_MIN_GC': level['MIN'],
            'PRIMER_OPT_GC_PERCENT': level['OPT'],
            'PRIMER_MAX_GC': level['MAX'],
        }
    if blockType == 'internal_tm_block':
        return {
            'PRIMER_INTERNAL_MIN_TM': level['MIN'],
            'PRIMER_INTERNAL_OPT_TM': level['OPT'],
            'PRIMER_INTERNAL_MAX_TM': level['MAX'],
        }
    if blockType == 'internal_size_block':
        return {
            'PRIMER_INTERNAL_MIN_SIZE': level['MIN'],
            'PRIMER_INTERNAL_OPT_SIZE': level['OPT'],
            'PRIMER_INTERNAL_MAX_SIZE': level['MAX'],
        }
    if blockType == 'internal_gc_block':
        return {
            'PRIMER_INTERNAL_MIN_GC': level['MIN'],
            'PRIMER_INTERNAL_OPT_GC_PERCENT': level['OPT'],
            'PRIMER_INTERNAL_MAX_GC': level['MAX'],
        }
    if blockType == 'leq':
        return {spec['primer3_key']: level}

    raise ValueError(f'Unknown tunable type for {tunableKey}: {blockType}')


def buildGlobalArgs(
    constants: dict[str, Any],
    tunables: dict[str, dict[str, Any]],
    levelMap: dict[str, int],
) -> dict[str, Any]:
    """Merge constants with the current tunable level selections."""
    globalArgs = dict(constants)

    for tunableKey, spec in tunables.items():
        levelIndex = levelMap.get(tunableKey, 0)
        updates = _resolveTunableUpdates(tunableKey, spec, levelIndex)

        for key, value in updates.items():
            if key in constants:
                continue
            if key in globalArgs:
                raise ValueError(f'Duplicate global arg key from tunables: {key}')
            globalArgs[key] = value

    return globalArgs


def runPrimer3(seqArgs: dict[str, Any], globalArgs: dict[str, Any]) -> dict[str, Any]:
    """Call Primer3 and return the parsed dictionary response."""
    return primer3.bindings.design_primers(seq_args=seqArgs, global_args=globalArgs)


def isSuccessfulResult(data: object) -> bool:
    """Return True when Primer3 produced at least one full left/right/internal trio."""
    if not isinstance(data, dict):
        return False

    requiredKeys = ('PRIMER_PAIR', 'PRIMER_LEFT', 'PRIMER_RIGHT', 'PRIMER_INTERNAL')
    for key in requiredKeys:
        value = data.get(key)
        if not isinstance(value, list) or not value:
            return False

    return True
