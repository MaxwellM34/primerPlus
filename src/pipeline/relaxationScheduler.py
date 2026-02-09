"""Probabilistic tunable-relaxation scheduler.

Lower-importance tunables are relaxed more often, but all tunables can be
selected until they exhaust their level list.
"""

from __future__ import annotations

import random
from collections.abc import Callable, Iterator
from typing import Any

_rng = random.Random(0)


def buildImportanceWeights(multiplierFn: Callable[[int], float], maxImportance: int) -> dict[int, float]:
    """Convert an importance multiplier into normalized-like descending weights."""
    weights = {1: 1.0}
    for importance in range(1, maxImportance):
        weights[importance + 1] = weights[importance] / float(multiplierFn(importance))
    return weights


def createRelaxationScheduler(
    tunables: dict[str, dict[str, Any]],
    levelMap: dict[str, int],
    *,
    multiplierFn: Callable[[int], float],
    maxImportance: int,
) -> Iterator[tuple[str, int]]:
    """Yield `(tunableKey, nextLevelIndex)` steps until no further relaxation exists."""
    weights = buildImportanceWeights(multiplierFn, maxImportance)

    groups: dict[int, list[str]] = {}
    for key, spec in tunables.items():
        groups.setdefault(spec['importance'], []).append(key)

    for importance in groups:
        groups[importance].sort()

    cursors = {importance: 0 for importance in groups}

    def hasRelaxable(importance: int) -> bool:
        return any(levelMap[key] < len(tunables[key]['levels']) - 1 for key in groups[importance])

    def nextKeyInGroup(importance: int) -> str | None:
        keys = groups[importance]
        if not keys:
            return None

        start = cursors[importance]
        total = len(keys)
        for offset in range(total):
            index = (start + offset) % total
            key = keys[index]
            if levelMap[key] < len(tunables[key]['levels']) - 1:
                cursors[importance] = (index + 1) % total
                return key
        return None

    while True:
        activeGroups = [(importance, weights[importance]) for importance in sorted(groups) if hasRelaxable(importance)]
        if not activeGroups:
            return

        totalWeight = sum(weight for _, weight in activeGroups)
        roll = _rng.random() * totalWeight
        chosenImportance = activeGroups[-1][0]

        for importance, weight in activeGroups:
            roll -= weight
            if roll <= 0:
                chosenImportance = importance
                break

        key = nextKeyInGroup(chosenImportance)
        if key is None:
            continue

        yield key, levelMap[key] + 1


def formatRelaxationSummary(tunables: dict[str, dict[str, Any]], levelMap: dict[str, int]) -> list[str]:
    """Render a stable text summary of active tunable levels."""
    lines: list[str] = []
    for key in sorted(tunables):
        levelIndex = levelMap[key]
        levelValue = tunables[key]['levels'][levelIndex]
        lines.append(f'- {key}: level {levelIndex} -> {levelValue}')
    return lines
