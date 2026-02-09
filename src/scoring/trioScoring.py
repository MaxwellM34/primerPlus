"""Quantile-based scoring for primer/probe trios.

The scorer computes deltas from configurable optimal values, maps each delta to
quantile bins, and combines per-metric scores into weighted group percentages.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ..io.jsonIo import readJson, writeJson


optimals = {
    'pair': {
        'COMPL_ANY_TH': 0.0,
        'COMPL_END_TH': 0.0,
        'PRODUCT_SIZE': 100,
        'PRODUCT_TM': 83,
    },
    'primer': {
        'length': 21,
        'BOUND': 30,
        'GC_PERCENT': 50.0,
        'SELF_ANY_TH': 0.0,
        'SELF_END_TH': 0.0,
        'HAIRPIN_TH': 0.0,
        'END_STABILITY': 3.5,
    },
    'internal': {
        'length': 23,
        'TM_DIFF': 10.5,
        'BOUND': 70,
        'GC_PERCENT': 50.0,
        'SELF_ANY_TH': 0.0,
        'SELF_END_TH': 0.0,
        'HAIRPIN_TH': 0.0,
        'END_STABILITY': 3.5,
    },
}


groupWeights = {
    'pair': 1.0 / 3.0,
    'primers': 1.0 / 3.0,
    'internal': 1.0 / 3.0,
}


quantileSteps = [1 / 6, 2 / 6, 3 / 6, 4 / 6, 5 / 6]


def asFloat(value: object) -> float | None:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def getLength(entry: dict[str, Any]) -> int | None:
    sequence = entry.get('SEQUENCE')
    if isinstance(sequence, str) and sequence:
        return len(sequence)

    coords = entry.get('COORDS')
    if isinstance(coords, (list, tuple)) and len(coords) >= 2:
        coordLength = asFloat(coords[1])
        if coordLength is not None:
            return int(coordLength)

    length = entry.get('length')
    if isinstance(length, (int, float)):
        return int(length)

    return None


def getMetric(entry: dict[str, Any], key: str) -> float | None:
    if key == 'length':
        length = getLength(entry)
        return float(length) if length is not None else None
    return asFloat(entry.get(key))


def numericDelta(value: object, target: object) -> float | None:
    observed = asFloat(value)
    expected = asFloat(target)
    if observed is None or expected is None:
        return None
    return abs(observed - expected)


def tmDiffDelta(
    left: dict[str, Any],
    right: dict[str, Any],
    internal: dict[str, Any],
    target: object,
) -> float | None:
    leftTm = getMetric(left, 'TM')
    rightTm = getMetric(right, 'TM')
    internalTm = getMetric(internal, 'TM')
    targetDiff = asFloat(target)

    if leftTm is None or rightTm is None or internalTm is None or targetDiff is None:
        return None

    primerAverage = (leftTm + rightTm) / 2.0
    return abs((internalTm - primerAverage) - targetDiff)


def computeDeltas(trio: dict[str, Any], optimalValues: dict[str, dict[str, object]]) -> dict[str, dict[str, float | None]]:
    left = trio.get('left', {})
    right = trio.get('right', {})
    internal = trio.get('internal', {})
    pair = trio.get('pair', {})

    deltas: dict[str, dict[str, float | None]] = {
        'pair': {},
        'left': {},
        'right': {},
        'internal': {},
    }

    for key, target in optimalValues['pair'].items():
        deltas['pair'][key] = numericDelta(getMetric(pair, key), target)

    for key, target in optimalValues['primer'].items():
        deltas['left'][key] = numericDelta(getMetric(left, key), target)
        deltas['right'][key] = numericDelta(getMetric(right, key), target)

    for key, target in optimalValues['internal'].items():
        if key == 'TM_DIFF':
            deltas['internal'][key] = tmDiffDelta(left, right, internal, target)
        else:
            deltas['internal'][key] = numericDelta(getMetric(internal, key), target)

    return deltas


def buildMetricList() -> list[tuple[str, str]]:
    metrics: list[tuple[str, str]] = []
    for key in optimals['pair']:
        metrics.append(('pair', key))
    for key in optimals['primer']:
        metrics.append(('left', key))
        metrics.append(('right', key))
    for key in optimals['internal']:
        metrics.append(('internal', key))
    return metrics


def computeBins(values: list[float]) -> np.ndarray | None:
    if not values:
        return None
    return np.quantile(np.array(values), quantileSteps)


def scoreFromBins(value: float | None, bins: np.ndarray | None) -> int:
    if value is None or bins is None or len(bins) == 0:
        return 0
    if value <= bins[0]:
        return 5
    if len(bins) > 1 and value <= bins[1]:
        return 4
    if len(bins) > 2 and value <= bins[2]:
        return 3
    if len(bins) > 3 and value <= bins[3]:
        return 2
    if len(bins) > 4 and value <= bins[4]:
        return 1
    return 0


def buildMetricBins(trios: list[dict[str, Any]], metrics: list[tuple[str, str]]) -> dict[str, np.ndarray]:
    binsByMetric: dict[str, np.ndarray] = {}

    for group, key in metrics:
        values: list[float] = []
        for trio in trios:
            deltas = trio.get('deltas', {})
            groupDeltas = deltas.get(group, {})
            value = groupDeltas.get(key)
            if isinstance(value, (int, float)):
                values.append(float(value))

        bins = computeBins(values)
        if bins is not None:
            binsByMetric[f'{group}.{key}'] = bins

    return binsByMetric


def normalizeWeights(weights: dict[str, float]) -> dict[str, float]:
    totalWeight = sum(weights.values())
    if totalWeight <= 0:
        count = len(weights) if weights else 1
        return {key: 1.0 / count for key in weights}
    return {key: value / totalWeight for key, value in weights.items()}


def averageScores(values: list[float]) -> float:
    if not values:
        return 0.0
    return sum(values) / len(values)


def computeGroupScores(trio: dict[str, Any], binsByMetric: dict[str, np.ndarray]) -> dict[str, float]:
    deltas = trio.get('deltas', {})

    pairScores: list[float] = []
    for key in optimals['pair']:
        metricId = f'pair.{key}'
        pairScores.append(float(scoreFromBins(deltas.get('pair', {}).get(key), binsByMetric.get(metricId))))

    primerScores: list[float] = []
    for key in optimals['primer']:
        leftMetricId = f'left.{key}'
        rightMetricId = f'right.{key}'
        leftScore = scoreFromBins(deltas.get('left', {}).get(key), binsByMetric.get(leftMetricId))
        rightScore = scoreFromBins(deltas.get('right', {}).get(key), binsByMetric.get(rightMetricId))
        primerScores.append((leftScore + rightScore) / 2.0)

    internalScores: list[float] = []
    for key in optimals['internal']:
        metricId = f'internal.{key}'
        internalScores.append(float(scoreFromBins(deltas.get('internal', {}).get(key), binsByMetric.get(metricId))))

    return {
        'pair': averageScores(pairScores),
        'primers': averageScores(primerScores),
        'internal': averageScores(internalScores),
    }


def computeWeightedPercent(scoresByGroup: dict[str, float], weights: dict[str, float]) -> float:
    normalizedWeights = normalizeWeights(weights)

    weightedScore = (
        scoresByGroup.get('pair', 0.0) * normalizedWeights.get('pair', 0.0)
        + scoresByGroup.get('primers', 0.0) * normalizedWeights.get('primers', 0.0)
        + scoresByGroup.get('internal', 0.0) * normalizedWeights.get('internal', 0.0)
    )
    return (weightedScore / 5.0) * 100.0


def scoreTrios(
    trios: list[dict[str, Any]],
    *,
    optimalValues: dict[str, dict[str, object]] = optimals,
    weights: dict[str, float] = groupWeights,
) -> list[dict[str, Any]]:
    """Attach delta metrics and quantile suitability scores to each trio."""
    for trio in trios:
        trio['deltas'] = computeDeltas(trio, optimalValues)

    metricList = buildMetricList()
    binsByMetric = buildMetricBins(trios, metricList)

    for trio in trios:
        scoreByGroup = computeGroupScores(trio, binsByMetric)
        trio['groupScores'] = scoreByGroup
        trio['groupWeights'] = normalizeWeights(weights)
        trio['quantilePercent'] = computeWeightedPercent(scoreByGroup, weights)

    return sorted(trios, key=lambda item: item.get('quantilePercent', 0.0), reverse=True)


def loadAndScoreTrios(
    inputPath: str = 'allTrios.json',
    outputPath: str | None = None,
) -> dict[str, Any]:
    """Load design output JSON, score trios, optionally persist scored payload."""
    payload = readJson(inputPath)
    if not isinstance(payload, dict):
        raise ValueError('Input JSON must be an object with a "trios" array.')

    trios = payload.get('trios')
    if not isinstance(trios, list):
        raise ValueError('Input JSON does not include a valid "trios" list.')

    normalizedTrios = [trio for trio in trios if isinstance(trio, dict)]
    ranked = scoreTrios(normalizedTrios)

    payload['trios'] = ranked
    payload['scoring'] = {
        'method': 'quantile',
        'groupWeights': normalizeWeights(groupWeights),
        'count': len(ranked),
    }

    if outputPath:
        writeJson(outputPath, payload)

    return payload
