"""Candidate filtering, ranking, and coordinate extraction utilities."""

from __future__ import annotations

from typing import Any


def _extractCoords(entry: Any) -> tuple[int, int] | None:
    if isinstance(entry, dict):
        coords = entry.get('COORDS')
        if isinstance(coords, (list, tuple)) and len(coords) >= 2:
            return int(coords[0]), int(coords[1])

    if isinstance(entry, (list, tuple)) and len(entry) >= 2:
        return int(entry[0]), int(entry[1])

    return None


def getPrimerPositions(data: dict[str, Any], *, index: int = 0) -> dict[str, int] | None:
    """Return positional metadata for left/right/internal oligos at a candidate index."""
    try:
        leftEntry = data['PRIMER_LEFT'][index]
        rightEntry = data['PRIMER_RIGHT'][index]
        internalEntry = data['PRIMER_INTERNAL'][index]
    except (KeyError, IndexError, TypeError):
        return None

    leftCoords = _extractCoords(leftEntry)
    rightCoords = _extractCoords(rightEntry)
    internalCoords = _extractCoords(internalEntry)

    if not leftCoords or not rightCoords or not internalCoords:
        return None

    leftStart, leftLen = leftCoords
    rightStart, rightLen = rightCoords
    internalStart, internalLen = internalCoords

    return {
        'left_start': leftStart,
        'left_len': leftLen,
        'right_start': rightStart,
        'right_len': rightLen,
        'internal_start': internalStart,
        'internal_len': internalLen,
    }


def getInternalSequence(
    data: dict[str, Any],
    sequence: str,
    *,
    internalStart: int | None = None,
    internalLen: int | None = None,
    index: int = 0,
) -> str | None:
    """Resolve the internal probe sequence from explicit fields or sequence slicing."""
    explicit = data.get(f'PRIMER_INTERNAL_{index}_SEQUENCE')
    if isinstance(explicit, str) and explicit:
        return explicit

    if internalStart is None or internalLen is None:
        positions = getPrimerPositions(data, index=index)
        if not positions:
            return None
        internalStart = positions['internal_start']
        internalLen = positions['internal_len']

    if internalStart < 0 or internalLen <= 0:
        return None

    end = internalStart + internalLen
    if end > len(sequence):
        return None

    return sequence[internalStart:end]


def passesKrakenSpecificRules(
    data: dict[str, Any],
    sequence: str,
    rules: dict[str, Any],
    *,
    index: int = 0,
) -> bool:
    """Apply spacing and sequence-prefix constraints required by the Kraken workflow."""
    if not rules.get('enabled', False):
        return True

    positions = getPrimerPositions(data, index=index)
    if not positions:
        return False

    leftEnd = positions['left_start'] + positions['left_len'] - 1
    internalStart = positions['internal_start']
    internalEnd = internalStart + positions['internal_len'] - 1
    rightThreePrime = positions['right_start']
    rightFivePrime = rightThreePrime - positions['right_len'] + 1

    leftGap = internalStart - leftEnd - 1
    rightGap = rightFivePrime - internalEnd - 1
    if leftGap < rules['min_left_probe_gap']:
        return False
    if rightGap < rules['min_right_probe_gap']:
        return False

    probePrefixLen = rules['probe_no_g_prefix_len']
    internalEntry = None
    internals = data.get('PRIMER_INTERNAL', [])
    if isinstance(internals, list) and index < len(internals):
        internalEntry = internals[index]

    internalSequence = None
    if isinstance(internalEntry, dict):
        internalSequence = internalEntry.get('SEQUENCE')

    if not internalSequence:
        internalSequence = getInternalSequence(
            data,
            sequence,
            internalStart=internalStart,
            internalLen=positions['internal_len'],
            index=index,
        )

    if not internalSequence or len(internalSequence) < probePrefixLen:
        return False

    if 'G' in internalSequence[:probePrefixLen].upper():
        return False

    return True


def _readPenalty(pair: Any, data: dict[str, Any], index: int) -> float:
    if isinstance(pair, dict) and pair.get('PENALTY') is not None:
        return float(pair['PENALTY'])

    legacyPenalty = data.get(f'PRIMER_PAIR_{index}_PENALTY')
    if legacyPenalty is not None:
        return float(legacyPenalty)

    return float('inf')


def pickBestPairIndex(data: dict[str, Any], sequence: str, rules: dict[str, Any]) -> int | None:
    """Select the lowest-penalty trio that also passes Kraken-specific rules."""
    ranked = rankPairsByPenalty(data, sequence, rules)
    return ranked[0][1] if ranked else None


def rankPairsByPenalty(
    data: dict[str, Any],
    sequence: str,
    rules: dict[str, Any],
    *,
    limit: int | None = None,
) -> list[tuple[float, int]]:
    """Return `(penalty, index)` tuples sorted from best to worst candidate."""
    pairs = data.get('PRIMER_PAIR', [])
    lefts = data.get('PRIMER_LEFT', [])
    rights = data.get('PRIMER_RIGHT', [])
    internals = data.get('PRIMER_INTERNAL', [])

    if not all(isinstance(item, list) for item in (pairs, lefts, rights, internals)):
        return []

    candidateCount = min(len(pairs), len(lefts), len(rights), len(internals))
    ranked: list[tuple[float, int]] = []

    for index in range(candidateCount):
        if not passesKrakenSpecificRules(data, sequence, rules, index=index):
            continue
        penalty = _readPenalty(pairs[index], data, index)
        ranked.append((penalty, index))

    ranked.sort(key=lambda item: (item[0], item[1]))
    if limit is not None:
        return ranked[:limit]
    return ranked
