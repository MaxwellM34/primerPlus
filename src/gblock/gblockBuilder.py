"""gBlock construction utilities for primer/probe candidates."""

from __future__ import annotations

import random
from typing import Any

from primer3 import bindings


_bases = ('A', 'C', 'G', 'T')


def reverseComplement(sequence: str) -> str:
    """Return reverse-complement DNA sequence."""
    complementTable = str.maketrans('ATCGatcg', 'TAGCtagc')
    return sequence.translate(complementTable)[::-1]


def _resolveLegacyCandidate(payload: dict[str, Any], candidateIndex: int) -> dict[str, Any]:
    data = payload.get('data', {})
    if not isinstance(data, dict):
        raise ValueError('Legacy payload must include a "data" object.')

    try:
        leftSeq = data[f'PRIMER_LEFT_{candidateIndex}_SEQUENCE']
        rightSeq = data[f'PRIMER_RIGHT_{candidateIndex}_SEQUENCE']
        internalSeq = data.get(f'PRIMER_INTERNAL_{candidateIndex}_SEQUENCE')

        leftCoords = data[f'PRIMER_LEFT_{candidateIndex}']
        rightCoords = data[f'PRIMER_RIGHT_{candidateIndex}']
        internalCoords = data.get(f'PRIMER_INTERNAL_{candidateIndex}')
    except KeyError as error:
        raise ValueError(f'Missing legacy key: {error}') from error

    return {
        'left': {'SEQUENCE': leftSeq, 'COORDS': leftCoords},
        'right': {'SEQUENCE': rightSeq, 'COORDS': rightCoords},
        'internal': {'SEQUENCE': internalSeq, 'COORDS': internalCoords},
    }


def _resolveTrioCandidate(payload: dict[str, Any], candidateIndex: int) -> dict[str, Any]:
    trios = payload.get('trios')
    if not isinstance(trios, list) or not trios:
        raise ValueError('Trio payload must include a non-empty "trios" list.')

    byIndex = next(
        (
            trio
            for trio in trios
            if isinstance(trio, dict) and trio.get('index') == candidateIndex
        ),
        None,
    )

    if byIndex is not None:
        candidate = byIndex
    elif 0 <= candidateIndex < len(trios) and isinstance(trios[candidateIndex], dict):
        candidate = trios[candidateIndex]
    else:
        raise ValueError(f'Candidate index {candidateIndex} not found in trios payload.')

    left = candidate.get('left')
    right = candidate.get('right')
    internal = candidate.get('internal')

    if not all(isinstance(section, dict) for section in (left, right, internal)):
        raise ValueError('Trio candidate must include left/right/internal sections.')

    return {
        'left': left,
        'right': right,
        'internal': internal,
    }


def _resolveCandidate(payload: dict[str, Any], candidateIndex: int) -> dict[str, Any]:
    if 'trios' in payload:
        return _resolveTrioCandidate(payload, candidateIndex)
    if 'data' in payload:
        return _resolveLegacyCandidate(payload, candidateIndex)
    raise ValueError('Unsupported payload format. Expected "trios" or "data".')


def _coords(section: dict[str, Any]) -> tuple[int, int]:
    coords = section.get('COORDS')
    if not isinstance(coords, (list, tuple)) or len(coords) < 2:
        raise ValueError('Missing COORDS values in candidate section.')
    return int(coords[0]), int(coords[1])


def _randomDna(length: int, randomizer: random.Random) -> str:
    if length <= 0:
        return ''
    return ''.join(randomizer.choices(_bases, k=length))


def _fillerCounts(candidate: dict[str, Any], *, includeProbe: bool) -> tuple[int, int, int]:
    leftStart, leftLen = _coords(candidate['left'])
    rightStart, rightLen = _coords(candidate['right'])

    leftEnd = leftStart + leftLen
    rightFivePrime = rightStart - rightLen + 1

    noProbeGap = max(0, rightFivePrime - leftEnd)
    if not includeProbe:
        return 0, 0, noProbeGap

    internalStart, internalLen = _coords(candidate['internal'])
    internalEnd = internalStart + internalLen

    leftProbeGap = max(0, internalStart - leftEnd)
    rightProbeGap = max(0, rightFivePrime - internalEnd)

    return leftProbeGap, rightProbeGap, noProbeGap


def buildGblock(
    payload: dict[str, Any],
    *,
    candidateIndex: int = 0,
    includeProbe: bool = True,
    upstreamLength: int = 30,
    downstreamLength: int = 30,
    randomSeed: int = 67,
) -> dict[str, Any]:
    """Build a synthetic gBlock sequence around a selected primer candidate."""
    randomizer = random.Random(randomSeed)
    candidate = _resolveCandidate(payload, candidateIndex)

    leftSequence = candidate['left'].get('SEQUENCE')
    rightSequence = candidate['right'].get('SEQUENCE')
    internalSequence = candidate['internal'].get('SEQUENCE')

    if not isinstance(leftSequence, str) or not leftSequence:
        raise ValueError('Candidate is missing left primer sequence.')
    if not isinstance(rightSequence, str) or not rightSequence:
        raise ValueError('Candidate is missing right primer sequence.')
    if includeProbe and (not isinstance(internalSequence, str) or not internalSequence):
        raise ValueError('Probe-inclusive gBlock requires an internal probe sequence.')

    leftProbeGap, rightProbeGap, noProbeGap = _fillerCounts(candidate, includeProbe=includeProbe)

    upstream = _randomDna(upstreamLength, randomizer)
    downstream = _randomDna(downstreamLength, randomizer)
    rightOnTemplate = reverseComplement(rightSequence)

    if includeProbe:
        leftFiller = _randomDna(leftProbeGap, randomizer)
        rightFiller = _randomDna(rightProbeGap, randomizer)
        insert = leftSequence + leftFiller + internalSequence + rightFiller + rightOnTemplate
    else:
        middleFiller = _randomDna(noProbeGap, randomizer)
        insert = leftSequence + middleFiller + rightOnTemplate

    gblockSequence = upstream + insert + downstream

    return {
        'candidateIndex': candidateIndex,
        'includeProbe': includeProbe,
        'randomSeed': randomSeed,
        'gblock': gblockSequence,
        'primers': {
            'left': leftSequence,
            'right': rightSequence,
            'internal': internalSequence,
        },
        'fillerLengths': {
            'leftProbeGap': leftProbeGap,
            'rightProbeGap': rightProbeGap,
            'noProbeGap': noProbeGap,
            'upstream': upstreamLength,
            'downstream': downstreamLength,
        },
    }


def _hairpinResultToDict(result: Any) -> dict[str, Any]:
    return {
        'tm': getattr(result, 'tm', None),
        'dg': getattr(result, 'dg', None),
        'dh': getattr(result, 'dh', None),
        'ds': getattr(result, 'ds', None),
        'structureFound': getattr(result, 'structure_found', None),
    }


def calculateHairpin(sequence: str, *, maxSequenceLength: int = 60) -> dict[str, Any]:
    """Run Primer3 hairpin analysis.

    Primer3 thermodynamic calculations require at least one sequence <= 60 bp.
    For longer sequences we scan fixed-size windows and report the highest-Tm
    window result.
    """
    if len(sequence) <= maxSequenceLength:
        return {
            'mode': 'full',
            'sequenceLength': len(sequence),
            'windowLength': len(sequence),
            'bestWindow': {
                'start': 0,
                'sequence': sequence,
                **_hairpinResultToDict(bindings.calc_hairpin(sequence)),
            },
        }

    bestWindow: dict[str, Any] | None = None
    bestTm = float('-inf')

    for start in range(0, len(sequence) - maxSequenceLength + 1):
        windowSequence = sequence[start:start + maxSequenceLength]
        thermo = bindings.calc_hairpin(windowSequence)
        result = _hairpinResultToDict(thermo)
        tm = result.get('tm')
        tmValue = float(tm) if isinstance(tm, (int, float)) else float('-inf')
        if tmValue > bestTm:
            bestTm = tmValue
            bestWindow = {
                'start': start,
                'sequence': windowSequence,
                **result,
            }

    return {
        'mode': 'windowed',
        'sequenceLength': len(sequence),
        'windowLength': maxSequenceLength,
        'bestWindow': bestWindow,
    }
