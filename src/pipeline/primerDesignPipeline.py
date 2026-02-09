"""Iterative Primer3 design pipeline with tunable-relaxation fallback."""

from __future__ import annotations

import itertools
from typing import Any

from ..config.designConfig import (
    krakenSpecificRules,
    maxImportance,
    multiplierFn,
    primer3Constants,
    tunableRules,
)
from ..io.fastaIo import readFasta
from ..io.jsonIo import writeJson
from ..primer3.primer3Engine import buildGlobalArgs, isSuccessfulResult, runPrimer3
from ..selection.primerSelection import (
    getPrimerPositions,
    passesKrakenSpecificRules,
    pickBestPairIndex,
)
from ..selection.sequenceView import formatSequenceWithHighlights, printSequenceWithHighlights
from .relaxationScheduler import createRelaxationScheduler, formatRelaxationSummary


def formatEntryLines(label: str, entry: Any, *, maxItems: int = 50) -> list[str]:
    """Format a dict/list object into aligned console-friendly text lines."""
    lines = [f'{label}:']

    if isinstance(entry, dict):
        items = list(itertools.islice(entry.items(), maxItems))
    elif isinstance(entry, (list, tuple)):
        items = list(enumerate(entry[:maxItems]))
    else:
        items = [(label, entry)]

    keyWidth = max((len(str(key)) for key, _ in items), default=0)
    for key, value in items:
        valueText = '' if value is None else str(value)
        valueLines = valueText.splitlines() or ['']

        prefix = f'  {str(key).ljust(keyWidth)}: '
        lines.append(prefix + valueLines[0])

        if len(valueLines) > 1:
            indent = '  ' + ' ' * keyWidth + '  '
            for continuation in valueLines[1:]:
                lines.append(indent + continuation)

    lines.append('')
    return lines


def printEntry(label: str, entry: Any, *, maxItems: int = 10) -> None:
    """Print a compact object preview for a primer entry."""
    for line in formatEntryLines(label, entry, maxItems=maxItems):
        print(line)


def _buildTrioRecord(
    data: dict[str, Any],
    sequence: str,
    *,
    index: int,
    firstBaseIndex: int,
) -> dict[str, Any]:
    trio = {
        'index': index,
        'pair': data['PRIMER_PAIR'][index],
        'left': data['PRIMER_LEFT'][index],
        'right': data['PRIMER_RIGHT'][index],
        'internal': data['PRIMER_INTERNAL'][index],
    }

    positions = getPrimerPositions(data, index=index)
    if positions:
        trio['positions'] = positions
        trio['sequenceView'] = {
            'format': 'html',
            'lineWidth': 50,
            'groupSize': 10,
            'firstBaseIndex': firstBaseIndex,
            'lines': formatSequenceWithHighlights(
                sequence,
                positions,
                firstBaseIndex=firstBaseIndex,
                outputFormat='html',
            ),
        }

    return trio


def runPrimerDesign(
    *,
    sequencePath: str = 'assets/templateSequence.fasta',
    outputPath: str = 'allTrios.json',
    sequenceId: str = 'MH1000',
    printBestEntries: bool = True,
) -> dict[str, Any]:
    """Run iterative Primer3 search and write all valid trios to JSON."""
    sequence = readFasta(sequencePath)

    levelMap = {key: 0 for key in tunableRules}
    scheduler = createRelaxationScheduler(
        tunableRules,
        levelMap,
        multiplierFn=multiplierFn,
        maxImportance=maxImportance,
    )

    seqArgs = {
        'SEQUENCE_ID': sequenceId,
        'SEQUENCE_TEMPLATE': sequence,
    }

    attempt = 0
    relaxedKey = 'none'
    relaxedLevel: int | None = None

    while True:
        globalArgs = buildGlobalArgs(primer3Constants, tunableRules, levelMap)
        data = runPrimer3(seqArgs, globalArgs)

        pairs = data.get('PRIMER_PAIR', []) if isinstance(data, dict) else []
        pairCount = len(pairs) if isinstance(pairs, list) else 0

        levelDisplay = '-' if relaxedLevel is None else str(relaxedLevel)
        print(f'Attempt {attempt} | relaxed: {relaxedKey} | level: {levelDisplay} | pairs: {pairCount}')

        bestIndex = pickBestPairIndex(data, sequence, krakenSpecificRules)
        if isSuccessfulResult(data) and bestIndex is not None:
            best = {
                'pair': data['PRIMER_PAIR'][bestIndex],
                'left': data['PRIMER_LEFT'][bestIndex],
                'right': data['PRIMER_RIGHT'][bestIndex],
                'internal': data['PRIMER_INTERNAL'][bestIndex],
            }

            if printBestEntries:
                printEntry('LEFT', best['left'])
                printEntry('RIGHT', best['right'])
                printEntry('INTERNAL', best['internal'])
                printEntry('PAIR', best['pair'])
                print('SEQUENCE VIEW:')
                printSequenceWithHighlights(
                    sequence,
                    data,
                    index=bestIndex,
                    firstBaseIndex=primer3Constants.get('PRIMER_FIRST_BASE_INDEX', 1),
                )

            lefts = data.get('PRIMER_LEFT', [])
            rights = data.get('PRIMER_RIGHT', [])
            internals = data.get('PRIMER_INTERNAL', [])
            candidateCount = min(len(pairs), len(lefts), len(rights), len(internals))
            firstBaseIndex = primer3Constants.get('PRIMER_FIRST_BASE_INDEX', 1)

            trios: list[dict[str, Any]] = []
            for index in range(candidateCount):
                if not passesKrakenSpecificRules(data, sequence, krakenSpecificRules, index=index):
                    continue
                trios.append(
                    _buildTrioRecord(
                        data,
                        sequence,
                        index=index,
                        firstBaseIndex=firstBaseIndex,
                    )
                )

            output = {
                'attempt': attempt,
                'relaxedKey': relaxedKey,
                'relaxedLevel': relaxedLevel,
                'pairCount': pairCount,
                'sequenceId': sequenceId,
                'sequence': sequence,
                'bestIndex': bestIndex,
                'levelMap': levelMap,
                'levelSummary': formatRelaxationSummary(tunableRules, levelMap),
                'totalAttempts': attempt + 1,
                'trios': trios,
            }
            writeJson(outputPath, output)
            print(f'Wrote {outputPath}')
            print('Final level map:')
            for line in output['levelSummary']:
                print(line)
            print(f'Total attempts: {output["totalAttempts"]}')
            return output

        try:
            nextKey, nextLevel = next(scheduler)
        except StopIteration:
            output = {
                'attempt': attempt,
                'relaxedKey': relaxedKey,
                'relaxedLevel': relaxedLevel,
                'pairCount': pairCount,
                'sequenceId': sequenceId,
                'sequence': sequence,
                'bestIndex': None,
                'levelMap': levelMap,
                'levelSummary': formatRelaxationSummary(tunableRules, levelMap),
                'totalAttempts': attempt + 1,
                'trios': [],
                'status': 'no_solution',
            }
            writeJson(outputPath, output)
            print('No valid primer pair with internal oligo found after exhausting tunables.')
            print(f'Wrote {outputPath}')
            print(f'Total attempts: {output["totalAttempts"]}')
            return output

        levelMap[nextKey] = nextLevel
        attempt += 1
        relaxedKey = nextKey
        relaxedLevel = nextLevel
