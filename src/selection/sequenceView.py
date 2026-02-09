"""Human-readable sequence visualization helpers."""

from __future__ import annotations

from typing import Any

from .primerSelection import getPrimerPositions


def formatSequenceWithHighlights(
    sequence: str,
    positions: dict[str, int] | None,
    *,
    lineWidth: int = 50,
    groupSize: int = 10,
    firstBaseIndex: int = 1,
    useAnsi: bool = True,
    outputFormat: str | None = None,
) -> list[str]:
    """Build a line-wrapped sequence view with left/internal/right regions highlighted."""
    if not positions:
        return []

    sequenceLength = len(sequence)
    styleMap: list[str | None] = [None] * sequenceLength

    def clampSpan(start: int | None, length: int | None, *, reverse: bool = False) -> tuple[int, int] | None:
        if start is None or length is None or length <= 0:
            return None

        zeroBased = start - firstBaseIndex
        if reverse:
            end = zeroBased + 1
            spanStart = end - length
        else:
            spanStart = zeroBased
            end = spanStart + length

        if spanStart >= sequenceLength or end <= 0:
            return None

        return max(spanStart, 0), min(end, sequenceLength)

    spans = {
        'left': clampSpan(positions['left_start'], positions['left_len']),
        'internal': clampSpan(positions['internal_start'], positions['internal_len']),
        'right': clampSpan(positions['right_start'], positions['right_len'], reverse=True),
    }

    for style, span in spans.items():
        if not span:
            continue
        spanStart, spanEnd = span
        for index in range(spanStart, spanEnd):
            styleMap[index] = style

    if outputFormat is None:
        outputFormat = 'ansi' if useAnsi else 'plain'

    if outputFormat == 'ansi':
        styles = {
            'left': '\033[30;42m',
            'internal': '\033[97;45m',
            'right': '\033[97;44m',
        }
        reset = '\033[0m'
    elif outputFormat == 'html':
        styles = {
            'left': '<span style="background-color:#22c55e;color:#000;">',
            'internal': '<span style="background-color:#a855f7;color:#fff;">',
            'right': '<span style="background-color:#2563eb;color:#fff;">',
        }
        reset = '</span>'
    else:
        styles = {}
        reset = ''

    def applyStyles(text: str, styleSlice: list[str | None]) -> str:
        if outputFormat == 'plain':
            return text

        chunks: list[str] = []
        currentStyle: str | None = None

        for character, style in zip(text, styleSlice):
            if style != currentStyle:
                if currentStyle is not None:
                    chunks.append(reset)
                if style is not None:
                    chunks.append(styles.get(style, ''))
                currentStyle = style
            chunks.append(character)

        if currentStyle is not None:
            chunks.append(reset)

        return ''.join(chunks)

    labelWidth = max(4, len(str(sequenceLength + firstBaseIndex)))
    lines: list[str] = []

    for lineStart in range(0, sequenceLength, lineWidth):
        lineEnd = min(lineStart + lineWidth, sequenceLength)
        lineSequence = sequence[lineStart:lineEnd]
        lineStyles = styleMap[lineStart:lineEnd]

        groups: list[str] = []
        for groupStart in range(0, len(lineSequence), groupSize):
            groupEnd = min(groupStart + groupSize, len(lineSequence))
            groupSequence = lineSequence[groupStart:groupEnd]
            groupStyles = lineStyles[groupStart:groupEnd]
            groups.append(applyStyles(groupSequence, groupStyles))

        label = str(lineStart + firstBaseIndex).rjust(labelWidth)
        lines.append(f'{label}  ' + ' '.join(groups))

    return lines


def printSequenceWithHighlights(
    sequence: str,
    data: dict[str, Any],
    *,
    index: int = 0,
    lineWidth: int = 50,
    groupSize: int = 10,
    firstBaseIndex: int = 1,
    useAnsi: bool = True,
    outputFormat: str | None = None,
) -> bool:
    """Render one candidate's highlighted sequence view to stdout."""
    positions = getPrimerPositions(data, index=index)
    if not positions:
        print('No primer positions available for sequence view.')
        return False

    lines = formatSequenceWithHighlights(
        sequence,
        positions,
        lineWidth=lineWidth,
        groupSize=groupSize,
        firstBaseIndex=firstBaseIndex,
        useAnsi=useAnsi,
        outputFormat=outputFormat,
    )
    for line in lines:
        print(line)
    return True
