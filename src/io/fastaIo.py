"""FASTA file helpers used across design and gBlock workflows."""

from __future__ import annotations

from pathlib import Path


def readFasta(filePath: str | Path) -> str:
    """Read a FASTA file and return the concatenated sequence lines."""
    path = Path(filePath)
    lines = path.read_text(encoding='utf-8').splitlines()

    sequenceParts: list[str] = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith('>'):
            continue
        sequenceParts.append(stripped)

    sequence = ''.join(sequenceParts)
    if not sequence:
        raise ValueError(f'No sequence data found in {path}')
    return sequence


def writeFasta(filePath: str | Path, sequence: str, *, header: str = 'sequence') -> None:
    """Write a sequence to FASTA with 80-character line wrapping."""
    path = Path(filePath)
    wrapped = [sequence[index:index + 80] for index in range(0, len(sequence), 80)]
    content = '\n'.join([f'>{header}', *wrapped]) + '\n'
    path.write_text(content, encoding='utf-8')
