"""JSON file helpers shared by scoring, pipeline, and gBlock tools."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def readJson(filePath: str | Path) -> Any:
    """Read and decode a JSON file."""
    path = Path(filePath)
    try:
        with path.open('r', encoding='utf-8') as handle:
            return json.load(handle)
    except FileNotFoundError as error:
        raise FileNotFoundError(f'JSON file not found: {path}') from error
    except json.JSONDecodeError as error:
        raise ValueError(f'Invalid JSON in file: {path}') from error


def writeJson(filePath: str | Path, payload: Any, *, indent: int = 2) -> None:
    """Encode and write JSON payload to disk."""
    path = Path(filePath)
    with path.open('w', encoding='utf-8') as handle:
        json.dump(payload, handle, indent=indent)
        handle.write('\n')
