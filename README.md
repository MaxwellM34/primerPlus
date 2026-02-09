# primerPlus (Refactored)

This repo is organized into focused modules and camelCase naming.

## Layout

- `src/config/designConfig.py`: Primer3 constants and tunable relaxation config.
- `src/io/`: FASTA and JSON helpers.
- `src/primer3/primer3Engine.py`: Primer3 argument building and execution.
- `src/pipeline/`: Design pipeline and relaxation scheduler.
- `src/selection/`: Candidate filtering and sequence highlighting.
- `src/scoring/`: Quantile scoring and GreenLight rule checks.
- `src/gblock/`: gBlock generation and hairpin analysis.
- `assets/`: input templates and sample data.

## CLI scripts

- `python runPrimerDesign.py`
- `python runTrioScoring.py`
- `python runGblockDesign.py --includeProbe`

Each script supports `--help` for full options.
