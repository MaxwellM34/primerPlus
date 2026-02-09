"""CLI entrypoint for iterative primer design."""

from __future__ import annotations

import argparse

from src.pipeline.primerDesignPipeline import runPrimerDesign


def main() -> None:
    parser = argparse.ArgumentParser(description='Run iterative Primer3 design with tunable relaxation.')
    parser.add_argument('--sequencePath', default='assets/templateSequence.fasta', help='Input FASTA path.')
    parser.add_argument('--outputPath', default='allTrios.json', help='Output JSON path.')
    parser.add_argument('--sequenceId', default='MH1000', help='Sequence identifier passed to Primer3.')
    parser.add_argument(
        '--noEntryPrint',
        action='store_true',
        help='Disable detailed printout of the best left/right/internal/pair entries.',
    )
    args = parser.parse_args()

    runPrimerDesign(
        sequencePath=args.sequencePath,
        outputPath=args.outputPath,
        sequenceId=args.sequenceId,
        printBestEntries=not args.noEntryPrint,
    )


if __name__ == '__main__':
    main()
