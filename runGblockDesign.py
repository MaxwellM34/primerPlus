"""CLI entrypoint for gBlock generation and hairpin analysis."""

from __future__ import annotations

import argparse

from src.gblock.gblockBuilder import buildGblock, calculateHairpin
from src.io.jsonIo import readJson, writeJson


def main() -> None:
    parser = argparse.ArgumentParser(description='Build a gBlock from a selected primer candidate.')
    parser.add_argument('--inputPath', default='allTrios.json', help='Input JSON containing trios or legacy data.')
    parser.add_argument('--outputPath', default='gblockResult.json', help='Output JSON path for gBlock payload.')
    parser.add_argument('--candidateIndex', type=int, default=0, help='Candidate index to use.')
    parser.add_argument('--includeProbe', action='store_true', help='Include internal probe and split filler gaps.')
    parser.add_argument('--upstreamLength', type=int, default=30, help='Random upstream flank length.')
    parser.add_argument('--downstreamLength', type=int, default=30, help='Random downstream flank length.')
    parser.add_argument('--randomSeed', type=int, default=67, help='Random seed for deterministic filler sequence.')
    args = parser.parse_args()

    payload = readJson(args.inputPath)
    if not isinstance(payload, dict):
        raise ValueError('Input JSON must be an object payload.')

    result = buildGblock(
        payload,
        candidateIndex=args.candidateIndex,
        includeProbe=args.includeProbe,
        upstreamLength=args.upstreamLength,
        downstreamLength=args.downstreamLength,
        randomSeed=args.randomSeed,
    )

    result['hairpin'] = calculateHairpin(result['gblock'])

    writeJson(args.outputPath, result)

    print(f'Wrote {args.outputPath}')
    print(f"gBlock length: {len(result['gblock'])}")


if __name__ == '__main__':
    main()
