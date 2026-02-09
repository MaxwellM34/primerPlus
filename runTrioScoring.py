"""CLI entrypoint for trio quantile scoring."""

from __future__ import annotations

import argparse

from src.scoring.trioScoring import loadAndScoreTrios


def main() -> None:
    parser = argparse.ArgumentParser(description='Score primer/probe trios with quantile-based suitability.')
    parser.add_argument('--inputPath', default='allTrios.json', help='Design output JSON to score.')
    parser.add_argument('--outputPath', default='scoredTrios.json', help='Where to write scored JSON.')
    parser.add_argument('--topN', type=int, default=5, help='How many top trios to print.')
    args = parser.parse_args()

    payload = loadAndScoreTrios(inputPath=args.inputPath, outputPath=args.outputPath)
    trios = payload.get('trios', []) if isinstance(payload, dict) else []

    print(f'Computed quantile scores for {len(trios)} trios.')
    for trio in trios[: max(0, args.topN)]:
        index = trio.get('index', 'n/a')
        percent = trio.get('quantilePercent', 0.0)
        print(f'{index}: suitability={percent:.2f}/100')

    print(f'Wrote {args.outputPath}')


if __name__ == '__main__':
    main()
