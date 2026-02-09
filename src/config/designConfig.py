"""Centralized configuration for PrimerPlus design and filtering.

This module intentionally keeps all tunable values in one place so design policy
can be updated without changing pipeline code.
"""

from __future__ import annotations


primer3Constants = {
    'PRIMER_PICK_LEFT_PRIMER': 1,
    'PRIMER_PICK_RIGHT_PRIMER': 1,
    'PRIMER_PICK_INTERNAL_OLIGO': 1,
    'PRIMER_NUM_RETURN': 5000,
    'PRIMER_INTERNAL_NUM_RETURN': 5000,
    'PRIMER_ANNEALING_TEMP': 61,
    'PRIMER_GC_CLAMP': 1,
    'PRIMER_MAX_POLY_X': 5,
    'PRIMER_MAX_NS_ACCEPTED': 1,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_SALT_DIVALENT': 1.5,
    'PRIMER_DNTP_CONC': 0.6,
    'PRIMER_DNA_CONC': 400.0,
    'PRIMER_INTERNAL_SALT_MONOVALENT': 50,
    'PRIMER_INTERNAL_SALT_DIVALENT': 1.5,
    'PRIMER_INTERNAL_DNTP_CONC': 0.0,
    'PRIMER_INTERNAL_DNA_CONC': 400.0,
    'PRIMER_INTERNAL_MAX_NS_ACCEPTED': 0,
    'PRIMER_MAX_SELF_END': 1000,
    'PRIMER_PAIR_MAX_COMPL_END': 1000,
    'PRIMER_MAX_SELF_ANY': 1000,
    'PRIMER_PAIR_MAX_COMPL_ANY': 1000,
    'PRIMER_INTERNAL_MAX_SELF_END': 1000,
    'PRIMER_INTERNAL_MAX_SELF_ANY': 1000,
    'PRIMER_SECONDARY_STRUCTURE_ALIGNMENT': 1,
    'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
    'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT': 1,
}


tunableRules = {
    'PRIMER_TM_BLOCK': {
        'importance': 7,
        'type': 'tm_block',
        'levels': [
            {'MIN': 58.0, 'OPT': 58.5, 'MAX': 59.5, 'DIFF': 1.0},
            {'MIN': 57.5, 'OPT': 58.5, 'MAX': 59.75, 'DIFF': 1.5},
            {'MIN': 57.0, 'OPT': 58.5, 'MAX': 60, 'DIFF': 2.0},
            {'MIN': 56.5, 'OPT': 58.5, 'MAX': 60.3, 'DIFF': 2.0},
        ],
    },
    'PRIMER_SIZE_BLOCK': {
        'importance': 3,
        'type': 'size_block',
        'levels': [
            {'MIN': 18, 'OPT': 22, 'MAX': 25},
            {'MIN': 17, 'OPT': 22, 'MAX': 26},
            {'MIN': 16, 'OPT': 22, 'MAX': 28},
        ],
    },
    'PRIMER_PRODUCT_SIZE_BLOCK': {
        'importance': 3,
        'type': 'product_size_block',
        'levels': [
            [[75, 150]],
            [[70, 160]],
            [[60, 180]],
        ],
    },
    'PRIMER_GC_BLOCK': {
        'importance': 4,
        'type': 'gc_block',
        'levels': [
            {'MIN': 40, 'OPT': 50, 'MAX': 60},
            {'MIN': 35, 'OPT': 50, 'MAX': 65},
        ],
    },
    'PRIMER_MAX_SELF_END_TH': {
        'importance': 6,
        'type': 'leq',
        'levels': [10, 20, 30],
        'primer3_key': 'PRIMER_MAX_SELF_END_TH',
    },
    'PRIMER_PAIR_MAX_COMPL_END_TH': {
        'importance': 6,
        'type': 'leq',
        'levels': [10, 20, 25],
        'primer3_key': 'PRIMER_PAIR_MAX_COMPL_END_TH',
    },
    'PRIMER_MAX_SELF_ANY_TH': {
        'importance': 5,
        'type': 'leq',
        'levels': [10, 20, 30, 40],
        'primer3_key': 'PRIMER_MAX_SELF_ANY_TH',
    },
    'PRIMER_PAIR_MAX_COMPL_ANY_TH': {
        'importance': 5,
        'type': 'leq',
        'levels': [20, 30, 35],
        'primer3_key': 'PRIMER_PAIR_MAX_COMPL_ANY_TH',
    },
    'PRIMER_MAX_HAIRPIN_TH': {
        'importance': 5,
        'type': 'leq',
        'levels': [10, 20, 30, 40],
        'primer3_key': 'PRIMER_MAX_HAIRPIN_TH',
    },
    'PRIMER_INTERNAL_TM_BLOCK': {
        'importance': 7,
        'type': 'internal_tm_block',
        'levels': [
            {'MIN': 67.0, 'OPT': 69.0, 'MAX': 72.0},
            {'MIN': 66.0, 'OPT': 69.0, 'MAX': 72.0},
            {'MIN': 65.0, 'OPT': 69.0, 'MAX': 72.0},
            {'MIN': 65.0, 'OPT': 69.0, 'MAX': 73.0},
            {'MIN': 65.0, 'OPT': 69.0, 'MAX': 73.0},
        ],
    },
    'PRIMER_INTERNAL_MAX_SELF_END_TH': {
        'importance': 6,
        'type': 'leq',
        'levels': [10, 15, 20, 25],
        'primer3_key': 'PRIMER_INTERNAL_MAX_SELF_END_TH',
    },
    'PRIMER_INTERNAL_MAX_SELF_ANY_TH': {
        'importance': 5,
        'type': 'leq',
        'levels': [10, 20, 30, 35],
        'primer3_key': 'PRIMER_INTERNAL_MAX_SELF_ANY_TH',
    },
    'PRIMER_INTERNAL_MAX_HAIRPIN_TH': {
        'importance': 5,
        'type': 'leq',
        'levels': [10, 20, 30, 40],
        'primer3_key': 'PRIMER_INTERNAL_MAX_HAIRPIN_TH',
    },
    'PRIMER_INTERNAL_SIZE_BLOCK': {
        'importance': 2,
        'type': 'internal_size_block',
        'levels': [
            {'MIN': 18, 'OPT': 22, 'MAX': 25},
            {'MIN': 17, 'OPT': 22, 'MAX': 28},
            {'MIN': 17, 'OPT': 22, 'MAX': 30},
        ],
    },
    'PRIMER_INTERNAL_GC_BLOCK': {
        'importance': 4,
        'type': 'internal_gc_block',
        'levels': [
            {'MIN': 45, 'OPT': 50, 'MAX': 55},
            {'MIN': 40, 'OPT': 50, 'MAX': 60},
            {'MIN': 35, 'OPT': 50, 'MAX': 65},
            {'MIN': 30, 'OPT': 50, 'MAX': 70},
        ],
    },
    'PRIMER_INTERNAL_MAX_POLY_X': {
        'importance': 1,
        'type': 'leq',
        'levels': [3, 4, 5],
        'primer3_key': 'PRIMER_INTERNAL_MAX_POLY_X',
    },
}


maxImportance = 7


# Edit this curve to adjust how quickly high-importance knobs become less likely
# to be relaxed in the scheduler.
def multiplierFn(importanceLevel: int) -> float:
    return 2 * importanceLevel


# Extra rules used after Primer3 generation. These enforce spacing and sequence
# constraints specific to the Kraken qPCR workflow.
krakenSpecificRules = {
    'enabled': True,
    'min_left_probe_gap': 5,
    'min_right_probe_gap': 10,
    'probe_no_g_prefix_len': 3,
}


__all__ = [
    'primer3Constants',
    'tunableRules',
    'maxImportance',
    'multiplierFn',
    'krakenSpecificRules',
]
