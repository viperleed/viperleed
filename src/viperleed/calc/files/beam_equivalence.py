"""Functions for reading and writing the BEAM_EQUIVALENCE file.

BEAM_EQUIVALENCE assigns calculated beams to experimental beams
explicitly, one beam assignment per line, for example
    (1 0) = (1 0) + 0.5*(-1 0)
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

from pathlib import Path
import re

import numpy as np
from quicktions import Fraction

from viperleed.calc.classes.beam import MAX_DENOMINATOR
from viperleed.calc.classes.beam_correspondence import BeamCorrespondence
from viperleed.calc.files.beams import snap_to_fraction
from viperleed.calc.lib.string_utils import strip_comments

# A non-negative number: an integer, a decimal, or a fraction like 1/2
_UNSIGNED = r'(?:\d+/\d+|\d+(?:\.\d*)?|\.\d+)'

# One beam index: a number, optionally with a sign
_INDEX = rf'[+-]?{_UNSIGNED}'

# A beam label: two indices in parentheses, separated by
# whitespace, a comma, or a vertical bar, e.g., (1/2 | -1)
_LABEL_RE = re.compile(
    rf'\(\s*(?P<first>{_INDEX})(?:\s*[,|]\s*|\s+)(?P<second>{_INDEX})\s*\)'
    )

# One term of the right-hand side of a beam assignment: a beam label,
# optionally preceded by a non-negative weight and '*', e.g., 1/3*(1 0)
_TERM_RE = re.compile(
    rf'(?:(?P<weight>{_UNSIGNED})\s*\*\s*)?'
    r'(?P<label>\([^()]*\))'
    )

# Relative tolerance within which the writer gives a weight as a
# fraction: well below that of BeamCorrespondence, so that only
# weights like 1/3 (rounded to a float) are written as fractions
_FRACTION_REL_TOL = 1e-12

_HEADER = ('# BEAM_EQUIVALENCE: experimental beam = (weight*) '
           'calculated beam + ...')


def parse_beam_label(label):
    """Return the (h, k) indices of a beam label as Fractions.

    Parameters
    ----------
    label : str
        A beam label such as '(1 0)', '(1/2, -1)', or '(0.5 | 0)'.
        Indices are integers, fractions, or decimals, separated by
        whitespace, a comma, or a vertical bar. Whitespace is allowed
        anywhere inside the parentheses and around them.

    Returns
    -------
    tuple
        The two indices (h, k) as Fractions. Decimals that are within
        the rounding to six decimals of a fraction with a small
        denominator are snapped to that fraction, as for IVBEAMS.

    Raises
    ------
    ValueError
        If `label` is not two indices in parentheses.
    """
    match = _LABEL_RE.fullmatch(label.strip())
    if not match:
        msg = (f'Invalid beam label {label.strip()!r}. Expected two '
               'indices in parentheses, e.g., (1/2 -1).')
        raise ValueError(msg)
    try:
        return tuple(snap_to_fraction(Fraction(match[group]))
                     for group in ('first', 'second'))
    except ZeroDivisionError:
        msg = f'Invalid beam label {label.strip()!r}: division by zero.'
        raise ValueError(msg) from None


def read_beam_equivalence(path):
    """Return the beam correspondence read from a BEAM_EQUIVALENCE file.

    Each non-empty line holds one beam assignment,
        <label> = <term> + <term> + ...
    with each term either a beam <label> (weight 1) or
    <weight>*<label>. Weights are non-negative integers, decimals,
    or fractions such as 1/3. Comments start with '!', '#', or '%'.

    Parameters
    ----------
    path : str or Path
        The file to read.

    Returns
    -------
    BeamCorrespondence
        One beam assignment per line, in the order of the file.

    Raises
    ------
    ValueError
        If a line cannot be parsed, if a beam is on the left-hand side
        of more than one line or more than once on one right-hand
        side, if a line has no non-zero weight, or if the file holds
        no beam assignment. The message names the file and, where
        applicable, the line number.
    """
    path = Path(path)
    assignments = {}
    first_line = {}  # {experimental beam: line number}
    with path.open('r', encoding='utf-8') as file:
        for line_nr, raw_line in enumerate(file, start=1):
            line = strip_comments(raw_line)
            if not line:
                continue
            try:
                exp_key, calc_weights = _parse_assignment(line)
            except ValueError as exc:
                msg = f'{path}, line {line_nr}: {exc}'
                raise ValueError(msg) from None
            if exp_key in first_line:
                msg = (f'{path}, line {line_nr}: Experimental beam '
                       f'{_format_label(exp_key)} is already assigned '
                       f'on line {first_line[exp_key]}.')
                raise ValueError(msg)
            first_line[exp_key] = line_nr
            assignments[exp_key] = calc_weights
    if not assignments:
        msg = f'{path}: No beam assignment found.'
        raise ValueError(msg)
    return BeamCorrespondence(assignments)


def _parse_assignment(line):
    """Return the experimental beam and {calculated beam: weight}."""
    try:
        left, right = line.split('=')
    except ValueError:
        msg = (f'Invalid line {line!r}. Expected one beam assignment, '
               'e.g., (1 0) = (1 0) + 0.5*(-1 0).')
        raise ValueError(msg) from None
    exp_key = parse_beam_label(left)
    calc_weights = {}
    for term in _split_terms(right):
        weight, calc_key = _parse_term(term)
        if calc_key in calc_weights:
            msg = (f'Calculated beam {_format_label(calc_key)} appears '
                   'more than once on the right-hand side.')
            raise ValueError(msg)
        calc_weights[calc_key] = weight
    # Validate the line on its own, for the line number
    try:
        BeamCorrespondence({exp_key: calc_weights})
    except OverflowError:
        msg = 'A weight is too large to be represented as a float.'
        raise ValueError(msg) from None
    return exp_key, calc_weights


def _split_terms(right):
    """Return the terms of the right-hand side of a beam assignment."""
    # Split at '+' outside parentheses: labels may contain a '+'
    terms = re.split(r'\+(?![^()]*\))', right)
    if any(not term.strip() for term in terms):
        msg = (f'Invalid right-hand side {right.strip()!r}. Expected '
               "terms like (1 0) or 0.5*(1 0), separated by '+'.")
        raise ValueError(msg)
    return terms


def _parse_term(term):
    """Return the weight and calculated beam of one term."""
    match = _TERM_RE.fullmatch(term.strip())
    if not match:
        msg = (f'Invalid term {term.strip()!r}. Expected a beam, e.g., '
               '(1 0), optionally with a non-negative weight, e.g., '
               '0.5*(1 0) or 1/3*(1 0).')
        raise ValueError(msg)
    weight = match['weight']
    try:
        weight = Fraction(weight) if weight else Fraction(1)
    except ZeroDivisionError:
        msg = f'Invalid weight {weight!r}: division by zero.'
        raise ValueError(msg) from None
    return weight, parse_beam_label(match['label'])


def write_beam_equivalence(correspondence, path):
    """Write a beam correspondence to a BEAM_EQUIVALENCE file.

    The file starts with a comment line, followed by one line per
    beam assignment, in the order of `correspondence`. Indices are
    written as integers or fractions, e.g., (1/2 -1). Weights are
    left out for beam assignments whose weights are all equal.

    Parameters
    ----------
    correspondence : BeamCorrespondence
        The beam correspondence to write.
    path : str or Path
        The file to write to. It is overwritten if it exists.

    Returns
    -------
    None.
    """
    lines = [_HEADER]
    lines.extend(format_assignment(exp_key, calc_weights)
                 for exp_key, calc_weights
                 in correspondence.assignments.items())
    Path(path).write_text('\n'.join(lines) + '\n', encoding='utf-8')


def format_assignment(exp_key, calc_weights):
    """Return the BEAM_EQUIVALENCE line of one beam assignment.

    Weights are left out if they are all equal.

    Parameters
    ----------
    exp_key : tuple
        The (h, k) Fractions of the experimental beam.
    calc_weights : dict
        Weight of each calculated beam, keyed by (h, k) Fractions.

    Returns
    -------
    str
        The line, without a trailing newline.
    """
    equal = len(set(calc_weights.values())) == 1
    terms = (_format_label(calc_key) if equal
             else f'{_format_weight(weight)}*{_format_label(calc_key)}'
             for calc_key, weight in calc_weights.items())
    joined = ' + '.join(terms)
    return f'{_format_label(exp_key)} = {joined}'


def _format_label(key):
    """Return the label of a beam from its (h, k) Fractions."""
    first, second = key
    return f'({first} {second})'


def _format_weight(weight):
    """Return a string for an assignment weight that reads back as it.

    A fraction with a small denominator, if it equals `weight` within
    a relative tolerance of _FRACTION_REL_TOL, otherwise the shortest
    decimal that reads back as `weight` (without exponent, which the
    reader does not accept).
    """
    fraction = Fraction(weight).limit_denominator(MAX_DENOMINATOR)
    if abs(fraction - weight) <= _FRACTION_REL_TOL * weight:
        return str(fraction)
    return np.format_float_positional(weight, trim='-')
