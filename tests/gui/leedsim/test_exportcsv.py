"""Tests for module exportcsv of viperleed.gui.leedsim."""

import re

import numpy as np
import pytest
from pytest_cases import fixture

from viperleed.gui.leedsim.exportcsv import _format_header_

_MODULE = 'viperleed.gui.leedsim.exportcsv'
_FRACTION_RE = r'-?\d+(/[1-9]\d+)?'
_FLOAT_RE = r'-?\d+(.\d*)?'
_DOMAIN_RE = r'\(?\d+\)?'
_COLUMN_LABELS_RE = re.compile(
    r'\(\s*h\s*[|]\s*k\s*\),'  # (h|k) fractional, with spaces around
    r'\s*h\s*,\s*k\s*,'        # h, k float, with spaces around
    r'\s*gx\s*,\s*gy\s*,'      # wave vectors, with spaces around
    # r'\s*domain\(s\)\s*,'      # domain labels  << does not match???
    )
_COLUMN_CONTENTS_RE = re.compile(
    rf'\(\s*{_FRACTION_RE}\s*[|]\s*{_FRACTION_RE}\s*\),'
    + rf'\s*{_FLOAT_RE},'*4                # h, k, gx, gy
    + r'\s*-?\d+,'                         # group index
    rf'\s*{_DOMAIN_RE}([+]{_DOMAIN_RE})*'  # domains
    )


class TestFormatHeader:
    """Tests for the _format_header_ helper."""

    @fixture(name='lengths')
    def fixture_lengths(self):
        """Return an example dictionary of column widths."""
        return {
            'numerator': 2,
            'denominator': 1,
            'hk_integer': 3,
            'g_integer': 4,
            'group': 5,
            'domains': 6,
            }

    @fixture(name='mock_impl')
    def fixture_mock_impl(self, mocker):
        """Replace implementation details with mocks."""
        return {
            'basis': mocker.patch(f'{_MODULE}._format_lattice_basis_',
                                  return_value='formatted basis'),
            'vectors': mocker.patch(f'{_MODULE}._format_basis_vectors_',
                                    return_value='formatted_vector')
            }

    @fixture(name='mock_leed')
    def fixture_mock_leed(self, mocker):
        """Return a fake LEEDPattern."""
        leed = mocker.MagicMock(
            n_domains=2,
            max_energy=123.4,
            superlattices=np.array([np.eye(2, dtype=int),
                                    2*np.eye(2, dtype=int)]),
            )
        leed.reciprocal_lattices = lattices = {
            'bulk': mocker.MagicMock(),
            'surf': mocker.MagicMock(),
            }
        lattices['bulk'].real_basis = np.eye(2)
        return leed

    @staticmethod
    def check_constant_contents(header):
        """Check the static contents of header lines."""
        *commented, column_labels = header
        assert not column_labels.startswith('#')
        # pylint: disable-next=magic-value-comparison
        assert 'ViPErLEED' in commented[0]
        assert _COLUMN_LABELS_RE.match(column_labels)

        # Other content that is always present
        header_joined = '\n'.join(header)
        expected = (
            'Max. LEED Energy',
            'Exporting beams from',
            )
        assert all(e in header_joined for e in expected)

        # Checks for Issue #86
        commented_lines = '\n'.join(commented).splitlines()
        assert all(line.startswith('#') for line in commented_lines)
        assert not any('"' in line for line in header)

    @pytest.mark.usefixtures('mock_impl')
    def test_basic(self, lengths, mock_leed):
        """Check the results of a call without optional arguments."""
        header = _format_header_(lengths, mock_leed)
        assert isinstance(header, list)
        self.check_constant_contents(header)

    def test_implementation(self, lengths, mock_leed, mock_impl, mocker):
        """Check expected calls to other helpers."""
        _format_header_(lengths, mock_leed)
        bulk = mock_leed.reciprocal_lattices['bulk']
        calls = {
            'basis': [mocker.call(latt.reciprocal_basis, latt.cell_shape)
                      for latt in mock_leed.reciprocal_lattices.values()],
            'vectors': [
                # This has two more calls for the superstructure bases.
                # We only check this in the n_calls below.
                mocker.call(bulk.reciprocal_basis),
                ],
            }
        n_calls = {'basis': 2, 'vectors': 3}
        for mock_name, expect_calls in calls.items():
            mock = mock_impl[mock_name]
            mock.assert_has_calls(expect_calls)
            assert mock.call_count == n_calls[mock_name]

    def test_missing_keys(self):
        """Check complaints when an the first argument is invalid."""
        bad_lengths = {'numerator': 2,
                       'denominator': 1}  # Missing others
        with pytest.raises(ValueError):
            _format_header_(bad_lengths, None)

    @pytest.mark.usefixtures('mock_impl')
    def test_subset_domains(self, lengths, mock_leed):
        """Check the outcome when exporting only selected domains."""
        header = _format_header_(lengths, mock_leed, domains=[1])
        self.check_constant_contents(header)
        header_joined = '\n'.join(header)
        assert 'Domain 1' not in header_joined
        assert 'Domain 2' in header_joined

    @pytest.mark.usefixtures('mock_impl')
    def test_with_name_source(self, lengths, mock_leed):
        """Check the results of a call with file and structure names."""
        fname = 'input.ini'
        struct_name = 'test_structure'
        header = _format_header_(lengths,
                                 mock_leed,
                                 source=fname,
                                 name=struct_name)
        assert fname in header[1]
        assert any(f'Structure: {struct_name}' in line for line in header)
        self.check_constant_contents(header)
