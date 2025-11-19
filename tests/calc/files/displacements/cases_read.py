"""Test cases for reading a DISPLACEMENTS file."""

__authors__ = ('Michele Riva (@michele-riva)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-06-08'
__license__ = 'GPLv3+'


class CasesEmptyFile:
    """Collection of cases of DISPLACEMENTS file with no valid content."""

    def case_comments_only(self, write_displacements):
        """Prepare a DISPLACEMENTS file with no contents."""
        contents = """
# A comment line with a hash character
! == SEARCH commented with an exclamation mark
% = GEO_DELTA   ! A commented tag
#    Fe* L(1-5) z = -0.5 0.5 0.05   ! And a commented valid line
"""
        return write_displacements(contents)

    def case_empty_lines(self, write_displacements):
        """Prepare a DISPLACEMENTS file with no contents."""
        return write_displacements('\n  \n      \n')

    def case_only_loop_tags(self, write_displacements):
        """Prepare a DISPLACEMENTS file with only <loop> tags."""
        contents = """
<loop>

</loop>
"""
        return write_displacements(contents)

    def case_no_lines(self, write_displacements):
        """Prepare a DISPLACEMENTS file with no contents."""
        return write_displacements('')


class CasesInvalidDomainDisplacements:
    """Invalid DISPLACEMENTS files for multi-domain calculations."""

    def case_all_domain_names_mismatched(self, domains, write_displacements):
        """Prepare a DISPLACEMENTS with invalid DOMAIN names."""
        contents = """
== DOMAIN not-among-the-names
  = GEO_DELTA
     * L(1-5) xy[1 0] = -0.5 0.5 0.05"""
        return domains, write_displacements(contents)

    def case_only_blocks_outside_domain(self, domains, write_displacements):
        """Prepare a single-domain DISPLACEMENTS for a DOMAINs calculation."""
        contents = """
= GEO_DELTA   # Not assigned to a DOMAIN, will be skipped
   * L(1-5) z = 0 0.2 0.05"""
        return domains, write_displacements(contents)

    def case_all_domain_blocks_empty(self, domains, write_displacements):
        """Prepare a DISPLACEMENTS with empty domain DISPLACEMENTS."""
        contents = """
= DOMAIN 0
= DOMAIN 1
    # comments only
= DOMAIN 2
    = GEO_DELTA  ! empty"""
        return domains, write_displacements(contents)
