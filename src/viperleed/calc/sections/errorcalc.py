"""Section Error calculation."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2021-03-18'
__license__ = 'GPLv3+'

import copy
import logging
import os
from pathlib import Path

from viperleed.calc.classes.r_error import R_Error
from viperleed.calc.classes.searchpar import SearchPar
from viperleed.calc.files import ioerrorcalc
from viperleed.calc.files.displacements import readDISPLACEMENTS_block
from viperleed.calc.sections.deltas import deltas as section_deltas
from viperleed.calc.sections.rfactor import rfactor as section_rfactor
from viperleed.calc.sections.superpos import superpos as section_superpos


logger = logging.getLogger(__name__)


def errorcalc(sl, rp):
    """Runs 1D error calculations based on DISPLACEMENTS."""

    if rp.domainParams:
        # TODO!! requested by Tilman!!
        # !!! Should be straightforward. Just take the beams from other domains
        # as constant (from refcalc), vary for the one domain as always
        logger.error("Error calculation not implemented for multiple domains.")
        return

    if rp.best_v0r is None:
        logger.info("Error calculation called without a stored inner "
                    "potential shift. Running R-factor calculation "
                    "from refcalc-fd.out to determine inner potential shift.")
        section_rfactor(sl, rp, index=11)
        logger.info("Finished R-factor pre-run, now starting error "
                    "calculation.\n")

    errors = []
    if rp.search_index > 0 and rp.search_index >= len(rp.disp_blocks):
        # running after one or more searches; read the last DISPLACEMENTS block
        rp.search_index -= 1
    seg_info = {"geo": "geometric", "vib": "vibration",
                "occ": "occupation"}
    for mode in "geo", "vib", "occ":
        sl.restoreOriState()  # reset positions, store any changes as offsets
        # read DISPLACEMENTS block - ONLY geo OR vib
        deltas_required = readDISPLACEMENTS_block(
            rp,
            sl,
            rp.disp_blocks[rp.search_index],
            only_mode=mode
            )
        rp.disp_block_read = True  # to prevent deltas segment from re-reading
        if not deltas_required:
            continue
        logger.info("\nStarting error calculations for "
                    + seg_info[mode] + " displacements.")
        # run delta calculations
        section_deltas(sl, rp)
        sl.deltas_initialized = True
        if rp.STOP:     # since this may stop the deltas, also check here
            return
        rp.generateSearchPars(sl)
        # group atoms under variation by linking
        atom_groups = [[at] for at in rp.search_atlist
                       if [sp for sp in rp.searchpars
                           if sp.atom == at and sp.steps > 1]]
        i = 0
        while i < len(atom_groups):
            found = None
            for at in atom_groups[i]:
                sps = [sp for sp in rp.searchpars if
                       sp.atom == at and sp.mode == mode and sp.el != "vac"
                       and ((isinstance(sp.linkedTo, SearchPar)
                             and sp.linkedTo.atom not in atom_groups[i]) or
                            (isinstance(sp.restrictTo, SearchPar)
                             and sp.restrictTo.atom not in atom_groups[i]))]
                if not sps:
                    continue
                if isinstance(sps[0].linkedTo, SearchPar):
                    found = [ag for ag in atom_groups
                             if sps[0].linkedTo.atom in ag][0]
                else:
                    found = [ag for ag in atom_groups
                             if sps[0].restrictTo.atom in ag][0]
                if found:
                    break
            if found:  # merge the two groups
                found.extend(atom_groups[i])
                atom_groups.pop(i)
            else:
                i += 1
        # run superpos and rfactor:
        for ag in atom_groups:
            # rp.search_atlist = ag  # limit variations to only that group
            only_vary = [sp for sp in rp.searchpars if sp.atom in ag]
            logger.info("\nNow calculating " + seg_info[mode] + " errors for "
                        "atom group: " + ", ".join(str(at) for at in ag))
            logger.info("Running superpos...")
            section_superpos(sl, rp, for_error=True, only_vary=only_vary)
            if rp.halt >= rp.HALTING:
                return
            if os.path.isfile("ROUTSHORT"):
                os.remove("ROUTSHORT")
            logger.info("Starting R-factor calculation...")
            rfaclist = section_rfactor(sl, rp, index=12, for_error=True,
                                       only_vary=only_vary)
            logger.info("Finished with " + seg_info[mode] + " errors for "
                        "atom group: " + ", ".join(str(at) for at in ag))
            error_disp_label = ag[0].displist[0].disp_labels[mode]
            error_lin_disps = ag[0].displist[0].disp_lin_steps[mode]
            errors.append(R_Error(rp.R_FACTOR_TYPE,
                                  ag, mode, rfaclist,
                                  error_disp_label,
                                  error_lin_disps,
                                  v0i=rp.V0_IMAG,
                                  energy_range=rp.total_energy_range()))
            if rp.halt >= rp.HALTING:
                return
    if len(errors) == 0:
        logger.info("Error calculation: Returning with no output.")
        return

    # Inform user that statistical error estimates are only available
    # for Pendry R-factor.
    if rp.R_FACTOR_TYPE != 1:
        logger.info("Estimates for statistical uncertainties "
                    "of parameters are only available for the Pendry "
                    "R-factor.")
    else:
        # Write var(R) to log as info.
        var_r_info = ioerrorcalc.extract_var_r(errors)
        if any(var_r_info.values()):
            var_str = []
            for mode, var_r in var_r_info.items():
                if not var_r:
                    continue
                var_str.append(
                    f"{mode}: {ioerrorcalc.format_col_content(var_r)}"
                    )
            var_str = ", ".join(var_str)
            logger.info(f"Found values for var(R): {var_str}")
        else:
            logger.info("Could not estimate var(R) for any error mode.")

    # Errors_summary.csv and Errors.zip
    (summary_csv_content,
     individual_files) = ioerrorcalc.generate_errors_csv(errors)
    ioerrorcalc.write_errors_archive(
        individual_files,
        compression_level=rp.ZIP_COMPRESSION_LEVEL,
        archive_fname='Errors.zip',
        )
    ioerrorcalc.write_errors_summary_csv(summary_csv_content,
                                         summary_fname='Errors_summary.csv')

    # Errors.pdf
    errors_figs = ioerrorcalc.make_errors_figs(errors, formatting=rp.PLOT_IV)
    ioerrorcalc.write_errors_pdf(errors_figs, filename='Errors.pdf')
