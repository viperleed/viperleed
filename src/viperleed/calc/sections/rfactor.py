"""Section Rfactor."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-08-11'
__license__ = 'GPLv3+'

import logging
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

from viperleed.calc.constants import DEFAULT_OUT, DEFAULT_SUPP, DEFAULT_TENSORS
from viperleed.calc.files import iorfactor, iotensors
from viperleed.calc.files.iorefcalc import readFdOut
from viperleed.calc.files.iorfactor import (beamlist_to_array,
                                            write_Rfactorpdf_from_splines)
from viperleed.calc.lib import fs_utils, leedbase, spline_interpolation
from viperleed.calc.lib import rfactor as rfactor_lib
from viperleed.calc.lib.checksums import validate_multiple_files
from viperleed.calc.lib.rfactor.utils import average_beam_array

logger = logging.getLogger(__name__)


def rfactor(
    sl, rp, index, for_error=False, only_vary=None
):  # TODO: Parameters __doc__
    """Runs the r-factor calculation for either the reference calculation
    (index 11) or the superpos (index 12)."""
    if rp.THEO_ENERGIES.n_energies < 2:
        logger.info(
            'Only one theoretical energy found: Cannot '
            'calculate a meaningful R-Factor. Stopping...'
        )
        return []
    if for_error and rp.best_v0r is None:
        # Cannot proceed with planned section
        err_msg = (
            'Cannot calculate R-factor for error without a stored'
            'V0r value. Execute a normal R-factor calculation first.'
        )
        rp.setHaltingLevel(3)
        logger.error(err_msg)
        raise RuntimeError(err_msg)

    # Both refcalc (11) and superpos should work fine with new R-factor
    if index == 11:
        name = 'refcalc'
    else:
        name = 'superpos'
    _need_spectra = (
        index == 11
        and not len(rp.refcalc_fdout)
        or index == 12
        and not len(rp.superpos_specout)
    )
    if _need_spectra:
        _fetch_and_check_spectra(rp, index, name)
    theobeams = rp.theobeams[name]
    expbeams = rp.expbeams

    # Branch off for new R factor calculation
    if not rp.R_FACTOR_LEGACY:
        return run_new_rfactor(sl, rp, for_error, name, theobeams, expbeams)
    return run_legacy_rfactor(
        sl, rp, for_error, name, theobeams, index, only_vary
    )


def _fetch_and_check_spectra(rp, index, name):
    """Look for full-dynamic spectra, then check consistency."""
    if index == 11:
        fn = Path('refcalc-fd.out')
    else:
        fn = Path('superpos-spec.out')
    directory = None
    path = None
    if fn.is_file():
        directory = Path.cwd().name
        path = fn
    elif (DEFAULT_OUT / fn).is_file():
        directory = DEFAULT_OUT
        path = DEFAULT_OUT / fn
    elif index == 11:
        # try getting from Tensors
        iotensors.fetch_unpacked_tensor(rp.TENSOR_INDEX)
        directory = Path(f'{DEFAULT_TENSORS}_{rp.TENSOR_INDEX:03d}')
        if (DEFAULT_TENSORS / directory / fn).is_file():
            path = DEFAULT_TENSORS / directory / fn

    if path:
        logger.warning(
            'R-factor calculation was called without stored '
            f'spectrum data. Reading from file {fn} in {directory} '
            'folder...'
        )
    else:
        logger.error(
            'Cannot execute R-factor calculation: no stored '
            f'spectrum data and no {fn} file was found.'
        )
        raise RuntimeError('No spectrum data found')
    try:
        theobeams, theospec = readFdOut(readfile=path)
    except Exception:
        logger.error(f'Failed to read {path}')
        raise

    # if we haven't returned, then it was read. check data vs ivbeams:
    if index == 11:
        rp.refcalc_fdout = theospec
    else:
        rp.superpos_specout = theospec
    rp.theobeams[name] = theobeams

    eps = 1e-3
    consistent = len(rp.ivbeams) == len(theobeams) and all(
        any(ivbeam.isEqual(theobeam, eps=eps) for theobeam in theobeams)
        for ivbeam in rp.ivbeams
    )
    if not consistent:
        logger.error(
            'The list of beams read from IVBEAMS is not equivalent '
            f'to the list of beams in {path}. R-Factor calculation '
            'cannot proceed.'
        )
        raise ValueError('Contradiction in beam sets')


def run_new_rfactor(sl, rp, for_error, name, theobeams, expbeams):
    logger.warning(
        'Using new R-factor calculation. This is still experimental!'
    )

    # select R-factor function by name
    r_func = rfactor_lib.select_rfactor(rp.R_FACTOR_TYPE.name)
    logger.info(f'Calculating R-factor type {rp.R_FACTOR_TYPE.name}.')

    (_, theo_range, iv_shift, intpol_step) = (
        iorfactor.prepare_rfactor_energy_ranges(
            rp, theobeams, for_error, n_expand=(rp.INTPOL_DEG - 1) // 2
        )
    )
    out_grid = np.arange(
        theo_range.min, theo_range.max + 0.1 * intpol_step, intpol_step
    )

    if rp.INTPOL_DEG != 3:
        logger.error(
            'New R-factor calculation only supports cubic '
            'interpolation (INTPOL_DEG=3).'
        )
        rp.setHaltingLevel(1)
        return []

    # find correspondence experimental to theoretical beams:
    beam_correspondence = leedbase.getBeamCorrespondence(sl, rp)
    beam_correspondence = tuple(beam_correspondence)
    logger.debug(f'Beam correspondence: {beam_correspondence}')

    # calculate splines from the spectra
    logger.debug('Constructing splines from data...')

    ## experiment
    exp_energies, _, _, exp_intensities = beamlist_to_array(rp.expbeams)
    exp_spline = spline_interpolation.interpolate_ragged_array(
        x=exp_energies,
        y=exp_intensities,
        bc_type='not-a-knot',  # TODO: from Rfactor parameter?
    )
    exp_spline = spline_interpolation.CachedSpline(exp_spline)

    ## theory
    (theo_grid, _, _, theo_intensities) = iorfactor.beamlist_to_array(theobeams)
    theo_intensities = average_beam_array(theo_intensities, beam_correspondence)
    theo_spline = spline_interpolation.interpolate_ragged_array(
        x=theo_grid,
        y=theo_intensities,
        bc_type='not-a-knot',  # TODO: from Rfactor parameter?
    )
    theo_spline = spline_interpolation.CachedSpline(theo_spline)

    logger.debug('Sweeping V0r...')
    # sweep V0r shifts - TODO: continuous shifts?
    shifts = np.arange(
        rp.IV_SHIFT_RANGE.start,
        rp.IV_SHIFT_RANGE.stop + rp.IV_SHIFT_RANGE.step,
        rp.IV_SHIFT_RANGE.step,
    )
    r_values = []
    for shift in shifts:
        r_values.append(
            r_func(
                theo_spline,
                rp.V0_IMAG,
                intpol_step,
                out_grid,
                data_spline_1=exp_spline,
                data_spline_2=theo_spline,
                shift_2nd_spline=shift,
                groups=None,
            )
        )

    # use best shift
    r_values = np.array(r_values)
    best_index = np.argmin(r_values)
    best_shift = shifts[best_index]
    logger.info(
        f'Best inner potential shift: {best_shift:.2f} eV '
        f'with R = {r_values[best_index]:.4f}'
    )
    rp.best_v0r = best_shift

    if best_shift != 0:
        # re-build theoretical splines one more time on the shifted grid
        theo_spline = spline_interpolation.interpolate_ragged_array(
            x=theo_grid + best_shift,
            y=theo_intensities,
            bc_type='not-a-knot',  # TODO: from Rfactor parameter?
            )
        theo_spline = spline_interpolation.CachedSpline(theo_spline)

    # calculate R-factors at best shift
    logger.debug('Calculating R-factors...')
    r_fac_overall = r_func(
        rp.V0_IMAG,
        intpol_step,
        out_grid,
        data_spline_1=exp_spline,
        data_spline_2=theo_spline,
        groups=None,
        )

    integer_fractional_mask = determine_integer_or_fractional(rp)
    r_fac_integer, r_fac_fractional = r_func(
        rp.V0_IMAG,
        intpol_step,
        out_grid,
        data_spline_1=exp_spline,
        data_spline_2=theo_spline,
        groups=integer_fractional_mask,
        num_groups=2,
    )

    r_fac_per_beam = r_func(
        rp.V0_IMAG,
        intpol_step,
        out_grid,
        data_spline_1=exp_spline,
        data_spline_2=theo_spline,
        groups='beam',
    )

    # plotting
    if rp.PLOT_IV['plot']:
        outname = f'Rfactor_plots_{name}.pdf'
        aname = f'Rfactor_analysis_{name}.pdf'
        write_Rfactorpdf_from_splines(rp,
                                      theo_spline,
                                      exp_spline,
                                      out_grid,
                                      r_fac_per_beam,
                                      outname=outname,
                                      analysis_file=aname,
                                      )

    # store R-factors
    rp.last_R = r_fac_overall
    rp.stored_R[name] = (r_fac_overall, r_fac_integer, r_fac_fractional)

    # return the per-beam R-factors
    return r_fac_per_beam


def run_legacy_rfactor(sl, rp, for_error, name, theobeams, index, only_vary):
    if index == 11:
        theospec = rp.refcalc_fdout
    elif index == 12:
        theospec = rp.superpos_specout
    # WEXPEL before PARAM, to make sure number of exp. beams is correct
    try:
        iorfactor.writeWEXPEL(sl, rp, theobeams, for_error=for_error)
    except Exception:
        logger.error('Exception during writeWEXPEL: ')
        raise
    try:
        iorfactor.writeRfactPARAM(
            rp, theobeams, for_error=for_error, only_vary=only_vary
        )
    except Exception:
        logger.error('Exception during writeRfactPARAM: ')
        raise

    # get fortran files and compile
    try:
        tl_source = rp.get_tenserleed_directory()
        tl_path = tl_source.path
        libpath = tl_path / 'lib'
        libname = next(libpath.glob('rfacsb*'))  # StopIteration??
        srcpath = tl_path / 'src'
        srcname = next(srcpath.glob('rfactor.*'))  # StopIteration??
        shutil.copy2(libname, libname.name)  # Copy here from source
        shutil.copy2(srcname, srcname.name)  # Copy here from source
    except Exception:
        logger.error(
            'Error getting TensErLEED files for r-factor calculation:'
        )
        raise
    if rp.SUPPRESS_EXECUTION:
        logger.warning(
            'SUPPRESS_EXECUTION parameter is on. R-factor '
            'calculation will not proceed. Stopping...'
        )
        rp.setHaltingLevel(3)
        return []
    # validate checksums
    if not rp.TL_IGNORE_CHECKSUM:
        files_to_check = (Path(libpath) / libname, Path(srcpath) / srcname)
        validate_multiple_files(
            files_to_check, logger, 'R-factor', rp.TL_VERSION
        )

    logger.info('Compiling fortran input files...')
    rfacname = f'rfactor-{rp.timestamp}'
    if rp.FORTRAN_COMP[0] == '':
        rp.getFortranComp()
    # compile
    _pre, _post = rp.FORTRAN_COMP
    ctasks = [
        (f'{_pre} -o rfacsb.o -c', libname.name, _post),
        (f'{_pre} -o main.o -c', srcname.name, _post),
        (f'{_pre} -o {rfacname}', 'rfacsb.o main.o', _post),
    ]

    # Compile
    compile_log = Path('compile-rfactor.log')
    try:
        leedbase.fortran_compile_batch(ctasks, logname=compile_log)
    except Exception:
        logger.error('Error compiling fortran files: ', exc_info=True)
        # move compile log to compile_logs directory
        leedbase.copy_compile_log(rp, compile_log, 'rfactor-compile')
        raise

    # move log file to supp
    try:
        fs_utils.move(compile_log, 'compile_logs' / compile_log)
    except OSError:
        logger.warning(f'Could not move {compile_log} to {DEFAULT_SUPP}')
    # run
    rfaclogname = Path(rfacname).with_suffix('.log')
    logger.info(
        'Starting R-factor calculation...\n'
        f'R-factor log will be written to file {rfaclogname}'
    )
    try:
        with rfaclogname.open('w') as log:
            subprocess.run(
                str(Path(rfacname).resolve()),
                input=theospec,
                encoding='ascii',
                stdout=log,
                stderr=log,
            )
    except Exception:
        logger.error(
            'Error during R-factor calculation. Also check R-factor log file.'
        )
        raise
    logger.info('Finished R-factor calculation. Processing files...')
    # rename and move files
    try:
        os.rename('WEXPEL', 'rfactor-WEXPEL')
    except OSError:
        logger.warning(
            'Failed to rename R-factor input file WEXPEL to rfactor-WEXPEL'
        )
    try:
        os.rename('PARAM', 'rfactor-PARAM')
    except OSError:
        logger.warning(
            'Failed to rename R-factor input file PARAM to rfactor-PARAM'
        )
    if not Path('ROUT').is_file():
        logger.error('No ROUT file was found after R-Factor calculation!')
        rp.setHaltingLevel(2)
        return []

    # read output
    if for_error:
        try:
            rfaclist = iorfactor.readROUTSHORT()
        except Exception:
            logger.error(
                'Error reading ROUTSHORT file', exc_info=rp.is_debug_mode
            )
            rp.setHaltingLevel(2)
        return rfaclist

    try:
        (rfac, r_int, r_frac), v0rshift, rfaclist = iorfactor.readROUT()
    except Exception:  # TODO catch correct exception
        logger.error('Error reading ROUT file', exc_info=rp.is_debug_mode)
        rp.setHaltingLevel(2)
        return []

    logger.info(
        'With inner potential shift of {:.2f} eV: R = {:.4f}\n'.format(
            v0rshift, rfac
        )
    )
    rp.best_v0r = v0rshift
    dir_list = [Path(), Path(DEFAULT_OUT)]
    for dir_name in dir_list:
        for f_name in dir_name.glob(f'R_OUT*'):
            if not f_name.is_file():
                continue
            try:  # delete old R_OUT files
                f_name.unlink()
            except Exception:
                pass
    if rfac <= 0.00001:
        logger.warning(
            'ROUT reports R-Factor as zero. This may mean '
            'something went wrong in the reference '
            'calculation or in the R-factor calculation. '
            'If you are comparing with pseudo-experiment data, '
            'you can ignore this warning.'
        )

    f_name = f'R_OUT_{name}_R={rfac:.4f}'
    rp.last_R = rfac
    rp.stored_R[name] = (rfac, r_int, r_frac)
    try:
        os.rename('ROUT', f_name)
    except OSError:
        logger.warning(
            f'Failed to rename R-factor output file ROUT to {f_name}'
        )
    if len(rfaclist) != len(rp.expbeams):
        logger.warning(
            'Failed to read R-Factors per beam from R-factor output file ROUT.'
        )
        rfaclist = [-1] * len(rp.expbeams)

    if rp.PLOT_IV['plot']:
        outname = f'Rfactor_plots_{name}.pdf'
        aname = f'Rfactor_analysis_{name}.pdf'
        labelstyle = 'overbar' if rp.PLOT_IV['overbar'] else 'minus'
        labelwidth = max(
            beam.getLabel(style=labelstyle)[1] for beam in rp.expbeams
        )
        try:
            iorfactor.writeRfactorPdf(
                [
                    (b.getLabel(lwidth=labelwidth, style=labelstyle)[0], r)
                    for b, r in zip(rp.expbeams, rfaclist)
                ],
                outName=outname,
                analysisFile=aname,
                v0i=rp.V0_IMAG,
                formatting=rp.PLOT_IV,
            )
        except Exception:  # TODO catch correct exception
            logger.warning('Error plotting R-factors.', exc_info=True)
    return rfaclist


def determine_integer_or_fractional(rp):
    """Determine whether beams are integer or fractional."""
    iorf = []
    for i, beam in enumerate(rp.expbeams):
        if beam.hk[0] % 1.0 != 0.0 or beam.hk[1] % 1.0 != 0.0:
            iorf.append(1)
        else:
            iorf.append(0)
    return tuple(iorf)
