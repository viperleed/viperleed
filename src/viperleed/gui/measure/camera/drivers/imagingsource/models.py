"""Module models of viperleed.gui.measure.camera.drivers.imagingsource.

Contains useful data relative to specific camera models
currently fabricated by The Imaging Source.
Data comes from https://www.theimagingsource.com/products/

Last updated: 2021-11-10 Only GigE monochrome industrial cameras
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-11-10'
__license__ = 'GPLv3+'

from enum import Enum


class _Sensors(Enum):
    """Camera sensors. For now only dynamic range in bits is stored."""

    IMX546 = (12,)
    IMX267 = (12,)
    IMX304 = (12,)
    ICX618 = (12,)
    IMX287 = (12,)
    MT9V024 = (10,)
    AR0134 = (12,)
    ICX445 = (12,)
    PYTHON_1300 = (10,)
    IMX273 = (12,)
    ICX274 = (12,)
    IMX290 = (12,)
    IMX174 = (12,)
    IMX236 = (12,)
    IMX249 = (12,)
    IMX265 = (12,)
    IMX264 = (12,)
    AR0521 = (12,)
    MT9P031 = (12,)
    IMX178 = (12,)
    IMX226 = (12,)
    IMX183 = (12,)
    PYTHON_2000 = (10,)
    PYTHON_5000 = (10,)
    MT9J003 = (12,)
    ICX618ALA = (12,)
    ICX445ALA = (12,)
    MT9M021 = (12,)
    ICX274AL = (12,)
    IMX236LLJ = (12,)

    @property
    def dynamic_range(self):
        """Return the dynamic range of pixels in bits."""
        return self.value[0]


class ISModels(Enum):
    """Available model numbers for Imaging Source cameras."""

    # pylint: disable=invalid-name
    # The invalid-name would be issued for camera models
    # that contain the small "e" (PoE models). Easier to
    # disable the invalid-name in this case than rewriting
    # the attribute getter for the class.

    DMK_33GX546 = _Sensors.IMX546
    DMK_38GX267 = _Sensors.IMX267
    DMK_38GX304 = _Sensors.IMX304
    DMK_33G618 = _Sensors.ICX618
    DMK_33GX287 = _Sensors.IMX287
    DMK_33GV024 = _Sensors.MT9V024
    DMK_33GR0134 = _Sensors.AR0134
    DMK_33G445 = _Sensors.ICX445
    DMK_33GP1300 = _Sensors.PYTHON_1300
    DMK_33GP1300e = _Sensors.PYTHON_1300
    DMK_33GX273 = _Sensors.IMX273
    DMK_33G274 = _Sensors.ICX274
    DMK_33GX290 = _Sensors.IMX290
    DMK_33GX290e = _Sensors.IMX290
    DMK_33GX174	= _Sensors.IMX174
    DMK_33GX174e = _Sensors.IMX174
    DMK_33GX236	= _Sensors.IMX236
    DMK_33GX249	= _Sensors.IMX249
    DMK_33GX249e = _Sensors.IMX249
    DMK_33GX265	= _Sensors.IMX265
    DMK_33GX265e = _Sensors.IMX265
    DMK_33GX264	= _Sensors.IMX264
    DMK_33GX264e = _Sensors.IMX264
    DMK_33GR0521 = _Sensors.AR0521
    DMK_33GP031	= _Sensors.MT9P031
    DMK_33GX178	= _Sensors.IMX178
    DMK_33GX178e = _Sensors.IMX178
    DMK_33GX226	= _Sensors.IMX226
    DMK_33GX183	= _Sensors.IMX183
    DMK_33GP2000e = _Sensors.PYTHON_2000
    DMK_33GP5000e = _Sensors.PYTHON_5000
    DMK_33GJ003e = _Sensors.MT9J003
    DMK_23G618 = _Sensors.ICX618ALA
    DMK_23GV024	= _Sensors.MT9V024
    DMK_23G445 = _Sensors.ICX445ALA
    DMK_23GM021 = _Sensors.MT9M021
    DMK_23G274 = _Sensors.ICX274AL
    DMK_23GX236 = _Sensors.IMX236LLJ
    DMK_23GP031	= _Sensors.MT9P031
    DMK_23G618_I = _Sensors.ICX618ALA
    DMK_23GV024_I = _Sensors.MT9V024
    DMK_23GM021_I = _Sensors.MT9M021
    DMK_23G445_I = _Sensors.ICX445ALA
    DMK_23G274_I = _Sensors.ICX274AL
    DMK_23GP031_I = _Sensors.MT9P031
    # pylint: enable=invalid-name

    @property
    def dynamic_range(self):
        """Return the dynamic range of pixels in bits."""
        return self.value.dynamic_range
