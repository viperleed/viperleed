"""Module l_max of viperleed.calc.classes.rparams.special.

Defines the LMax class, a convenience container for parameter LMAX.
"""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-11-11'
__license__ = 'GPLv3+'


from .base import SpecialParameter, SpecialParameterError


class RFactorTypeError(SpecialParameterError):
    """Exception for RFactorType-related errors."""


class RFactorType(SpecialParameter, param='R_FACTOR_TYPE'):
    """Enum-like parameter for selecting the R-factor type.

    This class provides a robust interface for specifying which
    R-factor is used. It accepts both integer codes and a range of
    string identifiers.

    Canonical mapping:
        1 → 'pendry'
        2 → 'r2'
        3 → 'zj'
        4 → 'smooth'

    Parameters
    ----------
    value : {int, str, RFactorType}
        The value to parse. Integers 1–4 map directly to the canonical
        names above. Strings are case-insensitive and may include any
        of the registered synonyms.

    Attributes
    ----------
    id : int
        Integer code identifying the selected R-factor type.
    name : str
        Canonical string name of the selected R-factor type.

    Notes
    -----
    This class is automatically registered under the parameter name
    ``'R_FACTOR_TYPE'`` via :meth:`SpecialParameter.__init_subclass__`.
    """

    # Define synonyms for matching
    _ALIASES = {
        'pendry': ['pendry_r', 'rp', 'pendry', 'pend', 'p'],
        'R2': ['r2', 'r squared'],
        'zj': ['zanazzi jona', 'zanazzi-jona', 'zanazzi_jona', 'zannazi', 'z'],
        'smooth': [
            'smoothed',
            'smooth pendry',
            'smooth-pendry',
            'schmid',
            's',
        ],
    }

    _CODE_TO_NAME = {
        1: 'pendry',
        2: 'r2',
        3: 'zj',
        4: 'smooth',
    }

    # Constant integer codes
    PENDRY = 1
    R2 = 2
    ZJ = 3
    SMOOTH = 4

    # Build reverse lookup including aliases
    _NAME_TO_CODE = None  # filled below

    def __init__(self, value):
        r_fac_id, name = self._coerce(value)
        self.id = r_fac_id
        self.name = name

    # --------- construction helpers ---------
    @classmethod
    def from_value(cls, value):
        """Create an instance from ``value``.

        Parameters
        ----------
        value : {int, str, RFactorType}
            Input to parse into an R-factor type.

        Returns
        -------
        RFactorType
            Parsed instance.
        """
        return cls(value)

    @classmethod
    def _build_reverse_maps(cls):
        """(Re)build the normalized name->code map including aliases."""
        name_to_code = {}
        # canonical names
        for code, canon in cls._CODE_TO_NAME.items():
            name_to_code[cls._norm_name(canon)] = code
        # aliases
        for canon, aliases in cls._ALIASES.items():
            code = next(
                k
                for k, v in cls._CODE_TO_NAME.items()
                if v.lower() == canon.lower()
            )
            for a in aliases:
                name_to_code[cls._norm_name(a)] = code
        cls._NAME_TO_CODE = name_to_code

    @classmethod
    def _coerce(cls, value):
        """Parse input into (code, name)."""
        if cls._NAME_TO_CODE is None:
            cls._build_reverse_maps()

        # Already an instance
        if isinstance(value, cls):
            return value.id, value.name

        # Int-like
        if isinstance(value, int):
            try:
                name = cls._CODE_TO_NAME[value]
            except KeyError:
                msg = f'Invalid R_FACTOR_TYPE code: {value}'
                raise RFactorTypeError(msg) from None
            return value, name

        # String-like
        if isinstance(value, str):
            v = value.strip()
            if v.isdigit():
                code = int(v)
                try:
                    name = cls._CODE_TO_NAME[code]
                except KeyError:
                    msg = f'Invalid R_FACTOR_TYPE code: {value}'
                    raise RFactorTypeError(msg) from None
                return code, name
            # name/alias
            key = cls._norm_name(v)
            try:
                code = cls._NAME_TO_CODE[key]
            except KeyError:
                allowed = ', '.join(cls._CODE_TO_NAME.values())
                msg = f"Invalid R_FACTOR_TYPE name: '{value}'. Allowed: {allowed}"
                raise RFactorTypeError(msg) from None
            return code, cls._CODE_TO_NAME[code]

        msg = f'Cannot parse R_FACTOR_TYPE from {type(value).__name__}: {value!r}'
        raise TypeError(msg)

    @staticmethod
    def _norm_name(input_str):
        """Normalize a name for robust matching."""
        # unify type & whitespace
        output_str = input_str.strip().casefold()
        # replace common separators
        return (
            output_str.replace('²', '2')
            .replace('^', '')
            .replace('_', '')
            .replace('-', '')
            .replace(' ', '')
        )

    def __repr__(self):
        """Return a representation string of this RFactorType."""
        return f"RFactorType(code={self.id}, name='{self.name}')"

    def __str__(self):
        """Return a string version of this RFactorType."""
        return self.name

    def __int__(self):
        """Return the integer code of this RFactorType."""
        return self.id

    def __eq__(self, other):
        """Compare with int, str, or RFactorType."""
        try:
            r_fac_id, _ = self._coerce(other)
        except (TypeError, ValueError):
            return NotImplemented
        return self.id == r_fac_id
