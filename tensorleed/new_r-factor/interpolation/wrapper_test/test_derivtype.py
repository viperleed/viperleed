from __future__ import print_function, absolute_import, division
import _test_derivtype
import f90wrap.runtime
import logging

class Interpolation(f90wrap.runtime.FortranModule):
    """
    Module interpolation
    
    
    Defined at test_deriv_types.f90 lines 14-68
    
    """
    @f90wrap.runtime.register_class("test_derivtype.info_values")
    class info_values(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=info_values)
        
        
        Defined at test_deriv_types.f90 lines 21-28
        
        """
        def __init__(self, handle=None):
            """
            self = Info_Values()
            
            
            Defined at test_deriv_types.f90 lines 21-28
            
            
            Returns
            -------
            this : Info_Values
            	Object to be constructed
            
            
            Automatically generated constructor for info_values
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _test_derivtype.f90wrap_info_values_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Info_Values
            
            
            Defined at test_deriv_types.f90 lines 21-28
            
            Parameters
            ----------
            this : Info_Values
            	Object to be destructed
            
            
            Automatically generated destructor for info_values
            """
            if self._alloc:
                _test_derivtype.f90wrap_info_values_finalise(this=self._handle)
        
        @property
        def deg(self):
            """
            Element deg ftype=integer  pytype=int
            
            
            Defined at test_deriv_types.f90 line 22
            
            """
            return _test_derivtype.f90wrap_info_values__get__deg(self._handle)
        
        @deg.setter
        def deg(self, deg):
            _test_derivtype.f90wrap_info_values__set__deg(self._handle, deg)
        
        @property
        def n_knots(self):
            """
            Element n_knots ftype=integer  pytype=int
            
            
            Defined at test_deriv_types.f90 line 23
            
            """
            return _test_derivtype.f90wrap_info_values__get__n_knots(self._handle)
        
        @n_knots.setter
        def n_knots(self, n_knots):
            _test_derivtype.f90wrap_info_values__set__n_knots(self._handle, n_knots)
        
        @property
        def nt(self):
            """
            Element nt ftype=integer  pytype=int
            
            
            Defined at test_deriv_types.f90 line 24
            
            """
            return _test_derivtype.f90wrap_info_values__get__nt(self._handle)
        
        @nt.setter
        def nt(self, nt):
            _test_derivtype.f90wrap_info_values__set__nt(self._handle, nt)
        
        @property
        def kl(self):
            """
            Element kl ftype=integer  pytype=int
            
            
            Defined at test_deriv_types.f90 line 25
            
            """
            return _test_derivtype.f90wrap_info_values__get__kl(self._handle)
        
        @kl.setter
        def kl(self, kl):
            _test_derivtype.f90wrap_info_values__set__kl(self._handle, kl)
        
        @property
        def ku(self):
            """
            Element ku ftype=integer  pytype=int
            
            
            Defined at test_deriv_types.f90 line 25
            
            """
            return _test_derivtype.f90wrap_info_values__get__ku(self._handle)
        
        @ku.setter
        def ku(self, ku):
            _test_derivtype.f90wrap_info_values__set__ku(self._handle, ku)
        
        @property
        def nleft(self):
            """
            Element nleft ftype=integer  pytype=int
            
            
            Defined at test_deriv_types.f90 line 26
            
            """
            return _test_derivtype.f90wrap_info_values__get__nleft(self._handle)
        
        @nleft.setter
        def nleft(self, nleft):
            _test_derivtype.f90wrap_info_values__set__nleft(self._handle, nleft)
        
        @property
        def nright(self):
            """
            Element nright ftype=integer  pytype=int
            
            
            Defined at test_deriv_types.f90 line 26
            
            """
            return _test_derivtype.f90wrap_info_values__get__nright(self._handle)
        
        @nright.setter
        def nright(self, nright):
            _test_derivtype.f90wrap_info_values__set__nright(self._handle, nright)
        
        @property
        def lhs_rows(self):
            """
            Element lhs_rows ftype=integer  pytype=int
            
            
            Defined at test_deriv_types.f90 line 27
            
            """
            return _test_derivtype.f90wrap_info_values__get__lhs_rows(self._handle)
        
        @lhs_rows.setter
        def lhs_rows(self, lhs_rows):
            _test_derivtype.f90wrap_info_values__set__lhs_rows(self._handle, lhs_rows)
        
        @property
        def lhs_cols(self):
            """
            Element lhs_cols ftype=integer  pytype=int
            
            
            Defined at test_deriv_types.f90 line 27
            
            """
            return _test_derivtype.f90wrap_info_values__get__lhs_cols(self._handle)
        
        @lhs_cols.setter
        def lhs_cols(self, lhs_cols):
            _test_derivtype.f90wrap_info_values__set__lhs_cols(self._handle, lhs_cols)
        
        @property
        def rhs_cols(self):
            """
            Element rhs_cols ftype=integer  pytype=int
            
            
            Defined at test_deriv_types.f90 line 28
            
            """
            return _test_derivtype.f90wrap_info_values__get__rhs_cols(self._handle)
        
        @rhs_cols.setter
        def rhs_cols(self, rhs_cols):
            _test_derivtype.f90wrap_info_values__set__rhs_cols(self._handle, rhs_cols)
        
        def __str__(self):
            ret = ['<info_values>{\n']
            ret.append('    deg : ')
            ret.append(repr(self.deg))
            ret.append(',\n    n_knots : ')
            ret.append(repr(self.n_knots))
            ret.append(',\n    nt : ')
            ret.append(repr(self.nt))
            ret.append(',\n    kl : ')
            ret.append(repr(self.kl))
            ret.append(',\n    ku : ')
            ret.append(repr(self.ku))
            ret.append(',\n    nleft : ')
            ret.append(repr(self.nleft))
            ret.append(',\n    nright : ')
            ret.append(repr(self.nright))
            ret.append(',\n    lhs_rows : ')
            ret.append(repr(self.lhs_rows))
            ret.append(',\n    lhs_cols : ')
            ret.append(repr(self.lhs_cols))
            ret.append(',\n    rhs_cols : ')
            ret.append(repr(self.rhs_cols))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("test_derivtype.Info_Values_Xinfo_Size_Array")
    class Info_Values_Xinfo_Size_Array(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=info_values_xinfo_size_array)
        
        
        Defined at test_deriv_types.f90 lines 21-28
        
        super-type
        Automatically generated to handle derived type arrays as a new derived type
        """
        def __init__(self, handle=None):
            """
            self = Info_Values_Xinfo_Size_Array()
            
            
            Defined at test_deriv_types.f90 lines 21-28
            
            
            Returns
            -------
            this : Info_Values_Xinfo_Size_Array
            	Object to be constructed
            
            
            Automatically generated constructor for info_values_xinfo_size_array
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _test_derivtype.f90wrap_info_values_xinfo_size_array_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Info_Values_Xinfo_Size_Array
            
            
            Defined at test_deriv_types.f90 lines 21-28
            
            Parameters
            ----------
            this : Info_Values_Xinfo_Size_Array
            	Object to be destructed
            
            
            Automatically generated destructor for info_values_xinfo_size_array
            """
            if self._alloc:
                _test_derivtype.f90wrap_info_values_xinfo_size_array_finalise(this=self._handle)
        
        def init_array_items(self):
            self.items = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _test_derivtype.f90wrap_info_values_xinfo_size_array__array_getitem__items,
                                            _test_derivtype.f90wrap_info_values_xinfo_size_array__array_setitem__items,
                                            _test_derivtype.f90wrap_info_values_xinfo_size_array__array_len__items,
                                            """
            Element items ftype=type(info_values) pytype=Info_Values
            
            
            Defined at  line 0
            
            """, Interpolation.info_values)
            return self.items
        
        _dt_array_initialisers = [init_array_items]
        
    
    @staticmethod
    def pack_values(deg, n_knots, nt, kl, ku, nleft, nright, lhs_rows, lhs_cols, \
        rhs_cols):
        """
        info = pack_values(deg, n_knots, nt, kl, ku, nleft, nright, lhs_rows, lhs_cols, \
            rhs_cols)
        
        
        Defined at test_deriv_types.f90 lines 31-49
        
        Parameters
        ----------
        deg : int
        n_knots : int
        nt : int
        kl : int
        ku : int
        nleft : int
        nright : int
        lhs_rows : int
        lhs_cols : int
        rhs_cols : int
        
        Returns
        -------
        info : Info_Values_Xinfo_Size_Array
        	super-type
        
        
        """
        info = _test_derivtype.f90wrap_pack_values(deg=deg, n_knots=n_knots, nt=nt, \
            kl=kl, ku=ku, nleft=nleft, nright=nright, lhs_rows=lhs_rows, \
            lhs_cols=lhs_cols, rhs_cols=rhs_cols)
        info = \
            f90wrap.runtime.lookup_class("test_derivtype.Info_Values_Xinfo_Size_Array").from_handle(info, \
            alloc=True)
        return info
    
    @staticmethod
    def unpack_values(self):
        """
        deg, n_knots, nt, kl, ku, nleft, nright, lhs_rows, lhs_cols, rhs_cols = \
            unpack_values(self)
        
        
        Defined at test_deriv_types.f90 lines 51-68
        
        Parameters
        ----------
        info : Info_Values
        
        Returns
        -------
        deg : int
        n_knots : int
        nt : int
        kl : int
        ku : int
        nleft : int
        nright : int
        lhs_rows : int
        lhs_cols : int
        rhs_cols : int
        
        """
        deg, n_knots, nt, kl, ku, nleft, nright, lhs_rows, lhs_cols, rhs_cols = \
            _test_derivtype.f90wrap_unpack_values(info=self._handle)
        return deg, n_knots, nt, kl, ku, nleft, nright, lhs_rows, lhs_cols, rhs_cols
    
    @property
    def sp(self):
        """
        Element sp ftype=integer pytype=int
        
        
        Defined at test_deriv_types.f90 line 17
        
        """
        return _test_derivtype.f90wrap_interpolation__get__sp()
    
    @property
    def dp(self):
        """
        Element dp ftype=integer pytype=int
        
        
        Defined at test_deriv_types.f90 line 18
        
        """
        return _test_derivtype.f90wrap_interpolation__get__dp()
    
    @property
    def qp(self):
        """
        Element qp ftype=integer pytype=int
        
        
        Defined at test_deriv_types.f90 line 19
        
        """
        return _test_derivtype.f90wrap_interpolation__get__qp()
    
    @property
    def info_size(self):
        """
        Element info_size ftype=integer pytype=int
        
        
        Defined at test_deriv_types.f90 line 20
        
        """
        return _test_derivtype.f90wrap_interpolation__get__info_size()
    
    def __str__(self):
        ret = ['<interpolation>{\n']
        ret.append('    sp : ')
        ret.append(repr(self.sp))
        ret.append(',\n    dp : ')
        ret.append(repr(self.dp))
        ret.append(',\n    qp : ')
        ret.append(repr(self.qp))
        ret.append(',\n    info_size : ')
        ret.append(repr(self.info_size))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

interpolation = Interpolation()

