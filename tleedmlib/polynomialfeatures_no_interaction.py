# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 10:22:22 2021

@author: Florian Kraushofer
"""

import numpy as np
from scipy import sparse

from sklearn.preprocessing import PolynomialFeatures
from sklearn.utils.validation import (check_is_fitted,
                                      _deprecate_positional_args,
                                      FLOAT_DTYPES, check_array)
from sklearn.preprocessing._csr_polynomial_expansion import (
                                                    _csr_polynomial_expansion)


class PolyFeatNoMix(PolynomialFeatures):
    """Child of the PolynomialFeatures class which excludes mixing terms
     like x1*x2, leaving only terms of type x1, x1^2, etc. and a constant bias
     term - i.e. basically the opposite of the 'interaction_only' keyword
     argument of PolynomialFeatures. See documentation for
     sklearn.preprocessing.PolynomialFeatures for details.
    ----------
    degree : int, default=2
        The degree of the polynomial features.
    include_bias : bool, default=True
        If True (default), then include a bias column, the feature in which
        all polynomial powers are zero (i.e. a column of ones - acts as an
        intercept term in a linear model).
    order : {'C', 'F'}, default='C'
        Order of output array in the dense case. 'F' order is faster to
        compute, but may slow down subsequent estimators.
        .. versionadded:: 0.21
    Attributes
    ----------
    powers_ : ndarray of shape (n_output_features, n_input_features)
        powers_[i, j] is the exponent of the jth input in the ith output.
    n_input_features_ : int
        The total number of input features.
    n_output_features_ : int
        The total number of polynomial output features. The number of output
        features is computed by iterating over all suitably sized combinations
        of input features.
    Notes
    -----
    Be aware that the number of features in the output array scales
    polynomially in the number of features of the input array, and
    exponentially in the degree. High degrees can cause overfitting.
    See :ref:`examples/linear_model/plot_polynomial_interpolation.py
    <sphx_glr_auto_examples_linear_model_plot_polynomial_interpolation.py>`
    """
    @_deprecate_positional_args
    def __init__(self, degree=2, *, include_bias=True, order='C'):
        self.degree = degree
        self.interaction_only = False
        self.include_bias = include_bias
        self.order = order

    @staticmethod
    def _combinations(n_features, degree, interaction_only, include_bias):
        # return ((i,)*j for j in range(n_features)
        #         for i in range(degree+1) if j >= i)
        return ((tuple(),) + tuple((i,)*j for j in range(1, degree+1)
                                          for i in range(n_features)))

    def transform(self, X):
        """Transform data to polynomial features
        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The data to transform, row by row.
            Prefer CSR over CSC for sparse input (for speed), but CSC is
            required if the degree is 4 or higher. If the degree is less than
            4 and the input format is CSC, it will be converted to CSR, have
            its polynomial features generated, then converted back to CSC.
            If the degree is 2 or 3, the method described in "Leveraging
            Sparsity to Speed Up Polynomial Feature Expansions of CSR Matrices
            Using K-Simplex Numbers" by Andrew Nystrom and John Hughes is
            used, which is much faster than the method used on CSC input. For
            this reason, a CSC input will be converted to CSR, and the output
            will be converted back to CSC prior to being returned, hence the
            preference of CSR.
        Returns
        -------
        XP : {ndarray, sparse matrix} of shape (n_samples, NP)
            The matrix of features, where NP is the number of polynomial
            features generated from the combination of inputs. If a sparse
            matrix is provided, it will be converted into a sparse
            ``csr_matrix``.
        """
        check_is_fitted(self)

        X = check_array(X, order='F', dtype=FLOAT_DTYPES,
                        accept_sparse=('csr', 'csc'))

        n_samples, n_features = X.shape

        if n_features != self.n_input_features_:
            raise ValueError("X shape does not match training shape")

        if sparse.isspmatrix_csr(X):
            if self.degree > 3:
                return self.transform(X.tocsc()).tocsr()
            to_stack = []
            if self.include_bias:
                to_stack.append(np.ones(shape=(n_samples, 1), dtype=X.dtype))
            to_stack.append(X)
            for deg in range(2, self.degree+1):
                Xp_next = _csr_polynomial_expansion(X.data, X.indices,
                                                    X.indptr, X.shape[1],
                                                    self.interaction_only,
                                                    deg)
                if Xp_next is None:
                    break
                to_stack.append(Xp_next)
            XP = sparse.hstack(to_stack, format='csr')
        elif sparse.isspmatrix_csc(X) and self.degree < 4:
            return self.transform(X.tocsr()).tocsc()
        else:
            if sparse.isspmatrix(X):
                combinations = self._combinations(n_features, self.degree,
                                                  self.interaction_only,
                                                  self.include_bias)
                columns = []
                for comb in combinations:
                    if comb:
                        out_col = 1
                        for col_idx in comb:
                            out_col = X[:, col_idx].multiply(out_col)
                        columns.append(out_col)
                    else:
                        bias = sparse.csc_matrix(np.ones((X.shape[0], 1)))
                        columns.append(bias)
                XP = sparse.hstack(columns, dtype=X.dtype).tocsc()
            else:
                XP = np.empty((n_samples, self.n_output_features_),
                              dtype=X.dtype, order=self.order)

                # What follows is a faster implementation of:
                # for i, comb in enumerate(combinations):
                #     XP[:, i] = X[:, comb].prod(1)
                # Second optimisation happens for degrees >= 3.
                # Xi^3 is computed reusing previous computation:
                # Xi^3 = Xi^2 * Xi.

                if self.include_bias:
                    XP[:, 0] = 1
                    base_col = 1
                else:
                    base_col = 0

                XP[:, base_col:base_col + n_features] = X
                for d in range(1, self.degree):
                    start = base_col + (d-1) * n_features
                    end = base_col + d * n_features
                    np.multiply(XP[:, start:end], X,
                                out=XP[:, start+n_features:end+n_features],
                                casting='no')

        return XP