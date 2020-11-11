import sympy
import numpy as np
from tleedmlib.leedbase import reduceUnitCell

def leftmult2(m, i0, i1, a, b, c, d):
    for j in range(m.cols):
        x, y = m[i0, j], m[i1, j]
        m[i0, j] = a * x + b * y
        m[i1, j] = c * x + d * y

def rightmult2(m, j0, j1, a, b, c, d):
    for i in range(m.rows):
        x, y = m[i, j0], m[i, j1]
        m[i, j0] = a * x + c * y
        m[i, j1] = b * x + d * y

def smithNormalForm(m, domain=sympy.ZZ):
    m = sympy.Matrix(m)
    s = sympy.eye(m.rows)
    t = sympy.eye(m.cols)
    last_j = -1
    for i in range(m.rows):
        for j in range(last_j+1, m.cols):
            if not m.col(j).is_zero:
                break
        else:
            break
        if m[i,j] == 0:
            for ii in range(m.rows):
                if m[ii,j] != 0:
                    break
            leftmult2(m, i, ii, 0, 1, 1, 0)
            rightmult2(s, i, ii, 0, 1, 1, 0)
        rightmult2(m, j, i, 0, 1, 1, 0)
        leftmult2(t, j, i, 0, 1, 1, 0)
        j = i
        upd = True
        while upd:
            upd = False
            for ii in range(i+1, m.rows):
                if m[ii, j] == 0:
                    continue
                upd = True
                if domain.rem(m[ii, j], m[i, j]) != 0:
                    coef1, coef2, g = domain.gcdex(int(m[i,j]), int(m[ii, j]))
                    coef3 = domain.quo(m[ii, j], g)
                    coef4 = domain.quo(m[i, j], g)
                    leftmult2(m, i, ii, coef1, coef2, -coef3, coef4)
                    rightmult2(s, i, ii, coef4, -coef2, coef3, coef1)
                coef5 = domain.quo(m[ii, j], m[i, j])
                leftmult2(m, i, ii, 1, 0, -coef5, 1)
                rightmult2(s, i, ii, 1, 0, coef5, 1)
            for jj in range(j+1, m.cols):
                if m[i, jj] == 0:
                    continue
                upd = True
                if domain.rem(m[i, jj], m[i, j]) != 0:
                    coef1, coef2, g = domain.gcdex(int(m[i,j]), int(m[i, jj]))
                    coef3 = domain.quo(m[i, jj], g)
                    coef4 = domain.quo(m[i, j], g)
                    rightmult2(m, j, jj, coef1, -coef3, coef2, coef4)
                    leftmult2(t, j, jj, coef4, coef3, -coef2, coef1)
                coef5 = domain.quo(m[i, jj], m[i, j])
                rightmult2(m, j, jj, 1, -coef5, 0, 1)
                leftmult2(t, j, jj, 1, coef5, 0, 1)
        last_j = j
    for i1 in range(min(m.rows, m.cols)):
        for i0 in reversed(range(i1)):
            coef1, coef2, g = domain.gcdex(int(m[i0, i0]), int(m[i1,i1]))
            if g == 0:
                continue
            coef3 = domain.quo(m[i1, i1], g)
            coef4 = domain.quo(m[i0, i0], g)
            leftmult2(m, i0, i1, 1, coef2, coef3, coef2*coef3-1)
            rightmult2(s, i0, i1, 1-coef2*coef3, coef2, coef3, -1)
            rightmult2(m, i0, i1, coef1, 1-coef1*coef4, 1, -coef4)
            leftmult2(t, i0, i1, coef4, 1-coef1*coef4, 1, -coef1)
    return (s, m, t)

def getSmallestSupercell(sl1, sl2):
    m1 = sympy.Matrix(sl1)
    m2 = sympy.Matrix(sl2)
    if all([(abs(v) >= 1 or v == 0) for v in m1*m2.inv()]):
        return np.array(m1).astype(int)
    if all([(abs(v) >= 1 or v == 0) for v in m2*m1.inv()]):
        return np.array(m2).astype(int)
    at = 0
    while at < 2:
        a = m2.inv()*m1
        print("a = {}".format(a))
        n = 1 / min(min([abs(v) for v in a if v != 0]), 1)
        n = n / sympy.igcd(*[v for v in a*n])
        a *= n      # smallest integer form of a
        (b, s, t) = smithNormalForm(a)
        print("s = {}".format(s))
        if abs(s[1,1]) > abs(s[0,0]):
            break
        else:
            m1 = sympy.Matrix([[1,0],[0,-1]])*m1
            at += 1
            print("reversing..")
    r = sympy.diag(*[n/np.gcd(s[i,i],n) for i in (0,1)])
    sc = m1*t.inv()*r
    ab, _, _ = reduceUnitCell(np.array(sc).astype(int))
    return ab
    
sl1 = np.array([[1,1],[1,-1]])
sl2 = np.array([[1,2],[1,-2]])

ab = getSmallestSupercell(sl1,sl2)
print(ab)
