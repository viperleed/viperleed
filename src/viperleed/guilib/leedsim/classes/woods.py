"""
======================================
  ViPErLEED Graphical User Interface
======================================
 *** module guilib.leedsim.classes.woods ***

Defines the Woods class, used for conversion of a Woods' notation string into
a matrix and vice versa, as well as for formatting

Author: Michele Riva
Created: 2021-03-13
"""

import re
import itertools

import numpy as np


# TODO: fix method names
# TODO: allow also "rect" special Woods notation for hexagonal lattices (ONLY??)

degrees = '\u00b0'

class Woods:
    common = {'p(1\u00d71)', 'p(1\u00d72)', 'p(2\u00d71)', 'p(2\u00d72)',
              'p(3\u00d71)', 'p(3\u00d72)', 'p(3\u00d73)', 'c(2\u00d72)',
              'c(4\u00d74)', 'c(6\u00d72)', 'c(8\u00d72)'}
    examples = {
        'Oblique': common,
        'Square': common | {'c(4\u00d72)',
                            'p(\u221a2\u00d7\u221a2)R45' + degrees,
                            'p(\u221a5\u00d7\u221a5)R26.6' + degrees,
                            'p(2\u221a2\u00d7\u221a2)R45' + degrees,
                            'c(3\u221a2\u00d7\u221a2)R45' + degrees,
                            'c(5\u221a2\u00d7\u221a2)R45' + degrees},
        'Rectangular': common,
        'Hexagonal': common | {'c(4\u00d72)',
                               'p(\u221a3\u00d7\u221a3)R30' + degrees,
                               'p(\u221a7\u00d7\u221a7)R19.1' + degrees,
                               'p(2\u221a3\u00d72\u221a3)R30' + degrees},
        'Rhombic': common,
        }
    
    def woodsToMatrix(self, woods, bulk):
        parsed = self.parseWoods(woods)
        if parsed is None:
            return None
        g1 = parsed[1]
        g2 = parsed[2]
        alpha = np.radians(parsed[3])
        R1 = np.linalg.norm(bulk[0])
        R2 = np.linalg.norm(bulk[1])
        q = R2/R1
        omega = np.arccos(np.dot(bulk[0], bulk[1])/(R1*R2))
        
        m = np.array([[g1*np.sin(omega-alpha), g1*np.sin(alpha)/q],
                      [-g2*q*np.sin(alpha), g2*np.sin(omega+alpha)]])
        m = m/np.sin(omega)
        
        if parsed[0] == 'c':
            m = np.dot([[1, 1], [-1, 1]], m)/2
        
        if self.isCommensurate(m):
            m = np.array([int(np.round(mij)) 
                          for mij in m.ravel()]).reshape(m.shape)
        else:
            m = None  # incommensurate
        return m
    
    def parseWoods(self, woods):
        if not isinstance(woods, str):
            return None
        
        reWoods = re.compile(
            r'''^(?P<prefix>[pc])                # * primitive or centered
            \(                                   # * open parenthesis
            (?P<g1>\d+)?                         # * direction1, integer part
            (\u221a(?P<g1rt>\d+))?               # * direction1, radical part
            \u00d7                               # * times
            (?P<g2>\d+)?                         # * direction2, integer part
            (\u221a(?P<g2rt>\d+))?               # * direction2, radical part
            \)                                   # * close parenthesis
            ((R(?P<alpha>\d+(\.\d+)?))\u00b0)?$  # * rotation angle
            ''', re.VERBOSE)
        # notice that the ^ and $ anchors make sure that the string is 
        # matched as a whole
        
        # check if it matches the full woods
        m = reWoods.match(woods)
        if m is not None:
            w = m.groupdict()
            g1 = 1
            g2 = 1
            if w['g1'] is not None:
                g1 *= float(w['g1'])
            if w['g1rt'] is not None:
                g1 *= np.sqrt(float(w['g1rt']))
            if w['g2'] is not None:
                g2 *= float(w['g2'])
            if w['g2rt'] is not None:
                g2 *= np.sqrt(float(w['g2rt']))
            if w['alpha'] is not None:
                alpha = float(w['alpha'])
            else:
                alpha = 0
            return [w['prefix'], g1, g2, alpha]
        return None
    
    def matrixToWoods(self, m, bulk):
        m = np.array(m)
        
        pc = self.primitiveOrCentered(m, bulk)
        if pc:  # matrix is woods representable as primitive or centered
            toFormat = []
            for gSq in [pc[1]**2, pc[2]**2]:
                gSq = np.round(gSq)
                (g2, grt) = self.squareToProdOfSquares(gSq)
                integer = str(int(round(np.sqrt(g2))))
                
                dirStr = '' #format the direction in here
                if grt > 1: # insert root part
                    dirStr = '\u221a' + str(int(grt))
                if dirStr == '': # if there is no root part, always insert 
                                 # the integer part
                    dirStr = integer
                else:
                    if integer != '1': # otherwise place it in only if it's 
                                       # not 1
                        dirStr = integer + dirStr
                
                toFormat.append(dirStr)
            woods = '{}({})'.format(pc[0], '\u00d7'.join(toFormat))
            cosAlpha = pc[3]
            if np.abs(cosAlpha) > 1e-3 and 1-np.abs(cosAlpha) > 1e-3:
                #angle not 0, 90 nor 180
                alpha = np.round(np.degrees(np.arccos(cosAlpha)),
                                 decimals = 1)
                woods += 'R' + str(alpha) + degrees
            return woods
        return None
    
    def isCommensurate(self, m, eps=1e-3):
        if m is None:
            return False
        
        mu = np.linalg.det(m)
        if np.round(mu) == 0.0:  # matrix is singular
            return False

        if np.abs(mu/np.round(mu) - 1) > eps:  # determinant is not integer
            return False
        
        # now check whether any element is non-integer
        for mij in m.ravel():
            if np.round(mij) == 0.0:
                if abs(mij)/np.sqrt(np.abs(mu)) > eps:
                    return False
            elif np.abs(mij/np.round(mij) - 1) > eps:
                return False
        return True
    
    def isRepresentable(self, m, basis):
        transform = np.dot(m, basis)
        nBasis = np.linalg.norm(basis, axis=1)
        nTransf = np.linalg.norm(transform, axis=1)
        g1 = nTransf[0]/nBasis[0]
        g2 = nTransf[1]/nBasis[1]
        mu = np.abs(np.linalg.det(m))
        
        return abs(mu/(g1*g2) - 1) < 1e-8
    
    def primitiveOrCentered(self, m, basis):
        invT = np.array([[1, -1], [1, 1]])
        
        primitive = self.isRepresentable(m, basis)
        ctrd = self.isRepresentable(np.dot(invT, m), basis)
        if primitive:
            prefix = 'p'
        elif ctrd:
            prefix = 'c'
            m = np.dot(invT, m)
        else:
            return False
        
        nBasis = np.linalg.norm(basis, axis=1)
        transform = np.dot(m, basis)
        nTransf = np.linalg.norm(transform, axis=1)
        g1 = nTransf[0]/nBasis[0]
        g2 = nTransf[1]/nBasis[1]
        cosAlpha = np.dot(transform[0], basis[0])/(nTransf[0]*nBasis[0])
        
        return (prefix, g1, g2, cosAlpha)
    
    def squareToProdOfSquares(self, number):
        # takes number, finds all prime factors, returns a tuple, first 
        # element is a product of all primes showing up an even number of 
        # times, the second one the rest. Useful to turn, e.g., sqrt(12) 
        # into 2sqrt(3)
        
        factors = list(self.primeFactors(number))
        if not factors:
            factors = [1]
        uniqueFact = sorted(list(set(factors)))
        countFact = [((factors.count(fac) // 2)*2, factors.count(fac) % 2)
                      for fac in uniqueFact]
        (pow2, restPow) = zip(*countFact)
        (squares, remainders) = zip(*[(fact**pow, fact**rem)
                                      for (fact, pow, rem)
                                      in zip(uniqueFact, pow2, restPow)
                                      ])
        
        return (np.prod(squares), np.prod(remainders))
    
    def primeFactors(self, n):  # TODO: get a reference for this code
        f = 2
        increments = itertools.chain(
                                  [1, 2, 2],
                                  itertools.cycle([4, 2, 4, 2, 4, 6, 2, 6])
                                  )
        for incr in increments:
            if f*f > n:
                break
            while n % f == 0:
                yield f
                n //= f
            f += incr
        if n > 1:
            yield n
