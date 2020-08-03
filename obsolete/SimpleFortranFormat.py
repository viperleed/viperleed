# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 11:59:13 2019

@author: FF

Contains basic fortran formatting routines. NOT FINISHED, AND NOT IN USE
"""

def Ax(s, x):
    """Takes a string and a field width; if the string fits, it will be padded
    to the field width with spaces on the left. If it doesn't fit, returns ***"""
    if len(s) > x:
        o = '*' * x
    else:
        o = s.rjust(x)
    return o

def Ix(i, x):
    """Takes an integer and a field width; if the integer fits, it will be padded
    to the field width with spaces on the left. If it doesn't fit, returns ***"""
    return Ax(str(i),x)
#
#def Fx_y(f,x,y):
#    """Takes a float, a total field width x and a number of digits y that should 
#    be printed after the comma. If the float fits, it will be padded to the 
#    field width with spaces on the left. If it doesn't fit, returns ***"""
    

def LAx(llist, x, elperline=0):
    """Generates a line of output from a list of strings and a field width per 
    string; for each string, if it fits, it will be padded to the field width 
    with spaces on the left. If it doesn't fit, *** are inserted instead. If
    elperline is given, a linebreak will be inserted after elperline elements."""
    count = 0
    out = ''
    for s in llist:
        if elperline > 0 and count >= elperline: 
            out += '\n'
            count = 0
        out += Ax(s,x)
        count += 1
    return out

def LIx(llist, x, elperline=0):
    """Generates a line of output from a list of integers and a field width per 
    integer; for each integer, if it fits, it will be padded to the field width 
    with spaces on the left. If it doesn't fit, *** are inserted instead. If
    elperline is given, a linebreak will be inserted after elperline elements."""
    slist = []
    for i in llist:
        slist.append(str(i))
    return LAx(slist, x, elperline)

#out = LAx(['f','as', 'b', 'c', '12', 'd'], 3)
llist = '   4   6   9  10  11  12  14  17  18  21  22  25  28  31  32 34  37  38  41  42  45  46  47  48  51  52  54  57  58  61 62  65  66  69'.split()
out = LIx(llist,4, 15)
#print("'"+out+"'")
print(out)