import math
from math import nan
import matplotlib.pyplot as plt
import numpy as np
from sympy import *

Y0 = Symbol('Y(x_0)')
y0 = Symbol('y(x_0)')
y1 = Symbol("y'(x_0)")
y2 = Symbol("y''(x_0)")
y3 = Symbol("y'''(x_0)")
y4 = Symbol("y''''(x_0)")
h = Symbol('h')
paso = 3/2 *h
y_c = Y0 \
    + paso*y0 \
    + (paso**2)/math.factorial(2)*y1 \
    + (paso**3)/math.factorial(3)*y2 \
    + (paso**4)/math.factorial(4)*y3 

paso2 = 1/2 * paso
y_c_2 = Y0 \
    + paso2*y0 \
    + (paso2**2)/math.factorial(2)*y1 \
    + (paso2**3)/math.factorial(3)*y2 \
    + (paso2**4)/math.factorial(4)*y3 
e = y_c - y_c_2
print(e)