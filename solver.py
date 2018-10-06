import parser
from mpmath import *

def euler_method(y0, t0, h, n, fxy):
    expr = parser.expr(fxy).compile()
    print(expr)
    for _ in range(0, n):
        #print("%.7f" %t0, " - %.7f" %y0)
        
        x = t0
        y = y0
        
        ans = y0 + h*eval(expr)
        print(ans)
        t0 = h + t0
        y0 = ans
        
        
with open("entrada.txt") as file:
    for line in file:
        strin = line.split(' ')
        
        if strin[0] == "euler":
            y0 = strin[1]
            t0 = strin[2]
            h = strin[3]
            n = strin[4]
            fxy = strin[5]
            
            euler_method(y0, t0, h, n, fxy)