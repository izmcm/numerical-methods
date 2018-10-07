import parser
from mpmath import *
import matplotlib.pyplot as plt

arr_t = []
arr_y = []

########################## EULER ##########################
def euler_method(y0, t0, h, n, fty):
    file_out = open("saida.txt", "a")
    file_out.write("Metodo de Euler\n")
    file_out.write("y( " + str(t0) + " ) = " + str(y0) + "\n")
    file_out.write("h = " + str(h) + "\n")

    expr = parser.expr(fty).compile()

    for i in range(0, n + 1):
        file_out.write(str(i) + " " + str(y0) + "\n")
        
        arr_t.append(t0)
        arr_y.append(y0)
        
        t = t0
        y = y0
        
        # yn+1 = yn + h*f(tn, yn)
        y0 += h*eval(expr)

        t0 += h
    
    file_out.write("\n")

    ## graphic
    plt.title("Metodo de Euler")
    plt.plot(arr_t, arr_y)
    plt.ylabel("y values")
    plt.xlabel("t values")
    plt.show()

########################## EULER INVERSO ##########################
def backward_euler_method(y0, t0, h, n, fty):
    file_out = open("saida.txt", "a")
    file_out.write("Metodo de Euler Inverso\n")
    file_out.write("y( " + str(t0) + " ) = " + str(y0) + "\n")
    file_out.write("h = " + str(h) + "\n")

    expr = parser.expr(fty).compile()

    for i in range(0, n + 1):
        file_out.write(str(i) + " " + str(y0) + "\n")
        
        arr_t.append(t0)
        arr_y.append(y0)
        
        t = t0
        y = y0
        
        # yn+1 = yn + h*f(tn, yn)
        y_euler = y0 + h*eval(expr) # descobre yn+1
        
        t0 += h
        t = t0
        y = y_euler
        
        # yn+1 = yn + h*f(tn+1, yn+1)
        y0 += h*eval(expr)
    
    file_out.write("\n")
    
    ## graphic
    plt.title("Metodo de Euler Inverso")
    plt.plot(arr_t, arr_y)
    plt.ylabel("y values")
    plt.xlabel("t values")
    plt.show()
    
########################## EULER APRIMORADO ##########################
def best_euler_method(y0, t0, h, n, fty):
    file_out = open("saida.txt", "a")
    file_out.write("Metodo de Euler Aprimorado\n")
    file_out.write("y( " + str(t0) + " ) = " + str(y0) + "\n")
    file_out.write("h = " + str(h) + "\n")

    expr = parser.expr(fty).compile()

    for i in range(0, n + 1):
        file_out.write(str(i) + " " + str(y0) + "\n")
        
        arr_t.append(t0)
        arr_y.append(y0)
        
        t = t0
        y = y0
        
        # f(tn, yn)
        f_tn_yn = eval(expr)
        
        # yn+1 = yn + h*f(tn, yn)
        y_euler = y0 + h*eval(expr) # descobre yn+1
        
        t0 += h
        t = t0
        y = y_euler
        
        # f(tn+1, yn+1)
        f_tn1_yn1 = eval(expr)
        
        # yn+1 = yn + (h/2)*(f(tn, yn) + f(tn+1, yn+1))
        y0 += (h/2)*(f_tn_yn + f_tn1_yn1)
    
    file_out.write("\n")

    ## graphic
    plt.title("Metodo de Euler Aprimorado")
    plt.plot(arr_t, arr_y)
    plt.ylabel("y values")
    plt.xlabel("t values")
    plt.show()

########################## MAIN ##########################
with open("entrada.txt") as file_in:
    for line in file_in:
        strin = line.split(' ')
        
        arr_y.clear()
        arr_t.clear()
        
        if strin[0] == "euler":
            y0 = strin[1]
            t0 = strin[2]
            h = strin[3]
            n = strin[4]
            fty = strin[5].split("\n ")
         
            euler_method(float(y0), float(t0), float(h), int(n), fty[0])
        
        if strin[0] == "euler_inverso":
            y0 = strin[1]
            t0 = strin[2]
            h = strin[3]
            n = strin[4]
            fty = strin[5].split("\n ")
         
            backward_euler_method(float(y0), float(t0), float(h), int(n), fty[0])
            
        if strin[0] == "euler_aprimorado":
            y0 = strin[1]
            t0 = strin[2]
            h = strin[3]
            n = strin[4]
            fty = strin[5].split("\n ")
         
            best_euler_method(float(y0), float(t0), float(h), int(n), fty[0])