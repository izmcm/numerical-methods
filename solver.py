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

########################## RUNGE KUTTA ##########################
def runge_kutta_method(y0, t0, h, n, fty):
    file_out = open("saida.txt", "a")
    file_out.write("Metodo de Runge-Kutta\n")
    file_out.write("y( " + str(t0) + " ) = " + str(y0) + "\n")
    file_out.write("h = " + str(h) + "\n")

    expr = parser.expr(fty).compile()

    for i in range(0, n + 1):
        file_out.write(str(i) + " " + str(y0) + "\n")
        
        arr_t.append(t0)
        arr_y.append(y0)
        
        ### Kn1 ###
        t = t0 # tn
        y = y0 # yn
        
        # f(tn, yn)
        kn1 = eval(expr)

        ### Kn2 ###
        t = t0 + h/2 # tn + h/2
        y = y0 + (h/2)*kn1 # yn + (h/2)*kn1

        # f(tn + h/2, yn + (h/2)*kn1)
        kn2 = eval(expr)

        ### Kn3 ###
        t = t0 + h/2 # tn + h/2
        y = y0 + (h/2)*kn2 # yn + (h/2)*kn2

        # f(tn + h/2, yn + (h/2)*kn2)
        kn3 = eval(expr)

        ### Kn4 ###
        t = t0 + h # tn + h
        y = y0 + h*kn3 # yn + h*kn3

        # f(tn + h, yn + h*kn3)
        kn4 = eval(expr)
        
        # yn+1 = yn + (h/6)*(kn1 + kn2 + kn3 + kn4)
        y0 += (h/6)*(kn1 + 2*kn2 + 2*kn3 + kn4)
        t0 += h
    
    file_out.write("\n")

    ## graphic
    plt.title("Metodo de Runge-Kutta")
    plt.plot(arr_t, arr_y)
    plt.ylabel("y values")
    plt.xlabel("t values")
    plt.show()

########################## ADAMS BASHFORTH ##########################
def adams_bashforth_method(y0, t0, h, n, fty, order):
    file_out = open("saida.txt", "a")
    file_out.write("Metodo de Adams Bashforth\n")
    file_out.write("y( " + str(t0) + " ) = " + str(y0) + "\n")
    file_out.write("h = " + str(h) + "\n")

    expr = parser.expr(fty).compile()

    # yn+1 = yn + h*((3/2)*f(tn, yn) - (1/2)*f(tn-1, yn-1))
    if order == 2:
    	# calculando pontos iniciais com euler
	    for i in range(0, 3):
	        file_out.write(str(i) + " " + str(y0) + "\n")

	        arr_t.append(t0)
	        arr_y.append(y0)
	        
	        t = t0
	        y = y0
	        
	        # yn+1 = yn + h*f(tn, yn) - euler
	        y0 += h*eval(expr)
	        t0 += h

	    # terminando com adams bashforth de ordem 2
	    for i in range(3, n + 1):
	        file_out.write(str(i) + " " + str(y0) + "\n")

	        arr_t.append(t0)
	        arr_y.append(y0)

	        # f(tn, yn)
	        t = t0
	        y = y0

	        f_tn_yn = eval(expr)

	        # f(tn-1, yn-1)
	        t = arr_t[len(arr_t) - 1]
	        y = arr_y[len(arr_t) - 1]

	        f_tn1_yn1 = eval(expr)

	        y0 += h*((3/2)*f_tn_yn - (1/2)*f_tn1_yn1)
	        t0 += h

    file_out.write("\n")

    ## graphic
    plt.title("Metodo de Euler Aprimorado")
    plt.plot(arr_t, arr_y)
    #plt.ylabel("y values")
    #plt.xlabel("t values")
    #plt.show()

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
        
        elif strin[0] == "euler_inverso":
            y0 = strin[1]
            t0 = strin[2]
            h = strin[3]
            n = strin[4]
            fty = strin[5].split("\n ")
         
            backward_euler_method(float(y0), float(t0), float(h), int(n), fty[0])
            
        elif strin[0] == "euler_aprimorado":
            y0 = strin[1]
            t0 = strin[2]
            h = strin[3]
            n = strin[4]
            fty = strin[5].split("\n ")
         
            best_euler_method(float(y0), float(t0), float(h), int(n), fty[0])

        elif strin[0] == "runge_kutta":
            y0 = strin[1]
            t0 = strin[2]
            h = strin[3]
            n = strin[4]
            fty = strin[5].split("\n ")
         
            runge_kutta_method(float(y0), float(t0), float(h), int(n), fty[0])

        elif strin[0] == "adam_bashforth_by_euler":
        	y0 = strin[1]
        	t0 = strin[2]
        	h = strin[3]
        	n = strin[4]
        	fty = strin[5]
        	order = strin[6].split("\n ")
         
        	adams_bashforth_method(float(y0), float(t0), float(h), int(n), fty, int(order[0]))







