import parser
from mpmath import *
import matplotlib.pyplot as plt

arr_t = []
arr_y = []

################################################################################### GERAL ###################################################################################
def show_graphic():
    plt.plot(arr_t, arr_y)
    plt.ylabel("y values")
    plt.xlabel("t values")
    plt.show()

def do_eval(t, y, expr):
    return eval(expr)

def print_in_file(method, h):
    file_out = open("saida.txt", "a")
    file_out.write("Metodo de " + str(method) + "\n")
    file_out.write("y( " + str(arr_t[0]) + " ) = " + str(arr_y[0]) + "\n")
    file_out.write("h = " + str(h) + "\n")

    for i in range(0, len(arr_y)):
        file_out.write(str(i) + " " + str(arr_y[i]) + "\n")

    file_out.write("\n")

################################################################################### EULER ###################################################################################
def euler_method(y0, t0, h, n, fty, flag_printar):
    expr = parser.expr(fty).compile()

    for i in range(0, n + 1):        
        arr_t.append(t0)
        arr_y.append(y0)

        # yn+1 = yn + h*f(tn, yn)
        y0 += h*do_eval(t0, y0, expr)
        t0 += h

    if flag_printar == 1:
        print_in_file("Euler", h)

        ## graphic
        plt.title("Metodo de Euler")
        show_graphic()

    else:
        return t0

############################################################################## EULER INVERSO ##############################################################################
def backward_euler_method(y0, t0, h, n, fty, flag_printar):
	expr = parser.expr(fty).compile()

	for i in range(0, n + 1):        
		arr_t.append(t0)
		arr_y.append(y0)
        
        # yn+1 = yn + h*f(tn, yn)
		y_euler = y0 + h*do_eval(t0, y0, expr) # descobre yn+1
        
		t0 += h
        
		# yn+1 = yn + h*f(tn+1, yn+1)
		y0 += h*do_eval(t0, y_euler, expr)
    
	if flag_printar == 1:
		print_in_file("Euler Inverso", h)

		## graphic
		plt.title("Metodo de Euler Inverso")
		show_graphic()
	
	else:
		return t0
    
############################################################################# EULER APRIMORADO #############################################################################
def best_euler_method(y0, t0, h, n, fty, flag_printar):
	expr = parser.expr(fty).compile()

	for i in range(0, n + 1):        
		arr_t.append(t0)
		arr_y.append(y0)
        
        # f(tn, yn)
		f_tn_yn = do_eval(t0, y0, expr)
        
        # yn+1 = yn + h*f(tn, yn)
		y_euler = y0 + h*f_tn_yn # descobre yn+1
        
		t0 += h
        
        # f(tn+1, yn+1)
		f_tn1_yn1 = do_eval(t0, y_euler, expr)
        
        # yn+1 = yn + (h/2)*(f(tn, yn) + f(tn+1, yn+1))
		y0 += (h/2)*(f_tn_yn + f_tn1_yn1)
	    
	if flag_printar == 1:
		print_in_file("Euler Aprimorado", h)

	    ## graphic
		plt.title("Metodo de Euler Aprimorado")
		show_graphic()

	else:
		return t0

############################################################################## RUNGE KUTTA ##############################################################################
def runge_kutta_method(y0, t0, h, n, fty, flag_printar):
	expr = parser.expr(fty).compile()

	for i in range(0, n + 1):        
		arr_t.append(t0)
		arr_y.append(y0)
        
        ### Kn1 ###
        # f(tn, yn)
		kn1 = do_eval(t0, y0, expr)

        ### Kn2 ###
        # f(tn + h/2, yn + (h/2)*kn1)
		kn2 = do_eval(t0 + h/2, y0 + (h/2)*kn1, expr)

        ### Kn3 ###
        # f(tn + h/2, yn + (h/2)*kn2)
		kn3 = do_eval(t0 + h/2, y0 + (h/2)*kn2, expr)

        ### Kn4 ###
        # f(tn + h, yn + h*kn3)
		kn4 = do_eval(t0 + h, y0 + h*kn3, expr)
        
        # yn+1 = yn + (h/6)*(kn1 + kn2 + kn3 + kn4)
		y0 += (h/6)*(kn1 + 2*kn2 + 2*kn3 + kn4)
		t0 += h

	if flag_printar == 1:
	    print_in_file("Runge-Kutta", h)

	    ## graphic
	    plt.title("Metodo de Runge-Kutta")
	    show_graphic()

	else:
		return t0

############################################################################## ADAMS BASHFORTH ##############################################################################
def adams_bashforth_method(t0, h, n, fty, order, method):
    
    # yn+1 = yn + h*((3/2)*f(tn, yn) - (1/2)*f(tn-1, yn-1))
    if order == 2:
        ab_2(t0, h, n, fty, method)

    # yn+1 = yn + h*((23/12)*f(tn, yn) - (4/3)*f(tn-1, yn-1) + (5/12)*f(tn-2, yn-2))
    elif order == 3:
        ab_3(t0, h, n, fty, method)

    # yn+1 = yn + h*((55/24)*f(tn, yn) - (59/24)*f(tn-1, yn-1) + (37/24)*f(tn-2, yn-2) - (3/8)*f(tn-3, yn-3))
    elif order == 4:
        ab_4(t0, h, n, fty, method)

    # yn+1 = yn + h*((1901/720)*f(tn, yn) - (1387/360)*f(tn-1, yn-1) + (109/30)*f(tn-2, yn-2) - (637/360)*f(tn-3, yn-3) + (251/720)*f(tn-4, yn-4))
    elif order == 5:
        ab_5(t0, h, n, fty, method)

    # yn+1 = yn + h*((4277/1440)*f(tn, yn) - (2641/480)*f(tn-1, yn-1) + (4991/720)*f(tn-2, yn-2) - (3649/720)*f(tn-3, yn-3) + (959/480)*f(tn-4, yn-4) - (95/288)*f(tn-5, yn-5))
    elif order == 6:
        ab_6(t0, h, n, fty, method)

    # yn+1 = yn + h*((198721/60480)*f(tn, yn) - (18637/2520)*f(tn-1, yn-1) + (235183/20160)*f(tn-2, yn-2) - (10754/945)*f(tn-3, yn-3) + (135713/20160)*f(tn-4, yn-4) - (5603/2520)*f(tn-5, yn-5) + (19087/60480)*f(tn-6, yn-6))
    elif order == 7:
        ab_7(t0, h, n, fty, method)


####################################################################### ADAMS BASHFORTH - ORDENS 1~8 ########################################################################
def ab_8(t0, h, n, fty, method):
    expr = parser.expr(fty).compile()

    y0 = arr_y[len(arr_y) - 1]

    for i in range(len(arr_y), n + 1):
    	arr_t.append(t0)

        # f(tn, yn)
    	f_tn_yn = do_eval(arr_t[len(arr_y) - 1], y0, expr)

        # f(tn-1, yn-1)
    	f_tn1_yn1 = do_eval(arr_t[len(arr_y) - 2], arr_y[len(arr_y) - 2], expr)

        # f(tn-2, yn-2)
    	f_tn2_yn2 = do_eval(arr_t[len(arr_y) - 3], arr_y[len(arr_y) - 3], expr)

        # f(tn-3, yn-3)
    	f_tn3_yn3 = do_eval(arr_t[len(arr_y) - 4], arr_y[len(arr_y) - 4], expr)

        # f(tn-4, yn-4)
    	f_tn4_yn4 = do_eval(arr_t[len(arr_y) - 5], arr_y[len(arr_y) - 5], expr)

        # f(tn-5, yn-5)
    	f_tn5_yn5 = do_eval(arr_t[len(arr_y) - 6], arr_y[len(arr_y) - 6], expr)

        # f(tn-6, yn-6)
    	f_tn6_yn6 = do_eval(arr_t[len(arr_y) - 7], arr_y[len(arr_y) - 7], expr)

    	# f(tn-6, yn-6)
    	f_tn7_yn7 = do_eval(arr_t[len(arr_y) - 8], arr_y[len(arr_y) - 8], expr)

		#yn+1 = yn + h*((16083/4480)*f(tn, yn) - (1152169/120960)*f(tn-1, yn-1) + (242653/13440)*f(tn-2, yn-2) - (296053/13440)*f(tn-3, yn-3) + (2102243/120960)*f(tn-4, yn-4) - (115747/13440)*f(tn-5, yn-5) + (32863/13440)*f(tn-6, yn-6) - (5257/17280)*f(tn-7, yn-7))
    	y0 += h*((16083/4480)*f_tn_yn - (1152169/120960)*f_tn1_yn1 + (242653/13440)*f_tn2_yn2 - (296053/13440)*f_tn3_yn3 + (2102243/120960)*f_tn4_yn4 - (115747/13440)*f_tn5_yn5 + (32863/13440)*f_tn6_yn6 - (5257/17280)*f_tn7_yn7)
    	t0 += h

    	arr_y.append(y0)

    if method == "nothing":
        print_in_file("Adams Bashforth de 8 Ordem", h)
        plt.title("Metodo de Adams Bashforth de 8 Ordem")
    
    elif method == "Euler":
        print_in_file("Adams Bashforth de 8 Ordem por Euler", h)
        plt.title("Metodo de Adams Bashforth de 8 Ordem por Euler")

    elif method == "Euler Inverso":
        print_in_file("Adams Bashforth de 8 Ordem por Euler Inverso", h)
        plt.title("Metodo de Adams Bashforth de 8 Ordem por Euler Inverso")

    elif method == "Euler Aprimorado":
        print_in_file("Adams Bashforth de 8 Ordem por Euler Aprimorado", h)
        plt.title("Metodo de Adams Bashforth de 8 Ordem por Euler Aprimorado")

    elif method == "Runge-Kutta":
        print_in_file("Adams Bashforth de 8 Ordem por Runge-Kutta", h)
        plt.title("Metodo de Adams Bashforth de 8 Ordem por Runge-Kutta")

    show_graphic()

def ab_7(t0, h, n, fty, method):
    expr = parser.expr(fty).compile()

    y0 = arr_y[len(arr_y) - 1]

    for i in range(len(arr_y), n + 1):
    	arr_t.append(t0)

        # f(tn, yn)
    	f_tn_yn = do_eval(arr_t[len(arr_y) - 1], y0, expr)

        # f(tn-1, yn-1)
    	f_tn1_yn1 = do_eval(arr_t[len(arr_y) - 2], arr_y[len(arr_y) - 2], expr)

        # f(tn-2, yn-2)
    	f_tn2_yn2 = do_eval(arr_t[len(arr_y) - 3], arr_y[len(arr_y) - 3], expr)

        # f(tn-3, yn-3)
    	f_tn3_yn3 = do_eval(arr_t[len(arr_y) - 4], arr_y[len(arr_y) - 4], expr)

        # f(tn-4, yn-4)
    	f_tn4_yn4 = do_eval(arr_t[len(arr_y) - 5], arr_y[len(arr_y) - 5], expr)

        # f(tn-5, yn-5)
    	f_tn5_yn5 = do_eval(arr_t[len(arr_y) - 6], arr_y[len(arr_y) - 6], expr)

        # f(tn-6, yn-6)
    	f_tn6_yn6 = do_eval(arr_t[len(arr_y) - 7], arr_y[len(arr_y) - 7], expr)

        # yn+1 = yn + h*((198721/60480)*f(tn, yn) - (18637/2520)*f(tn-1, yn-1) + (235183/20160)*f(tn-2, yn-2) - (10754/945)*f(tn-3, yn-3) + (135713/20160)*f(tn-4, yn-4) - (5603/2520)*f(tn-5, yn-5) + (19087/60480)*f(tn-6, yn-6))
    	y0 += h*((198721/60480)*f_tn_yn - (18637/2520)*f_tn1_yn1 + (235183/20160)*f_tn2_yn2 - (10754/945)*f_tn3_yn3 + (135713/20160)*f_tn4_yn4 - (5603/2520)*f_tn5_yn5 + (19087/60480)*f_tn6_yn6)
    	t0 += h

    	arr_y.append(y0)

    if method == "nothing":
        print_in_file("Adams Bashforth de 7 Ordem", h)
        plt.title("Metodo de Adams Bashforth de 7 Ordem")
    
    elif method == "Euler":
        print_in_file("Adams Bashforth de 7 Ordem por Euler", h)
        plt.title("Metodo de Adams Bashforth de 7 Ordem por Euler")

    elif method == "Euler Inverso":
        print_in_file("Adams Bashforth de 7 Ordem por Euler Inverso", h)
        plt.title("Metodo de Adams Bashforth de 7 Ordem por Euler Inverso")

    elif method == "Euler Aprimorado":
        print_in_file("Adams Bashforth de 7 Ordem por Euler Aprimorado", h)
        plt.title("Metodo de Adams Bashforth de 7 Ordem por Euler Aprimorado")

    elif method == "Runge-Kutta":
        print_in_file("Adams Bashforth de 7 Ordem por Runge-Kutta", h)
        plt.title("Metodo de Adams Bashforth de 7 Ordem por Runge-Kutta")

    show_graphic()

def ab_6(t0, h, n, fty, method):
	expr = parser.expr(fty).compile()

	y0 = arr_y[len(arr_y) - 1]

	for i in range(len(arr_y), n + 1):
		arr_t.append(t0)

        # f(tn, yn)
		f_tn_yn = do_eval(arr_t[len(arr_y) - 1], y0, expr)

        # f(tn-1, yn-1)
		f_tn1_yn1 = do_eval(arr_t[len(arr_y) - 2], arr_y[len(arr_y) - 2], expr)

        # f(tn-2, yn-2)
		f_tn2_yn2 = do_eval(arr_t[len(arr_y) - 3], arr_y[len(arr_y) - 3], expr)

        # f(tn-3, yn-3)
		f_tn3_yn3 = do_eval(arr_t[len(arr_y) - 4], arr_y[len(arr_y) - 4], expr)

        # f(tn-4, yn-4)
		f_tn4_yn4 = do_eval(arr_t[len(arr_y) - 5], arr_y[len(arr_y) - 5], expr)

        # f(tn-5, yn-5)
		f_tn5_yn5 = do_eval(arr_t[len(arr_y) - 6], arr_y[len(arr_y) - 6], expr)

        # yn+1 = yn + h*((4277/1440)*f(tn, yn) - (2641/480)*f(tn-1, yn-1) + (4991/720)*f(tn-2, yn-2) - (3649/720)*f(tn-3, yn-3) + (959/480)*f(tn-4, yn-4) - (95/288)*f(tn-5, yn-5))
		y0 += h*((4277/1440)*f_tn_yn - (2641/480)*f_tn1_yn1 + (4991/720)*f_tn2_yn2 - (3649/720)*f_tn3_yn3 + (959/480)*f_tn4_yn4 - (95/288)*f_tn5_yn5)
		t0 += h

		arr_y.append(y0)

	if method == "nothing":
		print_in_file("Adams Bashforth de 6 Ordem", h)
		plt.title("Metodo de Adams Bashforth de 6 Ordem")
    
	elif method == "Euler":
		print_in_file("Adams Bashforth de 6 Ordem por Euler", h)
		plt.title("Metodo de Adams Bashforth de 6 Ordem por Euler")

	elif method == "Euler Inverso":
		print_in_file("Adams Bashforth de 6 Ordem por Euler Inverso", h)
		plt.title("Metodo de Adams Bashforth de 6 Ordem por Euler Inverso")

	elif method == "Euler Aprimorado":
		print_in_file("Adams Bashforth de 6 Ordem por Euler Aprimorado", h)
		plt.title("Metodo de Adams Bashforth de 6 Ordem por Euler Aprimorado")

	elif method == "Runge-Kutta":
		print_in_file("Adams Bashforth de 6 Ordem por Runge-Kutta", h)
		plt.title("Metodo de Adams Bashforth de 6 Ordem por Runge-Kutta")

	show_graphic()

def ab_5(t0, h, n, fty, method):
	expr = parser.expr(fty).compile()

	y0 = arr_y[len(arr_y) - 1]

	for i in range(len(arr_y), n + 1):
		arr_t.append(t0)

		# f(tn, yn)
		f_tn_yn = do_eval(arr_t[len(arr_y) - 1], y0, expr)

        # f(tn-1, yn-1)
		f_tn1_yn1 = do_eval(arr_t[len(arr_y) - 2], arr_y[len(arr_y) - 2], expr)

        # f(tn-2, yn-2)
		f_tn2_yn2 = do_eval(arr_t[len(arr_y) - 3], arr_y[len(arr_y) - 3], expr)

        # f(tn-3, yn-3)
		f_tn3_yn3 = do_eval(arr_t[len(arr_y) - 4], arr_y[len(arr_y) - 4], expr)

        # f(tn-4, yn-4)
		f_tn4_yn4 = do_eval(arr_t[len(arr_y) - 5], arr_y[len(arr_y) - 5], expr)

        # yn+1 = yn + h*((1901/720)*f(tn, yn) - (1387/360)*f(tn-1, yn-1) + (109/30)*f(tn-2, yn-2) - (637/360)*f(tn-3, yn-3) + (251/720)*f(tn-4, yn-4))
		y0 += h*((1901/720)*f_tn_yn - (1387/360)*f_tn1_yn1 + (109/30)*f_tn2_yn2 - (637/360)*f_tn3_yn3 + (251/720)*f_tn4_yn4)
		t0 += h

		arr_y.append(y0)

	if method == "nothing":
		print_in_file("Adams Bashforth de 5 Ordem", h)
		plt.title("Metodo de Adams Bashforth de 5 Ordem")
    
	elif method == "Euler":
		print_in_file("Adams Bashforth de 5 Ordem por Euler", h)
		plt.title("Metodo de Adams Bashforth de 5 Ordem por Euler")

	elif method == "Euler Inverso":
		print_in_file("Adams Bashforth de 5 Ordem por Euler Inverso", h)
		plt.title("Metodo de Adams Bashforth de 5 Ordem por Euler Inverso")

	elif method == "Euler Aprimorado":
		print_in_file("Adams Bashforth de 5 Ordem por Euler Aprimorado", h)
		plt.title("Metodo de Adams Bashforth de 5 Ordem por Euler Aprimorado")

	elif method == "Runge-Kutta":
		print_in_file("Adams Bashforth de 5 Ordem por Runge-Kutta", h)
		plt.title("Metodo de Adams Bashforth de 5 Ordem por Runge-Kutta")

	show_graphic()

def ab_4(t0, h, n, fty, method):
	expr = parser.expr(fty).compile()

	y0 = arr_y[len(arr_y) - 1]

	for i in range(len(arr_y), n + 1):
		arr_t.append(t0)

		# f(tn, yn)
		f_tn_yn = do_eval(arr_t[len(arr_y) - 1], y0, expr)

		# f(tn-1, yn-1)
		f_tn1_yn1 = do_eval(arr_t[len(arr_y) - 2], arr_y[len(arr_y) - 2], expr)

		# f(tn-2, yn-2)
		f_tn2_yn2 = do_eval(arr_t[len(arr_y) - 3], arr_y[len(arr_y) - 3], expr)

		# f(tn-3, yn-3)
		f_tn3_yn3 = do_eval(arr_t[len(arr_y) - 4], arr_y[len(arr_y) - 4], expr)

		# yn+1 = yn + h*((55/24)*f(tn, yn) - (59/24)*f(tn-1, yn-1) + (37/24)*f(tn-2, yn-2) - (3/8)*f(tn-3, yn-3))
		y0 += h*((55/24)*f_tn_yn - (59/24)*f_tn1_yn1 + (37/24)*f_tn2_yn2 - (3/8)*f(tn-3, yn-3))
		t0 += h

		arr_y.append(y0)

	if method == "nothing":
		print_in_file("Adams Bashforth de 4 Ordem", h)
		plt.title("Metodo de Adams Bashforth de 4 Ordem")
    
	elif method == "Euler":
		print_in_file("Adams Bashforth de 4 Ordem por Euler", h)
		plt.title("Metodo de Adams Bashforth de 4 Ordem por Euler")

	elif method == "Euler Inverso":
		print_in_file("Adams Bashforth de 4 Ordem por Euler Inverso", h)
		plt.title("Metodo de Adams Bashforth de 4 Ordem por Euler Inverso")

	elif method == "Euler Aprimorado":
		print_in_file("Adams Bashforth de 4 Ordem por Euler Aprimorado", h)
		plt.title("Metodo de Adams Bashforth de 4 Ordem por Euler Aprimorado")

	elif method == "Runge-Kutta":
		print_in_file("Adams Bashforth de 4 Ordem por Runge-Kutta", h)
		plt.title("Metodo de Adams Bashforth de 4 Ordem por Runge-Kutta")

	show_graphic()

def ab_3(t0, h, n, fty, method):
	expr = parser.expr(fty).compile()

	y0 = arr_y[len(arr_y) - 1]

	for i in range(len(arr_y), n + 1):
		arr_t.append(t0)

		# f(tn, yn)
		f_tn_yn = do_eval(arr_t[len(arr_y) - 1], y0, expr)

		# f(tn-1, yn-1)
		f_tn1_yn1 = do_eval(arr_t[len(arr_y) - 2], arr_y[len(arr_y) - 2], expr)

		# f(tn-2, yn-2)
		f_tn2_yn2 = do_eval(arr_t[len(arr_y) - 3], arr_y[len(arr_y) - 3], expr)

		# yn+1 = yn + h*((23/12)*f(tn, yn) - (4/3)*f(tn-1, yn-1) + (5/12)*f(tn-2, yn-2))
		y0 += h*((23/12)*f_tn_yn - (4/3)*f_tn1_yn1 + (5/12)*f_tn2_yn2)
		t0 += h

		arr_y.append(y0)

	if method == "nothing":
		print_in_file("Adams Bashforth de 3 Ordem", h)
		plt.title("Metodo de Adams Bashforth de 3 Ordem")
    
	elif method == "Euler":
		print_in_file("Adams Bashforth de 3 Ordem por Euler", h)
		plt.title("Metodo de Adams Bashforth de 3 Ordem por Euler")

	elif method == "Euler Inverso":
		print_in_file("Adams Bashforth de 3 Ordem por Euler Inverso", h)
		plt.title("Metodo de Adams Bashforth de 3 Ordem por Euler Inverso")

	elif method == "Euler Aprimorado":
		print_in_file("Adams Bashforth de 3 Ordem por Euler Aprimorado", h)
		plt.title("Metodo de Adams Bashforth de 3 Ordem por Euler Aprimorado")

	elif method == "Runge-Kutta":
		print_in_file("Adams Bashforth de 3 Ordem por Runge-Kutta", h)
		plt.title("Metodo de Adams Bashforth de 3 Ordem por Runge-Kutta")

	show_graphic()

def ab_2(t0, h, n, fty, method):
	expr = parser.expr(fty).compile()

	y0 = arr_y[len(arr_y) - 1]

	for i in range(len(arr_y), n + 1):
		arr_t.append(t0)

		# f(tn, yn)
		f_tn_yn = do_eval(arr_t[len(arr_y) - 1], y0, expr)

		# f(tn-1, yn-1)
		f_tn1_yn1 = do_eval(arr_t[len(arr_y) - 2], arr_y[len(arr_y) - 2], expr)

		# yn+1 = yn + h*((3/2)*f(tn, yn) - (1/2)*f(tn-1, yn-1))
		y0 += h*((3/2)*f_tn_yn - (1/2)*f_tn1_yn1)
		t0 += h

		arr_y.append(y0)

	if method == "nothing":
		print_in_file("Adams Bashforth de 2 Ordem", h)
		plt.title("Metodo de Adams Bashforth de 2 Ordem")
    
	elif method == "Euler":
		print_in_file("Adams Bashforth de 2 Ordem por Euler", h)
		plt.title("Metodo de Adams Bashforth de 2 Ordem por Euler")

	elif method == "Euler Inverso":
		print_in_file("Adams Bashforth de 2 Ordem por Euler Inverso", h)
		plt.title("Metodo de Adams Bashforth de 2 Ordem por Euler Inverso")

	elif method == "Euler Aprimorado":
		print_in_file("Adams Bashforth de 2 Ordem por Euler Aprimorado", h)
		plt.title("Metodo de Adams Bashforth de 2 Ordem por Euler Aprimorado")

	elif method == "Runge-Kutta":
		print_in_file("Adams Bashforth de 2 Ordem por Runge-Kutta", h)
		plt.title("Metodo de Adams Bashforth de 2 Ordem por Runge-Kutta")

	show_graphic()


################################################################################### MAIN ###################################################################################
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
      			   
			euler_method(float(y0), float(t0), float(h), int(n), fty[0], 1)
        
		elif strin[0] == "euler_inverso":
			y0 = strin[1]
			t0 = strin[2]
			h = strin[3]
			n = strin[4]
			fty = strin[5].split("\n ")
         
			backward_euler_method(float(y0), float(t0), float(h), int(n), fty[0], 1)
            
		elif strin[0] == "euler_aprimorado":
			y0 = strin[1]
			t0 = strin[2]
			h = strin[3]
			n = strin[4]
			fty = strin[5].split("\n ")
         
			best_euler_method(float(y0), float(t0), float(h), int(n), fty[0], 1)

		elif strin[0] == "runge_kutta":
			y0 = strin[1]
			t0 = strin[2]
			h = strin[3]
			n = strin[4]
			fty = strin[5].split("\n ")
         
			runge_kutta_method(float(y0), float(t0), float(h), int(n), fty[0], 1)

		elif strin[0] == "adam_bashforth":
			order_str = strin[len(strin) - 1].split("\n ")
			order = int(order_str[0])

			arr_y = [float(a) for a in strin[1:order + 1]]
       			     
			t0 = strin[order + 1]
			h = strin[order + 2]
			n = strin[order + 3]
			fty = strin[order + 4]

			t0 = float(t0)
			h = float(h)

			for i in range(0, len(arr_y)):
				arr_t.append(t0)
				t0 += h

			adams_bashforth_method(t0, h, int(n), fty, int(order), "nothing")

		elif strin[0] == "adam_bashforth_by_euler":
			y0 = strin[1]
			t0 = strin[2]
			h = strin[3]
			n = strin[4]
			fty = strin[5]
			order = strin[6].split("\n ")

			t0 = euler_method(float(y0), float(t0), float(h), int(order[0]) - 1, fty, 0)
			adams_bashforth_method(float(t0), float(h), int(n), fty, int(order[0]), "Euler")

		elif strin[0] == "adam_bashforth_by_euler_inverso":
			y0 = strin[1]
			t0 = strin[2]
			h = strin[3]
			n = strin[4]
			fty = strin[5]
			order = strin[6].split("\n ")

			t0 = backward_euler_method(float(y0), float(t0), float(h), int(order[0]) - 1, fty, 0)
			adams_bashforth_method(float(t0), float(h), int(n), fty, int(order[0]), "Euler Inverso")

		elif strin[0] == "adam_bashforth_by_euler_aprimorado":
			y0 = strin[1]
			t0 = strin[2]
			h = strin[3]
			n = strin[4]
			fty = strin[5]
			order = strin[6].split("\n ")

			t0 = best_euler_method(float(y0), float(t0), float(h), int(order[0]) - 1, fty, 0)
			adams_bashforth_method(float(t0), float(h), int(n), fty, int(order[0]), "Euler Aprimorado")

		elif strin[0] == "adam_bashforth_by_runge_kutta":
			y0 = strin[1]
			t0 = strin[2]
			h = strin[3]
			n = strin[4]
			fty = strin[5]
			order = strin[6].split("\n ")

			t0 = runge_kutta_method(float(y0), float(t0), float(h), int(order[0]) - 1, fty, 0)
			adams_bashforth_method(float(t0), float(h), int(n), fty, int(order[0]), "Runge-Kutta")

################################################################################# FIM DA MAIN #################################################################################






