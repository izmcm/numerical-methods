# metodos-numericos
Implementação dos métodos numéricos dados na disciplina Métodos Numéricos e Computacionais do Centro de Informática - UFPE

Sistemas Operacionais testados: macOS 10 e Ubuntu 18

### Pré-requisitos e instalações
 * [Python3](https://www.python.org/download/releases/3.0/)

No macOS usando [Homebrew](https://brew.sh/)
```
$ brew install python3
```

No Ubuntu
```
$ sudo apt-get install python3
$ sudo apt-get install python3-pip
```

 * [MpMath](http://mpmath.org/) - para ler expressões
```
$ sudo pip3 install mpmath
```

 * [Matplotlib](https://matplotlib.org/) - para plotar gráficos
```
$ sudo apt-get install python3-matplotlib
```
### Execução

1. Caminhe até o diretório onde foi baixado o programa
2. No arquivo "entrada.txt", escreva as entradas no formato abaixo:

- Para euler, euler inverso, euler aprimorado e runge-kutta:  
nome-do-método y0 t0 h n f(t,y)

- Para adams-bashforth, adams-multon e fórmula inversa:  
método lista-de-valores-iniciais-de-y' t0 h n f(t,y) ordem"

- Para adams bashforth, adams multon e fórmula inversa por outros métodos:  
método y0 t0 h n f(t,y) ordem"

**observações**  
' a quantidade de valores iniciais é igual ao valor da ordem que se deseja calcular em bashforth e da ordem-1 em multon e fórmula inversa  

" adams bashforth e adams multon calculam ordens de 1 a 7  
" fórmula inversa calcula ordens de 2 a 6  

* Exemplo de entrada:
```
euler 0 0 0.1 20 1-t+4*y
euler_inverso 0 0 0.1 20 1-t+4*y
euler_aprimorado 0 0 0.1 20 1-t+4*y
runge_kutta 0 0 0.1 20 1-t+4*y
adam_bashforth 0.0 0.1 0.23 0.402 0.6328 0 0.1 20 1-t+4*y 5
adam_bashforth_by_euler 0 0 0.1 20 1-t+4*y 6
adam_bashforth_by_euler_inverso 0 0 0.1 20 1-t+4*y 6
adam_bashforth_by_euler_aprimorado 0 0 0.1 20 1-t+4*y 6
adam_bashforth_by_runge_kutta 0 0 0.1 20 1-t+4*y 6
adam_multon 0.0 0.1 0.23 0.402 0.6328 0 0.1 20 1-t+4*y 6
adam_multon_by_euler 0 0 0.1 20 1-t+4*y 6
adam_multon_by_euler_inverso 0 0 0.1 20 1-t+4*y 6
adam_multon_by_euler_aprimorado 0 0 0.1 20 1-t+4*y 6
adam_multon_by_runge_kutta 0 0 0.1 20 1-t+4*y 6
formula_inversa 0.0 0.1 0.23 0.402 0.6328 0 0.1 20 1-t+4*y 6
formula_inversa_by_euler 0 0 0.1 20 1-t+4*y 6
formula_inversa_by_euler_inverso 0 0 0.1 20 1-t+4*y 6
formula_inversa_by_euler_aprimorado 0 0 0.1 20 1-t+4*y 6
formula_inversa_by_runge_kutta 0 0 0.1 20 1-t+4*y 6
```

3. Execute o código com
```
$ python3 solver.py
```

4. Os gráficos ficarão disponíveis conforme a execução, mas os valores de y0 de cada passo deverão ser consultados no arquivo "saida.txt"

* Exemplo de saída no arquivo "saida.txt"


