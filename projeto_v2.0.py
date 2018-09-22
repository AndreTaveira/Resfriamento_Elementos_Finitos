'''
    Exercício Programa 3 - Método de Romberg para Integração de Funções

    Integrantes                |  NUSP      | Turma
    André Angelo Taveira       |  8586749   |  03
    Fernando Tramarim Trovões  |  7991642   |  03
'''

import math
n_prod=0
ksi=0
kal=0
L=0
d=0
t_amb = 0
Qzmais = 0
Qzmenos = 0
sigma = 0
theta = 0

def decomposicao_LU_tridiagonal(a, b, c):
    '''(vetor, vetor, vetor) -> (vetor, vetor)
    Recebe os vetores que compoe a matriz tridiagonal Anxn e devolve o vetor dos multiplcadores da matriz L e a diagonal principal da matriz U'''
    u = []
    l = []
    u.append(b[0])
    l.append(0)
    for i in range(1, len(a), 1):
        l.append(a[i]/u[i-1])
        u.append(b[i] - l[i]*c[i-1])
    return (l,u)

def resolucao_tridiagonal(d,l,u,c):
    '''(vetor, vetor, vetor, vetor) -> (vetor)
    Resolve o sistema linear Ax = d, onde A é uma matriz tridiagonal nxn
    Recebe o vetor d e a decomposicao LU de A e devolve o vetor x, solucao do sistema'''
    y = []
    x = []
    '''Ly = d'''
    y.append(d[0])
    for i in range(1, len(d), 1):
        y.append(d[i] - l[i]*y[i-1])
    '''Ux = y'''
    x.append(y[len(d)-1]/u[len(d)-1])
    for i in range(1, len(d),1):
        x.append((y[len(d)-1-i]-c[len(d)-1-i]*x[i-1])/u[len(d)-1-i])
    x.reverse()
    return x

def romb(a,b,n,e,ITMAX,funcao,i_prod):
    '''
    a => /float/ limite inferior de integração
    b => /float/ limite superior de integração
    n => /inteiro/ valor especificado que determina tamanho do triangulo do método de Romberg
    e => /float/ epsilon especificado para delimitação de erro
    ITMAX => /inteiro/ número máximo de iterações
    funcao => /função/ este argumento é uma função que calcula uma f(x) especificada
    '''

    '''
    Neste momento criamos a matrizT, que é um array de arrays que armazena um triangulo de Tik's.
    '''
    ite = 0
    Th0 = trapz(a,b,0,0,funcao,i_prod)
    matrizT = criaMatriz_romb(a,b,n,Th0,funcao,i_prod)

    '''
    O triangulo em matrizT será deslocado para baixo, sendo cada valor da coluna deslocado para uma posição anterior e
    acrescentado o novo valor na ultima posição (o primeiro valor de cada coluna é descartado). Uma iterada, portanto,
    consiste em deslocar esse triangulo para baixo.
    '''
    while( (ite<ITMAX) and (abs(matrizT[n][0]-matrizT[n-1][1]) >= (e*abs(matrizT[n][0]))) ):
        for p in range (0,n,1):
            shift_array(matrizT[p])
        matrizT[0][len(matrizT[0])-1] = trapz(a,b,matrizT[0][len(matrizT[0])-2],(n+1+ite),funcao,i_prod)
        
        for i in range (1,n+1,1):
            matrizT[i][len(matrizT[i])-1] = (matrizT[i-1][len(matrizT[i-1])-1]+((matrizT[i-1][len(matrizT[i-1])-1]-matrizT[i-1][len(matrizT[i-1])-2])/((4**i)-1)))
        ite += 1

    return matrizT[n][0]
    
'Esta função desloca todos valores de um vetor, menos o primeiro, para a posição anterior. Note que é usada no método da iterada em romb'
def shift_array(arr=[]):
    n=len(arr)
    for i in range (0,n-1,1):
        arr[i] = arr[i+1]

def criaMatriz_romb(a,b,n,Th0,funcao,i_prod):
    '''
    Nesta função, iremos criar um conjunto de arrays representantes das colunas do triangulo
    de Tik's com tamanho n (ou seja, ultimo T será Tnn). Futuramente a estrutura desse
    triangulo será utilizado para processar as futuras iterações.
    '''

    '''
    Declaramos um array que irá armazenar arrays como elementos. Cada elemento de
    matrizT será um array que representa uma coluna de Tik's.
    '''
    matrizT=[[]]
    '''
    Preencheremos então a primeira coluna, que deve ser gerada pela funcao trapz.
    '''
    matrizT[0].append(Th0)
    for i in range (1,n+1,1):
        matrizT[0].append(trapz(a,b,matrizT[0][i-1],i,funcao,i_prod))


    'Preencheremos aqui as outras colunas, que derivam da primeira'
    for k in range (1,n+1,1):
        matrizT.append([])
        for i in range (k,n+1,1):
            matrizT[k].append(matrizT[k-1][i-k+1]+((matrizT[k-1][i-k+1]-matrizT[k-1][i-k])/((4**k)-1)))
    return matrizT


'''
    Abaixo definimos as funções que serão integradas. Note que essas funções
    servirão de argumento no chamamento da função romb.
'''

def prod_int1(x,i_prod):
    Q = (12*x*(1-x))-2
    if (x <= i_prod/(n_prod+1) and (x >= (i_prod-1)/(n_prod+1))):
        Phi = (x-((i_prod-1)/(n_prod+1))) * (n_prod+1)
    else:
        if (x > i_prod/(n_prod+1) and (x <= (i_prod+1)/(n_prod+1))):
            Phi = (((i_prod+1)/(n_prod+1))-x) * (n_prod+1)
        else:
            Phi = 0
    v = Q*Phi

    return v

def prod_int2(x,i_prod):
    Q = Qzmais - Qzmenos
    if (x <= (L*i_prod)/(n_prod+1) and (x >= (L*(i_prod-1))/(n_prod+1))):
        Phi = (((x/L)-((i_prod-1)/(n_prod+1))) * (n_prod+1))
    else:
        if (x > (L*i_prod)/(n_prod+1) and (x <= (L*(i_prod+1))/(n_prod+1))):
            Phi = ((((i_prod+1)/(n_prod+1))-(x/L)) * (n_prod+1))
        else:
            Phi = 0
    v = Q*Phi
    return v

def prod_int3(x,i_prod):
    Qmais = Qzmais*(math.exp(((-1*((x-(L/2))**2))/(sigma**2))))
    Qmenos = Qzmenos
    Q = Qmais - Qmenos
    if (x <= (L*i_prod)/(n_prod+1) and (x >= (L*(i_prod-1))/(n_prod+1))):
        Phi = (((x/L)-((i_prod-1)/(n_prod+1))) * (n_prod+1))
    else:
        if (x > (L*i_prod)/(n_prod+1) and (x <= (L*(i_prod+1))/(n_prod+1))):
            Phi = ((((i_prod+1)/(n_prod+1))-(x/L)) * (n_prod+1))
        else:
            Phi = 0
    v = Q*Phi
    return v

def prod_int4(x,i_prod):
    Qmais = Qzmais*(math.exp(((-1*((x-(L/2))**2))/(sigma**2))))
    Qmenos = Qzmenos*(math.exp(((-1*((x)**2))/(theta**2))) + math.exp(((-1*((x-L)**2))/(theta**2))))
    Q = Qmais - Qmenos
    if (x <= (L*i_prod)/(n_prod+1) and (x >= (L*(i_prod-1))/(n_prod+1))):
        Phi = (((x/L)-((i_prod-1)/(n_prod+1))) * (n_prod+1))
    else:
        if (x > (L*i_prod)/(n_prod+1) and (x <= (L*(i_prod+1))/(n_prod+1))):
            Phi = ((((i_prod+1)/(n_prod+1))-(x/L)) * (n_prod+1))
        else:
            Phi = 0
    v = Q*Phi
    return v


'Esta função serve apenas como apoio para o calculo de uma somatória em trapz'
def soma(a,h,i,funcao,i_prod):
    s=0
    for j in range (1,(2**(i-1)+1),1):
        s += funcao(a+(2*j-1)*h,i_prod)
    return s

'Aqui há uma implementação simples do metodo dos trapézio, conforme implementação sugerida na apostila do exercício'
def trapz(a,b,t,i,funcao,i_prod):
    if (i>0):
        h_i= (b-a)/(2**i)
        t = 1/2*t+ h_i*soma(a,h_i,i,funcao,i_prod)
    if (i==0):
        t=(b-a)*(funcao(a,i_prod)+funcao(b,i_prod))/2

    return t

def preenche_vetor(v,a,n):
    ''' dado um vetor criado, sem elementos, preenche este vetor com um valor especificado "a" com "n" elementos'''
    ''' vetor, float, inteiro -> null '''
    for i in range (0,n,1):
        v.append(a)
    return

def produto_escalar(u,v):
    '''u e v devem ter o mesmo tamanho'''
    prod_escalar = 0
    for k in range (0,len(u),1):
        prod_escalar += u[k]*v[k]
    return prod_escalar

def tridiagonal_homogenea():
    a=[]
    b=[]
    c=[]
    a.append(0)
    preenche_vetor(a,((-ksi*(n_prod+1))/L),n_prod-1)
    preenche_vetor(b,((2*ksi*(n_prod+1))/L),n_prod)
    preenche_vetor(c,((-ksi*(n_prod+1))/L),n_prod-1)
    c.append(0)
    return a,b,c

def diagonais_nHomogeneas(): 
    b=[]
    for i in range (1,n_prod+1,1):
        if (((i-1)*(L/(n_prod+1))) < (L/2-d)):
            b.append((2*kal*(n_prod+1))/L)
        if (((i-1)*(L/(n_prod+1))) >= (L/2-d)) and (((i+1)*(L/(n_prod+1))) <= (L/2+d)) :
            b.append((2*ksi*(n_prod+1))/L)
        if (((i+1)*(L/(n_prod+1))) > (L/2+d)):
            b.append((2*kal*(n_prod+1))/L)
    a=[]
    a.append(0)
    for i in range (1,n_prod,1):
        if (((i-1)*(L/(n_prod+1))) < (L/2-d)):
            a.append((-kal*(n_prod+1))/L)
        if (((i-1)*(L/(n_prod+1))) >= (L/2-d)) and (((i+1)*(L/(n_prod+1))) <= (L/2+d)):
            a.append((-ksi*(n_prod+1))/L)
        if (((i+1)*(L/(n_prod+1))) > (L/2+d)):
            a.append((-kal*(n_prod+1))/L)
    c=[]
    for i in range (1,n_prod,1):
        if (((i-1)*(L/(n_prod+1))) < (L/2-d)):
            c.append((-kal*(n_prod+1))/L)
        if (((i-1)*(L/(n_prod+1))) >= (L/2-d)) and (((i+1)*(L/(n_prod+1))) <= (L/2+d)):
            c.append((-ksi*(n_prod+1))/L)
        if (((i+1)*(L/(n_prod+1))) > (L/2+d)):
            c.append((-kal*(n_prod+1))/L)
    c.append(0)

    return a,b,c

def d_tri_homogenea(n,e,ITMAX,funcao):
    d_th=[]
    for i_prod in range (1,n_prod+1,1):
        d_th.append(romb((L*(i_prod-1))/(n_prod+1),(L*(i_prod+1))/(n_prod+1),n,e,ITMAX,funcao,i_prod))
    return d_th

def gera_phi(x):
    phi_vetor=[]
    for i in range (1,n_prod+1,1):
        if ((x >= (L*(i-1))/(n_prod+1)) and (x <= (L*i)/(n_prod+1))):
            Phi = (x/L - ((i-1)/(n_prod+1)))*(n_prod+1)
            phi_vetor.append(Phi)
        else:
            if (x > (L*i)/(n_prod+1) and (x <= (L*(i+1))/(n_prod+1))):
                Phi = (((i+1)/(n_prod+1))-x/L) * (n_prod+1)
            else:
                Phi = 0
            phi_vetor.append(Phi)
    return phi_vetor

def resolve_matriz_eq(funcao, tridiagonal):
    
    n=4
    a,b,c = tridiagonal()
    Lt,Ut = decomposicao_LU_tridiagonal(a, b, c)
    d_th =  d_tri_homogenea(n,0.000001,15,funcao)
    alpha = resolucao_tridiagonal(d_th,Lt,Ut,c)

    return alpha

def resolve_eqDif(x, alpha):
    
    phi_vetor= gera_phi(x)
    un_barra = produto_escalar(alpha,phi_vetor)

    return un_barra + t_amb

def main():
    global L
    global n_prod
    global ksi
    global kal
    global d
    global t_amb
    global Qzmais
    global Qzmenos
    global sigma
    global theta
    encerrar = False
    while(not encerrar):
    
        print("")
        print("Modos: ")
        print("1 - Teste de convergência. ")
        print("2 - Teste nas condições do relatório. ")
        print("3 - Teste para parâmetros determinados pelo usuário. ")
        print("")
        modo = int(input("Digite 1, 2 ou 3 conforme seleção desejada: "))
        print("")
        
        if (modo == 1):

            L, ksi = 1, 1
            'Teste de convergência' 
            N_testes = [7,15,31,63,127,255]
            for j in N_testes:
                n_prod = j
                norm_dif = 0
                alpha = resolve_matriz_eq(prod_int1,tridiagonal_homogenea)
                h_quadrado = (1/(j+1))**2
                print("N é: %d"%j)
                print("h^2: %.8f "%h_quadrado)
                for k in range (0,10*j+1,1):
                    y = k/(10*j)
                    phi_vetor= gera_phi(y)
                    un_barrax = produto_escalar(alpha,phi_vetor)
                    un = (y**2)*((y-1)**2)
                    dif = abs(un_barrax - un)
                    if (norm_dif < dif):
                        norm_dif = dif
                print("Valor do erro ||u-un||max para n = %d é: %.8f"%(j,norm_dif))
                print("h^2/||u-un||max é: %.8f"%(h_quadrado/norm_dif))
                print("")
                
        elif (modo == 2):
            sair = False
            while(not sair):
                ksi = float(input("De um valor para k [W/(m*K)] do silício: "))
                kal = float(input("De um valor para k [W/(m*K)] do aluminio: "))
                print("ATENÇÃO: N deve ser alto para o teste com k não homogêneo ter boa aproximação")
                n_prod = int(input("De um valor inteiro para N, tamanho da matriz de solução do sistema: "))
                L = float(input("De o tamanho do comprimento do chip L em metros: "))
                t_amb = float(input("De a temperatura ambiente em kelvins: "))
                Qzmais = float(input("Dê o valor de Q0+: "))
                Qzmenos = float(input("Dê o valor de Q0-: "))
                sigma = float(input("Dê o valor de Sigma: "))
                theta = float(input("Dê o valor de Theta: "))
                d = float(input("De o tamanho de d (raio em torno do centro do chip que delimita a parte de silício) em metros: "))
                print("")


                # Teste para forçante de calor constante
                print("Chip de material homogêneo com Q+ e Q- constante")
                alpha = resolve_matriz_eq(prod_int2, tridiagonal_homogenea)
                x = 0
                for t in range (0,21,1):
                    result = resolve_eqDif(x, alpha)
                    print("%.2f"%result)
                    x += L/20
                print("")

                #Teste para forçante de calor variavel 1
                print("Chip de material homogêneo com Q+ gaussiana e Q- constante")
                alpha = resolve_matriz_eq(prod_int3, tridiagonal_homogenea)
                x = 0
                for t in range (0,21,1):
                    result = resolve_eqDif(x, alpha)
                    print("%.2f"%result)
                    x += L/20
                print("")


                #Teste para forçante de calor variavel 2
                print("Chip de material homogêneo com Q+ e Q- gaussiana")
                alpha = resolve_matriz_eq(prod_int4, tridiagonal_homogenea)
                x = 0
                for t in range (0,21,1):
                    result = resolve_eqDif(x, alpha)
                    print("%.2f"%result)
                    x += L/20
                print("")


                print("Teste para chips constituidos por materiais não homogêneos")
                #Teste para forçante de calor variavel 3
                alpha = resolve_matriz_eq(prod_int4, diagonais_nHomogeneas)
                
                x = 0
                for t in range (0,21,1):
                    result = resolve_eqDif(x, alpha)
                    print("%.2f"%result)
                    x += L/20
                print("")


                saida = int(input("Para novo teste digite 1, para sair digite 0: "))
                if (saida == 0):
                    sair = True
                               
        elif (modo == 3):
            sair = False
            while(not sair):
                ksi = float(input("De um valor para k [W/(m*K)] do silício: "))
                kal = float(input("De um valor para k [W/(m*K)] do aluminio: "))
                print("ATENÇÃO: N deve ser alto para o teste com k não homogêneo")
                n_prod = int(input("De um valor inteiro para N, tamanho da matriz de solução do sistema: "))
                L = float(input("De o tamanho do comprimento do chip L em metros: "))
                t_amb = float(input("De a temperatura ambiente em kelvins: "))
                x = float(input("De a posicao a ser testada em metros: "))
                Qzmais = float(input("Dê o valor de Q0+"))
                Qzmenos = float(input("Dê o valor de Q0-"))
                sigma = float(input("Dê o valor de Sigma"))
                theta = float(input("Dê o valor de Theta"))
                d = float(input("De o tamanho de d (raio em torno do centro do chip que delimita a parte de silício) em metros: "))
                
                # Teste para forçante de calor constante
                print("Chip de material homogêneo com Q+ e Q- constante")
                alpha = resolve_matriz_eq(prod_int2, tridiagonal_homogenea)
                result = resolve_eqDif(x, alpha)
                print("T(%f) = %f"%(x,result))
            
                #Teste para forçante de calor variavel 1
                print("Chip de material homogêneo com Q+ gaussiana e Q- constante")
                alpha = resolve_matriz_eq(prod_int4, tridiagonal_homogenea)
                result = resolve_eqDif(x, alpha)
                print("T(%f) = %f"%(x,result))

                #Teste para forçante de calor variavel 2
                print("Chip de material homogêneo com Q+ e Q- gaussiana")
                alpha = resolve_matriz_eq(prod_int4, tridiagonal_homogenea)
                result = resolve_eqDif(x, alpha)
                print("T(%f) = %f"%(x,result))

                print("Faremos agora um teste para chips constituidos por materiais não homogêneos")
                #Teste para forçante de calor variavel 3
                alpha = resolve_matriz_eq(prod_int2, diagonais_nHomogeneas)
                result = resolve_eqDif(x, alpha)
                print("T(%f) = %f"%(x,result))

                saida = int(input("Para novo teste digite 1, para sair digite 0: "))
                if (saida == 0):
                    sair = True
        else:
            print("Digite o modo corretamente")
            
        encerra = int(input("Para continuar digite 1, para encerrar programa digite 0: "))
        if (encerra == 0):
            encerrar = True
        
    return 0
    
main()
