#Imports 

import numpy as np
import time
from concurrent.futures import ProcessPoolExecutor
import json

mu_0 = 4 * np.pi * 1e-7    

a_leg = 0.0281  # Profundidade  (parte que está dentro da janela)
b_leg = 0.0260  # Largura       (parte que está fora da janela)

### Largura e Altura da janela
w_w = 0.0238    # Largura   (x)
h_w = 0.0519    # Altura    (y)

### PRIMÁRIO
N_1 = 19        # Número de voltas
I_1 = 1         # Corrente
a_1 = 0.0038    # Largura
h_1 = 0.0489    # Altura

J_1 = (N_1 * I_1) / (a_1 * h_1)  # Densidade de corrente


### SECUNDÁRIO
N_2 = 18          # Número de voltas
I_2 = -N_1/N_2    # Corrente (deve ter sentido oposto ao primário)
a_2 = 0.0045      # Largura
h_2 = 0.0379      # Altura

J_2 = (N_2 * I_2) / (a_2 * h_2)  # Densidade de corrente


### Disposição dos enrolamentos na janela

h_1_minus = 0.0015 
h_1_plus = 0.0015 + h_1

h_2_minus = 0.0015 + 0.0055
h_2_plus = 0.0015 + 0.0055 + h_2

a_1_minus = 0
a_1_plus = a_1

a_2_minus = a_1 + 0.007 
a_2_plus = a_1 + 0.007 + a_2

def calc_termAm0Jm0(m):

    J_1_m0 = ((2 * J_1) / (m * h_w * np.pi)) * (np.sin((m * np.pi * a_1_plus) / w_w) - np.sin((m * np.pi * a_1_minus) / w_w)) * (h_1_plus - h_1_minus)
    A_1_m0 = (mu_0 * J_1_m0) / (((m * np.pi) / (w_w))**2)

    J_2_m0 = ((2 * J_2) / (m * h_w * np.pi)) * (np.sin((m * np.pi * a_2_plus) / w_w) - np.sin((m * np.pi * a_2_minus) / w_w)) * (h_2_plus - h_2_minus)
    A_2_m0 = (mu_0 * J_2_m0) / (((m * np.pi) / (w_w))**2)
    
    A_m0 = A_1_m0 + A_2_m0
    J_m0 = J_1_m0 + J_2_m0

    return A_m0 * J_m0 

def calc_termA0nJ0n(n):

    J_1_0n = ((2 * J_1) / (n * w_w * np.pi)) * (np.sin((n * np.pi * h_1_plus) / h_w) - np.sin((n * np.pi * h_1_minus) / h_w)) * (a_1_plus - a_1_minus)
    A_1_0n = (mu_0 * J_1_0n) / (((n * np.pi) / (h_w))**2)
     
    J_2_0n = ((2 * J_2) / (n * w_w * np.pi)) * (np.sin((n * np.pi * h_2_plus) / h_w) - np.sin((n * np.pi * h_2_minus) / h_w)) * (a_2_plus - a_2_minus)
    A_2_0n = (mu_0 * J_2_0n) / (((n * np.pi) / (h_w))**2)
    
    A_0n = A_1_0n + A_2_0n
    J_0n = J_1_0n + J_2_0n


    return A_0n * J_0n 

def calc_termAmnJmn(t):

    m, n = t

    J_1_mn = ((4 * J_1) / (m * n * (np.pi)**2)) * (np.sin((m * np.pi * a_1_plus) / w_w) - np.sin((m * np.pi * a_1_minus) / w_w)) * (np.sin((n * np.pi * h_1_plus) / h_w) - np.sin((n * np.pi * h_1_minus) / h_w))
    A_1_mn = (mu_0 * J_1_mn) / (((m * np.pi) / (w_w))**2 + ((n * np.pi) / (h_w))**2)
    
    J_2_mn = ((4 * J_2) / (m * n * (np.pi)**2)) * (np.sin((m * np.pi * a_2_plus) / w_w) - np.sin((m * np.pi * a_2_minus) / w_w)) * (np.sin((n * np.pi * h_2_plus) / h_w) - np.sin((n * np.pi * h_2_minus) / h_w))
    A_2_mn = (mu_0 * J_2_mn) / (((m * np.pi) / (w_w))**2 + ((n * np.pi) / (h_w))**2)        

    A_mn = A_1_mn + A_2_mn
    J_mn = J_1_mn + J_2_mn


    return A_mn * J_mn


def nonCircular():
    
    global w_w, h_w
    M = N = 100

    #####   INSIDE WINDOW   #####
    # Número de iterações dos somatórios
   
    # Variáveis dos somatórios
    sum_m0 = sum_0n = sum_mn = 0

    # Somatório A_m0*J_m0
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(calc_termAm0Jm0, range(1, M)))
    sum_m0 = sum(results)

    # Somatório A_0n*J_0n
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(calc_termA0nJ0n, range(1, N)))
    sum_0n = sum(results)

    # Somatório A_mn*J_mn

    with ProcessPoolExecutor() as executor:
            results = list(executor.map(calc_termAmnJmn, ((m, n) for m in range(1, M) for n in range(1, N + 1))))
    sum_mn = sum(results)

    # Corrente de referência 
    I_ref = I_1 

    ### Indutância de Dispersão p.u.l. (Dentro da Janela)
    Lk_pul_IW = ((w_w*h_w) / (2 * I_ref**2)) * (sum_m0 + sum_0n + 0.5*sum_mn)   # (4.16)

    print('Leakage p.u.l. dentro da janela:                         ', Lk_pul_IW, 'H / p.u.l.')


    ### Contribuição de Indutância de Dispersão Dentro da Janela:
    print('Contribuição de Indutância de Dispersão dentro da janela:', Lk_pul_IW*a_leg*2, 'H')

    #####   OUTSIDE WINDOW   #####

    ### "JANELA INFINITA"
    
    w_w = (0.0238+0.0519)*5  # Largura   (x)
    h_w = (0.0238+0.0519)*5  # Altura    (y)


    ### OS ENROLAMENTOS DEVEM ESTAR NO "CENTRO" DA "JANELA INFINITA" 

    h_1_minus = -h_1/2 + h_w/2
    h_1_plus = h_1/2 + h_w/2

    h_2_minus = -h_2/2 + h_w/2
    h_2_plus = h_2/2 + h_w/2


    # Número de iterações dos somatórios
    M = N = 130 

    # Variáveis dos somatórios
    sum_m0 = sum_0n = sum_mn = 0 

    # Somatório A_m0*J_m0
    for m in range(1, M):
        J_1_m0 = ((2 * J_1) / (m * h_w * np.pi)) * (np.sin((m * np.pi * a_1_plus) / w_w) - np.sin((m * np.pi * a_1_minus) / w_w)) * (h_1_plus - h_1_minus)
        A_1_m0 = (mu_0 * J_1_m0) / (((m * np.pi) / (w_w))**2)

        J_2_m0 = ((2 * J_2) / (m * h_w * np.pi)) * (np.sin((m * np.pi * a_2_plus) / w_w) - np.sin((m * np.pi * a_2_minus) / w_w)) * (h_2_plus - h_2_minus)
        A_2_m0 = (mu_0 * J_2_m0) / (((m * np.pi) / (w_w))**2)
        
        A_m0 = A_1_m0 + A_2_m0
        J_m0 = J_1_m0 + J_2_m0

        sum_m0 += A_m0 * J_m0

    # Somatório A_0n*J_0n
    for n in range(1, N):
        J_1_0n = ((2 * J_1) / (n * w_w * np.pi)) * (np.sin((n * np.pi * h_1_plus) / h_w) - np.sin((n * np.pi * h_1_minus) / h_w)) * (a_1_plus - a_1_minus)
        A_1_0n = (mu_0 * J_1_0n) / (((n * np.pi) / (h_w))**2)
        
        J_2_0n = ((2 * J_2) / (n * w_w * np.pi)) * (np.sin((n * np.pi * h_2_plus) / h_w) - np.sin((n * np.pi * h_2_minus) / h_w)) * (a_2_plus - a_2_minus)
        A_2_0n = (mu_0 * J_2_0n) / (((n * np.pi) / (h_w))**2)
        
        A_0n = A_1_0n + A_2_0n
        J_0n = J_1_0n + J_2_0n

        sum_0n += A_0n * J_0n

    # Somatório A_mn*J_mn
    for m in range(1, M):
        for n in range(1, N + 1):
            J_1_mn = ((4 * J_1) / (m * n * (np.pi)**2)) * (np.sin((m * np.pi * a_1_plus) / w_w) - np.sin((m * np.pi * a_1_minus) / w_w)) * (np.sin((n * np.pi * h_1_plus) / h_w) - np.sin((n * np.pi * h_1_minus) / h_w))
            A_1_mn = (mu_0 * J_1_mn) / (((m * np.pi) / (w_w))**2 + ((n * np.pi) / (h_w))**2)
            
            J_2_mn = ((4 * J_2) / (m * n * (np.pi)**2)) * (np.sin((m * np.pi * a_2_plus) / w_w) - np.sin((m * np.pi * a_2_minus) / w_w)) * (np.sin((n * np.pi * h_2_plus) / h_w) - np.sin((n * np.pi * h_2_minus) / h_w))
            A_2_mn = (mu_0 * J_2_mn) / (((m * np.pi) / (w_w))**2 + ((n * np.pi) / (h_w))**2)        

            A_mn = A_1_mn + A_2_mn
            J_mn = J_1_mn + J_2_mn

            sum_mn += A_mn * J_mn



    ### Indutância de Dispersão p.u.l. (Fora da Janela)
    Lk_pul_OW = ((w_w*h_w) / (2 * I_ref**2)) * (sum_m0 + sum_0n + 0.5*sum_mn)   # (4.16)

    print('')
    print('Indutância de Dispersão p.u.l. fora da janela:           ', Lk_pul_OW, 'H / p.u.l.')

    ### Contribuição de Indutância de Dispersão Fora da Janela:
    print('Contribuição de Indutância de Dispersão fora da janela:  ', Lk_pul_OW*b_leg*2, 'H')


    #####   OUTSIDE WINDOW   #####
    #####       CANTOS       #####


    # Número de iterações dos somatórios
    M = N = 150 

    # Variáveis dos somatórios
    sum_m0 = sum_0n = sum_mn = 0 


    J_1_dict_OW = {}
    J_2_dict_OW = {}


    for m in range(1, M):


        if (m, 0) not in J_1_dict_OW:
            J_1_m0 = ((2 * J_1) / (m * h_w * np.pi)) * (np.sin((m * np.pi * a_1_plus) / w_w) - np.sin((m * np.pi * a_1_minus) / w_w)) * (h_1_plus - h_1_minus)
            J_1_dict_OW[(m, 0)] = J_1_m0  
        else:
            J_1_m0 = J_1_dict_OW[(m, 0)] 


        if (m, 0) not in J_2_dict_OW:
            J_2_m0 = ((2 * J_2) / (m * h_w * np.pi)) * (np.sin((m * np.pi * a_2_plus) / w_w) - np.sin((m * np.pi * a_2_minus) / w_w)) * (h_2_plus - h_2_minus)
            J_2_dict_OW[(m, 0)] = J_2_m0  
        else:
            J_2_m0 = J_2_dict_OW[(m, 0)]  



        A_1_m0 = (mu_0 * J_1_m0) / (((m * np.pi) / (w_w))**2)
        A_2_m0 = (mu_0 * J_2_m0) / (((m * np.pi) / (w_w))**2)

        A_m0 = A_1_m0 + A_2_m0
        J_m0 = J_1_m0 + J_2_m0

    
        J_til_m0_temp = 0
        J_til_m0 = (0 + 1)*J_m0 


        for p in range(1, M):
            if (m + p) % 2 == 1: 
                if (p, 0) not in J_1_dict_OW:
                    J_1_p0 = ((2 * J_1) / (p * h_w * np.pi)) * (np.sin((p * np.pi * a_1_plus) / w_w) - np.sin((p * np.pi * a_1_minus) / w_w)) * (h_1_plus - h_1_minus)
                else:
                    J_1_p0 = J_1_dict_OW[(p, 0)]
            
                if (p, 0) not in J_2_dict_OW:
                    J_2_p0 = ((2 * J_2) / (p * h_w * np.pi)) * (np.sin((p * np.pi * a_2_plus) / w_w) - np.sin((p * np.pi * a_2_minus) / w_w)) * (h_2_plus - h_2_minus)
                else:
                    J_2_p0 = J_2_dict_OW[(p, 0)]

                J_p0 = J_1_p0 + J_2_p0

                J_til_m0_temp += (8/(np.pi**2)) * (1/(2*p**2) - (m**2 + p**2)/((m**2 - p**2)**2)) * J_p0



        sum_m0 += A_m0 * (J_til_m0 + J_til_m0_temp)


    for n in range(1, N):
        
        if (0, n) not in J_1_dict_OW:
            J_1_0n = ((2 * J_1) / (n * w_w * np.pi)) * (np.sin((n * np.pi * h_1_plus) / h_w) - np.sin((n * np.pi * h_1_minus) / h_w)) * (a_1_plus - a_1_minus)
            J_1_dict_OW[(0, n)] = J_1_0n  
        else:
            J_1_0n = J_1_dict_OW[(0, n)]  
        
        if (0, n) not in J_2_dict_OW:
            J_2_0n = ((2 * J_2) / (n * w_w * np.pi)) * (np.sin((n * np.pi * h_2_plus) / h_w) - np.sin((n * np.pi * h_2_minus) / h_w)) * (a_2_plus - a_2_minus)
            J_2_dict_OW[(0, n)] = J_2_0n  
        else:
            J_2_0n = J_2_dict_OW[(0, n)]  
        

        A_1_0n = (mu_0 * J_1_0n) / (((n * np.pi) / (h_w))**2)
        A_2_0n = (mu_0 * J_2_0n) / (((n * np.pi) / (h_w))**2)
        
        A_0n = A_1_0n + A_2_0n
        J_0n = J_1_0n + J_2_0n
        

        J_til_0n_temp = 0
        J_til_0n = (0 + 1)*J_0n 
        
        
        for p in range(1, M):
            if p % 2 == 1:  
                    
                    if (p, n) not in J_1_dict_OW:
                        J_1_pn = ((4 * J_1) / (p * n * (np.pi)**2)) * (np.sin((p * np.pi * a_1_plus) / w_w) - np.sin((p * np.pi * a_1_minus) / w_w)) * (np.sin((n * np.pi * h_1_plus) / h_w) - np.sin((n * np.pi * h_1_minus) / h_w))
                        J_1_dict_OW[(p, n)] = J_1_pn  
                    else:
                        J_1_pn = J_1_dict_OW[(p, n)]  

                    
                    if (p, n) not in J_2_dict_OW:
                        J_2_pn = ((4 * J_2) / (p * n * (np.pi)**2)) * (np.sin((p * np.pi * a_2_plus) / w_w) - np.sin((p * np.pi * a_2_minus) / w_w)) * (np.sin((n * np.pi * h_2_plus) / h_w) - np.sin((n * np.pi * h_2_minus) / h_w))
                        J_2_dict_OW[(p, n)] = J_2_pn  
                    else:
                        J_2_pn = J_2_dict_OW[(p, n)]  

                    J_pn = J_1_pn + J_2_pn

                    J_til_0n_temp += (8/(np.pi**2)) * (h_w**2/(p**2 * h_w**2 + n**2 * w_w**2) - (1/p**2)*(1 + (n**2 * w_w**2) / (p**2 * h_w**2 + n**2 * w_w**2))) * J_pn* 0.5

        
        sum_0n += A_0n * (J_til_0n + J_til_0n_temp)


    for m in range(1, M):
        for n in range(1, N):

            if (m, n) not in J_1_dict_OW:
                J_1_mn = ((4 * J_1) / (m * n * (np.pi)**2)) * (np.sin((m * np.pi * a_1_plus) / w_w) - np.sin((m * np.pi * a_1_minus) / w_w)) * (np.sin((n * np.pi * h_1_plus) / h_w) - np.sin((n * np.pi * h_1_minus) / h_w))
                J_1_dict_OW[(m, n)] = J_1_mn  
            else:
                J_1_mn = J_1_dict_OW[(m, n)]  

            if (m, n) not in J_2_dict_OW:
                J_2_mn = ((4 * J_2) / (m * n * (np.pi)**2)) * (np.sin((m * np.pi * a_2_plus) / w_w) - np.sin((m * np.pi * a_2_minus) / w_w)) * (np.sin((n * np.pi * h_2_plus) / h_w) - np.sin((n * np.pi * h_2_minus) / h_w))
                J_2_dict_OW[(m, n)] = J_2_mn  
            else:
                J_2_mn = J_2_dict_OW[(m, n)]  



            A_1_mn = (mu_0 * J_1_mn) / (((m * np.pi) / (w_w))**2 + ((n * np.pi) / (h_w))**2)
            A_2_mn = (mu_0 * J_2_mn) / (((m * np.pi) / (w_w))**2 + ((n * np.pi) / (h_w))**2)
            A_mn = A_1_mn + A_2_mn
            J_mn = J_1_mn + J_2_mn

            
            J_til_mn_temp = 0
            J_til_mn = (0 + 1) * J_mn


            for p in range(1, M):
                if (m + p) % 2 == 1:  
                    
                    if (p, n) not in J_1_dict_OW:
                        J_1_pn = ((4 * J_1) / (p * n * (np.pi)**2)) * (np.sin((p * np.pi * a_1_plus) / w_w) - np.sin((p * np.pi * a_1_minus) / w_w)) * (np.sin((n * np.pi * h_1_plus) / h_w) - np.sin((n * np.pi * h_1_minus) / h_w))
                        J_1_dict_OW[(p, n)] = J_1_pn  
                    else:
                        J_1_pn = J_1_dict_OW[(p, n)]  

                    
                    if (p, n) not in J_2_dict_OW:
                        J_2_pn = ((4 * J_2) / (p * n * (np.pi)**2)) * (np.sin((p * np.pi * a_2_plus) / w_w) - np.sin((p * np.pi * a_2_minus) / w_w)) * (np.sin((n * np.pi * h_2_plus) / h_w) - np.sin((n * np.pi * h_2_minus) / h_w))
                        J_2_dict_OW[(p, n)] = J_2_pn  
                    else:
                        J_2_pn = J_2_dict_OW[(p, n)]

                    J_pn = J_1_pn + J_2_pn
                    J_til_mn_temp += (8 / (np.pi**2)) * (0.5 * (h_w**2 / (p**2 * h_w**2 + n**2 * w_w**2)) - (m**2 + p**2) / ((m**2 - p**2)**2)) * J_pn

            sum_mn += A_mn * (J_til_mn + J_til_mn_temp)


    ### Indutância de Dispersão p.u.a. (Fora da Janela)
    Lk_pua = ((h_w*(w_w**2)) / (4 * I_ref**2)) * (sum_m0 + sum_0n + 0.5*sum_mn)     # (4.38)

    print('')
    print('Indutância de Dispersão p.u.a. dos cantos                ', Lk_pua, 'H / p.u.a.')

    ### Contribuição de Indutância de Dispersão dos cantos:
    print('Contribuição de Indutância de Dispersão dos cantos:      ', Lk_pua*2*np.pi, 'H')

    print('')
    print('')
    print('Indutância de Dispersão total:                           ', Lk_pua*2*np.pi + Lk_pul_OW*b_leg*2 + Lk_pul_IW*a_leg*2, 'H')
        
    in_leak = {'Leakage p.u.l. fora da janela': Lk_pul_IW}
    in_contribuition = {'Contribuicao de Indutancia de Dispersao fora da janela': Lk_pul_IW*a_leg*2}
    out_leak = {'Leakage p.u.l. fora da janela': Lk_pul_OW}
    out_contribuition = {'Contribuicao de Indutancia de Dispersao fora da janela': Lk_pul_OW*b_leg*2}
    side_leak = {'Leakage p.u.a. dos cantos': Lk_pua}
    side_contribuition = {'Contribuicao de Indutancia de Dispersao dos cantos': Lk_pua*2*np.pi}
    total = {'Indutancia de Dispersao total': Lk_pua*2*np.pi + Lk_pul_OW*b_leg*2 + Lk_pul_IW*a_leg*2}

    with open('pyfile_results.json', 'w') as f:
        results = [
            out_leak,
            out_contribuition,
            in_leak,
            in_contribuition,
            side_leak,
            side_contribuition,
            total
        ]
        all = {'RNC': results}
        json.dump(all, f, indent=4)

def circular():
    #####   INSIDE WINDOW   #####
    # Diâmetro da perna central
    d_leg = 0.0269

    ### Largura e Altura da janela
    w_w = 0.0310    # Largura   (x)
    h_w = 0.1010    # Altura    (y)


    ### PRIMÁRIO
    N_1 = 20        # Número de voltas
    I_1 = 1         # Corrente
    a_1 = 0.0041    # Largura
    h_1 = 0.0720    # Altura

    J_1 = (N_1 * I_1) / (a_1 * h_1)  # Densidade de corrente


    ### SECUNDÁRIO
    N_2 = 1350        # Número de voltas
    I_2 = -N_1/N_2    # Corrente (deve ter sentido oposto ao primário)
    a_2 = 0.0043      # Largura
    h_2 = 0.0720      # Altura

    J_2 = (N_2 * I_2) / (a_2 * h_2)  # Densidade de corrente


    ### Disposição dos enrolamentos na janela

    h_1_minus = 0.016
    h_1_plus = 0.016 + h_1

    h_2_minus = 0.016
    h_2_plus = 0.016 + h_2

    a_1_minus = 0.0031
    a_1_plus = 0.0031 + a_1

    a_2_minus = 0.0031 + a_1 + 0.0077
    a_2_plus = 0.0031 + a_1 + 0.0077 + a_2

    # Número de iterações dos somatórios
    M = N = 150 

    # Variáveis dos somatórios
    sum_m0 = sum_0n = sum_mn = 0 


    J_1_dict_IW = {}
    J_2_dict_IW = {}


    for m in range(1, M):


        if (m, 0) not in J_1_dict_IW:
            J_1_m0 = ((2 * J_1) / (m * h_w * np.pi)) * (np.sin((m * np.pi * a_1_plus) / w_w) - np.sin((m * np.pi * a_1_minus) / w_w)) * (h_1_plus - h_1_minus)
            J_1_dict_IW[(m, 0)] = J_1_m0  
        else:
            J_1_m0 = J_1_dict_IW[(m, 0)] 


        if (m, 0) not in J_2_dict_IW:
            J_2_m0 = ((2 * J_2) / (m * h_w * np.pi)) * (np.sin((m * np.pi * a_2_plus) / w_w) - np.sin((m * np.pi * a_2_minus) / w_w)) * (h_2_plus - h_2_minus)
            J_2_dict_IW[(m, 0)] = J_2_m0  
        else:
            J_2_m0 = J_2_dict_IW[(m, 0)]  



        A_1_m0 = (mu_0 * J_1_m0) / (((m * np.pi) / (w_w))**2)
        A_2_m0 = (mu_0 * J_2_m0) / (((m * np.pi) / (w_w))**2)

        A_m0 = A_1_m0 + A_2_m0
        J_m0 = J_1_m0 + J_2_m0

    
        J_til_m0_temp = 0
        J_til_m0 = (d_leg/w_w + 1)*J_m0 


        for p in range(1, M):
            if (m + p) % 2 == 1: 
                if (p, 0) not in J_1_dict_IW:
                    J_1_p0 = ((2 * J_1) / (p * h_w * np.pi)) * (np.sin((p * np.pi * a_1_plus) / w_w) - np.sin((p * np.pi * a_1_minus) / w_w)) * (h_1_plus - h_1_minus)
                else:
                    J_1_p0 = J_1_dict_IW[(p, 0)]
            
                if (p, 0) not in J_2_dict_IW:
                    J_2_p0 = ((2 * J_2) / (p * h_w * np.pi)) * (np.sin((p * np.pi * a_2_plus) / w_w) - np.sin((p * np.pi * a_2_minus) / w_w)) * (h_2_plus - h_2_minus)
                else:
                    J_2_p0 = J_2_dict_IW[(p, 0)]

                J_p0 = J_1_p0 + J_2_p0

                J_til_m0_temp += (8/(np.pi**2)) * (1/(2*p**2) - (m**2 + p**2)/((m**2 - p**2)**2)) * J_p0



        sum_m0 += A_m0 * (J_til_m0 + J_til_m0_temp)


    for n in range(1, N):
        
        if (0, n) not in J_1_dict_IW:
            J_1_0n = ((2 * J_1) / (n * w_w * np.pi)) * (np.sin((n * np.pi * h_1_plus) / h_w) - np.sin((n * np.pi * h_1_minus) / h_w)) * (a_1_plus - a_1_minus)
            J_1_dict_IW[(0, n)] = J_1_0n  
        else:
            J_1_0n = J_1_dict_IW[(0, n)]  
        
        if (0, n) not in J_2_dict_IW:
            J_2_0n = ((2 * J_2) / (n * w_w * np.pi)) * (np.sin((n * np.pi * h_2_plus) / h_w) - np.sin((n * np.pi * h_2_minus) / h_w)) * (a_2_plus - a_2_minus)
            J_2_dict_IW[(0, n)] = J_2_0n  
        else:
            J_2_0n = J_2_dict_IW[(0, n)]  
        

        A_1_0n = (mu_0 * J_1_0n) / (((n * np.pi) / (h_w))**2)
        A_2_0n = (mu_0 * J_2_0n) / (((n * np.pi) / (h_w))**2)
        
        A_0n = A_1_0n + A_2_0n
        J_0n = J_1_0n + J_2_0n
        

        J_til_0n_temp = 0
        J_til_0n = (d_leg/w_w + 1)*J_0n 
        
        
        for p in range(1, M):
            if p % 2 == 1:  
                    
                    if (p, n) not in J_1_dict_IW:
                        J_1_pn = ((4 * J_1) / (p * n * (np.pi)**2)) * (np.sin((p * np.pi * a_1_plus) / w_w) - np.sin((p * np.pi * a_1_minus) / w_w)) * (np.sin((n * np.pi * h_1_plus) / h_w) - np.sin((n * np.pi * h_1_minus) / h_w))
                        J_1_dict_IW[(p, n)] = J_1_pn  
                    else:
                        J_1_pn = J_1_dict_IW[(p, n)]  

                    
                    if (p, n) not in J_2_dict_IW:
                        J_2_pn = ((4 * J_2) / (p * n * (np.pi)**2)) * (np.sin((p * np.pi * a_2_plus) / w_w) - np.sin((p * np.pi * a_2_minus) / w_w)) * (np.sin((n * np.pi * h_2_plus) / h_w) - np.sin((n * np.pi * h_2_minus) / h_w))
                        J_2_dict_IW[(p, n)] = J_2_pn  
                    else:
                        J_2_pn = J_2_dict_IW[(p, n)]  

                    J_pn = J_1_pn + J_2_pn

                    J_til_0n_temp += (8/(np.pi**2)) * (h_w**2/(p**2 * h_w**2 + n**2 * w_w**2) - (1/p**2)*(1 + (n**2 * w_w**2) / (p**2 * h_w**2 + n**2 * w_w**2))) * J_pn* 0.5

        
        sum_0n += A_0n * (J_til_0n + J_til_0n_temp)


    for m in range(1, M):
        for n in range(1, N):

            if (m, n) not in J_1_dict_IW:
                J_1_mn = ((4 * J_1) / (m * n * (np.pi)**2)) * (np.sin((m * np.pi * a_1_plus) / w_w) - np.sin((m * np.pi * a_1_minus) / w_w)) * (np.sin((n * np.pi * h_1_plus) / h_w) - np.sin((n * np.pi * h_1_minus) / h_w))
                J_1_dict_IW[(m, n)] = J_1_mn  
            else:
                J_1_mn = J_1_dict_IW[(m, n)]  

            if (m, n) not in J_2_dict_IW:
                J_2_mn = ((4 * J_2) / (m * n * (np.pi)**2)) * (np.sin((m * np.pi * a_2_plus) / w_w) - np.sin((m * np.pi * a_2_minus) / w_w)) * (np.sin((n * np.pi * h_2_plus) / h_w) - np.sin((n * np.pi * h_2_minus) / h_w))
                J_2_dict_IW[(m, n)] = J_2_mn  
            else:
                J_2_mn = J_2_dict_IW[(m, n)]  



            A_1_mn = (mu_0 * J_1_mn) / (((m * np.pi) / (w_w))**2 + ((n * np.pi) / (h_w))**2)
            A_2_mn = (mu_0 * J_2_mn) / (((m * np.pi) / (w_w))**2 + ((n * np.pi) / (h_w))**2)
            A_mn = A_1_mn + A_2_mn
            J_mn = J_1_mn + J_2_mn

            
            J_til_mn_temp = 0
            J_til_mn = (d_leg / w_w + 1) * J_mn


            for p in range(1, M):
                if (m + p) % 2 == 1:  
                    
                    if (p, n) not in J_1_dict_IW:
                        J_1_pn = ((4 * J_1) / (p * n * (np.pi)**2)) * (np.sin((p * np.pi * a_1_plus) / w_w) - np.sin((p * np.pi * a_1_minus) / w_w)) * (np.sin((n * np.pi * h_1_plus) / h_w) - np.sin((n * np.pi * h_1_minus) / h_w))
                        J_1_dict_IW[(p, n)] = J_1_pn  
                    else:
                        J_1_pn = J_1_dict_IW[(p, n)]  

                    
                    if (p, n) not in J_2_dict_IW:
                        J_2_pn = ((4 * J_2) / (p * n * (np.pi)**2)) * (np.sin((p * np.pi * a_2_plus) / w_w) - np.sin((p * np.pi * a_2_minus) / w_w)) * (np.sin((n * np.pi * h_2_plus) / h_w) - np.sin((n * np.pi * h_2_minus) / h_w))
                        J_2_dict_IW[(p, n)] = J_2_pn  
                    else:
                        J_2_pn = J_2_dict_IW[(p, n)]

                    J_pn = J_1_pn + J_2_pn
                    J_til_mn_temp += (8 / (np.pi**2)) * (0.5 * (h_w**2 / (p**2 * h_w**2 + n**2 * w_w**2)) - (m**2 + p**2) / ((m**2 - p**2)**2)) * J_pn

            sum_mn += A_mn * (J_til_mn + J_til_mn_temp)

    # Corrente de referência 
    I_ref = I_1 

    ### Indutância de Dispersão p.u.a. (Dentro da Janela)
    Lk_pua_IW = ((h_w*(w_w**2)) / (4 * I_ref**2)) * (sum_m0 + sum_0n + 0.5*sum_mn)  # (4.38)

    #  Seção angular que está dentro da janela (radianos)
    alpha = 1.69835/2

    ### Indutância de Dispersão p.u.a. (Dentro da Janela)
    print('Indutância de Dispersão p.u.a. dentro da janela:         ', Lk_pua_IW, 'H / p.u.a.')


    ### Contribuição de Indutância de Dispersão Dentro da Janela:
    print('Contribuição de Indutância de Dispersão dentro da janela:', Lk_pua_IW*alpha, 'H')
    print('')


    ##### OUTSIDE WINDOW #####

    ### "JANELA INFINITA"

    w_w = (0.031+0.101)*5   # Largura   (x)
    h_w = (0.031+0.101)*5   # Altura    (y)


    ### OS ENROLAMENTOS DEVEM ESTAR NO "CENTRO" DA "JANELA INFINITA" 

    h_1_minus = h_w/2 - h_1/2
    h_1_plus = h_w/2 + h_1/2

    h_2_minus = h_w/2 - h_2/2
    h_2_plus = h_w/2 + h_2/2

    # Número de iterações dos somatórios
    M = N = 150

    sum_m0 = sum_0n = sum_mn = 0

    J_1_dict_OW = {}
    J_2_dict_OW = {}


    for m in range(1, M):


        if (m, 0) not in J_1_dict_OW:
            J_1_m0 = ((2 * J_1) / (m * h_w * np.pi)) * (np.sin((m * np.pi * a_1_plus) / w_w) - np.sin((m * np.pi * a_1_minus) / w_w)) * (h_1_plus - h_1_minus)
            J_1_dict_OW[(m, 0)] = J_1_m0  
        else:
            J_1_m0 = J_1_dict_OW[(m, 0)] 


        if (m, 0) not in J_2_dict_OW:
            J_2_m0 = ((2 * J_2) / (m * h_w * np.pi)) * (np.sin((m * np.pi * a_2_plus) / w_w) - np.sin((m * np.pi * a_2_minus) / w_w)) * (h_2_plus - h_2_minus)
            J_2_dict_OW[(m, 0)] = J_2_m0  
        else:
            J_2_m0 = J_2_dict_OW[(m, 0)]  



        A_1_m0 = (mu_0 * J_1_m0) / (((m * np.pi) / (w_w))**2)
        A_2_m0 = (mu_0 * J_2_m0) / (((m * np.pi) / (w_w))**2)

        A_m0 = A_1_m0 + A_2_m0
        J_m0 = J_1_m0 + J_2_m0

    
        J_til_m0_temp = 0
        J_til_m0 = (d_leg/w_w + 1)*J_m0 


        for p in range(1, M):
            if (m + p) % 2 == 1: 
                if (p, 0) not in J_1_dict_OW:
                    J_1_p0 = ((2 * J_1) / (p * h_w * np.pi)) * (np.sin((p * np.pi * a_1_plus) / w_w) - np.sin((p * np.pi * a_1_minus) / w_w)) * (h_1_plus - h_1_minus)
                else:
                    J_1_p0 = J_1_dict_OW[(p, 0)]
            
                if (p, 0) not in J_2_dict_OW:
                    J_2_p0 = ((2 * J_2) / (p * h_w * np.pi)) * (np.sin((p * np.pi * a_2_plus) / w_w) - np.sin((p * np.pi * a_2_minus) / w_w)) * (h_2_plus - h_2_minus)
                else:
                    J_2_p0 = J_2_dict_OW[(p, 0)]

                J_p0 = J_1_p0 + J_2_p0

                J_til_m0_temp += (8/(np.pi**2)) * (1/(2*p**2) - (m**2 + p**2)/((m**2 - p**2)**2)) * J_p0



        sum_m0 += A_m0 * (J_til_m0 + J_til_m0_temp)


    for n in range(1, N):
        
        if (0, n) not in J_1_dict_OW:
            J_1_0n = ((2 * J_1) / (n * w_w * np.pi)) * (np.sin((n * np.pi * h_1_plus) / h_w) - np.sin((n * np.pi * h_1_minus) / h_w)) * (a_1_plus - a_1_minus)
            J_1_dict_OW[(0, n)] = J_1_0n  
        else:
            J_1_0n = J_1_dict_OW[(0, n)]  
        
        if (0, n) not in J_2_dict_OW:
            J_2_0n = ((2 * J_2) / (n * w_w * np.pi)) * (np.sin((n * np.pi * h_2_plus) / h_w) - np.sin((n * np.pi * h_2_minus) / h_w)) * (a_2_plus - a_2_minus)
            J_2_dict_OW[(0, n)] = J_2_0n  
        else:
            J_2_0n = J_2_dict_OW[(0, n)]  
        

        A_1_0n = (mu_0 * J_1_0n) / (((n * np.pi) / (h_w))**2)
        A_2_0n = (mu_0 * J_2_0n) / (((n * np.pi) / (h_w))**2)
        
        A_0n = A_1_0n + A_2_0n
        J_0n = J_1_0n + J_2_0n
        

        J_til_0n_temp = 0
        J_til_0n = (d_leg/w_w + 1)*J_0n 
        
        
        for p in range(1, M):
            if p % 2 == 1:  
                    
                    if (p, n) not in J_1_dict_OW:
                        J_1_pn = ((4 * J_1) / (p * n * (np.pi)**2)) * (np.sin((p * np.pi * a_1_plus) / w_w) - np.sin((p * np.pi * a_1_minus) / w_w)) * (np.sin((n * np.pi * h_1_plus) / h_w) - np.sin((n * np.pi * h_1_minus) / h_w))
                        J_1_dict_OW[(p, n)] = J_1_pn  
                    else:
                        J_1_pn = J_1_dict_OW[(p, n)]  

                    
                    if (p, n) not in J_2_dict_OW:
                        J_2_pn = ((4 * J_2) / (p * n * (np.pi)**2)) * (np.sin((p * np.pi * a_2_plus) / w_w) - np.sin((p * np.pi * a_2_minus) / w_w)) * (np.sin((n * np.pi * h_2_plus) / h_w) - np.sin((n * np.pi * h_2_minus) / h_w))
                        J_2_dict_OW[(p, n)] = J_2_pn  
                    else:
                        J_2_pn = J_2_dict_OW[(p, n)]  

                    J_pn = J_1_pn + J_2_pn

                    J_til_0n_temp += (8/(np.pi**2)) * (h_w**2/(p**2 * h_w**2 + n**2 * w_w**2) - (1/p**2)*(1 + (n**2 * w_w**2) / (p**2 * h_w**2 + n**2 * w_w**2))) * J_pn* 0.5

        
        sum_0n += A_0n * (J_til_0n + J_til_0n_temp)


    for m in range(1, M):
        for n in range(1, N):

            if (m, n) not in J_1_dict_OW:
                J_1_mn = ((4 * J_1) / (m * n * (np.pi)**2)) * (np.sin((m * np.pi * a_1_plus) / w_w) - np.sin((m * np.pi * a_1_minus) / w_w)) * (np.sin((n * np.pi * h_1_plus) / h_w) - np.sin((n * np.pi * h_1_minus) / h_w))
                J_1_dict_OW[(m, n)] = J_1_mn  
            else:
                J_1_mn = J_1_dict_OW[(m, n)]  

            if (m, n) not in J_2_dict_OW:
                J_2_mn = ((4 * J_2) / (m * n * (np.pi)**2)) * (np.sin((m * np.pi * a_2_plus) / w_w) - np.sin((m * np.pi * a_2_minus) / w_w)) * (np.sin((n * np.pi * h_2_plus) / h_w) - np.sin((n * np.pi * h_2_minus) / h_w))
                J_2_dict_OW[(m, n)] = J_2_mn  
            else:
                J_2_mn = J_2_dict_OW[(m, n)]  



            A_1_mn = (mu_0 * J_1_mn) / (((m * np.pi) / (w_w))**2 + ((n * np.pi) / (h_w))**2)
            A_2_mn = (mu_0 * J_2_mn) / (((m * np.pi) / (w_w))**2 + ((n * np.pi) / (h_w))**2)
            A_mn = A_1_mn + A_2_mn
            J_mn = J_1_mn + J_2_mn

            
            J_til_mn_temp = 0
            J_til_mn = (d_leg / w_w + 1) * J_mn


            for p in range(1, M):
                if (m + p) % 2 == 1:  
                    
                    if (p, n) not in J_1_dict_OW:
                        J_1_pn = ((4 * J_1) / (p * n * (np.pi)**2)) * (np.sin((p * np.pi * a_1_plus) / w_w) - np.sin((p * np.pi * a_1_minus) / w_w)) * (np.sin((n * np.pi * h_1_plus) / h_w) - np.sin((n * np.pi * h_1_minus) / h_w))
                        J_1_dict_OW[(p, n)] = J_1_pn  
                    else:
                        J_1_pn = J_1_dict_OW[(p, n)]  

                    
                    if (p, n) not in J_2_dict_OW:
                        J_2_pn = ((4 * J_2) / (p * n * (np.pi)**2)) * (np.sin((p * np.pi * a_2_plus) / w_w) - np.sin((p * np.pi * a_2_minus) / w_w)) * (np.sin((n * np.pi * h_2_plus) / h_w) - np.sin((n * np.pi * h_2_minus) / h_w))
                        J_2_dict_OW[(p, n)] = J_2_pn  
                    else:
                        J_2_pn = J_2_dict_OW[(p, n)]

                    J_pn = J_1_pn + J_2_pn
                    J_til_mn_temp += (8 / (np.pi**2)) * (0.5 * (h_w**2 / (p**2 * h_w**2 + n**2 * w_w**2)) - (m**2 + p**2) / ((m**2 - p**2)**2)) * J_pn

            sum_mn += A_mn * (J_til_mn + J_til_mn_temp)



    ### Indutância de Dispersão p.u.a. (Fora da Janela)
    Lk_pua_OW = ((h_w*(w_w**2)) / (4 * I_ref**2)) * (sum_m0 + sum_0n + 0.5*sum_mn)  # (4.38)


    ### Indutância de Dispersão p.u.a. (Fora da Janela)
    print('Indutância de Dispersão p.u.a. fora da janela:           ', Lk_pua_OW, 'H / p.u.a.')


    ### Contribuição de Indutância de Dispersão Dentro da Janela:
    print('Contribuição de Indutância de Dispersão fora da janela:  ', Lk_pua_OW*(2*np.pi-alpha), 'H')

    print('Indutância de Dispersão total:                           ', Lk_pua_OW*(2*np.pi-alpha) + Lk_pua_IW*alpha, 'H')

ini = time.time()
nonCircular()
fim = time.time()
print('Tempo de execução : ', fim-ini, 's')

with open('pyfile_results.json', 'r') as f:
    data = json.load(f)

with open('jupyter_results.json', 'r') as f:
    data2 = json.load(f)

def compare_dicts(list1, list2):
    differences = []
    
    if len(list1) != len(list2):
        differences.append(f"Lists have different lengths: {len(list1)} vs {len(list2)}")
    
    for i, (item1, item2) in enumerate(zip(list1, list2)):
        for key in item1:
            if key not in item2:
                differences.append(f"Key '{key}' is missing in item {i+1} of the second list")
            elif item1[key] != item2[key]:
                differences.append(f"Value for key '{key}' in item {i+1} differs: {item1[key]} vs {item2[key]}")
    
    return differences

differences = compare_dicts(data['RNC'], data2['RNC'])

if differences:
    print("\nDifferences found:")
    for diff in differences:
        print(diff)
else:
    print("No differences found.")
