using LinearAlgebra
using Plots

include("pre.jl")
include("barra.jl")
include("global.jl")

# mesma função da main desenvovlida em sala. Contudo, os apoios 
# e forças do problema são dadas em função do nó e grau de li-
# berdade local do nó

function main2()

    # Entrada de dados na mão!
    ne = 5
    nnos = 4
    conectividades = [2 1 ;
                      2 3 ;
                      4 3 ;
                      1 4 ;
                      1 3]

    coord = [0.0 0.0 ;
             1.0 0.0 ;
             1.0 1.0 ;
             0.0 1.0]                  

    VE = 210E9*ones(ne)
    VA = 1E-4*ones(ne)

    # Pré-processamento
    VL,Vtheta = Pre_processa(ne,coord,conectividades)

    # Apoios (cond. de contorno essenciais)
    # apoios = [1;2;4]

    # Fornecendo as informações de apoios com base no nó (coluna 1)
    # e grau de liberdade local do nó (coluna 2)
    apoios = [1 1;
              1 2;
              2 2]
          
    # Forças (cond. de contorno naturais)
    #forcas = [3  300.0 ;
              #5  600.0 ;
              #6 -500.0]
    
    # Matriz de forças com base no nó (coluna 1), grau de liberdade
    # local do nó (coluna 2) e valor (coluna 3)
    forcas = [2 1 300.0 ;
              3 1 600.0 ;
              3 2 -500.0]

    # Monta a matriz global do problema
    K = Global(ne,nnos,conectividades, VE, VA, VL, Vtheta)

    # Montar o vetor de forças global
    F = zeros(2*nnos)

    # Para cada informação em forças 
    # posiciona a força no vetor de forças
    # globais do problema
    for i=1:size(forcas,1)
        # Garante que gl é um inteiro
        # gl global é calculado pela equação de Localização
        gl = 2*(Int(forcas[i,1])-1) + Int(forcas[i,2])
        valor = forcas[i,3]
        F[gl] = valor
    end #i

    # Aplica as condições de contorno essenciais homogêneas
    for i=1:size(apoios,1)
        # Transforma gl local no gl global
        gl = 2*(Int(apoios[i,1])-1) + Int(apoios[i,2])

        # Zera linha e coluna da rigidez
        K[gl,:] .= 0.0
        K[:,gl] .= 0.0

        # Zera o vetor de forças nesse gl
        F[gl] = 0.0

        # Coloca 1.0 na diagonal
        K[gl,gl] = 1.0

    end #i gls

    # Soluciona o sistema de equações, obtendo os deslocamentos
    # KU = F
    U = K\F

end #main2()
