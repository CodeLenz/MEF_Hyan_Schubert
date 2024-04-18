#
# conectividades é uma matriz com ne linhas (número de elementos)
# e 2 colunas (nó inicial e nó final)
# VE, VA, VL e Vtheta são vetores de dimensão ne com 
# os E, A, L e theta de cada elemento
#
function Rigidez_global(ne,nnos,conectividades, VE, VA, VL, Vtheta)

    # Precisamos definir a matriz global
    K = zeros(2*nnos,2*nnos)

    # Loop nos elementos da malha
    for ele=1:ne

        # Recupera as informações do elemento
        Ee = VE[ele]
        Ae = VA[ele]
        Le = VL[ele]
        te = Vtheta[ele]

        # Monta a matriz de rigidez do elemento
        # no sistema local
        Kl = Rigidez_barra(Ee,Ae,Le)

        # Monta a matriz de transformação T
        T = Matriz_rotacao(te)
        
        # Passa Kl para o sistema global
        Kg = transpose(T)*Kl*T

        # Agora precisamos posicionar Kg na matriz
        # global do problema

        # Recupera os nós do elemento
        no1 = conectividades[ele,1]
        no2 = conectividades[ele,2]

        # Vetor com os gls GLOBAIS do elemento
        gls = [ 2*(no1-1)+1 ;2*(no1-1)+2; 2*(no2-1)+1; 2*(no2-1)+2 ]

        # Soma Kg nas posições gls em K
        K[gls,gls] .= K[gls,gls] .+ Kg

    end #ele

    # Retorna a matriz de rigidez do problema
    return K

end

# Montar o vetor de forças global
function Forca_global(nnos,forcas)

    # Vetor de forças global
    F = zeros(2*nnos)
    
    # Posiciona os valores de força no vetor 
    # de forças global
    for i=1:size(forcas,1)
        # Garante que gl é um inteiro
        # gl global é calculado pela equação de Localização
        gl = 2*(Int(forcas[i,1])-1) + Int(forcas[i,2])
        valor = forcas[i,3]
        F[gl] = valor
    end #i

    # Retorna a matriz de força global do problema
    return F
    
end