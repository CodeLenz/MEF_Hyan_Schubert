   # N1 = 0.25(1-r)(1-s)
   # N2 = 0.25(1+r)(1-s)
   # N3 = 0.25(1+r)(1+s)
   # N4 = 0.25(1-r)(1+s)

   # Elemento de bolha
   # N5 = 1-r^2
   # N6 + 1-s^2

#
# Calcula a matriz B e o determinante da matriz
# Jacobiana no ponto (r,s) 
#
function B_bolha(r,s,X,Y)

   # Vetores com as derivadas das funções de interpolação
   # em relação a r e a s                                         # Bolha
   dNr = [ -0.25*(1-s) ; 0.25*(1-s) ; 0.25*(1+s)  ; -0.25*(1+s)  ; -2*r ; 0 ]
   dNs = [ -0.25*(1-r) ;-0.25*(1+r) ; 0.25*(1+r)  ;  0.25*(1-r)  ; 0 ; -2*s ]

   # Aloca e calcula a matriz J
   J = zeros(2,2)

   # Loop do somatório para cada posição de J
   # Até o 4° termo pois o 5° e 6° são referentes à bolha
   for i=1:4
       J[1,1] += dNr[i]*X[i]
       J[1,2] += dNr[i]*Y[i]
       J[2,1] += dNs[i]*X[i]
       J[2,2] += dNs[i]*Y[i]

   end #i

   # Calcula o determinante
   dJ = det(J)

   # Calcula a inversa da J
   invJ = inv(J)

    # Aloca a matriz B original
    B = zeros(3,12)

   # Loop pelos blocos da B, com as correções de derivadas
   c = 1
   for i=1:6
     
       # Correção das derivadas para este bloco
       dNxy = invJ*[dNr[i];dNs[i]]

       # Posiciona na B
       B[1,c]   = dNxy[1]
       B[2,c+1] = dNxy[2]
       B[3,c]   = dNxy[2]
       B[3,c+1] = dNxy[1]

       # Atualiza a apontador de bloco
       c += 2

   end #i

   return dJ, B
end



#
# Monta a matriz de rigidez de um elemento
# bolha de 4 nós
#
include("guyan.jl")

function K_bolha(C::AbstractMatrix,X::Vector,Y::Vector,esp::Float64; bolha)

    # Aloca a matriz do elemento
    KA = zeros(12,12)

    # Define a quadratura 
    pg = [-1/sqrt(3) ; 1/sqrt(3)]
    W  = [1.0 ; 1.0]

    # Loop pelos pontos de Gauss
    # Ao inves de calcular para os infinitos pontos do elemento
    # podemos calcular apenas nos pontos de Gauss
    for i=1:2
        r  = pg[i]
        wr = W[i]
        for j=1:2
           s  = pg[j]  
           ws = W[j]

           # calcula o dJ e o B
           dJ, B = B_bolha(r,s,X,Y)

           # Acumula o produto nesses ptos
           KA = KA + B'*C*B*dJ*esp

        end #j
    end # i

    if bolha
        
        # Aplica Condensação Estática (Redução de Guyan)
        # Graus de liberdade mestres
        glm = Int[1;2;3;4;5;6;7;8]
        K = Guyan(KA, glm )

        # Retorna a matriz de rigidez do elmemento
        return K
    
    else
        return KA[1:8,1:8]
    end

end
