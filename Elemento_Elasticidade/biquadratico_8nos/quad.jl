###################################################################

#         ELEMENTO BI QUADRÁTICO DE 8 NÓS (SERENDIPITY)

###################################################################

# Funções de interpolação
# N1 = 0.25*(1-r)*(1-s)*(-s-r-1)
# N2 = 0.25*(r+1)*(1-s)*(-s+r-1)
# N3 = 0.25*(r+1)*(s+1)*(s+r-1)
# N4 = 0.25*(1-r)*(s+1)*(s-r-1)
# N5 = 0.5*(1-r^2)*(1-s)
# N6 = 0.5*(r+1)*(1-s^2)
# N7 = 0.5*(1-r^2)*(s+1)
# N8 = 0.5*(1-r)*(1-s^2)


#
# Calcula a matriz B e o determinante da matriz
# Jacobiana no ponto (r,s) 
#
function B_bi8n(r,s,X,Y)

   # Vetores com as derivadas das funções de interpolação
   # em relação a r e a s
   dNr = [ -(s^2+(2*r-1)*s-2*r)/4;  (s^2+(-2*r-1)*s+2*r)/4; 
            (s^2+(2*r+1)*s+2*r)/4;  -(s^2+(1-2*r)*s-2*r)/4;
            r*s-r;    (s^2-1)/2;    -r*s-r;      (s^2-1)/2]

   dNs = [  -((2*r-2)*s+r^2-r)/4;   ((2*r+2)*s-r^2-r)/4  ; 
             ((2*r+2)*s+r^2+r)/4;    -((2*r-2)*s-r^2+r)/4; 
             (r^2-1)/2;   (-r-1)*s;  -(r^2-1)/2;  (r-1)*s]

   # Aloca e calcula a matriz J
   J = zeros(2,2)

   # Loop do somatório para cada posição de J
   for i=1:8
       J[1,1] += dNr[i]*X[i]
       J[1,2] += dNr[i]*Y[i]
       J[2,1] += dNs[i]*X[i]
       J[2,2] += dNs[i]*Y[i]
   end #i

   # Calcula o determinante
   dJ = det(J)

   # Calcula a inversa da J
   invJ = inv(J)

   # Aloca a matriz B
   B = zeros(3,16)

   # Loop pelos blocos da B, com as correções de derivadas
   c = 1
   for i=1:8
     
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
# biquadrático isoparamétrico de 8 nós
#
function K_bi8n(C::AbstractMatrix,X::Vector,Y::Vector,esp::Float64)

    # Aloca a matriz do elemento
    K = zeros(16,16)

    # Define a quadratura 
    pg = [-sqrt(0.6) ; -sqrt(0.6); 0]
    W  = [5/9 ; 5/9; 8/9]

    # Loop pelos pontos de Gauss
    # Ao inves de calcular para os infinitos pontos do elemento
    # podemos calcular apenas nos pontos de Gauss
    for i=1:3
        r  = pg[i]
        wr = W[i]
        for j=1:3
           s  = pg[j]  
           ws = W[j]

           # calcula o dJ e o B
           dJ, B = B_bi8n(r,s,X,Y)

           # Acumula o produto nesses ptos
           K = K + B'*C*B*dJ*esp

        end #j
    end # i

    # Retorna a matriz de rigidez do elmemento
    return K

end
