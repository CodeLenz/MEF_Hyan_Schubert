using LinearAlgebra
using SparseArrays

include("gerador.jl")
include("material.jl")
include("quad_enriquecido.jl")
include("apoios.jl")
include("global.jl")

function exemplo_cis_parasita(nx,ny,bolha=true)


    Lx = 1.0
    Ly = 0.1
    espessura = 0.1
    force = 1000.0
    
    # Gera a malha da viga longa
    nnos, ne, coord, conectividades,apoios, forcas = Gera_malha(Lx, nx, Ly, ny, force)

    E=210E9
    VE = E*ones(ne)
    Vnuxy = (0.0)*ones(ne)
    
    # Deslocamento máximo ANALÍTICO
    Iz = (espessura*Ly^3)/12
    Ua_max = -(force*Lx^3)/(3*E*Iz)

    # Monta a matriz global do problema
    K = Global(ne,nnos,conectividades,coord, VE, Vnuxy, espessura,bolha=bolha)

    # Monta o vetor de forças globais concentradas
    F = Forca_global(nnos,forcas)

    # Modifica o sistema pela aplicação das condições de contorno
    # homogêneas
    # Aplica_CC_homo!(apoios,K,F)
    KA, FA = Aplica_CC_lagrange(nnos,apoios,K,F)

    # Soluciona o sistema de equações, obtendo os deslocamentos
    # KU = F
    UA = KA\FA

    # Só os deslocamentos que nos interessam
    U = UA[1:2*nnos]

    # Deflexão máxima NUMÉRICA
    Un_max=U[end]

    # Retorna U
    return Un_max, Ua_max

end


#
# Rotina de teste
#
function testa_bagaca(n=10)

    # Com bolha
    N = Int64[]
    V = Float64[] 
    for i=1:n
        num,ana=exemplo_cis_parasita(4*i,i,true)
        push!(N,i)
        push!(V,num/ana)
    end
    display(plot(N,V,label="Bolha",xlabel="Número de elementos",ylabel="Deslocamento/referência",linewidth=2))

    # Sem bolha
    N = Int64[]
    V = Float64[] 
    for i=1:n
        num,ana=exemplo_cis_parasita(4*i,i,false)
        push!(N,i)
        push!(V,num/ana)
    end
    display(plot!(N,V,label="Bilinear",xlabel="Número de elementos",ylabel="Deslocamento/referência",linewidth=2))

    # Grava a figura
    savefig("comparacao.png")

end