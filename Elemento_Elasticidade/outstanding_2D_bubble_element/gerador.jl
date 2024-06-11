function Gera_malha(Lx::Float64, nx::Int64, Ly::Float64, ny::Int64, force=1.0;
                     origin=(0.0,0.0))

    # Assertions
    @assert Lx>0  "Bmesh_solid_2D:: Lx must be > 0"
    @assert Ly>0  "Bmesh_solid_2D:: Ly must be > 0"
    @assert nx>=1 "Bmesh_solid_2D:: nx must be >= 1"
    @assert ny>=1 "Bmesh_solid_2D:: ny must be >= 1"

    # Primeiro geramos os nós de baixo para cima, esquerda para a direita
    nn = (nx+1)*(ny+1)
    coord = zeros(nn,2)

    # Dimensões de cada elemento
    dx = Lx / nx
    dy = Ly / ny

    # Gera a matrix de coordenadas nodais
    x = origin[1]-dx
    y = origin[2]-dy
    cont = 0
    for j=1:ny+1
        y += dy
        for i=1:nx+1
            x += dx
            cont += 1
            coord[cont,:] = [x y]
        end #i
        x = origin[1]-dx
    end #j

    # Gera a matrix de conectividades
    ne = (nx)*(ny)
    connect = zeros(Int64,ne,4)
    no1 = 1
    no2 = 2
    no3 = (nx+1)+2
    no4 = no3-1
    cont = 0
    for j=1:ny
        for i=1:nx
            cont += 1
            connect[cont,:] = [no1 no2 no3 no4]
            no1 += 1; no2 += 1; no3 += 1; no4 += 1
        end #i
        no1 += 1; no2 += 1; no3 += 1; no4 += 1
    end #j

    # Condições de contorno essenciais (engaste)
    nhebc = 2*(ny+1)
    hebc = Array{Float64}(undef,nhebc,3)
    node = 1
    pos = 0
    for i=1:(ny+1)
        pos += 1
        hebc[pos,:] = [node 1 0.0]
        pos += 1
        hebc[pos,:] = [node 2 0.0]
        node += (nx+1)
    end
    
    # Generate the load information
    no_forca = nx+1
    nbc = [no_forca 2 -force]

    return nn, ne, coord, connect, hebc, nbc

end
