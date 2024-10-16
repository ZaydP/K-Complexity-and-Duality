using LinearAlgebra
using DelimitedFiles
using SparseArrays

function comm(a,b)
    return a*b-b*a
end

function InnProd(a, b)
    flat_a = vec(a)
    flat_b = vec(b)
    c = dot(flat_a,flat_b)
    return c
end

function SpinStateMatrix(OccRep, Spin)
    # Given a binary occupation representation of the spin chain, return the associated matrix representation
    id=[1 0; 0 1]
    a = 1
    for i in 1:length(OccRep)
        if OccRep[i] == 0
            a = kron(a, id)
        end
        if OccRep[i] == 1
            a = kron(a, Spin)
        end
    end
    return a
end

function IsingSum(L, InteractionSpin, IntSpinNo, BSpin1, BSpin2)
    # Inputs: Length of Ising Chain, Interaction Spin type, Locality of the Interaction, magnetic field spin 1, and magnetic field spin 2
    # Outputs: Total matrix for the Ising chain
    
    # generate all the possible arrangements of the spins on the chain in binary occupation rep form:
    MatrixWidth = 2^L
    IntSpinPerms = zeros(Int, L-1, L)
    for i = 1:L-1
        IntSpinPerms[i, :] = cat(zeros(Int, i-1), ones(Int, IntSpinNo), zeros(Int, L-(i-1)-IntSpinNo); dims=1)
    end
    
    BSpinPerms = zeros(Int, L, L)
    for i = 1:L
        BSpinPerms[i, :] = cat(zeros(Int, i-1), ones(Int, 1), zeros(Int, L-(i-1)-1); dims=1)
    end
    
    # for each rep, send to SpinStateMatrix to get matrix back, and sum up
    IntAll = zeros(MatrixWidth, MatrixWidth)
    if InteractionSpin != 0
        for i = 1:L-1
            IntAll += SpinStateMatrix(IntSpinPerms[i, :], InteractionSpin)
        end
    end
    
    BSpin1All = zeros(MatrixWidth, MatrixWidth)
    if BSpin1 != 0
        for i = 1:L
            BSpin1All += SpinStateMatrix(BSpinPerms[i, :], BSpin1)
        end
    end
    
    BSpin2All = zeros(MatrixWidth, MatrixWidth)
    if BSpin2 != 0
        for i = 1:L
            BSpin2All += SpinStateMatrix(BSpinPerms[i, :], BSpin2)
        end
    end
    
    IsingMatrix = IntAll + BSpin1All + BSpin2All
    return IsingMatrix
end

function modified_gram_schmidt(V)
    if isempty(V)
        return []
    end

    n, m = size(V)
    Q = similar(V, typeof(V[1]))

    for j in 1:m
        q = copy(V[:, j])
        for i in 1:j-1
            q -= InnProd(Q[:, i], V[:, j]) * Q[:, i]
        end
        Q[:, j] = q / InnProd(q, q)  # Normalize
    end

    return Q
end

function modified_gram_schmidt(V, inner_product_function)
    n, m = size(V)
    Q = similar(V, typeof(V[1]))
    R = zeros(eltype(V), m, m)

    for j in 1:m
        v = copy(V[:, j])
        for i in 1:j-1
            R[i, j] = inner_product_function(Q[:, i], V[:, j])
            v -= R[i, j] * Q[:, i]
        end
        R[j, j] = sqrt(inner_product_function(v, v))
        Q[:, j] = v / R[j, j]
    end

    return Q, R
end


p1 = [0 1; 1 0]
p2 = [0 -1*im;1*im 0]
p3 = [1 0; 0 -1] 

modified_gram_schmidt((p1;p2))

#for L=3 use nmax>34,   and prec>18 . ratio = 0.529
#for L=4 use nmax>124,  and prec>75 . ratio = 0.605
#for L=5 use nmax>516,  and prec>361. ratio = 0.700
#for L=6 use nmax>2016, and prec>?
#for L=7 use nmax>8100?

L=8
nmax=32800
prec=27000
prec_base2=round(Int,prec*3.4)
btol=0.0001

H=sparse(IsingSum(L,p3,2,p3,p1))
Op=sparse(IsingSum(L,0,2,p3,0))

# setprecision(BigFloat,prec,base=10)
setprecision(BigFloat,prec_base2)
@time begin
    bs=LanczosDoc(H,Op,nmax,btol);
end

open("/home/pndzay001/julia_code/julia_outputs/julia_output_L8.txt", "w") do file
    writedlm(file,bs)
end

# bruh=readdlm("julia_output_L6.txt")
