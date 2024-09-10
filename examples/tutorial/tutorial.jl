using XDiag, LinearAlgebra, GLMakie, TOML

let
    # Tutorial displaying Hilbert space
    hspace = Spinhalf(2)
    for s in hspace
        println(s)
    end

    # Tutorial computing eigenvectors, eigenvalues
    ops = OpSum()
    ops += Op("HB", 1.0, [1, 2])

    H = matrix(ops, hspace)
    display(H)

    evals, evecs = eigen(H)
    display(evals)
    display(evecs)

    # Exercise triangular
    hspace = Spinhalf(3)
    ops = OpSum()
    ops += Op("HB", 1.0, [1, 2])
    ops += Op("HB", 1.0, [2, 3])
    ops += Op("HB", 1.0, [3, 1])
    H = matrix(ops, hspace)
    evals, evecs = eigen(H)
    display(evals)
    display(evecs)

    # Tutorial Sz conservation
    for nup in 0:3
        block = Spinhalf(3, nup)
        H = matrix(ops, block)
        evals, evecs = eigen(H)
        @show nup
        display(evals)
        display(evecs)
    end

    # Tutorial translational symmetry
    N = 4
    nup = 2
    ops = OpSum()
    for i in 1:N
        ops += Op("HB", 1.0, [i, mod1(i+1, N)])
    end
    
    id = Permutation([1, 2, 3, 4])
    T = Permutation([2, 3, 4, 1])
    group = PermutationGroup([id, T, T*T, T*T*T])

    # 0-momentum irrep
    irrep = Representation([1, 1, 1, 1])
    block = Spinhalf(N, nup, group, irrep)
    for s in block
        println(s)
    end
    H = matrix(ops, block)
    println("0 momentum")
    display(H)
    display(eigvals(H))

    # pi/2-momentum irrep
    irrep = Representation([1, 1im, -1, -1im])
    block = Spinhalf(N, nup, group, irrep)
    H = matrix(ops, block)
    println("pi/2 momentum")
    display(H)
    display(eigvals(Hermitian(H)))

    # pi-momentum irrep
    irrep = Representation([1, -1, 1, -1])
    block = Spinhalf(N, nup, group, irrep)
    H = matrix(ops, block)
    println("pi momentum")
    display(H)
    display(eigvals(H))

    # Exersize: Get full spectrum
    N = 4
    ops = OpSum()
    for i in 1:N
        ops += Op("HB", 1.0, [i, mod1(i+1, N)])
    end
    block = Spinhalf(N)
    H = matrix(ops, block)
    eigs = eigvals(Symmetric(H))

    id = Permutation([1, 2, 3, 4])
    T = Permutation([2, 3, 4, 1])
    group = PermutationGroup([id, T, T*T, T*T*T])

    ir1 = Representation([1, 1, 1, 1]) 
    ir2 = Representation([1, -1, 1, -1])    
    ir3 = Representation([1, im, -1, -im])    
    ir4 = Representation([1, -im, -1, im])    

    eigs_sym = []
    for nup in 0:N
        for irrep in [ir1, ir2, ir3, ir4]
            block = Spinhalf(N, nup, group, irrep)
            H = matrix(ops, block)
            append!(eigs_sym, eigvals(Hermitian(H)))
        end
    end
    sort!(eigs_sym)
    display(eigs)
    display(eigs_sym)
    @show isapprox(eigs, eigs_sym)

    # Tutorial: Sparse diagonalization and algebra
    N = 10
    ops = OpSum()
    for i in 1:N
        ops += Op("HB", 1.0, [i, mod1(i+1, N)])
    end
    block = Spinhalf(N)
    H = matrix(ops, block)
    eigs = eigvals(Symmetric(H))

    e0, gs = eig0(ops, block)

    Hgs = zeros(block)
    apply(ops, gs, Hgs)
    @show e0, eigs[1]
    @show norm(gs)
    @show dot(gs, Hgs)
    @show inner(ops, gs)

    rstate = rand(block)
    @show inner(ops, rstate)
    
    # Exercise: Implementing the Lanczos algorithm
    N = 10
    ops = OpSum()
    for i in 1:N
        ops += Op("HB", 1.0, [i, mod1(i+1, N)])
    end
    block = Spinhalf(N)
    
    function lanczos(ops, block; precision=1e-12, max_iterations=1000)
    
        alphas = Float64[]
        betas = Float64[]

        v1 = rand(block; normalized=true)
        v0 = zeros(block)
        w = zeros(block)

        alpha = 0.0
        beta = 0.0
            
        prev_energy = 0
        for iteration in 1:max_iterations

            # Perform basic Lanczos iteration step
            apply(ops, v1, w)
            alpha = dot(v1, w)
            w = w - alpha * v1 - beta * v0
            v0 = v1

            beta = norm(w)
            v1 = w / beta

            push!(alphas, alpha)

            # Build up T-matrix
            t_matrix = diagm(0 => alphas, 1 => betas, -1 => betas)
            t_eigs = eigvals(Symmetric(t_matrix))

            push!(betas, beta)
            @show iteration, t_eigs[1]

            # Return if converged
            if isapprox(t_eigs[1], prev_energy; rtol=precision)
                println("Lanczos converged in $iteration steps")
                return t_eigs, alphas, betas
            end
            prev_energy = t_eigs[1]
        end            
    end

    t_eigs, alphas, betas = lanczos(ops, block)
    @show t_eigs[1]
    @show eigval0(ops, block)


    # Tutorial: Plot convergence of Lanczos algorithm
    f = Figure()
    ax = Axis(f[1, 1])
    niter = length(alphas)
    for iter in 1:niter
        t_matrix = diagm(0 => alphas[1:iter], 1 => betas[1:iter-1], -1 => betas[1:iter-1])
        t_eigs = eigvals(Symmetric(t_matrix))
        scatter!(ax, iter*ones(size(t_eigs)), t_eigs)
    end
    ax.xlabel = "iteration"
    ax.ylabel = "eigenvalue"
    display(f)

    # Ground state correlations of triangular lattice Heisenberg model
    N = 12
    latfile = TOML.parsefile("triangular.$N.J1J2.toml")
    
    coords = latfile["Coordinates"]
    interactions = latfile["Interactions"]
    ops = OpSum()
    for interaction in interactions
        type = interaction[1]
        couplingname = interaction[2]
        s1 = interaction[3] + 1
        s2 = interaction[4] + 1
        ops += Op(type, couplingname, [s1, s2])
    end
    ops["J1"] = 1.0
    ops["J2"] = 0.0

    block = Spinhalf(N, N÷2)
    e0, gs = eig0(ops, block)
    @show e0
    
    for i in 2:N
        s1si = Op("HB", 1.0, [1, i])
        corr = inner(s1si, gs)
        println("$i $corr")
    end

    # Exercise: plot correlations on the lattice
    # ...


    # Tutorial: compute ground state correlations using symmetries
    permutation_vecs = latfile["Symmetries"]
    permutations = Permutation[]
    for permutation_vec in permutation_vecs
        push!(permutations, Permutation(permutation_vec .+ 1))
    end
    group = PermutationGroup(permutations)
    irrep = Representation(ones(size(group)))
    block = Spinhalf(N, N÷2, group, irrep)
    e0, gs = eig0(ops, block)
    @show e0

    for i in 2:N
        s1si = symmetrize(Op("HB", 1.0, [1, i]), group)
        corr = inner(s1si, gs)
        println("$i $corr")
    end
    
end
