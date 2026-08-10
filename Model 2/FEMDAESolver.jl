module FEMDAESolver

    using Ferrite, SparseArrays, WriteVTK, Plots, LinearAlgebra, DifferentialEquations

    export solve_fem_absorption, solve_fem_desorption, Constants_Model2

    struct Constants_Model2
        R::Float64
        T::Float64
        μ::Float64
        velocity::Float64
        ϵ::Float64
        MA::Float64
        rho_sat::Float64
        p_eqa::Float64
        ρ_e::Float64
    end

    struct RHSparams
        constants::Constants_Model2
        dh::DofHandler
        cellvalues_rhog::CellValues
        cellvalues_rhos::CellValues
        M::SparseMatrixCSC          # {Float64,Int} ?
        rhs::Vector{Float64}        # Can start as zeros, but we do not want to allocate each time step, so we will reuse it
        rhsL::SparseMatrixCSC       # Linear part of right-hand side
        assemble_rhs!::Function     # Complete assembly of the right-hand side, including the nonlinear part
    end

    function doassemble_M!(
        M_rho::SparseMatrixCSC, 
        cellvalues_rho_g::CellValues, 
        cellvalues_rho_s::CellValues,
        ϵ::Float64,
        dh::DofHandler
        )

        # number of basis functions per field, initialize Me
        n_basefuncs_rho_g  = getnbasefunctions(cellvalues_rho_g)
        n_basefuncs_rho_s  = getnbasefunctions(cellvalues_rho_s)
        first_cell = first(CellIterator(dh))
        ndofs_cell = length(celldofs(first_cell))
        Me = zeros(ndofs_cell, ndofs_cell)

        
        assembler = start_assemble(M_rho)
        rhog_dofs = dof_range(dh, :rho_g)
        rhos_dofs  = dof_range(dh, :rho_s)

        # loop over all cells
        for cell in CellIterator(dh)

            # get local dofs of this cell for rho
            cell_dofs = celldofs(cell)
            fill!(Me,0)
            
            # rho
            Ferrite.reinit!(cellvalues_rho_g, cell)
            for q_point in 1:getnquadpoints(cellvalues_rho_g)
                dx = getdetJdV(cellvalues_rho_g, q_point)
                for i in 1:n_basefuncs_rho_g
                    psi_i = shape_value(cellvalues_rho_g, q_point, i)
                    for j in 1:n_basefuncs_rho_g
                        psi_j = shape_value(cellvalues_rho_g, q_point, j)
                        # pim pam pet
                        Me[rhog_dofs[i], rhog_dofs[j]] += ϵ * psi_i * psi_j * dx
                    end
                end
            end

            # rhos
            Ferrite.reinit!(cellvalues_rho_s, cell)
            for q_point in 1:getnquadpoints(cellvalues_rho_s)
                dx = getdetJdV(cellvalues_rho_s, q_point)
                for i in 1:n_basefuncs_rho_s
                    psi_i = shape_value(cellvalues_rho_s, q_point, i)
                    for j in 1:n_basefuncs_rho_s
                        psi_j = shape_value(cellvalues_rho_s, q_point, j)
                        # pim pam pet
                        Me[rhos_dofs[i], rhos_dofs[j]] += (1-ϵ) * psi_i * psi_j * dx
                    end
                end
            end

            assemble!(assembler, cell_dofs, Me)
        end
        return M_rho
    end

    function doassemble_D!(
        D::SparseMatrixCSC, 
        cellvalues_rho_g::CellValues, 
        cellvalues_rho_s::CellValues, 
        μ::Float64,
        dh::DofHandler
        )

        assembler = start_assemble(D)

        # number of basis functions per field
        n_basefuncs_rho_g  = getnbasefunctions(cellvalues_rho_g)

        for cell in CellIterator(dh)

            # get global DOFs for this cell
            cell_dofs = celldofs(cell)
            ndofs_cell = length(cell_dofs)
            De = zeros(ndofs_cell, ndofs_cell)

            # rho is the only important bit since rho_s does not diffuse
            Ferrite.reinit!(cellvalues_rho_g, cell)

            for q_point in 1:getnquadpoints(cellvalues_rho_g)
                dx = getdetJdV(cellvalues_rho_g, q_point)

                for i in 1:n_basefuncs_rho_g
                    dpsi_i = shape_gradient(cellvalues_rho_g, q_point, i)

                    for j in 1:n_basefuncs_rho_g
                        dpsi_j = shape_gradient(cellvalues_rho_g, q_point, j)
                        # D∫ (∇ψ^i) ⋅ (∇ψ^j)dx
                        De[i, j] += μ * (dpsi_j ⋅ dpsi_i) * dx
                    end
                end
            end

            assemble!(assembler, cell_dofs, De)
        end

        return D
    end

    function doassemble_C!(
        C::SparseMatrixCSC, 
        cellvalues_rho::CellValues, 
        cellvalues_rhos::CellValues, 
        dh::DofHandler,
        velocity::Float64
        )
    
        n_basefuncs_rho = getnbasefunctions(cellvalues_rho)
        Ce = zeros(n_basefuncs_rho, n_basefuncs_rho)

        assembler = start_assemble(C)

        for cell in CellIterator(dh)
            
            # get global DOFs for this cell
            cell_dofs = celldofs(cell)
            ndofs_cell = length(cell_dofs)
            Ce = zeros(ndofs_cell, ndofs_cell)

            # rho is the only important bit since rho_s does not diffuse
            Ferrite.reinit!(cellvalues_rho, cell)
            fill!(Ce, 0)

            for q_point in 1:getnquadpoints(cellvalues_rho)
                dx = getdetJdV(cellvalues_rho, q_point)

                for i in 1:n_basefuncs_rho
                    psi_i = shape_value(cellvalues_rho, q_point, i)

                    for j in 1:n_basefuncs_rho
                        dpsi_j = shape_gradient(cellvalues_rho, q_point, j)
                        # -v \psi^i(\partial_x\psi^j)dx
                        Ce[i, j] += (velocity ⋅ dpsi_j) *  psi_i * dx
                    end
                end
            end
            assemble!(assembler, cell_dofs, Ce)
        end

        return C
    end

    function assemble_rhs_absorption!(
            N::Vector{Float64},
            u::Vector{Float64},
            params::RHSparams
        )

        N .= params.rhsL * u # This line computes the linear part of the right-hand side by multiplying the linear part of the right-hand side matrix (params.rhsL) with the solution vector (u). The result is stored in N, which represents the complete right-hand side vector for the system of equations.

        # now for the non-linear part:

        rhog_dofs = dof_range(params.dh, :rho_g)
        rhos_dofs = dof_range(params.dh, :rho_s)

        for cell in CellIterator(params.dh)

            cell_dofs = celldofs(cell)

            rhog_local = view(u, cell_dofs[rhog_dofs])
            rhos_local = view(u, cell_dofs[rhos_dofs])

            Ferrite.reinit!(params.cellvalues_rhog, cell)
            Ferrite.reinit!(params.cellvalues_rhos, cell)

            for q in 1:getnquadpoints(params.cellvalues_rhog)

                dx = getdetJdV(params.cellvalues_rhog,q)

                # calculate the values of rhog, rhos, v, and ∂v at the quadrature point
                rhog_q = sum(
                    rhog_local[j] *
                    shape_value(params.cellvalues_rhog,q,j)
                    for j in eachindex(rhog_local)
                )

                rhos_q = sum(
                    rhos_local[j] *
                    shape_value(params.cellvalues_rhos,q,j)
                    for j in eachindex(rhos_local)
                )

                # dot m = A*log(R*T*rhog/p_eqa)*(rho_sat-rhos)
                mdot_q = params.constants.MA*log(params.constants.R*params.constants.T*rhog_q/params.constants.p_eqa)*(params.constants.rho_sat-rhos_q)
                for i in 1:length(rhog_local)
                    ψ = shape_value(params.cellvalues_rhog,q,i)
                    N[cell_dofs[rhog_dofs[i]]] -= ψ*mdot_q*dx
                end

                for i in 1:length(rhos_local)
                    ψ = shape_value(params.cellvalues_rhos,q,i)
                    N[cell_dofs[rhos_dofs[i]]] += ψ*mdot_q*dx
                end
            end
        end
        return N
    end


    function solve_fem(
        cellvalues_rho_g,
        cellvalues_rho_s,
        ch,
        dh,
        constants,
        T_end,
        assemble_rhs!,
        t_saveat=2.0
        )
        # Assemble mass matrix
        M = allocate_matrix(dh)
        doassemble_M!(
            M,
            cellvalues_rho_g,
            cellvalues_rho_s,
            constants.ϵ,
            dh,
        )

        # Assemble linear RHS matrices
        D = allocate_matrix(dh)
        C = allocate_matrix(dh)

        doassemble_D!(
            D,
            cellvalues_rho_g,
            cellvalues_rho_s,
            constants.μ,
            dh,
        )

        doassemble_C!(
            C,
            cellvalues_rho_g,
            cellvalues_rho_s,
            dh,
            constants.velocity,
        )

        rhsL = -D - C

        # Construct parameters only after everything is ready
        params = RHSparams(
            constants,
            dh,
            cellvalues_rho_g,
            cellvalues_rho_s,
            M,
            zeros(Float64, ndofs(dh)),
            rhsL,
            assemble_rhs!,
        )


        # State vectors
        u      = zeros(ndofs(dh))  # solution vector
        du     = similar(u)        # time derivative of solution vector

        # Initial condition
        apply_analytical!(u, dh, :rho_g, x -> params.constants.p_eqa/(params.constants.R*params.constants.T))
        apply_analytical!(u, dh, :rho_s, x -> params.constants.ρ_e)

        # Initial BC application
        update!(ch, 0.0)
        apply!(u, ch);

        # residual function for DAE solver
        function res!(res, du, u, params, t)

            params.assemble_rhs!(params.rhs, u, params)

            mul!(res, params.M, du)      # res = M*du
            res .-= params.rhs

            update!(ch, t)
            for k in 1:length(ch.prescribed_dofs)
                dof = ch.prescribed_dofs[k]
                value = ch.inhomogeneities[k]
                res[dof] = u[dof] - value
            end
        end

        # enforce initial du
        update!(ch, 0.0)
        apply!(u, ch)          # enforce initial state

        dt = 1e-6

        params.assemble_rhs!(params.rhs, u, params)

        B = params.rhs * dt + (params.M * u)
        A = copy(params.M)
        apply!(A, B, ch)

        u_new = A \ B
        du = ( u_new - u) / dt

        # for k in eachindex(ch.prescribed_dofs)
        #     dof = ch.prescribed_dofs[k]
        #     println(
        #         "DOF ", dof,
        #         " u=", u[dof],
        #         " du=", du[dof]
        #     )
        # end

        differential_vars = trues(length(u))
        tspan = (0.0, T_end)


        prob = DAEProblem(
            res!,
            du,
            u,
            tspan,
            params;
            differential_vars=differential_vars
        )

        sol = solve(prob, saveat=t_saveat, reltol=1e-6, abstol=1e-6);

        return sol
    end

    # absorption and desorption solvers that call the generic solve_fem function with the appropriate assembly function for the right-hand side
    function solve_fem_absorption(
        cellvalues_rho_g,
        cellvalues_rho_s,
        ch,
        dh,
        constants,
        T_end,
        dt
        )

        return solve_fem(
            cellvalues_rho_g,
            cellvalues_rho_s,
            ch,
            dh,
            constants,
            T_end,
            assemble_rhs_absorption!,
            dt
        )
    end

    function solve_fem_desorption(
        cellvalues_rho_g,
        cellvalues_rho_s,
        ch,
        dh,
        constants,
        T_end
        )
        
        return solve_fem(
            cellvalues_rho_g,
            cellvalues_rho_s,
            ch,
            dh,
            constants,
            T_end,
            assemble_rhs_desorption!
        )
    end    


end