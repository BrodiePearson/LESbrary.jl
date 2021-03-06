function pressure(model)
    p_hyd, p_non = model.pressures
    return p_hyd + p_non
end

# Will need to be revised...
buoyancy(model) = model.tracers.b

function horizontally_averaged_velocities(model)

    # Extract short field names
    u, v, w = model.velocities

    # Define horizontal averages
    U = HorizontalAverage(u)
    V = HorizontalAverage(v)

    averages = Dict(
                    :U => model -> U(model),
                    :V => model -> V(model)
                    )

    return averages
end

function horizontally_averaged_tracers(model)

    averages = Dict()

    for tracer in keys(model.tracers)
        c = getproperty(model.tracers, tracer)
        C = HorizontalAverage(c)
        averages[tracer] = model -> C(model)
    end

    return averages
end

function velocity_covariances(model, scratch = CellField(model.architecture, model.grid);
                                     b = buoyancy(model))

    u, v, w = model.velocities

    uu = HorizontalAverage(u * u, scratch)
    vv = HorizontalAverage(v * v, scratch)
    ww = HorizontalAverage(w * w, scratch)

    uv = HorizontalAverage(u * v, scratch)
    wv = HorizontalAverage(w * v, scratch)
    wu = HorizontalAverage(w * u, scratch)

    ub = HorizontalAverage(u * b, scratch)
    vb = HorizontalAverage(v * b, scratch)
    wb = HorizontalAverage(w * b, scratch)

    covariances = Dict(
                       :uv => model -> uv(model),
                       :vv => model -> vv(model),
                       :ww => model -> ww(model),
                       :wu => model -> wu(model),
                       :wv => model -> wv(model),
                       :uv => model -> uv(model),
                       :ub => model -> ub(model),
                       :vb => model -> vb(model),
                       :wb => model -> wb(model),
                       )

    return covariances
end

function tracer_covariances(model, scratch = CellField(model.architecture, model.grid),
                                   b = buoyancy(model))

    u, v, w = model.velocities

    covariances = Dict()

    for tracer in keys(model.tracers)
        c = getproperty(model.tracers, tracer)

        cc = HorizontalAverage(c * c, scratch)
        uc = HorizontalAverage(u * c, scratch)
        vc = HorizontalAverage(v * c, scratch)
        wc = HorizontalAverage(w * c, scratch)

        covariances[Symbol(tracer, tracer)] = model -> cc(model)
        covariances[Symbol(:u, tracer)] = model -> uc(model)
        covariances[Symbol(:v, tracer)] = model -> vc(model)
        covariances[Symbol(:w, tracer)] = model -> wc(model)

        # Add covariance of tracer with buoyancy
        if tracer != :b
            bc = HorizontalAverage(b * c, scratch)
            covariances[Symbol(:b, tracer)] = model -> bc(model)
        end
    end

    return covariances
end

"""
    third_order_velocity_statistics(model, scratch = CellField(model.architecture, model.grid))

Returns a dictionary of functions that calculate horizontally-averaged third-order statistics
that involve the velocity field.

Includes statistics associated with pressure.
"""
function third_order_velocity_statistics(model, scratch = CellField(model.architecture, model.grid))

    u, v, w = model.velocities

    uuu = HorizontalAverage(u * u * u, scratch)
    vvv = HorizontalAverage(v * v * v, scratch)
    www = HorizontalAverage(w * w * w, scratch)

    uuv = HorizontalAverage(u * u * v, scratch)
    uvv = HorizontalAverage(u * v * v, scratch)

    uuw = HorizontalAverage(u * u * w, scratch)
    uww = HorizontalAverage(u * w * w, scratch)

    vvw = HorizontalAverage(v * v * w, scratch)
    vww = HorizontalAverage(v * w * w, scratch)

    wvu = HorizontalAverage(w * v * u, scratch)

    # Pressure-strain terms
    p = pressure(model)

    Σˣʸ = (∂y(u) + ∂x(v)) / 2 # FFC
    Σˣᶻ = (∂z(u) + ∂x(w)) / 2 # FCF
    Σʸᶻ = (∂z(v) + ∂y(w)) / 2 # CFF

    up = HorizontalAverage(u * p, scratch)
    vp = HorizontalAverage(v * p, scratch)
    wp = HorizontalAverage(w * p, scratch)

    pΣˣˣ = HorizontalAverage(p * ∂x(u), scratch)
    pΣʸʸ = HorizontalAverage(p * ∂y(v), scratch)
    pΣᶻᶻ = HorizontalAverage(p * ∂z(w), scratch)

    pΣˣʸ = HorizontalAverage(p * Σˣʸ, scratch)
    pΣʸᶻ = HorizontalAverage(p * Σʸᶻ, scratch)
    pΣˣᶻ = HorizontalAverage(p * Σˣᶻ, scratch)

    third_order_statistics = Dict(
                                   :uuu => model -> uuu(model),
                                   :vvv => model -> vvv(model),
                                   :www => model -> www(model),
                                   :uuv => model -> uuv(model),
                                   :uvv => model -> uvv(model),
                                   :uuw => model -> uuw(model),
                                   :uww => model -> uww(model),
                                   :vvw => model -> vvw(model),
                                   :vww => model -> vww(model),
                                   :wvu => model -> wvu(model),

                                    :up => model -> up(model),
                                    :vp => model -> vp(model),
                                    :wp => model -> wp(model),

                                  :pΣˣˣ => model -> pΣˣˣ(model),
                                  :pΣʸʸ => model -> pΣʸʸ(model),
                                  :pΣᶻᶻ => model -> pΣᶻᶻ(model),

                                  :pΣˣʸ => model -> pΣˣʸ(model),
                                  :pΣʸᶻ => model -> pΣʸᶻ(model),
                                  :pΣˣᶻ => model -> pΣˣᶻ(model),
                                  )

    return third_order_statistics
end

function third_order_tracer_statistics(model, scratch = CellField(model.architecture, model.grid))

    u, v, w = model.velocities

    # Pressure-strain terms
    p = pressure(model)

    third_order_statistics = Dict()

    for tracer in keys(model.tracers)
        c = getproperty(model.tracers, tracer)

        cwu = HorizontalAverage(c * w * u, scratch)
        wcc = HorizontalAverage(w * c * c, scratch)

        cpx = HorizontalAverage(c * ∂x(p), scratch)
        cpy = HorizontalAverage(c * ∂y(p), scratch)
        cpz = HorizontalAverage(c * ∂z(p), scratch)

        third_order_statistics[Symbol(:w, tracer, tracer)] = model -> wcc(model)
        third_order_statistics[Symbol(tracer, :wu)] = model -> cwu(model)

        third_order_statistics[Symbol(tracer, :px)] = model -> cpx(model)
        third_order_statistics[Symbol(tracer, :py)] = model -> cpy(model)
        third_order_statistics[Symbol(tracer, :pz)] = model -> cpz(model)
    end

    return third_order_statistics
end

function subfilter_viscous_dissipation(model)

    u, v, w = model.velocities

    νₑ = model.diffusivities.νₑ

    Σˣˣ = ∂x(u)
    Σʸʸ = ∂y(v)
    Σᶻᶻ = ∂z(w)
    Σˣʸ = (∂y(u) + ∂x(v)) / 2
    Σˣᶻ = (∂z(u) + ∂x(w)) / 2
    Σʸᶻ = (∂z(v) + ∂y(w)) / 2

    ϵ = 2 * νₑ * ( Σˣˣ^2 + Σʸʸ^2 + Σᶻᶻ^2 + Σˣʸ^2 + Σˣᶻ^2 + Σʸᶻ^2 )

    return ϵ
end

function first_order_statistics(model, scratch = CellField(model.architecture, model.grid))

    output = merge(
                   horizontally_averaged_velocities(model),
                   horizontally_averaged_tracers(model),
                   )

    p = pressure(model)
    P = HorizontalAverage(p, scratch)

    output[:P] = model -> P(model)

    return output
end

function second_order_statistics(model, scratch = CellField(model.architecture, model.grid))

    output = merge(
                   velocity_covariances(model, scratch),
                   tracer_covariances(model, scratch),
                   )

    return output
end

function third_order_statistics(model, scratch = CellField(model.architecture, model.grid))

    output = merge(
                   third_order_velocity_statistics(model, scratch),
                   third_order_tracer_statistics(model, scratch),
                   )

    ϵ = HorizontalAverage(subfilter_viscous_dissipation(model), scratch)
    output[:ϵ] = model -> ϵ(model)

    return output
end

function first_through_third_order(model, scratch = CellField(model.architecture, model.grid))

    output = merge(
                   first_order_statistics(model, scratch),
                   second_order_statistics(model, scratch),
                   third_order_statistics(model, scratch),
                   )

    return output
end

function horizontal_averages(model)

    # Create scratch space for calculations
    scratch = CellField(model.architecture, model.grid)

    # Extract short field names
    u, v, w = model.velocities

    # Define horizontal averages
    U = HorizontalAverage(u)
    V = HorizontalAverage(v)
    e = TurbulentKineticEnergy(model)

    W³ = HorizontalAverage(w^3, scratch)
    wu = HorizontalAverage(w*u, scratch)
    wv = HorizontalAverage(w*v, scratch)

    primitive_averages = (
                           U = model -> U(model),
                           V = model -> V(model),
                           E = model -> e(model),
                          W³ = model -> W³(model),
                          wu = model -> wu(model),
                          wv = model -> wv(model),
                         )

    vel_variances = velocity_variances(model; scratch=scratch)

    # Add subfilter stresses (if they exist)
    average_stresses = Dict()

    νₑ = model.diffusivities.νₑ

    NU = HorizontalAverage(νₑ)

    τ₁₃ = @at (Face, Cell, Face) νₑ * (-∂z(u) - ∂x(w))
    τ₂₃ = @at (Cell, Face, Face) νₑ * (-∂z(v) - ∂y(w))
    τ₃₃ = @at (Cell, Cell, Face) νₑ * (-∂z(w) - ∂z(w))

    T₁₃ = HorizontalAverage(τ₁₃, scratch)
    T₂₃ = HorizontalAverage(τ₂₃, scratch)
    T₃₃ = HorizontalAverage(τ₃₃, scratch)

    average_stresses[:τ₁₃] = model -> T₁₃(model)
    average_stresses[:τ₂₃] = model -> T₂₃(model)
    average_stresses[:τ₃₃] = model -> T₃₃(model)
    average_stresses[:νₑ] = model -> NU(model)

    average_stresses = (; zip(keys(average_stresses), values(average_stresses))...)

    average_tracers = Dict()
    average_fluxes = Dict()
    average_diffusivities = Dict()

    for tracer in keys(model.tracers)

        advective_flux_key = Symbol(:w, tracer)
        subfilter_flux_key = Symbol(:q₃_, tracer)
           diffusivity_key = Symbol(:κₑ_, tracer)

        w = model.velocities.w
        c = getproperty(model.tracers, tracer)

        # Average tracer
        average_tracer = HorizontalAverage(c)
        average_tracers[tracer] = model -> average_tracer(model)

        # Advective flux
        advective_flux = w * c
        average_advective_flux = HorizontalAverage(advective_flux, scratch)
        average_fluxes[advective_flux_key] = model -> average_advective_flux(model)

        # Subfilter diffusivity (if it exists)
        try
            κₑ = getproperty(model.diffusivities.κₑ, tracer)
            average_diffusivity = HorizontalAverage(κₑ)
            average_diffusivities[diffusivity_key] = model -> average_diffusivity(model)

            subfilter_flux = @at (Cell, Cell, Face) -∂z(c) * κₑ
            average_subfilter_flux = HorizontalAverage(subfilter_flux, scratch)
            average_fluxes[subfilter_flux_key] = model -> average_subfilter_flux(model)
        catch
        end
    end

    average_tracers = (; zip(keys(average_tracers), values(average_tracers))...)
    average_fluxes = (; zip(keys(average_fluxes), values(average_fluxes))...)
    average_diffusivities = (; zip(keys(average_diffusivities), values(average_diffusivities))...)

    return merge(primitive_averages,
                 vel_variances,
                 average_tracers,
                 average_stresses,
                 average_fluxes,
                 average_diffusivities)
end
