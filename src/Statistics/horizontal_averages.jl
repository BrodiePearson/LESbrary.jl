function velocity_variances(model; scratch = CellField(model.architecture, model.grid))
    u, v, w = model.velocities

    U² = HorizontalAverage(u^2, scratch)
    V² = HorizontalAverage(v^2, scratch)
    W² = HorizontalAverage(w^2, scratch)

    return U², V², W²
end

function horizontal_averages(model)

    # Create scratch space for calculations
    scratch = CellField(model.architecture, model.grid)

    # Extract short field names
    u, v, w = model.velocities
    U², V², W² = velocity_variances(model, scratch=scratch)
    # Extract Pressure (I'm not sure this is available in model)
    # p = model.pressure

    # Define horizontal averages
    U = HorizontalAverage(u)
    V = HorizontalAverage(v)
    e = TurbulentKineticEnergy(model)
    # When pressure is available calculate it's layer mean
    # P = HorizontalAverage(p)

    W³ = HorizontalAverage(w^3, scratch)
    wu = HorizontalAverage(w*u, scratch)
    wv = HorizontalAverage(w*v, scratch)
    uv = HorizontalAverage((u-U)*(v-V), scratch)
    wwu = HorizontalAverage(w*w*(u-U), scratch)
    wwv = HorizontalAverage(w*w*(v-V), scratch)
    wuu = HorizontalAverage(w*(u-U)*(u-U), scratch)
    wvv = HorizontalAverage(w*(v-V)*(v-V), scratch)
    wuv = HorizontalAverage(w*(u-U)*(v-V), scratch)
    
    # Diagnostics that involve pressure - need pressure passed
    #wp = HorizontalAverage(w*(p-P), scratch)
    #vp = HorizontalAverage((v-V)*(p-P), scratch)
    #up = HorizontalAverage((u-U)*(p-P), scratch)
    #pp = HorizontalAverage((p-P)*(p-P), scratch)
    #p∂u∂x = HorizontalAverage((p-P)*∂x(u), scratch)
    #p∂u∂y = HorizontalAverage((p-P)*∂y(u), scratch)
    #p∂v∂z = HorizontalAverage((p-P)*∂z(u), scratch)
    #p∂v∂x = HorizontalAverage((p-P)*∂x(v), scratch)
    #p∂v∂y = HorizontalAverage((p-P)*∂y(v), scratch)
    #p∂v∂z = HorizontalAverage((p-P)*∂z(v), scratch)
    #p∂w∂x = HorizontalAverage((p-P)*∂x(w), scratch)
    #p∂w∂y = HorizontalAverage((p-P)*∂y(w), scratch)
    #p∂w∂z = HorizontalAverage((p-P)*∂z(w), scratch)

    primitive_averages = (
                  U = model -> U(model),
                  V = model -> V(model),
                  E = model -> e(model),

                 U² = model -> U²(model),
                 V² = model -> V²(model),
                 W² = model -> W²(model),
                 W³ = model -> W³(model),

                 wu = model -> wu(model),
                 wv = model -> wv(model),
                 uv = model -> uv(model),

                 wwu = model -> wwu(model),
                 wwv = model -> wwv(model),
                 wuu = model -> wuu(model),
                 wvv = model -> wvv(model),
                 wuv = model -> wuv(model),

                 #wp = model -> wu(model),
                 #up = model -> wv(model),
                 #vp = model -> uv(model),
                 #pp = model -> wu(model),

                 # Syntax doesn't seem to work for the following
                 #p∂u∂x = model -> p∂u∂x(model),
                 #p∂u∂y = model -> p∂u∂y(model),
                 #p∂u∂z = model -> p∂u∂z(model),
                 #p∂v∂x = model -> p∂v∂x(model),
                 #p∂v∂y = model -> p∂v∂y(model),
                 #p∂v∂z = model -> p∂v∂z(model),
                 #p∂w∂x = model -> p∂w∂x(model),
                 #p∂w∂y = model -> p∂w∂y(model),
                 #p∂w∂z = model -> p∂w∂z(model),
               )

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

    return merge(primitive_averages, average_tracers, average_stresses, average_fluxes, average_diffusivities)
end
