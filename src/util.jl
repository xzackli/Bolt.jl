
# utilities for x ↔ scale factor ↔ redshift
a2z(a::T) where T = one(T)/a - one(T)
z2a(z::T) where T = one(T)/(z + one(T))
a2x(a) = log(a)
x2a(x) = exp(x)
z2x(z) = a2x(z2a(z))
x2z(x) = a2z(x2a(x))

# utility function for constructing an interpolator
spline(x, y) = scale(interpolate(y, BSpline(Cubic(Line(OnGrid())))), x)
spline_∂ₓ(f, x_grid) = spline(x_grid, [Interpolations.gradient(f, x)[1] for x in x_grid])
spline_∂ₓ²(f, x_grid) = spline(x_grid, [Interpolations.hessian(f, x)[1] for x in x_grid])

δ_kron(i, j) = (i == j) ? 1 : 0



"""

    @⌛ code ...
    @⌛ function_definition() = .... 

Label a section of code to be timed. The first form uses the code
itselfs as a label, the second uses the function name, and its the
body of the function which is timed. 

To run the timer and print output, returning the result of the
calculation, use

    @show⌛ run_code()

Timing uses `TimerOutputs.get_defaulttimer()`. 
"""
macro ⌛(ex)
    source_str = last(splitpath(string(__source__.file)))*":"*string(__source__.line)
    try
        # function definition
        sdef = splitdef(ex)
        sdef[:body] = quote
            $TimerOutputs.@timeit $("$(string(sdef[:name]))(…)  ($source_str)") $(sdef[:body])
        end
        esc(combinedef(sdef))
    catch
        # anything else
        :(@timeit $("$(Base._truncate_at_width_or_chars(string(prewalk(rmlines,ex)),26))  ($source_str)") $(esc(ex)))
    end
end


"""
See [`@⌛`](@ref)
"""
macro show⌛(ex)
    quote
        reset_timer!($TimerOutputs.get_defaulttimer())
        result = $(esc(ex))
        show($TimerOutputs.get_defaulttimer())
        result
    end
end
