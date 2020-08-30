module BesselsZero

export besseljzero, besselyzero
import Roots

#Raices por intervalos https://zenodo.org/badge/latestdoi/87007945 por Dr. D. Sanders y Dr. L. Bennet
import IntervalArithmetic
import IntervalRootFinding
import SpecialFunctions: bessely, besselj
include("Airyzero.jl")

##REPLICO EL CÓDIGO DE LIBRERÍA BOOST PARA C++
##VER boost/math/special_functions/detail/bessel_zero
function equation_nist_10_21_19(ν::Real,a::Real)

    # Obtiene el estimado inicial de la m-ésima raíz de Jν o Yν.
    # Esta función se usa para órdenes m > 1.
    # El orden m ha sido usado para crear al parámetro a. 

    # Esta es la ecuación 10.21.19 de NIST Handbook. http://dlmf.nist.gov/10.21.E19
    μ = ν^2 * 4
    μ_minus_one = μ-1
    eigth_a_inv = 1/(a*8)
    eight_a_inv_squared = eigth_a_inv^2
    
    term3 = ((μ_minus_one * 4)  *     ((μ *    7 ) - 31 ))/3
    term5 = ((μ_minus_one * 32) *   ((((μ *   83 ) - 982  )  * μ ) +3779 ))/15
    term7 = ((μ_minus_one * 64) * ((((((μ * 6949 ) -153855 ) * μ ) +1585743)*μ)-6277237))/105

    return a + ((((                       - term7
                    * eight_a_inv_squared - term5)
                    * eight_a_inv_squared - term3)
                    * eight_a_inv_squared - μ_minus_one)
                    * eigth_a_inv);
end

function equation_as_9_3_39(z::Real,ζ::Real)
    # Regresa la función de ζ definida implicitamente en
    # A&S Ec. 9.3.39 como función z. 

    zsq_minus_one_sqrt = sqrt((z*z)-1)

    the_function = zsq_minus_one_sqrt - (acos(1/z)+((2/3)*(ζ*sqrt(ζ))))

    return the_function
end

function equation_as_9_3_39_derivative(z::Real)
    # Regresa la función de ζ definida implicitamente en
    # A&S Ec. 9.3.39 como función z. La función se Regresa
    # junto con el valor de su derivada respecto de z.  

    zsq_minus_one_sqrt = sqrt((z*z)-1)

    its_derivative = zsq_minus_one_sqrt/z

    return its_derivative
end

function equation_as_9_5_26(ν::Real,ai_bi_root::Real)

    # Obtiene el estimado el m-ésimo cero de Jν o Yν.
    # El orden m ha sido usado apra crear el parámetro de entrada ai_bi_root.
    # Aquí , ν > 2.2. El estimado se calcula de Abramowitz and Stegun 
    # eqs. 9.5.22 y 9.5.26, página 371.
    
    # La inversión de z como función de ζ se menciona en el texto  A&S siguiendo
    # eq. 9.5.26. Aquí se logra la inversión al realizar una expansión de Taylor de
    # eq. 9.3.39  para z grande a orden 2 y resolviendo la resultante ecuación cuadrática,
    # luego, se toma la raíz positiva de la solución cuadrática.
    # Es decir: (2/3)(-ζ)^(3/2) ≈ z + 1/(2z) - π/2 .
    # Que lleva a : z^2 - [(2/3)(-ζ)^(3/2) + π/2 ] z + 1/2 = 0

    # Con este estimado inicial, se hace una iteración de Newton-Raphson para refinar 
    # el valor de la raiz de z como función de ζ

    v_pow_third = cbrt(ν)
    v_pow_minus_two_thirds = 1/(v_pow_third^2)
    
    # Se obtene ζ usando el orden ν combinado con la m-ésima raíz de
    # una función de airy, como se muestra en A&S eq. 9.5.22
    ζ = v_pow_minus_two_thirds * (- ai_bi_root)

    ζ_sqrt = sqrt(ζ)

    # Se establece una ecuación cuadrática basada en la serie de Taylor
    # que se mencionó arriba
    b = -((((ζ * ζ_sqrt) * 2)/3) + π/2)

    z_estimate = (-b + sqrt(b^2-2))/2

    # Se calcula la raíz de z como función de ζ
    z = Roots.find_zero((x->equation_as_9_3_39(x,ζ),equation_as_9_3_39_derivative),z_estimate,Roots.Newton())

    #Continúa con la implementación de A&s ec. 9.3.39
    zsq_minus_one = z^2 - 1
    zsq_minus_one_sqrt = sqrt(zsq_minus_one)

    # A&S ec. 9.3.42
    b0_term_5_24 = 5 / ((zsq_minus_one * zsq_minus_one)*24)
    b0_term_1_8  = 1 / ( zsq_minus_one_sqrt * 8)
    b0_term_5_48 = 5 / ( ζ^2 * 48)
    
    b0 = -b0_term_5_48 + ((b0_term_5_24+b0_term_1_8)/ζ_sqrt)

    # Segunda línea de A&s ec. 9.5.26 para fₖ con k = 1
    f1 = ((z * ζ_sqrt) * b0) / zsq_minus_one_sqrt

    # A&S 9.5.22 expandida para k=1 (i.e. hasta el primer término de la serie)
    return (ν * z) + (f1 / ν)
end

function equation_nist_10_21_40_a(ν::Real)
    v_pow_third = cbrt(ν)
    v_pow_minus_two_thirds = 1 / v_pow_third^2
    return ν * (((((                          + 0.043
                    * v_pow_minus_two_thirds  - 0.0908)
                    * v_pow_minus_two_thirds  - 0.00397)
                    * v_pow_minus_two_thirds  + 1.033150)
                    * v_pow_minus_two_thirds  + 1.8557571)
                    * v_pow_minus_two_thirds  + 1)
end

function jν(ν::Real,x::Real)
    # Obtiene Jν(x)
    # Si ν>154 usa número de presición doble.
    # Julia usa la definición de 10.2.2 del NIST
    if ν>151
        return besselj(ν,big(x))
    else
        return besselj(ν,x)
    end 
end

function jν_prime(ν::Real,x::Real)
    #Obtiene Jν'(x), ver 10.6.2 y 10.6.3 del NIST
    if ν==0
        return -jν(1,x)
    else
        j_ν = jν(ν,x)
        j_ν_m1 = jν(ν-1,x)
        return j_ν_m1-(ν * j_ν)/x
    end
end

function jν_prime2(ν::Real,x::Real)
    if ν == 0
        return jν(1,x)/x - jν(0,x)
    else
        inv_x = 1/x
        ν_m1 = ν-1
        return jν(ν_m1-1,x) + inv_x * (-(ν_m1+ν)*jν(ν_m1,x)
                +inv_x*(ν_m1)*ν*jν(ν,x))  
    end
end

function initial_guess_jν(ν,m)
    #Manejo especial para órdenes negativos
    if ν<0
        if (m == 1) && (ν > 0.5)
            # Para ν pequeñas y negativas, se usan los resultados del
            # ajuste emírico de curvas. Sesión de Mathematica(R) para los coeficientes
            # Table[{n, BesselJZero[n, 1]}, {n, -(1/2), 0, 1/10}];
            # N[%, 20];
            # Fit[%, {n^0, n^1, n^2, n^3, n^4, n^5, n^6}, n] 
            guess = ((((((   -0.2321156900728593063
                          *ν -0.1493247777488004455)
                          *ν -0.1520541916723949868)
                          *ν +0.0781493056124887375)
                          *ν -0.1775757353768819240)
                          *ν +1.542805677045662954)
                          *ν +2.404825557695772769)
            return guess
        end 

        # Se crea el orden positivo y se extrae su parte entera positivo
        νν = -ν
        νν_floor = floor(νν)
        
        # La raiz a ser encontrada está acotada por las raices de la función de bessel 
        # a quien reflejan, el orden enter positivo es menor que, pero más cercano a νν.
        # Ver NIST sec. 10.21(i) https://dlmf.nist.gov/10.2.1#i
        root_hi = initial_guess_jν(νν_floor,m)

        if m == 1
            #El estimado de la primera raiz para orden negativo se encuentra usando un
            #algoritmo de busqueda de rango adaprativo.
            root_lo = root_hi - 0.1
            hi_end_if_bracket_is_negative = jν(ν,root_hi)<0
            while root_lo > eps()
                lo_end_of_bracket_is_negative = jν(ν,root_lo)<0
                if hi_end_if_bracket_is_negative != lo_end_of_bracket_is_negative
                    break
                end

                root_hi = root_lo

                if root_lo > 0.5
                    root_lo -= 0.5
                else
                    root_lo *= 0.75
                end
                    
            end
        else
            root_lo = initial_guess_jν(νν_floor,m-1)
        end
        guess_range = IntervalArithmetic.@interval(root_lo,root_hi)
        guess_pair = IntervalRootFinding.roots(z_eval ->jν(ν,z_eval),guess_range,IntervalRootFinding.Bisection)
        guess_pair = IntervalArithmetic.interval.(guess_pair)[1]

        return IntervalArithmetic.mid(guess_pair)
    end

    if m == 1
        if ν < 2.2
            # Para ν pequeña, se usan los resultados de ajuste émpirico de curva.
            # Sessión de Mathematica(R) para los coeficientes
            #Table[{n, BesselJZero[n, 1]}, {n, 0, 22/10, 1/10}];
            # N[%, 20];
            # Fit[%, {n^0, n^1, n^2, n^3, n^4, n^5, n^6}, n]
            guess = ((((((   - 0.000834237904601039788
                          *ν - 0.00759003563741012167)
                          *ν - 0.03064091477201269534)
                          *ν + 0.0782320880201056074)
                          *ν - 0.1696687125906200374)
                          *ν + 1.542187960073749658)
                          *ν + 2.404835991525463429)
        else
            # Para ν ≥ 2.2, se usa la primera línea de las ecs. 10.21.40 del NIST Handbook
            guess = equation_nist_10_21_40_a(ν)
        end
    else
        if ν < 2.2
            # Se usan ec. 10.21.19 del NIST Handbook
            a = (ν+(m*2)-0.5)*π/2

            guess = equation_nist_10_21_19(ν,a)
        else
            # Se obtiene el estimado de la m-ésima raiz de airyai
            airy_ai_root = Airyzero.airyzero(m)

            guess = equation_as_9_5_26(ν,airy_ai_root)
        end
    end

    return guess
end

function equation_nist_10_21_40_b(ν::Real)
    v_pow_third = cbrt(ν)
    v_pow_minus_two_thirds = 1/(v_pow_third^2)

    return ν * (((((                         - 0.001  
                    * v_pow_minus_two_thirds - 0.0060)
                    * v_pow_minus_two_thirds + 0.01198)
                    * v_pow_minus_two_thirds + 0.260351)
                    * v_pow_minus_two_thirds + 0.9315768)
                    * v_pow_minus_two_thirds + 1)
end

function yν(ν::Real,x::Real)
    # Obtiene Jν(x)
    # Si ν>151 usa número de presición doble.
    # Julia usa la definición de 10.2.3 del NIST
    if ν>151
        return bessely(ν,big(x))
    else
        return bessely(ν,x)
    end 
end

function yν_prime(ν::Real,x::Real)
    #Obtiene Yν'(x), ver 10.6.2 y 10.6.3 del NIST
    if ν==0
        return -yν(1,x)
    else
        y_ν = yν(ν,x)
        y_ν_m1 = yν(ν-1,x)
        return y_ν_m1-(ν * y_ν)/x
    end
end

function yν_prime2(ν::Real,x::Real)
    if ν == 0
        return yν(1,x)/x - yν(0,x)
    else
        inv_x = 1/x
        ν_m1 = ν-1
        return yν(ν_m1-1,x) + inv_x * (-(ν_m1+ν)*yν(ν_m1,x)
                +inv_x*(ν_m1)*ν*yν(ν,x))  
    end
end

function initial_guess_yν(ν::Real,m::Integer)
    # Calcula un estimado de la m-ésima raiz de Yν
    if ν<0
        # Crea el orden positivo y extrae las partes enteras de su suelo y techo positivo
        νν = -ν
        νν_floor = floor(νν)
        # La raiz a ser encontrada está acotada por las raices de la función de bessel 
        # a quien reflejan, el orden enter positivo es menor que, pero más cercano a νν.
        # Ver NIST sec. 10.21(i) https://dlmf.nist.gov/10.2.1#i

        # El caso especial de un orden semientero y negativo, se usa la relación
        # entra Yν y las funciones de Bessel esféricas, a fin de obtener la 
        # raiz acorada.
        # En estos casos especiales, Yν(-n/2,x) = sph_bessel_j(n/2,x)
        # para ν = -n/2

        if m==1
            if (νν - νν_floor)< 0.5
                root_hi = initial_guess_yν(νν_floor,m)
            else
                root_hi = initial_guess_yν(νν_floor+0.5,m)
            end

            root_lo = root_hi - 0.1

            hi_end_if_bracket_is_negative = yν(ν,root_hi) < 0

            while root_lo > eps()
                lo_end_of_bracket_is_negative = yν(ν,root_lo) < 0
                if hi_end_if_bracket_is_negative != lo_end_of_bracket_is_negative
                    break
                end

                root_hi = root_lo

                # Se reduce la cóna inferior usando un algoritmo adaptitivo
                if root_lo > 0.5
                    root_lo -= 0.5
                else
                    root_lo *= 0.75
                end
            end
        else
            if (νν - νν_floor) < 0.5
                root_lo = initial_guess_yν(νν_floor,m-1)
                root_hi = initial_guess_yν(νν_floor,m)
                root_lo += 0.01
                root_hi += 0.01
            else
                root_lo = initial_guess_yν(νν_floor+0.5,m-1)
                root_hi = initial_guess_yν(νν_floor+0.5,m)
                root_lo += 0.01
                root_hi += 0.01
            end
        end

        guess_range = IntervalArithmetic.@interval(root_lo,root_hi)
        @show guess_range
        guess_pair = IntervalRootFinding.roots(z_eval ->jν(ν,z_eval),guess_range,IntervalRootFinding.Bisection,1e-4)
        guess_pair = IntervalArithmetic.interval.(guess_pair)[1]

        return IntervalArithmetic.mid(guess_pair)
    end

    if m == unsigned(1)
        # Obtiene el estimado dela primera raiz
        if ν < 2.2
            # Para ν pequeña, se usan los resultados de ajuste empírico de curva.
            # Sessión de Mathematica(R) para los coeficientes
            # Table[{n, BesselYZero[n, 1]}, {n, 0, 22/10, 1/10}];
            # N[%, 20];
            # Fit[%, {n^0, n^1, n^2, n^3, n^4, n^5, n^6}, n] 
            guess = ((((((   - 0.002509590923565226250
                          *ν + 0.02129188704905315827)
                          *ν - 0.0764877854865257993)
                          *ν + 0.1591102681153622068)
                          *ν - 0.2416816687651963761)
                          *ν + 1.443784631088524401)
                          *ν + 0.893621151902004904)
        else
            # Para ν>2.2, se usa la segunda línea de las ecs. 10.21.40 en el NIST Handbook
            guess = equation_nist_10_21_40_b(ν)
        end
    else
        if ν < 2.2
            # Se usa la ec. 10.21.19 en el NIST Handbook
            a = (ν+(m*2)-1.5)*π/2
            guess = equation_nist_10_21_19(ν,a)
        else
            # Obtiene un estimado de la m-ésima raiz de airybi
            airy_bi_root = Airyzero.initial_guess_bi(m)

            # Usa la ec. 9.5.26 de A&S Handbook
            guess = equation_as_9_5_26(ν,airy_bi_root)
        end
    end
    return guess
end

#using ForwardDiff
function besselyzero(ν::Real,m::Integer)
    guess = initial_guess_yν.(ν, m)
    return Roots.find_zero((x->yν(ν,x),x->yν_prime(ν,x),x->yν_prime2(ν,x)),guess,Roots.Halley(),maxevals=2000,verbose=true)
end

function besseljzero(ν::Real,m::Integer)
    guess = initial_guess_jν.(ν, m)
    return Roots.find_zero((x->jν(ν,x),x->jν_prime(ν,x),x->jν_prime2(ν,x)),guess,Roots.Halley(),maxevals=2000,verbose=true)
end

end