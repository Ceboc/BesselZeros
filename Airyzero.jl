module Airyzero

#Funciones que regresa los ceros de las funciones de Airy
##REPLICO EL CÓDIGO DE LIBRERÍA BOOST PARA C++
##VER boost/math/special_functions/detail/airy_ai_bi_zero.hpp

#using Cxx;
export airyzero
export airyzerobi
#cxx""" 
    #include <iostream>
    #include <boost/math/special_functions/airy.hpp>
    #include <boost/multiprecision/cpp_dec_float.hpp>
#    typedef boost::multiprecision::cpp_dec_float_50 float_type;
#    """

#cxx"""
#    double airyzero1(int y) {
#    return boost::math::airy_ai_zero<double>(y);
#    }
#    """

#    airyzero(ind) = @cxx airyzero1(ind)

import Roots.find_zero
import Roots.Newton
import SpecialFunctions.airyai
import SpecialFunctions.airyaiprime
import MLStyle: @match

function equation_as_10_4_105(z::Real)
    one_over_z = 1/z
    one_over_z_squared = one_over_z ^ 2

    z_pow_third = cbrt(z)
    z_pow_two_thirds = z_pow_third ^2

    fz = z_pow_two_thirds * (((((                          (162375596875 / 334430208 )
                                    * one_over_z_squared - (   108056875 /   6967296 ))
                                    * one_over_z_squared + (       77125 /     82944 ))
                                    * one_over_z_squared - (           5 /        36 ))
                                    * one_over_z_squared + (           5 /        48 ))
                                    * one_over_z_squared + 1)
    return fz 
end

function initial_guess_ai(m::Integer)
    if 0≤m≤10
    guess = @match m begin
        0  => 0
        1  => -2.33810741045976703849
        2  => -4.08794944413097061664
        3  => -5.52055982809555105913
        4  => -6.78670809007175899878
        5  => -7.94413358712085312314
        6  => -9.02265085334098038016
        7  => -10.0401743415580859306
        8  => -11.0085243037332628932
        9  => -11.9360155632362625170
        10 => -12.8287767528657572004
        end
    else
         t = (π*3)*(m*4-1)/8
         guess = -equation_as_10_4_105(t)
    end

    return guess
end

function initial_guess_bi(m::Integer)
    if 0≤m≤10
    guess = @match m begin
        0  => 0
        1  => -1.17371322270912792492
        2  => -3.27109330283635271568
        3  => -4.83073784166201593267
        4  => -6.16985212831025125983
        5  => -7.37676207936776371360
        6  => -8.49194884650938801345
        7  => -9.53819437934623888663
        8  => -10.5299135067053579244
        9  => -11.4769535512787794379
        10 => -12.3864171385827387456
        end
    else
         t = (π*3)*(m*4-3)/8
         guess = -equation_as_10_4_105(t)
    end

    return guess
end

function airyzero(n::Integer)
    an = initial_guess_ai(n)
    return find_zero((airyai,airyaiprime),an,Newton())
end

function airyzerobi(n::Integer)
    an = initial_guess_bi(n)
    return find_zero((airybi,airybiprime),an,Newton())
end


end