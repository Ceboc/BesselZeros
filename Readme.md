# Besselzero.jl

This module implements the code as shown in the Boost C++ library to calculate zeros of Bessel functions. [Zeros of Bessel Functions (Boost)](https://www.boost.org/doc/libs/1_74_0/libs/math/doc/html/math_toolkit/bessel/bessel_root.html)
    
    boost/math/special_functions/detail/bessel_zero.hpp

It's also include the module Airyzero.jl that implements in Julia the equivalent functions to calculate zeros of the Airy functions. [Finding Zeros of Airy Functions (Boost)](https://www.boost.org/doc/libs/1_74_0/libs/math/doc/html/math_toolkit/airy/airy_root.html)
    
    boost/math/special_functions/detail/airy_ai_bi_zero.hpp

### Warning!
- This is still a work in progress, there are still problems with negative indices for bessel's functions.
- Comments still in spanish.  
