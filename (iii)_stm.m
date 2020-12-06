## Copyright (C) 2020 Abhinav
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} stm (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Abhinav <Abhinav@DESKTOP-SNQ7GLK>
## Created: 2020-12-01

function retval = stm (V,E)
    ## V: Barrier's Potential, E: Particle's Energy
    
    a = 1:0.001:7;  #Set of different barrier width values in A 
    a=a*(10^(-10)); #coverting to metres
    
    m = 9.10938*(10^(-31));
    hbar = 1.05457*(10^(-34));
    ec = 1.60218*(10^(-19));
    
    k1 = sqrt((2*m*E*ec)/(hbar^2));
    k2 = sqrt((2*m*(E-V)*ec)/(hbar^2));
    
    # Calculating Transmittance for each value of barrier width
    A = e.^(i.*(1));
    E = ((4*k1*k2*(e.^(i.*(k2-k1).*a)))./((k1+k2)^2 .- e.^(i*2*k2.*a).*((k2-k1)^2)))*A;
    E_conj = ((4*k1*k2*(e.^(-i.*(k2-k1).*a)))./((k1+k2)^2 .- e.^(-i*2*k2.*a).*((k2-k1)^2)))*A;
    T = E .* E_conj;
    T = abs(T);
    
    # Plotting Transmittance as a function of barrier width
    figure(1);
    plot(a,T);
    xlabel('Width of the barrier');
    ylabel('Transmittance');
    title('Transmittance variation with barrier width');
   
endfunction

