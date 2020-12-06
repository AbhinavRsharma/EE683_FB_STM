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
## @deftypefn {} {@var{retval} =} fbarrier2 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Abhinav <Abhinav@DESKTOP-SNQ7GLK>
## Created: 2020-11-27

function retval = fbarrier2 (V , E, a) 
    ## V: Barrier's Potential, E: Particle's Energy, a: barrier width
    
    a=a*(10^(-10)); # coverting to metres
    m = 9.10938*(10^(-31));
    hbar = 1.05457*(10^(-34));
    ec = 1.60218*(10^(-19));
    
    k1 = sqrt((2*m*E*ec)/(hbar^2));
    k2 = sqrt((2*m*(E-V)*ec)/(hbar^2));
    
    A = e.^(i.*(1));
    B = (((k1+k2)*(k2-k1)*(e.^(i*2*k2*a)-1))/((k1+k2)^2 - e.^(i*2*k2*a)*((k2-k1)^2)))*A;
    C = (((k1+k2)*2*k1)/((k1+k2)^2 - e.^(i*2*k2*a)*((k2-k1)^2)))*A;
    D = (((k2-k1)*2*k1*(e.^(i*2*k2*a)))/((k1+k2)^2 - e.^(i*2*k2*a)*((k2-k1)^2)))*A;
    E = ((4*k1*k2*(e.^(i*(k2-k1)*a)))/((k1+k2)^2 - e.^(i*2*k2*a)*((k2-k1)^2)))*A;
    
    # set of position values before the barrier
    x_left = -20:0.001:0; 
    x_left = x_left.*(10^(-10));
    # set of position values in the barrier
    x_barrier = 0.001:0.001:a*(10^(10)); 
    x_barrier = x_barrier.*(10^(-10));
    # set of position values after the barrier
    x_right = a*(10^(10))+0.001:0.001:20; 
    x_right = x_right.*(10^(-10));
    
    # Calculating wavefunctions values at each position value above
    y_left = A.*(e.^(i.*k1.*x_left)) .+  B.*(e.^(-i.*k1.*x_left));
    y_barrier = C.*(e.^(i.*k2.*x_barrier)) .+  D.*(e.^(-i.*k2.*x_barrier));
    y_right = E.*(e.^(i.*k1.*x_right));
    
    # Calculating real part of the wavefunctions
    y_left_r = real(y_left);
    y_barrier_r = real(y_barrier);
    y_right_r = real(y_right);
    
    # Calculating imaginary part of the wavefunctions
    y_left_i = imag(y_left);
    y_barrier_i = imag(y_barrier);
    y_right_i = imag(y_right);
    
    # Calculating probability density at each position value 
    y_left_p = y_left_r.^2+y_left_i.^2;
    y_barrier_p = y_barrier_r.^2+y_barrier_i.^2;
    y_right_p = y_right_r.^2+y_right_i.^2;
    
    # Plotting real part of the wavefunctions
    y_left_r_mx = max(y_left_r);
    y_barrier_r_mx = max(y_barrier_r);
    y_right_r_mx = max(y_right_r);
    y_r_mx = max(y_right_r_mx,y_left_r_mx);
    y_r_mx = max(y_r_mx,y_barrier_r_mx);
    figure(1);
    plot(x_left,y_left_r,'r');
    xlabel('Position');
    ylabel('Real part');
    title('Real part of the wavefunction');
    hold on;
    plot(x_barrier,y_barrier_r,'b');
    hold on;
    plot(x_right,y_right_r,'k');
    legend('left of the barrier(V=0ev)','in the barrier(V=5ev)','right of the barrier(V=0ev)');
    axis([-20*(10^(-10)) 20*(10^(-10)) -(y_r_mx+0.5) (y_r_mx+0.5)]);
    
    # Plotting probability density of the wavefunctions
    y_left_p_mx = max(y_left_p);
    y_barrier_p_mx = max(y_barrier_p);
    y_right_p_mx = max(y_right_p);
    y_p_mx = max(y_right_p_mx,y_left_p_mx);
    y_p_mx = max(y_p_mx,y_barrier_p_mx);
    figure(2);
    plot(x_left,y_left_p,'r');
    xlabel('Position');
    ylabel('Probability Density');
    title('Probability Density of the wavefunction');
    hold on;
    plot(x_barrier,y_barrier_p,'b');
    hold on;
    plot(x_right,y_right_p,'k');
    legend('left of the barrier(V=0ev)','in the barrier(V=5ev)','right of the barrier(V=0ev)');
    axis([-20*(10^(-10)) 20*(10^(-10)) 0 (y_p_mx+0.5)]);
    
endfunction

