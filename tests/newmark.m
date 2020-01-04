%     Making a stupid test to see if newmark blows up


clear
clc
close all


M = 1.0;
C = 0.00;
K = 1.0;

omg = sqrt(K/M);
T   = 2*pi/omg;

a0 = 0;
v0 = 0;
x0 = 0;

xn = x0;
vn = v0;
an = a0;


dt = T/100;

nsteps = round(5*T/dt);

Fext = 10.0;
xhist = zeros(nsteps+1,1);
xhist(1) = xn;

for i=1:nsteps

%   if i>1
%     Fext = 0.;
%   end  

   f = 0;
   f = f + Fext;
   f = f + (4*M/(dt^2) + 2*C/dt)*xn;
   f = f + (4*M/dt     + C)*vn;
   f = f + M*an;

   A  = (4*M/(dt^2) + 2*C/dt + K);
   Ax = A*xn;

   xlag = xn;
%  Solve for increment   
   f = f - Ax;

   dx = f/A;

   xn = xn + dx;
   an = 4/(dt^2)*(xn-xlag) - 4/dt*vn - an; 
   vn = 2/dt*(xn - xlag) - vn;

   xhist(i+1) = xn;
end

plot(xhist)






