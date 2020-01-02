%     Testing givens rotations

clear
clc
close all


%A0 = [6 5 0 3;...
%      5 1 4 8;...
%      0 4 3 1;...
%      0 0 7 4;...
%      0 0 0 3];

%A0 = [2.1400 -9.4860 -1.1040  0.2779E+02;...
%      0.2088 -0.6699 -0.1638  0.4063E+01;...
%      0       25.180  0.1368 -0.8840E+02;...
%      0       0       0.5891  0.2340E+01;...
%      0       0       0       0.1791E+02];

A0 = [0.21404583E+01     -0.94858945E+01     -0.11036958E+01      0.27794887E+02;...
      0.20877078E+00     -0.66987236E+00     -0.16375599E+00      0.40632716E+01;...
      0.00000000E+00      0.25178930E+02      0.13681456E+02     -0.88403596E+02;...
      0.00000000E+00      0.00000000E+00      0.58909947E+00      0.23398410E+01;...
      0.00000000E+00      0.00000000E+00      0.00000000E+00      0.17906280E+02]


% Build Rotations for a Hessenberg Matrix

A1 = A0;          % Work with A1

[rw cl]=size(A1);


for i=1:cl

  r = A1(i:i+1,i);
  if (abs(A1(i+1,i))>0)
    rad = sqrt(r(1)^2 + r(2)^2);
    c = r(1)/rad;
    s = r(2)/rad;
  else
    c = 1;
    s = 0;
  end

  G          = eye(rw);
  G(i,i)     =  c;
  G(i,i+1)   =  s;
  G(i+1,i)   = -s;
  G(i+1,i+1) =  c;
  GR{i}      =  G;

  for j=i:cl
    r = A1(i:i+1,j);
    r2(1) =  c*r(1) + s*r(2);   
    r2(2) = -s*r(1) + c*r(2);
    A1(i:i+1,j)=r2';
  end

end

%A1

A2 = A0;
Q = eye(rw);
ngr = length(GR);
for i=1:ngr
  A2 = GR{i}*A2;
%  Q  = GR{ngr-i+1}'*Q;
   Q = Q*GR{i}';
end
A2;
%Q

%b0 = rand(rw,1);
b0 = zeros(rw,1);
b0(1) = 2.8166697378181341E-002;

b1 = b0;
for i=1:ngr
  b1 = GR{i}*b1;
end  
b1'

[u,s,v]=svd(A0);
soln1 = v*inv(s(1:ngr,1:ngr))*u(:,1:ngr)'*b0;
soln2 = inv(A1(1:ngr,1:ngr))*b1(1:ngr);

soln2'




