clear
clc
close all

fnamef = 'fluid/fluid0.f00001';
fnames = 'struc/corstruc0.f00001';

dsf = readnek(fnamef);
dss = readnek(fnames);

%fe = [6 12 25 26 35 29];

fe = [6 12];

for i=fe
  xf = dsf.flddata(1).data(:,:,:,i);
  yf = dsf.flddata(2).data(:,:,:,i);
  uf = dsf.flddata(3).data(:,:,:,i);

  surf(xf,yf,uf,'EdgeColor','interp'); hold on
end  

se = [];

for i=se
  xs = dss.flddata(1).data(:,:,:,i);
  ys = dss.flddata(2).data(:,:,:,i);
  us = dss.flddata(3).data(:,:,:,i);

  surf(xs,ys,us,'EdgeColor','interp'); hold on
end  

colorbar


fnamesg = 'struc/glostruc0.f00001';

dssg = readnek(fnamesg);

i=1;
v = dsf.flddata(3).data(:,:,:,i)'

i=7;
v = dsf.flddata(3).data(:,:,:,i)'















