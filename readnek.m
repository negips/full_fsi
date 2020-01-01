function [datastruc,data,lr1,nelt,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(fname)
%
% This function reads binary data from the nek5000 file format
%
%   [data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(fname)
%
%   INPUT
%   - fname:  name of the file 
%
%   OUTPUT
%   - data:   nek5000 data ordered as (iel,inode,[x|y|(z)|u|v|(w)|p|T|s_i])
%   - lr1:    element-size vector (lx1,ly1,lz1)
%   - elmap:  reading/writing map of the elements in the file
%   - time:   simulation time
%   - istep:  simulation step
%   - fields: fields saved in the file
%   - emode:  endian mode 'le' = little-endian, 'be' = big-endian
%   - wdsz:   single (4) or double (8) precision
%   - etag:   tag for endian indentification
%   - header: header of the file (string)
%   - status: status (< 0 something went wrong)
%
%
% Edit history: 20151028 Nicolo Fabbiane (nicolo@mech.kth.se)
%             : 20160707 Prabal Negi (negi@mech.kth.se): Added scalar reading and structure output.

%--------------------------------------------------------------------------
%  INITIALIZE OUTPUT
%--------------------------------------------------------------------------
data   = [];
etag   = [];
lr1    = [];
elmap  = [];
time   = [];
istep  = [];
fields = [];
wdsz   = [];
header = [];

disp(['Reading file: ', fname])

%--------------------------------------------------------------------------
%  OPEN THE FILE
%--------------------------------------------------------------------------
emode = 'le';
[infile,message] = fopen(fname,'r',['ieee-' emode]);
if infile == -1, disp(message), status = -1; return, end
%
% read header
header = fread(infile,132,'*char')';
%
% check endian
etag = round(fread(infile,1,'*float32')*1e5);
if (etag ~= 654321)
    fclose(infile);
    
    emode = 'be';
    [infile,message] = fopen(fname,'r',['ieee-' emode]);
    
    header = fread(infile,132,'*char')';
    etag = round(fread(infile,1,'*float32')*1e5);
    
    if (etag ~= 654321)
        disp('ERROR: could not interpret endianness.')
        status = -3; return
    end
end
etag = etag * 1e-5;

%--------------------------------------------------------------------------
% READ HEADER
%--------------------------------------------------------------------------
%
% word size
wdsz = str2double(header(6));
if (wdsz == 4)
    realtype = '*float32';
elseif (wdsz == 8)
    realtype = '*float64';
else
    fprintf('ERROR: could not interpret real type (wdsz = %i)',wdsz);
    status = -2; return
end
%
% element size
lr1 = [str2double(header(8:9))
       str2double(header(11:12))
       str2double(header(14:15))];
%
% compute the total number of points per element
npel = prod(lr1);
%
% compute number of active dimensions
ndim = 2 + (lr1(3)>1);
%
% number of elements
nel = str2double(header(17:26));
%
% number of elements in the file
nelt = str2double(header(28:37));
%
% time
time = str2double(header(39:58));
%
% istep
istep = str2double(header(60:68));
%
% get file id  
fid = str2double(header(70:75)); % TODO: multiple files not supported
%
% get tot number of files
nf = str2double(header(77:82)); % TODO: multiple files not supported
%
% getfields [XUPTS]
fields = strtrim(header(84:end));
var=zeros(1,5);
if sum(fields == 'X') > 0
  var(1) = ndim;
end
if sum(fields == 'U') > 0
  var(2) = ndim;
end
if sum(fields == 'P') > 0
  var(3) = 1;
end
if sum(fields == 'T') > 0
  var(4) = 1;
end
if sum(fields == 'S') > 0
      var(5) = str2double(fields(end-1:end)); % TODO: scalars not implemented
end
nfields = sum(var);
%
% read element map
elmap = fread(infile,nelt,'*int32').';

%--------------------------------------------------------------------------
% READ DATA
%--------------------------------------------------------------------------
dtemp = zeros(npel,nelt);
data = zeros(npel,nelt,nfields);
d1.data=dtemp;
datastruc = repmat(d1,1,nfields);

% For X,U           ! stored as vectors
for ivar = 1:2
    idim0 = sum(var(1:ivar-1));
    for iel = elmap
        for idim = (1:var(ivar))+idim0
            dtemp = fread(infile,npel,realtype);  
            data(:,iel,idim) = dtemp;
            datastruc(idim).data(:,iel)=dtemp;
        end
    end
end

% For P,T,S           ! stored as scalars 
idim0=sum(var(1:2));
nflds=idim0;
for ivar = 3:length(var)
  for idim = 1:var(ivar)
      nflds=nflds+1;  
      for iel = elmap
          dtemp = fread(infile,npel,realtype);  
          data(:,iel,idim+idim0) = dtemp;
          datastruc(nflds).data(:,iel)=dtemp;
      end
  end
end



%--------------------------------------------------------------------------
% CLOSE FILE
%--------------------------------------------------------------------------
status = fclose(infile);
