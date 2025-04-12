function writesu(outfile, seis, dt, domain)
%
%writesu   write data into su format file
%
%   written by Youshan Liu
%
%   seis:     array of the data
%   outfile:  filename
%   nt:       number of sampling points
%   nx:       number of traces
%   dt:       sampling interval
%   domain:   0 -> time domain (default); 1 -> depth domain
%
%   writesu(outfile, seis)
%   write the data stored in seis array into su file with dt = 1
%
%   writesu(outfile, seis, dt)
%   write the data stored in seis array into su file with dt as temporal sampling interval
%
%   writesu(outfile, seis, dt, 1)
%   write the data stored in seis array into su file with dt as spatial sampling interval
%


if (nargin < 2)
    fprintf('writesu: Not enough inputs ! \n');
    help writesu
end

if (nargin <= 2)
    dt = 1;
end

if (nargin <= 3)
    domain = 0;
end


[nt, nx] = size(seis);


head = zeros(120, 1);

head(58) = nt;
if (isempty(domain) || ~domain)
    head(59) = dt*1e6;
else
    head(59) = dt*1e3;
end


if (0 ~= exist(outfile, 'file'))
   delete(outfile);
end

fid = fopen(outfile, 'wb', 'n');
   for ix = 1:1:nx
       fwrite(fid, head(1:120), 'uint16');
       fwrite(fid, seis(1:nt,ix), 'float32');
   end
fclose(fid);

end
