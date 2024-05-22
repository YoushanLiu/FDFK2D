function [seis, nt, nx, dt, errflag] = readsu(infile, domain)
%
%readsu   read su format file
%
%   writen by Youshan Liu
%
%   infile:   filename
%   seis:     array of the data
%   nt:       number of sampling points
%   nx:       number of traces
%   dt:       sampling interval
%   domain:   0 -> time domain (default); 1 -> depth domain
%
%   [seis, nt, nx, dt, errflag] = readsu(infile)
%   read a su file in time-domain and return data in seis array
%
%   [seis, nt, nx, dt, errflag] = readsu(infile, 1)
%   read a su file in spatial-domain and return data in seis array
%


if (nargin < 1)
    fprintf('readsu: Not enough inputs ! \n');
    help readsu
end

if (nargin <= 1)
    domain = 0;
end



fid = fopen(infile, 'rb', 'n');

if (-1 == fid)

    fprintf('Error: Cannot open file %s ! \n', infile);
    errflag = 0;
    seis = [];
    nx = 0;
    nt = 0;
    dt = -1;
    return

else

    fseek(fid, 0, 'eof');
    nbyte = ftell(fid);
    if (0 == nbyte)
        fprintf('Error: %s file is NULL !', infile);
        errflag = 0;
        seis = [];
        nx = 0;
        nt = 0;
        dt = -1;
        return
    else
        errflag = 1;
    end

    frewind(fid);

    % get header infromation
	head = fread(fid, 120, 'uint16');
	nt = head(58);
	if (~domain || isempty(domain))
  		dt = head(59)*1e-6;
    else
  		dt = head(59)*1e-3;
	end

	nx = fix(nbyte/(240 + 4*nt));

	frewind(fid);

	seis = zeros(nt, nx);
    for ix = 1:1:nx
        fseek(fid, 240, 0);
        seis(1:nt,ix) = fread(fid, nt, 'float32');
    end

end

fclose(fid);

end
