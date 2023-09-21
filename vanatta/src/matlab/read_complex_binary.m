%
% Copyright 2001 Free Software Foundation, Inc.
%
% This file is part of GNU Radio
%
% GNU Radio is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3, or (at your option)
% any later version.
%
% GNU Radio is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with GNU Radio; see the file COPYING.  If not, write to
% the Free Software Foundation, Inc., 51 Franklin Street,
% Boston, MA 02110-1301, USA.
%

function v = read_complex_binary (filename, start, count)

%% usage: read_complex_binary (filename, [count])
%%
%%  open filename and return the contents as a column vector,
%%  treating them as 32 bit complex numbers
%%

narginchk (1,3);

if (nargin < 2)
    start = 0;
    count = Inf;
else % start = size(float)*start [basically start is in terms of bytes]
    start = start*4*length(size('float')); % 2: for real and imaginary, 4 for float
    %start = start*4*length(size('double')); % 2: for real and imaginary, 4 for float
end

f = fopen (filename, 'rb');
if (f < 0)
    v = 0;
else
    r=fseek(f, start, -1);
    if(r<0)
        fprintf(2, 'Error reading start\n');
        v=0;
        return;
    end
    t = fread (f, [2, count], 'float');
    fclose (f);
    v = t(1,:) + t(2,:)*1i;
    [r, c] = size (v);
    v = reshape (v, c, r);
end

if numel(v)<2
    error('File %s does not exist',filename);
end
end
