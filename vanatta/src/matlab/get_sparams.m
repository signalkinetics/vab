function sp = get_sparams(varargin)
% function accepts ABCD params of two port networks in front of ports 1 and
% 2
if mod(nargin,2) == 0
    error("Odd # of arguments required. Make sure # of nodes is even.")
end

N = nargin-1;
ports = varargin{nargin};
r1 = ports(1);
r2 = ports(2);

sp = zeros(size(varargin{1}),'like',varargin{1});

for n=1:N/2
    T1 = varargin{n};
    T2 = varargin{N-n+1};

    T = pagemtimes(T1,T2);

    Z = abcd2z(T);
    Y = z2y(Z);
    
    Zr = [r1 0; 0 r2];
    F = [1/(2*sqrt(r1)) 0; 0 1/(2*sqrt(r2))];
    
    i = eye(2);
    
%     sp(n:(N-2*n+1):(N-n+1),n:(N-2*n+1):(N-n+1),:) = F*(i-conj(Zr)*Y)*(i+Zr*Y)^(-1)*F^(-1);
    sp(n:(N-2*n+1):(N-n+1),n:(N-2*n+1):(N-n+1),:) = pagemtimes(F,pagemtimes(i-pagemtimes(conj(Zr),Y),pagemtimes(pageinv(i+pagemtimes(Zr,Y)),F^(-1))));
end


