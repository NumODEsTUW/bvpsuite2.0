function ret = trafo(varargin)
switch varargin{1}
    case 'xi'
        ret = xi(varargin{2},varargin{3},varargin{4});
    case 'xi_inv'
        ret = xi_inv(varargin{2},varargin{3},varargin{4});
    case 'trafomatrix'
        ret = trafomatrix(varargin{2},varargin{3},varargin{4});
    case 'trafomatrix_inv'
        ret = trafomatrix_inv(varargin{2},varargin{3},varargin{4});
    case 'splitting'        
        ret = 1;
end
end

function ret = xi(tau,k,a)
%is the k-th derivative of the transformation from [0,1] on [a,infty)
alpha=-1;
if(a>0)
    ret=a./tau.^(k+1).*(-1)^(k);
else
    if(k>0)        
        ret=prod(alpha-(0:k-1))*(1-tau).^(alpha-k);
    else
        ret=(1-tau).^alpha-1;
    end
end
end

function ret = xi_inv(t,k,a)
%is the k-th derivative of the transformation from [a,infty) on [0,1]
alpha=-1;
if(a>0)
    ret=a./t.^(k+1)*(-1)^(k);
else
    if(k>0)
        ret=prod(1/alpha-(0:k-1))*(1+t).^(1/alpha-k);
    else
        ret=1-(1+t).^(1/alpha);
    end
end
end

function [A] = trafomatrix(tau,n,a)
%trafo from z* defined on [0,1] on z defined on [a,infty)
global chainrulecoeffs;
A=zeros(n,n);
xitmp = zeros(n+1,1);
for i=1:n
    xitmp(i) = xi(tau,i-1,a);
end
for i=1:n     
    b=chainrulecoeffs{i}.b;    
    for j=1:length(b)      
       A(i,chainrulecoeffs{i}.b(j)+1) = chainrulecoeffs{i}.c(j)*prod(xitmp(1:i)'.^(chainrulecoeffs{i}.A(j,:)));
    end
end
end

function [A] = trafomatrix_inv(tau,n,a)
%trafo from z defined on [a,infty) on z* defined on [0,1]
B=trafomatrix(tau,n,a);
if(norm(B)<Inf)
    A=inv(trafomatrix(tau,n,a));
else
    s=size(B);
    A=zeros(s(1),s(2));
end
A(1,:)=0;
A(:,1)=0;
A(1,1)=1;
end

