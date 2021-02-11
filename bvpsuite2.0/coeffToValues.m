function [ z ] = coeffToValues( coeff, tau, l, rho, teval, d)
%COEFFTOVALUES Evalutes from the coefficients computed by bvpsuite values
%at teval
%   INPUT:  coeff       ... coefficient vector returned by bvpsuite
%           tau         ... grid on which the solution was computed
%           l           ... vector containing the orders of the components
%                           (i.e. [2 3], if the highest order of z1 is 2
%                           and of z2 is 3)
%           rho         ... collocation points given on [0,1]
%           teval       ... points at which solution should be evaluated
%                           (default: tau)
%           d           ... computes the d-th derivative of solution
%                           (default: 0)
%   OUTPUT: z           ... the values corresponding to teval
m = length(rho);
N=length(l);
if(nargin < 5)
    teval = tau;
end
if(nargin < 6)
    d = 0;
end

%Compute coefficents of RK-basis-polynomials in monom-basis
PHI = zeros(max(l),max(l));
for i=1:max(l)
    tmp = phi(i);
    PHI(i,max(l)-length(tmp)+1:max(l))=tmp;
end

PSI = zeros(m,max(l)+1,m+max(l)+1);
for k=1:m
    for i=0:max(l)
        tmp=psi(rho,k,i);
        PSI(k,i+1,end-length(tmp)+1:end)=tmp;
    end
end

z=zeros(N,length(teval));
for tpind = 1:length(teval)
    tp = teval(tpind);
    i=min(length(tau)-1,find(tp>=tau,1,'last'));
    for k=d+1:max(l)
        for j=1:N
            if(k<=l(j))
                z(j,tpind) = z(j,tpind) + coeff((sum(l)+m*N)*(i-1)+sum(min(l,k-1))+sum(l(1:j)>=k))*polyval(PHI(k-d,:),tp-tau(i));
            end
        end
    end
    for k=1:m
        for j=1:N
            if(l(j)-d+1>0)
                z(j,tpind) = z(j,tpind) + coeff((sum(l)+m*N)*(i-1)+sum(l)+N*(k-1)+j) * polyval(reshape(PSI(k,l(j)-d+1,:),m+max(l)+1,1),(tp-tau(i))/(tau(i+1)-tau(i))) * (tau(i+1)-tau(i))^(l(j)-d);
            else
                z(j,tpind) = z(j,tpind) + coeff((sum(l)+m*N)*(i-1)+sum(l)+N*(k-1)+j) * polyval(psi(rho,k,l(j)-d+1),(tp-tau(i))/(tau(i+1)-tau(i))) * (tau(i+1)-tau(i))^(l(j)-d);
            end
        end
    end
end

end


function poly = phi(k)
poly = [1/factorial(k-1),zeros(1,k-1)];
end

function poly = psi(rho,k,order)
poly = 1;
for j=1:length(rho)
    if(j~=k)
        poly = conv(poly,[1/(rho(k)-rho(j)),-rho(j)/(rho(k)-rho(j))]);
    end
end
for j=1:order
    poly = polyint(poly);
end
for j=1:(-order)
    poly=polyder(poly);
end
end

