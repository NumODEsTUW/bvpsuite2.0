function [x1tau,valx1tau,solStruct] = solveLinearProblem(problem,x1,rho)

m=length(rho);
n = feval_problem(problem,'n');
ordnung = feval_problem(problem,'orders');
parameter = feval_problem(problem,'parameters');
interval = feval_problem(problem,'interval');
c=feval_problem(problem,'c');
a=interval(1);
b=interval(2);
N=length(x1)-1;

[ctmp,sortind] = sort(c);
cint=zeros(length(ctmp),1);
for i=1:length(ctmp)
    j=0;
    if (ctmp(i)<a) || (ctmp(i)>b)
        err('equations_err1');
    end
    while ctmp(i)>=x1((j)+1) && j<N
        j=j+1;
    end
    cint(i)=j-1;
end
cint(sortind)=cint;

h(1:N)=x1(2:N+1)-x1(1:N);
tau(1:N,1:m)=x1(1:N).'*ones(1,m)+h(1:N).'*rho(1:m);
faktorielle=[1;1;2;6;24;120;720;5040;40320;362880;3628800;39916800;479001600;6227020800;];

p=zeros(parameter);

psi=zeros(m,max(ordnung)+m,max(ordnung));
for ord=1:max(ordnung)
    for i=1:m
        psi(i,1+max(ordnung)-ord:m+max(ordnung),ord)=Psi(i,rho,ord);
    end
end

psival=zeros(max(ordnung),m,m+2);
for ord=1:max(ordnung)
    for i=1:m
        %evaluation of psi
        psival(ord,i,1:m)=polyval(psi(i,:,ord),rho(1:m));
        psival(ord,i,m+1)=polyval(psi(i,:,ord),1);
        psival(ord,i,m+2)=polyval(psi(i,:,ord),0);
    end
end

%For 0-th order special definitions for psi are necessary

%if min(ordnung)==0 || (length(was)==7 && length(strfind(was,'wertabl'))>0) || (length(was)==8 && length(strfind(was,'residual'))>0)
%   for i=1:m
%s       psi0(i,:,1)=Psi(i,rho,0);
%   end
%end

for i=1:m
    psi0(i,:,1)=Psi(i,rho,0);
end

x0=zeros((length(x1)-1)*(n*m+sum(ordnung))+parameter,1);
for j=1:(length(x1)-1)*(n*m+sum(ordnung))+parameter
    x0(j,1)=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTE RESIDUUM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
residuum=functionFDF( 'F', problem,zeros(length(x0),1),x1,psival,psi,rho,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTE SYSTEMMATRIX%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
systemmatrix=functionFDF( 'DF', problem,zeros(length(x0),1),x1,psival,psi,rho,[]);

global M1
M1 = systemmatrix;

a=-systemmatrix\residuum;

j=0;
y=zeros(n,max(ordnung),N);
z=zeros(n,m,N);
for i=0:N-1
    for q=1:max(ordnung)
        for ni=1:n
            if q<=ordnung(ni)
                j=j+1;y(ni,q,i+1)=a(j);
            end
        end
    end
    z(1:n,1:m,i+1)=reshape(a(j+1:j+n*m),n,m);
    j=j+n*m;
end
p(1:parameter)=a(j+1:j+parameter).';
output=zeros(n,N+N*m+1);
for i=1:n
    stelle=0;
    for j=0:N-1
        if ordnung(i)~=0
            if j~=0
                help=Poptimiert(0,j,i,x1((j)+1),m,x1,h,y,z,psival,ordnung(i),m+2);
            else
                help=P(0,j,i,x1((j)+1),m,x1,h,y,z,psi,ordnung(i));
            end
        else
            help=P(0,j,i,x1((j)+1),m,x1,h,[],z,psi0,0);
        end
        stelle=stelle+1;output(i,stelle)=help;
        for k=1:m
            if ordnung(i)~=0
                help=Poptimiert(0,j,i,tau((j)+1,k),m,x1,h,y,z,psival,ordnung(i),k);
            else
                help=P(0,j,i,tau((j)+1,k),m,x1,h,[],z,psi0,0);
            end
            stelle=stelle+1;output(i,stelle)=help;
        end
    end
    if ordnung(i)~=0
        help=P(0,N-1,i,x1((N)+1),m,x1,h,y,z,psi,ordnung(i));
    else
        help=P(0,N-1,i,x1((N)+1),m,x1,h,[],z,psi0,0);
    end
    stelle=stelle+1;
    output(i,stelle)=help;
end
valx1tau=output;

stelle=0;
x1tau=zeros(1,N);
for j=0:N-1
    stelle=stelle+1;x1tau(stelle)=x1((j)+1);
    for k=1:m
        stelle=stelle+1;x1tau(stelle)=tau((j)+1,k);
    end
end
stelle=stelle+1;
x1tau(stelle)=x1((N)+1);

solStruct.x1 = x1tau(1:m+1:end);
solStruct.valx1 = valx1tau(:,1:m+1:end);
solStruct.x1tau = x1tau;
solStruct.valx1tau = valx1tau;
solStruct.parameters = p;
solStruct.coeff = a(1:end-length(p));
end




function ret=P(abl,i,komp,x,m,x1,h,y,z,psi,ord) %Polynomial of derivation abl

faktorielle=[1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800];
sum=0;
for j=1:ord
    if j-abl-1>=0
        sum=sum+((x-x1((i)+1))^(j-abl-1))/(faktorielle(j-abl))*y(komp,j,(i)+1);
    end
end
if ord==0
    for j=1:m
        sum=sum+h((i)+1)^(ord-abl)*(z(komp,j,(i)+1)*polyval(psi(j,:,ord-(abl+1)+1+1),(x-x1((i)+1))/h((i)+1)));
    end
else
    for j=1:m
        sum=sum+h((i)+1)^(ord-abl)*(z(komp,j,(i)+1)*polyval(psi(j,:,ord-(abl+1)+1),(x-x1((i)+1))/h((i)+1)));
    end
end
ret=sum;

end

function ret=Poptimiert(abl,i,komp,x,m,x1,h,y,z,psival,ord,k)
%Polynomial opitmized version
faktorielle=[1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800];

ret=sum(((x-x1((i)+1)).^((abl+1:ord)-abl-1))./(faktorielle((abl+1:ord)-abl)).*y(komp,abl+1:ord,(i)+1))+...
    sum(h((i)+1)^(ord-abl).*(z(komp,1:m,(i)+1).*psival(ord-(abl+1)+1,1:m,k)));

end

function Psireturn=Psi(n,rho,nr)
i=length(rho);
prod=1;
for s=1:i
    if (s~=n)
        prod=conv(prod,[1 -rho(s)])/(rho(n)-rho(s));
    end
end
for s=1:(nr)
    prod=polyint(prod);
end
Psireturn=prod;
end

