function [initProfiles,D,V,L,R] = computeEVPStart(problem,settings,numberofguesses)

evalMesh = feval(settings,'mesh');
collMethod = 'gauss';
collPoints = 4;
try
    pathfoll=feval(problem,'pathfollowing');
    lambda_p=pathfoll.start;
catch
    lambda_p=0;
end
interval = feval_problem(problem,'interval',[],[],[],[],[],[],0,0,lambda_p);
x1=linspace(interval(1),interval(2),ceil(5*numberofguesses/collPoints)+1);
evalMesh = (interval(2)-interval(1))/(evalMesh(end)-evalMesh(1))*(evalMesh+evalMesh(1))+interval(1);
m = collPoints;
if (strcmp(collMethod,'user'))
    [rho] = collPoints;
else
    [rho] = getStandardCollocationPoints(collMethod,m);
end
n = feval_problem(problem,'n',[],[],[],[],[],[],0,0,lambda_p);
ordnung = feval_problem(problem,'orders',[],[],[],[],[],[],0,0,lambda_p);
parameter = feval_problem(problem,'parameters',[],[],[],[],[],[],0,0,lambda_p);
interval = feval_problem(problem,'interval',[],[],[],[],[],[],0,0,lambda_p);
c=feval_problem(problem,'c',[],[],[],[],[],[],0,0,lambda_p);
%a=interval(1);
b=interval(2);
N=length(x1)-1;
h(1:N)=x1(2:N+1)-x1(1:N);
tau(1:N,1:m)=x1(1:N).'*ones(1,m)+h(1:N).'*rho(1:m);
faktorielle=[1;1;2;6;24;120;720;5040;40320;362880;3628800;39916800;479001600;6227020800;];


psi=[];
psival=[];
p=zeros(parameter,1);
y=zeros(n,max(ordnung),N);
z=zeros(n,m,N);

phishift =  [0; sum(tril(ones(max(ordnung),max(ordnung)))*(ones(max(ordnung),1)*ordnung>=(1:max(ordnung))'*ones(1,length(ordnung))),2)];

if isempty(psi)% is psi defined by user?
   for ord=1:max(ordnung)
       for i=1:m
           psi(i,1+max(ordnung)-ord:m+max(ordnung),ord)=Psi(i,rho,ord);
       end
   end
end
if isempty(psival)%is psi defined by user?
   for ord=1:max(ordnung)
       for i=1:m
           %evaluation of psi
           psival(ord,i,1:m)=polyval(psi(i,:,ord),rho(1:m));
           psival(ord,i,m+1)=polyval(psi(i,:,ord),1);
           psival(ord,i,m+2)=polyval(psi(i,:,ord),0);
       end
   end
end

%For 0-th order special definitions for psi are necessary

%if min(ordnung)==0 || (length(was)==7 && length(strfind(was,'wertabl'))>0) || (length(was)==8 && length(strfind(was,'residual'))>0)
%   for i=1:m
%s       psi0(i,:,1)=Psi(i,rho,0);
%   end
%end

for j=1:(length(x1)-1)*(n*m+sum(ordnung))+parameter
             x0(j,1)=1;
end
     
nonzeroent=(sum(ordnung)+parameter)*(max(ordnung)*n*max(2,length(c))+m*n+parameter)+...
   N*m*n*max(ordnung)*n+...
   N*m*n*m*n+...
   N*max(ordnung)*n*max(ordnung)*n+...
   N*max(ordnung)*n*m+...
   N*max(ordnung)*n*max(ordnung)*n;
ROWl=zeros(nonzeroent,1);
COLl=zeros(nonzeroent,1);
OUTPUTl=zeros(nonzeroent,1);
ROWr=zeros(nonzeroent,1);
COLr=zeros(nonzeroent,1);
OUTPUTr=zeros(nonzeroent,1);

rowcount = 1;
countl = 1;
countr = 1;

% Randbedingungen links
% Anmerkung: müssen homogen sein!
dBV = feval_problem(problem,'dBV',[],[],[],[],[],[],0,0,lambda_p);
[size_dBV(1),size_dBV(2),size_dBV(3),size_dBV(4)] = size(dBV);

for i=1:size_dBV(2)
   for j=1:n        
       for ll=1:size_dBV(4)   
           OUTPUTl(countl) = dBV(1,i,j,ll);
           ROWl(countl) = rowcount;
           COLl(countl) = phishift(ll)+j;
           countl = countl + 1;
           
           for k=ll:ordnung(j)
               OUTPUTl(countl) = h(end)^(k-ll)/faktorielle(1+k-ll)*dBV(2,i,j,ll);
               ROWl(countl) = rowcount;
               COLl(countl) = (N-1)*(sum(ordnung)+n*m)+phishift(k)+j;
               countl = countl + 1;
           end
           
           for k=1:m
               OUTPUTl(countl) = polyval(Psi(k,rho,ordnung(j)-ll+1),1)*h(end)^(ordnung(j)-ll+1)*dBV(2,i,j,ll);
               ROWl(countl)=rowcount;
               COLl(countl)=(N-1)*(sum(ordnung)+n*m)+sum(ordnung)+(k-1)*n+j;
               countl = countl + 1;
           end
       end
   end
   %fprintf('Randbed in Zeile %d \n',rowcount);
   rowcount = rowcount+1;
end

%Innere Bedingungen: Kollokation und Stetigkeit

for i=1:N % Bedingung in welchem Intervall?
   for k=1:m % Bedingung an welchem Kollokationspunkt?
       tik = tau(i,k);
       dFl = feval_problem(problem,'jacobian',[],[],[],[],tik,p,0,0,lambda_p); % here: lambda = 0
       dFr = dFl-feval_problem(problem,'jacobian',[],[],[],[],tik,p,1,0,lambda_p); % here: lambda = 1
       [size_dF(1),size_dF(2),size_dF(3)] = size(dFl);
       %Linke Seite       
       for ii=1:size_dF(1) % Nummer der Gleichung
          for jj=1:n % Welche Komponente
              %PHI-Koeffizienten
              for ll=1:ordnung(jj) % Bis zur welchen Ableitung
                  %Linke Seite
                 temp = 0; 
                 for kk=0:ll-1 % Welche Ableitung   
                     temp = temp + dFl(ii,jj,kk+1)*(tik-x1(i))^(ll-kk-1)/faktorielle(1+ll-kk-1);
                 end 
                 OUTPUTl(countl) = temp;
                 ROWl(countl)=rowcount;
                 COLl(countl)=(i-1)*(sum(ordnung)+n*m)+phishift(ll)+jj;                 
                 countl = countl + 1;                   
                 
                  %Rechte Seite
                 temp = 0; 
                 for kk=0:ll-1
                     temp = temp + dFr(ii,jj,kk+1)*(tik-x1(i))^(ll-kk-1)/faktorielle(1+ll-kk-1);
                 end 
                 OUTPUTr(countr) = temp;
                 ROWr(countr)=rowcount;
                 COLr(countr)=(i-1)*(sum(ordnung)+n*m)+phishift(ll)+jj;                 
                 countr = countr + 1;                    
              end
              %PSI-Koeffizienten
              for ll=1:m %PSI-Polynome
                  %Linke Seite
                  temp = 0;
                  for kk=0:ordnung(jj) %Ableitung
                      temp = temp + dFl(ii,jj,kk+1)*h(i)^(ordnung(jj)-kk)*polyval(Psi(ll,rho,ordnung(jj)-kk),(tik-x1(i))/h(i));
                  end                  
                  OUTPUTl(countl) = temp;
                  ROWl(countl)=rowcount;
                  COLl(countl)=(i-1)*(sum(ordnung)+n*m)+sum(ordnung)+(ll-1)*n+jj;
                  countl = countl + 1;
                  %Rechte Seite
                  temp = 0;
                  for kk=0:ordnung(jj) %Ableitung
                      temp = temp + dFr(ii,jj,kk+1)*h(i)^(ordnung(jj)-kk)*polyval(Psi(ll,rho,ordnung(jj)-kk),(tik-x1(i))/h(i));
                  end                  
                  OUTPUTr(countr) = temp;
                  ROWr(countr)=rowcount;
                  COLr(countr)=(i-1)*(sum(ordnung)+n*m)+sum(ordnung)+(ll-1)*n+jj;
                  countr = countr + 1;
              end              
          end
          %fprintf('Kollokationsbed in Zeile %d \n',rowcount);
          rowcount = rowcount + 1;
       end
   end
   
   %Stetigkeit
   if(i<N)
       for j=1:n %Komponente
           for k=1:ordnung(j) %(k-1). Ableitung
               %Von links kommend
               for ll=k:ordnung(j)
                   OUTPUTl(countl)=h(i)^(ll-k)/faktorielle(1+ll-k);
                   ROWl(countl)=rowcount;
                   COLl(countl)=(i-1)*(sum(ordnung)+n*m)+phishift(ll)+j;
                   countl = countl + 1;
               end
               
               for ll=1:m
                   OUTPUTl(countl) = polyval(Psi(ll,rho,ordnung(j)-k+1),1)*h(i)^(ordnung(j)-k+1);
                   ROWl(countl)=rowcount;
                   COLl(countl)=(i-1)*(sum(ordnung)+n*m)+sum(ordnung)+(ll-1)*n+j;
                   countl = countl + 1;
               end
               
               %Von rechts kommend
               OUTPUTl(countl)=-1;
               ROWl(countl)=rowcount;
               COLl(countl)=(i)*(sum(ordnung)+n*m)+phishift(k)+j;
               countl = countl + 1;
               
               %fprintf('Stetigkeitsbed in Zeile %d \n',rowcount);
               rowcount = rowcount + 1;
           end
       end
   end
end

countl = countl-1;
countr = countr-1;

ROWl = ROWl(1:countl);
ROWr = ROWr(1:countr);
COLl = COLl(1:countl);
COLr = COLr(1:countr);
OUTPUTl = OUTPUTl(1:countl);
OUTPUTr = OUTPUTr(1:countr);

L = sparse(ROWl,COLl,OUTPUTl,N*(sum(ordnung)+n*m),N*(sum(ordnung)+n*m));
R = sparse(ROWr,COLr,OUTPUTr,N*(sum(ordnung)+n*m),N*(sum(ordnung)+n*m));

[V,D] = eig(full(L),full(R));

D=diag(D);
[~,I] = sort(abs(D));


V=V(:,I(1:N*m));
D=D(I(1:N*m));


%EF = zeros(length(evalMesh),N*m);

for i = 1:N*m
   temp = coeffToValues(V(:,i),x1,ordnung,rho,evalMesh);
   % EF(:,i) = sqrt(length(evalMesh))*EF(:,i)/norm(EF(:,i),2);
   temp =  sqrt(n*length(evalMesh))*temp/norm(temp,2);   
   initProfile.initialMesh = evalMesh;
   initProfile.initialValues = temp;
   interval = feval(problem,'interval');
   if(interval(2) == Inf)
       if(interval(1) == 0 && trafo('splitting'))
           initProfile.initialMesh = unique([initProfile.initialMesh(1:end-1),1./initProfile.initialMesh(2:end)]);
           initProfile.initialValues = [initProfile.initialValues(1:end/2,1:(end-1)),initProfile.initialValues((end/2+1):end,end-1:-1:1)];
       else
           [initProfile.initialMesh,sortInd] = sort(trafo('xi',initProfile.initialMesh,0,interval(1)));
           initProfile.initialValues = initProfile.initialValues(:,sortInd);
       end
   end
   
   initProfile.lambda = D(i);
   initProfile.parameters = [];
   initProfiles{i} = initProfile;
end



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


function ret=Pablspezial(abl,i,komp,x,m,x1,h,y,z,psi,ord) %Polynom der Ableitung abl

faktorielle=[1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800];
sum=0;
for j=1:ord
   if j-abl-1>=0
       sum=sum+((x-x1((i)+1))^(j-abl-1))/(faktorielle(j-abl))*y(komp,j,(i)+1);
   end
end
if ord-abl==0
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

ret=sum(((x-x1((i)+1)).^([abl+1:ord]-abl-1))./(faktorielle([abl+1:ord]-abl)).*y(komp,abl+1:ord,(i)+1))+...
   sum(h((i)+1)^(ord-abl).*(z(komp,1:m,(i)+1).*psival(ord-(abl+1)+1,1:m,k)));

end

function ret=evalP0(i,komp,m,x1,h,y,z,psi,ord)
if (nargin<9) ord=2; end;

     sum=0;
     for j=1:ord

        poly=[1 -x1((i)+1)];
        for k=2:(j-1)
                poly=conv(poly,[1 -x1((i)+1)]);
        end
        if (j-1)==0
            poly=[1];
        end
        zwisch(m+ord-j+1:m+ord)=poly/(factorial(j-1))*y(komp,j,(i)+1);
        sum=sum+zwisch;
    end
for j=1:m
     if ord==0
         help=polyvalpoly(psi(j,:,1),[1 -x1((i)+1)]/h((i)+1));
     else
         help=polyvalpoly(psi(j,:,ord-1+1),[1 -x1((i)+1)]/h((i)+1));
     end
     %Attetion: help has to be of length m+ord - correction for
     %mixed order
     help=help(length(help)-(m+ord-1):length(help));
     zwisch(1:m+ord)=h((i)+1)^ord*(z(komp,j,(i)+1)*help);
     sum=sum+zwisch;
end
ret=sum;
end

function out=polyvalpoly(psi,E)
%polynomial as value in polynomial receiving polynomial

PsiK(1,1)=1;t=1;
for k=2:length(psi)
   t=conv(t,E);
   for m=1:length(t)
     PsiK(k,m)=t(m);
   end
end
PsiE(((length(E)-1))*(length(psi)-1)+1)=PsiK(1,1)*psi(length(psi));
PsiEout=PsiE;
for k=2:length(psi)
   PsiE(((length(E)-1))*(length(psi)-1)+1-(length(E)-1)*(k-1):((length(E)-1))*(length(psi)-1)+1)=...
                      psi(length(psi)-k+1)*PsiK(k,1:(length(E)-1)*(k-1)+1);
   PsiEout=PsiEout+PsiE;
end

out=PsiEout;
end

function Psireturn=Psi(n,rho,nr)
 i=length(rho);
 prod=[1];
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

function [A] = trafomatrix(t,n)
n=n-1;
A=diag([1,(-1).^(1:n)./t.^(2*(1:n))]);
if(n>1)
    A(1+(2:n),2) = ((-1).^(2:n).*factorial(2:n)./t.^(1+(2:n)))';
    for i=4:n+1
        A(i,3:i-1) = -(A(i-1,2:i-2)/t^2 + A(i-1,3:i-1)/t);
    end
end
end

