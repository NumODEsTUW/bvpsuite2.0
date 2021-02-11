function [ initProfile ] = initial_coefficients( problem,x1,initProfile,rho,transform )
% Compute polynomial coefficients for the polynomial interpolating the
% values in initProfile

m=length(rho);
n=feval_problem(problem,'n'); %number of solution components
ordnung=feval_problem(problem,'orders');
k=n;

if(nargin<5)
    transform = 1;
end

if(transform == 0)
    interval=feval_problem(problem,'interval');
else
    interval=feval(problem,'interval');
end

if(trafo('splitting') && interval(1)==0 && interval(2)==Inf)
    k=n/2;
end

N=length(x1)-1;
initialmesh = initProfile.initialMesh;
initialvalues = initProfile.initialValues;

parameter=feval_problem(problem,'parameters');

faktorielle=[1;1;2;6;24;120;720;5040;40320;362880;3628800;39916800;479001600;6227020800;];

if(parameter>0)
    par = initProfile.parameters;
else
    par=[];
end

if max(ordnung)>0
    for ord=1:max(ordnung)
        for i=1:m
            psi(i,1+max(ordnung)-ord:m+max(ordnung),ord)=Psi(i,rho,ord);
        end
    end
end
if min(ordnung)==0
    for i=1:m
        psi0(i,:,1)=Psi(i,rho,0);
    end
end
for i=0:N-1
    h(i+1)=x1(i+2)-x1(i+1);
end


for ord=1:max(ordnung)
    proj_intervallstellen=0:1/(m+ord-1):1;
    for i=1:m
        psival(ord,i,1:m+ord)=polyval(psi(i,:,ord),proj_intervallstellen);
    end
end

if ~isempty(initialvalues)
    
    %removal of entries that appear twice in "initialmesh"
    jj=0;
    for ii=2:length(initialmesh)
        if initialmesh(ii)~=initialmesh(ii-1)
            jj=jj+1;
            initialmesh2(jj)=initialmesh(ii-1);
            initialvalues2(1:k,jj)=initialvalues(1:k,jj);
        end
    end
    
    initialmesh2(jj+1)=initialmesh(end);
    
    initialvalues2(1:k,jj+1)=initialvalues(1:k,end);
end

for p=1:n
    ordp=ordnung(p);
    interpolationsstellen=zeros(1,(N-1)*(ordp+m-1)+1);
    
    if(interval(2)==Inf)
        if(trafo('splitting'))
            if p>k %components on [1,infty)
                
                ind_1=find(initialmesh2>=1);
                
                if length(ind_1)<=1
                    error('initial2');
                end
                
                %transformation to the interval [0,1]
                
                %********************************************************
                %define a computational grid according to the collocation
                %rules
                %********************************************************
                for i=0:N-1
                    %constructs grid (=mesh + coll.points)
                    %inserts coll. points between elements of x1
                    interpolationsstellen(i*(m+ordp-1)+1:(i+1)*(m+ordnung(p)-1)+1)=x1(i+1):h(i+1)/(m+ordp-1):x1(i+2);
                end
                
                %****************
                %Interpolation
                %***************
                
                if ind_1(1)-1 <=0
                    error('initial1');
                end
                
                if abs(initialmesh2(ind_1(1)-1)-0)< 10^-4
                    initialmesh2(ind_1(1)-1)=10^-4;
                end
                
                interpolationsstellen=[interpolationsstellen,1./initialmesh2(ind_1(1)-1)];
                interpolationsstellenwerte=interp1([1./initialmesh2(ind_1(1)-1),1./initialmesh2(ind_1)],[initialvalues2(p-k,ind_1(1)-1),initialvalues2(p-k,ind_1)],interpolationsstellen,'linear','extrap');
                
            else
                %if the number of the solution component is <=n/2
                %then no transformation has to be carried out
                
                %
                ind_1=find(initialmesh2<=1);
                
                if length(ind_1)<=1
                    error('initial1');
                end
                
                for i=0:N-1
                    %construction of the computational grid (=mesh + coll.points)
                    %inserts coll. points between elements of x1
                    interpolationsstellen(i*(m+ordp-1)+1:(i+1)*(m+ordnung(p)-1)+1)=x1(i+1):h(i+1)/(m+ordp-1):x1(i+2);
                end
                
                if ind_1(end)+1 > length(initialmesh2)
                    error('initial2');
                end
                
                interpolationsstellen=[interpolationsstellen,initialmesh2(ind_1(end)+1)];
                interpolationsstellenwerte=interp1([initialmesh2(ind_1),initialmesh2(ind_1(end)+1)],[initialvalues2(p,ind_1),initialvalues2(p,ind_1(end)+1)],interpolationsstellen,'spline');
                
                %figure
                %plot(interpolationsstellen,interpolationsstellenwerte,'r')
                
            end
        else
            %the given problem is posed on [a,infty) a ~= 0
            
            ind_1=find(initialmesh2>=1);
            
            %transformation to the interval [0,1]
            
            %********************************************************
            %define a computational grid according to the collocation
            %rules
            %********************************************************
            for i=0:N-1
                %constructs grid (=mesh + coll.points)
                %inserts coll. points between elements of x1
                interpolationsstellen(i*(m+ordp-1)+1:(i+1)*(m+ordnung(p)-1)+1)=x1(i+1):h(i+1)/(m+ordp-1):x1(i+2);
            end
            %****************
            %Interpolation
            %***************
            interpolationsstellenwerte=interp1(trafo('xi_inv',initialmesh2(ind_1),0,interval(1)),initialvalues2(p,ind_1),interpolationsstellen,'linear','extrap');
        end
        
    else
        
        for i=0:N-1
            
            % x1 ... new mesh
            %constructs grid (=mesh + coll.points)
            %the number of added collocation points varies with the order
            %of derivative of the solution component z_i
            interpolationsstellen(i*(m+ordp-1)+1:(i+1)*(m+ordnung(p)-1)+1)=x1(i+1):h(i+1)/(m+ordp-1):x1(i+2);
        end
        initialvalues2(p,:);
        interpolationsstellenwerte=interp1(initialmesh2,initialvalues2(p,:),interpolationsstellen,'spline');
    end
    
    
    %***********************************************
    %calculates coefficients in Runge-Kutta basis!!
    %***********************************************
    
    
    for i=0:N-1
        if (m==1 && ordnung(p)==0)
            intervallstellen(1)=x1((i+1));
        else
            intervallstellen=x1(i+1):h(i+1)/(m+ordnung(p)-1):x1(i+2);
        end
        % if length(koeff)==0
        %intervallwerte=interp1(stellen,werte(p,:),intervallstellen,'spline');
        intervallwerte=interpolationsstellenwerte(i*(m+ordnung(p)-1)+1:(i+1)*(m+ordnung(p)-1)+1);
        %  else
        %      intervallwerte=equations('wert',koeff,bvpfile,stellen,intervallstellen);
        %      intervallwerte=intervallwerte(p,:);
        %  end
        
        A=zeros(m+ordnung(p),ordnung(p)+m);
        for j=1:m+ordnung(p)
            for q=1:ordnung(p)
                A(j,q)=(intervallstellen(j)-x1((i)+1))^(q-1)/faktorielle(q);
            end
            
            for sp=1:m
                if ordp>0
                    %A(j,sp+ordnung(p))=h((i)+1)^ordnung(p)*polyval(psi(sp,:,ordnung(p)),(intervallstellen(j)-x1(i+1))/h(i+1));
                    A(j,sp+ordnung(p))=h((i)+1)^ordnung(p)*psival(ordnung(p),sp,j);
                else
                    A(j,sp)= polyval(psi0(sp,:),(intervallstellen(j)-x1((i)+1))/h((i)+1));
                end
            end
        end
        help=A\intervallwerte';
        intervallstart((i)+1,1:length(help),p)=help;
    end
end

j=0;
coeff=zeros(N*max(ordnung)*n+N*m*n+parameter,1);
for i=0:N-1
    for k=1:max(ordnung)
        for p=1:n
            if k<=ordnung(p)
                j=j+1;
                coeff(j,1)=intervallstart((i)+1,k,p);
            end
        end
    end
    for k=1:m
        for p=1:n
            j=j+1;
            coeff(j,1)=intervallstart((i)+1,k+ordnung(p),p);
        end
    end
end
for pii=1:parameter
    j=j+1;coeff(j,1)=par(pii);
end
initProfile.initialCoeff=coeff(1:j,1);


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

