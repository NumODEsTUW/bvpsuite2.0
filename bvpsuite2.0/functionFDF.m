function [ output ] = functionFDF( request, problem,a,x1,psival,psi,rho,predictor)
% request can be either 'F' to compute F(a)
%                    or 'DF' to compute DF(a)

m=length(rho);
n = feval_problem(problem,'n');
ordnung = feval_problem(problem,'orders');
parameter = feval_problem(problem,'parameters');
interval = feval_problem(problem,'interval');
interval_feval = feval(problem,'interval'); % needed for DF in case of pathfollowing
a1=interval(1);
b=interval(2);
c=feval_problem(problem,'c');
N=length(x1)-1;

[ctmp,sortind] = sort(c);
cint=zeros(length(ctmp),1);
for i=1:length(ctmp)
    j=0;
    if (ctmp(i)<a1) || (ctmp(i)>b)
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
faktorielle=[1;1;2;6;24;120;720;5040;40320;362880;3628800;39916800;479001600;6227020800];

%-----------------------------------------------
%entries of "a" will be transformed to the names of the documentation
%y,z,p see Runge-Kutta Basis in the manual!
j=0;
y=zeros(n,max(ordnung),N);
z=zeros(n,m,N);
for i=0:N-1
    for q=1:max(ordnung)
        for ni=1:n
            if q<=ordnung(ni)
                j=j+1; y(ni,q,(i)+1)=a(j);
            end
        end
    end
    z(1:n,1:m,(i)+1)=reshape(a(j+1:j+n*m),n,m);
    j=j+n*m;
end
p(1:parameter)=a(j+1:j+parameter).';

n_tmp = N*(sum(ordnung)+n*m)+parameter;
lambda = [];
if length(a)>n_tmp
    lambda_p=a(end);
    if feval(problem,'EVP')
        lambda = a(end-1);
    end
else
    try
        pathfoll=feval(problem,'pathfollowing');
        lambda_p=pathfoll.start;
    catch
        lambda_p=0;
    end
end

switch request
    % ***********************************************************************************
    case 'F' % *************************** Computation of F(X) ************************************
        % ******************************************************************************
        output=zeros(length(a),1);
        j=0;
        if isempty(c)
            fza=zeros(n,max(ordnung));
            fzb=zeros(n,max(ordnung));
            for oi=1:max(ordnung)
                for ni=1:n
                    if oi<=ordnung(ni)
                        fza(ni,oi)=y(ni,oi,(0)+1);
                        fzb(ni,oi)=Poptimiert(oi-1,N-1,ni,b,m,x1,h,y,z,psival,ordnung(ni),m+1);
                    end
                end
            end
            Rx=feval_problem(problem,'BV',[],fza,fzb,[],[],p,lambda,[],lambda_p);
            if sum(ordnung)+parameter>0
                output(1+j:sum(ordnung)+parameter+j,1)=Rx(1:sum(ordnung)+parameter);
                j=j+sum(ordnung)+parameter;
            end
        else
            fzc=zeros(n,max(ordnung),length(c));
            for oi=1:max(ordnung)
                for ni=1:n
                    if oi<=ordnung(ni)
                        for ci=1:length(c)
                            fzc(ni,oi,ci)=P(oi-1,cint(ci),ni,c(ci),m,x1,h,y,z,psi,ordnung(ni));
                        end
                    end
                end
            end
            Rx=feval_problem(problem,'BV',[],[],[],fzc,[],p,lambda,[],lambda_p);
            
            if sum(ordnung)+parameter>0
                output(1+j:sum(ordnung)+parameter+j,1)=Rx(1:sum(ordnung)+parameter);
                j=j+sum(ordnung)+parameter;
            end
        end
        
        %-----------------------------------------------
        %Evaluation of the equations for the solver
        fz=zeros(n,max(ordnung)+1);
        for i=0:N-1
            for k=1:m
                for oi=0:max(ordnung)-1
                    for ni=1:n
                        if oi<ordnung(ni)
                            fz(ni,oi+1)=Poptimiert(oi,i,ni,tau((i)+1,k),m,x1,h,y,z,psival,ordnung(ni),k); % psi wird immer von der Ordnung der entsprechenden Gleichung gewählt
                        end
                    end
                end
                for ni=1:n
                    fz(ni,ordnung(ni)+1)=z(ni,k,(i)+1);
                end
                if n>0
                    g_=feval_problem(problem,'problem',fz,[],[],[],tau((i)+1,k),p,lambda,[],lambda_p);
                    output(j+1:j+n,1)=g_(1:n);
                    j=j+n;
                end
            end
            if (i<N-1)
                for oi=0:max(ordnung)-1
                    for ni=1:n
                        if oi<=ordnung(ni)-1
                            j=j+1;u=Poptimiert(oi,i,ni,x1((i+1)+1),m,x1,h,y,z,psival,ordnung(ni),m+1)-y(ni,oi+1,(i+1)+1);
                            output(j,1)=u;
                        end
                    end
                end
            end
        end
        
        % In case of pathfollowing
        if length(a)>n_tmp && ~isempty(predictor)
            if ~predictor.lpfix
                output(N*(sum(ordnung)+n*m)+parameter+1,1)=(a.'-predictor.a_0)*predictor.tangent.'-predictor.steplength;
            else
                output(N*(sum(ordnung)+n*m)+parameter+1,1)=0;
            end
        end
        
        % ***********************************************************************************
    case 'DF' % *************************** Computation of DF(X) ************************************
        % ******************************************************************************
        
        %----Start Boundary conditions
        %R_0
        j=0;
        % Dimensioning of non-zero entries
        % Boundary conditions + Matrix
        nonzeroent=(sum(ordnung)+parameter)*(max(ordnung)*n*max(2,length(c))+m*n+parameter)+...
            N*m*n*max(ordnung)*n+...
            N*m*n*m*n+...
            N*max(ordnung)*n*max(ordnung)*n+...
            N*max(ordnung)*n*m+...
            N*max(ordnung)*n*max(ordnung)*n;
        ROW=zeros(nonzeroent,1);
        COL=zeros(nonzeroent,1);
        OUTPUT=zeros(nonzeroent,1);
        gen_c=0; %General count for sparse matrix efficiency
        if isempty(c)
            fza=zeros(n,max(ordnung));
            fzb=zeros(n,max(ordnung));
            for oi=1:max(ordnung)
                for ni=1:n
                    if oi<=ordnung(ni)
                        fza(ni,oi)=y(ni,oi,(0)+1);
                        fzb(ni,oi)=Poptimiert(oi-1,N-1,ni,b,m,x1,h,y,z,psival,ordnung(ni),m+1);
                    end
                end
            end
            % p1 index for boundary conditions of the equations, oi index for the
            % derivatives
            %p2 the "respect variables" of the differentiation have more than one
            %component
            DpR=feval_problem(problem,'dP_BV',[],fza,fzb,[],[],p,lambda,[],lambda_p);
            dBV=feval_problem(problem,'dBV',[],fza,fzb,[],[],p,lambda,[],lambda_p);
            [size_dBV(1),size_dBV(2),size_dBV(3),size_dBV(4)]=size(dBV);
            for p1=1:sum(ordnung)+parameter
                DRa_=zeros(max(ordnung),size_dBV(3));
                DRb_=zeros(max(ordnung),size_dBV(3));
                for k=1:max(ordnung)
                    temp = reshape(dBV(1,p1,:,k),size_dBV(3),1);
                    for l=1:size(temp)
                        DRa_(k,l)=temp(l);
                    end
                    temp = reshape(dBV(2,p1,:,k),size_dBV(3),1);
                    for l=1:size(temp)
                        DRb_(k,l)=temp(l);
                    end
                end
                
                if (N~=1)
                    j=j+1;
                    zaehler=0;
                    zaehler2=0;
                    %the derivatives with respect to y_0k follow
                    for k=1:max(ordnung)
                        help = dBV(1,p1,:,k);
                        for p2=1:n
                            if k<=ordnung(p2)
                                zaehler=zaehler+1;
                                gen_c=gen_c+1;
                                ROW(gen_c)=j;
                                COL(gen_c)=zaehler;
                                OUTPUT(gen_c)=help(p2);
                            end
                        end
                    end
                    %the derivatives with respect to y_(N-1)k follow
                    for k=1:max(ordnung)
                        for p2=1:n
                            if k<=ordnung(p2)
                                zaehler2=zaehler2+1;
                                help=0;
                                for abl=0:k-1
                                    help=help+DRb_(abl+1,p2)*((b-x1((N-1)+1))^(k-abl-1)/faktorielle(k-abl));
                                end
                                if (k-abl-1)>=0
                                    col=(N-1)*(sum(ordnung)+m*n)+zaehler2;
                                    gen_c=gen_c+1;
                                    ROW(gen_c)=j;
                                    COL(gen_c)=col;
                                    OUTPUT(gen_c)=help;
                                else
                                    col=(N-1)*(sum(ordnung)+m*n)+zaehler2;
                                    gen_c=gen_c+1;
                                    ROW(gen_c)=j;
                                    COL(gen_c)=col;
                                    OUTPUT(gen_c)=0;
                                end
                            end
                        end
                    end
                else
                    j=j+1;
                    zaehler2=0;
                    for k=1:max(ordnung)
                        for p2=1:n
                            if k<=ordnung(p2)
                                zaehler2=zaehler2+1;
                                help=0;
                                for abl=0:k-1
                                    help=help+DRb_(abl+1,p2)*((b-a1)^(k-abl-1)/faktorielle(k-abl));
                                end
                                if (k-abl-1)>=0
                                    gen_c=gen_c+1;
                                    ROW(gen_c)=j;
                                    COL(gen_c)=zaehler2;
                                    OUTPUT(gen_c)=DRa_(k,p2)+help;
                                else
                                    gen_c=gen_c+1;
                                    ROW(gen_c)=j;
                                    COL(gen_c)=zaehler2;
                                    OUTPUT(gen_c)=DRa_(k,p2);
                                end
                            end
                        end
                    end
                end
                %derivatives with respect to z_(N-1)q
                for q=1:m
                    for p2=1:n
                        zaehler2=zaehler2+1;
                        help=0;
                        for abl=0:ordnung(p2)-1
                            help=help+DRb_(abl+1,p2)*h((N-1)+1)^(ordnung(p2)-abl)*polyval(psi(q,:,ordnung(p2)-(abl+1)+1),1);
                        end
                        col=(N-1)*(sum(ordnung)+m*n)+zaehler2;
                        gen_c=gen_c+1;
                        ROW(gen_c)=j;
                        COL(gen_c)=col;
                        OUTPUT(gen_c)=help;
                    end
                end
                %derivatives with respect to the unknown parameters p_1 to
                %p_parameter
                if parameter>0
                    col=N*(n*m+sum(ordnung))+1:N*(n*m+sum(ordnung))+parameter;
                    gen_c_start=gen_c+1;
                    gen_c=gen_c+parameter+1;%length(col);
                    ROW(gen_c_start:gen_c)=j;
                    COL(gen_c_start:gen_c)=col;
                    OUTPUT(gen_c_start:gen_c)=DpR(p1,1:parameter);
                end
            end
        else   % for additional conditions
            fzc=zeros(n,max(ordnung),length(c));
            for oi=1:max(ordnung)
                for ni=1:n
                    if oi<=ordnung(ni)
                        for ci=1:length(c)
                            fzc(ni,oi,ci)=P(oi-1,cint(ci),ni,c(ci),m,x1,h,y,z,psi,ordnung(ni));
                        end
                    end
                end
            end
            % p1 index for boundary conditions of the equations, oi index for the
            % derivatives
            %p2 the "respect variables" of the differentiation have more than one
            %component
            %derivatives with respect to z_(N-1)q
            DpR=feval_problem(problem,'dP_BV',[],[],[],fzc,[],p,lambda,[],lambda_p);
            dBV=feval_problem(problem,'dBV',[],[],[],fzc,[],p,lambda,[],lambda_p);
            for p1=1:sum(ordnung)+parameter
                DR_=zeros(max(ordnung),n,length(c));
                for k=1:max(ordnung)
                    for ci=1:length(c)
                        temp = dBV(ci,p1,:,k);
                        for l=1:length(temp)
                            DR_(k,l,ci)=temp(l);
                        end
                    end
                end
                % Derivation with respect to y_c1_k, z_c1_k, y_c2_k, z_c2_k, ...
                % -----------------------------------------------------------
                j=j+1;
                for ci=1:length(c)
                    zaehler2=0;
                    for k=1:max(ordnung)
                        for p2=1:n
                            if k<=ordnung(p2)
                                zaehler2=zaehler2+1;
                                help2=0;
                                % Every polynomial containing the variable y_ci has
                                % to be differentiated and added
                                for cii=1:length(c)
                                    if cint(cii)==cint(ci)
                                        help=0;
                                        for abl=0:k-1
                                            help=help+DR_(abl+1,p2,cii)*((c(cii)-x1((  cint(cii)  )+1))^(k-abl-1)/faktorielle(k-abl));
                                        end
                                        help2=help2+help;
                                    end
                                end
                                col=cint(ci)*(sum(ordnung)+m*n)+zaehler2;
                                gen_c=gen_c+1;
                                ROW(gen_c)=j;
                                COL(gen_c)=col;
                                OUTPUT(gen_c)=help2;
                            end
                        end
                    end
                    %Derivation with respect to z_(N-1)q
                    for q=1:m
                        for p2=1:n
                            zaehler2=zaehler2+1;
                            help2=0;
                            for cii=1:length(c)
                                if cint(cii)==cint(ci)
                                    help=0;
                                    for abl=0:ordnung(p2)-1
                                        help=help+DR_(abl+1,p2,cii)*h((cint(cii))+1)^(ordnung(p2)-abl)*polyval(psi(q,:,ordnung(p2)-(abl+1)+1),(c(cii)-x1((  cint(cii)  )+1))/h( cint(cii) +1));
                                    end
                                    help2=help2+help;
                                end
                            end
                            col=cint(ci)*(sum(ordnung)+m*n)+zaehler2;
                            gen_c=gen_c+1;
                            ROW(gen_c)=j;
                            COL(gen_c)=col;
                            OUTPUT(gen_c)=help2;
                        end
                    end
                end
                
                %-----------------------------------------------------------
                %Derivatives with respect to unknown parameters p_1 bis p_parameter
                if parameter>0
                    col=N*(n*m+sum(ordnung))+1:N*(n*m+sum(ordnung))+parameter;
                    gen_c_start=gen_c+1;
                    gen_c=gen_c+length(col);
                    ROW(gen_c_start:gen_c)=j;
                    COL(gen_c_start:gen_c)=col;
                    OUTPUT(gen_c_start:gen_c)=DpR(p1,1:parameter);
                end
            end
        end
        
        %----------End Boundary Conditions
        
        %----------Start J_i with C_i (notation according to manual)
        fz=zeros(n,max(ordnung)+1);
        for i=0:N-1
            for k=1:m
                for oi=0:max(ordnung)-1
                    for ni=1:n
                        if oi<ordnung(ni)
                            fz(ni,oi+1)=Poptimiert(oi,i,ni,tau((i)+1,k),m,x1,h,y,z,psival,ordnung(ni),k);
                        end
                    end
                end
                for ni=1:n
                    fz(ni,ordnung(ni)+1)=z(ni,k,(i)+1);
                end
                %p1 indicates, which component of g_ik will be differentiated with respect
                %to y und z
                %        for oi=1:max(ordnung)+1
                %            % Dg_(:,:,oi)=feval_problem(bf,'Dg',oi,[],tau((i)+1,k),fz,[],[],[],[],p,lambda,[],lambda_p);
                %            temp = feval_problem(bf,'jacobian',fz,[],[],[],tau((i)+1,k),p,lambda,[],lambda_p);
                %            temp=temp(:,:,oi);
                %            for ind1=1:size(temp,2)
                %                for ind2=1:size(temp,1)
                %                    Dg_(ind1,ind2,oi)=temp(ind2,ind1);
                %                end
                %            end
                %        end %%% Da war glaub ich ein Fehler, deshalb gleich:
                Dg_=feval_problem(problem,'jacobian',fz,[],[],[],tau((i)+1,k),p,lambda,[],lambda_p);
                % Dpg=feval_problem(bf,'Dpg',[],[],tau((i)+1,k),fz,[],[],[],[],p,lambda,[],lambda_p);
                Dpg=feval_problem(problem,'dP',fz,[],[],[],tau((i)+1,k),p,lambda,[],lambda_p);
                for p1=1:n
                    j=j+1;
                    %Differentiation with respect to y_i, sequence: y_i1 to
                    %y_i(max(ordnung))
                    zaehler=0;
                    for oi2=1:max(ordnung)
                        %p2 indicates, which component of y_i(oi2) will be
                        %differentiated
                        for p2=1:n
                            %Differentiation with respect to y_i(oi2) only if
                            %maximum order of the equation equals oi2
                            if oi2<=ordnung(p2)
                                zaehler=zaehler+1;
                                help=0;
                                for r1=1:oi2
                                    help=help+Dg_(p1,p2,r1)*(((tau((i)+1,k)-x1((i)+1))^(oi2-r1))/(faktorielle(oi2-r1+1)));
                                end
                                col=i*m*n+sum(ordnung)*i+zaehler;
                                gen_c=gen_c+1;
                                ROW(gen_c)=j;
                                COL(gen_c)=col;
                                OUTPUT(gen_c)=help;
                            end
                        end
                    end
                    zaehler=zaehler+1;
                    for p2=1:n
                        help(p2,1:m)=0;
                        for r1=1:ordnung(p2)
                            help(p2,1:m)=help(p2,1:m)+Dg_(p1,p2,r1)*h((i)+1)^(ordnung(p2)-r1+1).*psival(ordnung(p2)-r1+1,1:m,k);
                        end
                        help(p2,1:m)=help(p2,1:m)+Dg_(p1,p2,ordnung(p2)+1).*(k==(1:m));
                    end
                    help2=help(:).';
                    col=i*m*n+sum(ordnung)*i+zaehler:i*m*n+sum(ordnung)*i+zaehler+n*m-1;
                    gen_c_start=gen_c+1;
                    gen_c=gen_c+length(col);
                    ROW(gen_c_start:gen_c)=j;
                    COL(gen_c_start:gen_c)=col;
                    OUTPUT(gen_c_start:gen_c)=help2;
                    if parameter>0
                        col=N*(n*m+sum(ordnung))+1:N*(n*m+sum(ordnung))+parameter;
                        gen_c_start=gen_c+1;
                        gen_c=gen_c+length(col);
                        ROW(gen_c_start:gen_c)=j;
                        COL(gen_c_start:gen_c)=col;
                        OUTPUT(gen_c_start:gen_c)=Dpg(p1,1:parameter);
                    end
                end
            end
            if i<N-1
                for abl=0:max(ordnung)-1
                    for p1=1:n
                        if abl<=ordnung(p1)-1
                            j=j+1;
                            %Differentiation with respect to y_i(oi2)
                            zaehler=0;
                            for oi2=1:max(ordnung)
                                for p2=1:n
                                    if oi2<=ordnung(p2)
                                        zaehler=zaehler+1;
                                        if p1==p2
                                            if oi2-abl-1>=0
                                                col=i*(m*n+sum(ordnung))+zaehler;
                                                gen_c=gen_c+1;
                                                ROW(gen_c)=j;
                                                COL(gen_c)=col;
                                                OUTPUT(gen_c)=((x1((i+1)+1)-x1((i)+1))^(oi2-abl-1))/(faktorielle(oi2-abl));
                                            else
                                                col=i*(m*n+sum(ordnung))+zaehler;
                                                gen_c=gen_c+1;
                                                ROW(gen_c)=j;
                                                COL(gen_c)=col;
                                                OUTPUT(gen_c)=0;
                                            end
                                        end
                                    end
                                end
                            end
                            %D with respect to z_iq
                            zaehler=0;
                            for q=1:m
                                zaehler=zaehler+1;
                                col=i*m*n+(i+1)*sum(ordnung)+(zaehler-1)*n+p1;
                                gen_c=gen_c+1;
                                ROW(gen_c)=j;
                                COL(gen_c)=col;
                                OUTPUT(gen_c)=h((i)+1)^(ordnung(p1)-abl)*psival(ordnung(p1)-(abl+1)+1,q,m+1);
                            end
                            zaehler=0;
                            for oi2=1:max(ordnung)
                                for p2=1:n
                                    if oi2<=ordnung(p2)
                                        zaehler=zaehler+1;
                                        if p1==p2
                                            col=(i+1)*m*n+(i+1)*sum(ordnung)+zaehler;
                                            gen_c=gen_c+1;
                                            ROW(gen_c)=j;
                                            COL(gen_c)=col;
                                            OUTPUT(gen_c)=-(oi2==abl+1);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        dimsparse=N*(n*m+sum(ordnung))+parameter;
        ROW=ROW(1:gen_c);
        COL=COL(1:gen_c);
        OUTPUT=OUTPUT(1:gen_c);
        output=sparse(ROW,COL,OUTPUT,dimsparse,dimsparse);
        
        % In case of pathfollowing
        if length(a)>n_tmp
            
            if isempty(c)
                fza=zeros(n,max(ordnung));
                fzb=zeros(n,max(ordnung));
                for oi=1:max(ordnung)
                    for ni=1:n
                        if oi<=ordnung(ni)
                            fza(ni,oi)=y(ni,oi,(0)+1);
                            fzb(ni,oi)=Poptimiert(oi-1,N-1,ni,b,m,x1,h,y,z,psival,ordnung(ni),m+1);
                        end
                    end
                end
                % If the problem is posed on [0,Inf)
                if interval_feval(2) == Inf && interval_feval(1) == 0
                    ret1 = feval(problem,'path_dBV',[],fza,fzb,[],[],p,lambda,lambda_p);
                    ret2 = zeros(sum(ordnung)/2,1);
                    p_dB = [ret1 ; ret2];
                else
                    p_dB=feval(problem,'path_dBV',[],fza,fzb,[],[],p,lambda,lambda_p);
                end
                % If the problem is an Eigenvalue problem
                if ~isempty(lambda)
                    p_dB = [p_dB ; 0; 0];
                end
            else   % for additional conditions
                fzc=zeros(n,max(ordnung),length(c));
                for oi=1:max(ordnung)
                    for ni=1:n
                        if oi<=ordnung(ni)
                            for ci=1:length(c)
                                fzc(ni,oi,ci)=P(oi-1,cint(ci),ni,c(ci),m,x1,h,y,z,psi,ordnung(ni));
                            end
                        end
                    end
                end
                % If the problem is posed on [0,Inf)
                if interval_feval(2) == Inf && interval_feval(1) == 0
                    ret1 = feval(problem,'path_dBV',[],[],[],fzc,[],p,lambda,lambda_p);
                    ret2 = zeros(sum(ordnung)/2,1);
                    p_dB = [ret1 ; ret2];
                else
                    p_dB=feval(problem,'path_dBV',[],[],[],fzc,[],p,lambda,lambda_p);
                end
                % If the problem is an Eigenvalue problem
                if ~isempty(lambda)
                    p_dB = [p_dB ; 0; 0];
                end
            end
            output(1:sum(ordnung)+parameter,N*(sum(ordnung)+n*m)+parameter+1)=p_dB;
            
            j=sum(ordnung)+parameter+1;
            
            fz=zeros(n,max(ordnung)+1);
            for i=0:N-1
                for k=1:m
                    for oi=0:max(ordnung)-1
                        for ni=1:n
                            if oi<ordnung(ni)
                                fz(ni,oi+1)=Poptimiert(oi,i,ni,tau((i)+1,k),m,x1,h,y,z,psival,ordnung(ni),k);
                            end
                        end
                    end
                    for ni=1:n
                        fz(ni,ordnung(ni)+1)=z(ni,k,(i)+1);
                    end
                    s=size(fz);
                    if interval_feval(2) == Inf && interval_feval(1) == 0
                        ret1 = feval(problem,'path_jac',fz(1:s(1)/2,:),[],[],[],tau((i)+1,k),p,lambda,lambda_p);
                        ret2 = feval(problem,'path_jac',fz((s(1)/2+1):s(1),:)*(trafomatrix(1,1/tau((i)+1,k),s(2)))',[],[],[],1/tau((i)+1,k),p,lambda,lambda_p);
                        stmp=size(ret2,1);
                        A=trafomatrix(1,1/tau((i)+1,k),1);
                        for iter=1:stmp
                            ret2(iter)=reshape(ret2(iter,:,:),1,1)*A;
                        end
                        p_jac = zeros(2*stmp,1);
                        p_jac(1:stmp,1)=ret1;
                        p_jac(stmp+(1:stmp),1)=ret2;
                        %p_jac = [ret1 ; ret2];
                    else
                        p_jac=feval(problem,'path_jac',fz,[],[],[],tau((i)+1,k),p,lambda,lambda_p);
                        % If the problem is an Eigenvalue problem
                        if ~isempty(lambda)
                            p_jac = [p_jac ; 0];
                        end
                    end
                    output(j:j+n-1,N*(sum(ordnung)+n*m)+parameter+1)=p_jac;
                    j=j+n;
                end
                j=j+sum(ordnung);
            end
            
            output(n_tmp+1,n_tmp+1) = 1;
            if ~isempty(predictor) && ~predictor.lpfix
                output(N*(sum(ordnung)+n*m)+parameter+1,:)=predictor.tangent;
            end
        end
        
end

end

% Local functions

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

function ret=Poptimiert(abl,i,komp,x,m,x1,h,y,z,psival,ord,k) %Polynomial opitmized version
faktorielle=[1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800];

ret=sum(((x-x1((i)+1)).^((abl+1:ord)-abl-1))./(faktorielle((abl+1:ord)-abl)).*y(komp,abl+1:ord,(i)+1))+...
    sum(h((i)+1)^(ord-abl).*(z(komp,1:m,(i)+1).*psival(ord-(abl+1)+1,1:m,k)));
end

function [A] = trafomatrix(a,t,n)
n=n-1;
A=diag([1,(-1).^(1:n).*a.^(1:n)./t.^(2*(1:n))]);
if(n>1)
    A(1+(2:n),2) = a*((-1).^(2:n).*factorial(2:n)./t.^(1+(2:n)))';
    for i=4:n+1
        A(i,3:i-1) = -((A(i-1,2:i-2)*a)/t^2 + (A(i-1,3:i-1).*((i-1)+(3:i-1)))/t);
    end
end
end

