function [ der ] = higherchainrule( n )
%CHAINRULE computes the chainrule up to the n-th derivative
der = cell(n+1,1);

der{1}.A=0; %is the 0-th derivate, i.e. f(tau(t))
der{1}.b=0;
der{1}.c=1;

for j=2:n+1 %compute the (j-1)-th derivative
    [M,N]=size(der{j-1}.A);
    [K,L]=find(der{j-1}.A);
    A = zeros(length(K),j);
    b = zeros(length(K),1);
    c = zeros(length(K),1);
    %first the tau-derivatives from the chain rule
    for ii = 1:length(K)
        A(ii,:)=[der{j-1}.A(K(ii),:),0]; %copy the row to the next A-matrix
        A(ii,L(ii)) = A(ii,L(ii)) - 1;
        A(ii,L(ii)+1) = A(ii,L(ii)+1) + 1;
        b(ii)=der{j-1}.b(K(ii)); %the tau-derivatives corresponds to the same f-derivative
        c(ii)=der{j-1}.c(K(ii))*der{j-1}.A(K(ii),L(ii)); %the new coefficient is the old coefficient times the power of the tau-derivative
    end
    %now the f-derivatives (note: the inner derivative gives one tau')
    c = [c;der{j-1}.c];
    b = [b;der{j-1}.b+1];
    A = [A;[der{j-1}.A,zeros(M,1)]+[zeros(M,1),ones(M,1),zeros(M,N-1)]];
    
    %Relaxation
    [U,IA,IU]=unique([A,b],'rows');
    der{j}.c = zeros(length(IA),1);
    for k = 1:length(IU)
       der{j}.c(IU(k))=der{j}.c(IU(k))+c(k);
    end
    der{j}.A=U(:,1:end-1);
    der{j}.b=U(:,end);
end
end

