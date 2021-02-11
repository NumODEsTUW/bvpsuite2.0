function [ret] = feval_problem(problem,request,varargin) %problem,request,z,za,zb,zc,t,p,lambda,EVPtrafo,lambda_p)

interval=feval(problem,'interval');
% if isempty(interval)
%     switch request
%         case 'interval'
%             ret = 0;
%             return
%         otherwise
%             ret=[0,0];
%     end
% end
EVPtrafo=1;
if(nargin>=10 && ~isempty(varargin{8}))
    EVPtrafo = varargin{8};
end
if(~strcmp(request,'n') && ~strcmp(request,'orders') && ~strcmp(request,'parameters'))
    z=zeros(feval_problem(problem,'n',[],[],[],[],[],[],[],EVPtrafo,[]),max(feval_problem(problem,'orders',[],[],[],[],[],[],[],EVPtrafo,[])));
    za=z;
    zb=z;
    zc=[];
    p=zeros(feval_problem(problem,'parameters',[],[],[],[],[],[],[],EVPtrafo,[]),1);
end
t=0;
lambda=0;
lambda_p=0;

if(nargin>=3 && ~isempty(varargin{1}))
    z = varargin{1};
end
if(nargin>=4 && ~isempty(varargin{2}))
    za = varargin{2};
end
if(nargin>=5 && ~isempty(varargin{3}))
    zb = varargin{3};
end
if(nargin>=6 && ~isempty(varargin{4}))
    zc = varargin{4};
end
if(nargin>=7 && ~isempty(varargin{5}))
    t = varargin{5};
end
if(nargin>=8 && ~isempty(varargin{6}))
    p = varargin{6};
end
if(nargin>=9 && ~isempty(varargin{7}))
    lambda = varargin{7};
end
if(nargin>=11 && ~isempty(varargin{9}))
    lambda_p = varargin{9};
end

if(interval(2)<Inf)
    if(feval(problem,'EVP')>0 && EVPtrafo)
        switch request
            case 'n'
                ret=feval(problem,request)+1;
            case 'orders'
                ret=feval(problem,request);
                ret=[ret 1];
            case 'problem'
                ret=feval(problem,request,z,[],[],[],t,p(1:end-1),p(end),lambda_p);
                ret=[ret; z(end,2)-norm(z(1:end-1,1),2)^2];
            case 'jacobian'
                ret_old=feval(problem,request,z,[],[],[],t,p(1:end-1),p(end),lambda_p);
                [s(1),s(2),s(3)]=size(ret_old);
                ret = zeros(s(1)+1,s(2)+1,max(2,s(3)));
                ret(1:s(1),1:s(2),1:s(3))=ret_old;
                ret(s(1)+1,s(2)+1,2)=1;
                for i=1:s(2)
                    ret(s(1)+1,i,1)=-2*z(i,1);
                end
            case 'interval'
                ret = feval(problem,request);
            case 'linear'
                ret = 0;
            case 'parameters'
                ret=feval(problem,request)+1;
            case 'c'
                ret = feval(problem,request);
                if(~isempty(ret))
                    ret = sort([interval(1) ret interval(2)]);
                end
            case 'BV'
                ret=feval(problem,request,[],za,zb,zc,[],p(1:end-1),p(end),lambda_p);
                ret=[ret; za(end,1); zb(end,1)-1];
            case 'dBV'
                ret_old=feval(problem,request,[],za,zb,zc,[],p(1:end-1),p(end),lambda_p);
                [s(1),s(2),s(3),s(4)]=size(ret_old);
                if(isempty(feval(problem,'c')))
                    ret = zeros(s(1),s(2)+2,s(3)+1,s(4));
                    ret(1:s(1),1:s(2),1:s(3),1:s(4))=ret_old;
                    ret(1,s(2)+1,s(3)+1,1)=1;
                    ret(s(1),s(2)+2,s(3)+1,1)=1;
                else
                    ret = zeros(s(1)+2,s(2)+2,s(3)+1,s(4));
                    ret(2:s(1)+1,1:s(2),1:s(3),1:s(4))=ret_old;
                    ret(1,s(2)+1,s(3)+1,1)=1;
                    ret(s(1)+2,s(2)+2,s(3)+1,1)=1;
                end
            case 'dP'
                ret=feval(problem,request,z,[],[],[],t,p(1:end-1),p(end),lambda_p);
                dLambda=feval(problem,'dLambda',z,[],[],[],t,p(1:end-1),p(end),lambda_p);
                ret = [ret,dLambda];
                ret = [ret;0];
            case 'dP_BV'
                ret=zeros(2+length(feval(problem,'BV',[],za,zb,zc,[],p(1:end-1),p(end),lambda_p)),1+feval(problem,'parameters'));
                ret(1:end-2,1:end-1)=feval(problem,request,[],za,zb,zc,[],p(1:end-1),p(end),lambda_p);
            case 'initProfile'
                ret=feval(problem,request);
            otherwise
                ret=feval(problem,request,z,za,zb,zc,t,p(1:end-1),p(end),lambda_p);
        end
    else
        if(strcmp(request,'n') || strcmp(request,'orders') || strcmp(request,'parameters'))
            ret = feval(problem,request);
        else
            ret = feval(problem,request,z,za,zb,zc,t,p,lambda,lambda_p);
        end
    end
else
    if(feval(problem,'EVP')>0 && EVPtrafo)
        switch request
            case 'n'
                ret=feval(problem,request)+1;
                if(interval(1) == 0 && interval(2) == Inf)
                    ret=2*ret;
                end
            case 'orders'
                ret=feval(problem,request);
                ret=[ret 1];
                if(interval(1) == 0 && interval(2) == Inf)
                    ret=[ret ret];
                end
            case 'problem'
                s=size(z);
                if(interval(1) == 0 && interval(2) == Inf)
                    dimz = s(1)/2-1;
                    ret=[feval(problem,request,z(1:dimz,:),[],[],[],t,p(1:end-1),p(end),lambda_p);...
                        z(dimz+1,2)-norm(z(1:dimz,1),2)^2;...
                        feval(problem,request,z((dimz+2):(s(1)-1),:)*(trafomatrix(1,1/t,s(2)))',[],[],[],1/t,p(1:end-1),p(end),lambda_p);...
                        z(end,2)*(-t^2)-norm(z(dimz+2:(s(1)-1),1),2)^2];
                else
                    zA=z*(trafomatrix(interval(1),interval(1)/t,s(2)))';
                    ret=feval(problem,request,zA,[],[],[],interval(1)/t,p(1:end-1),p(end),lambda_p);
                    ret=[ret; (-t^2/interval(1))*z(end,2)-norm(z(1:end-1,1),2)^2];
                end
            case 'jacobian'
                if(interval(1) == 0 && interval(2) == Inf)
                    s = size(z);
                    dimz = s(1)/2-1;
                    ret1=feval(problem,request,z(1:dimz,:),[],[],[],t,p(1:end-1),p(end),lambda_p);
                    ret2=feval(problem,request,z((dimz+2):(s(1)-1),:)*(trafomatrix(1,1/t,s(2)))',[],[],[],1/t,p(1:end-1),p(end),lambda_p);
                    [s(1),s(2),s(3)]=size(ret2);
                    A=trafomatrix(1,1/t,s(3));
                    for i=1:s(1)
                        ret2(i,:,:)=reshape(ret2(i,:,:),s(2),s(3))*A;
                    end
                    ret = zeros(2*s(1)+2,2*s(2)+2,s(3));
                    
                    ret(1:s(1),1:s(2),1:s(3))=ret1; %Die ursprünglichen Gleichungen
                    
                    ret(s(1)+1,s(2)+1,2)=1; %Die Normierungsbed auf [0,1]
                    for i=1:s(2)
                        ret(s(1)+1,i,1)=-2*z(i,1);
                    end
                    
                    ret(s(1)+1+(1:s(1)),s(2)+1+(1:s(2)),1:s(3))=ret2; %Die transf Glgen
                    
                    ret(2*s(1)+2,2*s(2)+2,2)=-t^2; %Die Normierungsbed transf von [1,Infty)
                    for i=s(2)+1+(1:s(2))
                        ret(2*s(1)+2,i,1)=-2*z(i,1);
                    end
                else
                    s = size(z);
                    A=trafomatrix(interval(1),interval(1)/t,s(2));
                    zA=z*(trafomatrix(interval(1),interval(1)/t,s(2)))';
                    ret_old=feval(problem,request,zA,[],[],[],interval(1)/t,p(1:end-1),p(end),lambda_p);
                    [s(1),s(2),s(3)]=size(ret_old);
                    ret = zeros(s(1)+1,s(2)+1,max(2,s(3)));
                    ret(1:s(1),1:s(2),1:s(3))=ret_old;
                    ret(s(1)+1,s(2)+1,2)=1;
                    for i=1:s(2)
                        ret(s(1)+1,i,1)=-2*zA(i,1);
                    end
                    [s(1),s(2),s(3)]=size(ret);
                    for i=1:s(1)
                        ret(i,:,:)=reshape(ret(i,:,:),s(2),s(3))*A;
                    end
                end
            case 'interval'
                ret = [0,1];
            case 'linear'
                ret = 0;
            case 'parameters'
                ret=feval(problem,request)+1;
            case 'c'
                ret = feval(problem,request);
                if(~isempty(ret))
                    ret = sort([interval(1) ret interval(2)]);
                end
            case 'BV'
                s=size(za);
                if(interval(1) == 0 && interval(2) == Inf)
                    orders = feval_problem(problem,'orders');
                    n=s(1)/2-1;
                    ret=[feval(problem,request,[],za(1:n,:),[za(n+2:(s(1)-1),1),zeros(n,s(2)-1)],[],[],p(1:end-1),p(end),lambda_p);...
                        za(n+1,1); za(end,1)-1];%zc hasn't been worked out yet
                    ret2 = zeros(sum(orders)/2,1);
                    for i=1:n+1
                        A=trafomatrix(1,1,orders(i));
                        ret2(sum(orders(1:i-1))+(1:orders(i)))=(zb(i,1:orders(i))-zb(i+n+1,1:orders(i))*A')';
                    end
                    ret=[ret;ret2];
                else
                    A=trafomatrix(interval(1),interval(1),s(2));
                    zaneu=[zb(:,1),zeros(s(1),s(2)-1)];
                    zbneu=za*A';
                    ret=feval(problem,request,z,zaneu,zbneu,[],[],p(1:end-1),p(end),lambda_p);
                    ret=[ret; zaneu(end,1); zbneu(end,1)-1];
                end
            case 'dBV'
                s=size(za);
                if(interval(1) == 0 && interval(2) == Inf)
                    orders = feval_problem(problem,'orders');
                    n=s(1)/2-1;
                    ret1=feval(problem,request,[],za(1:n,:),[za(n+2:(s(1)-1),1),zeros(n,s(2)-1)],[],[],p(1:end-1),p(end),lambda_p);
                    [s(1),s(2),s(3),s(4)]=size(ret1);
                    ret = zeros(2,s(2)+2+sum(orders)/2,2*n+2,s(4));
                    for i=1:s(2) %transform original BCs
                        ret(1,i,1:n,:) = ret1(1,i,:,:);
                        ret(1,i,(n+2):(2*n+1),1) = ret1(2,i,:,1);
                    end
                    ret(1,s(2)+1,n+1,1)=1;
                    ret(1,s(2)+2,2*(n+1),1)=1;
                    
                    for i=1:n+1 %smoothness conditions
                        A=trafomatrix(1,1,orders(i));
                        for j=1:orders(i)
                            ret(2,s(2)+2+sum(orders(1:i-1))+j,i,j) = 1;
                            ret(2,s(2)+2+sum(orders(1:i-1))+j,n+1+i,:) = -A(j,:);
                        end
                    end
                else
                    A=trafomatrix(interval(1),interval(1),s(2));
                    ret1=feval(problem,request,[],[zb(:,1),zeros(s(1),s(2)-1)],za*A',[],[],p(1:end-1),p(end),lambda_p);
                    [s(1),s(2),s(3),s(4)]=size(ret1);
                    ret2 = zeros(s(1),s(2)+2,s(3)+1,s(4));
                    ret2(1:s(1),1:s(2),1:s(3),1:s(4))=ret1;
                    ret2(1,s(2)+1,s(3)+1,1)=1;
                    ret2(s(1),s(2)+2,s(3)+1,1)=1;
                    [s(1),s(2),s(3),s(4)]=size(ret2);
                    ret=zeros(s);
                    A=trafomatrix(interval(1),interval(1),s(4));
                    for i=1:s(2)
                        ret(1,i,:,:)=reshape(ret2(2,i,:,:),s(3),s(4))*A';
                    end
                    A=trafomatrix(interval(1),interval(2),s(4));
                    for i=1:s(2)
                        ret(2,i,:,:)=reshape(ret2(1,i,:,:),s(3),s(4))*A';
                    end
                end
            case 'dP'
                s = size(z);
                if(interval(1) == 0 && interval(2) == Inf)
                    n=s(1)/2-1;
                    ret1=[feval(problem,request,z(1:n,:),[],[],[],t,p(1:end-1),p(end),lambda_p),...
                        feval(problem,'dLambda',z(1:n,:),[],[],[],t,p(1:end-1),p(end),lambda_p)];
                    A=trafomatrix(1,1/t,s(2));
                    zA=z((n+2):(2*n+1),:)*A';
                    ret2=[feval(problem,request,zA,[],[],[],1/t,p(1:end-1),p(end),lambda_p),...
                        feval(problem,'dLambda',zA,[],[],[],1/t,p(1:end-1),p(end),lambda_p)];
                    ret=[ret1;0;ret2;0];
                else
                    A=trafomatrix(interval(1),interval(1)/t,s(2));
                    zA=z*A';
                    ret_old=feval(problem,request,zA,[],[],[],interval(1)/t,p(1:end-1),p(end),lambda_p);
                    dLambda=feval(problem,'dLambda',zA,[],[],[],interval(1)/t,p(1:end-1),p(end),lambda_p);
                    ret_old = [ret_old,dLambda];
                    ret = [ret_old;0];
                    %                      [s(1),s(2),s(3)]=size(ret);
                    %                      A=trafomatrix(interval(1),interval(1)/t,s(2));
                    %                      for i=1:s(1)
                    %                          ret(i,:,:)=reshape(ret(i,:,:),s(2),s(3))*A;
                    %                      end
                end
            case 'dP_BV'
                s=size(za);
                if(interval(1) == 0 && interval(2) == Inf)
                    n=s(1)/2-1;
                    ret1 = feval(problem,'BV',[],za(1:n,:),[za(n+2:(s(1)-1),1),zeros(n,s(2)-1)],[],[],p(1:end-1),p(end),lambda_p);
                    [sBV(1),sBV(2),sBV(3),sBV(4)]=size(ret1);
                    ret = zeros(sum(feval_problem(problem,'orders'))+1+feval(problem,'parameters'),1+feval(problem,'parameters'));
                    retold = feval(problem,request,[],za(1:n,:),[za(n+2:(s(1)-1),1),zeros(n,s(2)-1)],[],[],p(1:end-1),p(end),lambda_p);
                    [sdP(1),sdP(2)]=size(retold);
                    ret(1:sdP(1),1:sdP(2))=retold;
                else
                    A=trafomatrix(interval(1),interval(1),s(2));
                    ret1=feval(problem,request,[],[zb(:,1),zeros(s(1),s(2)-1)],za*A',[],[],p(1:end-1),p(end),lambda_p);
                    [sdP(1),sdP(2)]=size(ret1);
                    ret = zeros(sdP(1)+2,1+feval(problem,'parameters'));
                    ret(1:sdP(1),1:sdP(2))=ret1;
                end
            case 'initProfile'
                ret=feval(problem,request);
            otherwise
                ret=feval(problem,request,z,za,zb,zc,t,p(1:end-1),p(end),lambda_p);
        end
    else
        switch request
            case 'n'
                ret=feval(problem,request);
                if(interval(1) == 0 && interval(2) == Inf)
                    ret=2*ret;
                end
            case 'orders'
                ret=feval(problem,request);
                if(interval(1) == 0 && interval(2) == Inf)
                    ret=[ret ret];
                end
            case 'problem'
                s=size(z);
                if(interval(1) == 0 && interval(2) == Inf)
                    ret=[feval(problem,request,z(1:s(1)/2,:),[],[],[],t,p,lambda,lambda_p);...
                        feval(problem,request,z((s(1)/2+1):s(1),:)*(trafomatrix(1,1/t,s(2)))',[],[],[],1/t,p,lambda,lambda_p)];
                else
                    ret=feval(problem,request,z*(trafomatrix(interval(1),interval(1)/t,s(2)))',[],[],[],interval(1)/t,p,lambda,lambda_p);
                end
            case 'jacobian'
                s = size(z);
                if(interval(1) == 0 && interval(2) == Inf)
                    ret1=feval(problem,request,z(1:s(1)/2,:),[],[],[],t,p,lambda,lambda_p);
                    ret2=feval(problem,request,z((s(1)/2+1):s(1),:)*(trafomatrix(1,1/t,s(2)))',[],[],[],1/t,p,lambda,lambda_p);
                    [s(1),s(2),s(3)]=size(ret2);
                    A=trafomatrix(1,1/t,s(3));
                    for i=1:s(1)
                        ret2(i,:,:)=reshape(ret2(i,:,:),s(2),s(3))*A;
                    end
                    ret = zeros(2*s(1),2*s(2),s(3));
                    ret(1:s(1),1:s(2),1:s(3))=ret1;
                    ret(s(1)+(1:s(1)),s(2)+(1:s(2)),1:s(3))=ret2;
                else
                    A=trafomatrix(interval(1),interval(1)/t,s(2));
                    ret=feval(problem,request,z*A',[],[],[],interval(1)/t,p,lambda,lambda_p);
                    [s(1),s(2),s(3)]=size(ret);
                    for i=1:s(1)
                        ret(i,:,:)=reshape(ret(i,:,:),s(2),s(3))*A;
                    end
                end
            case 'interval'
                ret = [0,1];
            case 'linear'
                ret = feval(problem,'linear');
            case 'parameters'
                ret= feval(problem,'parameters');%0;
            case 'c'
                ret = [];
            case 'BV'
                s=size(za);
                if(interval(1) == 0 && interval(2) == Inf)
                    orders = feval_problem(problem,'orders',[],[],[],[],[],[],[],EVPtrafo,lambda_p);
                    n=s(1)/2;
                    ret1=feval(problem,request,[],za(1:n,:),[za(n+1:2*n,1),zeros(n,s(2)-1)],[],[],p,lambda,lambda_p); %zc,p hasn't been worked out yet
                    ret2 = zeros(sum(orders)/2,1);
                    for i=1:n
                        A=trafomatrix(1,1,orders(i));
                        ret2(sum(orders(1:i-1))+(1:orders(i)))=(zb(i,1:orders(i))-zb(i+n,1:orders(i))*A')';
                    end
                    ret=[ret1;ret2];
                else
                    A=trafomatrix(interval(1),interval(1),s(2));
                    ret=feval(problem,request,z,[zb(:,1),zeros(s(1),s(2)-1)],za*A',[],[],p,lambda,lambda_p);
                end
            case 'dBV'
                if(isempty(feval(problem,'c')))
                    orders = feval_problem(problem,'orders',[],[],[],[],[],[],[],EVPtrafo,lambda_p);
                    s=size(za);
                    if(interval(1) == 0 && interval(2) == Inf)
                        n=s(1)/2;
                        ret1=feval(problem,request,[],za(1:n,:),[za(n+1:2*n,1),zeros(n,s(2)-1)],[],[],p,lambda,lambda_p);
                        [s(1),s(2),s(3),s(4)]=size(ret1);
                        ret = zeros(2,sum(orders),2*n,s(4));
                        for i=1:s(2) %transform original BCs
                            ret(1,i,1:n,:) = ret1(1,i,:,:);
                            ret(1,i,n+1:2*n,1) = ret1(2,i,:,1);
                        end
                        
                        for i=1:n %smoothness conditions
                            A=trafomatrix(1,1,orders(i));
                            for j=1:orders(i)
                                ret(2,s(2)+sum(orders(1:i-1))+j,i,j) = 1;
                                ret(2,s(2)+sum(orders(1:i-1))+j,n+i,:) = -A(j,:);
                            end
                        end
                        
                    else
                        A=trafomatrix(interval(1),interval(1),s(2));
                        ret1=feval(problem,request,[],[zb(:,1),zeros(s(1),s(2)-1)],za*A',[],[],p,lambda,lambda_p);
                        [s(1),s(2),s(3),s(4)]=size(ret1);
                        ret=zeros(s);
                        A=trafomatrix(interval(1),interval(1),s(4));
                        for i=1:s(2)
                            ret(1,i,:,:)=reshape(ret1(2,i,:,:),s(3),s(4))*A;
                        end
                        A=trafomatrix(interval(1),interval(2),s(4));
                        for i=1:s(2)
                            ret(2,i,:,:)=reshape(ret1(1,i,:,:),s(3),s(4))*A';
                        end
                    end
                end
            case 'initProfile'
                ret = feval(problem,request);
            case 'dP'
                s=size(z);
                if(interval(1) == 0 && interval(2) == Inf)
                    ret=[feval(problem,request,z(1:s(1)/2,:),[],[],[],t,p,lambda,lambda_p);...
                        feval(problem,request,z((s(1)/2+1):s(1),:)*(trafomatrix(1,1/t,s(2)))',[],[],[],1/t,p,lambda,lambda_p)];
                else
                    ret=feval(problem,request,z*(trafomatrix(interval(1),interval(1)/t,s(2)))',[],[],[],interval(1)/t,p,lambda,lambda_p);
                end
            case 'dP_BV'
                ret = zeros( 2*sum(feval(problem,'orders'))+feval(problem,'parameters'),1 );
                
            otherwise
                ret=feval(problem,request,z,za,zb,zc,t,p,lambda,lambda_p);
        end
    end
end

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

