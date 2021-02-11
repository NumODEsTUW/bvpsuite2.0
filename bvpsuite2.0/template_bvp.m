function [ret] = template_bvp(request,z,za,zb,zc,t,p,lambda,lambda_p)
switch request
    case 'n'
        ret=1;
    case 'orders'
        ret=[3];  
    case 'problem'
        %Problemstellung als System der F(z(i,j))=0 (jede Zeile entspricht einer Gleichung) mit z(i,j), wobei
        %i..Lösungskomponente
        %j..Ableitung
        %also z(i,j)=z_i^(j)(t)
        ret=[z(1,1)+z(1,3)-z(1,2)-z(1,4)];
    case 'jacobian'
        %DON'T CHANGE THIS LINE:
        ret = zeros(length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1),0)),feval(mfilename,'n'),max(feval(mfilename,'orders'))+1);
        %ret(i,j,k): Ableitung von 
        %i-ter Gleichung nach z(j,k)
         ret(1,1,1)=1;
         ret(1,1,2)=-1;
         ret(1,1,3)=1;
         ret(1,1,4)=-1;
    case 'interval'
        ret = [-1,1];  
    case 'linear'
        ret=1; 
    case 'parameters'
        ret=0;
    case 'c'
        ret = [];
    case 'BV'
        ret=[za(1,1)-exp(-1);
             za(1,2)-exp(-1);
             za(1,3)-exp(-1)];
    case 'dBV'
        %DON'T CHANGE THIS LINE:
        ret = zeros(max(length(feval(mfilename,'c')),2-length(feval(mfilename,'c'))),length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1),0)),feval(mfilename,'n'),max(feval(mfilename,'orders')));
        %ret(cInd,i,j,k)
        %cInd..spezifiziert Punkt, an dem Bedingung gestellt wird. (falls c=[] entspricht
        %1..a und 2..b, ansonsten Index von c)
        %i,j,k wie bei jacobian              
         ret(1,1,1,1)=1;
         ret(1,2,1,2)=1;
         ret(1,3,1,3)=1;
    case 'dP'
        ret=[0];
    case 'dP_BV'
        ret=[0;0;0];
    case 'initProfile'
        ret.initialMesh = linspace(-1,1,20);
        ret.initialValues = linspace(-1,1,20).^2;
    case 'EVP'
        ret=0;
    case 'dLambda'
        ret = 0;
    case 'pathfollowing'
        ret.activate=0;
    case 'path_jac'
        ret=zeros(feval(mfilename,'n'),1);
    case 'path_dBV'
        ret=zeros(sum(feval(mfilename,'orders'))+feval(mfilename,'parameters'),1);
    otherwise
        ret = 0;
                
end

end

