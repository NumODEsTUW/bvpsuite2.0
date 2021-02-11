function [ret] = bvps2_lin_probdef(request,z,za,zb,zc,t,p,lambda,lambda_p)
% Problem definition for linear problem in Manual

epsilon=10e-4;

switch request
    case 'n'
        ret=  1;
    case 'orders'
        ret=[ 2 ];  
    case 'problem'
        %Problemstellung als System der F(z(i,j))=0 (jede Zeile entspricht einer Gleichung) mit z(i,j), wobei
        %i..Lösungskomponente
        %j..Ableitung
        %also z(i,j)=z_i^(j-1)(t)
        ret=[
            epsilon*z(1,3)+z(1,2)-(1+epsilon)*z(1,1)
            ];
    case 'jacobian'
        %DON'T CHANGE THIS LINE:
        ret = zeros(length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders'))+1);
        %ret(i,j,k): Ableitung von 
        %i-ter Gleichung nach z(j,k)
         ret(1, 1, 1)= -(1+epsilon);
         ret(1, 1, 2)=1;
         ret(1, 1, 3)=epsilon;
    case 'interval'
        ret = [-1 ,1 ];  
    case 'linear'
        ret= 1; 
    case 'parameters'
        ret=0;
    case 'c'
        ret = [];
    case 'BV'
        ret=[za(1,1)-(1+exp(-2));
            zb(1,1)-(1+exp((-2*(1+epsilon))/epsilon))
            ];
    case 'dBV'
        %DON'T CHANGE THIS LINE:
        ret = zeros(max(length(feval(mfilename,'c')),2-length(feval(mfilename,'c'))),length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders')));
        %ret(cInd,i,j,k)
        %cInd..spezifiziert Punkt, an dem Bedingung gestellt wird. (falls c=[] entspricht
        %1..a und 2..b, ansonsten Index von c)
        %i,j,k wie bei jacobian              
         ret( 1,1 ,1 ,1 )= 1;
         ret( 2,2,1,1)=1;
    case 'dP'
        ret=[];
    case 'dP_BV'
        ret=[];
    case 'initProfile'
        ret.initialMesh = [ 1];
        ret.initialValues = [ 1];
    case 'EVP'
        ret = 0;
    case 'dLambda'
        ret = 0;
    case 'pathfollowing'
        ret.activate=0;
    otherwise
        ret = 0;
                
end

end