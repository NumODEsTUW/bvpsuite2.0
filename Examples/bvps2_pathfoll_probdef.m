function [ret] = bvps2_pathfoll_probdef(request,z,za,zb,zc,t,p,lambda,lambda_p)
% Problem definition for pathfollowing problem in Manual

switch request
    case 'n'
        ret= 1;
    case 'orders'
        ret=[ 3 ];
        
    case 'problem'
        %Problemstellung als System der F(z(i,j))=0 (jede Zeile entspricht einer Gleichung) mit z(i,j), wobei
        %i..Lösungskomponente
        %j..Ableitung
        %also z(i,j)=z_i^(j-1)(t)
        ret= t^2*z(1,4)-t*z(1,3)+z(1,2) - t^3*p - lambda_p*( t*z(1,2)^2-t*z(1,1)*z(1,3)+z(1,1)*z(1,2) );
        
    case 'jacobian'
        %DON'T CHANGE THIS LINE:
        ret = zeros(length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1),1,0)),feval(mfilename,'n'),max(feval(mfilename,'orders'))+1);
        %ret(i,j,k): Ableitung von
        %i-ter Gleichung nach z(j,k)
        ret(1,1,1)= -lambda_p*( -t*z(1,3)+z(1,2) );
        ret(1,1,2)= 1 - lambda_p*( 2*t*z(1,2) + z(1,1) );
        ret(1,1,3)= -t - lambda_p*( -t*z(1,1) );
        ret(1,1,4)= t^2;
    case 'interval'
        ret = [0 1];
    case 'linear'
        ret= 0;
    case 'parameters'
        ret=1;
    case 'c'
        ret = [];
    case 'BV'
        ret=[za(1,1);
            za(1,2);
            zb(1,1);
            zb(1,2)-1];
        
    case 'dBV'
        %DON'T CHANGE THIS LINE:
        ret = zeros(max(length(feval(mfilename,'c')),2-length(feval(mfilename,'c'))),length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1),1,0)),feval(mfilename,'n'),max(feval(mfilename,'orders')));
        %ret(cInd,i,j,k)
        %cInd..spezifiziert Punkt, an dem Bedingung gestellt wird. (falls c=[] entspricht
        %1..a und 2..b, ansonsten Index von c)
        %i,j,k wie bei jacobian
        ret(1,1,1,1)=1;
        ret(1,2,1,2)=1;
        ret(2,3,1,1)=1;
        ret(2,4,1,2)=1;
    case 'dP'
        ret= -t^3 ;
    case 'dP_BV'
        ret=[0;0;0;0];
    case 'initProfile'
        ret.initialMesh = [0 1];
        ret.initialValues = [1 1];
        ret.parameters=1;
    case 'EVP'
        ret = 0;
    case 'dLambda'
        ret = 0;
    case 'pathfollowing'
        ret.activate=1;
        ret.pathdata=@PathCharData;
        ret.startat = 'start';
        ret.start=0;
        ret.steplength=1;
        ret.counter=Inf;
        ret.pit_stop=[100,-20];
        ret.require_exact=[0 8 100];
        ret.dir='C:\Users\Merlin\Desktop\MA_pathfollowing\Matlab\reffig2\';
        ret.name='bvps2_pathfoll_ref';
    case 'path_jac'
        ret=zeros(feval(mfilename,'n'),1);
        ret(1)= -( t*z(1,2)^2-t*z(1,1)*z(1,3)+z(1,1)*z(1,2) );
    case 'path_dBV'
        ret=zeros(sum(feval(mfilename,'orders'))+feval(mfilename,'parameters'),1);
    otherwise
        ret = 0;
        
end

end

function ret = PathCharData(x1,coeff,ordnung,rho)

ret= coeff(end-1);

end