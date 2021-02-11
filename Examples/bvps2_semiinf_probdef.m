function [ret] = bvps2_semiinf_probdef(request,z,za,zb,zc,t,p,lambda,lambda_p)
% Problem definition for problem posed on a semi-infinite interval in Manual

switch request
    case 'n'
        ret= 1;
    case 'orders'
        ret=[ 2 ];  
    case 'problem'
        %Problemstellung als System der F(z(i,j))=0 (jede Zeile entspricht einer Gleichung) mit z(i,j), wobei
        %i..Lösungskomponente
        %j..Ableitung
        %also z(i,j)=z_i^(j-1)(t)
        ret=[
            z(1,3) + 2/t*z(1,2) - 4*(z(1,1)+1)*z(1,1)*(z(1,1)-0.1)
            ];
    case 'jacobian'
        %DON'T CHANGE THIS LINE:
        ret = zeros(length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders'))+1);
        %ret(i,j,k): Ableitung von 
        %i-ter Gleichung nach z(j,k)
         ret(1, 1, 1)= -4*z(1,1)*(z(1,1)-0.1) - 4*(z(1,1)+1)*(z(1,1)-0.1) - 4*(z(1,1)+1)*z(1,1);
         ret(1, 1, 2)= 2/t;
         ret(1, 1, 3)= 1;
    case 'interval'
        ret = [0,Inf];  
    case 'linear'
        ret= 0; 
    case 'parameters'
        ret=0;
    case 'c'
        ret = [];
    case 'BV'
        ret=[
            za(1,2);
            zb(1,1)-0.1];
    case 'dBV'
        %DON'T CHANGE THIS LINE:
        ret = zeros(max(length(feval(mfilename,'c')),2-length(feval(mfilename,'c'))),length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders')));
        %ret(cInd,i,j,k)
        %cInd..spezifiziert Punkt, an dem Bedingung gestellt wird. (falls c=[] entspricht
        %1..a und 2..b, ansonsten Index von c)
        %i,j,k wie bei jacobian              
         ret( 1,1 ,1 ,2 )= 1;
         ret( 2,2 ,1 ,1) = 1;
    case 'dP'
        ret=[];
    case 'dP_BV'
        ret=[];
    case 'initProfile'
        ret.initialMesh = [ 0.0225    0.1000    0.1775    0.2000    0.2225    0.3000    0.3775    0.4000    0.4225    0.5000    0.5775  0.6000    0.6225    0.7000    0.7775    0.8000    0.8225    0.9000    0.9775    1.0000    1.0231    1.1111      1.2157    1.2500    1.2862    1.4286    1.6063    1.6667    1.7317    2.0000    2.3666    2.5000    2.6493  3.3333    4.4936    5.0000    5.6351   10.0000 45];
        ret.initialValues = [-0.3042   -0.3037   -0.3024   -0.3020   -0.3014   -0.2991   -0.2962   -0.2952   -0.2942   -0.2902   -0.2857  -0.2842   -0.2828   -0.2773   -0.2713   -0.2694   -0.2675   -0.2607   -0.2535   -0.2513   -0.2490   -0.2400  -0.2286   -0.2248   -0.2207   -0.2041   -0.1825   -0.1750   -0.1669   -0.1336   -0.0899   -0.0751   -0.0592  0.0007    0.0593    0.0729    0.0834    0.0994 0.1];
    case 'EVP'
        ret=0;
    case 'dLambda'
        ret = 0;
    case 'pathfollowing'
        ret.activate=0;
    otherwise
        ret = 0;
        
end

end