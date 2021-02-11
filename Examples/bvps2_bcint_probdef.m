function [ret] = bvps2_bcint_probdef(request,z,za,zb,zc,t,p,lambda,lambda_p)

switch request
    case 'n'
        ret= 3 ;
    case 'orders'
        ret=[ 2 3 2 ];
    case 'problem'
        %Problemstellung als System der F(z(i,j))=0 (jede Zeile entspricht einer Gleichung) mit z(i,j), wobei
        %i..L�sungskomponente
        %j..Ableitung
        %also z(i,j)=z_i^(j-1)(t)
        
        C1=4e7;
        C2=5e-12;
        C3=-1.6e-3;
        C4=5.7e-2;
        C5=-5e2;
        C6=1e-15;
        
        ret=[
            C1*(z(1,1)-z(2,1))-C2*z(1,3)+z(3,2)*C3;
            C1*(z(1,2)-z(2,2))+C2*z(2,4)+C4*z(3,3)-C5;
            C6*(C1*(z(1,2)-z(2,2))+C2*z(2,4))+z(3,1)
            ];
    case 'jacobian'
        %DON'T CHANGE THIS LINE:
        ret = zeros(length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders'))+1);
        %ret(i,j,k): Ableitung von
        %i-ter Gleichung nach z(j,k)
        
        C1=4e7;
        C2=5e-12;
        C3=-1.6e-3;
        C4=5.7e-2;
        C6=1e-15;
        
        ret(1, 1, 1)= C1;
        ret(1, 1, 3)= -C2;
        ret(2, 1, 2)= C1;
        ret(3, 1, 2)= C6*C1;
        ret(1, 2, 1)= -C1;
        ret(2, 2, 2)= -C1;
        ret(2, 2, 4)= C2;
        ret(3, 2, 2)= -C6*C1;
        ret(3, 2, 4)= C6*C2;
        ret(1, 3, 2)= C3;
        ret(2, 3, 3)= C4;
        ret(3, 3, 1)= 1;
    case 'interval'
        D=1e-5;
        ret = [-D/2,D/2];
    case 'linear'
        ret=1;
    case 'parameters'
        ret=0;
    case 'c'
        D=1e-5;
        ret = [-D/2,0,D/2];
    case 'BV'
        B1=2e-1; B2=1.5e-1;
        
        ret=[
            zc(1,1,1)+B1;
            zc(1,1,3)-B1;
            zc(2,1,1)+B2;
            zc(2,1,3)-B2;
            zc(2,3,2);
            zc(3,1,1);
            zc(3,1,3);
            ];
    case 'dBV'
        %DON'T CHANGE THIS LINE:
        ret = zeros(max(length(feval(mfilename,'c')),2-length(feval(mfilename,'c'))),length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders')));
        %ret(cInd,i,j,k)
        %cInd..spezifiziert Punkt, an dem Bedingung gestellt wird. (falls c=[] entspricht
        %1..a und 2..b, ansonsten Index von c)
        %i,j,k wie bei jacobian
        ret( 1,1 ,1 ,1 )= 1;
        ret( 3,2 ,1 ,1 )= 1;
        ret( 1,3 ,2 ,1 )= 1;
        ret( 3,4 ,2 ,1 )= 1;
        ret( 2,5 ,2 ,3 )= 1;
        ret( 1,6 ,3 ,1 )= 1;
        ret( 3,7 ,3 ,1 )= 1;
    case 'dP'
        ret=[];
    case 'dP_BV'
        ret=[];
    case 'initProfile'
        ret.initialMesh=[];
        ret.initialValues=[];
    case 'EVP'
        ret = 0;
    case 'dLambda'
        ret = 0;
    case 'pathfollowing'
        ret.activate=0;
    otherwise
        ret = 0;
        
end