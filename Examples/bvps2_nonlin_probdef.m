function [ret] = bvps2_nonlin_probdef(request,z,za,zb,zc,t,p,lambda,lambda_p)
% Problem definition for nonlinear problem in Manual

nu=1/3;
mu=9;
gamma=1000;

switch request
    case 'n'
        ret= 2 ;
    case 'orders'
        ret=[ 2 2 ];
    case 'problem'
        ret=[
            z(1,3) + 3/t*z(1,2) + mu^2*z(2,1) + 2*gamma - z(1,1)*z(2,1);
            z(2,3) + 3/t*z(2,2) - mu^2*z(1,1) + 1/2*(z(1,1))^2];
    case 'jacobian'
        %DON'T CHANGE THIS LINE:
        ret = zeros(length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders'))+1);
        %Hier wurde das Format des return-Wertes festgelegt. Kann bei der Problemdefinition ignoriert werden.
        ret(1, 1, 1)= -z(2,1);
        ret(1, 1, 2)= 3/t;
        ret(1, 1, 3)= 1;
        ret(1, 2, 1)= mu^2 - z(1,1);
        ret(2, 1, 1)= -mu^2 + z(1,1);
        ret(2, 2, 2)= 3/t;
        ret(2, 2, 3)= 1;
    case 'interval'
        ret = [0,1];
    case 'linear'
        ret= 0;
    case 'parameters'
        ret=0;
    case 'c'
        ret = [];
    case 'BV'
        ret=[
            za(1,2);
            za(2,2);
            zb(1,1);
            zb(2,2) + (1-nu)*zb(2,1)];
    case 'dBV'
        %DON'T CHANGE THIS LINE:
        ret = zeros(max(length(feval(mfilename,'c')),2-length(feval(mfilename,'c'))),length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders')));
        %Hier wurde das Format des return-Wertes festgelegt. Kann bei der Problemdefinition ignoriert werden.
        ret(1,1,1,2)=1;
        ret(1,2,2,2)=1;
        ret(2,3,1,1)=1;
        ret(2,4,2,2)=1;
        ret(2,4,2,1)=(1-nu);
    case 'dP'
        ret=[];
    case 'dP_BV'
        ret=[];
    case 'initProfile'
        ret.initialMesh = linspace(0,1,6);
        ret.initialValues = ones(2,6);
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