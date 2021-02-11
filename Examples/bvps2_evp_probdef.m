function [ret] = bvps2_evp_probdef(request,z,za,zb,zc,t,p,lambda,lambda_p)
% Problem definition for EVP problem in Manual

c = 3;
switch request
    case 'n'
        ret = 1;
    case 'orders'
        ret = [2];
    case 'problem'
        ret = [
            -z(1,3)+c/t^2*z(1,1)-lambda*z(1,1)
            ];
    case 'jacobian'
        %DON'T CHANGE THIS LINE:
        ret = zeros(length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1),0)),feval(mfilename,'n'),max(feval(mfilename,'orders'))+1);
        ret(1,1,1) = -lambda+c/t^2;
        ret(1,1,3) = -1;
    case 'interval'
        ret = [0,pi];
    case 'linear'
        ret = 1;
    case 'parameters'
        ret = 0;
    case 'c'
        ret = [];
    case 'BV'
        ret = [
            za(1,1);
            zb(1,1) ];
    case 'dBV'
        %DON'T CHANGE THIS LINE:
        ret = zeros(max(length(feval(mfilename,'c')),2-length(feval(mfilename,'c'))),length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1),0)),feval(mfilename,'n'),max(feval(mfilename,'orders')));
        ret(1,1,1,1) = 1;
        ret(2,2,1,1) = 1;
    case 'dP'
        ret = [];
    case 'dP_BV'
        ret = [];
    case 'initProfile'
        ret.initialMesh = [0 pi];
        ret.initialValues = [1 1];
        ret.lambda = 1;
    case 'EVP'
        ret = 1;
    case 'dLambda'
        ret = -z(1,1);
    case 'pathfollowing'
        ret.activate=0;
    otherwise
        ret = 0;
        
end

end