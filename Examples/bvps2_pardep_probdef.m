function [ret] = bvps2_pardep_probdef(request,z,za,zb,zc,t,p,lambda,lambda_p)
% Problem definition for parameter dependent problem in Manual

switch request
    case 'n'
        ret = 2;
    case 'orders'
        ret = [ 1 1 ];
    case 'problem'
        ret = [
            (1-p)*z(1,2)+p*(sqrt(t)-z(1,1))/(t-z(1,1));
            z(2,2)-z(1,1)
            ];
    case 'jacobian'
        %DON'T CHANGE THIS LINE:
        ret = zeros(length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders'))+1);
        ret(1, 1, 1) = p*(sqrt(t)-t)/(t-z(1,1))^2;
        ret(1, 1, 2) = 1-p;
        ret(2, 1, 1) = -1;
        ret(2, 2, 2) = 1;
        if sum(isnan(ret))+sum(isinf(ret))>0
            disp('asd')
        end
    case 'interval'
        ret = [ 0, 1 ];
    case 'linear'
        ret = 0;
    case 'parameters'
        ret = 1;
    case 'c'
        ret = [];
    case 'BV'
        ret=[
            za(1,1);
            za(2,1);
            zb(2,1)-p];
    case 'dBV'
        %DON'T CHANGE THIS LINE:
        ret = zeros(max(length(feval(mfilename,'c')),2-length(feval(mfilename,'c'))),length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders')));
        ret(1,1,1,1) = 1;
        ret(1,2,2,1) = 1;
        ret(2,3,2,1) = 1;
    case 'dP'
        ret = [
            -z(1,2)+(sqrt(t)-z(1,1))/(t-z(1,1));
            0];
    case 'dP_BV'
        ret = [
            0;
            0;
            -1];
    case 'initProfile'
        ret.initialMesh = linspace(0,1,6);
        ret.initialValues = ones(2,6);
        ret.parameters = 1;
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
