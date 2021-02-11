function [ret] = bvps2_dae_probdef(request,z,za,zb,zc,t,p,lambda,lambda_p)
% Problem definition for DAE problem in Manual

J=1/2;
alpha=0;
rhobar=3;
b=10.3;

switch request
    case 'n'
        ret = 3;
    case 'orders'
        ret = [1 1 0];
    case 'problem'
        ret = [
            z(1,2)-z(2,1)*z(3,1)+alpha*J;
            z(2,2)-z(3,1)+1;
            z(1,1)-J^2/z(3,1)-z(3,1)];
    case 'jacobian'
        %DON'T CHANGE THIS LINE:
        ret = zeros(length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders'))+1);
        %Hier wurde das Format des return-Wertes festgelegt. Kann bei der Problemdefinition ignoriert werden.
        ret(1,1,2) = 1;
        ret(1,2,1) = -z(3,1);
        ret(1,3,1) = -z(2,1);
        ret(2,2,2) = 1;
        ret(2,3,1) = -1;
        ret(3,1,1) = 1;
        ret(3,3,1) = J^2/z(3,1)^2-1;
    case 'interval'
        ret = [0,b];
    case 'linear'
        ret = 0;
    case 'parameters'
        ret = 0;
    case 'c'
        ret = [];
    case 'BV'
        ret = [
            za(1,1)-J^2/rhobar-rhobar;
            zb(1,1)-J^2/rhobar-rhobar];
    case 'dBV'
        %DON'T CHANGE THIS LINE:
        ret = zeros(max(length(feval(mfilename,'c')),2-length(feval(mfilename,'c'))),length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders')));
        %Hier wurde das Format des return-Wertes festgelegt. Kann bei der Problemdefinition ignoriert werden.
        ret(1,1,1,1) = 1;
        ret(2,2,1,1) = 1;
    case 'dP'
        ret = [];
    case 'dP_BV'
        ret = [];
    case 'initProfile'
        ret.initialMesh = linspace(0,b,50);
        ret.initialValues = [
            rhobar*(1+J^2)/3*ones(1,length(ret.initialMesh));
            0*ones(1,length(ret.initialMesh));
            rhobar/3*ones(1,length(ret.initialMesh))];
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