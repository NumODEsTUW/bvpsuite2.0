function [sol_coarse,sol_fine,halve] = errorestimate(sol_coarse,problem,settings,predictor)

%Definiere Werte der zu schätzenden Lösung
N=length(sol_coarse.x1)-1;
%b=x1((N)+1);
%a1=x1((0)+1);
ordnung=feval_problem(problem,'orders');
n=length(ordnung);%number of equations
collMethod=feval(settings,'collMethod');
collPoints=feval(settings,'collPoints');%position of the collocation points transformed to the interval [0 1]
parameter=feval_problem(problem,'parameters');
linear=feval_problem(problem,'linear');
halve=0; % needed for pathfollowing

% prepare predictor
if ~exist('predictor','var')
    predictor=[];
end
if ~isempty(predictor) && predictor.display
    fprintf('  <strong>Error estimation...</strong>\n')
end
if isempty(predictor)
    predictor = [];
elseif ~isfield(predictor,'tangent')
    tmp_disp.display=predictor.display;
    predictor=tmp_disp;
    clear('tmp_disp')
end

m = collPoints;
if (strcmp(collMethod,'user'))
    [rho] = collPoints;
    m=length(rho);
else
    [rho] = getStandardCollocationPoints(collMethod,m);
end

x1tau_coarse = sol_coarse.x1tau;
valx1tau_coarse = sol_coarse.valx1tau;
coeff_coarse = sol_coarse.coeff;

x1 = x1tau_coarse(1:m+1:end);
valx1 = valx1tau_coarse(:,1:m+1:end);

% Add a mesh point in-between each two mesh points
x1_2(1:2:2*length(x1)) = x1;
x1_2(2:2:2*length(x1)-1) = (x1(2:length(x1))+x1(1:length(x1)-1))/2;

if(linear) && ~isfield(predictor,'tangent')
    [~,~,sol_fine] = solveLinearProblem(problem,x1_2,rho); 
else
    try
        sol_coarse.lambda_p;
        predictor.a_0;
        
        initP.initialMesh=x1_2;
        
        % n_tmp = N*(sum(ordnung)+n*m);
        % coeff_a_0 = predictor.a_0(1:n_tmp);
        % param_a_0 = predictor.a_0(n_tmp+1:n_tmp+parameter);
        % coeff_tang = predictor.tangent(1:n_tmp);
        % param_tang = predictor.tangent(n_tmp+1:n_tmp+parameter);
        %
        % initP.initialValues = coeffToValues(coeff_a_0, predictor.x1,ordnung,rho,x1_2);
        % initP.parameters = param_a_0;
        % initP = initial_coefficients(problem,x1_2,initP,rho,0);
        % predictor.a_0 = [ initP.initialCoeff ; predictor.lambda_p_0 ]' ;
        %
        % initP.initialValues = coeffToValues(coeff_tang, predictor.x1,ordnung,rho,x1_2);
        % initP.parameters = param_tang;
        % initP = initial_coefficients(problem,x1_2,initP,rho,0);
        % predictor.tangent = [ initP.initialCoeff ; 0 ]';
        
        initP.initialValues = coeffToValues(sol_coarse.coeff, sol_coarse.x1,ordnung,rho,x1_2);
        initP.parameters = sol_coarse.parameters;
        initP = initial_coefficients(problem,x1_2,initP,rho,0);
        initProfile = [ initP.initialCoeff ; sol_coarse.lambda_p ] ;
        
        % For the computation on the fine mesh keep lambda_p fixed
        predictor.lpfix = 1;
    catch
        initProfile.initialMesh = x1;
        initProfile.initialValues = valx1;
        initProfile.initialCoeff = coeff_coarse;
        initProfile.parameters=sol_coarse.parameters;
    end
    [~,~,sol_fine,halve] = solveNonLinearProblem(problem,settings,x1_2,rho,initProfile,0,predictor);
    if halve~=0
        sol_coarse=0; sol_fine=0;
        return
    end
end

sol_coarse.x1tau = x1tau_coarse;
sol_coarse.valx1tau = valx1tau_coarse;
valx1tau2_2 = coeffToValues(sol_fine.coeff,sol_fine.x1,ordnung,rho,sol_coarse.x1tau); %fine solution on coarse grid + collpoints
sol_coarse.errest = (2^m / (1-2^m)) * (valx1tau2_2 - sol_coarse.valx1tau);

valx1tau = coeffToValues(sol_coarse.coeff,sol_coarse.x1,ordnung,rho,sol_fine.x1tau); %coarse solution on fine grid + collpoints
sol_fine.errest = (1 / (1-2^m)) * (sol_fine.valx1tau - valx1tau);

end

function Psireturn=Psi(n,rho,nr)
i=length(rho);
prod=1;
for s=1:i
    if (s~=n)
        prod=conv(prod,[1 -rho(s)])/(rho(n)-rho(s));
    end
end
for s=1:(nr)
    prod=polyint(prod);
end
Psireturn=prod;
%tabulars for collocation points
end

