function [ residual ] = computeResidual( problem, coeff, parameters, tau, rho, teval, lambda_p )
% For a given (parameter-dependent) problem with approximate solution given
% by coeff on a mesh defined by tau and rho, the problem is evaluated at
% the points in teval.

orders = feval_problem(problem,'orders');

vals = zeros(length(orders),length(teval),max(orders)+1); %Pre-locate space for the values of the solution, the residual shall be computed of.

for oi = 0:max(orders)
    vals(:,:,oi+1) = coeffToValues(coeff,tau,orders,rho,teval,oi);
end

residual=zeros(length(orders),length(teval));
for k=1:length(teval)
   residual(:,k) =  feval_problem(problem,'problem',reshape(vals(:,k,:),length(orders),max(orders)+1),[],[],[],teval(k),parameters,[],[],lambda_p);
end

end

