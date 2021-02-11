function [new_x,halve,logstruct,itcount,fcount] = ...
    solve_nonlinear_sys(FDFhandle,x0,problem,tau,settings,psival,psi,rho,predictor)
%   SOLVE_NONLINEAR_SYS  Solve nonlinear systems obtained by SBVPCOL and SBVPERR
%
%   This routine is private and should not be accessed by other callers than SBVPCOL and SBVPERR

% ********** store some fields of bvpopt as variables for quick reference
AbsTol      = feval(settings,'absTolSolver');
RelTol      = feval(settings,'relTolSolver');
max_F_evals = feval(settings,'maxFunEvalsTRM');
max_iter    = feval(settings,'maxIterationsTRM');
if ~isempty(predictor) && isfield(predictor,'display')
    display     =  predictor.display; % Pathfollowing display messages
else
    display     = 1; % BVP solving display messages
end
log         = 0;
halve       = 0; % Needed in case of pathfollowing, otherwise not
TRM         = feval(settings,'allowTRM');

% set 1 to use lsqnonlin for the TRM. If set to 0, fsolve is used instead
lsq = 1;


% ********** some initializations
fcount = 1;   % number of function evaluations
itcount= 1;   % number of iterations
x = x0;
lambda = 1;
lambdamin = feval(settings,'lambdaMin');
nPreviousTRMIterates = 0;
logstruct = [];
updateJacFactor = feval(settings,'updateJacFactor');
switchToFFNFactor = feval(settings,'switchToFFNFactor');

if log
    logstruct.DOC = [];
    logstruct.tol_factor = [];
    logstruct.G0 = [];
    logstruct.G = [];
    logstruct.jac_update = [];
    logstruct.lambda = [];
    logstruct.nCorrections = [];
    
    logstruct.runtime.DFeval = [];
    logstruct.runtime.Feval = [];
    logstruct.runtime.lu = [];
    logstruct.runtime.resubstitute = [];
end

% Structure of the Algorithm:
% The Zerofinder consists of 3 different algorithms:
%    *) Fast Frozen Newton (cheap, but small domain of convergence (DOC) G1)
%    *) Predictor-Corrector Linesearch (expensive, but large DOC G2)
%    *) Trust Region Method (even more expensive, but even larger DOC G3)
%
% For the choice of the appropriate algorithm to use, it is necessary to determine the
% position of the current approximation w.r.t. the DOCs G1, G2, G3

% ********** Determine position of initial approximation

if display && isempty(predictor) % display information every iteration step
    fprintf('\n Determining zerofinder algorithm ... \n\n');
elseif display
    fprintf('  Performing Newton step...')
end

[TolFactor, DOC, new_x, U, L, G0, G, F, delta_x, simplified_delta_x, fcount, logstruct] = ...
    determine_position2(x,problem,tau,settings,psival,psi,lambdamin,fcount,logstruct,[],[],rho,predictor);


if TolFactor < 1
    if display && isempty(predictor)
        fprintf('\n\n Tolerances satisfied\n');
    elseif display
        fprintf(' Tolerances satisfied.\n')
    end
    return
end

if DOC == 1
    x = new_x; % Accept new approximation
end

% ********** initialize Zerofinder Display
if display  % display information every iteration step
    if ~isempty(predictor)
        fprintf('\n\n')
    end
    switch DOC
        case 1
            fprintf('  Fast Frozen Newton\n');
            fprintf('  STEP   F-COUNT   JAC_UPDATE   |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
        case 2
            fprintf('  Predictor-Corrector Line Search\n');
            fprintf('  STEP   F-COUNT     LAMBDA     |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
        case 3
            fprintf('  Trust Region Method\n');
    end
end

% ********** Choice-of-Algorithm Loop
while 1
    
    if log
        logstruct.DOC = [logstruct.DOC DOC];
    end
    
    if itcount > max_iter
        error(' Maximum number of iterations exceeded');
    end
    
    if fcount > max_F_evals
        error(' Maximum number of function evaluations exceeded');
    end
    
    
    
    switch DOC
        
        
        % **************************************************************************
        case 1 % ******************* Fast Frozen Newton ****************************
            % **************************************************************************
            % variables that have to be available: a (=x), G0, G, F (=F(x)), simplified_delta_x
            
            % ********** Determine if Jacobian has to be updated
            UpdateJac = G > updateJacFactor * G0;
            
            if log
                logstruct.lambda = [logstruct.lambda 1];
                logstruct.nCorrections = [logstruct.nCorrections 0];
                logstruct.jac_update = [logstruct.jac_update UpdateJac];
            end
            
            % ********** Determine Jacobian
            if UpdateJac % Update Jacobian
                cpt=cputime;
                DF = functionFDF('DF',problem,x,tau,psival,psi,rho,predictor);
                if log
                    logstruct.runtime.DFeval = [logstruct.runtime.DFeval cputime-cpt];
                    logstruct.condestDF = condest(DF);
                end
                
                cpt=cputime;
                [L,U]=lu(DF);                  % LU-decomposition of DF
                if log
                    logstruct.runtime.lu = [logstruct.runtime.lu cputime-cpt];
                end
                
                lastwarn('');                  % initialize warning state
                
                cpt=cputime;
                delta_x = U\(L\(- F));         % Newton correction
                if log
                    logstruct.runtime.resubstitute = ...
                        [logstruct.runtime.resubstitute cputime-cpt];
                end
                
                G0 = norm(delta_x);            % Norm of the Newton correction
            else % Iterate with frozen Jacobian
                delta_x = simplified_delta_x;
                G0 = G; % = norm(simplified_delta_x)
            end
            
            % ********** Perform Newton step
            new_x(:) = x(:) + delta_x;      % new approximation
            
            cpt = cputime;
            new_F = functionFDF('F',problem,new_x,tau,psival,psi,rho,predictor);   % new residual
            if log
                logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt];
            end
            
            fcount = fcount +1;             % increase number of function evaluations
            
            % ********** Determine Position of new solution approximation
            cpt = cputime;
            simplified_delta_x = U\(L\(-new_F));            % simplified delta_x
            if log
                logstruct.runtime.resubstitute = [logstruct.runtime.resubstitute cputime-cpt];
            end
            
            G = norm(simplified_delta_x);
            
            if log
                logstruct.G0 = [logstruct.G0 G0];
                logstruct.G =  [logstruct.G G];
            end
            
            if G < (1-lambdamin/2) * G0 % Approximation improved
                
                % ***** Check if Tolerances are satisfied
                % ***** TolFactor: Factor by which the correction had to be
                % ***** reduced in order to satisfy the tolerances
                TolFactor = check_tolerances(new_x,simplified_delta_x,AbsTol,RelTol);
                
                if log
                    logstruct.tol_factor = [logstruct.tol_factor TolFactor];
                end
                
                if TolFactor < 1
                    % Perform the last step towards the solution
                    new_x(:) = new_x(:) + simplified_delta_x;
                    
                    if display  % display information every iteration step
                        fprintf('  %3i   %5i    %8i     %13.2e    %.2e   %5.2e\n',...
                            itcount,fcount,UpdateJac,G/norm(x(:)),G0/G,TolFactor);
                    end
                    
                    if display && isempty(predictor) % off
                        fprintf('\n\n Tolerances satisfied\n');
                    elseif display
                        fprintf('\n  Tolerances satisfied.\n')
                    end
                    
                    return
                end
                
                DOC = 1;           % Keep iterating with the Fast Frozen Newton
                
                x = new_x;         % Accept new approximation
                F = new_F;
                
                if display  % display information every iteration step
                    fprintf('  %3i   %5i    %8i     %13.2e    %.2e   %5.2e\n',...
                        itcount,fcount,UpdateJac,G/norm(x(:)),G0/G,TolFactor);
                end
            elseif UpdateJac == 0 % Jacobian has not been updated in the previous step
                DOC = 1;           % -> Try FFN once again,
                %    but update Jacobian this time
                %    (We would have to evaluate the Jacobian
                %    anyway in the Line Search Algorithm)
                
                if display  % display information every iteration step
                    fprintf('  %3i   %5i    %8i     %13.2e    %.2e   %5.2e   NOT ACCEPTED\n',...
                        itcount,fcount,UpdateJac,G/norm(x(:)),G0/G,TolFactor);
                end
            else
                DOC = 2;           % -> Switch to Predictor Corrector Line Search
                
                % ********** Prepare for next step with the respective method
                lambda = 1;
                
                % Use Jacobian computed here
                if log
                    logstruct.jac_update  = [logstruct.jac_update  0];
                end
                
                if display
                    fprintf('  %3i   %5i    %8i     %13.2e    %.2e   %5.2e   NOT ACCEPTED\n',...
                        itcount,fcount,UpdateJac,G/norm(x(:)),G0/G,TolFactor);
                    
                    fprintf('\n  Switching to Predictor-Corrector Line Search\n');
                    fprintf('  STEP   F-COUNT     LAMBDA     |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
                end
            end
            
            itcount = itcount + 1;
            
            
            % **************************************************************************
        case 2 % *************** Predictor-Corrector Line Search *******************
            % **************************************************************************
            % Necessary variables at this point: a (=x), delta_x, simplified_delta_x,
            % lambda (=lambda_pred), G0, G (=G(lambda_pred))
            
            % Let G(lambda) = ||DF^(-1) * F(x + lambda*delta_x)||
            % At this point, x, delta_x, lambda, G(0), G'(0), G(lambda) are known
            % We check if lambda is accepted and correct lambda otherwise
            
            Accept_Lambda = G < (1-lambdamin/2)*G0;
            nCorrections = 0;
            
            while ~Accept_Lambda
                % ********** Repeatedly correct lambda until it is accepted or a termination condition is met
                
                if nCorrections == 0 % First correction, quadratical method
                    % We know G(0), G'(0) and G=G(lambda) by previous computations.
                    % To obtain a correction, we interpolate a quadratic polynomial through
                    % these values and define the minimum of this polynomial as the corrected lambda
                    % The respective polynomial reads G(x)=G(0) + x * G'(0) + x^2/lambda^2 *
                    % (G(lambda) - G(0) - lambda * G'(0)). Zeroing the first derivative yields:
                    
                    % store lambda (=lambda_pred) for possible further correction
                    l2 = lambda;
                    
                    Gprime0 = -G0;
                    lambda = - lambda^2 * Gprime0 / (2* (G - G0 - lambda * Gprime0));
                    
                    % Allow lambda_cor between lambda_pred / 10 and 1
                    lambda = max(l2/10 , min(lambda, 1));
                    
                    %***** Prepare for possible further correction (cubical, match notation with Num. Recipes)
                    l1 = lambda;
                    
                    nCorrections = 1;
                else % Subsequent corrections, cubical method
                    coeff = 1/(l1-l2) * [1/l1^2 -1/l2^2 ; -l2/l1^2 l1/l2^2] * ...
                        [G1 - Gprime0*l1 - G0 ; G2 - Gprime0*l2 - G0];
                    
                    ca = coeff(1);
                    cb = coeff(2);
                    
                    lambda = (-cb + sqrt(cb^2-3*ca*Gprime0))/(3*ca);
                    
                    % Allow lambda_cor between lambda_cor_old / 10 and lambda_cor_old /2
                    lambda = max(l1/10 , min(lambda, l1/2));
                    
                    % Prepare for possible further correction
                    l2 = l1;
                    l1 = lambda;
                    
                    nCorrections = nCorrections + 1;
                end
                
                if lambda < lambdamin
                    DOC = 3;  % switch to Trust Region Method
                    break;
                end
                
                new_x(:) = x(:) + lambda * delta_x;  % get new approximation
                
                cpt = cputime;
                F = functionFDF('F',problem,new_x,tau,psival,psi,rho,predictor);   % evaluate F at lambda_cor
                if log
                    logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt];
                end
                
                fcount = fcount + 1;
                
                cpt = cputime;
                simplified_delta_x = U\(L\-F);
                if log
                    logstruct.runtime.resubstitute = ...
                        [logstruct.runtime.resubstitute cputime-cpt];
                end
                
                G_old = G; % store old G for possible further correction
                G = norm(simplified_delta_x);
                
                Accept_Lambda = G < (1-lambdamin/2) * G0;
                
                % ***** Prepare for possible further correction
                G1 = G;
                G2 = G_old;
            end
            
            if Accept_Lambda % a lambda has been accepted
                x = new_x;    % Accept new solution
                
                if G < switchToFFNFactor * G0
                    DOC = 1;
                else
                    DOC = 2;
                end
            end
            
            itcount = itcount + 1;
            
            if log
                logstruct.lambda = [logstruct.lambda lambda];
                logstruct.nCorrections = [logstruct.nCorrections nCorrections];
                logstruct.G0 = [logstruct.G0 G0];
                logstruct.G =  [logstruct.G G];
            end
            
            % At this point, DOC has been determined. Check Termination conditions if DOC ~= 3
            if DOC ~= 3
                TolFactor = check_tolerances(x,simplified_delta_x,AbsTol,RelTol);
                
                if log
                    logstruct.tol_factor = [logstruct.tol_factor TolFactor];
                end
                
                if TolFactor < 1
                    new_x(:) = x(:) + simplified_delta_x;
                    
                    if display  % display information every iteration step
                        fprintf('  %3i   %5i        %5.4f    %12.2e    %.2e   %5.2e\n',...
                            itcount,fcount,lambda,G/norm(x(:)),G0/G,TolFactor);
                    end
                    
                    if display && isempty(predictor) % off
                        fprintf('\n\n Tolerances satisfied\n');
                    elseif display
                        fprintf('\n  Tolerances satisfied.\n')
                    end
                    
                    return
                end
            elseif log
                logstruct.tol_factor = [logstruct.tol_factor -1];
            end
            
            
            
            
            if display  % display information every iteration step
                fprintf('  %3i   %5i        %5.4f    %12.2e    %.2e   %5.2e\n',...
                    itcount,fcount,lambda,G/norm(x(:)),G0/G,TolFactor);
                
                if DOC == 1
                    fprintf('\n  Switching to Fast Frozen Newton\n');
                    fprintf('  STEP   F-COUNT   JAC_UPDATE   |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
                elseif DOC == 3
                    fprintf('\n  Switching to Trust Region Method\n');
                end
            end
            
            
            % If we keep up doing the Line Search, we have to predict lambda for the next step
            if DOC == 2
                
                % ********** Keep old values of delta_x for prediction
                G0_old = G0;
                
                % ********** Determine delta_x
                cpt = cputime;
                DF = functionFDF('DF',problem,x,tau,psival,psi,rho,predictor);
                if log
                    logstruct.runtime.DFeval = [logstruct.runtime.DFeval cputime-cpt];
                    logstruct.jac_update = [logstruct.jac_update 1];
                    logstruct.condestDF = condest(DF);
                end
                
                cpt=cputime;
                [L,U]=lu(DF);                  % LU-decomposition of DF
                if log
                    logstruct.runtime.lu = [logstruct.runtime.lu cputime-cpt];
                end
                
                lastwarn('');                  % initialize warning state
                
                cpt=cputime;
                delta_x = U\(L\(- F));         % Newton correction
                if log
                    logstruct.runtime.resubstitute = ...
                        [logstruct.runtime.resubstitute cputime-cpt];
                end
                
                G0 = norm(delta_x);
                
                % ********** Predict lambda [Deuflhardt74]
                %dpr%      warning off MATLAB:divideByZero  % If Jacobians are the same, delta_x = simplified_delta_x
                mu = lambda * G0_old / norm(simplified_delta_x - delta_x);
                %dpr%      warning on MATLAB:divideByZero
                
                if mu > 0.7
                    lambda = 1;
                else
                    lambda = mu;
                end
                
                % ********** Make sure all necessary variables are available for correction step
                new_x(:) = x(:) + lambda * delta_x;
                
                cpt = cputime;
                F = functionFDF('F',problem,new_x,tau,psival,psi,rho,predictor);
                if log
                    logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt];
                end
                
                fcount = fcount + 1;
                
                cpt = cputime;
                simplified_delta_x = U\(L\-F);
                if log
                    logstruct.runtime.resubstitute = ...
                        [logstruct.runtime.resubstitute cputime-cpt];
                end
                
                G = norm(simplified_delta_x);
            end
            
            
            % **************************************************************************
        case 3 % ********************* Trust Region Method *************************
            % **************************************************************************
            
            if isfield(predictor,'tangent')
                fprintf('  Trust Region Method is turned off during pathfollowing.\n')
                halve=3;
                return
            end
            
            if(~TRM)
                disp('tried to use TRM. Solution may be imprecise or wrong');
                return;
            end
            
            
            % Necessary variables at this point: a (=x), nPreviousTRMIterates
            % ********** Determine accuracy for Trust Region Method
            options = optimset('Display','off','MaxIter',max_iter,'MaxFunEvals',max_F_evals);
            
            options.TolX = sqrt(max(RelTol,AbsTol));
            options.TolFun = 0;        % Assemble necessary options for FSOLVE
            options.Jacobian = 'on';   % (Display, MaxIter, MaxFunEvals are user parameters)
            options.LargeScale = 'on';
            switch nPreviousTRMIterates
                case 0
                    options.TolX = sqrt(max(RelTol,AbsTol));
                    options.TolFun = 0;        % Assemble necessary options for FSOLVE
                    options.Jacobian = 'on';   % (Display, MaxIter, MaxFunEvals are user parameters)
                    options.LargeScale = 'on';
                case 1
                    options.TolX = max(RelTol,AbsTol);
                case 2
                    options.TolX = 1000*eps;
                case 3
                    options.Algorithm =  'levenberg-marquardt';
                case 4
                    fprintf('Zerofinder cannot fulfill the given tolerances');
                    return;
            end
            
            
            options.Display = 'iter';
            % *********** Limit the maximum number of iterations and FunEvals for TRM
            options.MaxIter = max_iter - itcount;
            options.MaxFunEvals = max_F_evals - fcount;
            
            % ********** Invoke Trust Region Method
            if lsq
                [x,~,F,ExitFlag,lsqLog,~,DF] = ...
                    lsqnonlin(FDFhandle,x,[],[],options,problem,tau,psival,psi,rho,predictor);
            else
                [x,F,ExitFlag,lsqLog,DF] = ...
                    fsolve(FDFhandle,x,options,problem,tau,psival,psi,rho,predictor);
            end
            
            nPreviousTRMIterates = nPreviousTRMIterates +1;
            fcount = fcount + lsqLog.funcCount;
            
            % ***** Determine whether anything has gone wrong
            if nPreviousTRMIterates==3
                switch ExitFlag
                    case 0
                        fprintf('Number of iterations or number of function evaluations exeed the given borders');
                        return;
                    case -1
                        fprintf('Zerofinding algorithm did not converge');
                        return;
                    case -2
                        fprintf('Zerofinding algorithm converged against a point that is no solution');
                        return;
                    case -3
                        fprintf('Zerofinding algorithm did not converge');
                        return;
                    case -4
                        fprintf('Zerofinding algorithm did not converge');
                        return;
                end
            end
            
            
            
            % ********** Determine position of new approximation
            [TolFactor, DOC, new_x, U, L, G0, G, F, delta_x, simplified_delta_x, fcount, logstruct] = ...
                determine_position2(x,problem,tau,settings,psival,psi,lambdamin,fcount,logstruct,F,DF,rho,predictor);
            
            if TolFactor < 1
                if display % off
                    fprintf('\n\n Tolerances satisfied\n');
                end
                return
            end
            
            % ***** We only log those quantities if we continue iterating
            if log
                logstruct.lambda = [logstruct.lambda -1];
                logstruct.nCorrections = [logstruct.nCorrections -1];
                logstruct.jac_update = [logstruct.jac_update -1];
                logstruct.G0 = [logstruct.G0 G0];
                logstruct.G =  [logstruct.G G];
                logstruct.tol_factor = [logstruct.tol_factor TolFactor];
            end
            
            
            if DOC == 1
                x = new_x; % Accept new approximation
            elseif DOC == 2  % use Jacobian computed here
                if log
                    logstruct.jac_update  = [logstruct.jac_update  0];
                end %dpr%
            end
            
            
            if display  % display information every iteration step
                switch DOC
                    case 1
                        fprintf('  Switching to Fast Frozen Newton\n');
                        fprintf('  STEP   F-COUNT   JAC_UPDATE   |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
                    case 2
                        fprintf('  Switching to Predictor-Corrector Line Search\n');
                        fprintf('  STEP   F-COUNT     LAMBDA     |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
                    case 3
                        fprintf('  Continuing with Trust Region Method\n');
                end
            end
    end
end
end



% *************************************************************************
% ********** DETERMINE POSITION OF INITIAL APPROXIMATION ******************
% *************************************************************************

function [TolFactor, DOC, new_x, U, L, G0, G, F, delta_x, simplified_delta_x, fcount, logstruct] = ...
    determine_position2(x,problem,tau,settings,psival,psi,lambdamin,fcount,logstruct,F,DF,rho,predictor)

log = 0;
AbsTol      = feval(settings,'absTolSolver');
RelTol      = feval(settings,'relTolSolver');
switchToFFNFactor = feval(settings,'switchToFFNFactor');

% ********** Determine Newton correction
% ***** Evaluate Jabobian, if not provided (by fsolve)
if isempty(F)
    cpt = cputime;
    DF=functionFDF('DF',problem,x,tau,psival,psi,rho,predictor);
    if log
        logstruct.runtime.DFeval = [logstruct.runtime.DFeval cputime-cpt];
        logstruct.condestDF = condest(DF);
    end
    
    cpt = cputime;
    F = functionFDF('F',problem,x,tau,psival,psi,rho,predictor);   % new residual
    if log
        logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt];
    end
    
    fcount = fcount +1;            % increase number of function evaluations
end

cpt=cputime;
[L,U]=lu(DF);                  % LU-decomposition of DF
if log
    logstruct.runtime.lu = [logstruct.runtime.lu cputime-cpt];
end

lastwarn('');                  % initialize warning state

cpt=cputime;

delta_x = U\(L\(- F));         % Newton correction
if log
    logstruct.runtime.resubstitute = [logstruct.runtime.resubstitute cputime-cpt];
end


G0 = norm(delta_x);            % Norm of the Newton correction

new_x = zeros(size(x));
new_x(:) = x(:) + delta_x;     % new approximation

% ***** In case we have very good initial approximation (e.g. from a previous mesh),
% ***** tolerances might be satisfied even now and we reduce computational
% ***** costs to the absolute minimum of one F- and one DF-evaluation
TolFactor = check_tolerances(new_x,delta_x,AbsTol,RelTol);
if TolFactor < 1
    DOC = [];  % dummy outputs
    G = [];
    simplified_delta_x = [];
    
    if log
        logstruct.tol_factor = [logstruct.tol_factor TolFactor];
        logstruct.G0 = [logstruct.G0 G0];
    end
    
    return
end

% ********** Check Monotonicity Condition for lambda = 1 and lambda = lambda_min
cpt = cputime;
F = functionFDF('F',problem,new_x,tau,psival,psi,rho,predictor);   % new residual
if log
    logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt];
end

fcount = fcount +1;             % increase number of function evaluations

cpt = cputime;
simplified_delta_x = U\(L\(-F));            % simplified delta_x
if log
    logstruct.runtime.resubstitute = [logstruct.runtime.resubstitute cputime-cpt];
end

G = norm(simplified_delta_x);

if G <(1-lambdamin/2) * G0 % Approximation improved
    if G < switchToFFNFactor * G0         % Approximation improved significantly
        DOC = 1;             % -> Fast Frozen Newton
        
        % See if tolerances are satisfied
        TolFactor = check_tolerances(new_x,simplified_delta_x,AbsTol,RelTol);
        if TolFactor < 1
            % Perform one last step towards the solution
            new_x(:) = new_x(:) + simplified_delta_x;
            
            if log
                logstruct.tol_factor = [logstruct.tol_factor TolFactor];
                logstruct.G0 = [logstruct.G0 G0];
                logstruct.G = [logstruct.G G];
            end
            
            return
        end
    else % Approximation improved slightly
        DOC = 2;             % -> Predictor Corrector Line Search
        if log
            logstruct.jac_update  = [logstruct.jac_update  0];   % Use Jacobian computed here
        end
    end
else % try another Newton step with smallest possible lambda
    
    new_x(:) = x(:) + lambdamin * delta_x;      % new approximation
    
    cpt = cputime;
    Fmin = functionFDF('F',problem,new_x,tau,psival,psi,rho,predictor);   % new residual
    if log
        logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt];
    end
    
    fcount = fcount +1;
    
    cpt = cputime;
    simplified_delta_x_min = U\(L\(-Fmin));
    if log
        logstruct.runtime.resubstitute = [logstruct.runtime.resubstitute cputime-cpt];
    end
    
    Gmin = norm(simplified_delta_x_min);
    
    if Gmin < (1-lambdamin/2) * G0 % Approximation improved
        DOC = 2; % Predictor Corrector Line Search
        if log
            logstruct.jac_update  = [logstruct.jac_update  0];   % Use Jacobian computed here
        end
    else
        DOC = 3; % -> Trust Region Method
    end
end
end


%********************************************************************************
%******************************* CHECK TOLERANCES *******************************
%********************************************************************************

function TolFactor = check_tolerances(x,delta,AbsTol,RelTol)
% Calculates the factor by which delta has to be reduced to satisfy the
% tolerances.
% out <= 1   => Tolerances satisfied
% out >  1   => Tolerances not satisfied

% ***** Determine Dimension of solution array x
[d , ~, N] = size(x);

% ***** Extract solution components at mesh points (last meshpoint (b) not included
InfNorm_y = max(abs(reshape(x(:,1,:),d,N)));

% ***** Extract last correction at mesh points
delta = reshape(delta,size(x));
InfNorm_delta = max(abs(reshape(delta(:,1,:),d,N)));

% ***** Componentwise check if tolerances are satisfied
TolFactor = max(InfNorm_delta ./ (AbsTol + InfNorm_y * RelTol));

if TolFactor<1
    pause(1)
end
end

