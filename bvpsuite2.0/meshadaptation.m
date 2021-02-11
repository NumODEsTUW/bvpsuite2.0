function [finalsol,toleriert,halve]=meshadaptation(problem,settings,x1,rhoColl,initProfile,predictor)

% meshadaptation ... calculate solution using adaptive meshes
%
% meshadaptation2(aTOL,rTOL,K,bvpfile,plotsol,x1,start,bvpopt,plotres,max
% iter) calculates the solution of a BVP, given in
% bvpfile using an adaptive mesh technique in order to solve the prescribed
% tolerance requirements given in aTOL (absolute tolerance) and rTOL
% (relative tolerance).
%
% K limits the maximal ratio between the shortest and the longest
% subinterval. plotsol (0,1) toggles the plot of the solution after the
% solution process. An initial mesh and an initial solution can be passed
% to the routine using x1 and start. If left empty, the mesh is determined
% from the bvpfile and the initial profile 0 is used. maxiter regulates the
% number of allowed adaptation steps.
%
% See also BVPSUITE, RUN
%

aTOL = feval(settings,'absTolMeshAdaptation');
rTOL = feval(settings,'relTolMeshAdaptation');
K = feval(settings,'K');
maxiter = feval(settings,'maxAdaptations');
finemesh = feval(settings,'finemesh');
halve=0; % Needed in case of pathfollowing, otherwise not

plotres = 0;

ordnung=feval_problem(problem,'orders');
n=length(ordnung);%number of equations
c=feval_problem(problem,'c');%boundary values inside the interval
collMethod=feval(settings,'collMethod');
collPoints=feval(settings,'collPoints');%position of the collocation points transformed to the interval [0 1]
parameter=feval_problem(problem,'parameters');
linear=feval_problem(problem,'linear');

m = length(rhoColl);

% prepare lambda_p and predictor
try
    pathfoll=feval(problem,'pathfollowing');
    lambda_p=pathfoll.start;
catch
    lambda_p=0;
end
if ~exist('predictor','var')
    predictor=[];
elseif ~isfield(predictor,'tangent')
    tmp_disp.display=predictor.display;
    predictor=tmp_disp;
    clear('tmp_disp')
else
    % The step-length might get halved, whenever N gets bigger then
    % meshFactorMax*N_orig
    % Before points are added to the mesh, the step-length will always get
    % halved first as a precautionary measure.
    meshFactorMax = max(1,predictor.meshFactorMax);
    N_orig = length(predictor.x1)-2;
end

% find initial mesh if not transmitted to meshadaptation
if length(x1) < 1
    x1 = feval(settings,'mesh');
    interval = feval_problem(problem,'interval');
    x1 = (interval(2)-interval(1))/(x1(end)-x1(1))*(x1+x1(1))+interval(1);
end
left_end = x1(1);
right_end = x1(end);

% Calculate mesh density rho from the mesh
rho=(1./diff(x1)/length(x1))';
update_mode=1; % 1 .. update rho, 2 .. update N

M=length(x1)-1; % number of intervals

N = M-1;
Ncompare = 1e8; % initialised with a high value
order_method = m;

% Variables
minimprove = 0.1; % 0.1 means at least 10 % points less to try further density adjustment
safetyN = 0; % 0.1 means at least 10 % more points than actually suggested - not accurate anymore - use safety_sigma instead
safety_sigma = 0.9; % safety factor used in the formula of N-preestimation
LPFilter_k = order_method*2.75; % controls the gain of the adjustment of the mesh density - should be coppled to the method order
res_mode = 'tcolmid'; % triggers the choice of points for the evaluation of the residual

numequ=length(ordnung);
if length(aTOL) <= numequ
    aTOL=ones(1,numequ)*aTOL(1);
end
if length(rTOL) <= numequ
    rTOL=ones(1,numequ)*rTOL(1);
end


% DEFINITIONS ==================================
% Number of INTERVALS is M
% Number of internal GRID POINTS is N = M-1
% STEP SIZE on the auxiliary grid is dxi = 1/M
% rho is a vector with M components
% ==============================================

% Main loop: solution of problem + grid refinement

densityUpdates = 0;

for j=1:maxiter
    % ***************************
    % Calculate next mesh from rho
    %*****************************
    
    % *************
    % SOLVE PROBLEM
    % *************
    
    if(feval_problem(problem,'linear')) && ( isempty(predictor) || ~isfield(predictor,'tangent') )
        [tcol,ycol,solstruct] = solveLinearProblem(problem,x1,rhoColl);
        tau = solstruct.x1;
        y = solstruct.valx1;
        coeff = solstruct.coeff;
    else
        [tcol,ycol,solstruct,halve] = solveNonLinearProblem(problem,settings,x1,rhoColl,initProfile,1-(j>1),predictor);
        if halve==0
            tau = solstruct.x1;
            y = solstruct.valx1;
            coeff = solstruct.coeff;
            if ( isempty(predictor) || ~isfield(predictor,'tangent') )
                initProfile.initialMesh = tau;
                initProfile.initialValues = y;
                initProfile.initialCoeff = solstruct.coeff;
                initProfile.parameters = solstruct.parameters;
            else
                initProfile = [solstruct.coeff', solstruct.parameters, solstruct.lambda_p ]';
                solstruct.initCoeff = initProfile;
                lambda_p = solstruct.lambda_p ;
            end
        else
            toleriert=0; finalsol=0;
            return
        end
    end
    
    par = solstruct.parameters;
    
    % **************
    % ESTIMATE ERROR
    % **************
    
    [sol_coarse,sol_fine,halve] = errorestimate(solstruct,problem,settings,predictor);
    if halve~=0
        toleriert=0; finalsol=0;
        return
    end
    ecol = sol_coarse.errest;
    ecoeff = sol_coarse.coeff;
    tcol2 = sol_fine.x1tau;
    ycol2 = sol_fine.valx1tau;
    tau2 = sol_fine.x1;
    y2 = sol_fine.valx1;
    ecol2 = sol_fine.errest;
    
    % Calculate error norms - global error ;
    %ycol_norm = max(abs(ycol),[],1)
    %err_norm = max(abs(ecol),[],1)
    ycol_norm = abs(ycol);
    err_norm = abs(ecol);
    
    %ycol_norm2 = max(abs(ycol2),[],1);
    %err_norm2 = max(abs(ecol2),[],1);
    ycol_norm2 = abs(ycol2);
    err_norm2 = abs(ecol2);
    
    qTol = zeros(size(err_norm,1),1);
    jmax = qTol;
    qTol2 = qTol;
    jmax2 = qTol;
    for tolt=1:size(err_norm,1)
        [qTol(tolt),jmax(tolt)] = max(err_norm(tolt,:) ./ (aTOL(tolt) + ycol_norm(tolt,:) .* rTOL(tolt))); % qTol is about tolfactor = errmax/tol
        [qTol2(tolt),jmax2(tolt)] = max(err_norm2(tolt,:) ./ (aTOL(tolt) + ycol_norm2(tolt,:) .* rTOL(tolt))); % qTol is about tolfactor = errmax/tol
    end
    
    if plotres
        disp(['Tolerance Factor: ' num2str(qTol)])
    end
    if finemesh
        if plotres
            disp(['Tolerance Factor on the fine grid: ' num2str(qTol2)])
        end
    end
    
    if qTol < 1
        if plotres
            disp('Tolerance satisfied')
        end
        x1=tau;
        valx1=y;
        x1tau=tcol;
        valx1tau=ycol;
        toleriert=0;
        break
    end
    if finemesh == 1
        if qTol2 < 1
            if plotres
                disp('Tolerance satisfied on the fine grid')
            end
            coeff=ecoeff;
            x1=tau2;
            valx1=y2;
            x1tau=tcol2;
            valx1tau=ycol2;
            toleriert=1;
            break
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%;
    % Update Procedure %;
    %%%%%%%%%%%%%%%%%%%%;
    
    if update_mode == 2; % update the number of mesh_points
        if densityUpdates <2
            update_mode = 1;
            densityUpdates = 1;
        else
            Nold = N;
            N = ceil((1+safetyN)*N*(max(max(err_norm,[],2)' ./ (safety_sigma*(aTOL))).^(1/(order_method+1))));
            %N = max(N,ceil(1.5^(1/order_method)*Nold)); % max increase
            % of number of points
        end
    elseif update_mode == 1; % update the mesh density
        % Decide about the next step - Update N or rho ;
        Ncompare_old = Ncompare;
        Ncompare = ceil((1+safetyN)*N*(max(max(err_norm,[],2)' ./ (safety_sigma*(aTOL))).^(1/(order_method+1))));
        if plotres
            disp(['N suggested: ' num2str(Ncompare)])
        end
        
        if Ncompare > (1-minimprove)*Ncompare_old % check if density is finally adjusted
            update_mode = 2;
            % here we have to do the following - it can happen, that
            % the suggested number of points increases by a huge amount
            % so it would be a bad idea to update the density in this
            % last step - we have to use the older (better) density
            % here - this may happen in rare cases
            if Ncompare <= Ncompare_old || (~exist('rhold','var'))
                N = Ncompare;
                Ncompare = 1e8;
            else
                N = Ncompare_old;
                Ncompare = 1e8;
                rho=rhold;
                x1=x1old;
            end
        else
            densityUpdates = densityUpdates + 1;
        end
    end
    
    % For pathfollowing a maximal number of allowed mesh points is given
    if isfield(predictor,'tangent') && N>=N_orig*meshFactorMax
        fprintf('\n  Suggested number of mesh points: %i\n Number of mesh points allowed in this step: %i\n', N+2, max(N_orig,N_orig*meshFactorMax))
        toleriert=0; finalsol=0;
        halve=1;
        return
    end
    
    % ******************
    % CALCULATE RESIDUAL
    % ******************
    if update_mode==1
        switch res_mode
            case 'mesh'
                resmesh = x1;
                %residual=abs(equations('residual',coeff,problem,tau,resmesh));
                residual=abs(computeResidual(problem,coeff,par,tau,rhoColl,resmesh,lambda_p));
                % to get the residual for adjusting the mesh density we need a single
                % vector, so we take the 2-norm of the integral of it
                
                residualfinal=zeros(1,size(residual,2)-1);
                for i=1:size(residual,1)
                    residualadd=(residual(i,1:end-1)+residual(i,2:end)).*diff(x1); % this means value on the lhs and rhs times the length
                    residualfinal=residualfinal+residualadd.^2;
                end
                residualfinal=sqrt(residualfinal);
            case 'tcol'
                resmesh = tcol;
                
                residual=abs(computeResidual(problem,coeff,par,tau,rhoColl,resmesh,lambda_p));
                % to get the residual for adjusting the mesh density we need a single
                % vector, so we take the 2-norm of the integral of it
                
                residualfinal=zeros(1,length(x1)-1);
                for i=1:size(residual,1)
                    residualadd=(residual(i,1:end-1)+residual(i,2:end)).*diff(tcol);
                    for k=1:length(x1)-1
                        residualadd2(k)=sum(residualadd(k*(order_method+1)-order_method:k*(order_method+1)));
                    end
                    residualfinal=residualfinal+residualadd2.^2;
                end
                residualfinal=sqrt(residualfinal);
                
            case 'midpoints'
                resmesh = x1(1:end-1)+diff(x1)/2;
                residual=abs(computeResidual(problem,coeff,par,tau,rhoColl,resmesh,lambda_p));
                residualfinal=zeros(1,size(residual,2));
                for i=1:size(residual,1)
                    residualadd=(residual(i,:)).*diff(x1); % this is the value at the midpoint times the length
                    residualfinal=residualfinal+residualadd.^2;
                end
                residualfinal=sqrt(residualfinal);
            case 'tcolmid'
                resmesh = tcol(1:end-1)+diff(tcol)/2;
                residual=abs(computeResidual(problem,coeff,par,tau,rhoColl,resmesh,lambda_p));
                % to get the residual for adjusting the mesh density we need a single
                % vector, so we take the 2-norm of the integral of it
                
                residualfinal=zeros(1,length(x1)-1);
                for i=1:size(residual,1)
                    residualadd=residual(i,:).*diff(tcol);
                    for k=1:length(x1)-1
                        residualadd2(k)=sum(residualadd(k*(order_method+1)-order_method:k*(order_method+1)))/(x1(k+1)-x1(k));
                    end
                    residualfinal=residualfinal+residualadd2.^2;
                end
                residualfinal=sqrt(residualfinal);
        end
    end
    
    Nnew = N;
    Mnew = Nnew+1;
    
    % Display stuff
    if plotres
        disp(['N used:      ' num2str(Nnew)])
        disp(['---------------'])
        disp(['Iteration ' num2str(j)])
    end
    if update_mode == 1
        if plotres
            disp(['Updating density'])
        end
        fprintf(1,'%s%d\n','  Density Update N=',N);
    elseif update_mode == 2
        if plotres
            disp('Updating N')
        end
        fprintf(1,'%s%d\n','N Update N=',N);
    end
    rhold = rho;
    x1old = x1;
    if update_mode == 1;
        % Update density function
        rho = LPfilter(residualfinal',rho,LPFilter_k);
    end;
    
    % check if the grid is within K-range
    maxrho = max(rho);
    rho=max(rho,maxrho/K);
    
    rho = bcTDFlogV4(rho);           % Smoothing
    rho = rho*sum(1./rho)/Mnew;      % Normalize
    
    x1 = (x1-x1(1))/(x1(end)-x1(1));
    I = cumtrapz(x1,([rho(1); rho]+[rho; rho(end)])'/2);
    x1 = pchip(I/I(end),x1,0:1/(Mnew):1);
    rho = 1./diff(x1)';
    x1 = left_end+(x1-0)*(right_end-left_end)/(1-0); %umrechnen von [0,1] auf [a,b]
    
    % Adapt the predictor at the end of each iteration
    fprintf('  Maximal error estimate value: %2.3e \n',max(max(abs(sol_coarse.errest))))
    if ~isempty(predictor) && isfield(predictor,'tangent')
        % update predictor
        initP.initialMesh=x1;
        
        % N_tmp = length(predictor.x1)-1;
        % n_tmp = N_tmp*(sum(ordnung)+n*m);
        % coeff_a_0 = predictor.a_0(1:n_tmp);
        % param_a_0 = predictor.a_0(n_tmp+1:n_tmp+parameter);
        %
        % initP.initialValues = coeffToValues(coeff_a_0, predictor.x1,ordnung,rhoColl,x1);
        % initP.parameters = param_a_0;
        % initP = initial_coefficients(problem,x1,initP,rhoColl,0);
        % predictor.a_0 = [ initP.initialCoeff ; predictor.lambda_p_0 ]' ;
        % predictor.tangent = update_tangent(problem,predictor.a_0,x1,rhoColl,ordnung).';
        % predictor.steplength = (predictor.lambda_p_p - predictor.a_0(end)) / predictor.tangent(end);
        
        % update initial Profile
        initP.initialValues = coeffToValues(initProfile, predictor.x1,ordnung,rhoColl,x1);
        initP.parameters = sol_coarse.parameters;
        initP = initial_coefficients(problem,x1,initP,rhoColl,0);
        initProfile = [ initP.initialCoeff ; lambda_p ] ;
        
        % For the computation on the new mesh keep lambda_p fixed
        predictor.lpfix = 1;
        
        predictor.x1 = x1;
    end
end

finalsol.x1 = x1;
finalsol.valx1 = valx1;
finalsol.x1tau = x1tau;
finalsol.valx1tau = valx1tau;
finalsol.coeff = coeff;
finalsol.parameters = par;
finalsol.errest = ecol;
if ~isempty(predictor)
    finalsol.lambda_p =  lambda_p;
end

end


% **********************************************************
% **********************************************************
% Functions for processing the residual ********************
% **********************************************************
% **********************************************************


function tosignal = resampleV4(fromgrid,fromsignal,togrid)
% Oversampling that guarantees that tosignal remains positive
% Written by GS, TU Wien, 24 August 2006
% V4 is original version

signal = abs(fromsignal);                   % For robustness only
meansig = norm(signal,1)/length(signal);
mag = max(signal);
signal = signal/mag;                        % Normalize to maximum of 1
signal = signal + 1e-10*exp(-10*signal);      % Limiter lifts values near zero 1e-3
signal = log(signal);                       % Transform to logarithmic scale
hi = max(signal);
lo = min(signal);
%signal = pchip(fromgrid,signal,togrid); % Resample on new grid
signal = max(lo,min(hi,signal));         % Make sure amplitude doesn't grow
signal = exp(signal);                       % Now mapped to new grid and positive
signal = signal + 1e-10*exp(-10*signal);      % Limiter lifts values near zero 1e-2
signal = bcTDFlogV4(signal);                % Smoothing
signal = mag*signal/max(signal);            % Restore signal amplitude
newmean = norm(signal,1)/length(signal);    % Correct the signal's mean value
signal = signal + meansig - newmean;
signal = max(0,signal);                     % Protect against negative values
signal = signal + 1e-10*exp(-10*signal);      % Limiter lifts values near zero 1e-3
tosignal = signal;                          % Output
end

function newrho = LPfilter(err,rho,LPFilter_k)
% Generate step density profile update
% For use with 2pBVP solver for 1st order systems
% Euler-Lagrange optimal grid is generated using
% using deadbeat control law. Local error is
% equidistributed over the interval.

% err is already on staggered grid in this version
k = LPFilter_k;       % k = 3 for DBC, k = 4 conv filter (root 1/3)
M = length(rho);

% Process input error
scalederr = M*abs(err).^(1/k);
% scalederr = processV4(bcTDFlogV4(scalederr)); % SWITCHED OFF

ratio = scalederr;
% ratio is suggested update of rho; compute next rho
rhonew = rho.*ratio;
%rhonew = bcTDFlogV4(rhonew);                   % SWITCHED OFF
newrho = rhonew;%*sum(1./rhonew)/M;
end


function out = bcTDFlogV4(in)
% Boundary corrected Topelitz digital filter, multiplicative (logarithmic)
% Applies 2nd order convolution LP filter repeatedly
% to input signal until output signal is sufficiently smooth

smooth = 50;     % Define smoothness requirement
N = length(in);
old = in;
signal = old;    % Mem alloc
snratio = 0;

while snratio < smooth
    % Apply filter until S/N ratio at this stage is at least = smooth
    for i=2:N-1
        signal(i) = (old(i-1)*old(i)^2*old(i+1))^(1/4);
    end
    % Boundary corrections
    signal(1) = (old(1)^3*old(2)^3/old(3)^2/old(4)*old(5))^(1/4);
    % signal(1) = (old(1)^2*old(2)^3/old(3))^(1/4);
    signal(N) = (old(N)^3*old(N-1)^3/old(N-2)^2/old(N-3)*old(N-4))^(1/4);
    % signal(N) = (old(N)^2*old(N-1)^3/old(N-2))^(1/4);
    
    % Compute noise
    noise = old - signal;
    s = norm(signal);
    n = norm(noise);
    snratio = s/(n + 1e-3*s);
    old = signal;
end

out = signal;
end

function grid = rhogrid(rho)
% Construct the staggered grid for the vector function rho

% Staggered grid
M = length(rho);
offset = 1/M/2;
grid = linspace(offset,1-offset,M)';
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

function new_tangent = update_tangent(problem,initProfile,x1_2,rho,ordnung)
m=length(rho);
psi=zeros(m,max(ordnung)+m,max(ordnung));
for ord=1:max(ordnung)
    for ii=1:m
        psi(ii,1+max(ordnung)-ord:m+max(ordnung),ord)=Psi(ii,rho,ord);
    end
end
psival=zeros(max(ordnung),m,m+2);
for ord=1:max(ordnung)
    for ii=1:m
        %evaluation of psi
        psival(ord,ii,1:m)=polyval(psi(ii,:,ord),rho(1:m));
        psival(ord,ii,m+1)=polyval(psi(ii,:,ord),1);
        psival(ord,ii,m+2)=polyval(psi(ii,:,ord),0);
    end
end
jac_F = functionFDF( 'DF', problem ,initProfile,x1_2,psival,psi,rho,[]);
help1=zeros(length(jac_F(:,1)),1);
help1(end)=1;
new_tangent =jac_F\help1;
new_tangent=new_tangent/sqrt(sum(new_tangent.*new_tangent));
end

