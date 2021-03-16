function [x1,valx1,solution] = bvpsuite2(problem,settings,initProfile,pathfoll)

% Matlab Code BVPSUITE2.0 (beta-version July 2018).
% Authors: Winfried Auzinger, Merlin Fallahpour, Georg Kitzhofer,
% Othmar Koch, Gernot Pulverer, Gustaf Söderlind, Ewa B. Weinmueller and
% Stefan Wurm.
% Institute for Analysis and Scientific Computing, Vienna University of
% Technology, Vienna, Austria.
%
% The package BVPSUITE2.0 has been developed at the Institute for Analysis
% and Scientific Computing, Vienna University of Technology, and can be
% used for the numerical solution of implicit boundary value problems
% (BVPs) in ordinary differential equations (ODEs) as well as eigenvalue
% problems (EVPs), and Index-1 Differential-Algebraic Equations (DAEs). The
% ODE system can have a general implicit form and be of mixed order subject
% to multi-point boundary conditions. Furthermore the interval can be
% finite or semi-infinite. The BVP may also include unknown parameters to
% be calculated together with the unknown solution z, in which case
% additional boundary conditions are required.
% BVPSUITE2.0 is a product of years of research at the Institute for
% Analysis and Scientific computing on the analysis and numerical solution
% of ODEs with time and space singularities. Through the years, the code
% evolved from SBVP1.0, a collocation code for singular ODEs of the first
% order, to BVPSUITE1.1, which can also handle explicit and implicit
% problems of arbitrary order, to BVPSUITE2.0, which improves on usability
% and readability of the code. The code in the package is now structured
% in a simpler modular form.
%
%   Copyright (c) July 2018 Winfried Auzinger
%                           Merlin Fallahpour
%                           Georg Kitzhofer
%                           Othmar Koch
%                           Gernot Pulverer
%                           Gustaf Söderlind
%                           Ewa B. Weinmüller
%                           Stefan Wurm
%                           Vienna University of Technology

global chainrulecoeffs;
chainrulecoeffs = higherchainrule(max(feval(problem,'orders')));
solution=[];

x1 = feval(settings,'mesh');
n = feval_problem(problem,'n');
collMethod = feval(settings,'collMethod');
collPoints = feval(settings,'collPoints');

try
    tmp_pathfoll=feval(problem,'pathfollowing');
    if exist('pathfoll','var')
        tmp_fn=fieldnames(pathfoll);
        for kk=1:numel(tmp_fn)
            tmp_pathfoll.(tmp_fn{kk})=pathfoll.(tmp_fn{kk});
        end
    end
    pathfoll=tmp_pathfoll;
catch
    pathfoll.activate=0;
end

showPlot=0;
debug=0;

if(nargin<3) || isempty(initProfile)
    try
        initProfile = feval_problem(problem,'initProfile');
    catch MException
        disp('Verwende Standardprofil');
        disp(MException);
        initProfile.initialMesh = x1;
        initProfile.initialValues = ones(n,length(x1));
    end
end

try
    initProfile.initialMesh = initProfile.x1;
    initProfile.initialValues = initProfile.valx1;
catch MException
    
end

if(isempty(initProfile.initialValues))
    disp('Verwende Standardprofil');
    disp(MException);
    initProfile.initialMesh = x1;
    initProfile.initialValues = ones(n,length(x1));
end

if(feval(problem,'EVP'))
    initProfile.initialValues(end+1,1)=norm(initProfile.initialValues(:,1),2)^2;
    for j=2:length(initProfile.initialMesh)
        initProfile.initialValues(end,j)=norm(initProfile.initialValues(1:end-1,1:j),2)^2;
    end
    if isfield(initProfile,'parameters')
        initProfile.parameters = [initProfile.parameters; initProfile.lambda];
    else
        initProfile.parameters = initProfile.lambda;
    end
end

%Gitter aus Settings wird auf Intervall aus bvp-file transformiert
interval = feval_problem(problem,'interval');
x1 = (interval(end)-interval(1))/(x1(end)-x1(1))*(x1+x1(1))+interval(1);

%Berechne Kollokationspunkte auf [0,1]
m = collPoints;
if (strcmp(collMethod,'user'))
    [rho] = collPoints;
else
    [rho] = getStandardCollocationPoints(collMethod,m);
end

if pathfoll.activate
    if ~strcmp(pathfoll.startat,'start')
        fprintf('\n<strong>Pathfollowing starting at %s</strong>\n',pathfoll.startat);
    else
        fprintf('\n<strong>Pathfollowing initializing ... </strong>\n');
        try
            pathfoll.initProfile=initProfile;
        catch
            fprintf('\nInitial profile not needed:')
        end
    end
    [p_save,p_exact,tur_pts]=pathfollowing(problem,settings,pathfoll,x1,rho);
    
    speichername=pathfoll.name;
    tmp=1;
    while exist( strcat(speichername,'_',num2str(tmp),'.mat'),'file' ) ==2
        tmp=tmp+1;
    end
    profilname=strcat(speichername,'_',num2str(tmp),'.mat');
    filename=strcat(pathfoll.dir,profilname);
    save(filename,'p_save','p_exact','tur_pts');
    
    fprintf('\nPathfollowing file saved as:\n%s\n',strcat(pathfoll.dir,profilname))
    
    fprintf('\nStarted run at lambda_p = %f.\n',p_save{2,1})
    if ~isempty(tur_pts)
        fprintf('Found <strong>TURNING POINT</strong> inbetween parameter values')
        for ii=1:size(tur_pts,2)
            fprintf('\n %i. <strong>***</strong>   %f and %f and %f  <strong>***</strong>',ii,tur_pts(1,ii),tur_pts(2,ii),tur_pts(3,ii))
        end
        fprintf('.\n')
    else
        fprintf('No turning points were encountered in this run.\n')
    end
    fprintf('Finished run after %i steps at lambda_p = %f.\n',size(p_save,2)-2,p_save{2,end-1})
    
    xfin=p_save{1,end-1}.x1;
    if length(x1)~=length(xfin) || max(abs(x1-xfin))>1e-12
        fprintf('Mesh adaptation was used during the run, ')
        if length(x1)~=length(xfin)
            fprintf('to augment the number of mesh points from %i to %i.\n',length(x1),length(xfin))
        else
            fprintf('to redistribute the %i mesh points fittingly.\n',length(x1))
        end
    elseif (feval(settings,'meshAdaptation'))
        fprintf('Mesh adaptation was not needed in this run.\n')
    end
    
    solution = p_save{1,end-1};
    solution.filename = filename;
end

%Löse lineares Problem
if (feval_problem(problem,'linear') == 1) && isempty(solution)
    dispDebug('Lineares Problem',debug);
    
    
    if(feval(settings,'meshAdaptation'))
        [solution] = meshadaptation(problem,settings,x1,rho);
    elseif(feval(settings,'errorEstimate'))
        [~,~,solution] = solveLinearProblem(problem,x1,rho);
        solution = errorestimate(solution,problem,settings);
    else
        [~,~,solution] = solveLinearProblem(problem,x1,rho);
    end
    
    
    
elseif isempty(solution)
    dispDebug('Nicht-Lineares Problem',debug);
    
    
    if(feval(settings,'meshAdaptation'))
        [solution] = meshadaptation(problem,settings,x1,rho,initProfile);
    elseif(feval(settings,'errorEstimate'))
        [~,~,solution] = solveNonLinearProblem(problem,settings,x1,rho,initProfile);
        solution = errorestimate(solution,problem,settings);
    else
        [~,~,solution] = solveNonLinearProblem(problem,settings,x1,rho,initProfile);
    end
    
end

interval = feval(problem,'interval');
if(interval(2) == Inf)
    if(trafo('splitting'))
        [x2,sortind] = sort(trafo('xi',solution.x1,0,1));
        solution.x1=[solution.x1(1:end-1),x2];
        solution.valx1 = [solution.valx1(1:end/2,1:(end-1)),solution.valx1((end/2+1):end,sortind)];
        [x2tau,sortind] = sort(trafo('xi',solution.x1tau,0,1));
        solution.x1tau=[solution.x1tau(1:end-1),x2tau];
        solution.valx1tau = [solution.valx1tau(1:end/2,1:(end-1)),solution.valx1tau((end/2+1):end,sortind)];
    else
        [solution.x1,sortind] = sort(trafo('xi',solution.x1,0,interval(1)));
        solution.valx1 = solution.valx1(:,sortind);
        [solution.x1tau,sortind] = sort(trafo('xi',solution.x1tau,0,interval(1)));
        solution.valx1tau = solution.valx1tau(:,sortind);
    end
end
x1 = solution.x1;
valx1 = solution.valx1;
if(feval(problem,'EVP')>0)
    valx1 = valx1(1:end-1,:);
    solution.lambda = solution.parameters(end);
    solution.parameters = solution.parameters(1:end-1);
    if isempty(solution.parameters)
        solution.parameters = [];
    end
end
if(showPlot>0)
    plot(x1,valx1);
end

end

