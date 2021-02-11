function [speicher,speicher_exact,tur_pts]=pathfollowing(problem,settings,pathfoll,x1,rho)
% Pathfollowing routine

%% Set options
eval_exact  = 1; % evaluate at the given parameter values
save_ws     = 1; % enable taking step(s) back if desired
dispexact   = 1; % plot at the given parameter values
dispres     = 1; % display the solution and parameter evolution plots
disppredcor = 1; % display with predictor and corrector step
log         = 1; % log some values
displog     = 1; % plot the logged values
dispmdep    = 1; % plot mesh density evolution
display     = 1; % display messages
% For secondary options, press Ctrl+F and search for 'Secondary'

%% Initializing
m           = length(rho);
xfin        = x1;
counter_mc  = 0; % Counter for mesh correction - DISABLED FOR NOW
ordnung     = feval_problem(problem,'orders');
n           = length(ordnung);
interval    = feval(problem,'interval');
infsplit    = 0;
se_nn       = 40;
% In case the right hand side of the interval is infinite, then x1 is splitted in [0,1] and [1,infty) when using splitting interval transformation
selim       = 0;
if interval(2) == Inf
    if interval(1) == 0 && trafo('splitting')
        infsplit = 1;
        n = n/2;
    end
    se_nn = 40;
    % Secondary option
    % The solution evolution variable 'selim' sets up to which value on the x-axis the solutions evolution should be drawn when the solution is computed for a problem posed on a semi-infinite interval
    selim = 10;
end
AbsTol     = feval(settings,'absTolSolver');
RelTol     = feval(settings,'relTolSolver');
minmeshpts = max(feval(settings,'minInitialMesh'),length(feval(settings,'mesh'))); % DISABLED for now

% Set the variables from user input or solver settings
tmp_fn = {'thetaMax','maxCorrSteps','maxSteplengthGrowth','angleMin','PredLengthFactor','CorrLengthGrowth','meshFactorMax'};
tmp_val = cell(1,numel(tmp_fn));
for ii=1:numel(tmp_fn)
    if isfield(pathfoll,tmp_fn{ii})
        tmp_val{ii}=pathfoll.(tmp_fn{ii});
    else
        tmp_val{ii}=feval(settings,tmp_fn{ii});
    end
end
[theta_max,maxCorrSteps,maxSteplengthGrowth,cos_min,pred_lf,corr_lf,meshFactorMax] = deal(tmp_val{:});

try
    data = pathfoll.pathdata;
catch
    msg = 'There is no ret.pathdata defined in the problem definition. Please refer to the manual of bvpsuite2.0 for some examples.';
    error(msg)
end
halve       = 0;
halve_nb    = 0; % times the step-length has been halved
corr_old    = Inf;
min_sl      = 1e-8; % Secondary option: minimal steplength
min_sl_flag = 0; % flag for minimal steplength

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

itnum = 0; % Number of completed steps
skip = 0; % Marker in case tagents angle is too steep
% Current number of steps until prompt
cter_flag = 0;
if isfield(pathfoll,'counter')
    counter = pathfoll.counter;
    if counter==Inf || (isfield(pathfoll,'only_counter') && pathfoll.only_counter)
        cter_flag = 1;
    end
else
    counter = 1;
end

prompt='\nType <strong>enter</strong> to carry out the chosen number of steps (ret.counter/1 when not modified),\n<strong>p</strong> and <strong>enter</strong> to change the maximal length of the predictor step,\n<strong>n</strong> and <strong>enter</strong> to change this chosen number,\n<strong>f</strong> and <strong>enter</strong> to plot&display solutions of computed steps, and \n<strong>s</strong> and <strong>enter</strong> to stop&save the computations\n';
subprompt1='Choose a number of steps: ';
subprompt2='Choose a maximal predictor step length: ';
subprompt3='Choose number of steps to go back: ';
subprompt4='Input a row vector [ . . . ] of steps: ';
subprompt5='Enter one of the options above to proceed: ';
f_size=8; % font size for the plots

% keep the previous interpreter settings and set all to latex
dtiii=get(0, 'DefaultTextInterpreter');
dliii=get(0, 'DefaultLegendInterpreter');
datliii=get(groot, 'DefaultAxesTickLabelInterpreter');
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex')

%% Prepare what is needed for the options
% load the exact values of the parameter that shall be computed
if isfield(pathfoll,'require_exact')
    exact_val   = sort(unique(pathfoll.require_exact));
    if strcmp(pathfoll.startat,'start')
        speicher_exact = {};
    end
else
    exact_val = [];
    speicher_exact = [];
end
if isfield(pathfoll,'max_pred_length')
    max_pred = pathfoll.max_pred_length;
    adapt=0;
else
    max_pred = [];
    adapt=0;
end
if isfield(pathfoll,'pit_stop')
    pit_stop = pathfoll.pit_stop;
else
    pit_stop = [];
end
% If only_exact is activated, then only speicher_exact is computed
if isfield(pathfoll,'only_exact') && ~isempty(exact_val)
    only_exact=pathfoll.only_exact;
    if only_exact
        eval_exact=1;
    end
else
    only_exact=0;
end
% to save the workspace in order to go back 'counter'-number of steps in the iteration
if save_ws
    prompt_ws='\nDo you want to go back some steps?\nIf yes, press <strong>b</strong> and <strong>enter</strong>.\nIf no just <strong>enter</strong>.\n';
    jump_cell = {};
end
% Which plots should be drawn?
if ~eval_exact
    dispexact = 0;
end
if isempty(data) % data should actually never be empty...
    dispres = 0;
end
if ~dispres
    disppredcor = 0;
end
if ~log && displog
    displog = 0;
end
if~(feval(settings,'meshAdaptation'))
    dispmdep = 0;
end
if dispmdep || dispres || displog
    pathfollplot=figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6],'PaperUnits','normalized');
    ax = gca;
    ax.FontSize=f_size;
    clf ;
    pathfollplot.Color = 'White' ;
end
if dispmdep
    % Mesh density evolution plot
    mdep_nn = 50;
    X_mdep  = linspace(0,1,mdep_nn);
    Y_mdep  = [];
    Z_mdep  = [];
    
    % Mesh evolution plot
    pts_nb = 20;
    if length(x1)<=pts_nb
        filler = interval(2)*ones(1,length(x1)-pts_nb);
        Pts_Mat=[x1 filler];
    else
        pts_step = ceil(length(x1)/pts_nb);
        filler = interval(2)*ones(1,pts_nb-length(x1(1:pts_step:end-1)));
        Pts_Mat = [x1(1:pts_step:end-1) filler];
    end
    
    % Call the figures
    if dispres
        path_pos = [0 0.5 1 0.5];
        mdep_pos = [0 0 1 0.5];
    else
        mdep_pos = [0 0 1 1];
    end
elseif dispres
    path_pos = [0 0 1 1];
end
if dispres
    figure(pathfollplot);
end
if log
    logstruct.theta_0         = [] ;
    logstruct.theta_max       = [] ;
    logstruct.steplength      = [] ;
    logstruct.steplength_pred = [pathfoll.steplength] ;
    logstruct.sl_adapt        = [];
    
    logstruct.norm_F          = [];
    logstruct.cpt             = [] ;
    logstruct.orthogonal      = [] ;
    logstruct.max_error       = [];
    logstruct.mesh_length     = [] ;
    
    logstruct.c_s             = Inf  ;
    logstruct.cos_ab          = Inf;
    logstruct.corr_dist       = Inf ;
    logstruct.norm_delta_0    = [] ;
    logstruct.norm_delta_1    = [] ;
    
    logstruct.mesha_halved    = [];
    logstruct.angle_halved    = [];
    logstruct.trm_halved      = [];
    logstruct.pdist_halved    = [];
    logstruct.cdist_halved    = [];
    
    % logstruct fields that are marks in the plot whenever steplength got halved for one of the reasons
    % 1. mesh adaptation was too long
    % 2. angle between consecutive steps was too big
    % 3. trust region method was activated
    % 4. corrector step was too long compared to predictor step
    % 5. corrector step was too long compared to previous corrector step
    fn_halve = fieldnames(logstruct);
    fn_halve = fn_halve(end-4:end);
    val_halve = [ 0 1 0 1 10 ];
    
    if displog
        if dispres && dispmdep
            path_pos = [0.5 0.5 0.5 0.5];
            mdep_pos = [0.5 0 0.5 0.5];
            log_pos  = [0 0 0.5 1];
        elseif dispres
            path_pos = [0.5 0 0.5 1];
            log_pos  = [0 0 0.5 1];
        elseif dispmdep
            mdep_pos = [0.5 0 0.5 1];
            log_pos  = [0 0 0.5 1];
        else
            log_pos  = [0 0 1 1];
        end
    end
else
    logstruct = [];
end
exact_lp = Inf;

%% Preprocess the two cases: Start from initial point or loaded data
% Start at the given initial point on the solution path
if strcmp(pathfoll.startat,'start')
    % Set the starting value for lambda_p
    lambda_p=pathfoll.start;
    
    % Before computing a predictor corrector step, the initial solution needs to be computed first
    itnum=itnum-1;
    
    % If the predictor and corrector steps should be drawn
    if disppredcor
        val = {};
    else
        val = [];
    end
    
    % Inital steplength and initalize struct for predictor, used to carry over information from the pathfollowing routine into the bvpsuite2.0 routine
    steplength = pathfoll.steplength;
    predictor  = struct;
    
    % Save the important values in predictor, that matter for the computation
    predictor.steplength = steplength;
    predictor.infsplit   = infsplit;
    predictor.lpfix      = 0;
    predictor.x1         = x1;
    predictor.lambda_p_0 = lambda_p;
    predictor.skip       = 0;
    
    % Also some solver settings (needed in PredCorrStrat and meshadaptation.m) and display setting
    predictor.maxCorrSteps        = maxCorrSteps;
    predictor.maxSteplengthGrowth = maxSteplengthGrowth;
    predictor.meshFactorMax       = meshFactorMax;
    predictor.display             = display;
    
    % Starting point data for the Newton iteration in case of nonlinear problem
    a_p=pathfoll.initProfile;
    
    % Search for the next given parameter value, for which the value should be computed, if required
    if ~isempty(exact_val)
        if steplength >0
            exact_lp = find(exact_val>=lambda_p,1);
        else
            exact_lp = find(exact_val<=lambda_p,1,'last');
        end
    end
    
    % Define some values needing to exist for the function call of Postprocess later, they do not have any specific values as of now
    tangent_new=0; theta_0=0; delta_0=0;
else
    try
        filename=strcat(pathfoll.dir,pathfoll.startat);
        load(filename);
    catch MException
        error(MException.message);
    end
    if ~only_exact
        speicher_exact = p_exact;
    else
        speicher_exact = {};
    end
    speicher = p_save(:,1:end-1); % savestruct column is not needed
    
    % Define x1 and a_p
    sol       = speicher{1,end-1}; % sol is the penultimate solution
    x1        = sol.x1;
    a_0       = [ sol.coeff ; sol.parameters ; sol.lambda_p ] ;
    % Define sol, predictor, xfin and a_c
    sol       = speicher{1,end}; % sol is the last computed solution
    predictor = sol.predictor;
    a_p       = a_0 + predictor.steplength.*predictor.tangent';
    xfin      = sol.x1;
    a_c       = [ sol.coeff ; sol.parameters ; sol.lambda_p ] ;
    
    savestruct = p_save{1,end};
    if log
        logstruct   = p_save{2,end};
    end
    
    jump_cell   = p_save{3,end};
    
    % Reassign all the variables from savestruct
    delta_0         = savestruct.delta_0;
    theta_0         = savestruct.theta_0;
    itnum           = savestruct.itnum;
    val             = savestruct.val;
    exact_lp        = savestruct.exact_lp;
    pred_tmp        = savestruct.pred_tmp;
    corr_old        = savestruct.corr_old;
    if isfield(savestruct,'halve_nb')
        halve_nb    = savestruct.halve_nb;
    else
        halve_nb    = 5;
    end
    if ~exist('pit_stop','var')
        pit_stop    = savestruct.pit_stop;
    end
    if dispmdep
        try
            Pts_Mat = savestruct.Pts_Mat;
            Y_mdep  = savestruct.y_mdep;
            Z_mdep  = savestruct.z_mdep;
        catch
            Y_mdep  = [];
        end
    end
    
    % If data has been changed, then recompute for the new data
    if max(data(xfin,a_c,ordnung,rho)~=speicher{3,end})
        [speicher,speicher_exact,val,corr_old] = correct_data(speicher,speicher_exact,val,data,ordnung,rho);
    end
    
    steplength = predictor.steplength;
    tangent_new = predictor.tangent';
    % Compute eval_exact if necessary
    if ~isempty(exact_val) && ~only_exact
        eval_exact = 1;
        % if sol.lambda_p<exact_val(1) and step-length<0 or
        tmp_bool1 = sol.lambda_p<exact_val(1) && steplength<0;
        % if sol.lambda_p>exact_val(end) and step-length>0
        tmp_bool2 = exact_val(end)<sol.lambda_p && steplength>0;
        % then for now, eval_exact is set to 0
        if tmp_bool1 || tmp_bool2
            eval_exact=0;
        end
    end
    
    % Update solver settings in case it was changed by user when reloading the run
    predictor.maxCorrSteps = maxCorrSteps;
    predictor.maxSteplengthGrowth = maxSteplengthGrowth;
    predictor.meshFactorMax = meshFactorMax;
    
    if ~only_exact
        % Prepare for the next step
        [a_0,a_p,tangent_new,steplength,predictor,delta_0,nda,nsd,theta_0,logstruct,eval_exact,exact_lp] = PredCorrStrat(problem,settings,predictor,a_c,a_p,x1,xfin,sol,tangent_new,ordnung,rho,psival,psi,steplength,exact_val,exact_lp,delta_0,theta_0,theta_max,log,logstruct,eval_exact,itnum,halve_nb,min_sl);
    else
        eval_exact = 1;
        exact_lp = Inf;
    end
    
    % Accept the new mesh for the next step
    x1 = xfin;
    
    lambda_p = speicher{2,1};
    
    itnum=itnum+1;
end

%% Pathfollowing routine
if itnum==-1 % Initial solution is computed first
    jj = 0;
else % Initial solution was already computed
    jj = 1;
end
% Secondary option: 
% After jj_max steps, the routine will stop & save. This number may be changed. This is a safety measure for the case where 'counter' is set to Inf.
jj_max = 1e4;
if counter==Inf
    counter=jj_max;
end

% Main loop: is only exited as soon as jj==counter and then 's' is chosen in the user prompt, either by the user or automatically.
while 1
    % Prechecks of the step-length
    if itnum~=-1 && ~only_exact
        % Check whether the step-length is shorter than the minimal steplength
        if abs(steplength)<min_sl
            min_sl_flag=1; % Set flag to stop later
            steplength=sign(steplength)*min_sl;
            predictor.steplength=steplength;
            a_p = a_0 + steplength.*tangent_new;
            predictor.lambda_p_p=a_p(end);
            if display
                fprintf('\n  The steplength was shorter than the minimal steplength %1.1e!\n  It has been augmented to the minimal steplength.\n',min_sl)
            end
        end
        % Check whether the user-defined maximal step-length is overreached, if it is, adapt the predictor to a closer step
        data_tmp=data(x1,a_p,ordnung,rho);
        pred_tmp=sqrt((a_p(end)-speicher{2,end})^2+(max(abs(speicher{3,end}-data_tmp)))^2);
        if ~isempty(max_pred) && pred_tmp>max_pred % pred_tmp~=max_pred
            div_tmp=max_pred/pred_tmp;
            steplength=steplength*div_tmp;
            predictor.steplength=steplength;
            a_p = a_0 + steplength.*tangent_new;
            predictor.lambda_p_p=a_p(end);
            
            pred_tmp = max_pred;
            
            if display
                fprintf('\n  Predictor step was too long. Steplength was reduced to %f.\n',steplength)
            end
            
            if log
                % save whether the step-length was dimished or not
                adapt = 1;
            end
        end
    end
    
    if only_exact
        fprintf('\nComputing the approximate solutions for lambda_p=[ require_exact ]:\n')
    elseif itnum~=-1 && display
        fprintf('\n  Starting at predictor value %f:\n\n',a_p(end))
    elseif display
        fprintf('\n  Starting at predictor value %f:\n\n',lambda_p)
    end
    
    % Corrector Step
    if halve==0
        cpt = cputime;
    end
    if log && itnum~=-1 && ~only_exact
        if skip==0
            for ii=1:length(fn_halve)
                logstruct.(fn_halve{ii}) = [ logstruct.(fn_halve{ii}) val_halve(ii)/(halve==ii) ];
            end
        else
            logstruct.(fn_halve{halve})(end) = val_halve(halve);
        end
    end
    if only_exact
        itnum = itnum-1;
    elseif (feval(settings,'meshAdaptation'))
        [sol,~,halve] = meshadaptation(problem,settings,x1,rho,a_p,predictor);
        % If points have been added in the mesh, then reset the counter, otherwise wait for ? number of iterations in which the mesh does not change, to halve the mesh length if needed
        % ..................... DISABLED FOR NOW .....................
        if length(x1)>length(xfin)
            counter_mc = 0;
        else
            counter_mc = counter_mc + 1;
        end
    elseif(feval(settings,'errorEstimate'))
        if ~isfield(predictor,'tangent') && feval(problem,'linear')
            [~,~,sol]=solveLinearProblem(problem,x1,rho);
        else
            [~,~,sol,halve] = solveNonLinearProblem(problem,settings,x1,rho,a_p,1,predictor);
        end
        if itnum==-1
            sol.lambda_p = lambda_p;
        end
        if halve==0
            initCoeff = [sol.coeff', sol.parameters, sol.lambda_p ]';
            sol.initCoeff = initCoeff;
            [sol,~,halve] = errorestimate(sol,problem,settings,predictor);
        end
    else
        if ~isfield(predictor,'tangent_new') && feval(problem,'linear')
            [~,~,sol]=solveLinearProblem(problem,x1,rho);
        else
            [~,~,sol,halve] = solveNonLinearProblem(problem,settings,x1,rho,a_p,1,predictor);
        end
        if log
            sol.errest = Inf;
        end
    end
    if itnum==-1 && ~only_exact && isstruct(sol)
        sol.lambda_p = lambda_p;
    end
    
    % If the number of mesh points does not need halving, then compute the angle between the tangents and the lengths of the predictor and corrector steps
    if halve==0 && ~only_exact
        % Save the predictor for possible future use
        sol.predictor=predictor;
        
        % Save the result of the corrector step on the solution path and the new mesh xfin
        a_c = [ sol.coeff ; sol.parameters ; sol.lambda_p ];
        predictor.lambda_p_0 = a_c(end);
        xfin = sol.x1;
        
        while itnum~=-1
            % Check whether the corrector step is short enough compared to the predictor step
            speicher_tmp=speicher;
            speicher_tmp{1,end+1}=sol;
            speicher_tmp{2,end}=a_c(end);
            speicher_tmp{3,end}=data(xfin,a_c,ordnung,rho);
            corr_tmp=sqrt((a_c(end)-a_p(end))^2+(max(abs(speicher_tmp{3,end}-data_tmp)))^2);
            factor1_tmp=corr_tmp/pred_tmp;
            factor2_tmp=corr_tmp/corr_old;
            if display
                fprintf('\n  Corrector step is %f-times as long as predictor step.\n',factor1_tmp)
            end
            if itnum~=0
                if display
                    fprintf('  Corrector step is %f-times as long as the previous corrector step.\n',factor2_tmp)
                end
            end
            if factor1_tmp>1/pred_lf
                halve=4;
                if display
                    fprintf('\n  Corrector step is at most allowed to be %f-times as long as predictor step!\n',1/pred_lf)
                end
                if ~min_sl_flag
                    break
                end
            end
            if factor2_tmp>corr_lf
                halve=5;
                if display
                    fprintf('\n  Corrector step is at most allowed to be %f-times as long as previous corrector step!\n',corr_lf)
                end
                if ~min_sl_flag
                    break
                end
            end
            
            if itnum~=0
                % Compute the cosine between the last computed step and the next tangent
                if (feval(settings,'meshAdaptation')) && (length(x1)~=length(xfin) || max(abs(x1-xfin))>1e-12)
                    initP.initialMesh=xfin;
                    initP.parameters = sol.parameters;
                    
                    initP.initialValues = coeffToValues(tangent_new, x1,ordnung,rho,xfin);
                    initP = initial_coefficients(problem,xfin,initP,rho,0);
                    tang_a = [ initP.initialCoeff ; tangent_new(end) ] ;
                else
                    tang_a = tangent_new;
                end
                jac_F = functionFDF( 'DF', problem ,a_c,xfin,psival,psi,rho,[]);
                tang_b = tangente_berechnen(jac_F);
                a_p_b = a_c + sign(steplength)*sign(tang_a.'*tang_b)*tang_b;
                data_tang=data(xfin,a_p_b,ordnung,rho);
                nn_data = length(data_tang);
                a=[(speicher_tmp{2,end}-speicher_tmp{2,end-1})*ones(nn_data,1) speicher_tmp{3,end}-speicher_tmp{3,end-1}];
                b=[(a_p_b(end)-speicher_tmp{2,end})*ones(nn_data,1) data_tang-speicher_tmp{3,end}];
                tmp_ab=zeros(1,nn_data);
                for ii=1:nn_data
                    n_a=a(ii,:)*a(ii,:)';
                    n_b=b(ii,:)*b(ii,:)';
                    tmp_ab(ii)=a(ii,:)*b(ii,:)'/sqrt(n_a*n_b);
                end
                cos_ab=min(tmp_ab);
                if display
                    fprintf('  Cosine of the angle between current step and next tangent: %1.1e\n',cos_ab)
                end
                % If the cosine smaller than the user set cosine, then the steplength is halved
                if cos_ab<cos_min
                    halve=2;
                    if display
                        fprintf('\n  The cosine is smaller than the minimal allowed cosine %1.1e\n',cos_min)
                    end
                end
            end
            break
        end
    end
    
    if itnum~=-1 && halve~=0 && ~min_sl_flag
        steplength = steplength/2;
        predictor.steplength=steplength;
        a_p = a_0 + steplength.*tangent_new;
        predictor.lambda_p_p=a_p(end);
        
        if display
            fprintf('  Steplength was halved from %f to %f!\n',2*steplength,steplength)
        end
        
        itnum=itnum-1; % itnum is augmented by 1 at the end of each loop
        skip = skip+1;
        predictor.skip = skip;
        counter = counter + 1;
    else % Proceed to finishing the step
        if ~only_exact
            halve_nb = skip;
            halve = 0;
            counter = counter-skip;
            jj = jj-skip;
            skip=0;
            predictor.skip = 0;
            if min_sl_flag && ~isstruct(sol) % pathfollowing is stuck
                msg1='The steplength is reduced too much! The pathfollowing routine has trouble following the path with the current solver settings! The current mininmal steplength (min_sl) is ';
                msg2='. To change this value look for min_sl at the beginning of the file pathfollowing.m in the bvpsuite2.0 package. Otherwise adapting the solver settings, especially using mesh adaptation (MA) and choosing the MA tolerances only slightly looser than the solver tolerances (i.e. solver tol: 10^-6, MA tol: 10^-4), may help.';
                msg=[msg1 num2str(min_sl) msg2];
                error(msg)
            end
            
            % Save the pathdata of interest in the variable speicher
            if itnum==-1
                speicher{1,1}=sol;
                speicher{2,end}=a_c(end);
                speicher{3,end}=data(xfin,a_c,ordnung,rho);
            else
                speicher=speicher_tmp;
                corr_old=corr_tmp;
            end
        end
        
        % Check if one of the required values of lambda_p has been passed if yes, compute the solution in this point and save the solution in speicher_exact
        if itnum~=-1 && only_exact==1 % only speicher_exact is computed in this run
            % Restart at the beginning & find the next exact value to compute
            lp_a0 = cell2mat(speicher(2,:));
            % Find the value exact_val(exact_lp) that is closest to the starting value of the parameter lambda_p
            for ii=2:length(lp_a0)
                tmp_sign = sign(lp_a0(ii)-lp_a0(ii-1));
                if tmp_sign==1
                    tmp_line = 1:length(exact_val);
                else
                    tmp_line = length(exact_val):-1:1;
                end
                for kk = tmp_line
                    tmp_bool1=tmp_sign*lp_a0(ii-1)<=tmp_sign*exact_val(kk) && tmp_sign*exact_val(kk)<tmp_sign*lp_a0(ii);
                    % If the required exact value is in-between two values of the path, then save ii and kk and break
                    if tmp_bool1
                        exact_lp=kk; % is the entry number in exact_val
                        cter_a0=ii-1; % is the step number in the path leading to a_0
                        sol = speicher{1,cter_a0+1}; % is the solution in the path leading to a_c exact_val(exact_lp) is in-between a_0(end) and a_c(end)
                        break
                    end
                end
                if exact_lp~=Inf
                    break
                end
            end
            % If the for-loop did not break
            if exact_lp==Inf
                msg='None of the required values in ret.require_exact is within the range of the computed path!';
                error(msg)
            end
        elseif itnum~=-1 && ~isempty(exact_val) % Compute lp_a0 and cter_a0 to check whether a required value has been passed or not
            lp_a0 = a_0(end);
            cter_a0 = 1;
        end
        % Check whether or not one of the values that should be computed at exactly, was already passed, if yes, compute at the value(s) that have been passed and are required
        if itnum~=-1 && ~isempty(exact_val) && eval_exact && abs(lp_a0(cter_a0)-sol.lambda_p)> abs(lp_a0(cter_a0)-exact_val(exact_lp))
            % More than one value can be in the interval that was computed, thus a while statement here
            while eval_exact && abs(lp_a0(cter_a0)-sol.lambda_p)> abs(lp_a0(cter_a0)-exact_val(exact_lp))
                if display
                    fprintf('\n<strong>Exact solution</strong> at %f computation:\n\n',exact_val(exact_lp))
                end
                
                % If by chance the value of lambda in the last computation was exactly the required point, then no computaiton is needed, otherwise it is
                if abs(lp_a0(cter_a0)-exact_val(exact_lp)) > 0
                    % Define sol_tmp and a_0_tmp whether only_exact is 1 or 0
                    if only_exact
                        sol_tmp = speicher{1,cter_a0};
                        a_0_tmp = [ sol_tmp.coeff ; sol_tmp.parameters ; sol_tmp.lambda_p ] ;
                    else
                        sol_tmp = sol;
                        a_0_tmp = a_0;
                    end
                    predictor = sol_tmp.predictor;
                    tangent = predictor.tangent';
                    sl_tmp = predictor.steplength;
                    predictor.lpfix = 1;
                    predictor.lambda_p_0 = exact_val(exact_lp);
                    a_p_exact = a_0_tmp + (exact_val(exact_lp)-lp_a0(cter_a0))/tangent(end)*tangent;
                    if(feval(settings,'meshAdaptation'))
                        [sol_exact] = meshadaptation(problem,settings,x1,rho,a_p_exact,predictor);
                    elseif(feval(settings,'errorEstimate'))
                        [~,~,sol_exact] = solveNonLinearProblem(problem,settings,x1,rho,a_p_exact,1,predictor);
                        initCoeff = [sol_exact.coeff', sol_exact.parameters, sol_exact.lambda_p ]';
                        sol_exact.initCoeff = initCoeff;
                        sol_exact = errorestimate(sol_exact,problem,settings,predictor);
                    else
                        [~,~,sol_exact] = solveNonLinearProblem(problem,settings,x1,rho,a_p_exact,1,predictor);
                    end
                    predictor.lpfix = 0;
                else % Solution was already computed
                    if ~only_exact
                        sol_exact = speicher{1,end-1};
                    else
                        sol_exact = speicher{1,cter_a0};
                    end
                    % Load the steplength of the last computed step
                    predictor = sol.predictor;
                    sl_tmp = predictor.steplength;
                end
                speicher_exact{1,end+1} = sol_exact;
                tmp_vct = [ sol_exact.coeff ; sol_exact.parameters ; sol_exact.lambda_p ] ;
                speicher_exact{2,end} = tmp_vct(end);
                speicher_exact{3,end} = data(sol_exact.x1,tmp_vct,ordnung,rho);
                
                if display
                    fprintf('\n<strong>Exact solution</strong> computed.\n')
                end
                predictor.lambda_p_0 = a_c(end); % In case it was changed
                
                if dispexact
                    plot_step(sol_exact,f_size,infsplit,n,selim)
                end
                
                % take the next value of the required values for next time
                if ~only_exact
                    exact_lp = exact_lp + sign(sl_tmp);
                    if exact_lp==0 || exact_lp>length(exact_val)
                        exact_lp = exact_lp - sign(sl_tmp);
                        eval_exact=0;
                    end
                else
                    % exact_lp = exact_lp+1?
                    tmp_bool1= sl_tmp>0 && exact_lp<length(exact_val) && (abs(lp_a0(cter_a0)-sol.lambda_p) > abs(lp_a0(cter_a0)-exact_val(exact_lp+1)));
                    if tmp_bool1
                        exact_lp = exact_lp+1;
                    end
                    % or exact_lp = exact_lp-1?
                    tmp_bool2= sl_tmp<0 && exact_lp>1 && (abs(lp_a0(cter_a0)-sol.lambda_p) > abs(lp_a0(cter_a0)-exact_val(exact_lp-1)));
                    if tmp_bool2
                        exact_lp = exact_lp-1;
                    end
                    % if neither one, then search for the next exact_val(exact_lp) along the path
                    if ~(tmp_bool1 || tmp_bool2)
                        exact_lp = Inf;
                        for ii=cter_a0+2:length(lp_a0)
                            for kk=1:length(exact_val)
                                tmp_bool1=exact_val(kk)<lp_a0(ii) && exact_val(kk)>lp_a0(ii-1);
                                tmp_bool2=exact_val(kk)>lp_a0(ii) && exact_val(kk)<lp_a0(ii-1);
                                % If the required exact value is in-between two values of the path, then save ii and kk and break
                                if tmp_bool1 || tmp_bool2
                                    exact_lp=kk;
                                    cter_a0=ii-1;
                                    sol = speicher{1,cter_a0+1};
                                    break
                                end
                            end
                            if exact_lp~=Inf
                                break
                            end
                        end
                        if exact_lp==Inf
                            if display
                                fprintf('\nAll the required exact solutions were computed and saved in p_exact.\n')
                            end
                            % Compute the next required value from here, which is not used in this run, since only_exact would lead straight to save&return
                            if speicher{1,end}.predictor.steplength >0
                                exact_lp = find(exact_val>=sol.lambda_p,1);
                                if isempty(exact_lp)
                                    exact_val=0;
                                    exact_lp=length(exact_val);
                                end
                            else
                                exact_lp = find(exact_val<=sol.lambda_p,1,'last');
                                if isempty(exact_lp)
                                    exact_val=0;
                                    exact_lp=1;
                                end
                            end
                            break
                        end
                    end
                end
            end
        end
        
        % When only_exact is activated, then stop here
        if only_exact==1
            jj=counter;
        else
            % Compute the predictor and corrector steps, if they are to be plotted
            if itnum~=-1 && isa(val,'cell')
                val{1,end+1} = a_p(end);
                val{2,end} = data(x1,a_p,ordnung,rho);
            end
            
            if log
                err = max(max(abs(sol.errest)));
                logstruct.max_error = [ logstruct.max_error err];
                logstruct.mesh_length = [logstruct.mesh_length [length(x1) ; length(xfin)]];
                logstruct.cpt = [logstruct.cpt cputime-cpt ] ;
                
                if itnum~=-1
                    logstruct.theta_0 = [ logstruct.theta_0 theta_0 ] ;
                    logstruct.theta_max = [ logstruct.theta_max theta_max ];
                    logstruct.steplength = [ logstruct.steplength steplength] ;
                    
                    logstruct.norm_delta_0 = [logstruct.norm_delta_0 nda ] ;
                    logstruct.norm_delta_1 = [logstruct.norm_delta_1 nsd ] ;
                    
                    if adapt
                        logstruct.sl_adapt = [ logstruct.sl_adapt 0];
                    else
                        logstruct.sl_adapt = [ logstruct.sl_adapt Inf];
                    end
                    
                    if length(a_p)~=length(a_c) || max(abs(xfin-x1))>1e-12 % mesh got adapted
                        p_ma.lpfix = 1;
                        F_tmp1 = functionFDF( 'F', problem ,a_p ,x1,psival,psi,rho,p_ma) ;
                        F_tmp2 = functionFDF( 'F', problem ,a_c ,xfin,psival,psi,rho,p_ma) ;
                    else
                        F_tmp1 = functionFDF( 'F', problem ,a_p ,x1,psival,psi,rho,predictor) ;
                        F_tmp2 = functionFDF( 'F', problem ,a_c ,xfin,psival,psi,rho,predictor) ;
                    end
                    logstruct.norm_F = [logstruct.norm_F [ norm(F_tmp1) ; norm(F_tmp2) ] ] ;
                    logstruct.orthogonal = [logstruct.orthogonal [abs(F_tmp1(end)) ; abs(F_tmp2(end)) ] ] ;
                    
                    if itnum~=0
                        logstruct.cos_ab = [logstruct.cos_ab cos_ab];
                    end
                end
            end
            
            % Compute the data for the plot of the mesh density evolution plot
            if dispmdep && itnum~=-1
                x_scaled = (xfin-xfin(1))/(xfin(end)-xfin(1)); % Scaling
                dens =1./(diff(x_scaled)*(length(x_scaled)-1)); % Mesh density
                x_int = x_scaled(1:end-1)+diff(x_scaled)/2;
                new_dens = pchip(x_int,dens,X_mdep);
                Y_mdep = [Y_mdep, itnum+1];
                Z_mdep = [Z_mdep; new_dens ];
                
                % Mesh evolution
                if length(xfin)<=pts_nb
                    filler = interval(2)*ones(1,length(xfin)-pts_nb);
                    Pts_Mat=[Pts_Mat; xfin filler];
                else
                    pts_step = ceil(length(xfin)/pts_nb);
                    filler = interval(2)*ones(1,pts_nb-length(xfin(1:pts_step:end-1)));
                    Pts_Mat = [Pts_Mat; xfin(1:pts_step:end-1) filler];
                end
            end
            
            % If the number of mesh points has been changed from the original number, then divide the number of mesh points by 2 if possible, otherwise reduce or keep the number of mesh points to minmeshpts
            % ..................... DISABLED FOR NOW .....................
            if counter_mc == 10^10 % Secondary option: Set a lower number
                if length(xfin)/2 > minmeshpts && mod(length(xfin),2)==1
                    xfin = xfin(1:2:end);
                elseif length(xfin)/2 > minmeshpts && mod(length(xfin),2)==0
                    xfin = [xfin(1:2:end),xfin(end)];
                elseif length(xfin)>minmeshpts
                    help = xfin(1:floor(length(xfin)/minmeshpts):end);
                    if help(end) == xfin(end)
                        xfin = help;
                    else
                        xfin = [help xfin(end)];
                    end
                end
                counter_mc = 0;
            end
            
            if halve_nb~=0
                fprintf('\nThe step-length was halved %i times in this step.', halve_nb)
            end
            
            if itnum==-1
                fprintf('\n<strong>Initial solution computed</strong> for value %f.\n', a_c(end) )
            elseif isempty(max_pred)
                fprintf('\n<strong>Step %i complete</strong> for value %f. (Step %i/%i)\n', itnum+1, a_c(end), jj, counter )
            else
                fprintf('\n<strong>Step %i complete</strong> for value %f with pl=%1.3f and max_pl=%1.3f. (Step %i/%i)\n', itnum+1, a_c(end), pred_tmp, max_pred, jj, counter )
            end
            
            % If min_sl was reached, then the user can choose whether to continue or not
            if min_sl_flag
                min_sl_flag = 0;
                halve = 0;
                jj=counter;
                if display
                    fprintf('\nThe minimal steplength %1.1e was reached!\nChoose whether to continue or stop here.\n',min_sl)
                end
            end
            
            % Check whether the pit_stop values have been reached
            if itnum~=-1 && ~isempty(pit_stop)
                tmp_bool1= a_0(end)<pit_stop(1) & pit_stop(1)<=a_c(end);
                tmp_bool2= a_0(end)>pit_stop(1) & pit_stop(1)>=a_c(end);
                if tmp_bool1 || tmp_bool2
                    jj=counter;
                    if display
                        fprintf('\nParameter value %f reached!!\nTime for a <strong>pit stop</strong>.\n',pit_stop(1))
                    end
                    pit_stop(1)=Inf;
                end
                tmp_bool1= speicher{3,end-1}<pit_stop(2) & pit_stop(2)<=speicher{3,end};
                tmp_bool2= speicher{3,end-1}>pit_stop(2) & pit_stop(2)>=speicher{3,end};
                if max(tmp_bool1) || max(tmp_bool2)
                    jj=counter;
                    if display
                        fprintf('\nData value %f reached!!\nTime for a <strong>pit stop</strong>.\n',pit_stop(2))
                    end
                    pit_stop(2)=Inf;
                end
            end
            
            % Save the values to jump back to this state from a later step
            if itnum>0
                savestruct = struct;
                
                savestruct.delta_0     = delta_0;
                savestruct.theta_0     = theta_0;
                savestruct.itnum       = itnum;
                savestruct.val         = val;
                savestruct.exact_lp    = exact_lp;
                savestruct.pred_tmp    = pred_tmp;
                savestruct.corr_old    = corr_old;
                savestruct.pit_stop    = pit_stop;
                savestruct.halve_nb    = halve_nb;
                if dispmdep
                    savestruct.Pts_Mat = Pts_Mat;
                    savestruct.y_mdep  = Y_mdep;
                    savestruct.z_mdep  = Z_mdep;
                end
                
                jump_cell{end+1}=savestruct;
            end
        end
    end
    
    % when the number of steps has been carried out
    if jj==counter
        jj=0;
        
        % Call the plot
        clf(pathfollplot) ;
        pathfollplot.Color = 'White' ;
        
        % Draw the path
        if dispres
            % Values of lambda_p
            X_se = cell2mat(speicher(2,:));
            
            % Mesh points
            if infsplit % Mesh was transformed from [0,infty) to [0,1]
                % In [0,1], only plot every 'DD'th value. This number is a divisor of se_nn.
                KK = 1:30;
                DD = KK(rem(se_nn,KK)==0);
                DD = min(DD(DD>2));
                Y_se = [ linspace(0,1,se_nn/DD+1), fliplr(1./linspace(1/se_nn,1-1/se_nn,se_nn-1)) ];
                selimpt = find(Y_se >= selim,1);
                Y_se = Y_se(1:selimpt);
            else
                Y_se = linspace(interval(1),interval(2),se_nn);
            end
            
            % Function values at X_se and Y_se
            Z_se = cell(1,n);
            for tt=1:n
                Z_se{tt} = zeros(length(X_se),length(Y_se));
                if infsplit
                    for kk=1:itnum+2
                        x1tmp = speicher{1,kk}.x1;
                        coefftmp = speicher{1,kk}.coeff;
                        help = coeffToValues( coefftmp, x1tmp,ordnung,rho,linspace(0,1,se_nn+1));
                        help1 = help(tt,1:DD:end); % values on [0,1]
                        help2 = help(n+tt,end-1:-1:size(help,2)+se_nn/DD+1-selimpt); % values on (1,selim]
                        Z_se{tt}(kk,:)=[help1,help2];
                    end
                else
                    for kk=1:itnum+2
                        help = coeffToValues( speicher{1,kk}.coeff, speicher{1,kk}.x1,ordnung,rho,Y_se);
                        Z_se{tt}(kk,:) = help(tt,:);
                    end
                end
                
            end
            
            % Draw the plots
            if disppredcor && itnum~=-1
                figdataplot(itnum+2,speicher,X_se,Y_se,Z_se,pathfollplot,path_pos,val,problem)
            elseif dispres && itnum~=-1
                figdataplot(itnum+2,speicher,X_se,Y_se,Z_se,pathfollplot,path_pos,zeros(length(speicher{3,1}),1),problem)
            end
        end
        
        if log && itnum~=-1 && displog % Draw log plot
            if only_exact
                p_ma.lpfix = 1;
                F_tmp2 = functionFDF( 'F', problem ,a_c ,xfin,psival,psi,rho,p_ma) ;
            end
            figlogdataplot(itnum+1,logstruct,pathfollplot,log_pos,AbsTol+max(F_tmp2)*RelTol)
        end
        
        if dispmdep && size(Z_mdep,1)>1 % Draw mesh and mesh density evolution plot
            figure(pathfollplot);
            
            % Mesh density evolution plot
            axes( 'Position', [ mdep_pos(1)+mdep_pos(3)/10, mdep_pos(2)+mdep_pos(4)/10, 5*mdep_pos(3)/10, 7*mdep_pos(4)/10] ) ;
            surf(X_mdep,Y_mdep,Z_mdep)
            zl = zlim;
            zlim([0 zl(2)])
            title('mesh density')
            ax = gca;
            ax.FontSize=f_size;
            ax.Title.FontWeight='normal';
            
            % Mesh evolution plot
            axes( 'Position', [ mdep_pos(1)+7*mdep_pos(3)/10, mdep_pos(2)+mdep_pos(4)/10, 2*mdep_pos(3)/10, 7*mdep_pos(4)/10] ) ;
            hold on
            for ii = Y_mdep(1)-1:itnum+1
                % Rescale in case meshadaptation was not used from the start but just from a loaded point
                oo=ii-Y_mdep(1)+1;
                plot(Pts_Mat(oo+1,:),ii,'k.','LineWidth',1.5);
            end
            title('mesh evolution')
            ax = gca;
            ax.FontSize=f_size;
            ax.Title.FontWeight='normal';
        end
        
        % If enabled, go back a step if wished
        if itnum>1 && save_ws && ~cter_flag && ~only_exact
            x_inp = input(prompt_ws,'s');
            switch x_inp
                case 'b' % Go back a number of steps
                    % Let the user input a number of steps to go back
                    while 1
                        y_inp=input(subprompt3);
                        if y_inp>=0 && y_inp>=0 && y_inp<=length(jump_cell)-2 && rem(y_inp,floor(y_inp))==0
                            back_cter=y_inp;
                            break
                        end
                        fprintf('Input should be a natural number between 0 and %i.\n',length(jump_cell)-2)
                    end
                    
                    % Adjust jump_cell and define savestruct
                    jump_cell = jump_cell(1:end-back_cter);
                    savestruct = jump_cell{end};
                    
                    % Load from savestruct
                    delta_0         = savestruct.delta_0;
                    theta_0         = savestruct.theta_0;
                    itnum           = savestruct.itnum;
                    val             = savestruct.val;
                    exact_lp        = savestruct.exact_lp;
                    pred_tmp        = savestruct.pred_tmp;
                    corr_old        = savestruct.corr_old;
                    pit_stop        = savestruct.pit_stop;
                    if isfield(savestruct,'halve_nb')
                        halve_nb    = savestruct.halve_nb;
                    else
                        halve_nb    = 5;
                    end
                    if dispmdep
                        try
                            Pts_Mat = savestruct.Pts_Mat;
                            Y_mdep  = savestruct.y_mdep;
                            Z_mdep  = savestruct.z_mdep;
                        catch
                            Y_mdep  = [];
                        end
                    end
                    
                    % Adjust speicher
                    speicher = speicher(:,1:end-back_cter);
                    
                    % Define x1 and a_p
                    sol       = speicher{1,end-1}; % sol is the penultimate solution
                    x1        = sol.x1;
                    a_0       = [ sol.coeff ; sol.parameters ; sol.lambda_p ] ;
                    % Define sol, predictor, xfin and a_c
                    sol       = speicher{1,end}; % sol is the last computed solution
                    predictor = sol.predictor;
                    a_p       = a_0 + predictor.steplength.*predictor.tangent';
                    xfin      = sol.x1;
                    a_c       = [ sol.coeff ; sol.parameters ; sol.lambda_p ] ;
                    
                    % Readjust all the fields in logstruct
                    if log && back_cter~=0
                        tmp_fn = fieldnames(logstruct);
                        for kk=1:numel(tmp_fn)
                            logstruct.(tmp_fn{kk})=logstruct.(tmp_fn{kk})(:,1:end-back_cter);
                        end
                    end
                    
                    predictor.maxCorrSteps        = maxCorrSteps;
                    predictor.maxSteplengthGrowth = maxSteplengthGrowth;
                    predictor.meshFactorMax       = meshFactorMax;
                    
                    % Load steplength and tangent_new from predictor
                    steplength = predictor.steplength;
                    tangent_new = predictor.tangent';
                    % Compute eval_exact if necessary
                    if ~isempty(exact_val)
                        tmp_bool1= exact_lp==1 && exact_val(exact_lp)<sol.lambda_p && steplength>0;
                        tmp_bool2= exact_lp==length(exact_val) && exact_val(exact_lp)>sol.lambda_p && steplength<0;
                        if tmp_bool1 || tmp_bool2
                            eval_exact=1;
                        else
                            eval_exact=0;
                        end
                    end
            end
        end
        
        % If only_exact is activated or the counter was set to Inf, then speicher and speicher_exact are saved, otherwise it depends on the user input
        if only_exact || cter_flag
            x_inp = 's';
        else
            x_inp = input(prompt,'s');
            while x_inp == 'f'
                while 1
                    y_inp=input(subprompt4);
                    if ~isempty(y_inp) && min(y_inp>=0) && min(y_inp<=itnum+1) && max(rem(y_inp+1,floor(y_inp+1)))==0
                        break
                    end
                    fprintf('Inputs must be natural numbers between 0 and %i!\n',itnum+1)
                end
                for ii=y_inp
                    plot_step(speicher{1,ii+1},f_size,infsplit,n,selim)
                end
                x_inp=input(subprompt5,'s');
            end
        end
        switch x_inp
            case 'p' % Change maximal steplength
                while 1
                    y_inp=input(subprompt2);
                    if ~isempty(y_inp) && y_inp>0
                        break
                    end
                    fprintf('Input must be a positive number!\n')
                end
                max_pred=y_inp;
            case 'n' % Change number of steps to carry out at once
                while 1
                    y_inp=input(subprompt1);
                    if ~isempty(y_inp) && y_inp>0 && rem(y_inp,floor(y_inp))==0
                        break
                    end
                    fprintf('Input must be a positive natural number!\n')
                end
                counter=y_inp;
            case 's' % Save & stop
                if ~only_exact
                    % Save
                    savestruct = struct;
                    
                    savestruct.delta_0     = delta_0;
                    savestruct.theta_0     = theta_0;
                    savestruct.itnum       = itnum;
                    savestruct.val         = val;
                    savestruct.exact_lp    = exact_lp;
                    savestruct.pred_tmp    = pred_tmp;
                    savestruct.corr_old    = corr_old;
                    savestruct.pit_stop    = pit_stop;
                    savestruct.halve_nb    = halve_nb;
                    if dispmdep
                        savestruct.Pts_Mat = Pts_Mat;
                        savestruct.y_mdep  = Y_mdep;
                        savestruct.z_mdep  = Z_mdep;
                    end
                end
                
                % savestruct is saved at the end of speicher, when the program is started from this point, it will access the last available data in this way
                speicher{1,end+1} = savestruct;
                if log
                    speicher{2,end} = logstruct;
                end
                if save_ws
                    speicher{3,end} = jump_cell;
                end
                
                lp_a0 = cell2mat(speicher(2,1:end-1));
                % Find the turning points in lp_a0
                tmp1 = find( lp_a0(1:end-1)>lp_a0(2:end) );
                tmp2 = find( lp_a0(1:end-1)<lp_a0(2:end) );
                if isempty(tmp1) || isempty(tmp2)
                    tur_pts=[];
                else
                    tur_pts=zeros(3,length(lp_a0));
                    kk=0;
                    for ii=2:length(lp_a0)
                        if (max(tmp1==ii) && max(tmp2==ii-1)) ||  (max(tmp2==ii) && max(tmp1==ii-1))
                            kk=kk+1;
                            tur_pts(1,kk)=lp_a0(ii-1);
                            tur_pts(2,kk)=lp_a0(ii);
                            tur_pts(3,kk)=lp_a0(ii+1);
                        end
                    end
                    tur_pts=tur_pts(:,1:kk);
                end
                
                % Reinstate the previous interpreter settings
                set(0, 'DefaultTextInterpreter', dtiii)
                set(0, 'DefaultLegendInterpreter', dliii)
                set(groot, 'DefaultAxesTickLabelInterpreter', datliii)
                
                return
        end
    end
    
    if skip==0 % If the step-length needs to be predicted
        % Prepare for the next step
        [a_0,a_p,tangent_new,steplength,predictor,delta_0,nda,nsd,theta_0,logstruct,eval_exact,exact_lp] = PredCorrStrat(problem,settings,predictor,a_c,a_p,x1,xfin,sol,tangent_new,ordnung,rho,psival,psi,steplength,exact_val,exact_lp,delta_0,theta_0,theta_max,log,logstruct,eval_exact,itnum,halve_nb,min_sl);
        
        % Accept the new mesh for the next step
        x1 = xfin;
    end
    
    itnum = itnum+1;
    jj = jj+1;
    
end

end

%% Local functions
function [tangente]=tangente_berechnen(DF)

help1=zeros(length(DF(:,1)),1);
help1(end)=1;

tangente=DF\help1;

tangente=tangente/sqrt(sum(tangente.*tangente));

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

function [] = figdataplot(ii,speicher,X_se,Y_se,Z_se,pathfollplot,path_pos,val,problem)

blue_rgb = [0 0.4470 0.7410];
red_rgb = [0.8500 0.3250 0.0980];
green_rgb = [0.4660 0.6740 0.1880];

figure(pathfollplot);
if isa(val,'cell')
    val1=cell2mat(val(1,:));
    val2=cell2mat(val(2,:));
else
    val2=val;
end
nn=size(val2,1);
n=length(Z_se);

f_size=8; % font size for plot

% data evolution plots
for jj=1:nn
    axes( 'Position', [ path_pos(1)+2*path_pos(3)/20, path_pos(2)+(1+(nn-jj)*10)*path_pos(4)/(nn*10), 7*path_pos(3)/20, 8*path_pos(4)/(nn*10)] ) ;
    hold on
    if ~isempty(val)
        for kk=1:ii-1
            if strcmp(problem,'testK_8')
                plot(1./[speicher{2,kk},val1(kk)],[speicher{jj+2,kk},val2(jj,kk)],'Color',green_rgb,'LineStyle','--') ;
                plot(1./[val1(kk),speicher{2,kk+1},],[val2(jj,kk),speicher{jj+2,kk+1}],'Color',red_rgb,'LineStyle',':','Marker','o','MarkerSize',6) ;
            else
                plot([speicher{2,kk},val1(kk)],[speicher{3,kk}(jj),val2(jj,kk)],'Color',green_rgb,'LineStyle','--') ;
                plot([val1(kk),speicher{2,kk+1},],[val2(jj,kk),speicher{3,kk+1}(jj)],'Color',red_rgb,'LineStyle',':','Marker','o','MarkerSize',6) ;
            end
        end
    end
    if strcmp(problem,'testK_8')
        a=speicher(2,1:ii);c=1./cell2mat(a);
    else
        a=speicher(2,1:ii);c=cell2mat(a);
    end
    b=speicher(3,1:ii);d=cell2mat(b);e=d(jj,:);
    plot(c,e,'Color',blue_rgb,'LineStyle','-','Marker','o','MarkerFaceColor',blue_rgb,'MarkerSize',3,'LineWidth',1.5)
    xlabel('pathfollowing parameter')
    if ii>50
        plot(c(50:50:ii),e(50:50:ii),'xw','MarkerSize',5)
    end
    ax = gca;
    ax.FontSize=f_size;
end

% solution evolution plot
if strcmp(problem,'testM_7')
    xstop=8;
    axes( 'Position', [ path_pos(1)+12*path_pos(3)/20, path_pos(2)+11*log_pos(4)/20, 7*path_pos(3)/20, 8*path_pos(4)/20] ) ;
    subplot(n*nn,2,2:2:2*3*nn)
    surf( [Y_se fliplr(1./Y_se(xstop:end-1))] ,X_se, [Z_se{1} fliplr(Z_se{3}(:,xstop:end-1))] )
    title('real part of the solution')
    xlabel('interval on mesh')
    ylabel('pathfollowing parameter')
    zlabel('function value')
    ax = gca;
    ax.FontSize=f_size;
    
    axes( 'Position', [ path_pos(1)+12*path_pos(3)/20, path_pos(2)+1*log_pos(4)/20, 7*path_pos(3)/20, 8*path_pos(4)/20] ) ;
    subplot(n*nn,2,2*3*nn+2:2:2*6*nn)
    surf( [Y_se fliplr(1./Y_se(xstop:end-1))] ,X_se, [Z_se{2} fliplr(Z_se{4}(:,xstop:end-1))] )
    title('real part of the solution')
    xlabel('interval on mesh')
    ylabel('pathfollowing parameter')
    zlabel('function value')
    ax = gca;
    ax.FontSize=f_size;
else
    for kk=1:n
        axes( 'Position', [ path_pos(1)+12*path_pos(3)/20, path_pos(2)+(1+(n-kk)*10)*path_pos(4)/(n*10), 7*path_pos(3)/20, 8*path_pos(4)/(n*10)] ) ;
        if strcmp(problem,'testK_8')
            surf(Y_se,1./X_se,Z_se{kk})
        else
            surf(Y_se,X_se,Z_se{kk})
        end
        name = strcat('solution function',num2str(kk));
        title(name)
        xlabel('interval on mesh')
        ylabel('pathfollowing parameter')
        zlabel('function value')
        ax = gca;
        ax.FontSize=f_size;
    end
end

end

function [] = figlogdataplot(ii,logstruct,pathfollplot,log_pos,tol)
nrow = 5;
ncol = 2;
f_size=8;
figure(pathfollplot);

axes( 'Position', [ log_pos(1)+2*log_pos(3)/20, log_pos(2)+42*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
plot(0:ii,logstruct.cpt)
hold on ;
plot(1:ii,logstruct.trm_halved(1:ii),'*r')
limsy=get(gca,'YLim');
set(gca,'Ylim',[0 limsy(2)+1]) ;
title('cputime')
sp(1) = gca;
sp(1).FontSize=f_size;
sp(1).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+12*log_pos(3)/20, log_pos(2)+42*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
plot(1:ii,logstruct.steplength)
hold on ;
plot(1:ii,logstruct.steplength_pred,'k')
plot(1:ii,zeros(1,ii),'k:')
plot(1:ii,logstruct.sl_adapt(1:ii),'*g')
set(gca,'Ylim',[-5/4*max(abs(logstruct.steplength)) 5/4*max(abs(logstruct.steplength))]) ;
legend('$s$','predicted $s$')
title('steplength $s$')
sp(2) = gca;
sp(2).FontSize=f_size;
sp(2).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+2*log_pos(3)/20, log_pos(2)+32*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
plot(1:ii,logstruct.theta_0)
hold on
plot(1:ii,ones(1,ii),'k--')
plot(1:ii,logstruct.theta_max,'r--')
%limsy=get(gca,'YLim');
set(gca,'Ylim',[0 3/2*max(abs(logstruct.theta_0))]);%max(1.1,limsy(2))]) ;
title('$\theta_0$')
sp(3) = gca;
sp(3).FontSize=f_size;
sp(3).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+12*log_pos(3)/20, log_pos(2)+32*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
semilogy(1:ii,logstruct.theta_max./logstruct.theta_0)
title('$\Theta_m / \Theta_0$')
sp(4) = gca;
sp(4).FontSize=f_size;
sp(4).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+2*log_pos(3)/20, log_pos(2)+22*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
plot(0:ii,logstruct.mesh_length(1,:))
hold on
plot(0:ii,logstruct.mesh_length(2,:),'k')
plot(1:ii,logstruct.mesha_halved(1:ii),'*r')
legend('at beginning','at end')
set(gca,'Ylim',[0 3/2*max(abs(logstruct.mesh_length(2,:)))]);
title('\# mesh points in step')
sp(5) = gca;
sp(5).FontSize=f_size;
sp(5).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+12*log_pos(3)/20, log_pos(2)+22*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
semilogy(1:ii,logstruct.norm_F(1,:))
hold on
semilogy(1:ii,logstruct.norm_F(2,:))
semilogy(1:ii,tol.*ones(1,ii),'k--')
legend('$||F(a_p)||$','$||F(a_c)||$')
title('$||F||$')
sp(6) = gca;
sp(6).FontSize=f_size;
sp(6).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+2*log_pos(3)/20, log_pos(2)+12*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
semilogy(1:ii,logstruct.orthogonal(1,:))
hold on
semilogy(1:ii,logstruct.orthogonal(2,:))
legend('$|F(a_p)(end)|$','$|F(a_c)(end)|$')
title('$|F(end)|$')
sp(7) = gca;
sp(7).FontSize=f_size;
sp(7).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+12*log_pos(3)/20, log_pos(2)+12*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
semilogy(1:ii,logstruct.norm_delta_0,'--')
hold on
semilogy(1:ii,logstruct.norm_delta_1,'--')
semilogy(1:ii,logstruct.corr_dist,'--')
semilogy(1:ii,logstruct.pdist_halved(1:ii),'*r')
semilogy(1:ii,logstruct.cdist_halved(1:ii),'*g')
legend('$||\Delta y_0||$','$||\Delta y_1||$','$||a_p-a_c||$')
title('$||...||$')
sp(8) = gca;
sp(8).FontSize=f_size;
sp(8).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+2*log_pos(3)/20, log_pos(2)+2*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
plot(1:ii,logstruct.c_s)
hold on
plot(1:ii,logstruct.cos_ab,'k')
plot(1:ii,logstruct.angle_halved(1:ii),'*r')
legend('$c_s$','cos($ab$)')
title('$c_s$ and cos($ab$)')
sp(9) = gca;
sp(9).FontSize=f_size;
sp(9).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+12*log_pos(3)/20, log_pos(2)+2*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
plot(0:ii,logstruct.max_error)
title('Maximal Error')
sp(10) = gca;
sp(10).FontSize=f_size;
sp(10).Title.FontWeight='normal';

for jj=1:nrow*ncol
    xlim(sp(jj),[0 ii+1])
end

end

function [a_0,a_p,tangent_new,steplength,predictor,delta_0,nda,nsd,theta_0,logstruct,eval_exact,exact_lp] = PredCorrStrat(problem,settings,predictor,a_c,a_p,x1,xfin,sol,tangent_new,ordnung,rho,psival,psi,steplength,exact_val,exact_lp,delta_0,theta_0,theta_max,log,logstruct,eval_exact,itnum,halve_nb,min_sl)

% In the case of mesh adaptation enabled, take special care to set the new starting point of the predictor step
if itnum~=-1 && (feval(settings,'meshAdaptation'))
    % If the mesh was adapted, then compute the tangent and a_0 for the new number of mesh points
    if length(x1)~=length(xfin) || max(abs(x1-xfin))>1e-12
        initP.initialMesh=xfin;
        initP.parameters = sol.parameters;
        
        initP.initialValues = coeffToValues(tangent_new, x1,ordnung,rho,xfin);
        initP = initial_coefficients(problem,xfin,initP,rho,0);
        tangent_new = [ initP.initialCoeff ; tangent_new(end) ] ;
        tangent_new(1:end-1) = tangent_new(1:end-1)/sqrt(sum(tangent_new.'*tangent_new)) ;
        
        initP.initialValues = coeffToValues(a_p, x1,ordnung,rho,xfin);
        initP = initial_coefficients(problem,xfin,initP,rho,0);
        a_p = [ initP.initialCoeff ; a_p(end) ] ;
        
        a_0 = [ sol.coeff ; sol.parameters ; a_c(end) ] ;
        
        predictor.x1 = xfin;
    else
        a_0 = a_c;
    end
else
    a_0 = a_c;
end

if itnum~=-1
    tangent_old = tangent_new;
end

% Compute the tangent at the new point
jac_F = functionFDF( 'DF', problem ,a_0,xfin,psival,psi,rho,[]);
tangent_new = tangente_berechnen(jac_F);

% Check for turning point
if itnum~=-1 && sign(tangent_old.'*tangent_new)<0
    fprintf('\n\n			 <strong>***</strong> Found <strong>TURNING POINT</strong> near %f! <strong>***</strong> \n\n',a_0(end));
    steplength=-steplength;
    % take the next value of the required values for next time
    if ~isempty(exact_val) && eval_exact==1
        if (exact_lp>1 && steplength<0) || (exact_lp<length(exact_val) && steplength>0)
            exact_lp = exact_lp + sign(steplength);
        else
            eval_exact=0;
        end
    elseif ~isempty(exact_val)
        eval_exact = 1;
    end
end

if itnum~=-1
    % Compute the new steplength
    c_s = abs(tangent_new.'*tangent_old) ;
    corr_dist = norm(a_p-a_0) ;
    tmp=2*norm(delta_0)/(c_s^2*corr_dist) ;
    beta = sqrt((theta_max/theta_0)*tmp);
    beta_max = max(min(1,predictor.maxSteplengthGrowth),predictor.maxSteplengthGrowth/(2^(0.5*halve_nb)));
    beta=min(beta,beta_max);
    steplength = beta*steplength;
    
    if predictor.display
        fprintf('\n  Steplength was %3.1e and is predicted to be %3.1e!\n',predictor.steplength,steplength);
    end
else
    c_s = 1;
end

if log && itnum~=-1
    logstruct.c_s = [ logstruct.c_s c_s ] ;
    
    logstruct.corr_dist = [logstruct.corr_dist corr_dist ] ;
    
    logstruct.steplength_pred = [logstruct.steplength_pred steplength];
end

% Update predictor
predictor.a_0 = a_0.';
predictor.tangent = tangent_new.';
predictor.steplength = steplength;
predictor.lambda_p_0 = a_0(end);

correct=0;
countercorr=0;
ncorrsteps=predictor.maxCorrSteps;
theta_0_min=Inf; steplength_min=steplength;
while correct~=1
    if abs(countercorr)==ncorrsteps
        steplength=steplength_min;
        predictor.steplength = steplength;
    end
    
    % Carry out predictor step
    a_p = a_0 + steplength.*tangent_new;
    
    % Update the value of the pathfollowing variable
    predictor.lambda_p_p = a_p(end);
    
    % Compute what will be needed to compute the new steplength after the corrector step
    F_a = functionFDF( 'F', problem ,a_p,xfin,psival,psi,rho,predictor);
    jac_F_a = functionFDF( 'DF', problem ,a_p,xfin,psival,psi,rho,predictor);
    delta_a = - jac_F_a\F_a ;
    new_a = a_p + delta_a ;
    F_new = functionFDF( 'F', problem ,new_a ,xfin,psival,psi,rho,predictor);
    simplified_delta = - jac_F_a\F_new ;
    
    nda = norm(delta_a);
    nsd = norm(simplified_delta);
    
    theta_0 = nsd / nda ;
    if theta_0<theta_0_min
        steplength_min=steplength;
    end
    
    delta_0=delta_a;
    beta = sqrt(1/2);%sqrt(theta_max/theta_0)
    predictor.steplength = steplength;
    
    % No Prediction-correction is needed if steplength is smaller than the minimally allowed step-length
    if abs(steplength)<min_sl
        break
    end
    
    % Prediction-correction is needed if theta_0 is too big
    if theta_0>theta_max*c_s^2 && abs(countercorr)~=ncorrsteps
        sl_tmp=steplength;
        steplength=beta*steplength;
        if predictor.display && countercorr==0
            fprintf('  <strong>Prediction correction:</strong>\n  Steplength is decreased from %1.1e to %1.1e',sl_tmp,steplength);
        elseif predictor.display
            fprintf(' to %1.1e',steplength)
        end
        countercorr=countercorr-1;
        % Procedure to increase steplength if needed -- disabled for now
        % elseif theta_0<theta_max*c0^2/4 && abs(countercorr)~=ncorrsteps
        %  sl_tmp=steplength;
        %  steplength=steplength/beta*0.99;
        %  countercorr=countercorr+1;
        %  fprintf('\n Steplength is increased from %f to %f!\n',sl_tmp,steplength);
    else
        if predictor.display && abs(countercorr)~=ncorrsteps
            fprintf('\n  Correction procedure should converge with steplength %3.1e!...\n',steplength);
        elseif predictor.display
            fprintf('\n  The steplength %f is used!...\n',steplength);
        end
        correct=1;
    end
end

end

function [speicher,speicher_exact,val,corr_old] = correct_data(speicher,speicher_exact,val,data,ordnung,rho)
sol = speicher{1,1};
a_0 = [ sol.coeff ; sol.parameters ; sol.lambda_p ] ;
x1 = sol.x1;
speicher{3,1} = data(x1,a_0,ordnung,rho);
for ii=2:size(speicher,2)
    predictor = speicher{1,ii}.predictor;
    a_p = a_0 + predictor.steplength.*predictor.tangent';
    if iscell(val)
        val{2,ii-1} = data(x1,a_p,ordnung,rho);
    end
    sol = speicher{1,ii};
    a_0 = [ sol.coeff ; sol.parameters ; sol.lambda_p ] ;
    x1 = sol.x1;
    speicher{3,ii} = data(x1,a_0,ordnung,rho);
end

data_tmp = data(speicher{1,end-1}.x1,a_p,ordnung,rho);
corr_old = sqrt((a_0(end)-a_p(end))^2+(max(abs(speicher{3,end}-data_tmp)))^2);

for ii=1:size(speicher_exact,2)
    sol = speicher_exact{1,ii};
    a_0 = [ sol.coeff ; sol.parameters ; sol.lambda_p ] ;
    x1 = sol.x1;
    speicher_exact{3,ii} = data(x1,a_0,ordnung,rho);
end

end

function [] = plot_step(sol,f_size,infsplit,n,selim)
exact_plot=figure('units','normalized','outerposition',[0.3 0.2 0.4 0.6],'PaperUnits','normalized');
ax = gca;
ax.FontSize=f_size;
clf ;
exact_plot.Color = 'White' ;
if infsplit
    selimpt = find((1./sol.x1)>=selim,1,'last');
end
for ii = 1:n
    axes( 'Position', [ 1/10, (1+(n-ii)*10)/(n*10), 8/10, 8/(n*10)] ) ;
    if infsplit
        plot([sol.x1 fliplr(1./sol.x1(selimpt:end-1))],[sol.valx1(ii,:) fliplr(sol.valx1(n+ii,selimpt:end-1))],'LineWidth',1.5)
    else
        plot(sol.x1,sol.valx1(ii,:),'LineWidth',1.5)
    end
    plottitle=['solution ',num2str(ii),' for $\lambda_p=$',num2str(sol.lambda_p)];
    title(plottitle)
    ax = gca;
    ax.FontSize=f_size;
    ax.Title.FontWeight='normal';
end
end

