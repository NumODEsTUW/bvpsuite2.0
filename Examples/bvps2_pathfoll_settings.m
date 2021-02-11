function [ret] = bvps2_pathfoll_settings(request)
% Solver settings for pathfollowing problem in Manual

switch request
    case 'mesh'
        ret=0:1/100:1;
    case 'collMethod'
        ret='gauss';
    case 'collPoints'
        ret=3;
    case 'meshAdaptation'
        ret=1;
    case 'errorEstimate'
        ret=1;
    case 'absTolSolver'
        ret=1e-6;
    case 'relTolSolver'
        ret=1e-6;
    case 'absTolMeshAdaptation'
        ret=1e-4;
    case 'relTolMeshAdaptation'
        ret=1e-4;
    case 'minInitialMesh'
        ret=50;
        
    case 'finemesh'
        ret = 0;
    
    case 'allowTRM'
        ret=1;
    case 'maxFunEvalsTRM'
        ret=90000000;
    case 'maxIterationsTRM'
        ret=90000000;
    case 'lambdaMin'
        ret = 0.001;
    case 'maxAdaptations'
        ret=18;
    case 'switchToFFNFactor'
        ret=0.5;
    case 'updateJacFactor'
        ret=0.5;
    case 'K'
        ret=200;
    
    case 'thetaMax' % Newly added, controls first Newton contraction factor 
        ret=0.01;
    case 'maxCorrSteps' % Maximal # of times the predicted step-length gets corrected if the first Newton contraction factor is not smaller than theta_max
        ret=5;
    case 'maxSteplengthGrowth' % Maximal growth of the steplength when predicting it for the next step
        ret=2;
    case 'angleMin' % Minimal cosine of angle between succesive tangents
        ret=0.75;
    case 'meshFactorMax' % Controls by which factor the number of points in the mesh can be augmented at once
        ret=2;
    case 'PredLengthFactor' % Maximal 1/factor difference between predictor and corrector length
        ret=1;
    case 'CorrLengthGrowth' % Maximal allowed growth of corrector step compared to corrector step of the previous step
        ret=8;
        
end

end

%{
- case thetaMax, angleMin, meshFactorMax, maxCorrSteps, PredLengthFactor, 
CorrLengthGrowth, maxSteplengthGrowth wurden hinzugefügt
%}