function [ret] = bvps2_semiinf_settings(request)
% Solver settings for problem posed on a semi-inifinite interval in Manual

switch request
    case 'mesh'
        ret=0:1/50:1;
    case 'collMethod'
        ret='gauss';
    case 'collPoints'
        ret=5;
    case 'meshAdaptation'
        ret=1;
    case 'errorEstimate'
        ret=1;
    case 'absTolSolver'
        ret=1e-10;
    case 'relTolSolver'
        ret=1e-10;
    case 'absTolMeshAdaptation'
        ret=1e-9;
    case 'relTolMeshAdaptation'
        ret=1e-9;
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
    
    case 'thetaMax'
        ret=0.01;
    case 'maxCorrSteps'
        ret=5;
    case 'maxSteplengthGrowth'
        ret=5;
    case 'angleMin'
        ret=0.75;
    case 'meshFactorMax'
        ret=1.5;
    case 'PredLengthFactor'
        ret=1;
    case 'CorrLengthGrowth'
        ret=8;
        
end

end