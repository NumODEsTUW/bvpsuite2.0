function [ret] = bvps2_bcint_settings(request)
switch request
    case 'mesh'
        ret=0:1/100:1;
    case 'collMethod'
        ret='gauss';
    case 'collPoints'
        ret=2;
    case 'meshAdaptation'
        ret=1;
    case 'errorEstimate'
        ret=0;
    case 'absTolSolver'
        ret=10^-6;
    case 'relTolSolver'
        ret=10^-6;
    case 'absTolMeshAdaptation'
        ret=10^-4;
    case 'relTolMeshAdaptation'
        ret=10^-4;
    case 'minInitialMesh'
        ret=50;
        
    case 'finemesh'
        ret = 0;
        
    case 'allowTRM'
        ret=1;
    case 'maxFunEvalsTRM'
        ret=100;
    case 'maxIterationsTRM'
        ret=100;
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