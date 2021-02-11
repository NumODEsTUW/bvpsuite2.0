function [ret] = default_settings(request)

switch request
    case 'mesh'
        ret=0:1/50:1;
    case 'collMethod'
        ret='gauss';
    case 'collPoints'
        ret=3;
    case 'meshAdaptation'
        ret=1;
    case 'errorEstimate'
        ret=0;
    case 'absTolSolver'
        ret=1e-12;
    case 'relTolSolver'
        ret=1e-12;
    case 'absTolMeshAdaptation'
        ret=1e-9;
    case 'relTolMeshAdaptation'
        ret=1e-9;
    case 'minInitialMesh'
        ret=50;
        
    case 'finemesh'
        ret = 0;
    
    % Solver settings -- see https://repositum.tuwien.at/handle/20.500.12708/9088
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
    
    % Pathfollowing settings -- see https://repositum.tuwien.at/handle/20.500.12708/1167
    case 'thetaMax'
        ret=0.05;
    case 'maxCorrSteps'
        ret=5;
    case 'maxSteplengthGrowth'
        ret=4;
    case 'angleMin'
        ret=0.75;
    case 'meshFactorMax'
        ret=2;
    case 'PredLengthFactor'
        ret=1;
    case 'CorrLengthGrowth'
        ret=4;
        
end

end

