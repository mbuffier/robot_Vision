function error_location = bundle_adjustment_test(pourcentageWanted, ITERATIONS, nbFeatures)

import gtsam.*

% Options
NUM_FRAMES = 60; % 0 for all
ADD_NOISE = 1;
blenddir = strcat(fileparts(mfilename('fullpath')), '/../blender/');

% Load data 
camera_gt = dlmread(strcat(blenddir, 'camera_poses.txt'));
features_gt = dlmread(strcat(blenddir, 'tracks_dist.txt'));
landmarks_gt = dlmread(strcat(blenddir, 'landmarks_3d.txt'));
landmarks_used = zeros(size(landmarks_gt,1),1);

% cut the number of measurement to the good number of frame 
features_gt = features_gt(1:NUM_FRAMES,:) ;

% camera position noisy 
camera_gt_noisy = zeros(size(camera_gt)) ;
landmarks_gt_noisy = zeros(size(landmarks_gt)) ;

if NUM_FRAMES < 1
    NUM_FRAMES = size(camera_gt, 1);
end

calib = Cal3_S2( ...
    634.8, ... % focal
    634.8, ... % focal
    0, ... % skew
    480,... % center
    270); % center

%% Setup noise
measurementNoiseSigma = 3;
pointNoiseSigma = 0.1;
rotationSigma = 0.2;
positionSigma = 3;
poseNoiseSigmas = [ positionSigma positionSigma positionSigma ...
                    rotationSigma rotationSigma rotationSigma]';
posePriorNoise  = noiseModel.Diagonal.Sigmas(poseNoiseSigmas);
pointPriorNoise  = noiseModel.Isotropic.Sigma(3,pointNoiseSigma);
measurementNoise = noiseModel.Isotropic.Sigma(2,measurementNoiseSigma);

%% Add noise to input data
if ADD_NOISE == 1
    disp('Adding noise...')
    
    for i=1:NUM_FRAMES
        
        % 2D features
        f = 1;
        while f < size(features_gt, 2) && features_gt(i,f) > 0
            features_gt(i,f+1:f+2) = features_gt(i,f+1:f+2) + measurementNoiseSigma * randn(1,2);
            f = f + 4;
        end
        
        % Camera Poses
        rot = Rot3.Quaternion( camera_gt(i,8), ...
                               camera_gt(i,5), ...
                               camera_gt(i,6), ...
                               camera_gt(i,7)).xyz;
        
        rot = rot + rotationSigma * randn(3,1);
        % I had a line here otherwise it gave an error
        rotq_inter = Rot3.RzRyRx(rot) ;
        rotq = rotq_inter.quaternion();
        camera_gt_noisy(i,8) = rotq(1) ;
        camera_gt_noisy(i,5) = rotq(2) ;
        camera_gt_noisy(i,6) = rotq(3) ;
        camera_gt_noisy(i,7) = rotq(4) ;         
        camera_gt_noisy(i,2:4) = camera_gt(i,2:4) + positionSigma * randn(1,3);       
    end
    
    % 3D landmarks
    landmarks_gt_noisy = landmarks_gt_noisy + randn(size(landmarks_gt_noisy));  
end

graph = NonlinearFactorGraph;
initialEstimate = Values;

%% Determine the number of measurement in pourcentage for this landmark 
% store the new number of measurements
newNumbMeas = zeros(nbFeatures,1) ;
% store the id of the landmark 
thisLandmarkNumber = zeros(nbFeatures,1) ;
% how many time it has been seen 
currentNumber = zeros(nbFeatures,1) ;
index = 1 ;

% go throught the number of landmark asked 
for i=1:4:nbFeatures*4
    thisLandmarkNumber(index,1) = features_gt(2,i) ;
    [whichFrames,~] = find(features_gt == thisLandmarkNumber(index,1)) ;
    
    % determine how many measurement there are in the system
    oldNumbMeas = size(whichFrames,1) ;
    
    % number of measurement wanted
    newNumbMeas(index,1) = floor(pourcentageWanted/100 * oldNumbMeas) ;
    currentNumber(index,1) = 0 ;
    index = index+1 ;
end

%% Add factors for all measurements
for i=1:NUM_FRAMES
    
    fprintf('Adding frame %d to graph...\n', i)
    
    cam_pose = Pose3(Rot3.Quaternion(camera_gt_noisy(i,8), camera_gt_noisy(i,5), camera_gt_noisy(i,6), camera_gt_noisy(i,7)), ...
        Point3(camera_gt_noisy(i,2), camera_gt_noisy(i,3), camera_gt_noisy(i,4)));
    
    % ****************************** TODO ****************************
    % Initialization of camera poses with positional priors
    graph.add(PriorFactorPose3(symbol('p', i), cam_pose, posePriorNoise));
    initialEstimate.insert(symbol('p', i), cam_pose) ;
    % ************************************************************
    
    f = 1; % column of current feature ID
    while f < size(features_gt, 2) && features_gt(i,f) > 0
        
        feature_id = features_gt(i,f);
        position = find(feature_id == thisLandmarkNumber) ;
               
        % if this landmark isn't concerned : add this measurement 
        if isempty(position)
            % find the position of this feature (noisy)
            feature_pos = Point3(landmarks_gt_noisy(feature_id,1),landmarks_gt_noisy(feature_id,2),landmarks_gt_noisy(feature_id,3));
            
            % ****************************** TODO ****************************
            % add the measurement between this feature and the camera
            thisMeasurement = Point2(features_gt(i,f+1), features_gt(i,f+2)) ;
            graph.add(GenericProjectionFactorCal3_S2(thisMeasurement, measurementNoise, symbol('p', i),  symbol('f', feature_id), calib));
            % ************************************************************
        
        % else we verify before that the number of measurement added is
        % below the threshold 
        elseif currentNumber(position,1) < newNumbMeas(position,1)
            % find the position of this feature (noisy)
            feature_pos = Point3(landmarks_gt_noisy(feature_id,1),landmarks_gt_noisy(feature_id,2),landmarks_gt_noisy(feature_id,3));
            % add the measurement
            thisMeasurement = Point2(features_gt(i,f+1), features_gt(i,f+2)) ;
            graph.add(GenericProjectionFactorCal3_S2(thisMeasurement, measurementNoise, symbol('p', i),  symbol('f', feature_id), calib));
            % ************************************************************
            currentNumber(position,1) = currentNumber(position,1) + 1;
        end
        
        % Initialise the point near ground-truth if we haven't seen it 
        if landmarks_used(feature_id,1) < 1
            
            % ****************************** TODO ****************************
            % Prior and initial estimate for the features point in the world
            graph.add(PriorFactorPoint3(symbol('f', feature_id), feature_pos, pointPriorNoise));
            initialEstimate.insert(symbol('f', feature_id), feature_pos) ;
            % ************************************************************
            
            landmarks_used(feature_id,1) = 1;
        end
        f = f + 4;
    end
end

optimizer = LevenbergMarquardtOptimizer(graph, initialEstimate);

% To evaluate the algorithm : we store the mean distance between the
% estimated and the groundtruth 
thisdist = 0 ;
for i=1:ITERATIONS
    fprintf('Starting iteration %d...\n', i);
    optimizer.iterate();
    
    result = optimizer.values();
    
    if i==ITERATIONS
        for j=1:size(thisLandmarkNumber,1)
            % to explore the distance to the ground-truth 
            feature_id = thisLandmarkNumber(j,1) ;
            thisNewpose = result.atPoint3(symbol('f', feature_id));
            pose_found = [thisNewpose.x, thisNewpose.y, thisNewpose.z] ;
            thisGT = landmarks_gt(thisLandmarkNumber(j,1),:) ;
            thisdist = thisdist+norm(pose_found - thisGT) ;
        end
    end
end
error_location = thisdist/size(thisLandmarkNumber,1) ;
end
