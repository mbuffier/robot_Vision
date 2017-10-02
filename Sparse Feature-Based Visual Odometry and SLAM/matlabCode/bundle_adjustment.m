function [groundTruth_position,initial_position,output_position, spacial_uncert, error_location, ATE_result, accuracyArray] = bundle_adjustment(uncertainty, ATE, accuracy, ITERATIONS)

import gtsam.*

% Options
NUM_FRAMES = 0; % 0 for all
ADD_NOISE = 1;
blenddir = strcat(fileparts(mfilename('fullpath')), '/../blender/');
outdir = strcat(fileparts(mfilename('fullpath')), '/output/');

% Load data 
camera_gt = dlmread(strcat(blenddir, 'camera_poses.txt'));
features_gt = dlmread(strcat(blenddir, 'tracks_dist.txt'));
landmarks_gt = dlmread(strcat(blenddir, 'landmarks_3d.txt'));
landmarks_used = zeros(size(landmarks_gt,1),1);

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

%% Add factors for all measurements
for i=1:NUM_FRAMES
    
    fprintf('Adding frame %d to graph...\n', i)
    
    cam_pose = Pose3(Rot3.Quaternion(camera_gt_noisy(i,8), camera_gt_noisy(i,5), camera_gt_noisy(i,6), camera_gt_noisy(i,7)), ...
        Point3(camera_gt_noisy(i,2), camera_gt_noisy(i,3), camera_gt_noisy(i,4)));
    
    % ****************************** TODO ****************************
    % Initialization of camera poses with positional priors
    % graph.add(PriorFactorPose3(symbol('p', i), cam_pose, posePriorNoise));
    initialEstimate.insert(symbol('p', i), cam_pose) ;
    % ************************************************************
    
    f = 1; % column of current feature ID
    while f < size(features_gt, 2) && features_gt(i,f) > 0
        
        feature_id = features_gt(i,f);
        
        % find the position of this feature (noisy) 
        feature_pos = Point3(landmarks_gt_noisy(feature_id,1),landmarks_gt_noisy(feature_id,2),landmarks_gt_noisy(feature_id,3));
        
        % Initialise the point near ground-truth
        if landmarks_used(feature_id,1) < 1
            
            % ****************************** TODO ****************************
            % Prior and initial estimate for the features point in the world
            graph.add(PriorFactorPoint3(symbol('f', feature_id), feature_pos, pointPriorNoise));
            initialEstimate.insert(symbol('f', feature_id), feature_pos) ;
            % ************************************************************
            
            landmarks_used(feature_id,1) = 1;
        end
        
        % ****************************** TODO ****************************
        % add the measurement between this feature and the camera
        thisMeasurement = Point2(features_gt(i,f+1), features_gt(i,f+2)) ;
        graph.add(GenericProjectionFactorCal3_S2(thisMeasurement, measurementNoise, symbol('p', i),  symbol('f', feature_id), calib));
        % ************************************************************
        
        f = f + 4;
    end
end

optimizer = LevenbergMarquardtOptimizer(graph, initialEstimate);

% To evaluate the algorithm 
spacial_uncert = zeros(ITERATIONS, 1) ;
error_location = zeros(ITERATIONS, 1) ;
camera_out = zeros(size(camera_gt));
ATE_result = zeros(ITERATIONS, 1) ;
accuracyArray = zeros(ITERATIONS, 1) ;

for i=1:ITERATIONS
    fprintf('Starting iteration %d...\n', i);
    optimizer.iterate();
    
    % provided with the information of covariances and distances to ground-truth 
    % for the landmarks to see the evolution 
    if uncertainty == 1
        result = optimizer.values();
        
        % **** I suppressed this call to marginal due to errors ****
        %marginals = Marginals(graph, result);        
        
        M = size(landmarks_used,1) ;
        %thisNorm = 0 ;
        distMean = 0 ;
        
        % for each landmarks 
        for k = 1:M
            if landmarks_used(k)==1
                % spacial uncertainty : the mean of the
                % magnitude of the standard deviation vector
%                 thisCov = marginals.marginalCovariance(symbol('f', k)) ;
%                 thisNorm = thisNorm+norm(diag(thisCov)) ;
                
                % error in the location : mean of the distance between
                % grund truth and the estimation 
                pose = result.atPoint3(symbol('f', k));
                pose_found = [pose.x, pose.y, pose.z] ;
                thisGT = landmarks_gt(k,:) ;
                
                thisdist = norm(pose_found - thisGT);
                distMean = distMean + thisdist ;
            end
        end
        % normalization 
%       spacial_uncert(i,1) = thisNorm./sum(landmarks_used,1) ;
        error_location(i,1) = distMean./sum(landmarks_used,1) ;
    end
    
    % use the ATE metric to observe the accuraty of the camera positions
    if ATE==1
        % create a matrix with the camera positions information  
        for j=1:size(camera_gt,1) 
            pose = result.atPose3(symbol('p', j));
            pos = pose.translation();
            quat = pose.rotation().quaternion();
            camera_out(j,:) = [camera_gt(j,1) pos.x pos.y pos.z quat(2) quat(3) quat(4) quat(1)];
        end
        % print it in a file 
        dlmwrite(strcat(outdir, 'vo_output_all_BA.txt'), camera_out, 'delimiter','\t','precision',6);
        
        % compute the ATE measurement for this iteration 
        [~,thisATE] = system('python ../external/tum-evaluation/evaluate_ate.py ../blender/camera_poses.txt output/vo_output_all_BA.txt --plot allCam_ate.pdf') ;
        ATE_result(i,1) = str2double(thisATE) ;
    end
    
    % store the error for this iteration 
    if accuracy==1
        result = optimizer.values();
        thisError = graph.error(result) ;
        accuracyArray(i,1) = thisError ;
    end
end

disp('Retrieving results...');
result = optimizer.values();
fprintf('Initial error: %f\n', graph.error(initialEstimate));
finalError = graph.error(result) ;
fprintf('Final error: %f\n',finalError);

% between 0 and 1 to see the evolution 
accuracyArray = (accuracyArray-min(accuracyArray(:)))./(max(accuracyArray(:))-min(accuracyArray(:))) ;

%% Output noisy input data and optimised result
disp('Exporting results...');
dlmwrite('input_poses.txt',camera_gt_noisy(1:NUM_FRAMES,:),'delimiter','\t','precision',6);
output_poses = zeros(NUM_FRAMES,size(camera_gt,2));

for p = 1:NUM_FRAMES
    pose = result.atPose3(symbol('p', p));
    pos = pose.translation();
    quat = pose.rotation().quaternion();
    output_poses(p,:) = [camera_gt_noisy(p,1) pos.x pos.y pos.z quat(4) quat(2) quat(3) quat(1)];
end

dlmwrite('output_poses.txt',output_poses,'delimiter','\t','precision',6);

%% ************* Todo : Generate a plot of the result *********************
groundTruth_position = camera_gt(:,2:4) ;
initial_position = camera_gt_noisy(:,2:4) ;
output_position = output_poses(:,2:4) ;
%*******************************************************************


end
