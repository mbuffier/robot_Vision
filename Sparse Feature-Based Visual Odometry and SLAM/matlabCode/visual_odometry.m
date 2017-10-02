function [rte_all, ate_all, rte_window, ate_window, stand_dev, meanMagSD, timeArray, covAfterBA] = visual_odometry(NUM_INITIALISE, WINDOW_SIZE, GPS, covariance, RTE_ATE, time)

import gtsam.*

%% Data Options
NUM_FRAMES = 0; % 0 for allus

% *************  TODO : find a nice value ************
% put into arguments 
% NUM_INITIALISE = 10 ; 
% WINDOW_SIZE = 10 ; 
% ************************

MAX_ITERATIONS = 30;
LOOP_THRESHOLD = 1000;

ADD_NOISE = true; % Try false for debugging (only)
SAVE_FIGURES = true; % Save plots to hard-drive?
PLOT_LANDMARKS = true;
PLOT_FEATURES = true;

outdir = strcat(fileparts(mfilename('fullpath')), '/output/');
blenddir = strcat(fileparts(mfilename('fullpath')), '/../blender/');
mkdir(outdir);

%% Load data
global features_gt;
global camera_gt;
global camera_out;
global landmarks_gt;
global landmarks_out;
global landmarks_used;
global calib;
camera_gt = dlmread(strcat(blenddir, 'camera_poses.txt')); % each line is: frame_id,x,y,z,qx,qy,qz,qw
camera_out = zeros(size(camera_gt)); % each line is: frame_id,x,y,z,qx,qy,qz,qw
features_gt = dlmread(strcat(blenddir, 'tracks_dist.txt')); % each line is list of: landmark_id,feature_x,feature_y,feature_d,...
landmarks_gt = dlmread(strcat(blenddir, 'landmarks_3d.txt')); % each line is: x,y,z
landmarks_out = zeros(size(landmarks_gt)); % each line is: x,y,z,last_seen_frame_id

% This table stores the last usage of each landmark (frame_id, whenever the
% landmark was visible somewhere in the window).
% The structure is used for two purposes:
% 1st: Use landmarks only once, when building a window
% 2nd: Decide whether a loop occured, based on stored frame-IDs
landmarks_used = zeros(size(landmarks_gt,1),1);

if NUM_FRAMES < 1
    NUM_FRAMES = size(camera_gt, 1);
end
calib = Cal3_S2( ...
    634.8, ... % focal
    634.8, ... % focal
    0, ... % skew
    480,... % center
    270); % center

% Setup noise
measurementNoiseSigma = 2;
pointNoiseSigma = 3;
GPSNoiseSigma = 2;

pointPriorNoise  = noiseModel.Isotropic.Sigma(3,pointNoiseSigma);
measurementNoise = noiseModel.Isotropic.Sigma(2,measurementNoiseSigma);
poseGPSNoise  = noiseModel.Isotropic.Sigma(3,GPSNoiseSigma);

%% Sanity checks
assert(WINDOW_SIZE <= NUM_INITIALISE);
assert(NUM_FRAMES <= size(camera_gt,1));


%% *************  TODO : initialization to evaluate ************
% RTE and ATE computed from the start of the run to the current time
ate_all = zeros(size(camera_gt,1), 1) ;
rte_all = zeros(size(camera_gt,1), 1) ;

% RTE and ATE computed over just the length of the sliding window
ate_window = zeros(size(camera_gt,1), 1) ;
rte_window = zeros(size(camera_gt,1), 1) ;

% for the standart deviation observation 
stand_dev = zeros(size(camera_gt,1),3) ;
meanMagSD = zeros(size(camera_gt,1),1) ;
covAfterBA = zeros(size(camera_gt,1),3) ;

% for the time of the computation 
timeArray = zeros(size(camera_gt,1),1) ;
firstBA = 0 ; 
% ************************************************

%% Setup figure
figure;
% Optionally choose size according to your needs
%set(gcf,'units','pixel');
%set(gcf,'papersize',[960,960]);
%set(gcf,'paperpositionmode','auto');
subplot(2,1,1); % Top for map, bottom for video
axis equal;
hold on;
% Plot ground-truth trajectory as green line
plot3(camera_gt(:,2), camera_gt(:,3), camera_gt(:,4), 'g');

%% Initialise with ground-truth data
% Of course, this can not be done in real-life(!)
for i=1:NUM_INITIALISE
    fprintf('Loading ground-truth data for frame %d...\n', i);
    camera_out(i,:) = camera_gt(i,:);
    f = 1; % column of current feature ID
    while f < size(features_gt, 2) && features_gt(i,f) > 0
        feature_id = features_gt(i,f);
        % Initialise the point near ground-truth
        if ~in_map(feature_id)
            landmarks_out(feature_id,1:3) = landmarks_gt(feature_id,1:3);
            landmarks_used(feature_id,1) = NUM_INITIALISE;
        end
        f = f + 4;
    end
end

%% Add noise to 2D observations
if ADD_NOISE
    disp('Adding noise...')
    for i=1:NUM_FRAMES
        % 2D features
        f = 1;
        while f < size(features_gt, 2) && features_gt(i,f) > 0
            features_gt(i,f+1:f+2) = features_gt(i,f+1:f+2) + measurementNoiseSigma * randn(1,2);
            f = f + 4;
        end
    end
end

%% Visual-Odometry part: Optimise a window for each new frame
for i=NUM_INITIALISE+1:NUM_FRAMES

    tic % measure time

    %% Load image (for plotting only)
    subplot(2,1,2);
    hold off;
    img = imread(strcat(blenddir, '/frames/',sprintf('%04d',i),'.jpg'));
    imshow(img);
    subplot(2,1,1);

    %% Adding new frame
    window_start = i - WINDOW_SIZE;
    fprintf('Building graph for frame %d (Starting at %d)...\n', i, window_start);
    graph = NonlinearFactorGraph;
    initialEstimate = Values;

    landmarks_in_window = 0;
    redetected_landmarks = 0;
    for j=window_start:i

        if j==i
            % Initialise assuming constant motion model
            lPose = get_pose(camera_out,j-1);
            llPose = get_pose(camera_out,j-2);
            rPose = llPose.between(lPose);
            cam_pose = lPose.compose(rPose);

            % Initialise from last pose (method above is better)
            % cam_pose = get_pose(camera_out,j-1);
        else
            cam_pose = get_pose(camera_out,j);
        end

        % GTSAM does not ensure that the rotation component of a pose is
        % orthogonal, hence do this manually (important!)
        rdet = det(cam_pose.rotation.matrix);
        if abs(rdet-1) > 0.0001
            fprintf('Correcting R-det: %f \n', rdet);
            [U,~,V]=svd(cam_pose.rotation.matrix);

            % ****************************** TODO ****************************
            %correct pose (single line missing)
            cam_pose = Pose3(Rot3(U*V.'), cam_pose.translation) ;
            % ****************************************************************
        end
        
        % ****************************** TODO ****************************
        % Initialization of camera poses without positional priors
        initialEstimate.insert(symbol('p', j), cam_pose) ;
        % ****************************************************************
        
        % ************************* TODO : add GPS observations ***********
        if GPS==1
            if j==1 || j==250 || j==500
                % add the standard deviation of 2 (from -1 to 1) to the
                % position 
                newTrans = camera_gt(j,2:4) + randn(1,3) ;
                % construct the new vector
                GPS_prior = [camera_gt(j,1) newTrans camera_gt(j,5:8)];
                % get the pose from it 
                cam_pose_GPS = get_pose(GPS_prior,1);
                
                % add the prior to the graph 
                graph.add(PoseTranslationPrior3D(symbol('p', j), cam_pose_GPS, poseGPSNoise));
            end
        end
        % **********************************************************
        
        % The following loop iterates over the (2D) features within frame 'j'
        f = 1; % column of current feature ID
        
        while f < size(features_gt, 2) && features_gt(j,f) > 0

            feature_id = features_gt(j,f);

            % Initialise the point near ground-truth
            if landmarks_used(feature_id) ~= i

                if in_map(feature_id)
                    % This landmark is known to the system
                    % Was it seen in the current window or do we have a
                    % potential loop-closure?
                    if landmarks_used(feature_id) < window_start 
                        redetected_landmarks = redetected_landmarks+1 ;
                        % Skip redetected landmarks and potentially close
                        % loop later
                        f = f + 4;
                        continue;
                    else
                        % The landmark was already seen, when building the window of frames
                        feature_pos = get_landmark(landmarks_out, feature_id);
                    end
                else
                    % "Triangulate"
                    rep = [(features_gt(j,f+1) - calib.px) / calib.fx, (features_gt(j,f+2) - calib.py) / calib.fy, 1] * features_gt(j,f+3);
                    rep = cam_pose.matrix * [rep'; 1];
                    feature_pos = Point3(rep(1),rep(2),rep(3));

                    % Ground-Truth position: Don't do this, unless you are debugging!
                    % feature_pos = Point3(landmarks_gt(feature_id,1),landmarks_gt(feature_id,2),landmarks_gt(feature_id,3));
                end

                % ****************************** TODO ****************************
                % Prior and initial estimate for the features point in the
                % world
                graph.add(PriorFactorPoint3(symbol('f', feature_id), feature_pos, pointPriorNoise));
                initialEstimate.insert(symbol('f', feature_id), feature_pos) ;
                % ****************************************************************

                landmarks_used(feature_id) = i;
                landmarks_in_window = landmarks_in_window + 1;
            end

            % ****************************** TODO ****************************
            % Insertion of measurements
            thisMeasurement = Point2(features_gt(j,f+1), features_gt(j,f+2)) ;
            graph.add(GenericProjectionFactorCal3_S2(thisMeasurement, measurementNoise, symbol('p', j),  symbol('f', feature_id), calib));
            % ****************************************************************

            f = f + 4;
        end
        subplot(2,1,1);
    end
    
    %% ****************************** TODO : Loop closing ****************************
    if redetected_landmarks > LOOP_THRESHOLD
        % Build a more complex graph if a loop was detected
        disp('Loop detected! Building large graph...');
        
        firstBA=firstBA+1 ;

        % add all the landmark and camera position found until then from
        % the first frame to window-1
        for j=1:window_start-1
            
            % the camera position 
            cam_pose = get_pose(camera_out,j);

            % corrected
            rdet = det(cam_pose.rotation.matrix);
            if abs(rdet-1) > 0.0001
                fprintf('Correcting R-det: %f \n', rdet);
                [U,~,V]=svd(cam_pose.rotation.matrix);
                cam_pose = Pose3(Rot3(U*V.'), cam_pose.translation) ;
            end
            
            % add to the initial estimate 
            initialEstimate.insert(symbol('p', j), cam_pose) ;
            
            % The following loop iterates over the (2D) features within frame 'j'
            f = 1; % column of current feature ID
            
            % add all the features for this camera 
            while f < size(features_gt, 2) && features_gt(j,f) > 0
                % find the id for this features 
                feature_id = features_gt(j,f);
                
                % if the landmark has been seen and not in the previous
                % window, we add it to the initial values and graph 
                if landmarks_used(feature_id) ~= i && in_map(feature_id)
                    feature_pos = get_landmark(landmarks_out, feature_id);
                    
                    graph.add(PriorFactorPoint3(symbol('f', feature_id), feature_pos, pointPriorNoise));
                    initialEstimate.insert(symbol('f', feature_id), feature_pos) ;
                    
                    landmarks_used(feature_id) = i;
                    landmarks_in_window = landmarks_in_window + 1;
                end
                
                % Insertion of measurements
                thisMeasurement = Point2(features_gt(j,f+1), features_gt(j,f+2)) ;
                graph.add(GenericProjectionFactorCal3_S2(thisMeasurement, measurementNoise, symbol('p', j),  symbol('f', feature_id), calib));
                
                f = f + 4;
            end
            
            % The following code has to be informed that the window is bigger now
            window_start = 1;
        end
        disp('Graph ready to be optimized..');
    end
    % ****************************************************************

    optimizer = LevenbergMarquardtOptimizer(graph, initialEstimate);
    initialError = graph.error(initialEstimate);
    
    %% Plot ground-truth and initial estimate of markers (2D)
    subplot(2,1,2);
    hold on;
    plot_features(i, 'g+'); % ground-truth
    %plot_features(i, 'cx', cam_pose, landmarks_gt); % initial estimate
    subplot(2,1,1);
    %% Plot ground-truth and initial estimate in 3D
    % plot3(camera_gt(i,2),camera_gt(i,3),camera_gt(i,4), 'g+');
    % plot3(cam_pose.translation().x,cam_pose.translation().y,cam_pose.translation().z,'cx');

    %optimizer.optimizeSafely;
    last_error = initialError;
    for m=1:MAX_ITERATIONS
        optimizer.iterate();

        %% Print iteration information
        result = optimizer.values();
        error = graph.error(result);
        fprintf('Initial error: %f, %d.-iteration error: %f (%3.3f %%)\n', initialError, m, error, 100 * error/ initialError);

        %% Optionally, plot the movement and orientation of the pose during optimisation
        % (slow)
%         pose = result.atPose3(symbol('p', i));
%         pos = pose.translation();
%         mat = pose.matrix;
%         h = mat * [0;0;1;1]; % This point indicates the orientation of the orientation of the pose
%         h = [h(1:3)'; pos.x,pos.y,pos.z];
%         plot3(pos.x,pos.y,pos.z, 'bo');
%         plot3(h(:,1),h(:,2),h(:,3), 'b.');

        %% Optionally, plot the movement of features during optimisation
        % (slow)
        % subplot(2,1,2);
        % hold on;
        % plot_features(i, 'r.', pose, landmarks_gt);
        % subplot(2,1,1);

        % ********************* TODO ***********************
        % i kept only this condition as the value for the error changed a
        % lot in function of the cases I thought this condition was more
        % pertinant
        if last_error - error < eps
            break
        end
        
        last_error = error;
    end

    result = optimizer.values();
    error = graph.error(result);
    fprintf('Initial error: %f, Final error: %f (%3.3f %%)\n', initialError, error, 100 * error/ initialError);
    fprintf('Landmarks in window: %d\n', landmarks_in_window);

    % Read all poses that were optimised, 'pose' represents the most recent
    % frame afterwards
    for j=window_start:i
        pose = result.atPose3(symbol('p', j));
        pos = pose.translation();
        quat = pose.rotation().quaternion();
        camera_out(j,:) = [camera_gt(j,1) pos.x pos.y pos.z quat(2) quat(3) quat(4) quat(1)];
    end

    % Read landmarks
    for l=1:length(landmarks_used)
        if landmarks_used(l) == i
            mark = result.atPoint3(symbol('f', l));
            landmarks_out(l,:) = [mark.x mark.y mark.z];
        end
    end

    %% Plotting

    % Plot markers on last frame
    if PLOT_FEATURES
        subplot(2,1,2);
        hold on;
        plot_features(i, 'rx', pose, landmarks_out);
        subplot(2,1,1);
    end

    % Plot trajectory
    try delete(tplt); end
    tplt = plot3(camera_out(1:i,2),camera_out(1:i,3),camera_out(1:i,4), 'r*');

    % Plot all landmarks
    if PLOT_LANDMARKS
        try delete(mplt); end
        mplt = plot3(landmarks_out(:,1),landmarks_out(:,2),landmarks_out(:,3), ...
                     'Marker','.', 'Markersize',1, 'Color','black','LineStyle','none');
    end

    % Update figure and optionally save
    drawnow
    if SAVE_FIGURES
        saveas(gcf, strcat(outdir, sprintf('plot-%04d',i),'.png'));
    end

    % ************ TODO : save the result for the evaluation ************
    % RTE and ATE evaluation
    if RTE_ATE==1
        disp('Exporting results for evaluation')
        dlmwrite(strcat(outdir, 'vo_output_thisWindow.txt'), camera_out(window_start:i,:), 'delimiter','\t','precision',6);
        dlmwrite(strcat(outdir, 'vo_output_all.txt'), camera_out(1:i,:), 'delimiter','\t','precision',6);
        
        [thisAteAll, thisRteAll,thisAteWindow, thisRteWindow] = evaluate_rte_ate(window_start, i) ;
        
        rte_all(i,1) = thisRteAll ;
        ate_all(i,1) = thisAteAll ;
        
        rte_window(i,1) = thisRteWindow ;
        ate_window(i,1) = thisAteWindow ;
    end
    % ******************************************************************
    
    % **************** TODO : observe the behaviors of covariance information **************
    % store the standard deviation at time i 
    if covariance==1 
        % Optionally, retrieve marginals for plotting
        marginals = Marginals(graph, result);
        
        thisCov = marginals.marginalCovariance(symbol('p', i)) ;
        perAxisSD = diag(thisCov(4:6,4:6))' ;
        
        stand_dev(i,:) = perAxisSD ;
        
        if i == NUM_INITIALISE+1
            meanMagSD(i,1) = norm(perAxisSD) ;
        else
            meanMagSD(i,1) = (meanMagSD(i-NUM_INITIALISE-1,1)+norm(perAxisSD))/2 ;
        end
    end
    
    % store the covariance of all poses after the first bundle adjustment
    if covariance == 1 && firstBA == 1 
        for k = NUM_INITIALISE+1:i
            thisCov = marginals.marginalCovariance(symbol('p', k)) ;
            perAxisSD = diag(thisCov(4:6,4:6))' ;

            covAfterBA(k,:) = perAxisSD ;
        end
    end
    % **********************************************************
    
    elapsedTime = toc ;
    
    % **************** TODO : store the time **************
    if time==1
        timeArray(i,1) = elapsedTime ;
    end
    % **********************************************************
end

% Export results to hard-drive
disp('Exporting final results...');
dlmwrite(strcat(outdir, 'vo_output_poses_final.txt'), camera_out, 'delimiter','\t','precision',6);

end



%% Some helper functions:
function p = get_pose(matrix, index)
    import gtsam.*;
    p = Pose3(Rot3.Quaternion(matrix(index,8), matrix(index,5), matrix(index,6), matrix(index,7)), ...
              Point3(matrix(index,2), matrix(index,3), matrix(index,4)));
end

function p = get_landmark(matrix, index)
    import gtsam.*;
    p = Point3(matrix(index,1),matrix(index,2),matrix(index,3));
end

% Returns true if a landmark was previously mapped, otherwise false
function t = in_map(landmark_index)
    import gtsam.*;
    global landmarks_used;
    if landmarks_used(landmark_index) > 0
        t = true;
    else
        t = false;
    end
end

% Project 'landmarks' into the fhe frame with ID 'frame_id', using the 'camera_pose'
% If no pose or landmarks are provided, ground-truth feature-locations are plotted.
function plot_features(frame_id, style, camera_pose, landmarks)
    import gtsam.*;
    global calib;
    global features_gt;

    % Plot ground-truth if no pose or landmarks are provided
    plot_gt = false;
    if nargin < 3
        plot_gt = true;
    else
        camera = SimpleCamera(camera_pose, calib);
    end

    % Iterate features
    f = 1; % column of current feature ID
    show_warning = false;
    while f < size(features_gt, 2) && features_gt(frame_id,f) > 0
        if plot_gt
            plot(features_gt(frame_id,f+1),features_gt(frame_id,f+2),'g+');
        else
            feature_id = features_gt(frame_id,f);
            try
                point2d = camera.project(Point3(landmarks(feature_id,1),landmarks(feature_id,2),landmarks(feature_id,3)));
                plot(point2d.x,point2d.y,style);
            catch
                show_warning = true;
            end
        end
        f = f+4;
    end
    if show_warning
        warning('Could not plot all features.');
    end
end
