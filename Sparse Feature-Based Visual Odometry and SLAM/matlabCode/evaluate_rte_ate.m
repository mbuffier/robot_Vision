function [ate_all, rte_all,ate_window, rte_window] = evaluate_rte_ate (startFrame, endFrame)

% extraction of the good value of camera poses 
createFile(startFrame, endFrame) ;

% RTE and ATE computed from the start of the run to the current time
[~,ate_all] = system('python ../external/tum-evaluation/evaluate_ate.py ../blender/camera_poses_all.txt output/vo_output_all.txt --plot allCam_ate.pdf') ;
[~,rte_all]  = system('python ../external/tum-evaluation/evaluate_rpe.py ../blender/camera_poses_all.txt output/vo_output_all.txt --plot allCam_rpe.pdf --fixed_delta') ;

ate_all = str2double(ate_all) ;
rte_all = str2double(rte_all) ;

% RTE and ATE computed over just the length of the sliding window
[~,ate_window] = system('python ../external/tum-evaluation/evaluate_ate.py ../blender/camera_poses_window.txt output/vo_output_thisWindow.txt --plot windowCam_ate.pdf') ;
[~,rte_window] = system('python ../external/tum-evaluation/evaluate_rpe.py ../blender/camera_poses_window.txt output/vo_output_thisWindow.txt --plot windowCam_rpe.pdf --fixed_delta') ;

ate_window = str2double(ate_window) ;
rte_window = str2double(rte_window) ;
end


function r = createFile(startFrame, endFrame)
% read the file with all ground truth positions 
allData = dlmread('../blender/camera_poses.txt') ;

% extract the wanted position information
matrixAll = allData(1:endFrame, :) ;
matrixThisWindow = allData(startFrame:endFrame,:) ;

% write files with the wanted positions 
dlmwrite('../blender/camera_poses_all.txt', matrixAll, 'delimiter','\t','precision',6);
dlmwrite('../blender/camera_poses_window.txt', matrixThisWindow, 'delimiter','\t','precision',6);
end

