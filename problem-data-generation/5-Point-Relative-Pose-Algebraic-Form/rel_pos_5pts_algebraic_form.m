%> Generate parameters and the corresponding solutions of
%  5-point relative pose estimation problem (Algebraic form)
%
%> Code Description: 
%     Use synthetic curve datasets, generate the target parameters and the 
%     ground truth data. These data are then used by the 
%     evaluation_5pt_rel_pose_algebraic.m function file for verifications.
%
%> (c) LEMS, Brown University
%> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
%> Nov. 16th, 2022

clear;
clc;
close all;

rng(0);

%> define path to synthetic curve dataset
src_dataset_path = '/home/chchien/datasets/synthcurves_dataset/spherical-ascii-100_views-perturb-radius_sigma10-normal_sigma0_01rad-minsep_15deg-no_two_cams_colinear_with_object/';

%> define experiment parameters
noise = 0;
view_indx = [80, 82];
selectCurveNum = 5;

%> define rotation variable representations
%rotation_representation = 'rotation_matrix';
rotation_representation = 'quaternion';
%rotation_representation = 'cayley';

%> valid curve ids only covers the ids from 3~39
validCurveIds = [];
for i = 3:39
    validCurveIds = [validCurveIds, i];
end
%sz = [size(validCurveIds, 2), 3];
%curve_RGB_color = unifrnd(0,1,sz);

%> pick random three curve ids
% rndCurveIds = round(unifrnd(3, 39, [1, selectCurveNum]));
% CurveIds = rndCurveIds;
assignCurveIds = [10, 12, 14, 16, 18];
CurveIds = assignCurveIds;

%> initialize: (dim1, dim2, dim3) is (points, point coordinates, views)
pickPoints2D   = zeros(5, 2, 2);

%> initialize 3D points and tangents
pickPoints3D   = zeros(5, 3);

%: stack all absolute camera extrinsic parameters from three views
abs_pose_R = zeros(3,3,2);
abs_pose_T = zeros(3,1,2);
abs_pose_C = zeros(3,1,2);

%: define a parameter array
parameters = zeros(20,1);

%: transform to metric representation
metric_pts = zeros(5, 2, 2);

%> change the point index to generate different set of target parameters
%pickPointsIndx = [90, 97, 56];
%pickPointsIndx = [10, 17, 55];
%pickPointsIndx = [90, 97, 56, 50, 29];
rndPickPointsIndx = [];

%> read image data, scene data, and camera parameters from the dataset
%  by looping over 2 views and select samples from 5 different curves
for v = 1:2
    for p = 1:5
        %> import image data, scene data, from a view and a curve
        [imgData_curve, sceneData_curve, Params] = readCurveSyntheticDataset(view_indx(v), src_dataset_path, selectCurveNum, CurveIds(p));
        
        %> from the imported data, rotation_representation = 'rotation_matrix';pick one sample randomly and store in
        %  arrametric_ys (tangents are already normalized)
        dataSizeCurve = size(imgData_curve.points, 1);
        
        %fprintf(string(rndPickPointsIndx));
        %fprintf('\t');
        %rndPickPointsIndx = pickPointsIndx(p);
        
        %: for 3D points and tangents and the 3D points along their tangents, only the first view is required
        if v == 1
            idx = round(unifrnd(1, dataSizeCurve, [1, 1]));
            rndPickPointsIndx = [rndPickPointsIndx, idx];
            pickPoints3D(p,:) = sceneData_curve.points(idx, :);
        end
        pickPoints2D(p,:,v) = imgData_curve.points(rndPickPointsIndx(p), :);

        %> the homogenous 2D image point should be padded with 1
        aug_pt2d = [pickPoints2D(p,:,v), 1];
        inv_K_pts2d = inv(Params.K) * aug_pt2d';
        metric_pts(p,:,v) = inv_K_pts2d(1:2,1)';
        
        %> pick up the depth of the picked 2D point, i.e., alpha
        unnorm_p = Params.K * [Params.R, Params.T] * [pickPoints3D(p,:), 1]';
        img_pt = unnorm_p;
        img_pt(:,1) = img_pt(:,1) ./ img_pt(3,1);
        if img_pt(1,1) - pickPoints2D(p,1,v) > 0.001 || img_pt(2,1) - pickPoints2D(p,2,v) > 0.001
            fprintf("invalid projection of point ");
            fprintf(string(p));
            fprintf(" in view ");
            fprintf(string(v));
            fprintf("\n");
        end
        
        %> store the depths
        alpha(p,v) = unnorm_p(3,1);        
    end
    
    abs_pose_R(:,:,v) = Params.R;
    abs_pose_T(:,:,v) = Params.T;
    abs_pose_C(:,:,v) = Params.C;
end

abs_T1 = abs_pose_T(:,:,1);
abs_T2 = abs_pose_T(:,:,2);
abs_R1 = abs_pose_R(:,:,1);
abs_R2 = abs_pose_R(:,:,2);

%: relative pose from view 1 to view 2
R12 = abs_R2 * abs_R1';
T12 = -abs_R2 * abs_R1' * abs_T1 + abs_T2;

%> TODO: if use Cayley, the input parameters have to be changed.
%> Convert from a general rotation matrix to a cayley
%> 1) (optional)define cayley representation
cay2R = @(r) [1+r(1,1)^2-(r(2,1)^2+r(3,1)^2), 2*(r(1,1)*r(2,1)-r(3,1)),       2*(r(1,1)*r(3,1)+r(2,1)); ...
              2*(r(1,1)*r(2,1)+r(3,1)),       1+r(2,1)^2-(r(1,1)^2+r(3,1)^2), 2*(r(2,1)*r(3,1)-r(1,1)); ...
              2*(r(1,1)*r(3,1)-r(2,1)),       2*(r(2,1)*r(3,1)+r(1,1)),       1+r(3,1)^2-(r(1,1)^2+r(2,1)^2)
             ];
%> 2) convert from rotation matrix to the exponential coordinates representation of the orientation
skew_r12 = inv(R12+eye(3))*(R12-eye(3));
r12 = [skew_r12(3,2); skew_r12(1,3); skew_r12(2,1)];
%denom_cayley_r12 = sqrt(1 + r12(1,1)^2 + r12(2,1)^2 + r12(3,1)^2);
%r12 = r12 ./ denom_cayley_r12;
%> (optional)
R12_cayley = cay2R(r12);

%> Quarternion parametrization
% quat2R = @(qx, qy, qz, qw)[1-2*qy*qy-2*qz*qz, 2*qx*qy-2*qw*qz,   2*qx*qz+2*qw*qy;
%                            2*qx*qy+2*qw*qz,   1-2*qx*qx-2*qz*qz, 2*(qy*qz-qw*qx);
%                            2*(qx*qz-qw*qy),   2*(qy*qz+qw*qx),   1-2*qx*qx-2*qy*qy];
% R12 = quat2R(q2_x, q2_y, q2_z, q2_w);
q_R12 = rotm2quat(R12);
%> (w, x, y, z)

unit_T12 = T12 ./ norm(T12);
normalize_T12 = T12 ./ T12(3);

%metric_pts(p,dim,v)
parameters(1:2,1)   = metric_pts(1,:,1).';
parameters(3:4,1)   = metric_pts(1,:,2).';
parameters(5:6,1)   = metric_pts(2,:,1).';
parameters(7:8,1)   = metric_pts(2,:,2).';
parameters(9:10,1)  = metric_pts(3,:,1).';
parameters(11:12,1) = metric_pts(3,:,2).';
parameters(13:14,1) = metric_pts(4,:,1).';
parameters(15:16,1) = metric_pts(4,:,2).';
parameters(17:18,1) = metric_pts(5,:,1).';
parameters(19:20,1) = metric_pts(5,:,2).';

if strcmp(rotation_representation, 'rotation_matrix')
    %> TODO
elseif strcmp(rotation_representation, 'quaternion')
    solutions = [normalize_T12(1), normalize_T12(2), ...
                 q_R12(1,2), q_R12(1,3), q_R12(1,4), q_R12(1,1)
                ];
elseif strcmp(rotation_representation, 'cayley')
    %> TODO
end

%> plug in parameters and variables into polynomials for correctness evaluations
%evaluation_5pt_rel_pose(solutions, parameters', rotation_representation, 1e-8);

evaluation_5pt_rel_pose_algebraic(solutions, parameters', rotation_representation, 1e-8);





