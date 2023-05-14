function evaluation_5pt_rel_pose_geometric(solutions, parameters, rotation_representation, tol)

    %: 5-point relative pose estimation
    %> 1) rotations
    if strcmp(rotation_representation, 'rotation_matrix')
        %> (9) variables for a general rotation matrix
        syms r2_11 r2_12 r2_13 r2_21 r2_22 r2_23 r2_31 r2_32 r2_33
    elseif strcmp(rotation_representation, 'quaternion')
        %> (4) varaibles for a quaternion
        syms q2_w q2_x q2_y q2_z
    elseif strcmp(rotation_representation, 'cayley')
        %> (3) variables for a cayley representation
        syms r12_x r12_y r12_z
    end
    
    %> 2) translations (3)
    syms t12_x t12_y
    assume(t12_x, 'real');
    assume(t12_y, 'real');

    %> polynomial parameters 5*2+5*2=20
    %> the 2D image point pairs (20), point(point,view)
    px = sym('px', [5,2], 'real')';
    py = sym('py', [5,2], 'real')';
    
%     syms point11_x point11_y %> 1st point of first view
%     syms point12_x point12_y %> 1st point of second view
%     syms point21_x point21_y %> 2nd point of first view
%     syms point22_x point22_y %> 2nd point of the second view
%     syms point31_x point31_y %> 3rd point of first view
%     syms point32_x point32_y %> 3rd point of second view
%     syms point41_x point41_y %> 4th point of first view
%     syms point42_x point42_y %> 4th point of second view
%     syms point51_x point51_y %> 5th point of first view
%     syms point52_x point52_y %> 5th point of second view
    
    %> 1) the charts (3)
    syms rho_p11
    %> points in homogeneous form
    points2D_v1 = [px(1,1), py(1,1), 1;     %> 1st point of first view
                   px(1,2), py(1,2), 1;     %> 2nd point of first view
                   px(1,3), py(1,3), 1;
                   px(1,4), py(1,4), 1;
                   px(1,5), py(1,5), 1];
    points2D_v2 = [px(2,1), py(2,1), 1;     %> 1st point of second view
                   px(2,2), py(2,2), 1;     %> 2nd point of second view
                   px(2,3), py(2,3), 1; 
                   px(2,4), py(2,4), 1; 
                   px(2,5), py(2,5), 1];

    %> rotations
    if strcmp(rotation_representation, 'rotation_matrix')
        R12 = [r2_11, r2_12, r2_13; r2_21, r2_22, r2_23; r2_31, r2_32, r2_33];
    elseif strcmp(rotation_representation, 'quaternion')
        quat2R = @(qx, qy, qz, qw)[1-2*qy*qy-2*qz*qz, 2*qx*qy-2*qw*qz,   2*qx*qz+2*qw*qy;
                                   2*qx*qy+2*qw*qz,   1-2*qx*qx-2*qz*qz, 2*(qy*qz-qw*qx);
                                   2*(qx*qz-qw*qy),   2*(qy*qz+qw*qx),   1-2*qx*qx-2*qy*qy];
        R12 = quat2R(q2_x, q2_y, q2_z, q2_w);
    elseif strcmp(rotation_representation, 'cayley')
        %> define caley representation
        cay2R = @(r) [1+r(1,1)^2-(r(2,1)^2+r(3,1)^2), 2*(r(1,1)*r(2,1)-r(3,1)),       2*(r(1,1)*r(3,1)+r(2,1)); ...
                      2*(r(1,1)*r(2,1)+r(3,1)),       1+r(2,1)^2-(r(1,1)^2+r(3,1)^2), 2*(r(2,1)*r(3,1)-r(1,1)); ...
                      2*(r(1,1)*r(3,1)-r(2,1)),       2*(r(2,1)*r(3,1)+r(1,1)),       1+r(3,1)^2-(r(1,1)^2+r(2,1)^2)
                     ];
        R12 = cay2R([r12_x; r12_y; r12_z]);
    end
    
    %> translation
    t12 = [t12_x; t12_y; 1];
    
    %> POINT EQUATIONS
    %> 1) 5 pairs of correspondences across two views
    Eq_p1 = dot(t12, cross(points2D_v2(1,:).', R12*points2D_v1(1,:).'));
    Eq_p2 = dot(t12, cross(points2D_v2(2,:).', R12*points2D_v1(2,:).'));
    Eq_p3 = dot(t12, cross(points2D_v2(3,:).', R12*points2D_v1(3,:).'));
    Eq_p4 = dot(t12, cross(points2D_v2(4,:).', R12*points2D_v1(4,:).'));
    Eq_p5 = dot(t12, cross(points2D_v2(5,:).', R12*points2D_v1(5,:).'));
    Eqs = [Eq_p1; Eq_p2; Eq_p3; Eq_p4; Eq_p5];
    
    quat_norm_constraint = q2_w^2 + q2_x^2 + q2_y^2 + q2_z^2 - 1;
    %Eqs = [Eqs; quat_norm_constraint];
    
    %: substitue values to the symbols
    %collect(Iterators.flatten([permutedims(x,[3,2,1]), permutedims(d,[3,2,1]), a[1,1], e[1,1], e[2,1]]))
    sym_parameters = [px(1,1), py(1,1), ...
                      px(2,1), py(2,1), ...
                      px(1,2), py(1,2), ...
                      px(2,2), py(2,2), ...
                      px(1,3), py(1,3), ...
                      px(2,3), py(2,3), ...
                      px(1,4), py(1,4), ...
                      px(2,4), py(2,4), ...
                      px(1,5), py(1,5), ...
                      px(2,5), py(2,5)
                     ];
    
    if strcmp(rotation_representation, 'rotation_matrix')
        %> TODO
    elseif strcmp(rotation_representation, 'quaternion')
        sym_variables = [t12_x, t12_y, ...
                         q2_x, q2_y, q2_z, q2_w
                        ];
    elseif strcmp(rotation_representation, 'cayley')
        %> TODO
    end

    invalid_flag = 0;
    %: Verify the RHS of the polynomials given the start parameters and the start solutions
    %> 1) point equations (5)
    for i = 1:5
        eq_val = subs(Eqs(i,1), [sym_parameters, sym_variables], [parameters, solutions]);
        if abs(eq_val) > tol 
            fprintf('point equations for v1 & v2 invalid solution at equation ');
            fprintf(string(i));
            fprintf('\n');
            invalid_flag = 1;
        end
    end
    
    %> 2) quaternion normalization constraint equation
    eq_val = subs(quat_norm_constraint, [sym_parameters, sym_variables], [parameters, solutions]);
    if abs(eq_val) > tol 
        fprintf('invalid solution in the quaternion normalization constraint\n');
        invalid_flag = 1;
    end
    
    %> if the evaluation is correct then print out the message
    if ~invalid_flag
        fprintf("Evaluation complete! Parameters and variables fit the polynomials!\n");
    end
end


