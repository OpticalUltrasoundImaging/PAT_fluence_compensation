function output = convex_cone(vt, p0_sNi)
% Normalization
vt = vt ./ norm(vt,2);
% error angle for all the sO2 values
ang = zeros(1,size(p0_sNi,1));
for i = 1:size(p0_sNi,1)
    cone_matrix = squeeze(p0_sNi(i,:,:));
    num_vectors = size(cone_matrix, 2);
    A = -eye(num_vectors); % Coefficients should be non-negative
    b = zeros(num_vectors, 1); % Right-hand side for the inequality constraints

    % Solve the linear system cone_matrix * alpha = vt
    options = optimoptions('linprog', 'Display', 'none');
    [alpha, ~, exitflag] = linprog([], A, b, cone_matrix, vt, [], [], options);

    if exitflag == 1
        ang(1,i) = 0;
    else
        options = optimoptions('quadprog', 'Display', 'none');
        % If vt is not in the cone, find the closest point on the cone
        % Solve the quadratic programming problem to minimize ||cone_matrix * alpha - vt||^2
        H = cone_matrix' * cone_matrix; % Quadratic term
        f = -cone_matrix' * vt; % Linear term
        [alpha_qp, ~, exitflag_qp] = quadprog(H, f, A, b, [], [], [], [], [], options);

        if exitflag_qp == 1
            closest_point = cone_matrix * alpha_qp;

            % Calculate the angle between vt and the closest point
            dot_product = dot(vt, closest_point);
            norm_vt = norm(vt);
            norm_closest_point = norm(closest_point);
            cosine_angle = dot_product / (norm_vt * norm_closest_point);
            angle_rad = acos(cosine_angle);
            angle_deg = rad2deg(angle_rad);
            ang(1,i) = angle_deg;
        else
            % disp('Failed to find the closest point on the cone.');
            % Need to be further fixed later (Jun 17, 2024)
            ang(1,i) = ang(1,i-1);
        end
    end
end


% zero elements
corrected_ang = ang;

% Assign corrected_ang to ang and make the prediction
[~, minIndex] = min(corrected_ang);
output = (minIndex-1) / 100;

end


