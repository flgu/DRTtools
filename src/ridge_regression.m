function handles = ridge_regression(handles)

    %   bounds ridge regression
    handles.lb = zeros(numel(handles.freq)+2,1);
    handles.ub = Inf*ones(numel(handles.freq)+2,1);
    handles.x_0 = ones(size(handles.lb));

    handles.options = optimset('algorithm','interior-point-convex','Display','off','TolFun',1e-15,'TolX',1e-10,'MaxFunEvals', 1E5);

    handles.b_re = real(handles.Z_exp);
    handles.b_im = imag(handles.Z_exp);
    
    %   compute epsilon
    handles.epsilon = compute_epsilon(handles.freq, handles.coeff, handles.rbf_type, handles.shape_control);

    %   calculate the A_matrix
    handles.A_re = assemble_A_re(handles.freq, handles.epsilon, handles.rbf_type);
    handles.A_im = assemble_A_im(handles.freq, handles.epsilon, handles.rbf_type);

    % up to here it works
end % fun def