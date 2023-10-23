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
    % below here not tested!
    %   adding the resistence column to the A_re_matrix
    handles.A_re(:,2) = 1;

    %   adding the inductance column to the A_im_matrix if necessary
    if  handles.inductance==2
        handles.A_im(:,1) = 2*pi*(handles.freq(:));
    end

    %   calculate the M_matrix
    switch handles.der_used
        case '1st-order'
            handles.M = assemble_M_1(handles.freq, handles.epsilon, handles.rbf_type);
        case '2nd-order'
            handles.M = assemble_M_2(handles.freq, handles.epsilon, handles.rbf_type);
    end

    %   Running ridge regression
    switch handles.data_used
        case 'Combined Re-Im Data'
            [H_combined,f_combined] = quad_format_combined(handles.A_re, handles.A_im, handles.b_re, handles.b_im, handles.M, handles.lambda);
            handles.x_ridge = quadprog(H_combined, f_combined, [], [], [], [], handles.lb, handles.ub, handles.x_0, handles.options);

            %prepare for HMC sampler
            handles.mu_Z_re = handles.A_re*handles.x_ridge;
            handles.mu_Z_im = handles.A_im*handles.x_ridge;

            handles.res_re = handles.mu_Z_re-handles.b_re;
            handles.res_im = handles.mu_Z_im-handles.b_im;

            sigma_re_im = std([handles.res_re;handles.res_im]);

            inv_V = 1/sigma_re_im^2*eye(numel(handles.freq));

            Sigma_inv = (handles.A_re'*inv_V*handles.A_re) + (handles.A_im'*inv_V*handles.A_im) + (handles.lambda/sigma_re_im^2)*handles.M;
            mu_numerator = handles.A_re'*inv_V*handles.b_re + handles.A_im'*inv_V*handles.b_im;
            
        case 'Im Data'
            [H_im,f_im] = quad_format(handles.A_im, handles.b_im, handles.M, handles.lambda);
            handles.x_ridge = quadprog(H_im, f_im, [], [], [], [], handles.lb, handles.ub, handles.x_0, handles.options);

            %prepare for HMC sampler
            handles.mu_Z_re = handles.A_re*handles.x_ridge;
            handles.mu_Z_im = handles.A_im*handles.x_ridge;

            handles.res_im = handles.mu_Z_im-handles.b_im;
            sigma_re_im = std(handles.res_im);

            inv_V = 1/sigma_re_im^2*eye(numel(handles.freq));

            Sigma_inv = (handles.A_im'*inv_V*handles.A_im) + (handles.lambda/sigma_re_im^2)*handles.M;
            mu_numerator = handles.A_im'*inv_V*handles.b_im;

        case 'Re Data'
            [H_re,f_re] = quad_format(handles.A_re, handles.b_re, handles.M, handles.lambda);
            handles.x_ridge = quadprog(H_re, f_re, [], [], [], [], handles.lb, handles.ub, handles.x_0, handles.options);

            %prepare for HMC sampler
            handles.mu_Z_re = handles.A_re*handles.x_ridge;
            handles.mu_Z_im = handles.A_im*handles.x_ridge;

            handles.res_re = handles.mu_Z_re-handles.b_re;
            sigma_re_im = std(handles.res_re);

            inv_V = 1/sigma_re_im^2*eye(numel(handles.freq));

            Sigma_inv = (handles.A_re'*inv_V*handles.A_re) + (handles.lambda/sigma_re_im^2)*handles.M;            
            mu_numerator = handles.A_re'*inv_V*handles.b_re;
    end % switch

    warning('off')
    handles.Sigma_inv = (Sigma_inv+Sigma_inv')/2;
    handles.mu = handles.Sigma_inv\mu_numerator; % linsolve
    warning('on')
    % map x to gamma

    [handles.gamma_ridge_fine,handles.freq_fine] = map_array_to_gamma(handles.freq_fine, handles.freq, handles.x_ridge(3:end), handles.epsilon, handles.rbf_type);
     handles.freq_fine = handles.freq_fine';
%   method_tag: 'none': havnt done any computation, 'simple': simple DRT,
%               'credit': Bayesian run, 'BHT': Bayesian Hilbert run
    handles.method_tag = 'simple';

    % deconvolved DRT callback displays DRT vs tau/f etc
    
end % fun def