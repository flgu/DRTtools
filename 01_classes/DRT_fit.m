classdef DRT_fit < handle

    properties

        inductance
        rbf_type % used for fitting
        data_used % flag if input data is available
        lambda % Thikonov regulation parameter
        coeff % used for fitting
        shape_control % used for fitting
        der_used
        method_tag
        data_exist

        % input data
        freq
        Z_prime_mat
        Z_double_prime_mat

        % original data
        freq_0
        Z_prime_mat_0
        Z_double_prime_mat_0
        Z_exp

        % tau axis related data
        taumax
        taumin
        freq_fine

        % fit results
        gamma_ridge_fine
        mu_Z_re % recalculated EIS - Real based on DRT
        mu_Z_im % recalculated EIS - Imag based on DRT
    end % props

    methods
        function this = DRT_fit()        
            % assign initial values
            this.inductance = 1; % previously assigned via set fct
            this.rbf_type = 'Gaussian';
            this.data_used = 'Combined Re-Im Data';
            this.lambda = 1e-3;
            this.coeff = 0.5;
            this.shape_control = 'FWHM Coefficient';
            this.der_used = '1st-order';    
        %   method_tag: 'none': havnt done any computation, 'simple': simple DRT,
        %               'credit': Bayesian run, 'BHT': Bayesian Hibert run
            this.method_tag = 'none'; 
            

            % initial values for data
            this.freq = [];
            this.Z_prime_mat = [];
            this.Z_double_prime_mat = [];
            this.freq_0 = [];
            this.Z_prime_mat_0 = [];
            this.Z_double_prime_mat_0 = [];
            this.Z_exp = [];
            this.data_exist = false;

            % tau data
            this.taumax = [];
            this.taumin = [];
            this.freq_fine = [];

            % fit results
            this.gamma_ridge_fine = [];
            this.mu_Z_re = [];
            this.mu_Z_im = [];
        end % constructor

        function add_data(this, freq, realPart, imagPart)

            A = [freq,realPart,imagPart];

            % set flag value true
            this.data_exist = true;

            % find incorrect rows with zero frequency
            index = find(A(:,1)==0); 
            A(index,:)=[];
            
            % flip freq, Z_prime and Z_double_prime so that data are in the desceding 
            % order of freq 
            if A(1,1) < A(end,1)
                A = fliplr(A')';
            end % if
            
            this.freq = A(:,1);
            this.Z_prime_mat = A(:,2);
            this.Z_double_prime_mat = A(:,3);
            
            % save original freq, Z_prime and Z_double_prime
            this.freq_0 = this.freq;
            this.Z_prime_mat_0 = this.Z_prime_mat;
            this.Z_double_prime_mat_0 = this.Z_double_prime_mat;
            
            this.Z_exp = this.Z_prime_mat(:)+ 1i*this.Z_double_prime_mat(:);
            
            this.method_tag = 'none'; % this is already set in the initialization, but if you load new data, the method is restored ok

        end % fun def

        function inductance_prep(this)

            if ~this.data_exist
                return
            end
        
            switch this.inductance
                case 1 %keep data fitting without inductance
                    this.freq = this.freq_0;
                    this.Z_prime_mat = this.Z_prime_mat_0;
                    this.Z_double_prime_mat = this.Z_double_prime_mat_0; 
        
                    this.Z_exp = this.Z_prime_mat(:)+ i*this.Z_double_prime_mat(:);
        
                case 2 %keep data fitting with inductance
                    this.freq = this.freq_0;
                    this.Z_prime_mat = this.Z_prime_mat_0;
                    this.Z_double_prime_mat = this.Z_double_prime_mat_0; 
        
                    this.Z_exp = this.Z_prime_mat(:)+ i*this.Z_double_prime_mat(:);
        
                case 3 %discard data
                    is_neg = -this.Z_double_prime_mat(:)<0;
                    index = find(is_neg==1);
                    this.Z_double_prime_mat(index) = [];
                    this.Z_prime_mat(index) = [];
                    this.freq(index) = [];
        
            end
              
            this.Z_exp = this.Z_prime_mat(:)+ i*this.Z_double_prime_mat(:);
            this.method_tag = 'none';
            
            this.taumax = ceil(max(log10(1./this.freq)))+0.5;    
            this.taumin = floor(min(log10(1./this.freq)))-0.5;
            this.freq_fine = logspace(-this.taumin, -this.taumax, 10*numel(this.freq));
        end % fun def

        function ridge_regression(this)
            
            %   bounds ridge regression
            lb = zeros(numel(this.freq)+2,1);
            ub = Inf*ones(numel(this.freq)+2,1);
            x_0 = ones(size(lb));

            % Fitting options for quadprog
            options = optimset('algorithm','interior-point-convex',...
                                'Display','off',...
                                'TolFun',1e-15,...
                                'TolX',1e-10,...
                                'MaxFunEvals', 1E5);

            % warum hier nicht gleich die input real & imag benutzen?
            handles.b_re = real(this.Z_exp);
            handles.b_im = imag(this.Z_exp);
            
            %   compute epsilon
            handles.epsilon = compute_epsilon(this.freq, this.coeff, this.rbf_type, this.shape_control);

            %   calculate the A_matrix
            handles.A_re = assemble_A_re(this.freq, handles.epsilon, this.rbf_type);
            handles.A_im = assemble_A_im(this.freq, handles.epsilon, this.rbf_type);

            % up to here it works
            % below here not tested!
            %   adding the resistence column to the A_re_matrix
            handles.A_re(:,2) = 1;

            %   adding the inductance column to the A_im_matrix if necessary
            if  this.inductance==2
                handles.A_im(:,1) = 2*pi*(handles.freq(:));
            end

            %   calculate the M_matrix
            switch this.der_used
                case '1st-order'
                    handles.M = assemble_M_1(this.freq, handles.epsilon, this.rbf_type);
                case '2nd-order'
                    handles.M = assemble_M_2(this.freq, handles.epsilon, this.rbf_type);
            end

            %   Running ridge regression
            switch this.data_used
                case 'Combined Re-Im Data'
                    [H_combined,f_combined] = quad_format_combined(handles.A_re,...
                                                                    handles.A_im,...
                                                                    handles.b_re,...
                                                                    handles.b_im,...
                                                                    handles.M,...
                                                                    this.lambda);
                    handles.x_ridge = quadprog(H_combined,...
                                                f_combined,...
                                                [], [], [], [],...
                                                lb,...
                                                ub,...
                                                x_0,...
                                                options);

                    %prepare for HMC sampler
                    this.mu_Z_re = handles.A_re*handles.x_ridge;
                    this.mu_Z_im = handles.A_im*handles.x_ridge;

                    handles.res_re = this.mu_Z_re-handles.b_re;
                    handles.res_im = this.mu_Z_im-handles.b_im;

                    sigma_re_im = std([handles.res_re;handles.res_im]);

                    inv_V = 1/sigma_re_im^2*eye(numel(this.freq));

                    Sigma_inv = (handles.A_re'*inv_V*handles.A_re) + (handles.A_im'*inv_V*handles.A_im) + (this.lambda/sigma_re_im^2)*handles.M;
                    mu_numerator = handles.A_re'*inv_V*handles.b_re + handles.A_im'*inv_V*handles.b_im;
                    
                case 'Im Data'
                    [H_im,f_im] = quad_format(handles.A_im,...
                                            handles.b_im,...
                                            handles.M,...
                                            this.lambda);
                    handles.x_ridge = quadprog(H_im,...
                                                f_im,...
                                                [], [], [], [],...
                                                lb,...
                                                ub,...
                                                x_0,...
                                                options);

                    %prepare for HMC sampler
                    this.mu_Z_re = handles.A_re*handles.x_ridge;
                    this.mu_Z_im = handles.A_im*handles.x_ridge;

                    handles.res_im = handles.mu_Z_im-handles.b_im;
                    sigma_re_im = std(handles.res_im);

                    inv_V = 1/sigma_re_im^2*eye(numel(this.freq));

                    Sigma_inv = (handles.A_im'*inv_V*handles.A_im) + (handles.lambda/sigma_re_im^2)*handles.M;
                    mu_numerator = handles.A_im'*inv_V*handles.b_im;

                case 'Re Data'
                    [H_re,f_re] = quad_format(handles.A_re,...
                                            handles.b_re,...
                                            handles.M,...
                                            this.lambda);
                    handles.x_ridge = quadprog(H_re,...
                                                f_re,...
                                                [], [], [], [],...
                                                lb,...
                                                ub,...
                                                x_0,...
                                                options);

                    %prepare for HMC sampler
                    this.mu_Z_re = handles.A_re*handles.x_ridge;
                    this.mu_Z_im = handles.A_im*handles.x_ridge;

                    handles.res_re = handles.mu_Z_re-handles.b_re;
                    sigma_re_im = std(handles.res_re);

                    inv_V = 1/sigma_re_im^2*eye(numel(this.freq));

                    Sigma_inv = (handles.A_re'*inv_V*handles.A_re) + (this.lambda/sigma_re_im^2)*handles.M;            
                    mu_numerator = handles.A_re'*inv_V*handles.b_re;
            end % switch

            warning('off')
            handles.Sigma_inv = (Sigma_inv+Sigma_inv')/2;
            handles.mu = handles.Sigma_inv\mu_numerator; % linsolve
            warning('on')
            % map x to gamma

            [this.gamma_ridge_fine,this.freq_fine] = map_array_to_gamma(this.freq_fine,...
                                                                            this.freq,...
                                                                            handles.x_ridge(3:end),...
                                                                            handles.epsilon,...
                                                                            this.rbf_type);
            this.freq_fine = this.freq_fine';
        %   method_tag: 'none': havnt done any computation, 'simple': simple DRT,
        %               'credit': Bayesian run, 'BHT': Bayesian Hilbert run
            this.method_tag = 'simple';

            % deconvolved DRT callback displays DRT vs tau/f etc
            
        end % fun def

        function varargout = Gamma_Tau_plot(this)

            fig = figure();
            ax_DRT = axes(fig);

            %   Running ridge regression first
        
            switch this.method_tag
                case 'simple'
                    plot(ax_DRT, 1./this.freq_fine, this.gamma_ridge_fine, '-k', 'LineWidth', 3);
    
                    y_min = 0; 
                    y_max = max(this.gamma_ridge_fine);
    
                case 'credit'
                    ciplot(handles.lower_bound_fine, handles.upper_bound_fine, 1./handles.freq_fine, 0.7*[1 1 1]);%plot CI
                    hold on
                    plot(1./handles.freq_fine, handles.gamma_mean_fine, '-b', 'LineWidth', 3);
                    plot(1./handles.freq_fine, handles.gamma_ridge_fine, '-k', 'LineWidth', 3);
                    
                    %             my_str = sprintf('%g \% CI', 99);
                    h = legend('99\% CI', 'Mean', 'MAP', 'Location','NorthWest');
                    set(h,'Interpreter', 'LaTex','Fontsize', 24)
                    legend boxoff
    
                    y_min = 0; 
                    y_max = max(handles.upper_bound_fine);
    
                case 'BHT'
                    plot(1./handles.freq_fine, handles.gamma_mean_fine_re, '-b', 'LineWidth', 3);
                    hold on
                    plot(1./handles.freq_fine, handles.gamma_mean_fine_im, '-k', 'LineWidth', 3);
    
                    h = legend('Mean Re', 'Mean Im', 'Location','NorthWest');
                    set(h,'Interpreter', 'LaTex','Fontsize', 24)
                    legend boxoff
    
                    y_min = min([handles.gamma_mean_fine_re(:);handles.gamma_mean_fine_im(:)]);
                    y_max = max([handles.gamma_mean_fine_re(:);handles.gamma_mean_fine_im(:)]);
                
                case 'peak'
                    plot(1./handles.freq_fine, handles.gamma_ridge_fine, '-k', 'LineWidth', 3);
                    hold on
                    
                    for i = 1:handles.N_peak
                        plot(1./handles.freq_fine, handles.gamma_gauss_mat(:,i), 'LineWidth', 3);
                    end 
                    
                    y_min = 0; 
                    y_max = max(handles.gamma_ridge_fine);
                    
            end
    
            %   adding labels
            xlabel(ax_DRT,'$\tau/s$', 'Interpreter', 'Latex','Fontsize',24)
            ylabel(ax_DRT,'$\gamma(\ln\tau)/\Omega$','Interpreter', 'Latex','Fontsize',24);
            set(ax_DRT,...
                'xscale','log',...
                'xlim',[min(1./this.freq_fine), max(1./this.freq_fine)],...
                'ylim',[y_min, 1.1*y_max],...
                'Fontsize',20,...
                'xtick',10.^[-10:2:10],...
                'TickLabelInterpreter','latex')

            grid(ax_DRT, "on")
            grid(ax_DRT, "minor")
            
            % add optional outputs
            varargout{1} = ax_DRT;
            varargout{2} = fig;
        end % fun def

        function varargout = Gamma_Freq_plot(this, options)
            %   Running ridge regression
            
            arguments
                this
                options.Color = "black";
                options.LineWidth double = 3;
                options.LineStyle string = "-";
                options.DisplayName string = "DRT Spectrum";
            end % args

            fig = figure();

            fig_width = 1280; % px
            fig_height = 720; % px

            % populate properties
            fig.Position = [100,100,fig_width, fig_height];


            ax_DRT = axes(fig);
                
            switch this.method_tag
                case 'simple'
                    plt = plot(ax_DRT, this.freq_fine, this.gamma_ridge_fine );
        
                    y_min = 0; 
                    y_max = max(this.gamma_ridge_fine);
        
                case 'credit'
                    ciplot(handles.lower_bound_fine, handles.upper_bound_fine, handles.freq_fine, 0.7*[1 1 1]);%plot CI
                    hold on
                    plot(handles.freq_fine, handles.gamma_mean_fine, '-b', 'LineWidth', 3);
                    plot(handles.freq_fine, handles.gamma_ridge_fine, '-k', 'LineWidth', 3);
                    
                %             my_str = sprintf('%g \% CI', 99);
                    h = legend('99\% CI', 'Mean', 'MAP', 'Location','NorthWest');
                    set(h,'Interpreter', 'LaTex','Fontsize', 24)
                    legend boxoff
        
                    y_min = 0; 
                    y_max = max(handles.upper_bound_fine);
        
                case 'BHT'
                    plot(handles.freq_fine, handles.gamma_mean_fine_re, '-b', 'LineWidth', 3);
                    hold on
                    plot(handles.freq_fine, handles.gamma_mean_fine_im, '-k', 'LineWidth', 3);
        
                    h = legend('Mean Re', 'Mean Im', 'Location','NorthWest');
                    set(h,'Interpreter', 'LaTex','Fontsize', 24)
                    legend boxoff
        
                    y_min = min([handles.gamma_mean_fine_re(:);handles.gamma_mean_fine_im(:)]);
                    y_max = max([handles.gamma_mean_fine_re(:);handles.gamma_mean_fine_im(:)]);
                    
                case 'peak'
                    plot(handles.freq_fine, handles.gamma_ridge_fine, '-k', 'LineWidth', 3);
                    hold on
                        
                    for i = 1:handles.N_peak
                        plot(handles.freq_fine, handles.gamma_gauss_mat(:,i), 'LineWidth', 3);
                    end 
                    
                    y_min = 0; 
                    y_max = max(handles.gamma_ridge_fine);
                    
            end
            
            plt.LineWidth = options.LineWidth;
            plt.LineStyle = options.LineStyle;
            plt.DisplayName = options.DisplayName;
            plt.Color = options.Color;

            %   adding labels
            xlabel(ax_DRT,'$f$/Hz', 'Interpreter', 'Latex','Fontsize',24)
            ylabel(ax_DRT,'$\gamma(\ln f)/\Omega$','Interpreter', 'Latex','Fontsize',24);
            set(ax_DRT,...
            'xscale','log',...
            'xlim',[min(this.freq_fine), max(this.freq_fine)],...
            'ylim',[y_min, 1.1*y_max],...
            'Fontsize',20,...
            'xtick',10.^[-10:2:10],...
            'TickLabelInterpreter','latex')
            hold off

            grid(ax_DRT, "on")
            grid(ax_DRT, "minor")
            
            % add optional outputs
            varargout{1} = ax_DRT;
            varargout{2} = fig;
        end % fun def


    end % methods
end % class def