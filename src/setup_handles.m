function handles = setup_handles()
    
    % set(handles.dis_button,'Value',1)
    % set(handles.plot_pop,'Value',1)
    % set(handles.derivative,'Value',1)
    % set(handles.shape,'Value',1)
    % set(handles.value,'String','1E-3')
    % set(handles.coef,'String','0.5')
    % set(handles.inductance,'Value',1)
    % set(handles.panel_drt, 'Visible', 'on');
    % set(handles.running_signal, 'Visible', 'off');


    handles.inductance = 1; % previously assigned via set fct
    handles.rbf_type = 'Gaussian';
    handles.data_used = 'Combined Re-Im Data';
    handles.lambda = 1e-3;
    handles.coeff = 0.5;
    handles.shape_control = 'FWHM Coefficient';
    handles.der_used = '1st-order';    
%   method_tag: 'none': havnt done any computation, 'simple': simple DRT,
%               'credit': Bayesian run, 'BHT': Bayesian Hibert run
    handles.method_tag = 'none'; 
    handles.data_exist = false;
end % fun def