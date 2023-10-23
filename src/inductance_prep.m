function handles = inductance_prep(handles)

    if ~handles.data_exist
        return
    end

    switch handles.inductance
        case 1 %keep data fitting without inductance
            handles.freq = handles.freq_0;
            handles.Z_prime_mat = handles.Z_prime_mat_0;
            handles.Z_double_prime_mat = handles.Z_double_prime_mat_0; 

            handles.Z_exp = handles.Z_prime_mat(:)+ i*handles.Z_double_prime_mat(:);

        case 2 %keep data fitting with inductance
            handles.freq = handles.freq_0;
            handles.Z_prime_mat = handles.Z_prime_mat_0;
            handles.Z_double_prime_mat = handles.Z_double_prime_mat_0; 

            handles.Z_exp = handles.Z_prime_mat(:)+ i*handles.Z_double_prime_mat(:);

        case 3 %discard data
            is_neg = -handles.Z_double_prime_mat(:)<0;
            index = find(is_neg==1);
            handles.Z_double_prime_mat(index) = [];
            handles.Z_prime_mat(index) = [];
            handles.freq(index) = [];

    end
      
    handles.Z_exp = handles.Z_prime_mat(:)+ i*handles.Z_double_prime_mat(:);
    handles.method_tag = 'none';
    
    handles.taumax = ceil(max(log10(1./handles.freq)))+0.5;    
    handles.taumin = floor(min(log10(1./handles.freq)))-0.5;
    handles.freq_fine = logspace(-handles.taumin, -handles.taumax, 10*numel(handles.freq));

end % fun def