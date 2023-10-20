function handles = prepare_A(A, handles)
    % find incorrect rows with zero frequency
    index = find(A(:,1)==0); 
    A(index,:)=[];
    
    % flip freq, Z_prime and Z_double_prime so that data are in the desceding 
    % order of freq 
    if A(1,1) < A(end,1)
       A = fliplr(A')';
    end % if
    
    handles.freq = A(:,1);
    handles.Z_prime_mat = A(:,2);
    handles.Z_double_prime_mat = A(:,3);
    
    % save original freq, Z_prime and Z_double_prime
    handles.freq_0 = handles.freq;
    handles.Z_prime_mat_0 = handles.Z_prime_mat;
    handles.Z_double_prime_mat_0 = handles.Z_double_prime_mat;
    
    handles.Z_exp = handles.Z_prime_mat(:)+ 1i*handles.Z_double_prime_mat(:);
    
    handles.method_tag = 'none';

end % fun def