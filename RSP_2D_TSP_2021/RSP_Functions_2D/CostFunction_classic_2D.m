


function [J_N] = CostFunction_classic_2D(params, x_K, x_N)


rho_a = params.rho_a;
sig_a = params.sig_a;

rho_n = params.rho_n;
sig_n = params.sig_n;

N = size(x_N,1);


%  Loop for caculating the cost function at each new postion (x_N)

J_N = zeros(N,1);


for i = 1 : N
    
    i;
    
    x_Ni = x_N(i,:);

    % To calculate the cost function
    
%     x_K_new = sort ([x_K ; x_Ni]);     % selected
    x_K_new = x_K;     % selected

    x_N_new =  x_Ni;      % rest
    
    
    C_NK1 = My_2D_Cov ( rho_a, sig_a, x_N_new, x_K_new );
    
    C_NN1 = My_2D_Cov ( rho_a, sig_a, x_N_new, x_N_new ) ;
    
    C_KK1 = My_2D_Cov ( rho_a, sig_a, x_K_new, x_K_new ) + My_2D_Cov ( rho_n, sig_n, x_K_new, x_K_new );
    L_C_KK1 = chol(C_KK1 + .5e-6*eye(size(C_KK1)), 'lower');
    
    
    C_hat_N1 = C_NN1 - C_NK1 * (  L_C_KK1.'   \   (  L_C_KK1 \ C_NK1' )   );
    
    J_N (i) = sum( diag(C_hat_N1) )/N;
   
end

1


end








