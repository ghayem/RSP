
function [J_y] = CostFunction_MI_2D(params, x_K, x_rest)

rho_a = params.rho_a;
sig_a = params.sig_a;
rho_n = params.rho_n;
sig_n = params.sig_n;
N = size(x_rest,1);
A = x_K;
%  Loop for caculating the cost function at each new postion (x_N)
J_y = zeros(N,1);

parfor i = 1 : N
    A_bar = x_rest;
    y_i = A_bar(i,:);
    A_bar(i,:)=[];
    
    C_yA = My_2D_Cov ( rho_a, sig_a, y_i, A );
    C_yy = My_2D_Cov ( rho_a, sig_a, y_i, y_i ) ;
    C_AA = My_2D_Cov ( rho_a, sig_a, A, A ) + My_2D_Cov ( rho_n, sig_n, A, A );
    L_C_AA = chol(C_AA + 5e-6*eye(size(C_AA)), 'lower');
    C1 = C_yy - C_yA * (  L_C_AA.'   \   (  L_C_AA \ C_yA' )   );
    %%%%
    C_yA_bar = My_2D_Cov ( rho_a, sig_a, y_i, A_bar );
    C_AA_bar = My_2D_Cov ( rho_a, sig_a, A_bar, A_bar ) + My_2D_Cov ( rho_n, sig_n, A_bar, A_bar ) ;
    L_C_AA_bar = chol(C_AA_bar + .5e-6*eye(size(C_AA_bar)), 'lower');
    C2 = C_yy - C_yA_bar * (  L_C_AA_bar.'   \   (  L_C_AA_bar \ C_yA_bar' )   );
    %%%%
    J_y (i) = C1 / C2;
   
end

% J_y=1

end








