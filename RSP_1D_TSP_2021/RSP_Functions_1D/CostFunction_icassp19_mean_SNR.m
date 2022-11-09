


function [J_N] = CostFunction_Proposed(params, a_K, m_K, x_K, x_N, m_N )

rho_n = params.rho_n;

sig_n = params.sig_n;

rho_a = params.rho_a;

sig_a = params.sig_a;


eps = 1e-10;


C_KK = Cov ( rho_a, sig_a, x_K, x_K ) ;
L_C_KK = chol(C_KK + eps*eye(size(C_KK)), 'lower');


G_KK = Cov ( rho_n, sig_n, x_K, x_K );
L_G_KK = chol(G_KK + eps*eye(size(G_KK)), 'lower');

N = size(x_N,1);

%%  Loop for caculating the cost function at each new postion (x_N)

J_N = zeros(N,1);

for i = 1:N
    
    x_Ni = x_N(i);
    m_Ni = m_N(i);
    
C_NK = Cov ( rho_a, sig_a, x_Ni, x_K );
G_NK = Cov ( rho_n, sig_n, x_Ni, x_K );


C_NN = Cov ( rho_a, sig_a, x_Ni, x_Ni ) ;
G_NN = Cov ( rho_n, sig_n, x_Ni, x_Ni ) ;


% To calculate R's

R_NN = pinv( ( G_NN - G_NK * (  L_G_KK'\ ( L_G_KK \ G_NK' )  ) ) ) ;
R_KN = -  L_G_KK.' \ ( L_G_KK \ (G_NK'  * R_NN) );
R_KK = (  (eye(size(G_KK)) - R_KN * G_NK ) / L_G_KK.'  ) / L_G_KK ;


% To calculate a_hat_N and C_hat_N

a_hat_N = m_Ni + C_NK * (  L_C_KK.'  \   ( L_C_KK \ ( a_K - m_K ) )  );

C_hat_N = C_NN - C_NK * (  L_C_KK.'   \   (  L_C_KK \ C_NK' )   );


% To calculate the cost function

J_N(i) = a_K' * R_KK * a_K + 2 * a_K' * R_KN * a_hat_N + trace ( R_NN * ( C_hat_N + a_hat_N * a_hat_N' ) );

end


end



