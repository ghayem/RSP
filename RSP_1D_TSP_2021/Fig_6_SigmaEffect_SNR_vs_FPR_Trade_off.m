
clear all;
close all;
clc;

%%%   Generate Grid
K = 100;  % # of possible total sensors
range = 1;  % Grid range
x_K = linspace(0, range, K)';   % places of possible total sensors

tol = 1e-6;

%%% spatial gain 'a'
rho_a = .1;
sig_a = 1;


%%% noise
%SNR
% SNR_dB = 1.5;
% SNR = (10^SNR_dB) / 10;
% sig_n = sig_a * S_t / sqrt(SNR);
sig_n = 1;
ratio_roh = 2e-1;
rho_n = ratio_roh * rho_a;


%%% Source
S_t = 1;

%%% Uncertainty
% nbSig_uncertainty = 15;
% s_uncertainty_set = logspace(-2, 0, nbSig_uncertainty).';
s_uncertainty_set = [0:0.1:1]';
nbSig_uncertainty=size(s_uncertainty_set,1);
l_uncertainty = rho_a;
flag_uncertainty = 1;

%%%% Bias
l_bias = rho_a;
s_bias = 0.3;
flag_bias = 0;


% Initial sensors
ind_Int = [50];
K_Int = size(ind_Int,1);
x_Int = x_K(ind_Int);

n = 1;  % # of mixture in Mixtured Gaussian Process
p = 5;  % number of sensors that we want to select greedy



%%% Delta parameter
jj_max=1;
j_step=20;
delta = j_step * [1:jj_max];

Kn = Cov ( rho_n, sig_n, x_K, x_K);
Ka = Cov ( rho_a, sig_a, x_K, x_K);

k = @(x,xP,s2,l) ( (s2^2) .* exp(-(x-xP).^2 / (2*l^2)));

mont_a_max = 10;
mont_m_max = 10;

for mont_sig_uncertainty = 1 : nbSig_uncertainty
    
s_uncertainty = s_uncertainty_set(mont_sig_uncertainty) ;


for mont_a = 1:mont_a_max
    
    
    a_K = Generate_GM(rho_a, sig_a, x_K, 0, n) ;
    
    
    for mont_m = 1 : mont_m_max
        
        mont_m;
        
        % Uncertainty
        uncertainty_a = Generate_GM(l_uncertainty, s_uncertainty, x_K, 0, 1);

        % Bias
        bias = Generate_GM(l_bias, s_bias, x_K, 0, 1);
        
        m_K = a_K +  (flag_uncertainty * uncertainty_a) +  (flag_bias * bias);
        
        a_hat_K  = m_K ;
        
        %%% Noise %%%
        n_K = Generate_GM(rho_n, sig_n, x_K, 0, n) ;
        
        %%% Measurements %%%
        y_K = a_K * S_t + n_K;
        
        %%% parameters
        % params_a
        params.roh_a = rho_a;
        params.sig_a = sig_a;
        % params_n
        params.roh_n = rho_n;
        params.sig_n = sig_n;
        
        params.sigmaS = S_t;
        
%% Data plot %%%
%%%%         params_plot
        fs = 17;   % font size
        lw = 3;    % line size
        f_yl = 30;  % font size label
        ms = 10;    % mark size
        % plot a_true
        figure(22)
        subplot(3,1,1);
        set(gca, 'fontsize', fs);
        plot(x_K, a_K,'b', 'linewidth',lw );hold on
        ylabel('$a_K$', 'Interpreter', 'latex', 'fontsize', f_yl)
        % plot mean
        subplot(3,1,1);
        set(gca, 'fontsize', fs);
        plot(x_K, m_K,'--b', 'linewidth', lw);
        ylabel('$m_k$', 'Interpreter', 'latex', 'fontsize', f_yl)
        xlabel('$x$', 'Interpreter', 'latex', 'fontsize', f_yl)
        legend('a','m')
        % plot noise
        subplot(3,1,2);
        set(gca, 'fontsize', fs);
        plot(x_K, n_K,'r', 'linewidth', lw);
        ylabel('$n_k$', 'Interpreter', 'latex', 'fontsize', f_yl)
        xlabel('$x$', 'Interpreter', 'latex', 'fontsize', f_yl)
        % plot y (measurements)
        subplot(3,1,3);
        set(gca, 'fontsize', fs);
        plot(x_K, y_K,'k', 'linewidth', lw);
        ylabel('$y_k$', 'Interpreter', 'latex', 'fontsize', f_yl)
        xlabel('$x$', 'Interpreter', 'latex', 'fontsize', f_yl)
                
        %% Assigning initial points (noisy situation)
        
        indperm = randperm(K);
        
        a_Int = a_K(ind_Int);
        a_hat_Int = a_hat_K(ind_Int);
        m_Int = m_K(ind_Int);
        y_Int = y_K(ind_Int);
        
        % figure(21)
        % subplot(3,1,1);
        % hold on
        % plot (x_Int, a_Int(:),'o', 'MarkerSize',ms, 'LineWidth',lw); hold on;
        
        %% Rest points
        Index = (1:K)';
        ind_rest = setdiff(Index,ind_Int);
        x_rest = x_K (ind_rest);
        m_rest = m_K (ind_rest);
        a_rest = a_K(ind_rest);
        a_hat_rest = a_hat_K(ind_rest);
        y_rest = y_K (ind_rest);
        
        
  %% True SNR(f_hat)
x = x_K;
a = a_K;
ma = m_K;
zK = m_K(ind_Int);
flagPrior = 1
xK   = x_Int;
ma_K = m_Int;
sigmaS = S_t;  % power of s
nbPts = K;

ln = rho_n;
sn2 = sig_n;


% SNR initial
Cn = Cov ( rho_n, sig_n, x_Int, x_Int);
f_hat = pinv(Cn)*m_Int;
SNR_True_Estimated_int = (S_t^2 * f_hat' * a_Int * a_Int' * f_hat ) / (f_hat' * Cn * f_hat)


SNR_True_Estimated = zeros(nbPts,1);
SNR_True_Estimated(ind_Int) = SNR_True_Estimated_int;

for i_rest = 1 : K-K_Int
    x_K_plus_N = [x_Int;x_rest(i_rest)];
    m_k_plus_N = [m_Int;m_rest(i_rest)];
    a_k_plus_N = [a_Int;a_rest(i_rest)];
    Cn = Cov ( rho_n, sig_n, x_K_plus_N, x_K_plus_N);
    f_hat = pinv(Cn)*m_k_plus_N;
    SNR_True_Estimated_rest(i_rest) = (S_t^2 * f_hat' * a_k_plus_N * a_k_plus_N' * f_hat ) / (f_hat' * Cn * f_hat);
end

SNR_True_Estimated(ind_rest) = SNR_True_Estimated_rest;

SNR_True_Estimated_norm_int = (SNR_True_Estimated - SNR_True_Estimated_int) / (max(SNR_True_Estimated) - SNR_True_Estimated_int+tol);

% figure(2)
% plot(SNR_True_Estimated_norm_int)
% grid on
% title('$SNR(\hat{f})$','interpreter','latex')

ind_failure = find(SNR_True_Estimated_norm_int<0);
failure_size(mont_a,mont_m) = size(ind_failure,1);

SNR_True_Estimated_norm_int_mont(mont_a,mont_m,:) = SNR_True_Estimated_norm_int;
SNR_True_Estimated_int_mont(mont_a,mont_m) = SNR_True_Estimated_int;

        
   %%     
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLASSICAL KRIGING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
%         ind_Int_classic = ind_Int;
%         x_Int_classic = x_Int;
%         m_Int_classic = m_Int;
%         a_Int_classic = a_Int;
%         a_hat_Int_classic = a_hat_Int;
%         y_Int_classic = y_Int;
%         
%         x_rest_classic = x_rest;
%         m_rest_classic = m_rest;
%         y_rest_classic = y_rest;
%         a_hat_rest_classic = a_hat_rest;
%         
%         for j = 1 : p
%             j_classic = j
%             mont_m
%             mont_a
%             % select new position
%             J_N = CostFunction_classic(params, x_Int_classic, x_rest_classic);
%             [J_max , Index_max] = max(J_N);  % pay attention, here you have to take the min of the cost_classic
%             x_new_classic = x_rest_classic(Index_max);
%             m_new_classic = m_rest_classic(Index_max);
%             y_new_classic = y_rest_classic(Index_max);
%             %         a_hat_new_classic = a_hat_rest_classic(Index_max);
%             
%             % 'a' estimation at new selected position
%             C_NK_new = Cov ( rho_a, sig_a, x_new_classic, x_Int_classic );
%             C_KK_new = Cov ( rho_a, sig_a, x_Int_classic, x_Int_classic ) ;
%             L_C_KK_new = chol(C_KK_new + tol*eye(size(C_KK_new)), 'lower');
%             a_new_classic = m_new_classic + C_NK_new * (  L_C_KK_new.'  \   ( L_C_KK_new \ ( a_Int_classic - m_Int_classic ) )  );
%             a_hat_new_classic = m_new_classic + C_NK_new * (  L_C_KK_new.'  \   ( L_C_KK_new \ ( a_hat_Int_classic - m_Int_classic ) )  );
%             
%             % Estimate 'S'
%             y_K_plus_N_classic =  [y_Int_classic; y_new_classic] ;
%             a_K_plus_N_classic =  [a_Int_classic; a_new_classic] ;
%             a_hat_K_plus_N_classic =  [a_hat_Int_classic; a_hat_new_classic] ;
%             ind_K_plus_N_classic =  [ind_Int_classic; Index_max] ;
%             
%             x_K_plus_N_classic =  [x_Int_classic; x_new_classic] ;
%             C_K_plus_N = Cov ( rho_n, sig_n, x_K_plus_N_classic, x_K_plus_N_classic );
%             L_C_K_plus_N = chol(C_K_plus_N + tol*eye(size(C_K_plus_N)), 'lower');
%             
%             %         S_hat_classic(j) = 1 / ( a_K_plus_N_classic' * (  L_C_K_plus_N.'  \   ( L_C_K_plus_N \ a_K_plus_N_classic ) ) ) * a_K_plus_N_classic' * ( L_C_K_plus_N.'  \ (L_C_K_plus_N \ y_K_plus_N_classic) );
%             
%             % SNR(f_hat)
%             a_True_K_plus_N_classic = a_K(ind_K_plus_N_classic);
%             f_hat_classic1 = pinv(C_K_plus_N) * a_K_plus_N_classic;
%             f_hat_classic2 = pinv(C_K_plus_N) * a_hat_K_plus_N_classic;
%             SNR_f_hat_classic1(j) = ( S_t^2 * (f_hat_classic1' * a_True_K_plus_N_classic)^2 ) / (f_hat_classic1' * C_K_plus_N * f_hat_classic1 );
%             SNR_f_hat_classic2(j) = ( S_t^2 * (f_hat_classic2' * a_True_K_plus_N_classic)^2 ) / (f_hat_classic2' * C_K_plus_N * f_hat_classic2 );
%             
%             % update data
%             J_classic(j).J = J_N;
%             
%             ind_Int_classic = ind_K_plus_N_classic;
%             x_Int_classic = x_K_plus_N_classic;
%             m_Int_classic = [m_Int_classic; m_new_classic] ;
%             a_Int_classic = a_K_plus_N_classic;
%             a_hat_Int_classic = a_hat_K_plus_N_classic;
%             y_Int_classic = y_K_plus_N_classic;
%             
%             x_rest_classic (Index_max) = [];
%             m_rest_classic (Index_max) = [];
%             y_rest_classic (Index_max) = [];
%             a_hat_rest_classic (Index_max) = [];
%             
%         end
%         
%         
%         %     S_hat_classic_mont(mont,:) = S_hat_classic;
%         SNR_f_hat_classic_mont_m1(mont_m,:) = SNR_f_hat_classic1;
%         SNR_f_hat_classic_mont_m2(mont_m,:) = SNR_f_hat_classic2;
%         
%         J_classic_mont.J(mont_m,:) = J_classic.J;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MI  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         %%
%         ind_Int_MI = ind_Int;
%         x_Int_MI = x_Int;
%         m_Int_MI = m_Int;
%         a_Int_MI = a_Int;
%         a_hat_Int_MI = a_hat_Int;
%         y_Int_MI = y_Int;
%         
%         x_rest_MI = x_rest;
%         m_rest_MI = m_rest;
%         y_rest_MI = y_rest;
%         a_hat_rest_MI = a_hat_rest;
%         
%         for j = 1 : p
%             j_MI = j
%             mont_m
%             mont_a
%             % select new position
%             J_N = CostFunction_MI(params, x_Int_MI, x_rest_MI);
%             [J_max , Index_max] = max(J_N);
%             x_new_MI = x_rest_MI(Index_max);
%             m_new_MI = m_rest_MI(Index_max);
%             y_new_MI = y_rest_MI(Index_max);
%             %         a_hat_new_MI = a_hat_rest_MI(Index_max);
%             
%             % 'a' estimation at new selected position
%             C_NK_new = Cov ( rho_a, sig_a, x_new_MI, x_Int_MI );
%             C_KK_new = Cov ( rho_a, sig_a, x_Int_MI, x_Int_MI ) ;
%             L_C_KK_new = chol(C_KK_new + tol*eye(size(C_KK_new)), 'lower');
%             a_new_MI = m_new_MI + C_NK_new * (  L_C_KK_new.'  \   ( L_C_KK_new \ ( a_Int_MI - m_Int_MI ) )  );
%             a_hat_new_MI = m_new_MI + C_NK_new * (  L_C_KK_new.'  \   ( L_C_KK_new \ ( a_hat_Int_MI - m_Int_MI ) )  );
%             
%             % Estimate 'S'
%             y_K_plus_N_MI =  [y_Int_MI; y_new_MI] ;
%             a_K_plus_N_MI =  [a_Int_MI; a_new_MI] ;
%             a_hat_K_plus_N_MI =  [a_hat_Int_MI; a_hat_new_MI] ;
%             ind_K_plus_N_MI =  [ind_Int_MI; Index_max] ;
%             
%             x_K_plus_N_MI =  [x_Int_MI; x_new_MI] ;
%             C_K_plus_N = Cov ( rho_n, sig_n, x_K_plus_N_MI, x_K_plus_N_MI );
%             L_C_K_plus_N = chol(C_K_plus_N + tol*eye(size(C_K_plus_N)), 'lower');
%             
%             %         S_hat_MI(j) = 1 / ( a_K_plus_N_MI' * (  L_C_K_plus_N.'  \   ( L_C_K_plus_N \ a_K_plus_N_MI ) ) ) * a_K_plus_N_MI' * ( L_C_K_plus_N.'  \ (L_C_K_plus_N \ y_K_plus_N_MI) );
%             
%             % SNR(f_hat)
%             a_True_K_plus_N_MI = a_K(ind_K_plus_N_MI);
%             f_hat_MI1 = pinv(C_K_plus_N) * a_K_plus_N_MI;
%             f_hat_MI2 = pinv(C_K_plus_N) * a_hat_K_plus_N_MI;
%             SNR_f_hat_MI1(j) = ( S_t^2 * (f_hat_MI1' * a_True_K_plus_N_MI)^2 ) / (f_hat_MI1' * C_K_plus_N * f_hat_MI1 );
%             SNR_f_hat_MI2(j) = ( S_t^2 * (f_hat_MI2' * a_True_K_plus_N_MI)^2 ) / (f_hat_MI2' * C_K_plus_N * f_hat_MI2 );
%             
%             % update data
%             J_MI(j).J = J_N;
%             
%             ind_Int_MI = ind_K_plus_N_MI;
%             x_Int_MI = x_K_plus_N_MI;
%             m_Int_MI = [m_Int_MI; m_new_MI] ;
%             a_Int_MI = a_K_plus_N_MI;
%             a_hat_Int_MI = a_hat_K_plus_N_MI;
%             y_Int_MI = y_K_plus_N_MI;
%             
%             x_rest_MI (Index_max) = [];
%             m_rest_MI (Index_max) = [];
%             y_rest_MI (Index_max) = [];
%             a_hat_rest_MI (Index_max) = [];
%             
%         end
%         
%         %     S_hat_MI_mont(mont,:) = S_hat_MI;
%         SNR_f_hat_MI_mont_m1(mont_m,:) = SNR_f_hat_MI1;
%         SNR_f_hat_MI_mont_m2(mont_m,:) = SNR_f_hat_MI2;
%         
%         
%         J_MI_mont.J(mont_m,:) = J_MI.J;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CostFunction_icassp19_mean_SNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%
%         ind_Int_proposed = ind_Int;
%         x_Int_proposed = x_Int;
%         m_Int_proposed = m_Int;
%         a_Int_proposed = a_Int;
%         a_hat_Int_proposed = a_hat_Int;
%         y_Int_proposed = y_Int;
%         
%         x_rest_proposed = x_rest;
%         m_rest_proposed = m_rest;
%         y_rest_proposed = y_rest;
%         a_hat_rest_proposed = a_hat_rest;
%         
%         for j = 1 : p
%             j_proposed = j
%             mont_m
%             mont_a
%             % select new position
%             J_N = CostFunction_icassp19_mean_SNR(params, a_Int_proposed, m_Int_proposed, x_Int_proposed, x_rest_proposed, m_rest_proposed);
%             [J_proposedn , Index_proposed] = max(J_N);
%             x_new_proposed = x_rest_proposed(Index_proposed);
%             m_new_proposed = m_rest_proposed(Index_proposed);
%             y_new_proposed = y_rest_proposed(Index_proposed);
%             %         a_hat_new_proposed = a_hat_rest_proposed(Index_max);
%             
%             % 'a' estimation at new selected position
%             C_NK_new = Cov ( rho_a, sig_a, x_new_proposed, x_Int_proposed );
%             C_KK_new = Cov ( rho_a, sig_a, x_Int_proposed, x_Int_proposed ) ;
%             L_C_KK_new = chol(C_KK_new + tol*eye(size(C_KK_new)), 'lower');
%             a_new_proposed = m_new_proposed + C_NK_new * (  L_C_KK_new.'  \   ( L_C_KK_new \ ( a_Int_proposed - m_Int_proposed ) )  );
%             a_hat_new_proposed = m_new_proposed + C_NK_new * (  L_C_KK_new.'  \   ( L_C_KK_new \ ( a_hat_Int_proposed - m_Int_proposed ) )  );
%             
%             % Estimate 'S'
%             y_K_plus_N_proposed =  [y_Int_proposed; y_new_proposed] ;
%             a_K_plus_N_proposed =  [a_Int_proposed; a_new_proposed] ;
%             a_hat_K_plus_N_proposed =  [a_hat_Int_proposed; a_hat_new_proposed] ;
%             ind_K_plus_N_proposed =  [ind_Int_proposed; Index_proposed] ;
%             
%             x_K_plus_N_proposed =  [x_Int_proposed; x_new_proposed] ;
%             C_K_plus_N = Cov ( rho_n, sig_n, x_K_plus_N_proposed, x_K_plus_N_proposed );
%             L_C_K_plus_N = chol(C_K_plus_N + tol*eye(size(C_K_plus_N)), 'lower');
%             
%             %         S_hat_proposed(j) = 1 / ( a_K_plus_N_proposed' * (  L_C_K_plus_N.'  \   ( L_C_K_plus_N \ a_K_plus_N_proposed ) ) ) * a_K_plus_N_proposed' * ( L_C_K_plus_N.'  \ (L_C_K_plus_N \ y_K_plus_N_proposed) );
%             
%             % SNR(f_hat)
%             a_True_K_plus_N_proposed = a_K(ind_K_plus_N_proposed);
%             f_hat_proposed1 = pinv(C_K_plus_N) * a_K_plus_N_proposed;
%             f_hat_proposed2 = pinv(C_K_plus_N) * a_hat_K_plus_N_proposed;
%             SNR_f_hat_proposed1(j) = ( S_t^2 * (f_hat_proposed1' * a_True_K_plus_N_proposed)^2 ) / (f_hat_proposed1' * C_K_plus_N * f_hat_proposed1 );
%             SNR_f_hat_proposed2(j) = ( S_t^2 * (f_hat_proposed2' * a_True_K_plus_N_proposed)^2 ) / (f_hat_proposed2' * C_K_plus_N * f_hat_proposed2 );
%             
%             % update data
%             J_proposed(j).J = J_N;
%             
%             ind_Int_proposed = ind_K_plus_N_proposed;
%             x_Int_proposed = x_K_plus_N_proposed;
%             m_Int_proposed = [m_Int_proposed; m_new_proposed] ;
%             a_Int_proposed = a_K_plus_N_proposed;
%             a_hat_Int_proposed = a_hat_K_plus_N_proposed;
%             y_Int_proposed = y_K_plus_N_proposed;
%             
%             x_rest_proposed (Index_proposed) = [];
%             m_rest_proposed (Index_proposed) = [];
%             y_rest_proposed (Index_proposed) = [];
%             a_hat_rest_proposed (Index_max) = [];
%             
%         end
%         
%         %     S_hat_proposed_mont(mont,:) = S_hat_proposed;
%         SNR_f_hat_proposed_mont_m1(mont_m,:) = SNR_f_hat_proposed1;
%         SNR_f_hat_proposed_mont_m2(mont_m,:) = SNR_f_hat_proposed2;
%         
%         J_proposed_mont.J(mont_m,:) = J_proposed.J;
%         
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROPOSED_PLUS: CostFunction_CDF_Complementary  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        for jj = 1:jj_max
            
        ind_Int_proposed_plus = ind_Int;
        x_Int_proposed_plus = x_Int;
        m_Int_proposed_plus = m_Int;
        a_Int_proposed_plus = a_Int;
        a_hat_Int_proposed_plus = a_hat_Int;
        y_Int_proposed_plus = y_Int;
        
        ind_rest_proposed_plus = ind_rest;
        x_rest_proposed_plus = x_rest;
        m_rest_proposed_plus = m_rest;
        a_rest_proposed_plus = a_rest;
        y_rest_proposed_plus = y_rest;
        %     a_hat_rest_proposed_plus = a_hat_rest;
        
        for j = 1 : p
            j_proposed_plus = j
            mont_m
            mont_a
            mont_sig_uncertainty
            jj
            
                % select new position
                J_N = CostFunction_CDF_Complementary(ind_Int_proposed_plus, ind_rest_proposed_plus, K, x_K, m_K, Ka,...
                    a_Int_proposed_plus, m_Int_proposed_plus, Kn, a_K,jj,j_step);
                
                
             %%%%%%%%%%%%%%% calculate FPR %%%%%%%%%%%%%%%
             
             J_N_norm = (J_N-min(J_N)) / (max(J_N) - min(J_N));
             index_select_CDF_Complementary = find(J_N_norm>0);
             
             if numel(ind_failure) == 0
                 FPR_CDF_Complementary(mont_m,jj) = 100 * numel(intersect(index_select_CDF_Complementary,ind_failure)) / (numel(ind_failure)+1);
             else
                 FPR_CDF_Complementary(mont_m,jj) = 100 * numel(intersect(index_select_CDF_Complementary,ind_failure)) / (numel(ind_failure));
             end
             
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                [J_proposed_plusn , Index_proposed_plus] = max(J_N);
                x_new_proposed_plus = x_rest_proposed_plus(Index_proposed_plus);
                m_new_proposed_plus = m_rest_proposed_plus(Index_proposed_plus);
                a_new_proposed_plus = a_rest_proposed_plus(Index_proposed_plus);
                y_new_proposed_plus = y_rest_proposed_plus(Index_proposed_plus);
                ind_new_proposed_plus = ind_rest_proposed_plus(Index_proposed_plus);
                %         a_hat_new_proposed_plus = a_hat_rest_proposed_plus(Index_max);
                
                % 'a' estimation at new selected position
                C_NK_new = Cov ( rho_a, sig_a, x_new_proposed_plus, x_Int_proposed_plus );
                C_KK_new = Cov ( rho_a, sig_a, x_Int_proposed_plus, x_Int_proposed_plus ) ;
                L_C_KK_new = chol(C_KK_new + tol*eye(size(C_KK_new)), 'lower');
%                 a_new_proposed_plus = m_new_proposed_plus + C_NK_new * (  L_C_KK_new.'  \   ( L_C_KK_new \ ( a_Int_proposed_plus - m_Int_proposed_plus ) )  );
%                 a_hat_new_proposed_plus = m_new_proposed_plus + C_NK_new * (  L_C_KK_new.'  \   ( L_C_KK_new \ ( a_hat_Int_proposed_plus - m_Int_proposed_plus ) )  );

                a_hat_new_proposed_plus = m_new_proposed_plus;

                % Estimate 'S'
                y_K_plus_N_proposed_plus =  [y_Int_proposed_plus; y_new_proposed_plus] ;
                a_K_plus_N_proposed_plus =  [a_Int_proposed_plus; a_new_proposed_plus] ;
                a_hat_K_plus_N_proposed_plus =  [a_hat_Int_proposed_plus; a_hat_new_proposed_plus] ;
                
                ind_K_plus_N_proposed_plus =  [ind_Int_proposed_plus; ind_new_proposed_plus] ;
                x_K_plus_N_proposed_plus =  [x_Int_proposed_plus; x_new_proposed_plus] ;
                C_K_plus_N = Cov ( rho_n, sig_n, x_K_plus_N_proposed_plus, x_K_plus_N_proposed_plus );
                L_C_K_plus_N = chol(C_K_plus_N + tol*eye(size(C_K_plus_N)), 'lower');
                
                %         S_hat_proposed_plus(j) = 1 / ( a_K_plus_N_proposed_plus' * (  L_C_K_plus_N.'  \   ( L_C_K_plus_N \ a_K_plus_N_proposed_plus ) ) ) * a_K_plus_N_proposed_plus' * ( L_C_K_plus_N.'  \ (L_C_K_plus_N \ y_K_plus_N_proposed_plus) );
                
                % SNR(f_hat)
                a_True_K_plus_N_proposed_plus = a_K(ind_K_plus_N_proposed_plus);
                f_hat_proposed_plus1 = pinv(C_K_plus_N) * a_K_plus_N_proposed_plus;
                f_hat_proposed_plus2 = pinv(C_K_plus_N) * a_hat_K_plus_N_proposed_plus;
                SNR_f_hat_proposed_plus1(j,jj) = ( S_t^2 * (f_hat_proposed_plus1' * a_True_K_plus_N_proposed_plus)^2 ) / (f_hat_proposed_plus1' * C_K_plus_N * f_hat_proposed_plus1 );
                SNR_f_hat_proposed_plus2(j,jj) = ( S_t^2 * (f_hat_proposed_plus2' * a_True_K_plus_N_proposed_plus)^2 ) / (f_hat_proposed_plus2' * C_K_plus_N * f_hat_proposed_plus2 );
                
                
                % update data
                J_proposed_plus(j,jj).J = J_N;
                
                ind_Int_proposed_plus = ind_K_plus_N_proposed_plus;
                x_Int_proposed_plus = x_K_plus_N_proposed_plus;
                m_Int_proposed_plus = [m_Int_proposed_plus; m_new_proposed_plus] ;
                a_Int_proposed_plus = a_K_plus_N_proposed_plus;
                a_hat_Int_proposed_plus = a_hat_K_plus_N_proposed_plus;
                y_Int_proposed_plus = y_K_plus_N_proposed_plus;
                
                ind_rest_proposed_plus (Index_proposed_plus) = [];
                x_rest_proposed_plus (Index_proposed_plus) = [];
                m_rest_proposed_plus (Index_proposed_plus) = [];
                a_rest_proposed_plus (Index_proposed_plus) = [];
                y_rest_proposed_plus (Index_proposed_plus) = [];
                %         a_hat_rest_proposed_plus (Index_max) = [];
                
        end
            
       
        end
        
         %     S_hat_proposed_plus_mont(mont,:) = S_hat_proposed_plus;
        SNR_f_hat_proposed_plus_mont_m1(mont_m,:) = SNR_f_hat_proposed_plus1(end,:);
        SNR_f_hat_proposed_plus_mont_m2(mont_m,:) = SNR_f_hat_proposed_plus2(end,:);
        
        SNR_improve(mont_a,mont_m,:) = SNR_f_hat_proposed_plus_mont_m2(mont_m,:) - SNR_True_Estimated_int;
        
        %         J_proposed_plus_mont.J(mont_m,jj) = J_proposed_plus.J;
        
        if numel(ind_failure) ==0
            mont_m=mont_m-1;
        end
        
    end
    %%
    
    FPR_CDF_Complementary_mont_m(mont_a,:) = mean(FPR_CDF_Complementary,1);
    
%     SNR_f_hat_classic_mean_mont_m1(mont_a,:) = mean(SNR_f_hat_classic_mont_m1,1);
%     SNR_f_hat_MI_mean_mont_m1(mont_a,:) = mean(SNR_f_hat_MI_mont_m1,1);
%     SNR_f_hat_proposed_mean_mont_m1(mont_a,:) = mean(SNR_f_hat_proposed_mont_m1,1);
    SNR_f_hat_proposed_plus_mean_mont_m1(mont_a,:) = mean(SNR_f_hat_proposed_plus_mont_m1,1);
    
%     SNR_f_hat_classic_mean_mont_m2(mont_a,:) = mean(SNR_f_hat_classic_mont_m2,1);
%     SNR_f_hat_MI_mean_mont_m2(mont_a,:) = mean(SNR_f_hat_MI_mont_m2,1);
%     SNR_f_hat_proposed_mean_mont_m2(mont_a,:) = mean(SNR_f_hat_proposed_mont_m2,1);
    SNR_f_hat_proposed_plus_mean_mont_m2(mont_a,:) = mean(SNR_f_hat_proposed_plus_mont_m2,1);
    
end
%%

SNR_improve_mean(mont_sig_uncertainty) = mean(mean(SNR_improve,1),2);

FPR_CDF_Complementary_mont_a(mont_sig_uncertainty) = mean(FPR_CDF_Complementary_mont_m,1);

% SNR_f_hat_classic_mean_mont_a1 = mean(SNR_f_hat_classic_mean_mont_m1,1);
% SNR_f_hat_MI_mean_mont_a1 = mean(SNR_f_hat_MI_mean_mont_m1,1);
% SNR_f_hat_proposed_mean_mont_a1 = mean(SNR_f_hat_proposed_mean_mont_m1,1);
SNR_f_hat_proposed_plus_mean_mont_a1(mont_sig_uncertainty) = mean(SNR_f_hat_proposed_plus_mean_mont_m1,1);

% SNR_f_hat_classic_mean_mont_a2 = mean(SNR_f_hat_classic_mean_mont_m2,1);
% SNR_f_hat_MI_mean_mont_a2 = mean(SNR_f_hat_MI_mean_mont_m2,1);
% SNR_f_hat_proposed_mean_mont_a2 = mean(SNR_f_hat_proposed_mean_mont_m2,1);
SNR_f_hat_proposed_plus_mean_mont_a2(mont_sig_uncertainty) = mean(SNR_f_hat_proposed_plus_mean_mont_m2,1);

SNR_True_Estimated_int_mont_sigma(mont_sig_uncertainty) = mean(mean(SNR_True_Estimated_int_mont))
end
%%
figure(123)
cla;
% plot(SNR_f_hat_proposed_plus_mean_mont_a2,'^','MarkerSize',12), hold on
lw=4;
fs=50;
colorstring = 'kbgry';
set(gca,'FontSize',fs)
set (gcf, 'Position' , [300, 300, 1500, 1500])

SNR_int_mean=mean(mean(SNR_True_Estimated_int_mont))

[hAx,hLine1,hLine2] =  plotyy(s_uncertainty_set,smooth(FPR_CDF_Complementary_mont_a),s_uncertainty_set, smooth(SNR_improve_mean));hold on

% pp1 = plot(s_uncertainty_set,SNR_True_Estimated_int_mont_sigma,'--k')
% legend([pp1], ['Initial SNR=',num2str(SNR_int_mean)],'FontSize',30)

% set([hLine1,hLine2,pp1],'linewidth',3)
set([hLine1,hLine2],'linewidth',lw)
set (hAx, 'FontSize' , fs)

xlabel('$\sigma_u$','interpreter','latex','FontSize',fs*1.5)

ylabel(hAx(1),'$\overline{FPR}~[\%]$','FontSize',fs*1.2,'interpreter','latex') % left y-axis 
ylabel(hAx(2),'$\overline{SNR(\hat{f})-SNR(\hat{f})^{(int)}}$~[dB]','interpreter','latex','FontSize',fs*1.2) % right y-axis
% title(['$\rho_a=$ ',num2str(rho_a),' , $K_{int}=$',num2str(K_Int), '$, P=$',num2str(p),...
%     ' , $K=$',num2str(K), '$, \delta=$',num2str(delta), '$, \sigma_m=$',num2str(flag_bias*s_bias), '$, \sigma_n=$',num2str(sig_n),'$dB$',...
%     ', $\beta=$',num2str(ratio_roh)],'interpreter','latex','FontSize',fs)
% title(['$ \delta=$',num2str(delta)],'interpreter','latex','FontSize',fs)

grid on 

set(hAx(1,1),'YLim',[0,100],'YTick',[0:20:100])
set(hAx(1,2),'YLim',[0,15],'YTick',[0:3:15])


SNR_True_Estimated_int_mont_sigma

%%

figure(1)
plot(a_K)
hold on
plot(n_K,'r')
hold on
plot(m_K,'k-')
grid on
legend('a','n','m')

%%
% 
% figure(2)
% range = 1;  % Grid range
% xx = linspace(0, range, p)';
% lw = 2;
% plot(1:p,SNR_f_hat_classic_mean_mont_a1,'LineWidth', lw); hold on;
% plot(1:p,SNR_f_hat_MI_mean_mont_a1, 'r','LineWidth', lw); hold on;
% plot(1:p,SNR_f_hat_proposed_mean_mont_a1, 'g','LineWidth', lw); hold on;
% plot(1:p,SNR_f_hat_proposed_plus_mean_mont_a1, 'c','LineWidth', lw);
% plot(1:p,SNR_f_hat_proposed_mean_mont_a2, '--g','LineWidth', lw); hold on;
% plot(1:p,SNR_f_hat_proposed_plus_mean_mont_a2, '--c','LineWidth', lw);
% set(gca,'FontSize', 18);
% ylabel('$SNR \hat{f}$','interpreter','latex','FontSize', 35);
% xlabel('Number of sensors', 'FontSize', 35,'interpreter','latex');
% legend({'Entropy', 'MI','$\hat{g}(X_N|a^*_K)$', '$1-G_N(W_M;\delta=0.5|a^*_K)$','$\hat{g}(X_N|\hat{a}_K)$', '$1-G_N(W_M;\delta=0.5|\hat{a}_K)$'},'interpreter','latex','FontSize',30);
% grid on
% 
% figure(3)
% range = 1;  % Grid range
% xx = linspace(0, range, p)';
% lw = 2;
% plot(1:p,SNR_f_hat_classic_mean_mont_a2,'LineWidth', lw); hold on;
% plot(1:p,SNR_f_hat_MI_mean_mont_a2, 'r','LineWidth', lw); hold on;
% plot(1:p,SNR_f_hat_proposed_mean_mont_a2, 'g','LineWidth', lw); hold on;
% plot(1:p,SNR_f_hat_proposed_plus_mean_mont_a2, 'c','LineWidth', lw);
% legend('classic', 'MI', 'proposed', 'proposed-plus');
% % legend('icassp 19', ' 1 - CDF', 'FontSize', 18);
% ylabel('$SNR( \hat{f})$', 'FontSize', 18,'interpreter','latex');
% xlabel('Number of sensors', 'FontSize', 18);
% 
% 
