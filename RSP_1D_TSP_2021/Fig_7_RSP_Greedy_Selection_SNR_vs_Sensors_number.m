
clear all; 
close all; 
clc;

%%%   Generate Grid
K = 200;  % # of possible total sensors
range = 1;  % Grid range
x_K = linspace(0, range, K)';   % places of possible total sensors
tol = 1e-8;
tol_failure = 0;

rho_a = .01;
sig_a = 1;
%%% Source
S_t = 1;

%%% Noise
% SNR_dB = 1.2;
% SNR = (10^SNR_dB) / 10;
% sig_n = sig_a * S_t / sqrt(SNR);
sig_n = 1;
ratio_roh = .5;
rho_n = ratio_roh * rho_a;

Kn = Cov ( rho_n, sig_n, x_K, x_K);
Ka = Cov ( rho_a, sig_a, x_K, x_K);

% ''a''
n = 1;  % # of mixture in Mixtured Gaussian Process

%%%% uncertainty
l_uncertainty = rho_a;
% SNR_dB_uncertainty = 1;
% SNR_uncertainty = (10^SNR_dB_uncertainty) / 10;
% s_uncertainty = sig_a / sqrt(SNR_uncertainty);
s_uncertainty = .1;
flag_uncertainty = 1;

%%%% Bias
l_bias = rho_a;
% SNR_dB_bias = 5;
% SNR_bias = (10^SNR_dB_bias) / 10;
% s_bias = sig_a / sqrt(SNR_bias);
s_bias = 0;
flag_bias = 0;

% Assigning initial points (noisy situation)
ind_Int = [100];
K_Int = size(ind_Int,1);
x_Int = x_K(ind_Int);

p = 10;  % number of sensors that we want to select greedy

% setting delta parameter
jj_max=1;
j_step=13;
delta = j_step * [1:jj_max]

%%% parameters
% params_a
params.rho_a = rho_a;
params.sig_a = sig_a;
% params_n
params.rho_n = rho_n;
params.sig_n = sig_n;
params.sigmaS = S_t;

mont_a_max = 10;
mont_m_max = 10;

for mont_a = 1:mont_a_max
    
    
    a_K = Generate_GM(rho_a, sig_a, x_K, 0, n) ;
    
    
    for mont_m = 1 : mont_m_max
        
        
        
        %%% Generate ''a'' %%%
        
        % Mean
        %     ma1 = Generate_GM(roh_a, sig_a, x_K, 0, 1) ;
        %     ma1 = 0;
        %     % ''a''
        %     n = 1;  % # of mixture in Mixtured Gaussian Process
        %     l_bias = roh_a;
        %     s_bias = 0.01;
        %     bias = Generate_GM(l_bias, s_bias, x_K, 0, n);
        %     m_K = a_K + bias;
        %
        
        uncertainty_a = Generate_GM(l_uncertainty, s_uncertainty, x_K, 0, n);
        bias = Generate_GM(l_bias, s_bias, x_K, 0, n);
        m_K = a_K + (flag_uncertainty * uncertainty_a) + (flag_bias * bias);
        
        
        %     la = roh_a;
        %     sa2 = 3*sig_a;      % uncertainty on spatial gain 'a'
        a_hat_K  = m_K;
        
        %%% Noise %%%
        
        n_K = Generate_GM(rho_n, sig_n, x_K, 0, n) ;
        
        %%% Measurements %%%
        y_K = a_K * S_t + n_K;
        
        
        
        %% Data plot %%%
%%%%         params_plot
%         fs = 17;   % font size
%         lw = 3;    % line size
%         f_yl = 30;  % font size label
%         ms = 10;    % mark size
%         % plot a_true
%         figure(22)
%         subplot(3,1,1);
%         set(gca, 'fontsize', fs);
%         plot(x_K, a_K,'b', 'linewidth',lw );hold on
%         ylabel('$a_K$', 'Interpreter', 'latex', 'fontsize', f_yl)
%         % plot mean
%         subplot(3,1,1);
%         set(gca, 'fontsize', fs);
%         plot(x_K, m_K,'--b', 'linewidth', lw);
%         ylabel('$m_k$', 'Interpreter', 'latex', 'fontsize', f_yl)
%         xlabel('$x$', 'Interpreter', 'latex', 'fontsize', f_yl)
%         legend('a','m')
%         % plot noise
%         subplot(3,1,2);
%         set(gca, 'fontsize', fs);
%         plot(x_K, n_K,'r', 'linewidth', lw);
%         ylabel('$n_k$', 'Interpreter', 'latex', 'fontsize', f_yl)
%         xlabel('$x$', 'Interpreter', 'latex', 'fontsize', f_yl)
%         % plot y (measurements)
%         subplot(3,1,3);
%         set(gca, 'fontsize', fs);
%         plot(x_K, y_K,'k', 'linewidth', lw);
%         ylabel('$y_k$', 'Interpreter', 'latex', 'fontsize', f_yl)
%         xlabel('$x$', 'Interpreter', 'latex', 'fontsize', f_yl)
        
%%        
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
        a_rest = a_K (ind_rest);
        a_hat_rest = a_hat_K(ind_rest);
        y_rest = y_K (ind_rest);
        %     N = (K - K_Int);
        
        
        
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
        ind_Int_classic = ind_Int;
        x_Int_classic = x_Int;
        m_Int_classic = m_Int;
        a_Int_classic = a_Int;
        a_hat_Int_classic = a_hat_Int;
        y_Int_classic = y_Int;
        
        x_rest_classic = x_rest;
        m_rest_classic = m_rest;
        a_rest_classic = a_rest;
        y_rest_classic = y_rest;
        a_hat_rest_classic = a_hat_rest;
        
        for j_p = 1 : p
            j_classic = j_p
            mont_m
            mont_a
            % select new position
            J_N = CostFunction_classic(params, x_Int_classic, x_rest_classic);
            [J_max , Index_max] = max(J_N);  % pay attention, here you have to take the min of the cost_classic
            x_new_classic = x_rest_classic(Index_max);
            m_new_classic = m_rest_classic(Index_max);
            a_new_classic = a_rest_classic(Index_max);
            y_new_classic = y_rest_classic(Index_max);
            %         a_hat_new_classic = a_hat_rest_classic(Index_max);
            
            % 'a' estimation at new selected position
            C_NK_new = Cov ( rho_a, sig_a, x_new_classic, x_Int_classic );
            C_KK_new = Cov ( rho_a, sig_a, x_Int_classic, x_Int_classic ) ;
            L_C_KK_new = chol(C_KK_new + 1e-10*eye(size(C_KK_new)), 'lower');
            %         a_new_classic = m_new_classic + C_NK_new * (  L_C_KK_new.'  \   ( L_C_KK_new \ ( a_Int_classic - m_Int_classic ) )  );
            %         a_hat_new_classic = m_new_classic + C_NK_new * (  L_C_KK_new.'  \   ( L_C_KK_new \ ( a_hat_Int_classic - m_Int_classic ) )  );
            
            %         a_new_classic = m_new_classic;
            a_hat_new_classic = m_new_classic;
            
            % Estimate 'S'
            y_K_plus_N_classic =  [y_Int_classic; y_new_classic] ;
            a_K_plus_N_classic =  [a_Int_classic; a_new_classic] ;
            a_hat_K_plus_N_classic =  [a_hat_Int_classic; a_hat_new_classic] ;
            ind_K_plus_N_classic =  [ind_Int_classic; Index_max] ;
            
            x_K_plus_N_classic =  [x_Int_classic; x_new_classic] ;
            C_K_plus_N = Cov ( rho_n, sig_n, x_K_plus_N_classic, x_K_plus_N_classic );
            L_C_K_plus_N = chol(C_K_plus_N + 1e-10*eye(size(C_K_plus_N)), 'lower');
            
            %         S_hat_classic(j) = 1 / ( a_K_plus_N_classic' * (  L_C_K_plus_N.'  \   ( L_C_K_plus_N \ a_K_plus_N_classic ) ) ) * a_K_plus_N_classic' * ( L_C_K_plus_N.'  \ (L_C_K_plus_N \ y_K_plus_N_classic) );
            
            % SNR(f_hat)
            a_True_K_plus_N_classic = a_K(ind_K_plus_N_classic);
            f_hat_classic1 = pinv(C_K_plus_N) * a_K_plus_N_classic;
            f_hat_classic2 = pinv(C_K_plus_N) * a_hat_K_plus_N_classic;
            SNR_f_hat_classic1(j_p) = ( S_t^2 * (f_hat_classic1' * a_True_K_plus_N_classic)^2 ) / (f_hat_classic1' * C_K_plus_N * f_hat_classic1 );
            SNR_f_hat_classic2(j_p) = ( S_t^2 * (f_hat_classic2' * a_True_K_plus_N_classic)^2 ) / (f_hat_classic2' * C_K_plus_N * f_hat_classic2 );
            
            % update data
            J_classic(j_p).J = J_N;
            
            ind_Int_classic = ind_K_plus_N_classic;
            x_Int_classic = x_K_plus_N_classic;
            m_Int_classic = [m_Int_classic; m_new_classic] ;
            a_Int_classic = a_K_plus_N_classic;
            a_hat_Int_classic = a_hat_K_plus_N_classic;
            y_Int_classic = y_K_plus_N_classic;
            
            x_rest_classic (Index_max) = [];
            m_rest_classic (Index_max) = [];
            y_rest_classic (Index_max) = [];
            a_hat_rest_classic (Index_max) = [];
            a_rest_classic (Index_max) = [];
            
        %%%% FPR %%%%
            J_N_norm = (J_N-min(J_N) ) /( max(J_N) - min(J_N) );
            index_select_classic = find(J_N_norm>tol_failure);
            if numel(ind_failure) == 0
                FPR_CDF_classic(mont_m,j_p) = 0;
            else
                FPR_CDF_classic(mont_m,j_p) = 100 * numel(intersect(index_select_classic,ind_failure)) / (numel(ind_failure));
            end
        %%%%%%%%%%%%%            
            
        end
        
        %     S_hat_classic_mont(mont,:) = S_hat_classic;
        SNR_f_hat_classic_mont_m1(mont_m,:) = SNR_f_hat_classic1;
        SNR_f_hat_classic_mont_m2(mont_m,:) = SNR_f_hat_classic2;
        
        J_classic_mont.J(mont_m,:) = J_classic.J;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MI  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%
        ind_Int_MI = ind_Int;
        x_Int_MI = x_Int;
        m_Int_MI = m_Int;
        a_Int_MI = a_Int;
        a_hat_Int_MI = a_hat_Int;
        y_Int_MI = y_Int;
        
        x_rest_MI = x_rest;
        m_rest_MI = m_rest;
        y_rest_MI = y_rest;
        a_rest_MI = a_rest;
        a_hat_rest_MI = a_hat_rest;
        
        for j_p = 1 : p
            j_MI = j_p
            mont_m
            mont_a
            % select new position
            J_N = CostFunction_MI(params, x_Int_MI, x_rest_MI);
            [J_max , Index_max] = max(J_N);
            x_new_MI = x_rest_MI(Index_max);
            m_new_MI = m_rest_MI(Index_max);
            a_new_MI = a_rest_MI(Index_max);
            y_new_MI = y_rest_MI(Index_max);
            %         a_hat_new_MI = a_hat_rest_MI(Index_max);
            
            % 'a' estimation at new selected position
            C_NK_new = Cov ( rho_a, sig_a, x_new_MI, x_Int_MI );
            C_KK_new = Cov ( rho_a, sig_a, x_Int_MI, x_Int_MI ) ;
            L_C_KK_new = chol(C_KK_new + 1e-10*eye(size(C_KK_new)), 'lower');
            %         a_new_MI = m_new_MI + C_NK_new * (  L_C_KK_new.'  \   ( L_C_KK_new \ ( a_Int_MI - m_Int_MI ) )  );
            %         a_hat_new_MI = m_new_MI + C_NK_new * (  L_C_KK_new.'  \   ( L_C_KK_new \ ( a_hat_Int_MI - m_Int_MI ) )  );
            
            %         a_new_MI = a_new_MI;
            a_hat_new_MI = m_new_MI;
            
            % Estimate 'S'
            y_K_plus_N_MI =  [y_Int_MI; y_new_MI] ;
            a_K_plus_N_MI =  [a_Int_MI; a_new_MI] ;
            a_hat_K_plus_N_MI =  [a_hat_Int_MI; a_hat_new_MI] ;
            ind_K_plus_N_MI =  [ind_Int_MI; Index_max] ;
            
            x_K_plus_N_MI =  [x_Int_MI; x_new_MI] ;
            C_K_plus_N = Cov ( rho_n, sig_n, x_K_plus_N_MI, x_K_plus_N_MI );
            L_C_K_plus_N = chol(C_K_plus_N + 1e-10*eye(size(C_K_plus_N)), 'lower');
            
            %         S_hat_MI(j) = 1 / ( a_K_plus_N_MI' * (  L_C_K_plus_N.'  \   ( L_C_K_plus_N \ a_K_plus_N_MI ) ) ) * a_K_plus_N_MI' * ( L_C_K_plus_N.'  \ (L_C_K_plus_N \ y_K_plus_N_MI) );
            
            % SNR(f_hat)
            a_True_K_plus_N_MI = a_K(ind_K_plus_N_MI);
            f_hat_MI1 = pinv(C_K_plus_N) * a_K_plus_N_MI;
            f_hat_MI2 = pinv(C_K_plus_N) * a_hat_K_plus_N_MI;
            SNR_f_hat_MI1(j_p) = ( S_t^2 * (f_hat_MI1' * a_True_K_plus_N_MI)^2 ) / (f_hat_MI1' * C_K_plus_N * f_hat_MI1 );
            SNR_f_hat_MI2(j_p) = ( S_t^2 * (f_hat_MI2' * a_True_K_plus_N_MI)^2 ) / (f_hat_MI2' * C_K_plus_N * f_hat_MI2 );
            
            % update data
            J_MI(j_p).J = J_N;
            
            ind_Int_MI = ind_K_plus_N_MI;
            x_Int_MI = x_K_plus_N_MI;
            m_Int_MI = [m_Int_MI; m_new_MI] ;
            a_Int_MI = a_K_plus_N_MI;
            a_hat_Int_MI = a_hat_K_plus_N_MI;
            y_Int_MI = y_K_plus_N_MI;
            
            x_rest_MI (Index_max) = [];
            m_rest_MI (Index_max) = [];
            y_rest_MI (Index_max) = [];
            a_hat_rest_MI (Index_max) = [];
            a_rest_MI (Index_max) = [];
            
            
        %%%% FPR %%%%
            J_N_norm = (J_N-min(J_N) ) /( max(J_N) - min(J_N) );
            index_select_MI = find(J_N_norm>tol_failure);
            if numel(ind_failure) == 0
                FPR_CDF_MI(mont_m,j_p) = 0;
            else
                FPR_CDF_MI(mont_m,j_p) = 100 * numel(intersect(index_select_MI,ind_failure)) / (numel(ind_failure));
            end
        %%%%%%%%%%%%%            
            
        end
        
        %     S_hat_MI_mont(mont,:) = S_hat_MI;
        SNR_f_hat_MI_mont_m1(mont_m,:) = SNR_f_hat_MI1;
        SNR_f_hat_MI_mont_m2(mont_m,:) = SNR_f_hat_MI2;
        
        
        J_MI_mont.J(mont_m,:) = J_MI.J;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CostFunction_icassp19_mean_SNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        ind_Int_proposed = ind_Int;
        x_Int_proposed = x_Int;
        m_Int_proposed = m_Int;
        a_Int_proposed = a_Int;
        a_hat_Int_proposed = a_hat_Int;
        y_Int_proposed = y_Int;
        
        x_rest_proposed = x_rest;
        m_rest_proposed = m_rest;
        y_rest_proposed = y_rest;
        a_hat_rest_proposed = a_hat_rest;
        a_rest_proposed = a_rest;
        
        
        for j_p = 1 : p
            j_proposed = j_p
            mont_m
            mont_a
            % select new position
            J_N = CostFunction_icassp19_mean_SNR(params, a_Int_proposed, m_Int_proposed, x_Int_proposed, x_rest_proposed, m_rest_proposed);
            [J_proposedn , Index_proposed] = max(J_N);
            x_new_proposed = x_rest_proposed(Index_proposed);
            m_new_proposed = m_rest_proposed(Index_proposed);
            a_new_proposed = a_rest_proposed(Index_proposed);
            y_new_proposed = y_rest_proposed(Index_proposed);
            %         a_hat_new_proposed = a_hat_rest_proposed(Index_max);
            
            % 'a' estimation at new selected position
            C_NK_new = Cov ( rho_a, sig_a, x_new_proposed, x_Int_proposed );
            C_KK_new = Cov ( rho_a, sig_a, x_Int_proposed, x_Int_proposed ) ;
            L_C_KK_new = chol(C_KK_new + 1e-10*eye(size(C_KK_new)), 'lower');
            %         a_new_proposed = m_new_proposed + C_NK_new * (  L_C_KK_new.'  \   ( L_C_KK_new \ ( a_Int_proposed - m_Int_proposed ) )  );
            %         a_hat_new_proposed = m_new_proposed + C_NK_new * (  L_C_KK_new.'  \   ( L_C_KK_new \ ( a_hat_Int_proposed - m_Int_proposed ) )  );
            
            a_hat_new_proposed = m_new_proposed;
            
            % Estimate 'S'
            y_K_plus_N_proposed =  [y_Int_proposed; y_new_proposed] ;
            a_K_plus_N_proposed =  [a_Int_proposed; a_new_proposed] ;
            a_hat_K_plus_N_proposed =  [a_hat_Int_proposed; a_hat_new_proposed] ;
            ind_K_plus_N_proposed =  [ind_Int_proposed; Index_proposed] ;
            
            x_K_plus_N_proposed =  [x_Int_proposed; x_new_proposed] ;
            C_K_plus_N = Cov ( rho_n, sig_n, x_K_plus_N_proposed, x_K_plus_N_proposed );
            L_C_K_plus_N = chol(C_K_plus_N + 1e-10*eye(size(C_K_plus_N)), 'lower');
            
            %         S_hat_proposed(j) = 1 / ( a_K_plus_N_proposed' * (  L_C_K_plus_N.'  \   ( L_C_K_plus_N \ a_K_plus_N_proposed ) ) ) * a_K_plus_N_proposed' * ( L_C_K_plus_N.'  \ (L_C_K_plus_N \ y_K_plus_N_proposed) );
            
            % SNR(f_hat)
            a_True_K_plus_N_proposed = a_K(ind_K_plus_N_proposed);
            f_hat_proposed1 = pinv(C_K_plus_N) * a_K_plus_N_proposed;
            f_hat_proposed2 = pinv(C_K_plus_N) * a_hat_K_plus_N_proposed;
            SNR_f_hat_proposed1(j_p) = ( S_t^2 * (f_hat_proposed1' * a_True_K_plus_N_proposed)^2 ) / (f_hat_proposed1' * C_K_plus_N * f_hat_proposed1 );
            SNR_f_hat_proposed2(j_p) = ( S_t^2 * (f_hat_proposed2' * a_True_K_plus_N_proposed)^2 ) / (f_hat_proposed2' * C_K_plus_N * f_hat_proposed2 );
            
            % update data
            J_proposed(j_p).J = J_N;
            
            ind_Int_proposed = ind_K_plus_N_proposed;
            x_Int_proposed = x_K_plus_N_proposed;
            m_Int_proposed = [m_Int_proposed; m_new_proposed] ;
            a_Int_proposed = a_K_plus_N_proposed;
            a_hat_Int_proposed = a_hat_K_plus_N_proposed;
            y_Int_proposed = y_K_plus_N_proposed;
            
            x_rest_proposed (Index_proposed) = [];
            m_rest_proposed (Index_proposed) = [];
            y_rest_proposed (Index_proposed) = [];
            a_hat_rest_proposed (Index_max) = [];
            a_rest_proposed (Index_max) = [];
            
            %%%%% FPR %%%%%
            J_N_norm = (J_N-min(J_N) ) /( max(J_N) - min(J_N) );
            index_select_proposed = find(J_N_norm>tol_failure);
            if numel(ind_failure) == 0
                FPR_proposed(mont_m,j_p) = 0;
            else
                FPR_proposed(mont_m,j_p) = 100 * numel(intersect(index_select_proposed,ind_failure)) / (numel(ind_failure));
            end
            %%%%%%%%%%%%%%%%%
                        
        end
        
        %     S_hat_proposed_mont(mont,:) = S_hat_proposed;
        SNR_f_hat_proposed_mont_m1(mont_m,:) = SNR_f_hat_proposed1;
        SNR_f_hat_proposed_mont_m2(mont_m,:) = SNR_f_hat_proposed2;
        
        J_proposed_mont.J(mont_m,:) = J_proposed.J;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROPOSED_PLUS: CostFunction_CDF_Complementary  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
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
        
        for j_p = 1 : p
            j_proposed_plus = j_p
            mont_m
            mont_a
            
            for jj = 1 : jj_max
                % select new position
                J_N = CostFunction_CDF_Complementary(ind_Int_proposed_plus, ind_rest_proposed_plus, K, x_K, m_K, Ka,...
                    a_Int_proposed_plus, m_Int_proposed_plus, Kn, a_K,jj,j_step);
                
            end
            
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
            L_C_KK_new = chol(C_KK_new + 1e-10*eye(size(C_KK_new)), 'lower');
            %         a_new_proposed_plus = m_new_proposed_plus + C_NK_new * (  L_C_KK_new.'  \   ( L_C_KK_new \ ( a_Int_proposed_plus - m_Int_proposed_plus ) )  );
            %         a_hat_new_proposed_plus = m_new_proposed_plus + C_NK_new * (  L_C_KK_new.'  \   ( L_C_KK_new \ ( a_hat_Int_proposed_plus - m_Int_proposed_plus ) )  );
            
            a_hat_new_proposed_plus = m_new_proposed_plus;
            
            % Estimate 'S'
            y_K_plus_N_proposed_plus =  [y_Int_proposed_plus; y_new_proposed_plus] ;
            a_K_plus_N_proposed_plus =  [a_Int_proposed_plus; a_new_proposed_plus] ;
            a_hat_K_plus_N_proposed_plus =  [a_hat_Int_proposed_plus; a_hat_new_proposed_plus] ;
            
            ind_K_plus_N_proposed_plus =  [ind_Int_proposed_plus; ind_new_proposed_plus] ;
            x_K_plus_N_proposed_plus =  [x_Int_proposed_plus; x_new_proposed_plus] ;
            C_K_plus_N = Cov ( rho_n, sig_n, x_K_plus_N_proposed_plus, x_K_plus_N_proposed_plus );
            L_C_K_plus_N = chol(C_K_plus_N + 1e-10*eye(size(C_K_plus_N)), 'lower');
            
            %         S_hat_proposed_plus(j) = 1 / ( a_K_plus_N_proposed_plus' * (  L_C_K_plus_N.'  \   ( L_C_K_plus_N \ a_K_plus_N_proposed_plus ) ) ) * a_K_plus_N_proposed_plus' * ( L_C_K_plus_N.'  \ (L_C_K_plus_N \ y_K_plus_N_proposed_plus) );
            
            % SNR(f_hat)
            a_True_K_plus_N_proposed_plus = a_K(ind_K_plus_N_proposed_plus);
            f_hat_proposed_plus1 = pinv(C_K_plus_N) * a_K_plus_N_proposed_plus;
            f_hat_proposed_plus2 = pinv(C_K_plus_N) * a_hat_K_plus_N_proposed_plus;
            SNR_f_hat_proposed_plus1(j_p) = ( S_t^2 * (f_hat_proposed_plus1' * a_True_K_plus_N_proposed_plus)^2 ) / (f_hat_proposed_plus1' * C_K_plus_N * f_hat_proposed_plus1 );
            SNR_f_hat_proposed_plus2(j_p) = ( S_t^2 * (f_hat_proposed_plus2' * a_True_K_plus_N_proposed_plus)^2 ) / (f_hat_proposed_plus2' * C_K_plus_N * f_hat_proposed_plus2 );
            
            
            % update data
            J_proposed_plus(j_p).J = J_N;
            
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
            
            
            %%%%% FPR %%%%%
            J_N_norm = (J_N-min(J_N) ) /( max(J_N) - min(J_N) );
            index_select_proposed_plus = find(J_N_norm>tol_failure);
            if numel(ind_failure) == 0
                FPR_CDF_proposed_plus(mont_m,j_p) = 0;
            else
                FPR_CDF_proposed_plus(mont_m,j_p) = 100 * numel(intersect(index_select_proposed_plus,ind_failure)) / (numel(ind_failure));
            end
            %%%%%%%%%%%%%%%%%
            
        end
        
        %     S_hat_proposed_plus_mont(mont,:) = S_hat_proposed_plus;
        SNR_f_hat_proposed_plus_mont_m1(mont_m,:) = SNR_f_hat_proposed_plus1;
        SNR_f_hat_proposed_plus_mont_m2(mont_m,:) = SNR_f_hat_proposed_plus2;
        
        J_proposed_plus_mont.J(mont_m,:) = J_proposed_plus.J;
        
    end
    %%
    %%%%% FPR average on mont_m
    FPR_CDF_classic_mean_mont_m(mont_a,:) = mean(FPR_CDF_classic,1);
    FPR_CDF_MI_mean_mont_m(mont_a,:) = mean(FPR_CDF_MI,1);
    FPR_CDF_proposed_mean_mont_m(mont_a,:) = mean(FPR_proposed,1);
    FPR_CDF_proposed_plus_mean_mont_m(mont_a,:) = mean(FPR_CDF_proposed_plus,1);
    
    
    
    SNR_f_hat_classic_mean_mont_m1(mont_a,:) = mean(SNR_f_hat_classic_mont_m1,1);
    SNR_f_hat_MI_mean_mont_m1(mont_a,:) = mean(SNR_f_hat_MI_mont_m1,1);
    SNR_f_hat_proposed_mean_mont_m1(mont_a,:) = mean(SNR_f_hat_proposed_mont_m1,1);
    SNR_f_hat_proposed_plus_mean_mont_m1(mont_a,:) = mean(SNR_f_hat_proposed_plus_mont_m1,1);
    
    SNR_f_hat_classic_mean_mont_m2(mont_a,:) = mean(SNR_f_hat_classic_mont_m2,1);
    SNR_f_hat_MI_mean_mont_m2(mont_a,:) = mean(SNR_f_hat_MI_mont_m2,1);
    SNR_f_hat_proposed_mean_mont_m2(mont_a,:) = mean(SNR_f_hat_proposed_mont_m2,1);
    SNR_f_hat_proposed_plus_mean_mont_m2(mont_a,:) = mean(SNR_f_hat_proposed_plus_mont_m2,1);
    
end

%%
%%%%% FPR average on mont_a
FPR_CDF_classic_mean_mont_a = mean(FPR_CDF_classic_mean_mont_m,1);
FPR_CDF_MI_mean_mont_a = mean(FPR_CDF_MI_mean_mont_m,1);
FPR_CDF_proposed_mean_mont_a = mean(FPR_CDF_proposed_mean_mont_m,1);
FPR_CDF_proposed_plus_mean_mont_a = mean(FPR_CDF_proposed_plus_mean_mont_m,1);

%%%%% SNR(f_hat) average on mont_a
SNR_f_hat_classic_mean_mont_a1 = mean(SNR_f_hat_classic_mean_mont_m1,1);
SNR_f_hat_MI_mean_mont_a1 = mean(SNR_f_hat_MI_mean_mont_m1,1);
SNR_f_hat_proposed_mean_mont_a1 = mean(SNR_f_hat_proposed_mean_mont_m1,1);
SNR_f_hat_proposed_plus_mean_mont_a1 = mean(SNR_f_hat_proposed_plus_mean_mont_m1,1);

SNR_f_hat_classic_mean_mont_a2 = mean(SNR_f_hat_classic_mean_mont_m2,1);
SNR_f_hat_MI_mean_mont_a2 = mean(SNR_f_hat_MI_mean_mont_m2,1);
SNR_f_hat_proposed_mean_mont_a2 = mean(SNR_f_hat_proposed_mean_mont_m2,1);
SNR_f_hat_proposed_plus_mean_mont_a2 = mean(SNR_f_hat_proposed_plus_mean_mont_m2,1);


%%

figure(2)
range = 1;  % Grid range
xx = linspace(0, range, p)';
lw = 4;
fs = 55;
plot(1:p,SNR_f_hat_classic_mean_mont_a1,'LineWidth', lw,'color',[0.4660 0.6740 0.1880]); hold on;
plot(1:p,SNR_f_hat_MI_mean_mont_a1,'LineWidth', lw,'Color',[0.8500 0.3250 0.0980]); hold on;
% plot(1:p,SNR_f_hat_proposed_mean_mont_a1, 'g','LineWidth', lw); hold on;
% plot(1:p,SNR_f_hat_proposed_plus_mean_mont_a1, 'c','LineWidth', lw);
plot(1:p,SNR_f_hat_proposed_mean_mont_a2,'LineWidth', lw,'Color',[0.6350 0.0780 0.1840]); hold on;
plot(1:p,SNR_f_hat_proposed_plus_mean_mont_a2,'LineWidth', lw,'Color',[0 0.4470 0.7410]);
% plot(1:p,SNR_f_hat_proposed_plus_mean_mont_a2,'--','LineWidth', lw,'Color',[0 0.4470 0.7410]);
set(gca,'FontSize', fs);
ylabel('$\overline{SNR (\hat{f})}~$[dB]','interpreter','latex','FontSize', fs);
xlabel('Number of sensors', 'FontSize', fs,'interpreter','latex');
% title(['$\sigma_{uncertainty}=$ ',num2str(s_uncertainty)],'interpreter','latex','FontSize', fs)
% legend({'Entropy', 'MI','$\hat{g}(X_N|a^*_K)$', '$1-G_N(W_M;\delta=0.5|a^*_K)$','$\hat{g}(X_N|\hat{a}_K)$', '$1-G_N(W_M;\delta=0.5|\hat{a}_K)$'},'interpreter','latex','FontSize',30);
legend({'$J_{H}(X_N|\hat{a}_K)$', '$J_{MI}(X_N|\hat{a}_K)$','$J_{E}(X_N|\hat{a}_K)$', ['$J_{P}(X_N,\delta=$', num2str(j_step) , '$|\hat{a}_K)$']...
%     ,'$J_{P}(X_N,\delta=5|\hat{a}_K)$','$J_{P}(X_N,\delta=0|\hat{a}_K)$'...,'$J_{P}(X_N,\delta=1.5|\hat{a}_K)$'...
    },'interpreter','latex','FontSize',fs, 'Location','best');
ylim([0,40])
xticks ([0:2:10])
grid on
set(gca,'FontSize',fs)
set (gcf, 'Position' , [50, 150, 1200, 1000])

%%%%%%% plot FPR %%%%%%
figure(3)
range = 1;  % Grid range
xx = linspace(0, range, p)';
lw = 4;
plot(1:p,smooth(FPR_CDF_classic_mean_mont_a),'LineWidth', lw,'Color',[0.4660 0.6740 0.1880]); hold on;
plot(1:p,smooth(FPR_CDF_MI_mean_mont_a), 'LineWidth', lw,'Color',[0.8500 0.3250 0.0980]); hold on;
plot(1:p,smooth(FPR_CDF_proposed_mean_mont_a),'LineWidth', lw,'Color',[0.6350 0.0780 0.1840]); hold on;
plot(1:p,smooth(FPR_CDF_proposed_plus_mean_mont_a),'LineWidth', lw,'Color',[0 0.4470 0.7410]);
% plot(1:p,smooth(FPR_CDF_proposed_plus_mean_mont_a),'--','LineWidth', lw,'Color',[0 0.4470 0.7410]);
legend({'$J_{H}(X_N|\hat{a}_K)$', '$J_{MI}(X_N|\hat{a}_K)$','$J_{E}(X_N|\hat{a}_K)$', ['$J_{P}(X_N,\delta=$', num2str(j_step) , '$|\hat{a}_K)$']...
%     ,'$J_{P}(X_N,\delta=5|\hat{a}_K)$','$J_{P}(X_N,\delta=0|\hat{a}_K)$','$J_{P}(X_N,\delta=1.5|\hat{a}_K)$'...
    },'interpreter','latex','FontSize',fs, 'Location','best');
ylabel('$\overline{FPR}~[\%]$','interpreter','latex','FontSize', fs);
xlabel('Number of sensors', 'FontSize', fs,'interpreter','latex');
% title(['$\sigma_{uncertainty}=$ ',num2str(s_uncertainty)],'interpreter','latex','FontSize', fs)
ylim([0,100])
xticks ([0:2:10])
grid on
set(gca,'FontSize',fs)
set (gcf, 'Position' , [1200, 150, 1200, 1000])


