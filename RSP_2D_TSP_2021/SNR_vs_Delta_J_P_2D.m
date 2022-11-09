%%%%%%%%%%%%%%%%%%% 2D %%%%%%%%%%%%%%%%
%%%%%%%%%%% Proposed Criterion based on the pdf of the SNR J_P
%%%%%%%%%%% Here, we will see that there is a single optimum value for the
%%%%%%%%%%% parameter Delta in a 2D setting
%%%%%%%%%%% We plot average SNR v.s. different Delta values

clear all;
close all;
clc;

%%%   Generate Grid
K = 20;  % # of possible total sensors
Points_Number = K^2;
range = 1;  % Grid range
xx_K = linspace(0, range, K)';   % places of possible total sensors
[x_K1,x_K2] = meshgrid(xx_K,xx_K);

x_K = [x_K1(:),x_K2(:)];
figure(10)
plot(x_K(:,1),x_K(:,2),'o')
title('2-D Grid')

tol = 1e-10;
tol_failure = 1e-8;

rho_a = .2;
sig_a = 1;
%%% Source
S_t = 2;

%%% Noise
% SNR_dB = 1.2;
% SNR = (10^SNR_dB) / 10;
% sig_n = sig_a * S_t / sqrt(SNR);
sig_n = 1.5;
ratio_rho = .5;
rho_n = ratio_rho * rho_a;

Kn = My_2D_Cov ( rho_n, sig_n, x_K, x_K);
Ka = My_2D_Cov ( rho_a, sig_a, x_K, x_K);

% ''a''
n = 1;  % # of mixture in Mixtured Gaussian Process

%%%% uncertainty
l_uncertainty = rho_a/3;
% SNR_dB_uncertainty = 1;
% SNR_uncertainty = (10^SNR_dB_uncertainty) / 10;
% s_uncertainty = sig_a / sqrt(SNR_uncertainty);
s_uncertainty = .8;
flag_uncertainty = 1;

%%%% Bias
l_bias = rho_a;
% SNR_dB_bias = 5;
% SNR_bias = (10^SNR_dB_bias) / 10;
% s_bias = sig_a / sqrt(SNR_bias);
s_bias = 0;
flag_bias = 0;

% Assigning initial points (noisy situation)
ind_Int = ceil(K^2/2-K/2);
% ind_Int = [ceil(K^2/2-K),ceil(K^2/2-K/2),ceil(K^2/2+K)]';
K_Int = size(ind_Int,1);
x_Int = x_K(ind_Int,:);

p = 5;  % number of sensors that we want to select greedy

% setting delta parameter
jj_max=1;
% j_step= 5;
delta = [0:15];
Del_size = size(delta,2);
% delta = j_step * [1:jj_max];
% jj_delta = 1;
%%% parameters
% params_a
params.rho_a = rho_a;
params.sig_a = sig_a;
% params_n
params.rho_n = rho_n;
params.sig_n = sig_n;
params.sigmaS = S_t;

mont_a_max = 3;
mont_m_max = 3;
%%

for mont_a = 1:mont_a_max
    
    
    aa = Generate_GM_2D(rho_a, sig_a, x_K, zeros(size(x_K,1),1), n) ;
    a_K_2D = reshape(aa,[K,K]);
    a_K = aa;
    
    figure(987)
    subplot(2,1,1)
    surf(x_K1,x_K2,a_K_2D)
    
    for mont_m = 1 : mont_m_max
        
        
        uncertainty_a = Generate_GM_2D(l_uncertainty, s_uncertainty, x_K, 0, n);
        bias = Generate_GM_2D(l_bias, s_bias, x_K, 0, n);
        m_K = a_K + (flag_uncertainty * uncertainty_a) + (flag_bias * bias);
        
        m_K_2D = reshape(m_K,[K,K]);
        figure(987)
        subplot(2,1,2)
        surf(x_K1,x_K2,m_K_2D)
        
        %     la = roh_a;
        %     sa2 = 3*sig_a;      % uncertainty on spatial gain 'a'
        a_hat_K  = m_K;
        
        %%% Noise %%%
        
        n_K = Generate_GM_2D(rho_n, sig_n, x_K, 0, n) ;
        
        %%% Measurements %%%
        y_K = a_K * S_t + n_K;
        
        
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
        Index = (1:K^2)';
        ind_rest = setdiff(Index,ind_Int);
        x_rest = x_K (ind_rest,:);
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
        Cn = My_2D_Cov ( rho_n, sig_n, x_Int, x_Int);
        f_hat = pinv(Cn)*m_Int;
        SNR_True_Estimated_int = (S_t^2 * f_hat' * a_Int * a_Int' * f_hat ) / (f_hat' * Cn * f_hat)
        
        
        SNR_True_Estimated = zeros(nbPts,1);
        SNR_True_Estimated(ind_Int) = SNR_True_Estimated_int;
        
        for i_rest = 1 : K^2-K_Int
            x_K_plus_N = [x_Int;x_rest(i_rest,:)];
            m_k_plus_N = [m_Int;m_rest(i_rest)];
            a_k_plus_N = [a_Int;a_rest(i_rest)];
            Cn = My_2D_Cov ( rho_n, sig_n, x_K_plus_N, x_K_plus_N);
            f_hat = pinv(Cn)*m_k_plus_N;
            SNR_True_Estimated_rest(i_rest,1) = (S_t^2 * f_hat' * a_k_plus_N * a_k_plus_N' * f_hat ) / (f_hat' * Cn * f_hat);
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROPOSED_PLUS: CostFunction_CDF_Complementary  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        for Del = 1:Del_size
            j_step = delta(Del);
            
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
                    
                    
                    tic
                    
                    % select new position
                    J_N = CostFunction_CDF_Complementary_2D(ind_Int_proposed_plus, ind_rest_proposed_plus, K^2, x_K, m_K, Ka,...
                        a_Int_proposed_plus, m_Int_proposed_plus, Kn, a_K,jj,j_step);
                    
                    t_proposed_plus(j_p,mont_a,mont_m)=toc;
                end
                
                [J_proposed_plusn , Index_proposed_plus] = max(J_N);
                x_new_proposed_plus = x_rest_proposed_plus(Index_proposed_plus,:);
                m_new_proposed_plus = m_rest_proposed_plus(Index_proposed_plus);
                a_new_proposed_plus = a_rest_proposed_plus(Index_proposed_plus);
                y_new_proposed_plus = y_rest_proposed_plus(Index_proposed_plus);
                ind_new_proposed_plus = ind_rest_proposed_plus(Index_proposed_plus);
                %         a_hat_new_proposed_plus = a_hat_rest_proposed_plus(Index_max);
                
                % 'a' estimation at new selected position
                C_NK_new = My_2D_Cov ( rho_a, sig_a, x_new_proposed_plus, x_Int_proposed_plus );
                C_KK_new = My_2D_Cov ( rho_a, sig_a, x_Int_proposed_plus, x_Int_proposed_plus ) ;
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
                C_K_plus_N = My_2D_Cov ( rho_n, sig_n, x_K_plus_N_proposed_plus, x_K_plus_N_proposed_plus );
                L_C_K_plus_N = chol(C_K_plus_N + 1e-10*eye(size(C_K_plus_N)), 'lower');
                
                %         S_hat_proposed_plus(j) = 1 / ( a_K_plus_N_proposed_plus' * (  L_C_K_plus_N.'  \   ( L_C_K_plus_N \ a_K_plus_N_proposed_plus ) ) ) * a_K_plus_N_proposed_plus' * ( L_C_K_plus_N.'  \ (L_C_K_plus_N \ y_K_plus_N_proposed_plus) );
                
                % SNR(f_hat)
                a_True_K_plus_N_proposed_plus = a_K(ind_K_plus_N_proposed_plus);
                f_hat_proposed_plus1 = pinv(C_K_plus_N) * a_K_plus_N_proposed_plus;
                f_hat_proposed_plus2 = pinv(C_K_plus_N) * a_hat_K_plus_N_proposed_plus;
                SNR_f_hat_proposed_plus1(j_p,Del) = ( S_t^2 * (f_hat_proposed_plus1' * a_True_K_plus_N_proposed_plus)^2 ) / (f_hat_proposed_plus1' * C_K_plus_N * f_hat_proposed_plus1 );
                SNR_f_hat_proposed_plus2(j_p,Del) = ( S_t^2 * (f_hat_proposed_plus2' * a_True_K_plus_N_proposed_plus)^2 ) / (f_hat_proposed_plus2' * C_K_plus_N * f_hat_proposed_plus2 );
                
                
                % update data
                
                ind_Int_proposed_plus = ind_K_plus_N_proposed_plus;
                x_Int_proposed_plus = x_K_plus_N_proposed_plus;
                m_Int_proposed_plus = [m_Int_proposed_plus; m_new_proposed_plus] ;
                a_Int_proposed_plus = a_K_plus_N_proposed_plus;
                a_hat_Int_proposed_plus = a_hat_K_plus_N_proposed_plus;
                y_Int_proposed_plus = y_K_plus_N_proposed_plus;
                
                ind_rest_proposed_plus (Index_proposed_plus) = [];
                x_rest_proposed_plus (Index_proposed_plus,:) = [];
                m_rest_proposed_plus (Index_proposed_plus) = [];
                a_rest_proposed_plus (Index_proposed_plus) = [];
                y_rest_proposed_plus (Index_proposed_plus) = [];
                %             a_hat_rest_proposed_plus (Index_proposed_plus) = [];
                
                
            end
            
            
        end
        
        %     S_hat_proposed_plus_mont(mont,:) = S_hat_proposed_plus;
        SNR_f_hat_proposed_plus_mont_m1(mont_m,:) = SNR_f_hat_proposed_plus1(end,:);
        SNR_f_hat_proposed_plus_mont_m2(mont_m,:) = SNR_f_hat_proposed_plus2(end,:);
        
        
        
    end
    %%
    
    
    SNR_f_hat_proposed_plus_mean_mont_m1(mont_a,:) = mean(SNR_f_hat_proposed_plus_mont_m1,1);
    
    SNR_f_hat_proposed_plus_mean_mont_m2(mont_a,:) = mean(SNR_f_hat_proposed_plus_mont_m2,1);
    
end

SNR_Del = mean(SNR_f_hat_proposed_plus_mean_mont_m2,1);


%%
figure(13)
lw = 3;
fs = 30;
plot(delta,smooth(SNR_Del),'LineWidth', lw,'color',[0 0.4470 0.7410]); hold on;
ylabel('$\overline{SNR (\hat{f})}~$[dB]','interpreter','latex','FontSize', fs);
xlabel('$\delta$', 'FontSize', fs,'interpreter','latex');
set(gca,'FontSize', fs);
grid on
