clear all
close all
clc

tic

%% Parameters
nbPts = 3e2;
x = linspace(0, 1, nbPts).';

% Parameters of a
la  = .2;
sa2 = 1;

% Parameters of the mean function
lma  = la;
sma2 = sa2*3;

% Parameters of noise
SNR_dB = 2;
ratio_rho = 1e-1;
SNR = (10^SNR_dB) / 10;

ln = ratio_rho * la;
sn2 = sa2 / sqrt(SNR);

% Covariance function
kSquare = @(x,s2,l)(s2 * exp(-squareform(pdist(x,'squaredeuclidean')) / (2*l^2)));
k = @(x,xP,s2,l) (s2 * exp(-(x-xP).^2 / (2*l^2)));



%%% SNR %%%


%% Selection X_(N-1) = [x_1, ... , x_{N-1}]
nbX = 3;

% indSelected = fix([1:nbX]/(nbX+1)*nbPts);
indSelected = [5,ceil(nbPts/2),nbPts-5];
%indSelected = randperm(nbPts);
indSelected = sort(indSelected(1:nbX));
xSelected   = x(indSelected);

for mont11 = 1 : 1
    
    
    
%% Generate noise
Kn = kSquare(x, sn2, ln);

%% Generate ma
Kma = kSquare(x, sma2, lma);
Lma = chol(Kma + 1e-8*eye(nbPts), 'lower');
ma  = Lma * randn(nbPts, 1);

maSelected = ma(indSelected);


%% Generate a
Ka = kSquare(x, sa2, la);
La = chol(Ka + 1e-8*eye(nbPts), 'lower');
a  = ma + La * randn(nbPts, 1);

% aSelected = ma(indSelected);
aSelected = a(indSelected);


%% Generate a_N | vec(a)_N-1)
indNonSelected = setdiff(1:nbPts, indSelected);

nbPtsTest = 100;

for indN = 1:nbPts
    xN(indN) = x(indN);
    
    if ismember(indN, indNonSelected)
        %indNonSelected   %(randomOrder(1));
        
        mu_aN(indN,1) = ma(indN) + Ka(indSelected, indN).' * ((Ka(indSelected,indSelected)+1e-8*eye(nbX)) \ (aSelected - maSelected));
        s2_aN(indN,1) = Ka(indN,indN) - Ka(indSelected, indN).' / (Ka(indSelected,indSelected)+1e-8*eye(nbX)) * Ka(indSelected, indN);
        
        aN(indN, 1:nbPtsTest) = mu_aN(indN) + sqrt(s2_aN(indN)) * randn(1, nbPtsTest);
        
    else
        mu_aN(indN,1) = a(indN);
        s2_aN(indN,1) = 0;
        
        aN(indN, 1:nbPtsTest) = a(indN);
        
    end
    
end


%% True SNR

SNR_true = zeros(nbPts,1);
SNR_true(indSelected) = aSelected.' / (Kn(indSelected,indSelected)+1e-8*eye(nbX)) * aSelected;
aSelected_True = a(indSelected);

for ind = 1:length(indNonSelected)
    indN = indNonSelected(ind);      %indNonSelected(randomOrder(1));
    a_total = [aSelected_True;a(indN)];
    x_total = [x(indSelected) ; x(indN)];
    Rn1 = pinv( kSquare(x_total, sn2, ln) );
    SNR_true(ind,1) = a_total' * Rn1 * a_total;
end

SNR_true_norm = (SNR_true - min(SNR_true)) / (max(SNR_true)-min(SNR_true));

%% True SNR with estimated f

% True SNR(f_hat)

zK = ma(indSelected);
tol = 1e-10;
flagPrior = 1;
xK   = x(indSelected);
ma_K = ma(indSelected);
sigmaS = 2;  % power of s


% 1. Estimation of a(x) from z(x)
if flagPrior == 1
    aEst = ma;
    
    [xx_, xy_] = meshgrid(x, x);
    Sigma_aEst = k(xx_, xy_, la, sa2);
    
else
    [xx_kk, xy_kk] = meshgrid(xK, xK);
    SigmaA_kk = k(xx_kk, xy_kk, la, sa2);
    
    [xx_k, xy_k] = meshgrid(xK, x);
    SigmaA_k  = k(xx_k, xy_k, la, sa2);
    
    aEst = ma + SigmaA_k * (SigmaA_kk \ (zK - ma_K));
    
    [xx_, xy_] = meshgrid(x, x);
    SigmaA_k  = k(xx_, xy_, la, sa2);
    
    Sigma_aEst = SigmaA_k - SigmaA_k * (SigmaA_kk \ SigmaA_k.');
end


% Computation of the true SNR with estimated a(x)

% SNR_True_Estimated_f = zeros(nbPts);

for indSensor1 = 1:nbPts
    for indSensor2 = 1:indSensor1
        indexExtendedSensors = unique([indSelected' ; indSensor1 ; indSensor2]);
        
        xExtended     = x(indexExtendedSensors);
        aExtended     = a(indexExtendedSensors);
        aEst_Extended = ma(indexExtendedSensors);
        
        xDist                   = squareform(pdist(xExtended, 'squaredeuclidean'));
        Rn_Extended_computation = sn2^2 * (exp( -xDist / (2*ln^2) ) + tol*eye(length(xExtended)));
        SigmaA_Extended         = Sigma_aEst(indexExtendedSensors, indexExtendedSensors);
        
        
        aTmp = Rn_Extended_computation \ aExtended;
        SNR_True_Estimated_f(indSensor1, indSensor2) = ...
            sigmaS^2 * ( (aTmp.' * aEst_Extended)^2 + aTmp.' * SigmaA_Extended * aTmp ) / ...
            ( (aEst_Extended.' / Rn_Extended_computation) * aEst_Extended + ...
            trace(Rn_Extended_computation \ SigmaA_Extended) );
    end
end

tmp = tril(SNR_True_Estimated_f, -1);   % Complete the upper part of SNR_True (should be symmetric matrix)
SNR_True_Estimated_f = SNR_True_Estimated_f + tmp.';
SNR_True_Estimated_diag = diag(SNR_True_Estimated_f);


toc

SNR_True_Estimated_int = SNR_True_Estimated_diag(indSelected(1));

SNR_True_Estimated_diag_norm_int = (SNR_True_Estimated_diag - SNR_True_Estimated_int) / (max(SNR_True_Estimated_diag) - SNR_True_Estimated_int);

index_failure = find(SNR_True_Estimated_diag_norm_int<0);

%% Compute SNR mean (icassp)

%   SNRExpectationTh = zeros(nbPts,1);

for ind = 1:length(indNonSelected)
    indN = indNonSelected(ind);      %indNonSelected(randomOrder(1));
    
    Rn     = eye(nbX+1) / (Kn([indSelected , indN], [indSelected , indN]) + 1e-8*eye(nbX+1));
    SNRExpectationTh(ind) = aSelected.' * Rn(1:nbX,1:nbX) * aSelected ...
        + 2 * aSelected.' * Rn(1:nbX,end) * mu_aN(indN) ...
        + Rn(end,end) * (s2_aN(indN) + mu_aN(indN)^2);
end

% SNRExpectationTh(indSelected) = aSelected.' / (Kn(indSelected,indSelected)+1e-8*eye(nbX)) * aSelected;



%% Compute SNR CDF
randomOrder = randperm(nbPts-nbX);

SNR_Initial     = aSelected.' / (Kn(indSelected,indSelected)+1e-8*eye(nbX)) * aSelected;
[SNR_ExpectedMin, indexMax] = max(SNRExpectationTh);


for ind = 1:length(indNonSelected)
    indN = indNonSelected(ind);      %indNonSelected(randomOrder(1));
    
    for indTest = 1:nbPtsTest
        SNR_N(indTest) = [aSelected ; aN(indN,indTest)].' / (Kn([indSelected indN],[indSelected indN]) + 1e-8*eye(nbX+1)) * [aSelected ; aN(indN,indTest)];
    end
    
    
    [SNRHist(:,ind), SNRHistAbsc(:,ind)] = hist(SNR_N, 35);
    SNRHist(:,ind) = SNRHist(:,ind) / (sum(SNRHist(:,ind)) * (SNRHistAbsc(2,ind)-SNRHistAbsc(1,ind)));
    
    
    
    Rn     = eye(nbX+1) / (Kn([indSelected , indN], [indSelected , indN]) + 1e-8*eye(nbX+1));
    alpha  = aSelected.' * Rn(1:nbX,1:nbX) * aSelected;
    beta   = aSelected.' * Rn(1:nbX,end);
    gamma  = Rn(end,end);
    lambda = (mu_aN(indN) + beta/gamma)^2 / s2_aN(indN);
    
    
    xSNRTh(:,ind) = linspace(min(SNRHistAbsc(:,ind)), max(SNRHistAbsc(:,ind)), 1e3).';
    SNR_pdfTh(:,ind) = 1/(gamma*s2_aN(indN)) * ncx2pdf(1/(gamma*s2_aN(indN))*(xSNRTh(:,ind) - (alpha - beta^2/gamma)), 1, lambda);
    
    
    probaImprovement(ind,1) = 1-ncx2cdf(1/(gamma*s2_aN(indN))*(SNR_ExpectedMin - (alpha - beta^2/gamma)), 1, lambda);
    probaDegradate(ind,1) = ncx2cdf(1/(gamma*s2_aN(indN))*(SNR_Initial - (alpha - beta^2/gamma)), 1, lambda);
    
    for g = 1:40
        thr(g) = .05 * g;
        threshold(g) = .05 * g * (SNR_ExpectedMin-SNR_Initial);
        probaDegradate_BIS(ind,g) = ncx2cdf(1/(gamma*s2_aN(indN))*(SNR_Initial + threshold(g) - (alpha - beta^2/gamma)), 1, lambda);
    end
    
    
    SNRExpectationTh_BIS(ind,1) = gamma*s2_aN(indN) * (1 + lambda) + (alpha - beta^2 / gamma);
    SNRVarTh_BIS(ind,1)         = (gamma*s2_aN(indN))^2 * 2*(1+2*lambda);
    
    
    
    SNRMean(ind) = mean(SNR_N);
    SNRVar(ind)  = var(SNR_N);
    
    
end



%%

indN = indexMax;%randomOrder(1);
indexPlot = [100 nbPts-100];

figure(1); clf
fill([x ; flipud(x)], [mu_aN+2*(s2_aN).^(1/2) ; flipud(mu_aN-2*(s2_aN).^(1/2))], .8*ones(1,3), 'edgeColor', 'none'); hold on

plot(x, a,'--','LineWidth',2.5); hold on
plot(xSelected, aSelected, 'o','MarkerSize',12,'LineWidth',2.5)
plot(xN, mu_aN, 'k','LineWidth',2.5);
plot(x(indN), aN(indN,:), '+r')
plot(x(indexPlot), aN(indexPlot,:), '+c')
grid on
xlabel('x')
ylabel('m_a(x)')
% legend('Uncertainty of a','a','intial points','mean(a)','estimated a at new position','FontSize',16);




figure(2); clf
set(gcf, 'Position',  [500, 500, 1200, 900])
subplot(131)
bar(SNRHistAbsc(:,indexPlot(1)), SNRHist(:,indexPlot(1))); hold on
plot(xSNRTh(:,indexPlot(1)), SNR_pdfTh(:,indexPlot(1)),'r','LineWidth',3)
title(['(a): Position index:', num2str(indexPlot(1))],'FontSize',14)
xlabel('SNR','FontSize',16)
ylabel('Frequency','FontSize',16)
set(gca, 'FontSize', 18)

subplot(132)
bar(SNRHistAbsc(:,indexPlot(2)), SNRHist(:,indexPlot(2))); hold on
plot(xSNRTh(:,indexPlot(2)), SNR_pdfTh(:,indexPlot(2)),'r','LineWidth',3)
title(['(b): Position index:',num2str(indexPlot(2))],'FontSize',14)
xlabel('SNR','FontSize',16)
ylabel('Frequency','FontSize',16)
set(gca, 'FontSize', 18)

subplot(133)
bar(SNRHistAbsc(:,indN), SNRHist(:,indN)); hold on
plot(xSNRTh(:,indN), SNR_pdfTh(:,indN),'r','LineWidth',3)
title(['(c): Position index:',num2str(indN)],'FontSize',14)
xlabel('SNR','FontSize',16)
ylabel('Frequency','FontSize',16)
set(gca, 'FontSize', 18)


figure(3); clf
axs3(1) = subplot(211);
fill([x(indNonSelected) ; flipud(x(indNonSelected))], ...
    [SNRExpectationTh(:)+2*(SNRVarTh_BIS(:)).^(1/2) ; flipud(SNRExpectationTh(:)-2*(SNRVarTh_BIS(:)).^(1/2))], .8*ones(1,3), 'edgeColor', 'none'); hold on
plot(x(indNonSelected), SNRExpectationTh_BIS,'k');
%     plot(x(indNonSelected), SNRExpectationTh);
plot(x(indNonSelected(indexMax)), SNRExpectationTh_BIS(indexMax), 'or');

grid on
xlabel('x')
ylabel('J(x) = mean(SNR(x))')
title(['SNR_{initial} = ' num2str(SNR_Initial)])

axs3(2) = subplot(212);
plot(x(indNonSelected), probaImprovement,'LineWidth',2.5); hold on
plt = plot(x(indNonSelected(indexMax)), probaImprovement(indexMax), 'o');
plt.HandleVisibility = 'off';

plot(x(indNonSelected), 1-probaDegradate_BIS,'LineWidth',2.5);
plt = plot(x(indNonSelected(indexMax)), 1-probaDegradate_BIS(indexMax), 'o');
plt.HandleVisibility = 'off';

plot(x(indNonSelected), probaDegradate,'LineWidth',2.5);
plt = plot(x(indNonSelected(indexMax)), probaDegradate(indexMax), 'o');
plt.HandleVisibility = 'off';

grid on
xlabel('x','FontSize',16)
ylabel('Proba improvement > max(average(SNR))','FontSize',16)
% legend('1-CDF(SNR_{max})', 'CDF(SNR''int+thr'')','CDF(SNR_{int})','FontSize',16)

%     linkaxes(axs3, 'x')


plotSNRTh = SNRExpectationTh - min(SNRExpectationTh);
plotSNRTh = plotSNRTh / max(plotSNRTh);

plotPdfTh = probaImprovement - min(probaImprovement);
plotPdfTh = plotPdfTh / max(plotPdfTh);

plot_SNR_true = SNR_true - min(SNR_true);
plot_SNR_true = plot_SNR_true / max(plot_SNR_true);

for mont = 1 : g
    plotPdfDegTh(:,mont) = smooth(smooth(1-probaDegradate_BIS(:,mont) - min(1-probaDegradate_BIS(:,mont))));
    plotPdfDegTh(:,mont) = plotPdfDegTh(:,mont) / max(plotPdfDegTh(:,mont));
end

%%
%%%%%% icassp_2019

% figure(5)
grid on
plot_SNR_Expectation = SNRExpectationTh - min(SNRExpectationTh);
plot_SNR_Expectation = plot_SNR_Expectation/ max(plot_SNR_Expectation);

p98 = plot(x(indNonSelected), smooth(smooth(plot_SNR_Expectation)),'m','LineWidth',3); hold on

% set(gca, 'FontSize', 20)
% xlabel('$x_N$','interpreter','latex','FontSize',28)
% ylabel('$\hat{g}(X_N|\mathbf{X}_K)$','interpreter','latex','FontSize',28)

legend({'$\delta=0.1$','$\hat{g}(X_N|\mathbf{X}_K)$','$SNR(\hat{f})$','$\theta$'},'Interpreter','latex')

% 
% theta1 = .75*ones(1,nbPts);
% theta2 = .5*ones(1,nbPts);
% theta3 = .25*ones(1,nbPts);
% theta4 = 0*ones(1,nbPts);
% 
% plot(x, theta1,'--k','LineWidth',2)
% plot(x, theta2,'--k','LineWidth',2)
% plot(x, theta3,'--k','LineWidth',2)
% plot(x, theta4,'--k','LineWidth',2)
%%
figure(101)

My_colors = [0 0.4470 0.7410;0.6350 0.0780 0.1840;0.4660 0.6740 0.1880;0.4940 0.1840 0.5560];
My_marker = {'o','x','>','s','d','^','v','<','p'};

theta1 = .7*ones(1,nbPts);
theta2 = .5*ones(1,nbPts);
theta3 = .2*ones(1,nbPts);
theta4 = 0*ones(1,nbPts);

f=10;
thr(f)
p01 = plot(x(indNonSelected), plotPdfDegTh(:,f),'Marker','o','MarkerIndices',1:40:1000,'Color',My_colors(1,:),'LineWidth',3,'MarkerSize',10); hold on
p02 = plot(x(indNonSelected), smooth(smooth(plot_SNR_Expectation)),'Color',My_colors(2,:),'LineWidth',3); hold on
p03 = plot(x, smooth(smooth(SNR_True_Estimated_diag_norm_int)),'--k','LineWidth',3,'Color',.9*[.5 .5 .5]); hold on
grid on

p04 = plot(x, theta2,'-.k','LineWidth',2);
p05 = plot(x, theta3,'-.k','LineWidth',2);
p06 = plot(x, theta4,'-.k','LineWidth',2);

set(gca, 'FontSize', 35)
xlabel('$x_N$','interpreter','latex','FontSize',40)
ylabel('Criterion','interpreter','latex','FontSize',35)

set(gcf, 'Position',  [500, 500, 1500, 1000])

legend([p01,p02,p03,p04],{'$1-G_N(W_M;\delta=0.0.75)$','$\hat{g}(X_N|\mathbf{X}_K)$','$SNR(\hat{f})$','$\theta$'},'Interpreter','latex')


%%  true SNR with estimated f
% 
% grid on
% % plot_True_SNR_Expectation_f = SNR_True_Estimated_diag - min(SNRExpectationTh);
% % plot_True_SNR_Expectation_f = SNR_True_Estimated_diag/ max(SNR_True_Estimated_diag);
% 
% plot(x, smooth(smooth(SNR_True_Estimated_diag_norm_int)),'b--','LineWidth',3); hold on
% grid on
% set(gca, 'FontSize', 20)
% xlabel('$x_N$','interpreter','latex','FontSize',28)
% ylabel('$SNR(\hat{f})$','interpreter','latex','FontSize',28)


% theta1 = .75*ones(1,nbPts);
% theta2 = .5*ones(1,nbPts);
% theta3 = .25*ones(1,nbPts);
% theta4 = 0*ones(1,nbPts);
% 
% plot(x, theta1,'--k','LineWidth',2)
% plot(x, theta2,'--k','LineWidth',2)
% plot(x, theta3,'--k','LineWidth',2)
% plot(x, theta4,'--k','LineWidth',2)
%% plot threshold


% theta1 = .95*ones(1,nbPts);
% theta2 = .8*ones(1,nbPts);
% theta3 = .55*ones(1,nbPts);
% 
% uiopen('/Users/ghayyemf/Desktop/FatemeGHAYYEM/Project/MyPhD/OSP2/Delta1.fig',1)
% plot(x, theta1,'k','LineWidth',2)
% plot(x, theta2,'k','LineWidth',2)
% plot(x, theta3,'k','LineWidth',2)
% legend('$SNR(\mathbf{f}^*)$','$G_N(\delta=0.25)$','$G_N(\delta=0.5)$','$G_N(\delta=0.75)$','$G_N(\delta=1)$','interpreter','latex','FontSize',22)
% figure(1);


%% Failure
clc
threshold = [1e-5 .2 .5 .75 .95];

for hh = 1 : size(threshold,2)
failur_coeficcient = threshold(hh);
% the true SNR
ROI = find(SNR_True_Estimated_diag_norm_int > failur_coeficcient);
Total_actual_positive = numel(ROI);

% Proposed
G1_1_cdf =  plotPdfDegTh(:,1);
Positive_G1 = find(G1_1_cdf > failur_coeficcient);
TP_G1(mont11,hh) = numel ( intersect(Positive_G1,ROI) );
FP_G1(mont11,hh) = numel ( setdiff(Positive_G1,ROI) );
Negative_G1 = find(G1_1_cdf < failur_coeficcient);
FN_G1(mont11,hh) = numel ( intersect(Negative_G1,ROI) );
TN_G1(mont11,hh) = numel ( setdiff(Negative_G1,ROI) );
Rcl_G1(mont11,hh) = 100 * TP_G1(mont11,hh) / Total_actual_positive;
Prc_G1(mont11,hh) = 100 * TP_G1(mont11,hh) / (TP_G1(mont11,hh)+FP_G1(mont11,hh));
Spc_G1(mont11,hh) = 100 * TN_G1(mont11,hh) / (TN_G1(mont11,hh)+FP_G1(mont11,hh));
F_Score_G1(mont11,hh) = 2 * (Rcl_G1(mont11,hh) * Prc_G1(mont11,hh)) / (Rcl_G1(mont11,hh) + Prc_G1(mont11,hh));
Hrm_G1(mont11,hh) = 2 * (Rcl_G1(mont11,hh) * Spc_G1(mont11,hh)) / (Rcl_G1(mont11,hh) + Spc_G1(mont11,hh));

G2_1_cdf =  plotPdfDegTh(:,4);
Positive_G2 = find(G2_1_cdf > failur_coeficcient);
TP_G2(mont11,hh) = numel ( intersect(Positive_G2,ROI) );
FP_G2(mont11,hh) = numel ( setdiff(Positive_G2,ROI) );
Negative_G2 = find(G2_1_cdf < failur_coeficcient);
FN_G2(mont11,hh) = numel ( intersect(Negative_G2,ROI) );
TN_G2(mont11,hh) = numel ( setdiff(Negative_G2,ROI) );
Rcl_G2(mont11,hh) = 100 * TP_G2(mont11,hh) / Total_actual_positive;
Prc_G2(mont11,hh) = 100 * TP_G2(mont11,hh) / (TP_G2(mont11,hh)+FP_G2(mont11,hh));
Spc_G2(mont11,hh) = 100 * TN_G2(mont11,hh) / (TN_G2(mont11,hh)+FP_G2(mont11,hh));
F_Score_G2(mont11,hh) = 2 * (Rcl_G2(mont11,hh) * Prc_G2(mont11,hh)) / (Rcl_G2(mont11,hh) + Prc_G2(mont11,hh));
Hrm_G2(mont11,hh) = 2 * (Rcl_G2(mont11,hh) * Spc_G2(mont11,hh)) / (Rcl_G2(mont11,hh) + Spc_G2(mont11,hh));

G3_1_cdf =  plotPdfDegTh(:,8);
Positive_G3 = find(G3_1_cdf > failur_coeficcient);
TP_G3(mont11,hh) = numel ( intersect(Positive_G3,ROI) );
FP_G3(mont11,hh) = numel ( setdiff(Positive_G3,ROI) );
Negative_G3 = find(G3_1_cdf < failur_coeficcient);
FN_G3(mont11,hh) = numel ( intersect(Negative_G3,ROI) );
TN_G3(mont11,hh) = numel ( setdiff(Negative_G3,ROI) );
Rcl_G3(mont11,hh)= 100 * TP_G3(mont11,hh) / Total_actual_positive;
Prc_G3(mont11,hh) = 100 * TP_G3(mont11,hh) / (TP_G3(mont11,hh)+FP_G3(mont11,hh));
Spc_G3(mont11,hh) = 100 * TN_G3(mont11,hh) / (TN_G3(mont11,hh)+FP_G3(mont11,hh));
F_Score_G3(mont11,hh) = 2 * (Rcl_G3(mont11,hh) * Prc_G3(mont11,hh)) / (Rcl_G3(mont11,hh) + Prc_G3(mont11,hh));
Hrm_G3(mont11,hh) = 2 * (Rcl_G3(mont11,hh) * Spc_G3(mont11,hh)) / (Rcl_G3(mont11,hh) + Spc_G3(mont11,hh));

G4_1_cdf =  plotPdfDegTh(:,10);
Positive_G4 = find(G4_1_cdf > failur_coeficcient);
TP_G4(mont11,hh) = numel ( intersect(Positive_G4,ROI) );
FP_G4(mont11,hh) = numel ( setdiff(Positive_G4,ROI) );
Negative_G4 = find(G4_1_cdf < failur_coeficcient);
FN_G4(mont11,hh) = numel ( intersect(Negative_G4,ROI) );
TN_G4(mont11,hh) = numel ( setdiff(Negative_G4,ROI) );
Rcl_G4(mont11,hh) = 100 * TP_G4(mont11,hh) / Total_actual_positive;
Prc_G4(mont11,hh) = 100 * TP_G4(mont11,hh) / (TP_G4(mont11,hh)+FP_G4(mont11,hh));
Spc_G4(mont11,hh) = 100 * TN_G4(mont11,hh) / (TN_G4(mont11,hh)+FP_G4(mont11,hh));
F_Score_G4(mont11,hh) = 2 * (Rcl_G4(mont11,hh) * Prc_G4(mont11,hh)) / (Rcl_G4(mont11,hh) + Prc_G4(mont11,hh));
Hrm_G4(mont11,hh) = 2 * (Rcl_G4(mont11,hh) * Spc_G4(mont11,hh)) / (Rcl_G4(mont11,hh) + Spc_G4(mont11,hh));

% icassp'19 SNR
G4_1_cdf = plot_SNR_Expectation;
Positive_Mean_SNR = find(plot_SNR_Expectation > failur_coeficcient);
TP_Mean_SNR(mont11,hh) = numel ( intersect(Positive_Mean_SNR,ROI) );
FP_Mean_SNR(mont11,hh) = numel ( setdiff(Positive_Mean_SNR,ROI) );
Negative_Mean_SNR = find(plot_SNR_Expectation < failur_coeficcient);
FN_Mean_SNR(mont11,hh) = numel ( intersect(Negative_Mean_SNR,ROI) );
TN_Mean_SNR(mont11,hh) = numel ( setdiff(Negative_Mean_SNR,ROI) );
Rcl_Mean_SNR(mont11,hh) = 100 * TP_Mean_SNR(mont11,hh) / Total_actual_positive;
Prc_Mean_SNR(mont11,hh) = 100 * TP_Mean_SNR(mont11,hh) / (TP_Mean_SNR(mont11,hh)+FP_Mean_SNR(mont11,hh));
Spc_Mean_SNR(mont11,hh) = 100 * TN_Mean_SNR(mont11,hh) / (TN_Mean_SNR(mont11,hh)+FP_Mean_SNR(mont11,hh));
F_Score_Mean_SNR(mont11,hh) = 2 * (Rcl_Mean_SNR(mont11,hh) * Prc_Mean_SNR(mont11,hh)) / (Rcl_Mean_SNR(mont11,hh) + Prc_Mean_SNR(mont11,hh));
Hrm_Mean_SNR(mont11,hh) = 2 * (Rcl_Mean_SNR(mont11,hh) * Spc_Mean_SNR(mont11,hh)) / (Rcl_Mean_SNR(mont11,hh) + Spc_Mean_SNR(mont11,hh));

end

mont11

end

%%
Rcl_G1_avg = nanmean(Rcl_G1,1)
Rcl_G2_avg = nanmean(Rcl_G2,1)
Rcl_G3_avg = nanmean(Rcl_G3,1)
Rcl_G4_avg = nanmean(Rcl_G4,1)
Rcl_Mean_SNR_avg = nanmean(Rcl_Mean_SNR,1)

Prc_G1_avg = nanmean(Prc_G1,1)
Prc_G2_avg = nanmean(Prc_G2,1)
Prc_G3_avg = nanmean(Prc_G3,1)
Prc_G4_avg = nanmean(Prc_G4,1)
Prc_Mean_SNR_avg = nanmean(Prc_Mean_SNR,1)


Spc_G1_avg =  nanmean(Spc_G1,1)
Spc_G2_avg =  nanmean(Spc_G2,1)
Spc_G3_avg =  nanmean(Spc_G3,1)
Spc_G4_avg =  nanmean(Spc_G4,1)
Spc_Mean_SNR_avg =  nanmean(Spc_Mean_SNR,1)


F_Score_G1_avg =  nanmean(F_Score_G1,1);
F_Score_G2_avg =  nanmean(F_Score_G2,1);
F_Score_G3_avg =  nanmean(F_Score_G3,1);
F_Score_G4_avg =  nanmean(F_Score_G4,1);
F_Score_Mean_SNR_avg =  nanmean(F_Score_Mean_SNR,1);


Hrm_G1_avg =  nanmean(Hrm_G1,1)
Hrm_G2_avg =  nanmean(Hrm_G2,1)
Hrm_G3_avg =  nanmean(Hrm_G3,1)
Hrm_G4_avg =  nanmean(Hrm_G4,1)
Hrm_Mean_SNR_avg =  nanmean(Hrm_Mean_SNR,1)
%%
 close all
figure(400); clf
dd = [5,10,19];
thr(dd)
My_colors = [0 0.4470 0.7410;0.6350 0.0780 0.1840;0.4660 0.6740 0.1880;0.4940 0.1840 0.5560];
My_marker = {'o','x','>','s','d','^','v','<','p'};

for f = 1:size(dd,2)
    lm =  My_marker{f};
    J_p = zeros(nbPts,1);
    ind_rest=setdiff(1:nbPts,indSelected);
    J_p(ind_rest)=plotPdfDegTh(:,dd(f));
    
        p(f) = plot(x, J_p,'Color',My_colors(f,:),'LineWidth',3,'MarkerSize',10); hold on
        
end
grid on
set(gca, 'FontSize', 35)
% legend('SNR_{true}','SNR(x_{N})','thr=0.25','thr=0.5','thr=0.75','thr=1','FontSize',16)
% legend('$SNR(\mathbf{f}^*)$','$G_N(\delta=0.25)$','$G_N(\delta=0.5)$','$G_N(\delta=0.75)$','$G_N(\delta=1)$','interpreter','latex','FontSize',22)
xlabel('$x_N$','interpreter','latex','FontSize',40)
ylabel('$\tilde{J}_P(\mathbf{X}_N,\delta | \mathbf{a}_K)$','interpreter','latex','FontSize',35)


%%%%%% True SNR with true f
% plot(x(indNonSelected), smooth(smooth(plot_SNR_true)),':','LineWidth',3); hold on
% plot(x(indNonSelected), plotSNRTh,':','LineWidth',2.5); hold on
%     plot(x(indNonSelected), plotPdfTh,'LineWidth',2.5); hold on


%%%%%%  true SNR with estimated f
hold on
grid on
%plot_True_SNR_Expectation_f = SNR_True_Estimated_diag - min(SNRExpectationTh);
% plot_True_SNR_Expectation_f = SNR_True_Estimated_diag/ max(SNR_True_Estimated_diag);

p222 = plot(x, SNR_True_Estimated_diag_norm_int,'--k','LineWidth',3,'Color',.9*[.5 .5 .5]); hold on
% grid on
% set(gca, 'FontSize', 20)
% xlabel('$x_N$','interpreter','latex','FontSize',28)
% ylabel('$SNR(\hat{f})$','interpreter','latex','FontSize',28)
J_E=zeros(nbPts,1);
plot_SNR_Expectation = SNRExpectationTh - min(SNRExpectationTh);
plot_SNR_Expectation = plot_SNR_Expectation/ max(plot_SNR_Expectation);
J_E(indNonSelected) = plot_SNR_Expectation;

p02 = plot(x, J_E,'Color',.9*[.5 .5 .5],'LineWidth',2); hold on

theta1 = .7*ones(1,nbPts);
theta2 = .5*ones(1,nbPts);
theta3 = .2*ones(1,nbPts);
theta4 = 0*ones(1,nbPts);
% % 
% p77 = plot(x, theta2,'-k','LineWidth',.5);
% % plot(x, theta2,'-k','LineWidth',.5);
% p22 = plot(x, theta3,'-k','LineWidth',.5);
plot(x, theta4,'-k','LineWidth',.75);

p23 = plot(x(indSelected), [0,0,0],'o','Color',.95*[.5 .5 .5],'LineWidth',3,'MarkerSize',25);

legend([p,p02,p222,p23],{['$\delta=$',num2str(thr(dd(1)))],['$\delta=$',num2str(thr(dd(2)))],['$\delta=$',num2str(thr(dd(3)))],'$\tilde{J}_E(X_N|\mathbf{X}_K)$','$\widetilde{SNR}(\hat{f})$','Initial positions'},'Interpreter','latex')
