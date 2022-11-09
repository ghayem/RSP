clear all
close all
clc

tic
%%
nbMcMc = 1;
for Mont = 1:nbMcMc
%% Parameters
nbPts = 100;
x = linspace(0, 1, nbPts).';

% Parameters of a
la  = .2;
sa2 = .15;

% Parameters of the mean function
lma  = la;
sma2 = sa2;

% Parameters of noise
ln  = .1;
sn2 = .5;

% Covariance function
kSquare = @(x,s2,l)(s2 * exp(-squareform(pdist(x,'squaredeuclidean')) / (2*l^2)));
k = @(x,xP,s2,l) (s2 * exp(-(x-xP).^2 / (2*l^2)));



%% Selection X_(N-1) = [x_1, ... , x_{N-1}]
nbX = 3;

indSelected = fix([1:nbX]/(nbX+1)*nbPts);
%indSelected = randperm(nbPts);
indSelected = sort(indSelected(1:nbX));
xSelected   = x(indSelected);


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

aSelected_True = a(indSelected);
for ind = 1:length(indNonSelected)
    indN = indNonSelected(ind);      %indNonSelected(randomOrder(1));
    a_total = [aSelected_True;a(indN)];
    x_total = [x(indSelected) ; x(indN)];
    Rn1 = pinv( kSquare(x_total, sn2, ln) );
    SNR_true(ind,1) = a_total' * Rn1 * a_total;
end


%% True SNR(f_hat)

zK = ma(indSelected);
tol = 1e-10;
flagPrior = 1
xK   = x(indSelected);
ma_K = ma(indSelected);
sigmaS = 2;  % power of s


% 1. Estimation of a(x) from z(x)
if flagPrior == 1
    aEst = ma;
    
    [xx_, xy_] = meshgrid(x, x);
    Sigma_aEst = k(xx_, xy_, sa2, la);
    
else
    [xx_kk, xy_kk] = meshgrid(xK, xK);
    SigmaA_kk = k(xx_kk, xy_kk, sa2, la);
    
    [xx_k, xy_k] = meshgrid(xK, x);
    SigmaA_k  = k(xx_k, xy_k, sa2, la);
    
    aEst = ma + SigmaA_k * (SigmaA_kk \ (zK - ma_K));
    
    [xx_, xy_] = meshgrid(x, x);
    SigmaA_k  = k(xx_, xy_, sa2, la);
    
    Sigma_aEst = SigmaA_k - SigmaA_k * (SigmaA_kk \ SigmaA_k.');
end


% Computation of the true SNR with estimated a(x)

SNR_True_Estimated = zeros(nbPts);

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
        SNR_True_Estimated(indSensor1, indSensor2) = ...
            sigmaS^2 * ( (aTmp.' * aEst_Extended)^2 + aTmp.' * SigmaA_Extended * aTmp ) / ...
            ( (aEst_Extended.' / Rn_Extended_computation) * aEst_Extended + ...
            trace(Rn_Extended_computation \ SigmaA_Extended) );
    end
end

tmp = tril(SNR_True_Estimated, -1);   % Complete the upper part of SNR_True (should be symmetric matrix)
SNR_True_Estimated = SNR_True_Estimated + tmp.';
SNR_True_Estimated_diag = diag(SNR_True_Estimated);


toc

SNR_True_Estimated_int = SNR_True_Estimated_diag(indSelected(1));

SNR_True_Estimated_diag_norm_int = (SNR_True_Estimated_diag - SNR_True_Estimated_int) / (max(SNR_True_Estimated_diag) - SNR_True_Estimated_int);

plot(SNR_True_Estimated_diag_norm_int, 'LineWidth', 2.5); hold on

index_failure = find(SNR_True_Estimated_diag_norm_int<0);

%% Compute SNR mean

% SNRExpectationTh = zeros(nbPts,1);

Rn_KK     = eye(nbX) / (Kn(indSelected , indSelected) + 1e-8*eye(nbX));
Ka_KK = Ka(indSelected , indSelected);

% SNRExpectationTh_int = sigmaS^2 * ( trace(Rn_KK * Ka_KK) + ma_K' * Rn_KK * ma_K);
% SNRExpectationTh_int = sigmaS^2 * ( ma_K' * Rn_KK * ma_K);

% SNRExpectationTh(indSelected,1) = repmat(SNRExpectationTh_int, nbX, 1);

for ind = 1:length(indNonSelected)
    indN = indNonSelected(ind);      %indNonSelected(randomOrder(1));
    
    Rn_KN     = eye(nbX+1) / (Kn([indSelected , indN], [indSelected , indN]) + 1e-8*eye(nbX+1));
    SNRExpectationTh(ind,1) = sigmaS^2 * ( aSelected.' * Rn_KN(1:nbX,1:nbX) * aSelected ...
        + 2 * aSelected.' * Rn_KN(1:nbX,end) * mu_aN(indN) ...
        + Rn_KN(end,end) * (s2_aN(indN) + mu_aN(indN)^2) );
end

SNRExpectationTh_norm_int= zeros(nbPts,1);

SNRExpectationTh_norm_int(indNonSelected,1) = (SNRExpectationTh - min(SNRExpectationTh)) / (max(SNRExpectationTh) - min(SNRExpectationTh));


plot(SNRExpectationTh_norm_int, 'LineWidth', 2.5); hold on






%% Compute SNR CDF

SNR_Initial = sigmaS^2 * ( aSelected.' / (Kn(indSelected,indSelected)+1e-8*eye(nbX)) * aSelected);
                
[SNR_ExpectedMax, indexMax] = max(SNRExpectationTh);

% probaDegradate_BIS = zeros(nbPts,1);


for ind = 1:length(indNonSelected)
    indN = indNonSelected(ind);      %indNonSelected(randomOrder(1));
    
    for indTest = 1:nbPtsTest
        SNR_N(indTest) = [aSelected ; aN(indN,indTest)].' / (Kn([indSelected indN],[indSelected indN]) + 1e-8*eye(nbX+1)) * [aSelected ; aN(indN,indTest)];
    end
    
    
    [SNRHist(:,ind), SNRHistAbsc(:,ind)] = hist(SNR_N, 35);
    SNRHist(:,ind) = SNRHist(:,ind) / (sum(SNRHist(:,ind)) * (SNRHistAbsc(2,ind)-SNRHistAbsc(1,ind)));
    
    
    
    Rn_KN     = eye(nbX+1) / (Kn([indSelected , indN], [indSelected , indN]) + 1e-8*eye(nbX+1));
    alpha  = aSelected.' * Rn_KN(1:nbX,1:nbX) * aSelected;
    beta   = aSelected.' * Rn_KN(1:nbX,end);
    gamma  = Rn_KN(end,end);
    lambda = (mu_aN(indN) + beta/gamma)^2 / s2_aN(indN);
    
    
    xSNRTh(:,ind) = linspace(min(SNRHistAbsc(:,ind)), max(SNRHistAbsc(:,ind)), 1e3).';
    SNR_pdfTh(:,ind) = 1/(gamma*s2_aN(indN)) * ncx2pdf(1/(gamma*s2_aN(indN))*(xSNRTh(:,ind) - (alpha - beta^2/gamma)), 1, lambda);
    
    
    probaImprovement(ind,1) = 1-ncx2cdf(1/(gamma*s2_aN(indN))*(SNR_ExpectedMax - (alpha - beta^2/gamma)), 1, lambda);
    probaDegradate(ind,1) = ncx2cdf(1/(gamma*s2_aN(indN))*(SNR_Initial - (alpha - beta^2/gamma)), 1, lambda);
    
    for g = 1:4
        threshold(g) = 1 * g * (SNR_ExpectedMax-SNR_Initial);
        probaDegradate_BIS(ind,g) = ncx2cdf(1/(gamma*s2_aN(indN))*(SNR_Initial + threshold(g) - (alpha - beta^2/gamma)), 1, lambda);
    end
    

end


CDF_Complementary = 1-probaDegradate_BIS(:,2);
CDF_min = min(CDF_Complementary);

CDF_Complementary_norm= zeros(nbPts,1);
CDF_Complementary_norm(indNonSelected,1) = (CDF_Complementary-CDF_min) / (max(CDF_Complementary) - CDF_min);

plot(CDF_Complementary_norm, 'LineWidth', 2.5),hold on

lgd = legend('$SNR(\hat{f})$' , '$J_E(X_N|\mathbf{X}_K) $' , '$J_P(X_N,\delta|\mathbf{X}_K)$','Interpreter' , 'latex');
lgd.Interpreter = 'latex';
lgd.FontSize = 24
grid on
set(gca,'FontSize', 18 )
xlabel('$X$','Interpreter' , 'latex', 'FontSize', 24);
ylabel('Criterion', 'FontSize', 24);

plot(0:nbPts,zeros(nbPts+1,1),'k', 'LineWidth', 1, 'HandleVisibility','off')
%%

% normalization _ Oracle
plot_SNR_true = SNR_true - min(SNR_true);
plot_SNR_true = plot_SNR_true / max(plot_SNR_true);

% normalization _ proposed
plotPdfTh = probaImprovement - min(probaImprovement);
plotPdfTh = plotPdfTh / max(plotPdfTh);


for h = 1 : 4
    plotPdfDegTh(:,h) = smooth(smooth(probaDegradate_BIS(:,h) - min(probaDegradate_BIS(:,h))));
    plotPdfDegTh(:,h) = plotPdfDegTh(:,h) / max(plotPdfDegTh(:,h));
end

% normalization _ Mean SNR _ (icassp 2019)
plot_SNR_Expectation = SNRExpectationTh - min(SNRExpectationTh);
plot_SNR_Expectation = plot_SNR_Expectation/ max(plot_SNR_Expectation);


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

eps = 1e-2; 

failur_coeficcient = .5;
% the true SNR
Positive_oracle = find(plot_SNR_true > failur_coeficcient);
Total_actual_positive = numel(Positive_oracle);

% Proposed
G1_1_cdf = 1 - plotPdfDegTh(:,1);
Positive_G1 = find(G1_1_cdf > failur_coeficcient);
TP_G1 = numel ( intersect(Positive_G1,Positive_oracle) );
FP_G1 = numel ( setdiff(Positive_G1,Positive_oracle) );
Rcl_G1= 100 * TP_G1 / Total_actual_positive;
Prc_G1 = 100 * TP_G1 / (TP_G1+FP_G1);
F_Score_G1(Mont) = eps + 2 * (Rcl_G1 * Prc_G1) / (Rcl_G1 + Prc_G1);

G2_1_cdf = 1 - plotPdfDegTh(:,2);
Positive_G2 = find(G2_1_cdf > failur_coeficcient);
TP_G2 = numel ( intersect(Positive_G2,Positive_oracle) );
FP_G2 = numel ( setdiff(Positive_G2,Positive_oracle) );
Rcl_G2 = 100 * TP_G2 / Total_actual_positive;
Prc_G2 = 100 * TP_G2 / (TP_G2+FP_G2);
F_Score_G2(Mont) = eps + 2 * (Rcl_G2 * Prc_G2) / (Rcl_G2 + Prc_G2);

G3_1_cdf = 1 - plotPdfDegTh(:,3);
Positive_G3 = find(G3_1_cdf > failur_coeficcient);
TP_G3 = numel ( intersect(Positive_G3,Positive_oracle) );
FP_G3 = numel ( setdiff(Positive_G3,Positive_oracle) );
Rcl_G3= 100 * TP_G3 / Total_actual_positive;
Prc_G3 = 100 * TP_G3 / (TP_G2+FP_G3);
F_Score_G3(Mont) = eps + 2 * (Rcl_G3 * Prc_G3) / (Rcl_G3 + Prc_G3);

G4_1_cdf = 1 - plotPdfDegTh(:,4);
Positive_G4 = find(G4_1_cdf > failur_coeficcient);
TP_G4 = numel ( intersect(Positive_G4,Positive_oracle) );
FP_G4 = numel ( setdiff(Positive_G4,Positive_oracle) );
Rcl_G4= 100 * TP_G4 / Total_actual_positive;
Prc_G4 = 100 * TP_G4 / (TP_G2+FP_G4);
F_Score_G4(Mont) = eps + 2 * (Rcl_G4 * Prc_G4) / (Rcl_G4 + Prc_G4);

% icassp'19 SNR
Positive_Mean_SNR = find(plot_SNR_Expectation > failur_coeficcient);
TP_Mean_SNR = numel ( intersect(Positive_Mean_SNR,Positive_oracle) );
FP_Mean_SNR = numel ( setdiff(Positive_Mean_SNR,Positive_oracle) );
Rcl_Mean_SNR= 100 * TP_Mean_SNR / Total_actual_positive;
Prc_Mean_SNR = 100 * TP_Mean_SNR / (TP_G2+FP_Mean_SNR);
F_Score_Mean_SNR(Mont) = 2 * (Rcl_Mean_SNR * Prc_Mean_SNR) / (Rcl_Mean_SNR + Prc_Mean_SNR);
%%

end

% F_Score_G1_average = nanmean(F_Score_G1,'all')
% F_Score_G2_average = nanmean(F_Score_G2,'all')
% F_Score_G3_average = nanmean(F_Score_G3,'all')
% F_Score_G4_average = nanmean(F_Score_G4,'all')
% F_Score_Mean_SNR_average = nanmean(F_Score_Mean_SNR,'all')
