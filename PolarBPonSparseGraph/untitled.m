%converts .mat to .txt H matrix

clear all;
clc;

% k = size(H,1)
% n = size(H,2)
% load('H256.mat')
% H2=zeros(379,507);

% load('H64.mat')
% H2=zeros(76,108);

% load('H64_Optimized.mat')
% H2=zeros(84,116);

%y=load('H_RPrev_Optimizedpoly1.mat')
% load('H_RPrev_Optimizedpoly1.mat')
% H2=zeros(85,117);

% y=load('N64_K32_H_RP_pop344.mat')
% load('N64_K32_H_RP_pop344.mat')
% H2=zeros(77,109);

% y=load('RP_H_MaxCN8_Pop300_N128_K64.mat')
% load('RP_H_MaxCN8_Pop300_N128_K64.mat')
% H2=zeros(191,255);

% y=load('RP_H_MaxCN8_Pop161_N64_K32.mat')
% load('RP_H_MaxCN8_Pop161_N64_K32.mat')
% H2=zeros(77,109);

% y=load('RP_H_MaxCN8_Polar_N64_K32.mat')
% load('RP_H_MaxCN8_Polar_N64_K32.mat')
% H2=zeros(84,116);

% load('H64_NOcheckComb.mat')
% H2=zeros(104,136);

% [I,J,S] = find(H);


% load('OptimizedRPPop200_BP_PRP_N128_K64_maxBP_iter200_MaxC8.mat')
% [I,J,S] = find(res.H);
% size(res.H)
% H2=zeros(191,255);

% load('RP_H_Girth_BER2.mat') %length 64
% load('RP_H_Girth_BER.mat') %length 128
% load('RP_H_Girth_BER3.mat') %length 256
% H=data.H_MATRIX{53}

% load('N128_K64_HgeneticBest.mat')
load('N128_K64_Hpolar.mat')
[I,J,S]=find(H);
for ii = 1:length(I)
    H2(I(ii),J(ii))=1;
end

fileID = fopen('H2.txt','w');
fprintf(fileID,'%d ',H2');
fclose(fileID);

size(H2)