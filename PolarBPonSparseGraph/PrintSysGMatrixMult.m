%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Code to print Generator matrox to file
%
%    Author: David Mitchell, dgmm@nmsu.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [gens] = PrintGMatrix(G);
clear all;
clc;

G=load('H64.mat')

k = 76;
n = 108;

%--------------------------------------
% save G-matrix to file


str=sprintf('H64.mat');
disp(str);

fp = fopen( str,'wt');

fprintf(fp,'n = %1d\n',n);
fprintf(fp,'k = %1d\n',k);

fprintf(fp,'G = [\n');

for i=1:k-1
    fprintf(fp,'[');
    for j=1:n-1
       fprintf(fp,'%1d, ',G(i,j));
    end
    fprintf(fp,'%1d],\n ',G(i,n));
end

fprintf(fp,'[');
for j=1:n-1
      fprintf(fp,'%1d, ',G(k,j));
end
fprintf(fp,'%1d]] ',G(k,n));
fclose(fp);

 