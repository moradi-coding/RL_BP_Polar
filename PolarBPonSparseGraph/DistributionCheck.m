
clear all;
clc;

m=249; n=498; 
% m=68; n=510; 
%m=132; n=990;
% m=93; n=186;

%---------------------- Distribution for G matrix ------------------------%

%Kelley/Kliewer irregular sub-code 1 distribution
% Landa=[0 4/16 5/16 4/16 3/16]; %row_wt distribution (wt. 1,2...)
% Ru=[0 0 0 2/10 2/10 6/10]; %col_wt distribution (wt. 1,2...)

% regular codes
% Landa=[0 1];  Ru=[0 0 1]; % (2,3) regular code
Landa=[0 0 1];  Ru=[0 0 0 0 0 1]; % (3,6) regular code
%Landa=[0 0 1]; Ru=[0 0 0 0 0 0 0 0 0 0 0 1]; % (3,12) regular code
% Landa=[0 1]; Ru=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]; % (2,15) regular code

G=degreeDist(n,m,Ru,Landa);
size(G)
%load('G.mat');

%     totalOnes = sum(sum(G));
%     degree = sum(G,1);
%     degree=sort(degree);
%     [~,j] = find(degree==degree(1));
%     k=1;
%     percentage = zeros(2,length(Landa));
%     while(j<=length(degree))
%         percentage(1,k) =  degree(j(end));
%         percentage(2,k) =  length(j)*degree(j(end))/totalOnes;
%         k = k+1;
%         if (j(end)+1 <= length(degree))
%             [~,j] = find(degree==degree(j(end)+1));
%         else
%             break;
%         end
%     end
%     str2='';
%     for i=1:k-1
%         if(i==1)
%             str2=strcat(str2,'Degree=');
%         else
%             str2=strcat(str2,',  Degree=');
%         end
%         str1 = num2str(percentage(1,i));
%         str1=strcat(str1,':');
%         str1=strcat(str1,num2str(percentage(2,i)));
%         str1=strcat(str1,'%');
%         str2=strcat(str2,str1);
%     end
%     hist(degree);
% %     text(1,-40,str2);
%     title(str2);
    save('G','G');
%     
    fileID = fopen('mat.txt','w');
    fprintf(fileID,'%d ',G');
    fclose(fileID);
%     G
     
    for i=1:m
        rw(i)=sum(G(i,:));
    end
    
    for i=1:n
        cw(i)=sum(G(:,i));
    end 
    
    cw
    rw
    
%     max(cw)
%     max(rw)
%     m
%     n
    
%  makeSCfromBlockV2();
%  fileID = fopen('col_wt_kk.txt','w');
%  fprintf(fileID,'%d ',max(cw));
%  fclose(fileID);
%  
%  fileID = fopen('row_wt_kk.txt','w');
%  fprintf(fileID,'%d ',max(rw));
%  fclose(fileID);
%  
%  fileID = fopen('mval.txt','w');
%  fprintf(fileID,'%d ',m);
%  fclose(fileID);
%  
%  fileID = fopen('nval.txt','w');
%  fprintf(fileID,'%d ',n);
%  fclose(fileID);
 