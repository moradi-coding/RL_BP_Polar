function G=degreeDist(n, m, Ru, Landa) 

    [variable,E1] = getNoNodeDegree(Landa, n);
    [checkk,E2] = getNoNodeDegree(Ru, m);
    len = max(variable);
    done = ones(length(variable),1);
    link = zeros(length(variable), len);
    true=1;
    count=1;
    while(true==1)
        var = variable;
        che = checkk;
        for i=1:length(var)
            [~, c] = find(che);
            if (length(c)>=var(i))
                indx = randperm(length(c));
                %c(indx)=c;
                c=c(indx);
                link(i,1:var(i)) = c(1:var(i));
                che(c(1:var(i))) = che(c(1:var(i)))-1;
                done(i)=0;
            end
        end
        if (sum(done,1)==0)
            true=0;%if everything is O.K.!!!!!
        else
            count=count+1;
            str='Try for ';
            str=strcat(str,num2str(count));
            indx = randperm(length(variable));
            variable(indx) = variable;
            link = zeros(length(variable), len);
%             disp(str);
        end
    end
    G = zeros(length(variable),length(checkk));
    for i=1:length(variable)
        [~,c]=find(link(i,:));
        G(i,link(i,c))=1;
    end
    indx = randperm(length(variable));
    G(indx,:) = G;
    indx = randperm(length(checkk));
    G(:,indx) = G;
    G=G';
end


function [var,Edge] = getNoNodeDegree(Ru, n)

    [r,c] = find(Ru);
    A= zeros(length(c),length(c));
    for i=1:length(c)
        A(i,:) = Ru(c(i)); 
    end
    I = eye(length(c));
    A= A-I;
    lastRow = zeros(1,length(c));
    for i=1:length(c)
        lastRow(i) = 1/c(i);
    end
    A = [A;lastRow];
    B=zeros(length(c)+1,1);
    B(end) = n;
    X = linsolve(A,B);
    for i=1:length(X)
        X2(i) = (X(i)/c(i));
    end
    X2 = round(X2);% number of variables with specified degrees!
    X2(end) = n - sum(X2(1:end-1),2);
    Edge = sum(X2.*c);
    degDist= zeros(1,length(Ru));
    degDist(c) = (X2.*c)/Edge;
    var = zeros(1,n);
    j=1;
    for i=1:length(X2)
        var(j:j+X2(i)-1)= c(i);
        j= j+X2(i);
    end

end