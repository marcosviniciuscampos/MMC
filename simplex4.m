function [X z  I iter tiposolucao] = simplex4(A,b,c,I)

%busca a quantidade de linhas m e colunas n
[m n] = size(A);

%gera vetor J 
J = [];
for i = 1:n,
    existe = 0;
    for j = 1:m
        if I(j) == i
            existe = 1;  
            break;
        end
    end
    if existe == 0
        J = [J, i];
    end       
end

c = c';

%Enquanto não chegar a solução ótima, prossiga
termina = 0;
iter = 0;
if (rank([A b]) > rank(A))
        tiposolucao = -1;        
else
    while (termina == 0)
    
    
        iter = iter + 1;
    
        %preenche a matriz I e o vetor CI
        AI = [];
        cI = [];
        for cont = 1:m,    
            AI = [AI A(:,I(cont))];
            cI = [cI c(I(cont))];
        end

        %preenche a matriz J e o vetor cJ
        AJ = []; 
        cJ = [];
        for cont = 1:n-m,
            AJ = [AJ A(:,J(cont))];
            cJ = [cJ c(J(cont))];        
        end

        %cálculo xI
        xI = inv(AI) * b;
        
        %cálculo xJ
        xJ = [];
        for i=1: n-m
            xJ = [xJ 0]; 
        end
    
        %cálculo pí
        pi = cI * inv(AI);

        %cálculo função objetivo
        z = pi*b; 

        %cálculo cchapeuJ
        cchapeuJ = (pi*AJ) - cJ;
       
        inv_AIAJ = inv(AI)*AJ;
        
        %Escolhe quem sai
        [menor r] = min(xI);
        
        if (menor >= 0)
            disp('solucao otima determinada encontrada');
            termina = 1;
            tiposolucao = 0;
            break
        end      
        idx_r = r;
        r = I(r);
        
        
        menor = inf;
        for i=1:n-m,
            if (inv_AIAJ(idx_r,i) < 0)
                valor = cchapeuJ(i) / inv_AIAJ(idx_r,i);
                if (valor < menor)
                    menor = valor;
                    k = i;
                end;
            end;
        end;
        
        if (menor == inf)
            termina = 1;
            tiposolucao = -1;
            continue;
        end;
       
        %coloca k no lugar de r em I, e o r no lugar de k em J
        aux = I(idx_r);
        I(idx_r) = J(k);
        J(k) = aux;
    end
end
%montar X
X = [];
for i = 1: m,
    X(I(i)) = xI(i);    
end
for i = 1: n-m,
    X(J(i)) = xJ(i);
end

 
end


