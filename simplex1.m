function [X z  I iter tiposolucao] = simplex1(A,b,c,I)


%busca a quantidade de linhas m e colunas n
[m n] = size(A);


%gera vetor J inicial
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
       
        
        %calcula o máximo de cchapeuJ e o respectivo indice k
        [maior k] = max(cchapeuJ);

        %se cchapeuJ for menor que 0, chegamos a solução ótima
        if (maior < 0)        
            termina = 1;
            tiposolucao = 0;
            break;    
        elseif (maior == 0)        
            termina = 1;
            tiposolucao = 1;
            break;
        end

        yk = inv(AI)*A(:,J(k));

        if (max(yk) <= 0)
            termina = 1;
            tiposolucao = 2;
            continue;
        end
        
        menor = 0;
        for i=1:m,    
            if yk(i) > 0                
                aux = xI(i)/yk(i);                
                if menor == 0
                    menor = aux;
                    r = i;
                elseif (aux < menor)
                    menor = aux;            
                    r = i;
                end               
            end
        end
        
        %coloca k no lugar de r em I, e o r no lugar de k em J
        aux = I(r);
        I(r) = J(k);
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

