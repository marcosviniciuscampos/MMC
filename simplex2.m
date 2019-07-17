function [X z  I iter tiposolucao] = simplex2(A,b,c)


%busca a quantidade de linhas m e colunas n
[m n] = size(A);

%cria identidade m por m
Identidade = eye(m);


%cria a matriz A Ficticia
AF = [A Identidade];

[m n] = size(AF);

%gera vetor I inicial
I = [];
for i = (n-m)+1:n,    
    I = [I i];
end

%gera vetor J inicial
J = [];
for i = 1:n,
    existe = 0
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

%guarda o valor de c do problema original
cOriginal = c;

%gera o valor de problema fictício
c = [];
for i=1 : n,
    if (i<=n-m)
        c = [c 0];
    else
        c = [c 1];
    end
end


%Resolve o problema fictício
termina = 0;
iter = 0;
disp(AF);
while (termina == 0)    
    iter = iter + 1;
    
    %preenche a matriz I e o vetor CI
    AI = [];
    cI = [];
    for cont = 1:m,    
        AI = [AI AF(:,I(cont))];
        cI = [cI c(I(cont))];
    end

    %preenche a matriz J e o vetor cJ
    AJ = []; 
    cJ = [];
    for cont = 1:n-m,
        AJ = [AJ AF(:,J(cont))];
        cJ = [cJ c(J(cont))];        
    end    
    
    %cálculo xI
    xI = inv(AI) * b;

    disp('xI');
    disp(xI);
    
    %cálculo xJ
    xJ = [];
    for i=1: n-m
        xJ = [xJ 0]; 
    end
    
    %cálculo pí
    pi = cI * inv(AI);

    %cálculo função objetivo
    q = pi*b;             
        
    var_artificial_base = 0;
    if (q == 0)
        for i = 1:m,
            for j = n-m+1:n,
                %verifica se existe variavel artificial na base
                if (I(i) == j)
                    var_artificial_base = 1;
                    break;
                end
            end
        end
        if (var_artificial_base == 0)
            break;
        end
    end
    
    %cálculo cchapeuJ
    cchapeuJ = pi*AJ - cJ;    
    
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
    
    %se cchapeuJ for menor que 0, chegamos a solução ótima
    if (maior < 0)        
        termina = 1;  
        continue;
    end
    
    negativo = 0;
    for i=1:m,
        if (xI(i) < 0)
            negativo = 1;
        end
    end
    if (maior == 0 && negativo == 0)
        termina = 1;
        continue;
    end
        
    
    yk = inv(AI)*AF(:,J(k));

    
    menor = 0;
    for i=1:m,    
        if yk(i) > 0
            termina = 0;
            aux = xI(i)/yk(i);
            if (aux > 0 || q < 0)
                if menor == 0
                    menor = aux;
                    r = i;
                elseif (aux < menor)
                    menor = aux;            
                    r = i;
                end
            end
        end
    end
    
    if (menor == 0)        
        termina = 1;        
    end

    %coloca k no lugar de r em I, e o r no lugar de k em J
    aux = I(r);
    I(r) = J(k);
    J(k) = aux;
    
    
end

%siginifca que não existe função factível
if (q > 0)
    tiposolucao = -1;
    X = [];
    z = [];
    disp('não tem solucao');
    return;
end
c = cOriginal;

%verifica se as variáveis artificiais estão na base retornada
linhas_excedentes = [];
for i = n-m+1:n,
   for j = 1: m,
       %A variável artificial i é uma das que fazem parte da base?
       if (I(j) == i)
           linhas_excedentes = [linhas_excedentes i];
           break;
       end
   end                
end    

[m2 n2] = size(linhas_excedentes);
%Caso estejam na base, são eliminadas as linhas correspondentes da matriz principal
if (n2 > 0)
    aux = [];
    for (i = n-m+1:n) 
        excedente = 0;
        for (j = 1:n2),
            if (i == linhas_excedentes(j))
                excedente = 1;     
                I(i-m) = [];
                b(i-m) = [];                
                break;
            end
        end
        if (excedente == 0)
            aux = [aux;A(i-m,:)];
        end
    end
    A = aux;
end

c = c';
%exclui variáveis artificiais
i = 1;
tamanhoJ = n-m;
while i < tamanhoJ 
    tamanhoAnt = tamanhoJ;
    for j = n-m+1:n,
        %verifica se existe variavel artificial na base
        if (J(i) == j)
            J(i) = [];
            tamanhoJ = tamanhoJ - 1;
            break;
        end
    end
    if tamanhoAnt == tamanhoJ
        i = i+1;
    end
end

%busca a quantidade de linhas m e colunas n
[m n] = size(A);

%Enquanto não chegar a solução ótima, prossiga
termina = 0;

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

%montar X
X = [];
for i = 1: m,
X(I(i)) = xI(i);    
end
for i = 1: n-m,
X(J(i)) = xJ(i);
end



disp('fim');