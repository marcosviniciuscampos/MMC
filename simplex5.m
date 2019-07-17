function [X z  I iter tiposolucao] = simplex5(A,b,c,PI_barra)

%busca a quantidade de linhas m e colunas n
[m n] = size(A);

c = c';
iter = 0;
termina = 0;
while (termina == 0)
    iter = iter + 1;
    %monta o vetor Q
    cchapeu = PI_barra*A-c;
    Q = [];
    for i = 1:n,
        if cchapeu(i) == 0
            Q = [Q i];
        end;        
    end
    
    [nada n_Q] = size(Q);
    
    
    %monta problema primal restrito
    c_primal = [];
    A_primal = [];       
        
    for i=1:n_Q,
        A_primal = [A_primal A(:,Q(i))];
        
    end;
    A_primal = [A_primal zeros(m-n_Q,m)'];
    for i = 1:m,
        c_primal = [c_primal 0];
    end;
    for i = 1:m,
        c_primal = [c_primal 1];
    end;    
    A_primal = [A_primal eye(m)];
    
    [m_primal n_primal] = size(A_primal);
    I = [];
    for i = 1:m_primal,
       I = [I n_primal - m_primal + i];       
    end;
        
        
        
        
    %resolve problema primal restrito
    %gera vetor J inicial
    J = [];
    for i = 1:n_primal,
        existe = 0;
        for j = 1:m_primal
            if I(j) == i
                existe = 1;  
                break;
            end
        end
        if existe == 0
            J = [J, i];
        end       
    end
    termina_primal = 0;
    while (termina_primal == 0)
        %preenche a matriz I e o vetor CI
        A_primal_I = [];
        c_primal_I = [];
        for cont = 1:m_primal,    
            A_primal_I = [A_primal_I A_primal(:,I(cont))];
            c_primal_I = [c_primal_I c_primal(I(cont))];
        end
        %preenche a matriz J e o vetor cJ
        A_primal_J = []; 
        c_primal_J = [];
        for cont = 1:n_primal-m_primal,
            A_primal_J = [A_primal_J A_primal(:,J(cont))];
            c_primal_J = [c_primal_J c_primal(J(cont))];        
        end
        %cálculo xI        
        xI = inv(A_primal_I) * b;
        %cálculo xJ
        xJ = [];
        for i=1: n_primal-m_primal
            xJ = [xJ 0]; 
        end
        %cálculo pí
        pi = c_primal_I * inv(A_primal_I);
        %cálculo função objetivo
        z_zero = pi*b; 
        %cálculo cchapeuJ
        if (n_primal-m_primal == 0)
            termina = 1;
            break;
        end;
    
        cchapeuJ = (pi*A_primal_J) - c_primal_J;
        %calcula o máximo de cchapeuJ e o respectivo indice k
        [maior k] = max(cchapeuJ);
        %se cchapeuJ for menor que 0, chegamos a solução ótima
        if (maior < 0)        
            termina_primal = 1;
            tiposolucao = 0;
            break;    
        elseif (maior == 0)        
            termina_primal = 1;
            tiposolucao = 1;
            break;
        end
        yk = inv(A_primal_I)*A_primal(:,J(k));

        if (max(yk) <= 0)
            termina_primal = 1;            
            continue;
        end        
        menor = 0;
        for i=1:m_primal,    
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
    X_primal = [];
    for i = 1: m,
        X_primal(I(i)) = xI(i);    
    end
    for i = 1: n_primal-m,
        X_primal(J(i)) = xJ(i);
    end

    
    % -----------********************--------TERMINA RESOLUÇÃO PROBLEMA
    % PRIMAL
                               
    
    if (z_zero == 0)
        termina = 1;               
    else                
        A_dual = A_primal';
        A_dual = [A_dual eye(n_primal)];
       
        
        
        % resolve dual através de folgas complementares
        [nada tamanho_X_primal] = size(X_primal);
        for i=1:m_primal+n_primal-tamanho_X_primal,            
            X_primal = [X_primal 0];
        end;        
        X_dual = [];        
        for i=1:n_primal,
            if (X_primal(i)>0)
                X_dual(i+m_primal) = 0;
            else
                X_dual(i+m_primal) = inf;
            end;
        end;
        for i=1:m_primal,
            if (X_primal(n_primal+i)>0)
                X_dual(i) = 0;
            else
                X_dual(i) = inf;
            end;
        end;                       
        A_dual_I = [];
        for i=1:m_primal+n_primal,
            if (X_dual(i) > 0)
                A_dual_I = [A_dual_I A_dual(:,i)]; 
            end;
        end;        
        %Resolve
        v = inv(A_dual_I)*c_primal';        
        v_star = [];
        for i=1: m,
            v_star = [v_star v(i)];
        end;


        
        
        
        
        %verifica se é primal infactível
        for i=1:n,
            primal_infactivel = 1;
            if (v_star*A(:,i) > 0)   
                primal_infactivel = 0;
                break;                
            end;
            if (primal_infactivel == 1)
                tipo_solucao = -1;
                termina = 1;
                break;
            end
        end
        
        if (primal_infactivel == 1)
            break;
        end;
        
        %calcula teta
        maior = 0;
        teta = [];
        for i=1:n,
            pertence = 0;
            for (j =1:n_Q)
                if (Q(j) == i)
                    pertence = 1;
                    break;
                end;
            end;
            if (pertence == 0)
                if (v_star * A(:,i) > 0)
                    teta = [teta  -1*(PI_barra*A(:,i) - c(i)) / (v_star * A(:,i))];
                end;
            end;
        end;
        teta = min(teta);
        PI_barra = PI_barra + (teta * v_star);
    end
    
   
end

%montar X
X = [];
for i = 1: m,
    X(I(i)) = xI(i);    
end

%gera vetor J do problema primal
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

for i =1: n-m,
    X(J(i)) = 0;
end

z = PI_barra *b;

end