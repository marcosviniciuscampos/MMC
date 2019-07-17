function [X z  I iter tiposolucao] = simplex3(A,b,c, L, U)

%busca a quantidade de linhas m e colunas n
[m n] = size(A);

%cria identidade m por m
Identidade = eye(m);

%completa os limites
for i = 1 : m,
    L = [L;0];
    U = [U;inf];
end

%cria a matriz A Ficticia
AF = [A Identidade];



[mF nF] = size(AF);





%gera vetor I inicial
I = [];
for i = (nF-mF)+1:nF,    
    I = [I i];
end


%gera vetor J1 inicial
J1 = [];
for i = 1:nF,
    existe = 0
    for j = 1:mF
        if I(j) == i
            existe = 1;  
            break;
        end
    end
    if (existe == 0 && L(i) ~= -Inf)
        J1 = [J1, i];
    end       
end

[nada,tamanho_J1] = size(J1);

%cálculo xJ1
xJ1 = [];
for i=1: tamanho_J1,
    xJ1 = [xJ1;L(J1(i))]; 
end


%gera vetor J2 inicial
J2 = [];
for i = 1:nF-mF,
    existe = 0;
    for j = 1:tamanho_J1
        if J1(j) == i
            existe = 1;  
            break;
        end
    end
    if (existe == 0)
        J2 = [J2, i];
    end       
end

[nada,tamanho_J2] = size(J2);

xJ2 = [];
for i = 1:tamanho_J2,
    if (U(J2(i)) == inf)
        xJ2 = [xJ2;0];
    else
        xJ2 = [xJ2;U(J2(i))]; 
    end;
end
%xJ2 = [0];

c = c';

%guarda o valor de c do problema original
cOriginal = c;

%gera o valor de problema fictício
c = [];
for i=1 : nF,
    if (i<=nF-mF)
        c = [c 0];
    else
        c = [c 1];
    end
end

%Enquanto não chegar a solução ótima, prossiga
termina = 0;
iter = 0;
xI = b;
while (termina == 0)
    if (rank([AF b]) > rank(AF))
        tiposolucao = -1;
        break;
    end
    
    %obtem a quantidade de variaveis em J1 e J2
    [nada,tamanho_J1] = size(J1);
    [nada,tamanho_J2] = size(J2);
    
    iter = iter + 1;
    
    %preenche a matriz I e o vetor CI
    AI = [];
    cI = [];
    for cont = 1:mF,    
        AI = [AI AF(:,I(cont))];
        cI = [cI c(I(cont))];
    end

    %preenche a matriz J1 e o vetor cJ1
    AJ1 = []; 
    cJ1 = [];
   % xJ1 = [];
    for cont = 1:tamanho_J1,
        AJ1 = [AJ1 AF(:,J1(cont))];
        cJ1 = [cJ1 c(J1(cont))];  
        %xJ1 = [xJ1;L(J1(cont))];
    end
    
    %preenche a matriz J2 e o vetor cJ2
    AJ2 = []; 
    cJ2 = [];
   % xJ2 = [];
    for cont = 1:tamanho_J2,
        AJ2 = [AJ2 AF(:,J2(cont))];
        cJ2 = [cJ2 c(J2(cont))];    
     %   xJ2 = [xJ2;U(J2(cont))];
    end
    
    
    
    %calculo pi
    pi = cI * inv(AI);
    
    
    
    %cálculo xI    
    if (tamanho_J1 == 0)
        xI = (inv(AI) * b) - (inv(AI)*AJ2*xJ2);
        q = cI*xI + cJ2*xJ2; 
    elseif (tamanho_J2== 0)
        xI = (inv(AI) * b) - (inv(AI)*AJ1*xJ1);
        q = cI*xI + cJ1*xJ1;
         q = pi*xI;
    else
        xI = (inv(AI) * b) - (inv(AI)*AJ1*xJ1) - (inv(AI)*AJ2*xJ2);
        q = cI*xI + cJ1*xJ1 + cJ2*xJ2; 
    end;

    % Encontrando q==0, verifica-se se existe variável artificial na base
    if (q == 0)
        existe_artificial = 0;
        for i = 1:m,
            for j = n+1 : nF,
                if (I(i) == j)
                    existe_artificial = 1;
                    break;
                end
            end
            if (existe_artificial == 1)
                break;
            end
        end
        if (existe_artificial == 0)
            termina = 1;
            continue;
        end
    end
    
        
        
   

    %cálculo do maior cchapeuJ1
    cchapeuJ1 = [];
    kJ1 = 0;
    maiorJ1 = 0;
    if (tamanho_J1 ~= 0)
        cchapeuJ1 = (pi*AJ1) - cJ1;
        [maiorJ1 kJ1] = max(cchapeuJ1);
        if (maiorJ1 <= 0)
            kJ1 = 0;       
        end
        
    end
    
    %cálculo do menor cchapeuJ2
    cchapeuJ2 = [];
    kJ2 = 0;
    menorJ2 = 0;
    if (tamanho_J2 ~= 0)
        cchapeuJ2 = (pi*AJ2) - cJ2;
        [menorJ2 kJ2] = min(cchapeuJ2);
        if (menorJ2 >= 0)
            kJ2 = 0;
        end
        
    end
    
    %significa que J1 e J2 são vazios
    if (kJ1 == 0 && kJ2 == 0)
        %solucao ótima única
        tiposolucao = 0;
        termina = 1;
        continue;
    else
        if (maiorJ1 >= (-1*menorJ2))
            k = J1(kJ1);
        else
            k = J2(kJ2);
        end
    end
    %Variável k entra na base
    
    %variável que entra vem de J1
    if (kJ1 ~= 0 && k == J1(kJ1))
                
        yk = inv(AI)*AF(:,J1(kJ1));
        
        %1º caso        
        %verifica se existe algum elemento positivo em yk. Caso não, gama 1
        %recebe infinito
        if (max(yk) <= 0)
            gama1 = inf;
        else
            r1 = 0;
            %havendo valor > 0, encontra-se a menor divisão
            for i = 1:mF
                if (yk(i) > 0)
                    aux = (xI(i) - L(I(i))) / yk(i);                    
                    if (r1 == 0)
                        r1 = I(i);
                        menor = aux;
                    else
                        if (aux < menor)
                            menor = aux;
                            r1 = I(i);
                        end
                    end
                end
            end
            gama1 = menor;
        end
            
        
        %2º caso
        %verifica se existe algum elemento negativo em yk. Caso não, gama 2
        %recebe infinito
        if (min(yk) >= 0)
            gama2 = inf;
        else            
            r2 = 0;
            %havendo valor < 0, encontra-se a menor divisão
            for i = 1:mF
                if (yk(i) < 0)
                    aux = ( U(I(i)) - xI(i) ) / (-1*yk(i));                    
                    if (r2 == 0)
                        r2 = I(i);
                        menor = aux;
                    else
                        if (aux < menor)
                            menor = aux;
                            r2 = I(i);
                        end
                    end
                end
            end
            gama2 = menor;
        end
        
        
        %3º caso
        gama3 = U(k) - L(k);
        r3 = k;
        
        [deltaK pos] = min([gama1 gama2 gama3]);
       
        switch pos
            case 1 
                r = r1
            case 2 
                r = r2
            case 3 
                r = r3
            otherwise
                disp('problema')
        end;
        
        if (deltaK == inf)
            %solução ilimitada
            termina = 1;
            tiposolucao = 2;
            continue;
        else
            xJ1(kJ1) = L(k) + deltaK;        
            xI = xI-(yk*deltaK);
            if (r ~= k)                 
                for i = 1:mF,
                    if I(i) == r
                        aux = i;
                    end
                end
                aux_val = xI(aux);
                xI(aux) = xJ1(kJ1);
                xJ1(kJ1) = aux_val;
                I(aux) = k;
                if (U(r)== xJ1(kJ1))
                    J2 = [J2 r];
                    J2(tamanho_J2+1) = r;                
                    xJ2(tamanho_J2+1) = U(r);
                else
                    J1(kJ1) = r;  
                end
            else
                aux = J1(kJ1);
                J1(kJ1) = [];
                aux_val = xJ1(kJ1);
                xJ1(kJ1) = [];
                J2 = [J2 aux];
                %xJ2(tamanho_J2+1) = aux_val;
                xJ2 = [xJ2;aux_val];
            end;
        end; 
    %quando k vem de J2
    elseif (k == J2(kJ2))
		yk = inv(AI)*AF(:,J2(kJ2));
		
        %1º caso
        %verifica se existe algum elemento negativo em yk. Caso não, gama 2
        %recebe infinito
        if (min(yk) >= 0)
            gama1 = inf;
        else            
            r1 = 0;
            %havendo valor < 0, encontra-se a menor divisão
            for i = 1:mF
                if (yk(i) < 0)
                    aux = ( xI(i) - L(I(i)) ) / (-1*yk(i));                    
                    if (r1 == 0)
                        r1 = I(i);
                        menor = aux;
                    else
                        if (aux < menor)
                            menor = aux;
                            r1 = I(i);
                        end
                    end
                end
            end
            gama1 = menor;
        end
		
		%2º caso
		if (max(yk) <= 0)
            gama2 = inf;
        else            
            r2 = 0;            
            for i = 1:mF
                if (yk(i) > 0)
                    aux = ( U(I(i)) - xI(i)) / (yk(i));                    
                    if (r2 == 0)
                        r2 = I(i);
                        menor = aux;
                    else
                        if (aux < menor)
                            menor = aux;
                            r2 = I(i);
                        end
                    end
                end
            end
            gama2 = menor;
        end
		
		%3º caso
		gama3 = U(k) - L(k);
        r3 = k;
		
		[deltaK pos] = min([gama1 gama2 gama3]);
       
        switch pos
            case 1 
                r = r1
            case 2 
                r = r2
            case 3 
                r = r3
            otherwise
                disp('problema')
        end;
        
        if (deltaK == inf)
            %solução ilimitada
            tiposolucao = 2;
            termina = 1;
            continue;
        else
            xJ2(kJ2) = U(k) - deltaK;    
            xI = xI+(yk*deltaK);
            
            if (r ~= k)
                for i = 1:mF,
                    if I(i) == r
                        aux = i;
                    end
                end
                aux_val = xI(aux);
                xI(aux) = xJ2(kJ2);
                xJ2(kJ2) = aux_val;
                I(aux) = k;   
                if (L(r)== xJ2(kJ2))
                    J1 = [J1 r];
                    J1(tamanho_J1+1) = r;                
                    xJ1(tamanho_J1+1) = L(r);
                else
                    J2(kJ2) = r;
                end
               
            else
                
                aux = J2(kJ2);
                J2(kJ2) = [];
                aux_val = xJ2(kJ2);
                xJ2(kJ2) = [];
                J1 = [J1 aux];
                xJ1(tamanho_J1+1) = aux_val;                
                %xJ2 = [xJ2;aux_val];
                xJ1 = [xJ1;aux_val];
            end;
        end; 
		
    end
    
         
end

%siginifca que não existe função factível
if (q > 0)
    tiposolucao = -1;
    disp('não tem solucao');
    return;
end
c = cOriginal;

%verifica se as variáveis artificiais estão na base retornada
linhas_excedentes = [];
for i = nF-mF+1:nF,
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
    for (i = nF-mF+1:nF) 
        excedente = 0;
        for (j = 1:n2),
            if (i == linhas_excedentes(j))
                excedente = 1;     
                I(i-mF) = [];
                b(i-mF) = [];                
                break;
            end
        end
        if (excedente == 0)
            aux = [aux;A(i-mF,:)];
        end
    end
    A = aux;
end


%FASE 2

%elimina variáveis artificiais
[nada,tamanho_J2] = size(J2);
[nada,tamanho_J1] = size(J1);
for (i = nF-mF+1:nF)
    j = 1;
    while (j <= tamanho_J1)
        if (i == J1(j))
            J1(j) = [];
            xJ1(j) = [];
            tamanho_J1 = tamanho_J1 - 1;           
        else
            j = j + 1;
        end
    end;
        
    j = 1;
    while (j <= tamanho_J2)
        if (i == J2(j))
            J2(j) = [];
            xJ2(j) = [];
            tamanho_J2 = tamanho_J2 - 1;            
        else
            j = j + 1;
        end;
    end;   
end;

%gera vetor J1 inicial
%J1 = [];
%for i = 1:n,
 %   existe = 0
  %  for j = 1:m
   %     if I(j) == i
    %        existe = 1;  
     %       break;
      %  end
    %end
    %if existe == 0
    %    J1 = [J1, i];
    %end       
%end



%cálculo xJ1
xJ1 = [];
for i=1: tamanho_J1,
    xJ1 = [xJ1;L(J1(i))]; 
end

%cálculo xJ2
xJ2 = [];
for i=1: tamanho_J2,
    xJ2 = [xJ2;U(J2(i))]; 
end


%vetor J2 começa vazio
%J2 = [];
%xJ2 = [];


%Enquanto não chegar a solução ótima, prossiga
termina = 0;
while (termina == 0)    
    %obtem a quantidade de variaveis em J1 e J2
    [nada,tamanho_J1] = size(J1);
    [nada,tamanho_J2] = size(J2);
    
    iter = iter + 1;
    
    %preenche a matriz I e o vetor CI
    AI = [];
    cI = [];
    for cont = 1:m,    
        AI = [AI A(:,I(cont))];
        cI = [cI c(I(cont))];
    end

    %preenche a matriz J1 e o vetor cJ1
    AJ1 = []; 
    cJ1 = [];
    %xJ1 = [];
    for cont = 1:tamanho_J1,
        AJ1 = [AJ1 A(:,J1(cont))];
        cJ1 = [cJ1 c(J1(cont))];  
        %xJ1 = [xJ1;L(J1(cont))];
    end
    
    %preenche a matriz J2 e o vetor cJ2
    AJ2 = []; 
    cJ2 = [];
    %xJ2 = [];
    for cont = 1:tamanho_J2,
        AJ2 = [AJ2 A(:,J2(cont))];
        cJ2 = [cJ2 c(J2(cont))];    
        %xJ2 = [xJ2;U(J2(cont))];
    end
    
    
    
    %calculo pi
    pi = cI * inv(AI);
    
    %cálculo xI    
    if (tamanho_J1 == 0)
        xI = (inv(AI) * b) - (inv(AI)*AJ2*xJ2);
        z = cI*xI + cJ2*xJ2; 
    elseif (tamanho_J2== 0)
        xI = (inv(AI) * b) - (inv(AI)*AJ1*xJ1);
        z = cI*xI + cJ1*xJ1;
    else
        xI = (inv(AI) * b) - (inv(AI)*AJ1*xJ1) - (inv(AI)*AJ2*xJ2);
        z = cI*xI + cJ1*xJ1 + cJ2*xJ2; 
    end;


   

    %cálculo do maior cchapeuJ1
    cchapeuJ1 = [];
    kJ1 = 0;
    maiorJ1 = 0;
    if (tamanho_J1 ~= 0)
        cchapeuJ1 = (pi*AJ1) - cJ1;
        [maiorJ1 kJ1] = max(cchapeuJ1);
        if (maiorJ1 <= 0)
            kJ1 = 0;       
        end
        
    end
    
    %cálculo do menor cchapeuJ2
    cchapeuJ2 = [];
    kJ2 = 0;
    menorJ2 = 0;
    if (tamanho_J2 ~= 0)
        cchapeuJ2 = (pi*AJ2) - cJ2;
        [menorJ2 kJ2] = min(cchapeuJ2);
        if (menorJ2 >= 0)
            kJ2 = 0;
        end
        
    end
    
    %significa que J1 e J2 são vazios
    if (kJ1 == 0 && kJ2 == 0)
        %solucao ótima única
        tiposolucao = 0;
        termina = 1;
        continue;
    else
        if (maiorJ1 >= (-1*menorJ2))
            k = J1(kJ1);
        else
            k = J2(kJ2);
        end
    end
    %Variável k entra na base
    
    %variável que entra vem de J1
    if (kJ1 ~= 0 && k == J1(kJ1))
                
        yk = inv(AI)*A(:,J1(kJ1));
        
        %1º caso        
        %verifica se existe algum elemento positivo em yk. Caso não, gama 1
        %recebe infinito
        if (max(yk) <= 0)
            gama1 = inf;
        else
            r1 = 0;
            %havendo valor > 0, encontra-se a menor divisão
            for i = 1:m
                if (yk(i) > 0)
                    aux = (xI(i) - L(I(i))) / yk(i);                    
                    if (r1 == 0)
                        r1 = I(i);
                        menor = aux;
                    else
                        if (aux < menor)
                            menor = aux;
                            r1 = I(i);
                        end
                    end
                end
            end
            gama1 = menor;
        end
            
        
        %2º caso
        %verifica se existe algum elemento negativo em yk. Caso não, gama 2
        %recebe infinito
        if (min(yk) >= 0)
            gama2 = inf;
        else            
            r2 = 0;
            %havendo valor < 0, encontra-se a menor divisão
            for i = 1:m
                if (yk(i) < 0)
                    aux = ( U(I(i)) - xI(i) ) / (-1*yk(i));                    
                    if (r2 == 0)
                        r2 = I(i);
                        menor = aux;
                    else
                        if (aux < menor)
                            menor = aux;
                            r2 = I(i);
                        end
                    end
                end
            end
            gama2 = menor;
        end
        
        
        %3º caso
        gama3 = U(k) - L(k);
        r3 = k;
        
        [deltaK pos] = min([gama1 gama2 gama3]);
       
        switch pos
            case 1 
                r = r1
            case 2 
                r = r2
            case 3 
                r = r3
            otherwise
                disp('problema')
        end;
        
        if (deltaK == inf)
            %solução ilimitada
            termina = 1;
            tiposolucao = 2;
            continue;
        else
            xJ1(kJ1) = L(k) + deltaK;        
            xI = xI-(yk*deltaK);
            if (r ~= k)                 
                for i = 1:mF,
                    if I(i) == r
                        aux = i;
                    end
                end
                aux_val = xI(aux);
                xI(aux) = xJ1(kJ1);
                xJ1(kJ1) = aux_val;
                I(aux) = k;
                %if (U(r)== xJ1(kJ1))
                 %   J2 = [J2 r];
                  %  J2(tamanho_J2+1) = r;                
                    %xJ2(tamanho_J2+1) = U(r);
                   % xJ2 = [xJ2;U(r)];
                %else
                    J1(kJ1) = r;  
                %end
            else
                aux = J1(kJ1);
                J1(kJ1) = [];
                aux_val = xJ1(kJ1);
                xJ1(kJ1) = [];
                J2 = [J2 aux];
                %xJ2(tamanho_J2+1) = aux_val;
                xJ2 = [xJ2;aux_val];
            end;
        end; 
    %quando k vem de J2
    elseif (k == J2(kJ2))
		yk = inv(AI)*A(:,J2(kJ2));
		
        %1º caso
        %verifica se existe algum elemento negativo em yk. Caso não, gama 2
        %recebe infinito
        if (min(yk) >= 0)
            gama1 = inf;
        else            
            r1 = 0;
            %havendo valor < 0, encontra-se a menor divisão
            for i = 1:m
                if (yk(i) < 0)
                    aux = ( xI(i) - L(I(i)) ) / (-1*yk(i));                    
                    if (r1 == 0)
                        r1 = I(i);
                        menor = aux;
                    else
                        if (aux < menor)
                            menor = aux;
                            r1 = I(i);
                        end
                    end
                end
            end
            gama1 = menor;
        end
		
		%2º caso
		if (max(yk) <= 0)
            gama2 = inf;
        else            
            r2 = 0;            
            for i = 1:m
                if (yk(i) > 0)
                    aux = ( U(I(i)) - xI(i)) / (yk(i));                    
                    if (r2 == 0)
                        r2 = I(i);
                        menor = aux;
                    else
                        if (aux < menor)
                            menor = aux;
                            r2 = I(i);
                        end
                    end
                end
            end
            gama2 = menor;
        end
		
		%3º caso
		gama3 = U(k) - L(k);
        r3 = k;
		
		[deltaK pos] = min([gama1 gama2 gama3]);
       
        switch pos
            case 1 
                r = r1
            case 2 
                r = r2
            case 3 
                r = r3
            otherwise
                disp('problema')
        end;
        
        if (deltaK == inf)
            %solução ilimitada
            tiposolucao = 2;
            termina = 1;
            continue;
        else            
            xJ2(kJ2) = U(k) - deltaK;    
            xI = xI+(yk*deltaK);
            
            if (r ~= k)
                for i = 1:mF,
                    if I(i) == r
                        aux = i;
                    end
                end
                aux_val = xI(aux);
                xI(aux) = xJ2(kJ2);
                xJ2(kJ2) = aux_val;
                I(aux) = k;   
                %if (L(r)== xJ2(kJ2))
                 %   J1 = [J1 r];
                  %  J1(tamanho_J1+1) = r;                
                  %  xJ1(tamanho_J1+1) = L(r);
                %else
                    J2(kJ2) = r;
                %end
               
            else
                aux = J2(kJ2);
                J2(kJ2) = [];
                aux_val = xJ2(kJ2);
                xJ2(kJ2) = [];
                J1 = [J1 aux];
                %xJ1(tamanho_J1+1) = aux_val;
                xJ1 = [xJ1;aux_val];
            end;        
        end; 
		
    end
    
         
end

%montar X
X = [];
for i = 1: m,
    X(I(i)) = xI(i);    
end
[nada,tamanho_J1] = size(J1);
[nada,tamanho_J2] = size(J2);
for i = 1: tamanho_J1,
    X(J1(i)) = xJ1(i);
end
for i = 1: tamanho_J2,
    X(J2(i)) = xJ2(i);
end

