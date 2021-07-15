clc
clear all
% MONORRESTITUIÇÃO PELAS EQUAÇÕES DE COLINEARIDADE

%% 1)DEFINIÇÕES INICIAIS 
%Altura aproximada de voo (metros)
hm=300;
%Distância Focal (mm)
f=15.8893;
%Tamanho do píxel(milímetros)
pixel=0.007;
%Número de colunas
nc = 3344;
%Número de linhas
nl = 2224;
%Coordenadas do Ponto Principal (milímetros)
xpp = 0.1746;
ypp = -0.0947;
%Numero de pontos de apoio
n=6;

%% 2)ESPAÇO IMAGEM

%Leitura de coordenadas dos pontos na imagem. Estas medidas (observações)
%devem ser efetuadas através de algum software (ex. multispec)
% Vetores coluna e linha referentes aos pontos de apoio (1:6) e de verificação(7:).

pts_digital = [
                129  317  %57
                2686 1900 %84
                1338 1477 %94
                1075 89   %108
                3057 503  %127
                85   1111 %138
                865  1048 %02 
                2205 456  %123 
                2370 1575 %129 
                1816 1432 %133 
                280 1530  %136 
                1692 1199]; %145 
               
               
c=pts_digital(:, 1);
l=pts_digital(:, 2);

%transformação: Ref digital (c,l) para Ref Digital em Milímetros com origem
%no centro da imagem
xmm = pixel*(c(:)-(nc-1)/2);
ymm = -pixel*(l(:)-(nl-1)/2);

%transformação: Ref Dig em Milímetros para o Ref Fotogramétrico com distorções
%Translação para coincidir origem com o Ponto Principal
xpd = xmm - xpp;
ypd = ymm - ypp;

%CORREÇÕES DAS DISTORÇÕES--------------------------------------------------
%Corrigir os efeitos que afastam da condição de colinearidade
%Atualmente são efetuados apenas as correções das distorções de lentes.

%Correção da Distorção Radial Simétrica
%Parâmetros de distorção radial simétrica (Obtidos por processo de calibração da câmera)
%k1 (mm-2)
k1= -2.54848278*10^-4;  
%k1 (mm-4)
k2= 1.48901365*10^-6;   
k3= 0;
%Distância dos pontos ao centro da imagem (radial)
r= sqrt(xpd.^2 + ypd.^2);
%Distorção radial simétrica para a componente x
drsx = (k1*r.^2+k2*r.^4).*xpd;
%Distorção radial simétrica para a componente y
drsy = (k1*r.^2+k2*r.^4).*ypd;


%Correção da Distorção Descentrada
%Parâmetros de distorção descentrada (mm-2)Calculados por processo de
%calibração
p1 = -8.92107195*10^-5;
p2 = -8.00528299*10^-5;
%Distorção descentrada para a componente x
ddcx = p1*(r.^2+2*xpd.^2)+(2*p2*xpd.*ypd);
%Distorção descentrada para a componente y
ddcy = p2*(r.^2+2*ypd.^2)+(2*p1*xpd.*ypd);

%Correção da refreação fotogramétrica
%Não está sendo efetuada.

%Correção das deformações do sensor
%Não está sendo efetuada.

%transformação: Ref Fotogramétrico COM distorções para o Ref Fotogramétrico SEM distorções
%Subtrair os efeitos de distorções do referentcial fotogramétric com distorções.
xp = xpd-drsx-ddcx;
yp = ypd-drsy-ddcy;

%vetor Lb (Vetor das Observações) 6 Pontos de apoio com coordenadas no ref
%fotogramétrico sem distorção para estimar os POE
Lb=[xp(1);yp(1);xp(2);yp(2);xp(3);yp(3);xp(4);yp(4);xp(5);yp(5);xp(6);yp(6)];
%% 2)ESPAÇO OBJETO
%espaço objeto
%Coordenadas Geodésicas dos seis pontos de APOIO  


 XYZ = [
        444.191; 669.765; 21.42;%57
        675.587; 374.946; 1.432;%84
        525.366; 492.427; 7.269;%94
        554.357;  662.66; 12.16;%108
        770.173; 537.567; 8.069;%127
        395.067; 588.307; 8.878];%138
       

%Coordenadas dos pontos de apoio para entrarem nas equações de injunção
XYZ1=  [
        444.191 669.765 21.42;%57
        675.587 374.946 1.432;%84
        525.366 492.427 7.269;%94
        554.357  662.66 12.16;%108
        770.173 537.567 8.069;%127
        395.067 588.307 8.878];%138

% %Coordenadas Geodésicas dos seis pontos de VERIFICAÇÃO 
 XYZV = [488.402; 563.665; 8.267;%02
         669.969; 577.359; 9.741;%123
         648.682; 432.527; 2.117;%129
         585.134; 476.260; 3.493;%133
         395.148; 530.693; 6.973;%136
         579.864; 510.927; 5.448];%145
         

%Precisão dos pontos de apoio (metros)
sXYZ=0.1;
%% AJUSTAMENTO
% Valores Iniciais dos parâmetros (incógnitas) que são os seis parâmetros de
% orientação exterior da imagem (X0, Y0, Z0, omega, phi, kapa)
X0=[500;500;hm;0;0;0;XYZ];

%Matriz Peso das medidas (espaço imagem)
P=eye(12)/(pixel)^2;

%Matriz Peso dos pontos de apoio (espaço objeto)
Pc=eye(18)/sXYZ^2;

%Matriz A - Definição da estrutura da matriz (facilita os processamentos)
A=zeros(12,24);

%Matriz C- Definição da estrutura da matriz
C=zeros(18,24);

itera=0;
fim=false;

while(~fim)
    
    
%Matriz M (Rotações)
%Está relacionada a atitude do sensor através dos parâmetros omega, phi e kappa.
M11=(cos(X0(5))*cos(X0(6)));
M12=(cos(X0(4))*sin(X0(6)))+(sin(X0(4))*sin(X0(5))*cos(X0(6)));
M13=(sin(X0(4))*sin(X0(6)))-(cos(X0(4))*sin(X0(5))*cos(X0(6)));
M21=(-cos(X0(5))*sin(X0(6)));
M22=(cos(X0(4))*cos(X0(6)))-(sin(X0(4))*sin(X0(5))*sin(X0(6)));
M23=(sin(X0(4))*cos(X0(6)))+(cos(X0(4))*sin(X0(5))*sin(X0(6)));
M31=(sin(X0(5)));
M32=(-sin(X0(4))*cos(X0(5)));
M33=(cos(X0(4))*cos(X0(5)));
    
   M=[M11 M12 M13;
      M21 M22 M23;
      M31 M32 M33];
   
 %Modelo Matemático- Equações de colinearidade DIRETA
 %xp=-f*(Qx/Qz)
 %yp=-f*(Qy/Qz)
  
%Enlace com o número de pontos utilizados 
 for i=1:n
     
     %Equações de Colinearidade
     Qx(i)=M11*(X0(6+3*i-2,1)-X0(1))+M12*(X0(6+3*i-1,1)-X0(2))+M13*(X0(6+3*i,1)-X0(3));
     Qy(i)=M21*(X0(6+3*i-2,1)-X0(1))+M22*(X0(6+3*i-1,1)-X0(2))+M23*(X0(6+3*i,1)-X0(3));
     Qz(i)=M31*(X0(6+3*i-2,1)-X0(1))+M32*(X0(6+3*i-1,1)-X0(2))+M33*(X0(6+3*i,1)-X0(3));
     
     %Vetor L0 (f(x0))
     L0(2*i-1,1)=-f*(Qx(i)/Qz(i));
     L0(2*i,1)=-f*(Qy(i)/Qz(i));

     %Derivadas parciais das equações com respeito aos parâmetros irão compor a A)
     %Derivadas em x:
   
     DXX0(i)=-f*(-Qz(i)*M11+Qx(i)*M31)/Qz(i)^2;
     DXY0(i)=-f*(-Qz(i)*M12+Qx(i)*M32)/Qz(i)^2;
     DXZ0(i)=-f*(-Qz(i)*M13+Qx(i)*M33)/Qz(i)^2;
     DXW(i)=-f*(Qz(i)*((X0(6+3*i-1,1)-X0(2))*(-M13)+(X0(6+3*i,1)-X0(3))*M12)-Qx(i)*((X0(6+3*i-1,1)-X0(2))*(-M33)+(X0(6+3*i,1)-X0(3))*M32))/Qz(i)^2;
          
     phi1=(X0(6+3*i-2,1)-X0(1))*(-sin(X0(5)))*cos(X0(6));
     phi2=(X0(6+3*i-1,1)-X0(2))*(sin(X0(4)))*cos(X0(5))*cos(X0(6));
     phi3=(X0(6+3*i,1)-X0(3))*(-cos(X0(4)))*cos(X0(5))*cos(X0(6));
     phi4=(X0(6+3*i-2,1)-X0(1))*cos(X0(5));
     phi5=(X0(6+3*i-1,1)-X0(2))*(sin(X0(4)))*sin(X0(5));
     phi6=(X0(6+3*i,1)-X0(3))*(-cos(X0(4)))*sin(X0(5));
    
     DXFI(i)=-f*(Qz(i)*(phi1+phi2+phi3)-Qx(i)*(phi4+phi5+phi6))/Qz(i)^2;
     DXK(i)=-f*Qy(i)/Qz(i);
   
    %Derivadas em y
     
    DYX0(i)=-f*(-Qz(i)*M21+Qy(i)*M31)/Qz(i)^2;
    DYY0(i)=-f*(-Qz(i)*M22+Qy(i)*M32)/Qz(i)^2;
    DYZ0(i)=-f*(-Qz(i)*M23+Qy(i)*M33)/Qz(i)^2;
    DYW(i)=-f*(Qz(i)*((X0(6+3*i-1,1)-X0(2))*(-M23)+(X0(6+3*i,1)-X0(3))*M22)-Qy(i)*((X0(6+3*i-1,1)-X0(2))*(-M33)+(X0(6+3*i,1)-X0(3))*M32))/Qz(i)^2;
    
    phiy1=(X0(6+3*i-2,1)-X0(1))*(sin(X0(5)))*sin(X0(6));
    phiy2=(X0(6+3*i-1,1)-X0(2))*(-sin(X0(4)))*cos(X0(5))*sin(X0(6));
    phiy3=(X0(6+3*i,1)-X0(3))*(cos(X0(4)))*cos(X0(5))*sin(X0(6));
    phiy4=(X0(6+3*i-2,1)-X0(1))*cos(X0(5));
    phiy5=(X0(6+3*i-1,1)-X0(2))*sin(X0(4))*sin(X0(5));
    phiy6=(X0(6+3*i,1)-X0(3))*(-cos(X0(4)))*sin(X0(5));
    
    DYFI(i)=-f*(Qz(i)*(phiy1+phiy2+phiy3)-Qy(i)*(phi4+phi5+phi6))/Qz(i)^2;
    DYK(i)=f*Qx(i)/Qz(i);
      
   %Preenchimento da Matriz A
   A(2*i-1,1:6)=[DXX0(i) DXY0(i) DXZ0(i) DXW(i) DXFI(i) DXK(i)];
   A(2*i,1:6)=[DYX0(i) DYY0(i) DYZ0(i) DYW(i) DYFI(i) DYK(i)];
  
  
    %INJUNÇÕES [G(Xa)=0]
    %Injunções de posição, portanto serão adicionadas 3*número de pontos de apoio(3 equações por ponto) equações
    %Matriz A - Aumentada com as injunções
    A(2*i-1,6+3*i-2)=-f*(Qz(i)*M11-Qx(i)*M31)/Qz(i)^2;
    A(2*i-1,6+3*i-1)=-f*(Qz(i)*M12-Qx(i)*M32)/Qz(i)^2;
    A(2*i-1,6+3*i)=-f*(Qz(i)*M13-Qx(i)*M33)/Qz(i)^2;
    A(2*i,6+3*i-2)=-f*(Qz(i)*M21-Qy(i)*M31)/Qz(i)^2;
    A(2*i,6+3*i-1)=-f*(Qz(i)*M22-Qy(i)*M32)/Qz(i)^2;
    A(2*i,6+3*i)=-f*(Qz(i)*M23-Qy(i)*M33)/Qz(i)^2;

    %Matriz C (zeros(18x6)eye(18)) das injunções
    C(3*i-2,6+3*i-2)=1;
    C(3*i-1,6+3*i-1)=1;
    C(3*i,6+3*i)=1;

    %Erro de fechamento das injunções
    E(3*i-2,1)=X0(6+3*i-2,1)-XYZ1(i,1);
    E(3*i-1,1)=X0(6+3*i-1,1)-XYZ1(i,2);
    E(3*i,1)=X0(6+3*i,1)-XYZ1(i,3);
 end;
 
 %% AJUSTAMENTO DAS OBSERVAÇÕES
  L=L0-Lb;
  N=A'*P*A;
  U=A'*P*L;
  NC=C'*Pc*C;%Injunção
  Uc=C'*Pc*E; %Injunção
  X=-inv(N+NC)*(U+Uc);
  fim = (max(abs(X(1:6,1))))<1.0E-6 | itera>30;
  Xa=X0+X;
  X0=Xa;
  
end

%Valores ajustados dos POE
Xa;
C;
V=A*X+L;%rESÍDUOS

%Avaliação da qualidade do ajustamento

%Graus de Liberdade (número de equações - número de parâmetros)+ número de
%equações de injunção (3*número de pontos)
GL=12-6+18;
%Variância a posteriori
varpost=(V'*P*V+E'*Pc*E)/(GL);
varXa=varpost*inv(N+NC);
dpXa=sqrt(diag(varXa));
%Teste do Chi^2
qui_cal = varpost*(GL);
% Chi^2 tabelado (Verificado na tabela em função do GL e nível de confiança)
qui1 = 12.401; qui2 = 39.364;
OBS_AJUST=[Xa dpXa];
disp('===================================================================')
disp('RESULTADO APÓS 6 ITERAÇÕES')
disp('===================================================================')

fprintf('Parâmetro de Orientação Exterior X=  %f +- %f \n', OBS_AJUST(1,1), OBS_AJUST(1,2));
fprintf('Parâmetro de Orientação Exterior Y=  %f +- %f \n', OBS_AJUST(2,1), OBS_AJUST(2,2));
fprintf('Parâmetro de Orientação Exterior Z=  %f +- %f \n', OBS_AJUST(3,1), OBS_AJUST(3,2));
fprintf('Parâmetro de Orientação Exterior w=  %f +- %f \n', OBS_AJUST(4,1), OBS_AJUST(4,2));
fprintf('Parâmetro de Orientação Exterior phi=  %f +- %f \n', OBS_AJUST(5,1), OBS_AJUST(5,2));
fprintf('Parâmetro de Orientação Exterior k=  %f +- %f \n', OBS_AJUST(6,1), OBS_AJUST(6,2));
    
for i = 7:24
    fprintf(' Injunções e precisões(i,1) =  %f +- %f \n', OBS_AJUST(i,1), OBS_AJUST(i,2));
end

fprintf('%f<=%f<=%f \n',qui1,qui_cal,qui2);
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
disp('Hipótese aceita no teste do Chi^2');
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')

%%  VERIFICAÇÃO DA MODELAGEM
%Aplicação das equações de coplinearidade inversa (R2 para R3)para a obtenção das
%coordenas geodésicas de pontos observados na imagem.

%Matriz de rotação  com omega, phi e kappa ajustados (obtidos na etapa anterior)
   r11 = cos(Xa(5))*cos(Xa(6));
   r12 = cos(Xa(4))*sin(Xa(6))+sin(Xa(4))*sin(Xa(5))*cos(Xa(6));
   r13 = sin(Xa(4))*sin(Xa(6))-cos(Xa(4))*sin(Xa(5))*cos(Xa(6));
   r21 = -cos(Xa(5))*sin(Xa(6));
   r22 = cos(Xa(4))*cos(Xa(6))-sin(Xa(4))*sin(Xa(5))*sin(Xa(6));
   r23 = sin(Xa(4))*cos(Xa(6))+cos(Xa(4))*sin(Xa(5))*sin(Xa(6));
   r31 = sin(Xa(5));
   r32 = -sin(Xa(4))*cos(Xa(5));
   r33 = cos(Xa(4))*cos(Xa(5));
 
%Equações de colinearidade Inversa (R2 para R3, com Z proveniente de fonte externa)

%Como não temos fonte externa para buscar o Z (MDT),os valores de Z
%buscados na componente Z dos pontos de apoi, já que estes não serão
%avaliados nesta transformação.

   for k =1:6
       
    Xv(k) = Xa(1)+(XYZV(3*k)-Xa(3))*((r11*xp(6+k)+r21*yp(6+k)+r31*(-f))/(r13*xp(6+k)+r23*yp(6+k)+r33*(-f)));
    Yv(k) = Xa(2)+(XYZV(3*k)-Xa(3))*((r12*xp(6+k)+r22*yp(6+k)+r32*(-f))/(r13*xp(6+k)+r23*yp(6+k)+r33*(-f)));
  
    %Discrepâncias entre os pontos calculados pela colinearidade e os
    %pontos de verificação
  
    dx(k)=Xv(k)-XYZV(3*k-2);
    
    dy(k)=Yv(k)-XYZV(3*k-1);
   end
     
fprintf('Discrepância entre as coordenadas X dos pontos calculados e dos pontos de verificação');
dx'
fprintf('Discrepância entre as coordenadas Y dos pontos calculados e dos pontos de verificação');
dy'
fprintf('Variância a posteriori');
varpost

%Média das discrepâncias na componente X
fprintf('Média dos resíduos na componente X');
mdx=mean(dx)
fprintf('Média dos resíduos na componente Y');
mdy=mean(dy)



%% Avaliação de tendências na modelagem pela colinearidade
%TESTE T Student
%Hipósete 0= modelagem livre de tendência (discrepâncias com media nula a um determinado nível de significância)
tx=(mdx/(std(dx)))*sqrt(n);
%tx (em módulo)deve ser menor que o t tabelado, dado em função do intervalo de confiança e do tamanho da amostra.
ty=(mdy/(std(dy)))*sqrt(n);






