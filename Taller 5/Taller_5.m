%Taller de series de Tiempo
%Manuel Andrés Niño Silva
%Métodos estocasticos en recursos hidraulicos
a =length(Caudales_diarios);
disp(a);
b =length(Caudales_diarios)/2;
disp(b);
%%  Estadisticos y momentos de Caudales Diarios

archivo_diario = 'Caudales Diarios.txt';
Caudales_diarios = load(archivo_diario);

%Estadisticos y momentos
Media_Diaria = mean(Caudales_diarios); %Media
S_Diaria = std(Caudales_diarios); %Desviación Estandar
S2_Diaria = S_Diaria^2; %Varianza
CV_diario = S_Diaria/Media_Diaria; %Coeficiente de variación
CS_diario = skewness(Caudales_diarios); %Coeficiente de asimetría
CC_diario = kurtosis(Caudales_diarios); %Coeficiente de curtosis
Cov_diaria = cov(Caudales_diarios(1:365),Caudales_diarios(366:730)); %Covarianza
Corr_diaria = corrcoef(Caudales_diarios(1:365),Caudales_diarios(366:730)); %Coeficiente de Correlación
[acf_diario,lag1]= autocorr(Caudales_diarios); %Autocorrelación ACF
[pacf_diario,lag2]= parcorr(Caudales_diarios); %Autocorrelación parcial PACF
[psd_diario,w_diario]= periodogram(Caudales_diarios); %Función de densidad espectral PSD

%Respuesta
fprintf('La media de los caudales diarios es: %.4f\n',Media_Diaria);disp(' ');
fprintf('La desviación estandar de los caudales diarios es: %.4f\n',S_Diaria);disp(' ');
fprintf('La variaanza de los caudales diarios es: %.4f\n',S2_Diaria);disp(' ');
fprintf('El coeficiente de variación de los caudales diarios es: %.4f\n',CV_diario);disp(' ');
fprintf('El coeficiente de asimetría de los caudales diarios es: %.4f\n lo que significa que está sesgada a la izquierda la distribución de datos' ,CS_diario);disp(' ');
fprintf('El coeficiente de Curtosis los caudales diarios es: %.4f\n',CC_diario);disp(' ');
fprintf('La covarianza de los caudales diarios del primer año y del segundo año es: %.4f\n',Cov_diaria(1,2));disp(' ');
fprintf('El coeficiente de correlación entre los caudales diarios del primer año y el segundo años es: %.4f\n',Corr_diaria(1,2));disp(' ');

%%  Estadisticos y momentos de Caudales mensuales
archivo_mensual = 'Caudales mensuales.txt';%Para crear los caudales mensuales se promediaron grupos de 30 datos diarios
Caudales_mensual = load(archivo_mensual);

%Estadisticos y momentos
Media_Mensual = mean(Caudales_mensual); %Media
S_Mensual = std(Caudales_mensual); %Desviación Estandar
S2_Mensual = S_Mensual^2; %Varianza
CV_Mensual= S_Mensual/Media_Mensual; %Coeficiente de variación
CS_Mensual = skewness(Caudales_mensual); %Coeficiente de asimetría
CC_Mensual = kurtosis(Caudales_mensual); %Coeficiente de curtosis
Cov_Mensual = cov(Caudales_mensual(1:12),Caudales_mensual(13:24)); %Covarianza
Corr_mensual = corrcoef(Caudales_mensual(1:12),Caudales_mensual(13:24)); %Coeficiente de Correlación
[acf_mensual,lag3]= autocorr(Caudales_mensual); %Autocorrelación ACF
[pacf_mensual,lag4]= parcorr(Caudales_mensual); %Autocorrelación parcial PACF
[psd_mensual,w_mensual]= periodogram(Caudales_mensual); %Función de densidad espectral PSD

%Respuesta
fprintf('La media de los caudales mensuales es: %.4f\n',Media_Mensual);disp(' ');
fprintf('La desviación estandar de los caudales mensuales es: %.4f\n',S_Mensual);disp(' ');
fprintf('La varianza de los caudales mensuales es: %.4f\n',S2_Mensual);disp(' ');
fprintf('El coeficiente de variación de los caudales mensuales es: %.4f\n',CV_Mensual);disp(' ');
fprintf('El coeficiente de asimetría de los caudales mensuales es: %.4f\n',CS_Mensual);disp(' ')
fprintf('El coeficiente de curtosis de los caudales mensuales es: %.4f\n',CC_Mensual);disp(' ');
fprintf('La covarianza de los caudales mensuales del primer y segundo año es: %.4f\n',Cov_Mensual(1,2));disp(' ');
fprintf('El coeficiente de correlación de los caudales mensuales del primer año y segundo año es: %.4f\n',Corr_mensual(1,2));disp(' ');

%%  Estadisticos y momentos de Cuadales Anuales

archivo_anual = 'Caudales Anuales.txt';%Los caudales anuales se obtuvieron promediando grupos de 12 datos de los caudales mensuales
Caudales_anuales = load(archivo_anual);

%Estadisticos y momentos
Media_Anual = mean(Caudales_anuales); %Media
S_Anual = std(Caudales_anuales); %Desviación Estandar
S2_Anual = S_Anual^2; %Varianza
CV_Anual = S_Anual/Media_Anual; %Coeficiente de variación
CS_Anual = skewness(Caudales_anuales); %Coeficiente de asimetría
CC_Anual = kurtosis(Caudales_anuales); %Coeficiente de curtosis
Cov_Anual = cov(Caudales_anuales(1:10),Caudales_anuales(11:20)); %Covarianza
Corr_anual = corrcoef(Caudales_anuales(1:10),Caudales_anuales(11:20)); %Coeficiente de Correlación
[acf_anual,lag5]= autocorr(Caudales_anuales); %Autocorrelación ACF
[pacf_anual,lag6]= parcorr(Caudales_anuales); %Autocorrelación parcial PACF
[psd_anual,w_anual]= periodogram(Caudales_anuales); %Función de densidad espectral PSD

%Respuesta
fprintf('La media de los caudales anuales es: %.4f\n',Media_Anual);disp(' ');
fprintf('La desviación estandar de los caudales anuales es: %.4f\n',S_Anual);disp(' ');
fprintf('La varianza de los caudales anuales es: %.4f\n',S2_Anual);disp(' ');
fprintf('El coeficiente de variación de los caudales anuales es: %.4f\n',CV_Anual);disp(' ');
fprintf('El coeficiente de asimetría de los caudales anuales es: %.4f\n',CS_Anual);disp(' ');
fprintf('El coeficiente de curtosis de los caudales anuales es: %.4f\n',CC_Anual);disp(' ');
fprintf('La covarianza de los caudales anuales entre los primeros 10 años y los siguientes 10 es: %.4f\n',Cov_Anual(1,2));disp(' ');
fprintf('El coeficiente de correlación de los caudales anuales entre los primeros 10 años y los siguientes 10 es: %.4f\n',Corr_anual(1,2));disp(' ');

%% Estadisticos y momentos de cuadales mensuales multianuales

archivo_mmulti = 'Caudal Mensual Multianual.txt';
Caudales_mmulti = load(archivo_mmulti);

%Estadisticos y momentos
Media_mmulti = mean(Caudales_mmulti); %Media
S_mmulti = std(Caudales_mmulti); %Desviación Estandar
S2_mmulti = S_mmulti^2; %Varianza
CV_mmulti = S_mmulti/Media_mmulti; %Coeficiente de variación
CS_mmulti = skewness(Caudales_mmulti); %Coeficiente de asimetría
CC_mmulti = kurtosis(Caudales_mmulti); %Coeficiente de curtosis
Cov_mmulti = cov(Caudales_mmulti); %Covarianza
Corr_mmulti = corrcoef(Caudales_mmulti); %Coeficiente de Correlación
[acf_mmulti,lag7]= autocorr(Caudales_mmulti); %Autocorrelación ACF
[pacf_mmulti,lag8]= parcorr(Caudales_mmulti); %Autocorrelación parcial PACF
[psd_mmulti,w_mmulti]= periodogram(Caudales_mmulti); %Función de densidad espectral PSD

%Respuesta
fprintf('La media de los caudales mensuales multianuales es: %.4f\n',Media_mmulti);disp(' ');
fprintf('La desviación estandar de los caudales mensuales multianuales es: %.4f\n',S_mmulti);disp(' ');
fprintf('La varianza de los caudales mensuales multianuales es: %.4f\n',S2_mmulti);disp(' ');
fprintf('El coeficiente de variación de los caudales mensuales multianuales es: %.4f\n',CV_mmulti);disp(' ');
fprintf('El coeficiebte de asimetría de los caudales mensuales multianuales es: %.4f\n',CS_mmulti);disp(' ');
fprintf('El coeficiente de curtosis de los caudales mensuales multianuales es: %.4f\n',CC_mmulti);disp(' ');
fprintf('La covarianza de los caudales mensuales multianuales es: %.4f\n Vemos que es igual que la varianza, ya que, esta serie se comparó consigo misma',Cov_mmulti);disp(' ');
fprintf('El coeficiente de correlación de los caudales mensuales multianuales es: %.0f\n Vemos que es igual a 1, ya que, esta serie se comparó consigo misma',Corr_mmulti);disp(' ');

%% Gráficas

%Serie de tiempo diaria
figure;
subplot(2,2,1);
plot(Caudales_diarios,'Color','red');
ylabel('Caudal (m3/s)');
xlabel('Día')
title('Serie de Caudales Diarios');
grid on; % Mostrar cuadrícula en el gráfico

%Serie de tiempo mensual
subplot(2,2,2);
plot(Caudales_mensual,'Color','blue');
ylabel('Caudal (m3/s)');
xlabel('Mes')
title('Serie de Caudales Mensuales');
grid on; % Mostrar cuadrícula en el gráfico

%Serie de tiempo anual
subplot(2,2,3);
plot(Caudales_anuales,'Color','Green');
ylabel('Caudal (m3/s)');
xlabel('Año')
title('Serie de Caudales Anuales');
grid on; % Mostrar cuadrícula en el gráfico

%Serie de tiempo mensual multianual
subplot(2,2,4);
meses = {'Enero', 'Febrero', 'Marzo', 'Abril', 'Mayo', 'Junio', 'Julio', 'Agosto', 'Septiembre', 'Octubre', 'Noviembre', 'Diciembre'};
plot(1:12,Caudales_mmulti,'Color','#6600a1');
xticks(1:12);
xticklabels(meses); 
xtickangle(90);
ylabel('Caudal (m3/s)');
xlabel('Mes')
title('Serie de Caudales Mensual Multianual');
grid on; % Mostrar cuadrícula en el gráfico

%% Correlograma

%Correlograma diario
figure;
subplot(2,2,1);
stem(lag1,acf_diario,"filled",'Color','red');
title('Correlograma Diario');
ylabel('Autocorrelación');
xlabel('Lag');

%Correlograma mensual
subplot(2,2,2);
stem(lag3,acf_mensual,"filled",'Color','blue');
title('Correlograma Mensual');
ylabel('Autocorrelación');
xlabel('Lag');

%Correlograma anual
subplot(2,2,3);
stem(lag5,acf_anual,"filled",'Color','green');
title('Correlograma Anual');
ylabel('Autocorrelación');
xlabel('Lag');

%Correlograma mensual multianual
subplot(2,2,4);
stem(lag7,acf_mmulti,"filled",'Color','#6600a1');
title('Correlograma Mensual Multianual');
ylabel('Autocorrelación');
xlabel('Lag');

%% Correlograma parcial

%Correlograma parcial diario
figure;
subplot(2,2,1);
stem(lag2,pacf_diario,"filled",'Color','red');
title('Correlograma Parcial Diario');
ylabel('Autocorrelación parcial');
xlabel('Lag');

%Correlograma parcial mensual
subplot(2,2,2);
stem(lag4,pacf_mensual,"filled",'Color','blue');
title('Correlograma Parcial Mensual');
ylabel('Autocorrelación parcial');
xlabel('Lag');

%Correlograma parcial anual
subplot(2,2,3);
stem(lag6,pacf_anual,"filled",'Color','green');
title('Correlograma Parcial Anual');
ylabel('Autocorrelación parcial');
xlabel('Lag');

%Correlograma parcial mensual multianual
subplot(2,2,4);
stem(lag8,pacf_mmulti,"filled",'Color','#6600a1');
title('Correlograma Parcial Mensual Multianual');
ylabel('Autocorrelación parcial');
xlabel('Lag');

%% Gráficas de los estadisticos

%Desviación estandar, Coeficiente de asimetria y Coeficiente de Curtosis
Y_S = [S_Diaria,S_Mensual,S_Anual,S_mmulti];
Y_CS =[CS_diario,CS_Mensual,CS_Anual,CS_mmulti];
Y_CC =[CC_diario,CC_Mensual,CC_Anual,CC_mmulti];
X_S = {'Diario','Mensual','Anual','M. Multianual'};

figure;
subplot(3,1,1);
plot(1:4,Y_S,'Color','red');
xticks(1:4);
xticklabels(X_S); 
xtickangle(90);
ylabel('Desviación Estandar');
title('Gráfica de la desviación estandar de las 4 Series');

subplot(3,1,2);
plot(1:4,Y_CS,'Color','blue');
xticks(1:4);
xticklabels(X_S); 
xtickangle(90);
ylabel('C. de Asimetría');
title('Gráfica del Coeficiente de Asimetría de las 4 Series');

subplot(3,1,3);
plot(1:4,Y_CC,'Color','green');
xticks(1:4);
xticklabels(X_S); 
xtickangle(90);
ylabel('C. de Curtosis');
title('Gráfica del Coeficiente de Curtosis de las 4 Series');

%% Peridogramas

%Periodogrma diario
figure;
subplot(2,2,1)
semilogy(w_diario,psd_diario,'Color','red');%Eje Y en escala Log
xlabel('Frecuencia')
ylabel('Densidad Espectral');
title('Periodograma Diario');
grid on;

%Periodogrma mensual
subplot(2,2,2)
semilogy(w_mensual,psd_mensual,'Color','blue');%Eje Y en escala Log
xlabel('Frecuencia')
ylabel('Densidad Espectral');
title('Periodograma Mensual');
grid on;

%Periodogrma anual
subplot(2,2,3)
semilogy(w_anual,psd_anual,'Color','green');%Eje Y en escala Log
xlabel('Frecuencia')
ylabel('Densidad Espectral');
title('Periodograma Anual');
grid on;

%Periodogrma mensual multianual
subplot(2,2,4)
semilogy(w_mmulti,psd_mmulti,'Color','#6600a1');%Eje Y en escala Log
xlabel('Frecuencia')
ylabel('Densidad Espectral');
title('Periodograma Mensual Multianual');
grid on;

%% Curvas IMD

%Curvas IMF diarias
[imf_diario,residual_diario,info_diaria] = emd(Caudales_diarios);
figure;
subplot(6,2,1);
plot(imf_diario(:, 1),'LineWidth',0.5);
title('Curva diaria imf 1')
subplot(6,2,2);
plot(imf_diario(:, 2),'LineWidth',0.5);
title('Curva diaria imf 2')
subplot(6,2,3);
plot(imf_diario(:, 3),'LineWidth',0.5);
title('Curva diaria imf 3')
subplot(6,2,4);
plot(imf_diario(:, 4),'LineWidth',1);
title('Curva diaria imf 4')
subplot(6,2,5);
plot(imf_diario(:, 5),'LineWidth',1);
title('Curva diaria imf 5')
subplot(6,2,6);
plot(imf_diario(:, 6),'LineWidth',1);
title('Curva diaria imf 6')
subplot(6,2,7);
plot(imf_diario(:, 7),'LineWidth',1);
title('Curva diaria imf 7')
subplot(6,2,8);
plot(imf_diario(:, 8),'LineWidth',1);
title('Curva diaria imf 8')
subplot(6,2,9);
plot(imf_diario(:, 9),'LineWidth',1);
title('Curva diaria imf 9')
subplot(6,2,10);
plot(imf_diario(:, 10),'LineWidth',1);
title('Curva diaria imf 10')
subplot(6,2,11:12);
plot(residual_diario,'LineWidth',1);
title('Residual Diario')

%Curvas IMF mensuales
[imf_mensual,residual_mensual,info_mensual] = emd(Caudales_mensual);
figure;
subplot(4,2,1);
plot(imf_mensual(:, 1),'LineWidth',0.5);
title('Curva mensual imf 1')
subplot(4,2,2);
plot(imf_mensual(:, 2),'LineWidth',0.5);
title('Curva mensual imf 2')
subplot(4,2,3);
plot(imf_mensual(:, 3),'LineWidth',0.5);
title('Curva mensual imf 3')
subplot(4,2,4);
plot(imf_mensual(:, 4),'LineWidth',1);
title('Curva mensual imf 4')
subplot(4,2,5);
plot(imf_mensual(:, 5),'LineWidth',1);
title('Curva mensual imf 5')
subplot(4,2,6);
plot(imf_mensual(:, 6),'LineWidth',1);
title('Curva mensual imf 6')
subplot(4,2,7);
plot(imf_mensual(:, 7),'LineWidth',1);
title('Curva mensual imf 7')
subplot(4,2,8);
plot(residual_mensual,'LineWidth',1);
title('Residual mensual')

%Curvas IMF anuales
[imf_anual,residual_anual,info_anual] = emd(Caudales_anuales);
figure;
subplot(2,2,1);
plot(imf_anual(:, 1),'LineWidth',0.5);
title('Curva anual imf 1')
subplot(2,2,2);
plot(imf_anual(:, 2),'LineWidth',0.5);
title('Curva anual imf 2')
subplot(2,2,3);
plot(imf_anual(:, 3),'LineWidth',0.5);
title('Curva anual imf 3')
subplot(2,2,4);
plot(residual_anual,'LineWidth',1);
title('Residual anual')

%Curvas IMF mensuales multianuales
[imf_mensual_multi,residual_mensual_multi,info_mensual_multi] = emd(Caudales_mmulti);
figure;
subplot(2,2,1);
plot(imf_mensual_multi(:, 1),'LineWidth',0.5);
title('Curva mensual multianual imf 1')
subplot(2,2,2);
plot(imf_mensual_multi(:, 2),'LineWidth',0.5);
title('Curva mensual multianual imf 2')
subplot(2,2,3:4);
plot(residual_mensual_multi,'LineWidth',1);
title('Residual mensual multianual')

%% Histogramas

%Diario
figure;
subplot(2,2,1)
histogram(Caudales_diarios,10,'FaceColor','red')
title('Caudales Diarios')
xlabel('Caudal (m3/s) ')
ylabel('Número de Datos')

%Mensual
subplot(2,2,2)
histogram(Caudales_mensual,10,'FaceColor','blue')
title('Caudales Mensuales')
xlabel('Caudal (m3/s) ')
ylabel('Número de Datos')

%Anual
subplot(2,2,3)
histogram(Caudales_anuales,5,'FaceColor','green')
title('Caudales Anuales')
xlabel('Caudal (m3/s)')
ylabel('Número de Datos')

%Mensual multianual
subplot(2,2,4)
histogram(Caudales_mmulti,3,'FaceColor','#6600a1')
title('Caudales Mensuales Multianuales')
xlabel('Caudal (m3/s)')
ylabel('Número de Datos')

%% Diagramas de cajas y bigotes

%Diario
figure;
subplot(4,1,1)
boxplot(Caudales_diarios,'Orientation','horizontal','OutlierSize',2);
title('Caudales Diarios')

%Mensual
subplot(4,1,2)
boxplot(Caudales_mensual,'Orientation','horizontal');
title('Caudales Mensuales')

%Anual
subplot(4,1,3)
boxplot(Caudales_anuales,'Orientation','horizontal');
title('Caudales Anuales')

%Mensual Multianual
subplot(4,1,4)
boxplot(Caudales_mmulti,'Orientation','horizontal');
title('Caudales Mensuales Multianuales')

%% Pruebas de hipotesis

% Pruebas de analisis de Saltos

% Prueba de pettitt
pettitt_diario = pettitt(Caudales_diarios(1:1000));
disp(['Ubicación del cambio: ', num2str(pettitt_diario(1))]);
disp(['Estadístico de prueba K: ', num2str(pettitt_diario(2))]);
disp(['Valor p: ', num2str(pettitt_diario(3))]);disp(' ');

pettitt_mensual = pettitt(Caudales_mensual);
disp(['Ubicación del cambio: ', num2str(pettitt_mensual(1))]);
disp(['Estadístico de prueba K: ', num2str(pettitt_mensual(2))]);
disp(['Valor p: ', num2str(pettitt_mensual(3))]);disp(' ');

pettitt_anual = pettitt(Caudales_anuales(1:20));
disp(['Ubicación del cambio: ', num2str(pettitt_anual(1))]);
disp(['Estadístico de prueba K: ', num2str(pettitt_anual(2))]);
disp(['Valor p: ', num2str(pettitt_anual(3))]);disp(' ');

pettitt_mmulti = pettitt(Caudales_mmulti);
disp(['Ubicación del cambio: ', num2str(pettitt_mmulti(1))]);
disp(['Estadístico de prueba K: ', num2str(pettitt_mmulti(2))]);
disp(['Valor p: ', num2str(pettitt_mmulti(3))]);disp(' ');

%%
%Prueba de suma de rangos
[p_dia_mes,h_dia_mes] = ranksum(Caudales_diarios,Caudales_mensual);
[p_dia_year,h_dia_year] = ranksum(Caudales_diarios,Caudales_anuales);
[p_dia_mmulti,h_dia_mmulti] = ranksum(Caudales_diarios,Caudales_mmulti);
[p_mes_year,h_mes_year] = ranksum(Caudales_mensual,Caudales_anuales);
[p_mes_mmulti,h_mes_mmulti] = ranksum(Caudales_mensual,Caudales_mmulti);
[p_year_mmulti,h_year_mmulti] = ranksum(Caudales_anuales,Caudales_mmulti);

fprintf('Loos valores de P y H entre las series diaria y mensual son: %.6f y %d.\n',p_dia_mes, h_dia_mes);
fprintf('Loos valores de P y H entre las series diaria y anual son: %.6f y %d.\n',p_dia_year, h_dia_year);
fprintf('Loos valores de P y H entre las series diaria y mensual multianual son: %.6f y %d.\n',p_dia_mmulti, h_dia_mmulti);
fprintf('Loos valores de P y H entre las series mensual y anual son: %.6f y %d.\n',p_mes_year, h_mes_year);
fprintf('Loos valores de P y H entre las series mensual y mensual multianual son: %.6f y %d.\n',p_mes_mmulti, h_mes_mmulti);
fprintf('Loos valores de P y H entre las series anual y mensual multianual son: %.6f y %d.\n',p_year_mmulti, h_year_mmulti);


%%
%Prueba de desviaciones acumuladas

%%
%Prueba de Kruskal-Wallis
p1_diario = kruskalwallis(Caudales_diarios);
p1_mensual = kruskalwallis(Caudales_mensual);
p1_anual = kruskalwallis(Caudales_anuales);
p1_mmulti = kruskalwallis(Caudales_mmulti);

%%
%Prueba CUSUM
h2_diario = cusumtest(Caudales_diarios(1:3716),Caudales_diarios(3717:7432));
h2_mensual = cusumtest(Caudales_mensual(1:124),Caudales_mensual(125:248));
h2_anual = cusumtest(Caudales_anuales(1:10),Caudales_anuales(11:20));
h2_mmulti = cusumtest(Caudales_mmulti(1:6),Caudales_mmulti(7:12));


fprintf('El valor de H para la CUSUM de los caudales diarios es: %.d',h2_diario );disp(' ');
fprintf('El valor de H para la CUSUM de los mensuales diarios es: %.d',h2_mensual);disp(' ');
fprintf('El valor de H para la CUSUM de los caudales anuales es: %.f',h2_anual);disp(' ');
fprintf('El valor de H para la CUSUM de los caudales mensuales multianuales es: %.f',h2_mmulti );disp(' ');
%%
%Prueba t de Student
h3_diario = ttest(Caudales_diarios,Media_Diaria);
h3_mensual = ttest(Caudales_mensual,Media_Mensual);
h3_anual = ttest(Caudales_anuales,Media_Anual);
h3_mmulti = ttest(Caudales_mmulti,Media_mmulti);

fprintf('El valor de  H para la serie diaria es: %d.\n', h3_diario);disp(' ');
fprintf('El valor de  H para la serie mensual es: %d.\n',h3_mensual);disp(' ');
fprintf('El valor de  H para la serie anual es: %d.\n',h3_anual);disp(' ');
fprintf('El valor de  H para la serie mensual multianual es: %d.\n',h3_mmulti);disp(' ');

%Prueba de la Relación de Verosimilitud de Worsley
%Prueba de Siegel-Tukey

%% Pruebas de hipótesis de tendencias

%Prueba de coeficiente de correlación ρ de Spearman
rho_diario = corr(Caudales_diarios(1:365), Caudales_diarios(366:730), 'type', 'Spearman');
rho_mensual = corr(Caudales_mensual(1:12),Caudales_mensual(13:24), 'type', 'Spearman');
rho_anual = corr(Caudales_anuales(1:10),Caudales_anuales(11:20), 'type', 'Spearman');
rho_mmulti = corr(Caudales_mmulti(1:6),Caudales_mmulti(7:12), 'type', 'Spearman');

fprintf('El valor del coeficiente de Spearman para la serie diaria es: %4f.\n', rho_diario);disp(' ');
fprintf('El valor del coeficiente de Spearman para la serie mensual es: %4f.\n', rho_mensual);disp(' ');
fprintf('El valor del coeficiente de Spearman para la serie anual es: %4f.\n', rho_anual);disp(' ');
fprintf('El valor del coeficiente de Spearman para la serie mensual multianual es: %4f.\n', rho_mmulti);disp(' ');

%%
%Prueba de Mann-Kendall
[h4_diario,p4_diario]=Mann_Kendall(Caudales_diarios,0.05);
[h4_mensual,p4_mensual]=Mann_Kendall(Caudales_mensual,0.05);
[h4_anual,p4_anual]=Mann_Kendall(Caudales_anuales,0.05);
[h4_mmulti,p4_mmulti]=Mann_Kendall(Caudales_mmulti,0.05);

fprintf('Los valores de P y H para la serie diaria son: %.6f y %d.\n',p4_diario, h4_diario);
fprintf('Los valores de P y H para la serie mensual son: %.6f y %d.\n',p4_mensual, h4_mensual);
fprintf('Los valores de P y H para la serie anual son: %.6f y %d.\n',p4_anual, h4_anual);
fprintf('Los valores de P y H para la serie mensual multianual son: %.6f y %d.\n',p4_mmulti, h4_mmulti);
%%
%Prueba de Mann-Kendall modificada
[h5_diario,p5_diario]=Mann_Kendall_Modified(Caudales_diarios,0.05);
[h5_mensual,p5_mensual]=Mann_Kendall_Modified(Caudales_mensual,0.05);
[h5_anual,p5_anual]=Mann_Kendall_Modified(Caudales_anuales,0.05);
[h5_mmulti,p5_mmulti]=Mann_Kendall_Modified(Caudales_mmulti,0.05);

fprintf('Los valores de P y H para la serie diaria son: %.6f y %d.\n',p5_diario, h5_diario);
fprintf('Los valores de P y H para la serie mensual son: %.6f y %d.\n',p5_mensual, h5_mensual);
fprintf('Los valores de P y H para la serie anual son: %.6f y %d.\n',p5_anual, h5_anual);
fprintf('Los valores de P y H para la serie mensual multianual son: %.6f y %d.\n',p5_mmulti, h5_mmulti);

%%
%Prueba de Regresión Lineal
Reg_diaria = fitlm(Caudales_diarios(731:1095),Caudales_diarios(1096:1460));
Reg_mensual = fitlm(Caudales_mensual(25:36),Caudales_mensual(37:48));
Reg_anual = fitlm(Caudales_anuales(1:10),Caudales_mensual(11:20));
Reg_mmulti = fitlm(Caudales_mmulti(1:6),Caudales_mmulti(7:12));

R2_diario = Reg_diaria.Rsquared.Ordinary;
R2_mensual = Reg_mensual.Rsquared.Ordinary;
R2_anual = Reg_anual.Rsquared.Ordinary;
R2_mmulti = Reg_mmulti.Rsquared.Ordinary;

fprintf('El valor R cuadrado para la serie diaria es: %.6f.\n',R2_diario);
fprintf('El valor R cuadrado para la serie mensual es: %.6f.\n',R2_mensual);
fprintf('El valor R cuadrado para la serie anual es: %.6f.\n',R2_anual);
fprintf('El valor R cuadrado para la serie mensual multianual es: %.6f.\n',R2_mmulti);

%% Pruebas de hipotesis de cambio de distribución

%Prueba de Kologorov-Smirnov
[h6_diario,p6_diario] = kstest((Caudales_diarios-Media_Diaria)/S_Diaria); %Normalizamos, ya que, la prueba se aplica a una distribución normal estandar 
[h6_mensual,p6_mensual] = kstest((Caudales_mensual-Media_Mensual)/S_Mensual);
[h6_anual,p6_anual] = kstest((Caudales_anuales-Media_Anual)/S_Anual);
[h6_mmulti,p6_mmulti] = kstest((Caudales_mmulti-Media_mmulti)/S_mmulti);

fprintf('Los valores de P y H para la serie diaria son: %.6f y %d.\n',p6_diario, h6_diario);
fprintf('Los valores de P y H para la serie mensual son: %.6f y %d.\n',p6_mensual, h6_mensual);
fprintf('Los valores de P y H para la serie anual son: %.6f y %d.\n',p6_anual, h6_anual);
fprintf('Los valores de P y H para la serie mensual multianual son: %.6f y %d.\n',p6_mmulti, h6_mmulti);

%%
%Prueba de Wald-Wolfowitz
[h7_diario,p7_diario,~] = runstest(Caudales_diarios);
[h7_mensual,p7_mensual,~] = runstest(Caudales_mensual);
[h7_anual,p7_anual,~] = runstest(Caudales_anuales);
[h7_mmulti,p7_mmulti,~] = runstest(Caudales_mmulti);

fprintf('Los valores de P y H para la serie diaria son: %.6f y %d.\n',p7_diario, h7_diario);
fprintf('Los valores de P y H para la serie mensual son: %.6f y %d.\n',p7_mensual, h7_mensual);
fprintf('Los valores de P y H para la serie anual son: %.6f y %d.\n',p7_anual, h7_anual);
fprintf('Los valores de P y H para la serie mensual multianual son: %.6f y %d.\n',p7_mmulti, h7_mmulti);

%Prueba de Terry

%% Pruebas de hipótesis de normalidad

%Prueba Chi Cuadrado
h9_diario = chi2gof(Caudales_diarios);
h9_mensual = chi2gof(Caudales_mensual);
h9_anual = chi2gof(Caudales_anuales);
h9_mmulti = chi2gof(Caudales_mmulti);

fprintf('El valor de H para la serie diaria es: %d.\n', h9_diario);
fprintf('El valor de H para la serie mensual es: %d.\n', h9_mensual);
fprintf('El valor de H para la serie anual es: %d.\n', h9_anual);
fprintf('El valor de H para la serie mensual multianual es: %d.\n',h9_mmulti);

%%
%Prueba de Anderson-Darling
[h10_diario,p10_diario] = adtest(Caudales_diarios);
[h10_mensual,p10_mensual] = adtest(Caudales_mensual);
[h10_anual,p10_anual] = adtest(Caudales_anuales);
[h10_mmulti,p10_mmulti] = adtest(Caudales_mmulti);

fprintf('Los valores de P y H para la serie diaria son: %.6f y %d.\n',p10_diario, h10_diario);
fprintf('Los valores de P y H para la serie mensual son: %.6f y %d.\n',p10_mensual, h10_mensual);
fprintf('Los valores de P y H para la serie anual son: %.6f y %d.\n',p10_anual, h10_anual);
fprintf('Los valores de P y H para la serie mensual multianual son: %.6f y %d.\n',p10_mmulti, h10_mmulti);

%% Pruebas de hipótesis de aleatoriedad e independencia

%Prueba de cruces de Mediana
%%
%Prueba de puntos de cambio

media_movil_diaria = movmean(Caudales_diarios, 10); % Utiliza una ventana de media móvil de tamaño 10. Calcular la media móvil de los datos
diferencia_diaria = abs(Caudales_diarios - media_movil_diaria);
umbral_diario = 2*S_Diaria; % Umbral para la detección de cambios
puntos_cambio = find(diferencia_diaria > umbral_diario);
disp('Puntos de cambio detectados en la serie diaria:');
disp(puntos_cambio);% Mostrar los puntos de cambio

%%
%Prueba de diferencia de rangos
[p13_diario,h13_diario] = signrank(Caudales_diarios);
[p13_mensual,h13_mensual] = signrank(Caudales_mensual);
[p13_anual,h13_anual] = signrank(Caudales_anuales);
[p13_mmulti,h13_mmulti] = signrank(Caudales_mmulti);

fprintf('Los valores de P y H para la serie diaria son: %.6f y %d.\n',p13_diario, h13_diario);
fprintf('Los valores de P y H para la serie mensual son: %.6f y %d.\n',p13_mensual, h13_mensual);
fprintf('Los valores de P y H para la serie anual son: %.6f y %d.\n',p13_anual, h13_anual);
fprintf('Los valores de P y H para la serie mensual multianual son: %.6f y %d.\n',p13_mmulti, h13_mmulti);

%%
%Prueba de autocorrelación

%%
%Prueba de autocorrelación de Barlett

%%
%Prueba de von Neumann
[q_prima_diario,grafica_diario,Yi,Xi,a] = VonNewman(Caudales_diarios);
[q_prima_mensual,grafica_mensual] = VonNewman(Caudales_mensual);
[q_prima_anual,grafica_anual] = VonNewman(Caudales_anuales);
[q_prima_mmulti,grafica_mmulti] = VonNewman(Caudales_mmulti);

fprintf('El valor de Q_prima diario es: %.8f ',q_prima_diario);disp(' ');
fprintf('El valor de Q_prima mensual es: %.8f ',q_prima_mensual);disp(' ');
fprintf('El valor de Q_prima anual es: %.8f ',q_prima_anual);disp(' ');
fprintf('El valor de Q_prima mensual multianual es: %.8f ',q_prima_mmulti);disp(' ');

%% Funciones
%Calculo de las sumatorias a partir de las que se obtiene Q
function [q_prima,grafica,Promedio,Suma_1,Suma_2] = VonNewman(Serie)

Promedio = mean(Serie);
Yi_menos_Yi_1 = zeros(length(Serie)-1);
Yi_menos_Y_promedio = zeros(length(Serie)-1);

for i = 2:length(Serie)
    X = Serie(i)-Serie(i-1);
    Yi_menos_Yi_1(i-1) = X^2;
end

for i = 1:length(Serie)
    Z = Serie(i)-Promedio;
    Yi_menos_Y_promedio(i)=Z^2;
end

Suma_1 = sum(Yi_menos_Yi_1);
Suma_2 = sum(Yi_menos_Y_promedio);

Q = Suma_1/Suma_2; %Calculo del valor de Q

n = length(Serie); %Numero de datos
q_prima = (Q-2)/(2*sqrt((n-2)*(n^2-1))); %Valor de Q'

eje_x = linspace(-3,3,200); %Valores del eje X
I_confianza = 1.96; %Valor de Z para un intervalo de confianza del 95%
dist_normal = (1/sqrt(2*pi)) * exp(-0.5*eje_x.^2); %Ecuación de la distribución normal estandar, como se distribuye Q prima

%Crear la gráfica de la distribución normal
grafica = figure;
title('Ubicación de Q prima')
plot(eje_x, dist_normal, 'b', 'LineWidth', 2);
hold on;

%Sombrear la zona en donde la hipotesis nula se cumple
x_sombra = linspace(-I_confianza,I_confianza,100);
y_sombra = (1/sqrt(2*pi)) * exp(-0.5*x_sombra.^2);
fill([x_sombra, fliplr(x_sombra)],[y_sombra, zeros(size(y_sombra))], 'r', 'FaceAlpha', 0.3)

%Ubicación del valor de Q prima en el gráfico de la distribución normal
Ubicacion_Q_prima_Y = (1/sqrt(2*pi)) * exp(-0.5*q_prima.^2);
plot([q_prima, q_prima], [0, Ubicacion_Q_prima_Y], 'g', 'LineWidth', 2);
hold off;
end