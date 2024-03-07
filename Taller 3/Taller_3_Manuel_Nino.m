%%
%Taller valores Extremos
%Manuel Andrés Niño Silva
%Codigo de Matlab

Caudales = [1266,1492,1862,861,715,1367,1837,1429,1429,1261,1607,2132,1652,1537,1155,1899,1956,1596,1380,745,2181,955,1007,824,1271,1044,597,1216,1061,1450,2016,1270,1341,1075,2482,874,1689,1554,1831,1149,1444,1791,1207,3050,1578,2817,792,1143,1698,2076,1204,1835];

Promedio = mean(Caudales);
Yi_menos_Yi_1 = zeros(size(Caudales)-1);
Yi_menos_Y_promedio = zeros(size(Caudales)-1);

for i = 2:size(Caudales,2)
    X = Caudales(i)-Caudales(i-1);
    Yi_menos_Yi_1(i-1) = X^2;
end

for i = 1:size(Caudales,2)
    Z = Caudales(i)-Promedio;
    Yi_menos_Y_promedio(i)=Z^2;
end

%Calculo de las sumatorias a partir de las que se obtiene Q
Suma_1 = sum(Yi_menos_Yi_1);
Suma_2 = sum(Yi_menos_Y_promedio);

Q = Suma_1/Suma_2; %Calculo del valor de Q

n = size(Caudales,2); %Numero de datos
Q_prima = (Q-2)/(2*sqrt((n-2)*(n^2-1))); %Valor de Q'
fprintf('El valor de Qprima es: %.6f\n',Q_prima)
fprintf('Y se ubica tal y como muestra el gráfico de la distribución normal estandar')

x = linspace(-3,3,200); %Valores del eje X
I_confianza = 1.96; %Valor de Z para un intervalo de confianza del 95%
dist_normal = (1/sqrt(2*pi)) * exp(-0.5*x.^2); %Ecuación de la distribución normal estandar, como se distribuye Q prima

%Crear la gráfica de la distribución normal
figure;
plot(x, dist_normal, 'b', 'LineWidth', 2);
hold on;

%Sombrear la zona en donde la hipotesis nula se cumple
x_sombra = linspace(-1.96,1.96,100);
y_sombra = (1/sqrt(2*pi)) * exp(-0.5*x_sombra.^2);
fill([x_sombra, fliplr(x_sombra)],[y_sombra, zeros(size(y_sombra))], 'r', 'FaceAlpha', 0.3)

%Ubicación del valor de Q prima en el gráfico de la distribución normal
Ubicacion_Q_prima_Y = (1/sqrt(2*pi)) * exp(-0.5*Q_prima.^2);
plot([Q_prima, Q_prima], [0, Ubicacion_Q_prima_Y], 'g', 'LineWidth', 2);
hold off;

%Etiquetas y titulo
xlabel('Q_p_r_i_m_a');
ylabel('P');
title('Comprobación de la hipotesis nula');

% Mostrar leyenda
legend('Distribución de Q_p_r_i_m_a', 'Área Hipotesis nula', 'Ubicación Q_p_r_i_m_a', 'Location', 'northwest');

%Respuesta
disp(' ')
fprintf('Ya que Qprima se ubica dentro del área sombreada, concluimos que los datos son independientes')
disp(' ')
fprintf('Esto se comprueba en la gráfica de Qi+1 vs Qi')

%Grafica de Q_i+1 vs Q_i
Caudales_i1 = Caudales(2:end);
Caudales_i = Caudales(1:end-1);
figure;
scatter(Caudales_i,Caudales_i1)
xlabel('Qi')
ylabel('Qi+1')
title('Qi vs Qi+1')

%%
%Ejercicio 2
Yi_menos_Yj = []; %Crear una lista para introducir el valor de la diferencia entre Y_i-Y_j
for i = 2:size(Caudales,2)
    Resta = Caudales(i)-Caudales(i-1);
    Yi_menos_Yj (i-1) = Resta;
end

%Crear un vector en el que se guarde el valor de funcion signo para cada
%valor de la resta de Y_i-Y_j
Vector_signo = [];

%Aplicar la función signo a cada valor en la lista Y_i_menos_Y_j y se
%guarda en la lista Vector_signo
for i = 1:size(Yi_menos_Yj,2)
    Signo = sign(Yi_menos_Yj(i));
    Vector_signo(i)=Signo;
end

T = sum(Vector_signo); %Sumar todos los valores obtenidos al aplicar la función signo, valor T
n = size(Caudales,2); %Número de datos
Var_T_prima = 1/18*(n*(n-1)*(2*n+5)); %Calculo de la varianza de T prima
media_T_prima = 0; %Media de te prima
T_prima = 18*T/(n*(n-1)*(2*n+5)); %Calculo de T_prima

%Funcion de distribución de T prima con la varianza calculada y media igual a 0
x = linspace(-600,600,10000);
Funt_prob_T_prima = (1/(sqrt(Var_T_prima)*sqrt(2*pi))) * exp(-0.5*(x.^2/Var_T_prima));
figure;
plot(x,Funt_prob_T_prima,'b','LineWidth',2);
hold on;

z = 1.96; %Valor Z para intervalo de confianza del 95% en distribución normal estandar
Z = media_T_prima+z*sqrt(Var_T_prima)/sqrt(n); %Calculo del valor límite en el eje X para que se cumpla la hipotesis nula

%Sombrear la zona en la que se cumple la hipotesis nula
x_sombra = linspace(-Z,Z,1000);
y_sombra = (1/(sqrt(Var_T_prima)*sqrt(2*pi))) * exp(-0.5*(x_sombra.^2/Var_T_prima));
fill([x_sombra, fliplr(x_sombra)],[y_sombra, zeros(size(y_sombra))], 'r', 'FaceAlpha', 0.3);

%Ubicar el valor de T_prima en la gráfica de la función de distribución
Ubicacion_T_prima_Y =(1/(sqrt(Var_T_prima)*sqrt(2*pi))) * exp(-0.5*(T_prima.^2/Var_T_prima)) ;
plot([T_prima, T_prima], [0, Ubicacion_T_prima_Y], 'g', 'LineWidth', 2);
hold off

fprintf('El valor de Tprima es: %.6f\n',T_prima)

%Etiquetas y titulo
xlabel('T_p_r_i_m_a');
ylabel('P');
title('Comprobación de la hipotesis nula');

% Mostrar leyenda
legend('Distribución de T_p_r_i_m_a', 'Área Hipotesis nula', 'Ubicación T_p_r_i_m_a', 'Location', 'northwest');
disp(' ')
fprintf('Ya que Tprima se ubica dentro del área sombreada, concluimos que los datos no tienen tendencia')

%%
%Ejercicio 3
%A partir de una función creada por matlab, se hace la gráfica del ajuste
%de probabilidad acumulada y se hace la curva en papel gumbel aunque
%muestra la probabilidad de que NO se supere un determinado valor de caudal
figure;
createFitGumbel(Caudales);

%Probabilidad de la inundación de 1000 años ya que P = 1/T, en donde T es el periodo de retorno
P_inundacion_1000_years = 1/1000; 

%Parametros de la distribucíon
k = -0.05138129257569;
sigma = 420.924472877;
mu = 1252.04669596;

%Variables para la constucción de la g´rafica del ajuste de distribución
%gumbel
Rango_Caudal = linspace(0,4000,40000);
P_gumbel=1-gevcdf(Rango_Caudal,k,sigma,mu); %Se corrige la distribución para que muestre la probabilidad de que sea excedido un valor de caudal

%Crea la gráfica
figure;
plot(Rango_Caudal,P_gumbel)
hold on

%Estimación del caudal de mil años hasta lograr encontrar el valor para el
%que la probabilidad de excedencia sea de 0.001
Caudal_estimado_1000_years = 3699;
P_gumbel_1000_years=1-gevcdf(Caudal_estimado_1000_years,k,sigma,mu);%Se halló el caudal para que esto de 0.001

%Gráfica el valor de caudal de la inundación de mil años
plot([Caudal_estimado_1000_years,Caudal_estimado_1000_years],[0,1])

%Etiquetas
xlabel('Caudal (m3/s)')
ylabel('Probabilidad de ocurrencia')
title('Probabilidad de inundación')

%Leyenda
legend('Distribución gumbel','Ubicación del evento con T = 1000 años','Location','northwest')

fprintf('La probabilidad de que ocurra el cuadal estimado para la inundación de 1000 años es: %.3f\n',P_gumbel_1000_years)
disp(['El valor de caudal para la inundación con perido de retorno igual a 1000 años es:',num2str(Caudal_estimado_1000_years),' m3/s'])
disp('Ya que su probabilidad de ocurrencia es de 0.001 o 0.1%')

%%
%Ejercicio 4
Caudales_ordenados = sort(Caudales,2,"ascend"); %Ordenar caudales de mayor a menor
Eje_Y = []; %Lista para guardar los valores de -ln(-ln(i/(n+1)))
n = size(Caudales_ordenados,2);

%Calculo de -ln(-ln(i/(n+1))) para cada caudal
for i=1:size(Caudales_ordenados,2)
    V_Gumbel = -log(-log(i/(n+1)));
    Eje_Y(i) = V_Gumbel;
end

%Regresión lineal a la gráfica de Caudales vs -ln(-ln(i/(n+1)))
regresion = polyfit(Caudales_ordenados,Eje_Y, 1);
x_regresion = min(Caudales_ordenados):0.1:max(Caudales_ordenados);
y_regresion = polyval(regresion, x_regresion);

%Crea la figura de la regresión lineal y datos originales
figure;
plot(Caudales_ordenados, Eje_Y, 'o', 'MarkerSize', 5);  % Grafica los datos originales
hold on;
plot(x_regresion, y_regresion, 'r', 'LineWidth', 2);  % Grafica la línea de regresión lineal
hold off;

xlabel('Caudal (m3/s)')
ylabel('Regresión')
title('Regresión lineal de Q vs -ln(-ln(i/(n+1)))')
legend('Datos originales', 'Regresión lineal','Location', 'northwest');

%Oncluir la ecuación de regresión
text(2000, -1, sprintf('Y = %.4fX + %.4f', regresion(1), regresion(2)), 'FontSize', 8);

%La función de distribución acumulada de gumbel es F(X) = exp(-exp(-(x-μ)/β)
%μ = mu en matlab y β = sigma en matlab
%Al hacer aplicar el doble logaritmo a la ecuación de gumbel nos queda (x-μ)/β
%y al hacer la regresión, obtenemos una ecuación de la forma ax+b = x/β-μ/β por lo tanto:
%ax = x/β => β = 1/a 
%b = -μ/β => μ =-bβ = -b/a
a = regresion(1);
b = regresion(2);
mu_regresion = -b/a;
Sigma_regresion = 1/a;

%Gráfica de la función de distribución con los parametros obtenidos
Rango_Caudal = linspace(0,5000,40000);
figure;
P_gumbel=1-gevcdf(Rango_Caudal,0,Sigma_regresion,mu_regresion);
plot(Rango_Caudal,P_gumbel)
hold on

%Estimación del caudal para el cual la probabilidad de superarlo sea 0.001
%teniendo en cuenta los parametros establecidos en la regresión
Caudal_estimado_1000_years_regresion = 4250;
Probabilidad_caudal_1000_years = 1-gevcdf(Caudal_estimado_1000_years_regresion,0,Sigma_regresion,mu_regresion);
disp(['El valor de la probablidad de ocurrencia es ',num2str(Probabilidad_caudal_1000_years,'%.4f\n'),' para una caudal de ',num2str(Caudal_estimado_1000_years_regresion),' m3/s']);

%Intervalos de confianza
%Valor de K para un T = 1000 años
K_T_1000 = -sqrt(6)/pi*(0.577+log(log(1/gevcdf(Caudal_estimado_1000_years_regresion,0,Sigma_regresion,mu_regresion))));
delta = sqrt(1+1.1396*K_T_1000+1.1*K_T_1000^2); 
d_estandar_gumbel = Sigma_regresion*pi/sqrt(6); %Se calcula la desviación estandar de la función de distribución usando los parametros establecidos
n = 52;
Error_estandar = d_estandar_gumbel*delta/sqrt(n); %E_estandar = d_estandar*δ/sqrt(n)
Z_alpha_sobre_2 = 1.960; %Valor de Z para nivel de significancia del 5%

Limite_confianza_superior = Caudal_estimado_1000_years_regresion+Z_alpha_sobre_2*Error_estandar; %Limite superior = Q_T + Z_a/2*E_estandar
Limite_confianza_inferior = Caudal_estimado_1000_years_regresion-Z_alpha_sobre_2*Error_estandar; %Limite superior = Q_T - Z_a/2*E_estandar

disp(['El limite de confianza inferior para la inundación de 1000 años es: ',num2str(Limite_confianza_inferior,'%.4f\n'),' m3/s'])
disp(['El limite de confianza superior para la inundación de 1000 años es: ',num2str(Limite_confianza_superior,'%.4f\n'),' m3/s'])

%Gráficar el caudal estimado de 1000 años y sus límites de confianza
plot([Caudal_estimado_1000_years_regresion,Caudal_estimado_1000_years_regresion],[0,1])
plot([Limite_confianza_inferior,Limite_confianza_inferior],[0,1])
plot([Limite_confianza_superior,Limite_confianza_superior],[0,1])
hold off

xlabel('Caudal (m3/s)')
ylabel('Probabilidad')
title('Ajuste de probabilidad gumbel usando regresión')
legend('Ajuste de distribución', 'Caudal 1000 años', 'Limite inferior','Limite Superior', 'Location', 'northwest');

%%
%Ejercicio 5

% Calcular media y varianza muestrales (Momentos)
media_caudales = mean(Caudales);
varianza_caudales = var(Caudales);

% Constante de Euler-Mascheroni
gamma = 0.5772156649;

% Resolver ecuaciones para los parámetros de la distribución Gumbel
Sigma_momentos = sqrt(6 * varianza_caudales / pi^2);
mu_momentos = media_caudales - gamma * Sigma_momentos;

% Mostrar los parámetros
fprintf('Parámetros de la distribución Gumbel:\n');
fprintf('Ubicación (mu): %.4f\n', mu_momentos);
fprintf('Escala (Sigma): %.4f\n', Sigma_momentos);

Caudal_estimado_1000_years_momentos = 4000; %Estimación de caudal con probabilidad de excedencia de 0.001 teniendo en cuenta los parametros calculados

%Gráfica del ajuste de probabilidad por el metodo de momentos
Rango_Caudal_momentos = linspace(0,4000,40000);
figure;
P_gumbel_momentos=1-gevcdf(Rango_Caudal_momentos,0,Sigma_momentos,mu_momentos);
plot(Rango_Caudal_momentos,P_gumbel_momentos)
hold on

Probabilidad_caudal_1000_years_momentos = 1-gevcdf(Caudal_estimado_1000_years_momentos,0,Sigma_momentos,mu_momentos);
disp(['El valor de caudal con probabilidad de excedencia de 0.001 (Tr = 1000 años) es: ', num2str(Caudal_estimado_1000_years_momentos),' m3/s']);

%Intervalos de confianza
K_T_1000 = -sqrt(6)/pi*(0.577+log(log(1/gevcdf(Caudal_estimado_1000_years_momentos,0,Sigma_momentos,mu_momentos))));
delta_momentos = sqrt(1+1.1396*K_T_1000+1.1*K_T_1000^2);
d_estandar_gumbel_momentos = Sigma_momentos*pi/sqrt(6);
n = 52;
Error_estandar_momentos = d_estandar_gumbel_momentos*delta_momentos/sqrt(n);
Z_alpha_sobre_2 = 1.960; %Nivel de significancia del 5%

Limite_confianza_superior_momentos = Caudal_estimado_1000_years_momentos+Z_alpha_sobre_2*Error_estandar_momentos;
Limite_confianza_inferior_momentos = Caudal_estimado_1000_years_momentos-Z_alpha_sobre_2*Error_estandar_momentos;

disp(['El limite de confianza inferior para la inundación de 1000 años es: ',num2str(Limite_confianza_inferior_momentos,'%.4f\n'),' m3/s'])
disp(['El limite de confianza superior para la inundación de 1000 años es: ',num2str(Limite_confianza_superior_momentos,'%.4f\n'),' m3/s'])

%%Gráficar el caudal estimado de 1000 años y sus límites de confianza
plot([Caudal_estimado_1000_years_momentos,Caudal_estimado_1000_years_momentos],[0,1])
plot([Limite_confianza_inferior_momentos,Limite_confianza_inferior_momentos],[0,1])
plot([Limite_confianza_superior_momentos,Limite_confianza_superior_momentos],[0,1])
hold off

xlabel('Caudal (m3/s)')
ylabel('Probabilidad')
title('Ajuste de probabilidad gumbel usando regresión')
legend('Ajuste de distribución', 'Caudal 1000 años', 'Limite inferior','Limite Superior', 'Location', 'northwest');

%%
%Ejercicio 6
%Función que ajusta los datos a una distribución LogNormal y lo gráfica 
createFitLogNormal(Caudales)

%Parametros de la distribución log normal
mu_lognormal = 7.23817;
sigma_lognormal = 0.348027;

%Grafíca de la función normal teniendo en cuenta la probabilidad de no
%excedencia
Caudal_lognormal = linspace(0,4800,40000);
P_lognormal = 1-logncdf(Caudal_lognormal,mu_lognormal,sigma_lognormal);
plot(Caudal_lognormal,P_lognormal)
hold on

%Estimación del caudal de 1000 años
Caudal_estimado_1000_years_lognormal = 4100;
Probabilidad_caudal_1000_years_lognormal = 1-logncdf(Caudal_estimado_1000_years_lognormal,mu_lognormal,sigma_lognormal);
disp(['El caudal con periodo de retorno de 1000 años de acuerdo con este ajuste es: ',num2str(Caudal_estimado_1000_years_lognormal),' m3/s'])

%Gráfica el caudal de 1000 años estimado
plot([Caudal_estimado_1000_years_lognormal,Caudal_estimado_1000_years_lognormal],[0,1])

legend('Ajuste LogNormal','Caudal de 1000 años','Location','northwest')
xlabel('Caudal (m3/s)')
ylabel('Probabilidad')
title('Ajuste LogNormal')

%%
%Ejercicio 7

%Aplicamos la prueba de bondad de ajuste de Kolmogorov-Smirnov
k = -0.05138129257569;
sigma = 420.924472877;
mu = 1252.04669596;

Caudales_vertical = Caudales';

Prueba_bondad_Ajuste_Gumbel_1 = [Caudales_vertical,gevcdf(Caudales_vertical,k,sigma,mu)];
[h,p] = kstest(Caudales_vertical,'CDF',Prueba_bondad_Ajuste_Gumbel_1,'Alpha',0.05);

% Mostrar resultados
fprintf('Valor p para el primer ajuste: %.4f\n', p);
fprintf('Valor h para el primer ajuste: %.0f\n', h);

%Respuesta
if h
    disp('Los datos no siguen una distribución gumbel.');
else
    disp('Los datos parecen seguir una distribución gumbel ya que h = 0 y el p valor es cercano a 1.');
end

Prueba_bondad_Ajuste_Gumbel_2 = [Caudales_vertical,gevcdf(Caudales_vertical,0,Sigma_regresion,mu_regresion)];
[h,p] = kstest(Caudales_vertical,'CDF',Prueba_bondad_Ajuste_Gumbel_2,'Alpha',0.05);

% Mostrar resultados
disp(' ')
fprintf('Valor p para el ajuste con regresión: %.4f\n', p);
fprintf('Valor h para el ajuste con regresión: %.0f\n', h);

%Respuesta
if h
    disp('Los datos no siguen una distribución gumbel.');
else
    disp('Los datos parecen seguir una distribución gumbel ya que h = 0 y el p valor es cercano a 1.');
end

Prueba_bondad_Ajuste_Gumbel_3 = [Caudales_vertical,gevcdf(Caudales_vertical,0,Sigma_momentos,mu_momentos)];
[h,p] = kstest(Caudales_vertical,'CDF',Prueba_bondad_Ajuste_Gumbel_3,'Alpha',0.05);

% Mostrar resultados
disp(' ')
fprintf('Valor p para el ajuste con momentos: %.4f\n', p);
fprintf('Valor h para el ajuste con momentos: %.0f\n', h);

%Respuesta
if h
    disp('Los datos no siguen una distribución gumbel.');
else
    disp('Los datos parecen seguir una distribución gumbel ya que h = 0 y el p valor es cercano a 1.');
end
%%
%Ejercicio 8

Prueba_bondad_Ajuste_lognormal = [Caudales_vertical,logncdf(Caudales_vertical,mu_lognormal,sigma_lognormal)];
[h,p] = kstest(Caudales_vertical,'CDF',Prueba_bondad_Ajuste_lognormal,'Alpha',0.05);

% Mostrar resultados
fprintf('Valor p para el ajuste lognormal: %.4f\n', p);
fprintf('Valor h para el ajuste lognormal: %.0f\n', h);

%Respuesta
if h
    disp('Los datos no siguen una distribución lognormal.');
else
    disp('Los datos parecen seguir una distribución lognormal ya que h = 0 y el p valor es cercano a 1.');
end

%%
%Ejercicio 9
%Cada año este caudal tiene una probabilidad de ocurrencia de 0.001
%entonces, asumiendo que cada evento es independiente, la probabilidad de 
%no ocurrencia de este evento en 40 años es:
P_inundacion_1000_years = 0.001;
P_no_ocurrencia = (1-P_inundacion_1000_years)^40;
P_ocurrencia = 1-P_no_ocurrencia;
fprintf('La probabilidad de que ocurra una inundación de 1000 años en los proximos 40 años es de %.4f\n',P_ocurrencia)

%%
%Ejercicio 10
%Debemos usar p(k) = n!/(k!(n-k)!)*(P(X≥x))^k*(1-P(X≥x))^(n-k) en donde n
%es el periodo considerado y k es la cantidad de veces que ocurre el
%evento. En este caso n = 100 años y k = 2 veces

n = 100;
k = 2;
p_inundacion_1000_years = 0.001;
P_ocurrencia_2_veces = factorial(n)/(factorial(k)*factorial(n-k))*(p_inundacion_1000_years)^2*(1-p_inundacion_1000_years)^(n-k);
fprintf('La probabilidad de que la inundación de 1000 años ocurra 2 veces en 100 años es: %.4f\n',P_ocurrencia_2_veces)
%%
%Funciones
function pd1 = createFitGumbel(Caudales)
%CREATEFIT    Create plot of datasets and fits
%   PD1 = CREATEFIT(CAUDALES)
%   Creates a plot, similar to the plot in the main distribution fitter
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with distributionFitter
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1
%
%   See also FITDIST.

% This function was automatically generated on 06-Mar-2024 19:16:05

% Output fitted probablility distribution: PD1

% Data from dataset "Caudales data":
%    Y = Caudales

% Force all inputs to be column vectors
Caudales = Caudales(:);

% Prepare figure
clf;
hold on;
LegHandles = []; LegText = {};


% --- Plot data originally in dataset "Caudales data"
hLine = probplot('extreme value',Caudales,[],[],'noref');
set(hLine,'Color',[0.333333 0 0.666667],'Marker','o', 'MarkerSize',6);
xlabel('Data');
ylabel('Probability')
LegHandles(end+1) = hLine;
LegText{end+1} = 'Caudales data';


% --- Create fit "Ajuste_gumbel"

% Fit this distribution to get parameter values
% To use parameter estimates from the original fit:
%     pd1 = ProbDistUnivParam('generalized extreme value',[ -0.05138129257569, 420.924472877, 1252.04669596])
pd1 = fitdist(Caudales, 'generalized extreme value');
hLine = probplot(gca,pd1);
set(hLine,'Color',[1 0 0],'LineStyle','-', 'LineWidth',2);
LegHandles(end+1) = hLine;
LegText{end+1} = 'Ajuste_gumbel';

% Adjust figure
box on;
hold off;

% Create legend from accumulated handles and labels
hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'FontSize', 9, 'Location', 'northwest');
set(hLegend,'Interpreter','none');
end

function pd1 = createFitLogNormal(Caudales)
%CREATEFIT    Create plot of datasets and fits
%   PD1 = CREATEFIT(CAUDALES)
%   Creates a plot, similar to the plot in the main distribution fitter
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with distributionFitter
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1
%
%   See also FITDIST.

% This function was automatically generated on 04-Mar-2024 21:36:55

% Output fitted probablility distribution: PD1

% Data from dataset "Caudales data":
%    Y = Caudales

% Force all inputs to be column vectors
Caudales = Caudales(:);

% Prepare figure
clf;
hold on;
LegHandles = []; LegText = {};


% --- Plot data originally in dataset "Caudales data"
[CdfY,CdfX] = ecdf(Caudales,'Function','cdf');  % compute empirical function
hLine = stairs(CdfX,CdfY,'Color',[0.333333 0 0.666667],'LineStyle','-', 'LineWidth',1);
xlabel('Data');
ylabel('Cumulative probability')
LegHandles(end+1) = hLine;
LegText{end+1} = 'Caudales data';

% Create grid where function will be computed
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);


% --- Create fit "Lognormal"

% Fit this distribution to get parameter values
% To use parameter estimates from the original fit:
%     pd1 = ProbDistUnivParam('lognormal',[ 7.238173103989, 0.3480267897421])
pd1 = fitdist(Caudales, 'lognormal');
[YPlot,YLower,YUpper] = cdf(pd1,XGrid,0.05);
hLine = plot(XGrid,YPlot,'Color',[1 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
if ~isempty(YLower)
    hBounds = plot([XGrid(:); NaN; XGrid(:)], [YLower(:); NaN; YUpper(:)],'Color',[1 0 0],...
        'LineStyle',':', 'LineWidth',1,...
        'Marker','none');
end
LegHandles(end+1) = hLine;
LegText{end+1} = 'Lognormal';
LegHandles(end+1) = hBounds;
LegText{end+1} = '95% confidence bounds';

% Adjust figure
box on;
grid on;
hold off;

% Create legend from accumulated handles and labels
hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'FontSize', 9, 'Location', 'northwest');
set(hLegend,'Interpreter','none');
end