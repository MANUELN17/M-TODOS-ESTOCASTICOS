%Manuel Andrés Niño Silva
%Taller 1 Métodos estocasticos en recursos hidraulicos

%Ejercicio 1
datos_A = [1.70, 6.26, 7.56, 7.92, 0.96, 2.47, 2.55, 0.28, 1.34, 0.71, 1.66, 2.99, 8.71, 0.09, 0.62, 0.99, 10.27, 2.96, 5.54, 3.61];

% Histograma con ancho de clase 1
figure; %Crea una ventana nueva con el histograma
histogram(datos_A, 'BinWidth', 1); %Asignación de datos que apareceran en el gráfico y su ancho de clase
title('Histograma con ancho de clase 1');
xlabel('Concentración (mg/kg)');
ylabel('Frecuencia');

% Histograma con ancho de clase 2
figure;
histogram(datos_A, 'BinWidth', 2);
title('Histograma con ancho de clase 2');
xlabel('Concentración (mg/kg)');
ylabel('Frecuencia');

% Histograma con ancho de clase 5
figure;
histogram(datos_A, 'BinWidth', 5);
title('Histograma con ancho de clase 5');
xlabel('Concentración (mg/kg)');
ylabel('Frecuencia');

%Datos entre 5 y 10
datos_cumplen = datos_A(datos_A >= 5 & datos_A <= 10);
cantidad_datos_cumplen = length(datos_cumplen);
fprintf('Un total de %d valores estan entre 5 y 10 en el conjunto A',cantidad_datos_cumplen)
%% 

%Ejercicio 2
datos_B = [1.30, 17.02, 19.74, 12.01, 0.66, 1.80, 15.91, 0.62, 2.15, 2.07, 4.68, 2.74, 11.72, 0.24, 2.30, 0.52, 5.67, 3.17, 5.92, 5.03];

% Histograma con ancho de clase 1
figure;
histogram(datos_B, 'BinWidth', 1);
title('Histograma con ancho de clase 1');
xlabel('Concentración (mg/kg)');
ylabel('Frecuencia');

% Histograma con ancho de clase 2
figure;
histogram(datos_B, 'BinWidth', 2);
title('Histograma con ancho de clase 2');
xlabel('Concentración (mg/kg)');
ylabel('Frecuencia');

% Histograma con ancho de clase 5
figure;
histogram(datos_B, 'BinWidth', 5);
title('Histograma con ancho de clase 5');
xlabel('Concentración (mg/kg)');
ylabel('Frecuencia');

%Datos entre 5 y 15
datos_cumplen = datos_B(datos_B >= 5 & datos_B <= 15);
cantidad_datos_cumplen = length(datos_cumplen);
fprintf('Un total de %d valores estan entre 5 y 15 en el conjunto B',cantidad_datos_cumplen)
%%
%Ejercicio 3
% Calcular la distribución acumulativa
acum_A = [0, cumsum(datos_A)]; %Crea una lista con los valores acumulados de las muestras A
acum_B = [0, cumsum(datos_B)]; %Crea una lista con los valores acumulados de las muestras B

% Trazar la distribución acumulativa
figure;
plot(0:length(acum_A)-1, acum_A, 'b', 'LineWidth', 2); %Asigna valores al gráfico, color y grozor de línea.
hold on; %Sigue trazando el gráfico de la distribución acumulada de B sin crear un gráfico nuevo
plot(0:length(acum_B)-1, acum_B, 'r', 'LineWidth', 2);
hold off; %Detiene el trazado de gráficos

% Añadir etiquetas y título
xlabel('Número de datos');
ylabel('Suma acumulativa (mg/kg)');
title('Distribución acumulativa de las concentraciones en A y B');
legend('Datos A', 'Datos B');

% Mostrar la leyenda
legend('show');

%%
%Ejercicio 4
datos_A = [1.70, 6.26, 7.56, 7.92, 0.96, 2.47, 2.55, 0.28, 1.34, 0.71, 1.66, 2.99, 8.71, 0.09, 0.62, 0.99, 10.27, 2.96, 5.54, 3.61];
       
datos_B = [1.30, 17.02, 19.74, 12.01, 0.66, 1.80, 15.91, 0.62, 2.15, 2.07, 4.68, 2.74, 11.72, 0.24, 2.30, 0.52, 5.67, 3.17, 5.92, 5.03];

% Media
media_A = mean(datos_A);
media_B = mean(datos_B);

% Varianza
varianza_A = var(datos_A);
varianza_B = var(datos_B);

% Asimetría
asimetria_A = calcular_asimetria(datos_A); %La función calcular_asimetria aparece al final del documento
asimetria_B = calcular_asimetria(datos_B);

% Curtosis
curtosis_A = calcular_curtosis(datos_A);
curtosis_B = calcular_curtosis(datos_B);

% Cuantiles
cuantiles_A = quantile(datos_A, [0.25, 0.50, 0.75]); %muestra los diferentes cuantiles para cada grupo de datos en una misma linea
cuantiles_B = quantile(datos_B, [0.25, 0.50, 0.75]);

% Mediana
mediana_A = median(datos_A);
mediana_B = median(datos_B);

% Rango intercuartílico
rango_intercuartilico_A = iqr(datos_A);
rango_intercuartilico_B = iqr(datos_B);

% Mostrar resultados en forma organizada
disp('Parametros estadisticos para datos A:');
disp(['Media: ', num2str(media_A)]);
disp(['Varianza: ', num2str(varianza_A)]);
disp(['Asimetría: ', num2str(asimetria_A)]);
disp(['Curtosis: ', num2str(curtosis_A)]);
disp(['Cuantiles (25%, 50%, 75%): ', num2str(cuantiles_A)]);
disp(['Mediana: ', num2str(mediana_A)]);
disp(['Rango intercuartílico: ', num2str(rango_intercuartilico_A)]);

disp(' ');

disp('Parametros estadisticos para datos B:');
disp(['Media: ', num2str(media_B)]);
disp(['Varianza: ', num2str(varianza_B)]);
disp(['Asimetría: ', num2str(asimetria_B)]);
disp(['Curtosis: ', num2str(curtosis_B)]);
disp(['Cuantiles (25%, 50%, 75%): ', num2str(cuantiles_B)]);
disp(['Mediana: ', num2str(mediana_B)]);
disp(['Rango intercuartílico: ', num2str(rango_intercuartilico_B)]);

%%
%Ejercicio 5
datos_A = [1.70, 6.26, 7.56, 7.92, 0.96, 2.47, 2.55, 0.28, 1.34, 0.71, 1.66, 2.99, 8.71, 0.09, 0.62, 0.99, 10.27, 2.96, 5.54, 3.61];
       
datos_B = [1.30, 17.02, 19.74, 12.01, 0.66, 1.80, 15.91, 0.62, 2.15, 2.07, 4.68, 2.74, 11.72, 0.24, 2.30, 0.52, 5.67, 3.17, 5.92, 5.03];

% Crear diagrama de cajas para datos A
figure;
boxplot(datos_A);
xlabel('Datos A');
ylabel('Valores (mg/kg)');
title('Diagrama de Cajas y Bigotes para Datos A');

% Crear diagrama de cajas para datos B
figure;
boxplot(datos_B);
xlabel('Datos B');
ylabel('Valores (mg/kg)');
title('Diagrama de Cajas y Bigotes para Datos B');

%Determinar cantidad de datos anomalos de A
Lim_sup_A = quantile(datos_A, 0.75)+iqr(datos_A)*1.5;
Lim_inf_A = quantile(datos_A, 0.25)-iqr(datos_A)*1.5;
Datos_anomalos_A = datos_A(datos_A < Lim_inf_A | datos_A > Lim_sup_A);
Cant_Datos_Anomalos_A = length(Datos_anomalos_A);
fprintf('Hay un total de %d datos anomalos en A y son:',Cant_Datos_Anomalos_A)
disp(Datos_anomalos_A);

disp('  ')

%Determinar cantidad de datos anomalos de B
Lim_sup_B = quantile(datos_B, 0.75)+iqr(datos_B)*1.5;
Lim_inf_B = quantile(datos_B, 0.25)-iqr(datos_B)*1.5;
Datos_anomalos_B = datos_B(datos_B < Lim_inf_B | datos_B > Lim_sup_B);
Cant_Datos_Anomalos_B = length(Datos_anomalos_B);
fprintf('Hay un total de %d datos anomalos en B y son:',Cant_Datos_Anomalos_B)
disp(Datos_anomalos_B);

disp('  ')
%%
%Ejercicio 6

% Datos de las concentraciones
datos_A = [1.70, 6.26, 7.56, 7.92, 0.96, 2.47, 2.55, 0.28, 1.34, 0.71,1.66, 2.99, 8.71, 0.09, 0.62, 0.99, 10.27, 2.96, 5.54, 3.61];

% Concentración crítica
concentracion_critica = 5; % mg/kg

% Contar muestras que superan la concentración crítica
muestras_exceden_critica = sum(datos_A > concentracion_critica);

% Área del sitio
area_sitio = 8000; % El área esta en m^2

% Calcular el área aproximada que debe ser limpiada en proporción a los datos criticos
area_limpiada_aproximada = (muestras_exceden_critica / length(datos_A)) * area_sitio;

fprintf('Área aproximada a limpiar: %.2f m^2\n', area_limpiada_aproximada);


%%
%Ejercicio 7

% Datos de las concentraciones
datos_B = [1.30, 17.02, 19.74, 12.01, 0.66, 1.80, 15.91, 0.62, 2.15, 2.07, 4.68, 2.74, 11.72, 0.24, 2.30, 0.52, 5.67, 3.17, 5.92, 5.03];
% Concentración crítica
concentracion_critica = 10; % mg/kg

% Contar muestras que superan la concentración crítica
muestras_exceden_critica = sum(datos_B > concentracion_critica);

% Área del sitio
area_sitio = 8000; % m^2

% Calcular el área aproximada que debe ser limpiada en proporción a los datos criticos
area_limpiada_aproximada = (muestras_exceden_critica / length(datos_B)) * area_sitio;

fprintf('Área aproximada a limpiar: %.2f m^2\n', area_limpiada_aproximada);

%%
%Ejercicio 8

% Datos A y B
datos_A = [1.70, 6.26, 7.56, 7.92, 0.96, 2.47, 2.55, 0.28, 1.34, 0.71,1.66, 2.99, 8.71, 0.09, 0.62, 0.99, 10.27, 2.96, 5.54, 3.61];
datos_B = [1.30, 17.02, 19.74, 12.01, 0.66, 1.80, 15.91, 0.62, 2.15, 2.07, 4.68, 2.74, 11.72, 0.24, 2.30, 0.52, 5.67, 3.17, 5.92, 5.03];

% Calcular el coeficiente de correlación
coeficiente_correlacion = corrcoef(datos_A, datos_B);

% Mostrar el coeficiente de correlación
fprintf('El coeficiente de correlación entre A y B es: %.4f\n', coeficiente_correlacion(1, 2));%Muestra el dato ubicado en la posición 1,2 de matriz que resulta de calcular el coeficiente

%%

%Ejercicio 9

% Datos A y B
datos_A = [1.70, 6.26, 7.56, 7.92, 0.96, 2.47, 2.55, 0.28, 1.34, 0.71, 1.66, 2.99, 8.71, 0.09, 0.62, 0.99, 10.27, 2.96, 5.54, 3.61];
       
datos_B = [1.30, 17.02, 19.74, 12.01, 0.66, 1.80, 15.91, 0.62, 2.15, 2.07, 4.68, 2.74, 11.72, 0.24, 2.30, 0.52, 5.67, 3.17, 5.92, 5.03];

% Valor crítico
valor_A = 5;
valor_B = 10;

% Contar muestras que cumplen ambas condiciones
muestras_cumplen_condicion = sum(datos_A < valor_A & datos_B < valor_B); %suma los datos que cumplen con las condiciones dadas

% Calcular la fracción de muestras que cumplen ambas condiciones
fraccion_muestras_cumplen_condicion = muestras_cumplen_condicion / (length(datos_A))*100; %Calcula la fracción de datos (En decimal) que cumple

fprintf('La fracción de los datos que cumplen ambas condiciones es: %d%%.\n',fraccion_muestras_cumplen_condicion);
fprintf('Es decir, de las %d muestras, %d cumplen ambas condiciones\n', length(datos_A), muestras_cumplen_condicion);
%%

%Ejercicio 10

% Datos A y B
datos_A = [1.70, 6.26, 7.56, 7.92, 0.96, 2.47, 2.55, 0.28, 1.34, 0.71, 1.66, 2.99, 8.71, 0.09, 0.62, 0.99, 10.27, 2.96, 5.54, 3.61];
       
datos_B = [1.30, 17.02, 19.74, 12.01, 0.66, 1.80, 15.91, 0.62, 2.15, 2.07, 4.68, 2.74, 11.72, 0.24, 2.30, 0.52, 5.67, 3.17, 5.92, 5.03];

% Valor crítico
valor_critico_A = 5;
valor_critico_B = 10;

% Contar muestras que cumplen al menos una de las condiciones
muestras_cumplen_condicion = sum(datos_A < valor_critico_A | datos_B < valor_critico_B);

% Calcular la fracción de muestras que cumplen al menos una de las condiciones
fraccion_muestras_cumplen_condicion = muestras_cumplen_condicion / length(datos_A)*100;

fprintf('La fracción de los datos que cumplen al menos una de las condiciones es: %d%%.\n', fraccion_muestras_cumplen_condicion);
fprintf('Es decir, de las %d muestras, %d cumplen alguna de las condiciones\n', length(datos_A), muestras_cumplen_condicion);
%%
%Funciones
function asimetria = calcular_asimetria(datos)
    media = mean(datos);
    desviacion_estandar = std(datos);
    n = length(datos);
    asimetria = (1/n) * sum((datos - media).^3)/(desviacion_estandar).^3;
end

function curtosis = calcular_curtosis(datos)
    media = mean(datos);
    desviacion_estandar = std(datos);
    n = length(datos);
    curtosis =((1/n) * sum((datos - media).^4)/(desviacion_estandar).^4)-3;
end