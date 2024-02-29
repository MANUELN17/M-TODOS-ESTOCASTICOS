%Manuel Andrés Niño Silva
%Taller 2
%Codigo en matlab
%%
%Ejercicio 1
Media_Y = 2;
Varianza_Y = 1.5;
Media_K = exp(Media_Y+Varianza_Y/2); %Haciendo uso de la formula μ_k=exp[μ_y+(σ_y^2)/2]
Varianza_K = exp(2*Media_Y+Varianza_Y)*(exp(Varianza_Y)-1); %Usando la formula σ_k^2=exp(2μ_y+σ_y^2 )*[exp(σ_y^2 )-1]
fprintf('La media de la variable K es %f',Media_K)
disp(' ')
fprintf('La varianza de la variable K es %f',Varianza_K)

%%
%Ejercicio 2
Media = 10;
Varianza = 200;
D_Estandar = sqrt(Varianza);
Conductividad_limite = 30;
Z = (log(Conductividad_limite)-Media)/D_Estandar; %Estadarizamos para poder calcular la probabilidad como en una distribución normal estandar
P_conductividad_menor_30 = 1-normcdf(Z, 0, 1); %Calculo de la probabilidad de que la conductividad se menor a 30 con el valor estandarizado
disp(' ')
fprintf('La probabilidad de que la conductibidad en este acuifero sea mayor a 30 m/d es de %f',P_conductividad_menor_30);

%%
%Ejercicio 3
Probabilidad_conductividad_para_cada_textura =[0,0,0,0.1,0.4,0.3,0.1,0.1;
0.3,0.4,0.2,0.1,0,0,0,0;
0.1,0.3,0.3,0.2,0.1,0,0,0]; %Matriz de probabilidades de conductividad dada una cierta textura del suelo

%Probabilidad de cada textura
P_arena = 0.7;
P_arcilla = 0.1;
P_turba = 0.1;

%Vector de probabilidades de cada textura de suelo
Probabilidad_de_cada_textura =[P_arena,P_arcilla,P_turba];

%Crear una matriz vacía con el mismo tamaño de la matriz de probabilidades de conductividad para cada textura
Distribucion_probabilidades = zeros(size(Probabilidad_conductividad_para_cada_textura));

%Realiza el calculo de la probabilidad de cada combinación de textura y conductividad en todo el acuifero y lo pone en una matriz, la suma de todos los valores de la matriz es 1
%Se usa P[A⋂B]=P[A]*P[B|A] En donde A sería la textura y B es la conductividad
for i = 1:size(Probabilidad_conductividad_para_cada_textura,1)
    Distribucion_probabilidades (i,:) = Probabilidad_conductividad_para_cada_textura(i,:)*Probabilidad_de_cada_textura(i);
end

% Títulos y subtítulos de la tabla con la distribución de probabilidades
titulo_resultado = 'Distribución de probabilidades';
subtitulos_filas = {'Arena', 'Arcilla', 'Turba'};
subtitulos_columnas = {'1x10^-3', '1x10^-2', '1x10^-1', '1x10^0','1x10^1','2x10^1','5x10^1','1x10^2'};

% Crear una figura
Tabla_Probabilidades = figure;

% Crear la tabla
tabla = uitable(Tabla_Probabilidades);

% Configurar la posición y tamaño de la tabla
tabla.Position = [10 80 600 200];

% Definir los datos de la tabla (la matriz junto con los títulos y subtítulos)
datos_tabla = num2cell(Distribucion_probabilidades);
datos_tabla = [[{titulo_resultado}; subtitulos_filas'], [subtitulos_columnas; datos_tabla]];

% Establecer los datos en la tabla
tabla.Data = datos_tabla;

% Configurar el estilo de la tabla
tabla.ColumnName = {};
tabla.RowName = {};
%%
% Ejercicio 4

%Parametros de cada variable
Med_Z1 = 10;
Var_Z1 = 300;
Med_Z2 = 25;
Var_Z2 = 450;
c_correlacion = 0.7;

%Punto a
SD_Z1 = sqrt(Var_Z1);%Desviación estandar de Z1
SD_Z2 = sqrt(Var_Z2);%Desviación estandar de Z2
Covarianza = SD_Z2*SD_Z1*c_correlacion; %Se calcula usando C_(z_1 z_2 )=r_(z_1 z_2 )*S_z1*S_(z_2 )
fprintf('La covarianza entre Z1 y Z2 es %f',Covarianza);
disp(' ');

%Punto b
Valor_esperado_Y = Med_Z1+Med_Z2; %El valor esperado de la suma de de 2 variables es la suma de sus valores esperados
fprintf('El valor esperado de Y = Z1 + Z2 es la suma de las medias de estas variables, es decir %d',Valor_esperado_Y);
disp(' ');

%Punto c
Varianza_Y = Var_Z1+Var_Z2+2*Covarianza; %La varianza de la suma de 2 variables dependientes es: Var(X+Y) = Var(X)+Var(Y)+2*Cov(X,Y)
fprintf('La varianza de Y = Z1 + Z2 es igual a %f teniendo en cuenta la covarianza entre las variables',Varianza_Y);

%%
%Punto 5

%Parametros de las variables
Med_Z1 = 10;
Var_Z1 = 300;
Med_Z2 = 25;
Var_Z2 = 450;
c_correlacion = 0.7;

% a. Calcular Pr[Z1 < 30]
Z1_menor_30 = normcdf(30, Med_Z1, sqrt(Var_Z1)); %30 es el valor que queremos evaluar y luego asignamos los parametros media y desviación estandar a la función
disp(['La probabilidad Pr[Z1 < 30] es: ', num2str(Z1_menor_30)]) %La distribuciones marginales de una distribución normal bivariada son distribuciones normales monovariadas

%b. Calcular Pr[Z2 < 40]
Z2_menor_40 = normcdf(40,Med_Z2, sqrt(Var_Z2));
disp(['La probabilidad Pr[Z2 < 40] es: ', num2str(Z2_menor_40)]) %La distribuciones marginales de una distribución normal bivariada son distribuciones normales monovariadas

% c. Calcular la probabilidad Pr[Z1 + Z2 < 50]
SD_Z1 = sqrt(Var_Z1);
SD_Z2 = sqrt(Var_Z2);
Covarianza = SD_Z2*SD_Z1*c_correlacion;
Varianza_Z1_mas_Z2 = Var_Z1+Var_Z2+2*Covarianza; %Calculo de la varianza de la suma de estas variables
Media_Z1_Z2 = Med_Z2+Med_Z1; %Calculo del valor esperado de la suma
Z1_mas_Z2_menor_50 = normcdf(50,Media_Z1_Z2, sqrt(Varianza_Z1_mas_Z2)); %Se calcula como si fuera una funcion normal monovariada con los parametros de la suma de las variables
disp(['La probabilidad Pr[Z1 + Z2 < 50]  es: ', num2str(Z1_mas_Z2_menor_50)])

% d. Calcular la probabilidad Pr[Z1 < 30 ⋂ Z2 < 40]

%Para el calculo se usa la ecuación de distribución normal bivariada 
%f(x1,x2 )=e^(-Q/2)/(2πσ1σ2 √(1-r(x1x2)^2 ))
%En donde Q=1/(1-rx1x2)* [(x1-μ_1 )^2/(σ_1^2 )-2rx1x2 (x1-μ_1 )(x2-μ_2 )/(σ_1 σ_2 )+((x2-μ_2 ))/(σ_2^2 )]
% Definir los parámetros de la distribución normal bivariada
mu = [Med_Z1 Med_Z2]; % vector de medias
sigma = [Var_Z1 Covarianza; Covarianza Var_Z2]; % matriz de covarianza

% Definir los límites de integración 
Zmin = [-Inf -Inf];
Zmax = [30 40];

% Calcular la probabilidad utilizando mvncdf
Z1_menor_30_Z2_menor_40 = mvncdf(Zmin, Zmax, mu, sigma);

disp(['La probabilidad de P[Z1 < 30 ⋂ Z2 < 40] es: ', num2str(Z1_menor_30_Z2_menor_40)]);

% e. Calcular la probabilidad de Pr[Z1 < 30 ⋃ Z2 < 40]
Union_Z1_menor_30_Z2_menor_40 = Z1_menor_30 + Z2_menor_40-Z1_menor_30_Z2_menor_40;% Usamos P[AUB] = P[A]+P[B]-P[A⋂B]
disp(['La probabilidad de P[Z1 < 30 U Z2 < 40] es: ', num2str(Union_Z1_menor_30_Z2_menor_40)])