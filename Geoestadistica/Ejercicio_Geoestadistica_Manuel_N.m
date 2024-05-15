%% Geoestadistica
% Manuel Andrés Niño Silva
% Métodos Estocasticos en recursos hídricos

%Se debe tener en cuenta que para cada conjunto que se quiera probar se
%debe cambiar el variograma y la ecuación para la varianza asociada al
%proceso de Kriging de acuerdo con lo especificado en cada problema.

%El Kriging ordinario de bloque aplicado sobre el área V en este código, se hace
%solamente sobre los vertices de esta región. Es por esto por lo que pueden
%variar los resultados con respecto a los presentados en el articulo, ya
%que, ellos no indican explicitamente como evaluaron la semivarianza sobre
%esta área y existen diversas formas de promediar y estimar este valor.

%% Definición de puntos

x = [60.6,	36.2; %Puntos fijos con mediciones
60.9,	17.7;];

t_original = [60.4, 37.2; %Puntos con posibilidad de nuevas medidas
60.9,   36.1;
56.3,	34.6;
76.4,	23.3;
63.9,	21.5;
77.4,   25.3;
73.9,   34.8;
69.4,	35.6;];

t_1 = [55.08, 39.21;
    69.75, 21.85;
    56.65, 39.53;
    63.12, 39.02;
    68.29, 26.52;
    65.21, 32.47;
    64.23, 32.71;
    72.79, 24.92;];

t_2 = [77.06, 21.71;
69.91, 39.21;
60.1, 38.31;
55.55, 28.39;
71.95, 33.19;
70.47, 36.19;
68.73, 26.17;
57.51, 22.75;];

t_3 = [71.42, 21.18;
74.09, 38.85;
71.56, 36.2;
73.24, 23.62;
74.34, 35.97;
60.96, 20.52;
69.91, 33.91;
77.12, 31.62;];

t_4 = [63.47,	24.88;
59.77,	33.92;
58.04,	28.12;
66.99,	28.13;
79.51,	30.62;
66.82,	20.86;
56.56,	38.47;
64.87,	27.04;];

t_5 = [74.51, 25.56;
64.48, 29.41;
79.01, 36.19;
65.98, 25.18;
65.32, 36.19;
77.61, 21.74;
65.7, 24.73;
71.08, 22.09;];

t_6 = [65.73, 28.42;
60.64, 27.25;
76.16, 37.4;
76.93, 34.41;
56.2, 29.53;
59.52, 34.31;
59.77, 27.73;
62.47, 33.17;];

t_7 = [77.3, 33.95;
65.73, 20.26;
74.79, 33.56;
66.98, 22.01;
56.63, 39.79;
74.74, 33.28;
73.51, 29.41;
62.63, 32.18;];

t_8 = [74.96, 30;
77.9, 20.92;
65.12, 22.04;
76.64, 35.12;
70.06, 26.36;
55.93, 29.02;
59.39, 32.81;
67.86, 26.02;];

t_9 = [66.04, 26.49;
77.38, 30.92;
57, 25.5;
59.38, 31.76;
77.65, 23.7;
63.36, 30.96;
62.32, 22.1;
72.07, 24.52;];

t_10 = [72.15, 26.16;
58.86, 37.31;
76.18, 39.33;
69.22, 36.92;
76.62, 36.5;
77.59, 38.12;
79.69, 36.98;
59.98, 37.59;];

vertices_cuadro = [57.5, 22.5; %Vertices del bloque V
57.5,32.5;
72.5,22.5;
72.5,32.5];

%% Puntos Originales
disp('_________________________________________________________')
disp('Puntos Originales')
disp(' ')
disp('Método de ennumeración total:')
b = kriging_1(x,t_original,vertices_cuadro,4); %Kriging con el método de ennumeración total (#1)
disp('Método de inclusión secuencial:')
r = kriging_2(x,t_original,vertices_cuadro,4); %Kriging con el método de inlcusión secuencial (#3)

%% Ejemplo 1
disp('_________________________________________________________')
disp('Ejemplo 1')
disp(' ')
disp('Método de ennumeración total:')
b_1 = kriging_1(x,t_1,vertices_cuadro,4);
disp('Método de inclusión secuencial:')
r_1 = kriging_2(x,t_1,vertices_cuadro,4);

%% Ejemplo 2
disp('_________________________________________________________')
disp('Ejemplo 2')
disp(' ')
disp('Método de ennumeración total:')
b_2 = kriging_1(x,t_2,vertices_cuadro,4);
disp('Método de inclusión secuencial:')
r_2 = kriging_2(x,t_2,vertices_cuadro,4);

%% Ejemplo 3
disp('_________________________________________________________')
disp('Ejemplo 3')
disp(' ')
disp('Método de ennumeración total:')
b_3 = kriging_1(x,t_3,vertices_cuadro,4);
disp('Método de inclusión secuencial:')
r_3 = kriging_2(x,t_3,vertices_cuadro,4);

%% Ejemplo 4
disp('_________________________________________________________')
disp('Ejemplo 4')
disp(' ')
disp('Método de ennumeración total:')
b_4 = kriging_1(x,t_4,vertices_cuadro,4);
disp('Método de inclusión secuencial:')
r_4 = kriging_2(x,t_4,vertices_cuadro,4);

%% Ejemplo 5
disp('_________________________________________________________')
disp('Ejemplo 5')
disp(' ')
disp('Método de ennumeración total:')
b_5 = kriging_1(x,t_5,vertices_cuadro,4);
disp('Método de inclusión secuencial:')
r_5 = kriging_2(x,t_5,vertices_cuadro,4);

%% Ejemplo 6
disp('_________________________________________________________')
disp('Ejemplo 6')
disp(' ')
disp('Método de ennumeración total:')
b_6 = kriging_1(x,t_6,vertices_cuadro,4);
disp('Método de inclusión secuencial:')
r_6 = kriging_2(x,t_6,vertices_cuadro,4);

%% Ejemplo 7
disp('_________________________________________________________')
disp('Ejemplo 7')
disp(' ')
disp('Método de ennumeración total:')
b_7 = kriging_1(x,t_7,vertices_cuadro,4);
disp('Método de inclusión secuencial:')
r_7 = kriging_2(x,t_7,vertices_cuadro,4);

%% Ejemplo 8
disp('_________________________________________________________')
disp('Ejemplo 8')
disp(' ')
disp('Método de ennumeración total:')
b_8 = kriging_1(x,t_8,vertices_cuadro,4);
disp('Método de inclusión secuencial:')
r_8 = kriging_2(x,t_8,vertices_cuadro,4);

%% Ejemplo 9
disp('_________________________________________________________')
disp('Ejemplo 9')
disp(' ')
disp('Método de ennumeración total:')
b_9 = kriging_1(x,t_9,vertices_cuadro,4);
disp('Método de inclusión secuencial:')
r_9 = kriging_2(x,t_9,vertices_cuadro,4);

%% Ejemplo 10
disp('_________________________________________________________')
disp('Ejemplo 10')
disp(' ')
disp('Método de ennumeración total:')
b_10 = kriging_1(x,t_10,vertices_cuadro,4);
disp('Método de inclusión secuencial:')
r_10 = kriging_2(x,t_10,vertices_cuadro,4);

%%
function [varianzas,mejor_comb,menor_varianza] = kriging_1(datos_1, datos_variables,vertices,n)
       
       %Función del método de eenumeración total
       % Se crea la función con 3 parametros de salida (Las varianzas de
       % todas las posibles combinaciones, la mejor combinación y la menor
       % varianza que es posible ontener) y 4 parametros de entrada (Los
       % puntos fijos, los puntos variables, los vertices del área sobre la
       % que se hace Kriging y el número n de puntos que se quieren añadir
       % de los N candidatos.

       comb_posibles = generar_indices(datos_variables,n);

       %La función generar indices crea todas las posibles combinaciones de
       %n puntos a partir de las N posibles candidatos es decir genera N
       %combinado n posibilidades, cada posibilidad teniendo las posiciones
       %que ocupan los n puntos probados, dentro del grupo de N candidatos.

       comb_n = datos_1; %Genera la primera combinación que solo incluye los puntos con mediciones.
       
       for i_1 = 1:length(comb_posibles) %Crea un loop que va de 1 hasta el número de combinaciones posibles de n puntos
           for j_1 = 1:size(comb_posibles,2) %Crea un segundo loop pero solo de 1 hasta el número de elementos que tiene cada combinación (Igual a n)
               a = comb_posibles(i_1,j_1);  %Toma uno de los n números dentro la combinación que se esté probando
               b = datos_variables(a,1); %Asigna la primera componenete de las coordenadas del punto que ocupa la posición a dentro del grupo de N candidatos 
               c = datos_variables(a,2); %Asigna la segunda componenete de las coordenadas del punto que ocupa la posición a dentro del grupo de N candidatos 
               comb_n(length(datos_1)+j_1,1) = b; %Añade la primera componente "b" a la combinación que será probada
               comb_n(length(datos_1)+j_1,2) = c; %Añade la segunda componente "c" a la combinación que será probada
           end
       m_d =calcular_matriz_distancias(comb_n); %Al usar esta función cálcula la matriz de distancias de la combinación que incluye los puntos fijos y los n
       %añadidos en una determinada combinación
       m_d_v = calcular_distancia_vertices(comb_n,vertices);%Cálcula la matriz de distancias entre los puntos dentro de la combinación de puntos a probar y
       %los vertices del cuadrado V de interes.
       m_s_v = calcula_semivarianza(m_d); %Cálcula la matriz de semivarianzas para la matriz de distancias teniendo en cuenta el semivariograma especificado 
       % en cada caso
       m_s_v_v = calcula_semivarianza(m_d_v); %Calcula la matriz de semivarianzas para la matriz de distancias a los vertices
       v_s_v_p_v = zeros(size(m_s_v_v,2),1); %Crea un vector de ceros
       for i_2 = 1:size(m_s_v_v,2)
           v_s_v_p_v(i_2,1) = mean(m_s_v_v(:,i_2)); %Se calcula la media de las semivarianzas de cada punto t con cada vertice, para obtener una valor 
       %medio que represente la semivarianza entre un determinado punto y el área V de ineteres para hacer Kriging
       end
       m_kriging = matriz_final(m_s_v); %Crea la matriz final para resolver el sistema de ecuaciones del Kriging ordinario
       v_kriging = vertcat(v_s_v_p_v,1); %Crea el vector de resultados del sistema de ecuaciones del Kriging Ordinario
       pesos = m_kriging\v_kriging; %Despeja el vector con los pesos (Lamda) y el termino adicional del Kriging ordinario
       multiplicacion = zeros(size(v_s_v_p_v,1),1);%Crea un vector  
       for i = 1:size(v_s_v_p_v,1)
           multiplicacion(i,1) = v_s_v_p_v(i,1)* pesos(i,1);%Multiplica cada peso con la semivarianza entre cada punto con o sin medición y el área sobre la 
       % que se quiere aplicar el Kriging
       end
       varianza_kringing = 0.18-sum(multiplicacion,1)-pesos(size(pesos,1),1); %Calcula la varianza asociada al proceso de Kriging realizado (sigma^2_k = sigma^2-sumatoria(Lamda*Semivarianza(x_i-x))-nu)
       %En donde la varianza es 0.18 dado que es el mayor valor que adopta
       %el semivariograma. Esta ecuación varia con cada problema.
       varianzas(i_1) = varianza_kringing; %Almacena las varianzas de Kriging calculadas
       comb_n = datos_1;%Resetea la combinación de puntos
       end 
     [min_valor, indice] = min(varianzas); %Delvuelve la varianza mínima obtenida y la posición que ocupa esta dentro del vector de las varianzas calculadas, esta posición coincide con la
     %posición de la combinación en la matriz de las posibles combinaciones
     menor_varianza = min_valor; %El valor de varianza de Kriging mínimo
     mejor_comb = [comb_posibles(indice,:)];%A partir de la posición de la mínima varianza de Kriging busca la combinación que ocupa dicha posición
     disp('La menor varianza posible que se puede lograr con los puntos extras es:')
     disp(menor_varianza);
     disp('Usando la combinación de puntos:')
     disp(mejor_comb)
end


%%
function [varianzas,mejor_comb,menor_varianza] = kriging_2(datos_1, datos_variables,vertices,n)

      %Función con el método de inlcusión secuencial
      %Tiene 3 parametros de salida, las varianzas calculadas, la mejor
      %combinación hallada y la varianza asociada a esa combinación. Cuenta
      %con los mismos 4 datos de entrada también.

      datos_adicionales = 1:length(datos_variables); %Crea un vector con los número de 1 hasta N (Número de puntos candidatos)
      mejor_comb = []; %Crea un vector donde se irá almacenando la combinación que resulte en la menor varianza
      for i = 1:n %Crea un bucle que controlrá la cantidad de puntos n que se quieren utilizar dentro del grupo de N candidatos.
          for j = 1:length(datos_adicionales)%Crea un bucle para probar los N puntos disponibles y ver cual tiene la menor varianza para luego
      %añadir este punto al vector datos_1 y repetir el bucle con los puntos restantes.
             combo = [datos_1;datos_variables(datos_adicionales(j),:)]; %Añade el j-esimo punto a la matriz de coordenadas sobre la que se harpa Kriging Ordinario
             m_d = calcular_matriz_distancias(combo); %Crea la matriz de distancias y realiza el proceso de Kriging ordinario como en el anterior caso.
             m_d_v = calcular_distancia_vertices(combo,vertices);
             m_s_v = calcula_semivarianza(m_d);
             m_s_v_v = calcula_semivarianza(m_d_v);
             v_s_v_p_v = zeros(size(m_s_v_v,2),1);
             for i_2 = 1:size(m_s_v_v,2)
                v_s_v_p_v(i_2,1) = mean(m_s_v_v(:,i_2));
             end
             m_kriging = matriz_final(m_s_v); 
             v_kriging = vertcat(v_s_v_p_v,1); %Añade un  1 al final del vector de los promedios de las semivarianzas entre los puntos y los vertices
             %para completar el sistema de ecuaciones
             pesos = m_kriging\v_kriging;
             multiplicacion = zeros(size(v_s_v_p_v,1),1);
             for i_1 = 1:size(v_s_v_p_v,1)
                multiplicacion(i_1,1) = v_s_v_p_v(i_1,1)* pesos(i_1,1);
             end
             varianza_kringing = 0.18-sum(multiplicacion,1)-pesos(size(pesos,1),1);
             varianzas(j) = varianza_kringing;
          end
          indice_minimo = find(varianzas == min(varianzas)); %Encuentra la posición ocupada por el menor valor de varianza calculado añadiendo un punto más,
          %esto nos llevará a su vez a elegir el mejor punto dentro de las N posibilidades y añadirlo al conjunto datos_1
          if i == n
              menor_varianza = min(varianzas); %Cuando ya se hayan añadido los n puntos deseados al conjunto datos_1 guardará la menor varianza en esta iteración
          %que es el valor que al final nos interesa.
          end
          varianzas = [];%Resetea el vector de varianzas para empezar la siguiente iteración
          datos_1 = [datos_1;datos_variables(datos_adicionales(indice_minimo),:)];%Añade el punto que presenta la menor varianza al conjunto datos_1, el cual incluye
          %los puntos fijos. Este proceso se repite hasta añadir los n puntos deseados
          mejor_comb(i) = datos_adicionales(indice_minimo);%Añade la posición del mejor punto de la iteración para ir formando la mejor combinación de n puntos posible
          datos_adicionales(indice_minimo) = [];% Elimina la posición del punto elegido para que no se tenga en cuenta en la siguiente iteración
      end
     disp('La menor varianza posible que se puede lograr con los puntos extras es:')
     disp(menor_varianza);
     disp('Usando la combinación de puntos:')
     disp(mejor_comb)
end

%%
function gamma = variograma(h) 
    % Define el variograma asociado al problema abordado en el articulo.
    % Este se debe cambiar de acuerdo al problema abordado.
    if h == 0
        gamma = 0; 
    elseif h > 0 && h < 40
        gamma = 0.08+0.1*((3/2)*(h/40)-(1/2)*(h/40)^3);
    elseif h >= 40
        gamma =0.18 ; 
    end
end

%%
function combinaciones = generar_indices(puntos,n)
        % Esta función se encarga de crear todos los grupos posibles con n
        % elementos usando los N puntos candidatos.
        N = length(puntos);
        if n <= 0 || n > N
        error('El valor de n debe ser un entero positivo menor o igual que la longitud del vector de datos.');
        end
        indice = combnk(1:N, n);
        combinaciones = flipud(indice);
        %disp(combinaciones)
end   

%%
function matriz_distancias = calcular_matriz_distancias(datos)
    %Esta función se encarga de calcular la matriz de distancias asociada a
    %una matriz con las coordenadas de un determinado conjunto de puntos.
    num_puntos = size(datos, 1);  % Número de puntos que se estan evaluando 
    matriz_distancias = zeros(num_puntos);  % Crea la matriz de distancias

    % Calcular distancias entre cada par de puntos
    for i = 1:num_puntos
        for j = 1:num_puntos
            % Calcular la distancia euclidiana entre los puntos i y j
            matriz_distancias(i, j) = sqrt((datos(i, 1) - datos(j, 1))^2 + (datos(i, 2) - datos(j, 2))^2);
        end
    end

    % Mostrar la matriz de distancias
    %disp('Matriz de Distancias:')
    %disp(matriz_distancias);
end

%%
function matriz_distancia_vertices = calcular_distancia_vertices(puntos,vertices)

         %Esta función crea una matriz con las distancias entre cada punto
         %candidato y los vertices del área sobre la que se desea realizar
         %el Kriging

         num_vertices = size(vertices,1);
         num_puntos = size(puntos,1);
         matriz_distancia_vertices = zeros(num_vertices,num_puntos); %Crea una matriz con un número de filas igual al número vertices y con un número de 
         %columna igual al número de puntos evaluados.

         for i = 1:num_vertices %Crea un bucle que cálcula la distancia entre puntos y vertices y los almacena en la matriz dispuesta para esto.
             for j = 1:num_puntos
                 matriz_distancia_vertices(i,j) = sqrt((vertices(i,1)-puntos(j,1))^2+(vertices(i,2)-puntos(j,2))^2);
             end
         end
%disp('Matriz de Distancias a vertices:')         
%disp(matriz_distancia_vertices);         
end

%%
function matriz_semivarianza = calcula_semivarianza(x)
%Crea una matriz de tamaño identico al de la matriz de distancias y cálcula
%la semivarianza asociada a cada distancia de acuerdo con el varograma que
%este definido.
 for i = 1:size(x,1)
    for j = 1:size(x,2)
        matriz_semivarianza(i,j) = variograma(x(i,j));
    end
 end
end

%%
function matriz_kriging = matriz_final(matriz_semivarianza)

%Esta función agrega algunos ceros y unos a la matriz de semivarianzas, con
%el fin de completar la matriz asociada al sistema de ecuaciones del
%Kriging ordinario.

for i = 1:size(matriz_semivarianza,1)+1
    for j = 1:size(matriz_semivarianza,2)+1
        if i <= size(matriz_semivarianza,1) && j<= size(matriz_semivarianza,2)
            matriz_kriging(i,j) = matriz_semivarianza(i,j);
        elseif i <= size(matriz_semivarianza,1) && j > size(matriz_semivarianza,2)
            matriz_kriging(i,j) = 1;
        elseif i > size(matriz_semivarianza,1) && j <= size(matriz_semivarianza,2)
            matriz_kriging(i,j) = 1;
        elseif i > size(matriz_semivarianza,1) && j > size(matriz_semivarianza,2)
            matriz_kriging(i,j) = 0;
        end
    end
end
end

      