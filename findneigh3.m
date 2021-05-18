function [ neighbors,sizeneigh,X,checked  ] = findneigh3(neighbors,sizeneigh,indexpart,connect,listaparticles,checked,X,indexagglo)
% FINDNEIGH  busca los vecinos de los vecinos de una particula inicial  dentro de un pool de celulas dado    
%   Detailed explanation goes here

while sizeneigh~=0
    checked(indexpart,1)=boolean(0);                         % etiquetarlo como chequeado (ETIQUETA 0) de la lista particular checked
    X(indexpart,11)= indexagglo;                           % asignarle en la matrix X el numero de agregado al que corresponde
    neighbors(1,sizeneigh)=0;                                % borrar el elemento chequeado de la lista de vecinos POR EL FINAL DEL VECTOR
    sizeneigh=sizeneigh-1;                                 % reducir en una unidad la variable "tamaño del vector vecinos"
    new=listaparticles(connect(indexpart,:)==1);           % almacenar los vecinos detectados en un vector
    diam=max(X(:,9));
 
    if ~isempty(new) && any(checked(new))                                      % si la particula considerada tiene vecinos
        E=[indexpart*ones(size(new,2),1) new'];  % Creo la matriz de índices E, donde el primer indice es el índice del elemento considerado y el segundo el indice de cada uno de sus vecinos  
        [ ~,~,y1,y2] = buildcontactvectors4(X(:,1:3),X(:,[4:6 7]),E); % Obtengo los vectores que me dan los puntos a la mínima distancia entre los esferocilindros
        notconnected=sqrt(sum((y1-y2).^2,2))>1.15*diam;  % Si la distancia entre estos puntos es mayor que un 15% el mayor diámetro, los rechazo como unidos
        if any(notconnected)
            connect(sub2ind(size(connect),E(notconnected,1),E(notconnected,2)))=0;
            new(notconnected)=[]; % <<<< ===== 
        end
        if ~isempty(new)
            if numel(new)>1
                neighbors(1,sizeneigh+1:sizeneigh+size(new,2))= listaparticles(connect(indexpart,:)==1);  % añadir a la lista de "vecinos" los vecinos de la particula por detras
            else
                neighbors(1,sizeneigh+1)= listaparticles(connect(indexpart,:)==1);  % añadir a la lista de "vecinos" los vecinos de la particula por detras
            end
            sizeneigh=sizeneigh+size(new,2);                                            % incrementar el contador de vecinos pendiente por leer con los nuevos vecinos detectados
            indexpart=neighbors(1,sizeneigh);                                             % asigno un nuevo target parar hacer recurrencia con la funcion
            [ neighbors,sizeneigh,X,checked  ] = findneigh3(neighbors,sizeneigh,indexpart,connect,listaparticles,checked,X,indexagglo);
        end
        
    else                                                   % si la particula no tiene más vecinos, entonces paso al siguiente de la lista de vecinos del agregado si esta no esta vacía
        if sizeneigh~=0                                     % si la lista de vecinos generales NO ESTÁ VACÍA
            indexpart=neighbors(1,sizeneigh);
            [ neighbors,sizeneigh,X,checked  ] = findneigh3(neighbors,sizeneigh,indexpart,connect,listaparticles,checked,X,indexagglo);
        else
            break    % si está vacío, salir
        end
    end
end


end

