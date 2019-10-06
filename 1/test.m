    constellation = 5;
    BitInSym = 4;                 % колличество бит на точку
    ConstName = '16APSK';

    
    data = (0:2^BitInSym-1)    % массив десятичных чисел
    bits = de2bi(data, BitInSym);   % матрица, строки которой являются
    % бинарным представлением вышеупомянутых чисел

    % .' - транспонирование
    % reshape - изменение размеров матрицы с сохранением элементов
    % reshape(..., 1, []) - хотим матрицу с 1 строкой и не паримся о
    % количестве столбцов
    
    %bits = reshape(bits(:,end:-1:1).', 1, []);
    bits = bits(:,end:-1:1);
    bits = bits.';
    bits = reshape(bits, 1, []);
    % в итоге получили поток битов, который состоит из наших бинарных чисел
   
       
    modData = mapping(bits, constellation);
    
    %> Визуализация
    scatterplot(modData)
    text(real(modData)+0.1, imag(modData), dec2bin(data))
    title(ConstName)
    axis([-2 2 -2 2])