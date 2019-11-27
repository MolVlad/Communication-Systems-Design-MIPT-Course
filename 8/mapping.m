% Отображение бит на совездие
%> @file mapping.m
% =========================================================================
%> @brief Отображение бит на созвездие
%> @param bits биты, отображаемые на созвездие (должно быть кратно числу бит
%> на точку созвездия)
%> @param constellation условное обозначение для типа модуляции:
%> [1] - BPSK, [2] - QPSK, [3] - 8PSK, [4] - 16APSK, [5] - 16QAM
%> @return symbols точки в IQ пространстве
%> @warning Созвездие должно быть нормированно по мощности на 1.
%> Расположение точек созвездия должно совпадать с заданным в задании.
%> Полезными при выполнении будут такие функции, как bi2de и reshape 
% =========================================================================
function symbols = mapping (bits, constellation)
    
    % Инициализация созвездия
    switch (constellation)
            case 1 % BPSK  
                BitInSym = 1;                 % колличество бит на точку
                points = [1 -1];
            case 2 % QPSK
                BitInSym = 2;
                points(1:2^BitInSym) = 0;
                points(1) = complex(cos(pi/4), sin(pi/4));
                points(2) = complex(cos(3*pi/4), sin(3*pi/4));
                points(3) = complex(cos(7*pi/4), sin(7*pi/4));
                points(4) = complex(cos(5*pi/4), sin(5*pi/4));
            case 3 % 8PSK
                BitInSym = 3;
                points(1:2^BitInSym) = 0;
                points(1) = complex(cos(pi/8), sin(pi/8));
                points(2) = complex(cos(3*pi/8), sin(3*pi/8));
                points(3) = complex(cos(7*pi/8), sin(7*pi/8));
                points(4) = complex(cos(5*pi/8), sin(5*pi/8));
                points(5) = complex(cos(15*pi/8), sin(15*pi/8));
                points(6) = complex(cos(13*pi/8), sin(13*pi/8));
                points(7) = complex(cos(9*pi/8), sin(9*pi/8));
                points(8) = complex(cos(11*pi/8), sin(11*pi/8));
            case 4 % 16APSK
                BitInSym = 4;
                points(1:2^BitInSym) = 0;
                points(1) = complex(2*cos(pi/4), 2*sin(pi/4));
                points(2) = complex(2*cos(7*pi/4), 2*sin(7*pi/4));
                points(3) = complex(2*cos(3*pi/4), 2*sin(3*pi/4));
                points(4) = complex(2*cos(5*pi/4), 2*sin(5*pi/4));
                points(5) = complex(2*cos(pi/12), 2*sin(pi/12));
                points(6) = complex(2*cos(23*pi/12), 2*sin(13*pi/12));
                points(7) = complex(2*cos(11*pi/12), 2*sin(11*pi/12));
                points(8) = complex(2*cos(13*pi/12), 2*sin(13*pi/12));
                points(9) = complex(2*cos(5*pi/12), 2*sin(5*pi/12));
                points(10) = complex(2*cos(19*pi/12), 2*sin(19*pi/12));
                points(11) = complex(2*cos(7*pi/12), 2*sin(7*pi/12));
                points(12) = complex(2*cos(17*pi/12), 2*sin(17*pi/12));
                points(13) = complex(cos(pi/4), sin(pi/4));
                points(14) = complex(cos(7*pi/4), sin(7*pi/4));
                points(15) = complex(cos(3*pi/4), sin(3*pi/4));
                points(16) = complex(cos(5*pi/4), sin(5*pi/4));
            case 5 % 16QAM
                BitInSym = 4;
                points(1:2^BitInSym) = 0;
                points(1) = complex(-3, 3);
                points(2) = complex(-3, 1);
                points(3) = complex(-3, -3);
                points(4) = complex(-3, -1);
                points(5) = complex(-1, 3);
                points(6) = complex(-1, 1);
                points(7) = complex(-1, -3);
                points(8) = complex(-1, -1);
                points(9) = complex(3, 3);
                points(10) = complex(3, 1);
                points(11) = complex(3, -3);
                points(12) = complex(3, -1);
                points(13) = complex(1, 3);
                points(14) = complex(1, 1);
                points(15) = complex(1, -3);
                points(16) = complex(1, -1);
    end
    
    % Нормировка по мощности
    n = sqrt(sum(points.*conj(points))/length(points));
    points=points/n;
      
    % Перевод поток битов в массив десятичных чисел
    data = reshape(bits, BitInSym, []).';
    data = data(:,end:-1:1);
    data = bi2de(data).';
    
    % Демодуляция
    l = length(data);
    symbols = zeros(1,l);
    for i=1:l
        symbols(i) = points(data(i)+1);
    end

end
