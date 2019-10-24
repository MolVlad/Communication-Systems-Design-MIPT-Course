% Обратное отображение точек созвездия
%> @file demapping.m
% =========================================================================
%> @brief Отображение бит на созвездие
%> @param symbols точки в IQ пространстве
%> @param constellation условное обозначение для типа модуляции:
%> [1] - BPSK, [2] - QPSK, [3] - 8PSK, [4] - 16APSK, [5] - 16QAM
%> @param softFlag [0] - жесткое решение, [1] - мягкое решение
%> @param SNR отношение сигнал/шум (необходимо для мягкого решения)
%> @return out биты при жестком решении, LLR при мягком решении 
% =========================================================================
function out = demapping (symbols, constellation, softFlag, SNR)
    
    % Инициализация созвездия
    switch (constellation)
            case 1 % BPSK  
                BitInSym = 1;                 % колличество бит на точку
                numberOfPoints = 2^BitInSym;
                points = [1 -1];
            case 2 % QPSK
                BitInSym = 2;
                numberOfPoints = 2^BitInSym;
                points(1:numberOfPoints) = 0;
                points(1) = complex(cos(pi/4), sin(pi/4));
                points(2) = complex(cos(3*pi/4), sin(3*pi/4));
                points(3) = complex(cos(7*pi/4), sin(7*pi/4));
                points(4) = complex(cos(5*pi/4), sin(5*pi/4));
            case 3 % 8PSK
                BitInSym = 3;
                numberOfPoints = 2^BitInSym;
                points(1:numberOfPoints) = 0;
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
                numberOfPoints = 2^BitInSym;
                points(1:numberOfPoints) = 0;
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
                numberOfPoints = 2^BitInSym;
                points(1:numberOfPoints) = 0;
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
    
    % Нормировка созвездия по мощности
    n = sqrt(sum(points.*conj(points))/length(points));
    points=points/n;
        
    if softFlag == 0
        out(1:length(symbols)) = 0;
        for i=1:length(symbols)
            minDistance = 999;
            nearestIndex = 0;
            for j=1:numberOfPoints
                distance = abs(symbols(i) - points(j));
                if  distance < minDistance
                    nearestIndex = j;
                    minDistance = distance;
                end
            end
            out(i) = nearestIndex - 1;
        end
        
        out = de2bi(out, BitInSym);
        out = out(:,end:-1:1);
        out = out.';
        out = reshape(out, 1, []);
        out = out';
    else        
        out = zeros(BitInSym,length(symbols));
        N0 = -1 / (10 ^ (SNR / 10));
        for i=1:length(symbols)
           upperSum(1:BitInSym) = 0;
           lowerSum(1:BitInSym) = 0;
           
           for j=1:numberOfPoints
               dif = symbols(i) - points(j);
               
               n = de2bi(j - 1, BitInSym);
               for k=1:BitInSym
                   if n(k) == 1
                       lowerSum(BitInSym + 1 - k) = lowerSum(BitInSym + 1 - k) + exp(abs(dif.*conj(dif)) / N0);
                   else
                       upperSum(BitInSym + 1 - k) = upperSum(BitInSym + 1 - k) + exp(abs(dif.*conj(dif)) / N0);
                   end
                end
            end
           
            div(1:BitInSym) = 0;
            for j=1:BitInSym
                div(j) = upperSum(j) / lowerSum(j);
            end
        
            out(:,i) = log(div)';
        end
        
        out = reshape(out, [], 1);
    end
end
