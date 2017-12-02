% Для выполнения задания за основу взят программный пакет "Robotics
% Toolbox" (http://petercorke.com)

% Задание параметров манипулятора
l = ([0.4 0.4])';
m = ([0.5 0.5])';
r = ([0.2 0.2])';
rl = ([0.005 0.005])';

% Создание объекта манипулятора
planar2l = make_robot('planar2l', l, m, r, rl);

% Задание граничных условий для расчета траектории
qs = ([0; 0]);
qf = ([pi/4; pi/2]);

% Интерполяция траектории
trjl2 = traj5(2, qs, qf);

qt = trjl2(:,:,1);
dqt = trjl2(:,:,2);
ddqt = trjl2(:,:,3);

% Определение обобщенных моментов, элементов вектора С, G и D
tau = rne(planar2l, qt, dqt, ddqt);
C = rne(planar2l, qt, dqt, zeros(numrows(tau),planar2l.n), zeros(1,3));
G = rne(planar2l, qt, zeros(numrows(tau),planar2l.n), zeros(numrows(tau),planar2l.n));
D = tau - (C + G);

function [R] = make_robot(name, l, m, r, rl)
    n = numrows(l);
    if numrows(m) ~= n || numrows(r) ~= n || numrows(rl) ~= n
            error('Mismatch of input calues dimensions');
    end
    
    % Создание объектов звеньев, составляющих манипулятор
    for i = 1:n
        L(i) = Revolute('d', 0, 'a', l(i), 'alpha', 0);
        L(i).m = m(i);
        L(i).r = ([r(i) 0 0]);
        L(i).I = [[m(i)*rl(i)^2/2 0 0];...
                  [0 m(i)*(3*rl(i)^2+l(i)^2)/12+m(i)*rl(i)^2 0];...
                  [0 0 m(i)*(3*rl(i)^2+l(i)^2)/12+m(i)*rl(i)^2]];
        %if i == 1
        %    R = SerialLink(L(i), 'name', name);
        %else
        %    r = SerialLink(L(i));
        %    R = SerialLink([R r]);
        %end
        R = SerialLink(L, 'name', name);
    end
    
    %R = SerialLink(L, 'name', name); 
    %R.base = ([[1 0 0 0];[0 0 -1 0];[0 1 0 0];[0 0 0 1]]);
end

function [tr] = traj5(in_t, in_qs, in_qf, in_dqs, in_dqf, in_ddqs, in_ddqf)
    
    n = numrows(in_qs);

    % Обнуление не заданных граничных условий
    if nargin == 3 
        in_dqs  = zeros(1,n)'; in_dqf  = zeros(1,n)';
        in_ddqs = zeros(1,n)'; in_ddqf = zeros(1,n)';
    end
    
    % интерполяция с интервалом в 10мс
    itr_n = in_t*100;
    t_step = 1/itr_n;
    
    % Соответствующая матрица из выражения (26)
    Am = [0  0  0  0  0  1;...
          0  0  0  0  1  0;...
          0  0  0  2  0  0;...
          1  1  1  1  1  1;...
          5  4  3  2  1  0;...
          20 12 6  2  0  0];
    Ami = inv(Am);
    
    % Определение полиномиальных коэффициентов для граничных условий
    % всех обобщенных координат A, скоростей dA и ускорений ddA
    for i=1:n
        A(i,:)   = Ami*[in_qs(i); in_dqs(i); in_ddqs(i); in_qf(i); in_dqf(i); in_ddqf(i)];
        dA(i,:)  = polyder(A(i,:));
        ddA(i,:) = polyder(dA(i,:));
    end
    
    % Вычисление значений интерполирующих полиномов по всей траектории
    for i=1:n
        itr = 1;
        for t=0:t_step:1
            tr(itr, i, 1) = polyval(A(i,:),t);
            tr(itr, i, 2) = polyval(dA(i,:),t);
            tr(itr, i, 3) = polyval(ddA(i,:),t);
            itr = itr + 1;
        end
    end          
end