signal = [];

for i=1:552
    signal = [signal yout.signals.values(:,1,i)'];
end

%x=signal(1183601:1183601+1000);
%x=signal(1209800-1800:1209800);
x=signal(40000-1000:40000+1000);


%scatterplot(x);

y = demapping(x,1,0);

% code = [-1 1 -1 1 1 1 -1 1 1 -1 -1 -1 1 1 1 1];
% 
% in = [zeros(1,length(code)) y zeros(1,length(code))];
% 
% r = [];
% 
% for j=1:length(in)-length(code)
%     
%     sum = 0;
%     
%     for i=1:length(code)
%         sum = sum + in(j+i)*code(i);
%     end
%     
%     r = [r sum];
% end
% 
% x = zeros(1, length(r));
% for i=1:length(r)
%     x(i)=i;
% end
% 
% stem(x,r)
% 
% data = [];
% 
% for i=1:length(r)
%     if r(i) > 8
%         data = [data 1];
%     elseif r(i) < -5
%         data = [data 0];
%     end
% end
