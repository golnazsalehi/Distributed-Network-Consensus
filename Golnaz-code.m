%% Setting the parameters
clear
close all

NUM_NODES = 10;
MAX_TIMESTEPS = 350;

x = rand(1,NUM_NODES) * 200; % x-coordinate
y = rand(1,NUM_NODES) * 100; % y-coordinate
dist = zeros(NUM_NODES);

GRAPH = [ ...
	0	0	1	1	0	0	0	0	0	0;
	0	0	1	1	0	0	0	0	0	0;
	1	1	0	0	1	1	0	0	0	0;
	1	1	0	0	1	1	0	0	0	0;
	0	0	1	1	0	0	1	1	0	0;
	0	0	1	1	0	0	1	1	0	0;
	0	0	0	0	1	1	0	0	1	1;
	0	0	0	0	1	1	0	0	1	1;
	0	0	0	0	0	0	1	1	0	0;
	0	0	0	0	0	0	1	1	0	0; ];
figure;
s = [1 1 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 7 8 8 8 8 9 9 10 10];
t = [3 4 3 4 1 2 5 6 1 2 5 6 3 4 7 8 3 4 7 8 5 6 9 10 5 6 9 10 8 7 8 7];
G = digraph(s,t);
plot(G,'Layout','layered')


for i=1:1:NUM_NODES
    for j=1:1:NUM_NODES
        if GRAPH(i, j) ~= 0 % link is found where i RXs from j
            dist(i, j) = sqrt( (x(i)-x(j))^2 + (y(i)-y(j))^2 );
        end
    end
end

offset = [2 3 8 1 4 5 6 7 1 0];
frequency = [1.1 .9 .75 1.3 .8 .1 .5 .43 .6 .33]*10;
%% Applying the algorithm 
c = 3*10^8; % speed of light
threshold = 0.6;

alfa = zeros(MAX_TIMESTEPS, NUM_NODES);
alfa(1, :) = ones(1, NUM_NODES);
beta = zeros(MAX_TIMESTEPS, NUM_NODES);


time_loc = zeros(MAX_TIMESTEPS, NUM_NODES); 
time_vir = zeros(MAX_TIMESTEPS, NUM_NODES); 

time_loc(1, :) = offset;
time_vir(1, :) = time_loc(1,:);

Allow2Recieve = -1*ones(1,NUM_NODES);
wait_time = zeros(1,NUM_NODES);


S = zeros(1,NUM_NODES);
A = cell(NUM_NODES,1);

l = 0;

tic
for t = 2:1:MAX_TIMESTEPS
    l = l + 1;

    time_loc(t, :) = l*frequency + offset;
    
    % Go through the graph for links to check if 'i' can send
    List = randperm(NUM_NODES);
    for i = List
      if Allow2Recieve(i) ~= -1 %check if i is receiving any messages
          continue
      end
      FLAG = -1;
      for j=1:NUM_NODES
          if i ~= j 
              if GRAPH(i, j) ~= 0 % link is found where i RXs from j
                   if (Allow2Recieve(j) == -1) && wait_time(j) < time_loc(t,j)
                       Allow2Recieve(j) = i;
                       FLAG = -2;
                   end
    
              end
         end
      end
      Allow2Recieve(i) = FLAG;
    end
    % Go through the graph for links send whatever that has been assigned
    for i=1:NUM_NODES
        if Allow2Recieve(i) > 0
            j = Allow2Recieve(i);
            d = dist(j,i);
            delay = d/c;
            contention = randi([0 ,2*15*20])/10^6;

            alfa_i = alfa(t-1,i);
            beta_i = beta(t-1,i);

            alfa_j = alfa(t-1,j);
            beta_j = beta(t-1,j);

            C_j = alfa_j * ((contention + l)*frequency(j) + offset(j)) + beta_j ; %sending time
            C_i = alfa_i * ((delay+ contention +l)*frequency(i) + offset(i)) + beta_i ; %receiving time
            wait_time(i) = C_i;
            if abs(C_i - C_j) > threshold
                A_i = A{i};
                Result = conditions_satisfied(A_i,j,S(j));
                if Result
                    [Alfa,Beta] = Complete_Update(A_i,j,C_j,C_i,alfa_i,beta_i);
                    A{i}  = [j,S(j),C_j,C_i];
                else
                    [Alfa,Beta] = Partial_Update(C_j,C_i,alfa_i,beta_i);
                    A{i,end} = [A{i,end};[j,S(j),C_j,C_i]];
                end
                alfa(t,i) = Alfa;
                beta(t,i) = Beta;
                S(i) = S(i) + 1;
            else
                alfa(t,i) = alfa(t-1,i);
                beta(t,i) = beta(t-1,i);
            end
        else
            alfa(t,i) = alfa(t-1,i);
            beta(t,i) = beta(t-1,i);
        end

    end
    time_vir(t, :) = alfa(t, :).*time_loc(t, :) + beta(t,:);
    Allow2Recieve = -1*ones(1,NUM_NODES);
end
time = toc;
fprintf('The proposed method: %f seconds.\n',time);
%%--Plotting the data
figure;
for i=1:1:NUM_NODES
    plot(1:MAX_TIMESTEPS, time_vir(:, i),'color', rand(1,3)), hold on
end
title('Virtual Time within Nodes');
xlabel('Iterations');
ylabel('Time (s)');
hold off;
%% Periodically Strongly connected
c = 3*10^8; % speed of light
threshold = 1;

alfa = zeros(MAX_TIMESTEPS, NUM_NODES);
alfa(1, :) = ones(1, NUM_NODES);
beta = zeros(MAX_TIMESTEPS, NUM_NODES);


time_loc = zeros(MAX_TIMESTEPS, NUM_NODES); 
time_vir = zeros(MAX_TIMESTEPS, NUM_NODES); 

time_loc(1, :) = offset;
time_vir(1, :) = time_loc(1,:);

Allow2Recieve = -1*ones(1,NUM_NODES);
wait_time = zeros(1,NUM_NODES);


S = zeros(1,NUM_NODES);
A = cell(NUM_NODES,1);

l = 0;

tic
for t = 2:1:MAX_TIMESTEPS
    l = l + 1;

    time_loc(t, :) = l*frequency + offset;
    
    % Go through the graph for links to check if 'i' can send
    del_links = randi([1 10]);
    TEMP_GRAPH = GRAPH; TEMP_GRAPH(del_links,:) = zeros(1,10);


    List = randperm(NUM_NODES);
    for i = List
      if Allow2Recieve(i) ~= -1 %check if i is receiving any messages
          continue
      end
      FLAG = -1;
      for j=1:NUM_NODES
          if i ~= j 
              if TEMP_GRAPH(i, j) ~= 0 % link is found where i RXs from j
                   if (Allow2Recieve(j) == -1) && wait_time(j) < time_loc(t,j)
                       Allow2Recieve(j) = i;
                       FLAG = -2;
                   end
    
              end
         end
      end
      Allow2Recieve(i) = FLAG;
    end
    % Go through the graph for links send whatever that has been assigned
    for i=1:NUM_NODES
        if Allow2Recieve(i) > 0
            j = Allow2Recieve(i);
            d = dist(j,i);
            delay = d/c;
            contention = randi([0 ,2*15*20])/10^6;

            alfa_i = alfa(t-1,i);
            beta_i = beta(t-1,i);

            alfa_j = alfa(t-1,j);
            beta_j = beta(t-1,j);

            C_j = alfa_j * ((contention + l)*frequency(j) + offset(j)) + beta_j ; %sending time
            C_i = alfa_i * ((delay+ contention +l)*frequency(i) + offset(i)) + beta_i ; %receiving time
            wait_time(i) = C_i;
            if abs(C_i - C_j) > threshold
                A_i = A{i};
                Result = conditions_satisfied(A_i,j,S(j));
                if Result
                    [Alfa,Beta] = Complete_Update(A_i,j,C_j,C_i,alfa_i,beta_i);
                    A{i}  = [j,S(j),C_j,C_i];
                else
                    [Alfa,Beta] = Partial_Update(C_j,C_i,alfa_i,beta_i);
                    A{i,end} = [A{i,end};[j,S(j),C_j,C_i]];
                end
                alfa(t,i) = Alfa;
                beta(t,i) = Beta;
                S(i) = S(i) + 1;
            else
                alfa(t,i) = alfa(t-1,i);
                beta(t,i) = beta(t-1,i);
            end
        else
            alfa(t,i) = alfa(t-1,i);
            beta(t,i) = beta(t-1,i);
        end

    end
    time_vir(t, :) = alfa(t, :).*time_loc(t, :) + beta(t,:);
    Allow2Recieve = -1*ones(1,NUM_NODES);
end
time = toc;

fprintf('periodically strongly connected: %f seconds.\n',time);
%%--Plotting the data
figure;
for i=1:1:NUM_NODES
    plot(1:MAX_TIMESTEPS, time_vir(:, i),'color', rand(1,3)), hold on
end
title('Virtual Time within Nodes');
xlabel('Iterations');
ylabel('Time (s)');
hold off;
%% Asynchronous Gossip Scheme
c = 3*10^8; % speed of light
threshold = 2;

alfa = zeros(MAX_TIMESTEPS, NUM_NODES);
alfa(1, :) = ones(1, NUM_NODES);
beta = zeros(MAX_TIMESTEPS, NUM_NODES);


time_loc = zeros(MAX_TIMESTEPS, NUM_NODES); 
time_vir = zeros(MAX_TIMESTEPS, NUM_NODES); 

time_loc(1, :) = offset;
time_vir(1, :) = time_loc(1,:);


S = zeros(1,NUM_NODES);
A = cell(NUM_NODES,1);
l = 0;

tic
for t = 2:1:MAX_TIMESTEPS
    l = l + 1;
    time_loc(t, :) = l*frequency + offset;

    % Go through the graph for links to check if 'i' can send
    List = randperm(NUM_NODES);
    for i = List
        Neighbors = logical(GRAPH(i,:));
        rows = 1:NUM_NODES;
        Neighbors = rows(Neighbors);
        L = length(Neighbors);
        Permutation = randperm(L);
        Per_of_Neighbors = Neighbors(Permutation);
        j = Per_of_Neighbors(1);

        % Sending the data

        d = dist(j,i);
        delay = d/c;
        contention = randi([0 ,2*15*20])/10^6;

        alfa_i = alfa(t-1,i);
        beta_i = beta(t-1,i);
    
        alfa_j = alfa(t-1,j);
        beta_j = beta(t-1,j);
    
         C_j = alfa_j * ((contention + l)*frequency(j) + offset(j)) + beta_j ; %sending time
         C_i = alfa_i * ((delay+ contention +l)*frequency(i) + offset(i)) + beta_i ; %receiving time
        if abs(C_i - C_j) > threshold
    
            A_i = A{i};
            Result = conditions_satisfied(A_i,j,S(j));
            if Result
                [Alfa,Beta] = Complete_Update(A_i,j,C_j,C_i,alfa_i,beta_i);
                A{i,end} = [A{i,end};[j,S(j),C_j,C_i]];
            else
                [Alfa,Beta] = Partial_Update(C_j,C_i,alfa_i,beta_i);
                A{i,end} = [A{i,end};[j,S(j),C_j,C_i]];
            end
            alfa(t,i) = Alfa;
            beta(t,i) = Beta;
            S(i) = S(i) + 1;
        else
            alfa(t,i) = alfa(t-1,i);
            beta(t,i) = beta(t-1,i);
        end

    end

    time_vir(t, :) = alfa(t, :).*time_loc(t, :) + beta(t,:);
end
time = toc;
fprintf('Asynchronous Gossip Scheme method: %f seconds.\n',time);
%%--Plotting the data
figure;
for i=1:1:NUM_NODES
    plot(1:MAX_TIMESTEPS, time_vir(:, i),'color', rand(1,3)), hold on
end
title('Virtual Time within Nodes');
xlabel('Iterations');
ylabel('Time (s)');
hold off;
%% Centralized method
tic
Mean_Freq = mean(frequency);
Mean_Offset = mean(offset);
alfa = Mean_Freq./frequency;
beta = Mean_Offset - alfa.*offset;
time = toc;
fprintf('Centralized method: %f seconds.\n',time);
%%--Plotting the data
t = 1:0.01:10;
figure;
for j=1:1:NUM_NODES
    time_vir = alfa(j) *(frequency(j)*t + offset(j)) + beta(j);
    plot(t, time_vir,'color', rand(1,3)), hold on
end
title('Virtual Time within Nodes');
xlabel('Iterations');
ylabel('Time (s)');
hold off;
%% Stochastic matrix
weights = randfixedsum(10,10,1,0,1)';
c = 3*10^8; % speed of light
threshold = 1;

alfa = zeros(MAX_TIMESTEPS, NUM_NODES);
alfa(1, :) = ones(1, NUM_NODES);
beta = zeros(MAX_TIMESTEPS, NUM_NODES);


time_loc = zeros(MAX_TIMESTEPS, NUM_NODES); 
time_vir = zeros(MAX_TIMESTEPS, NUM_NODES); 

time_loc(1, :) = offset;
time_vir(1, :) = time_loc(1,:);

Allow2Recieve = -1*ones(1,NUM_NODES);
wait_time = zeros(1,NUM_NODES);


S = zeros(1,NUM_NODES);
A = cell(NUM_NODES,1);

l = 0;

tic
for t = 2:1:MAX_TIMESTEPS
    l = l + 1;

    time_loc(t, :) = l*frequency + offset;
    
    % Go through the graph for links to check if 'i' can send
    List = randperm(NUM_NODES);
    for i = List
      if Allow2Recieve(i) ~= -1 %check if i is receiving any messages
          continue
      end
      FLAG = -1;
      for j=1:NUM_NODES
          if i ~= j 
              if GRAPH(i, j) ~= 0 % link is found where i RXs from j
                   if (Allow2Recieve(j) == -1) && wait_time(j) < time_loc(t,j)
                       Allow2Recieve(j) = i;
                       FLAG = -2;
                   end
    
              end
         end
      end
      Allow2Recieve(i) = FLAG;
    end
    % Go through the graph for links send whatever that has been assigned
    for i=1:NUM_NODES
        if Allow2Recieve(i) > 0
            j = Allow2Recieve(i);
            d = dist(j,i);
            delay = d/c;
            contention = randi([0 ,2*15*20])/10^6;

            alfa_i = alfa(t-1,i);
            beta_i = beta(t-1,i);

            alfa_j = alfa(t-1,j);
            beta_j = beta(t-1,j);

            C_j = alfa_j * ((contention + l)*frequency(j) + offset(j)) + beta_j ; %sending time
            C_j = C_j * weights(i,j);
            C_i = alfa_i * ((delay+ contention +l)*frequency(i) + offset(i)) + beta_i ; %receiving time
            wait_time(i) = C_i;
            if abs(C_i - C_j) > threshold
                A_i = A{i};
                Result = conditions_satisfied(A_i,j,S(j));
                if Result
                    [Alfa,Beta] = Complete_Update(A_i,j,C_j,C_i,alfa_i,beta_i);
                    A{i}  = [j,S(j),C_j,C_i];
                else
                    [Alfa,Beta] = Partial_Update(C_j,C_i,alfa_i,beta_i);
                    A{i,end} = [A{i,end};[j,S(j),C_j,C_i]];
                end
                alfa(t,i) = Alfa;
                beta(t,i) = Beta;
                S(i) = S(i) + 1;
            else
                alfa(t,i) = alfa(t-1,i);
                beta(t,i) = beta(t-1,i);
            end
        else
            alfa(t,i) = alfa(t-1,i);
            beta(t,i) = beta(t-1,i);
        end

    end
    time_vir(t, :) = alfa(t, :).*time_loc(t, :) + beta(t,:);
    Allow2Recieve = -1*ones(1,NUM_NODES);
end
time = toc;
fprintf('Stochastic matrix %f seconds.\n',time);
%%--Plotting the data
figure;
for i=1:1:NUM_NODES
    plot(1:MAX_TIMESTEPS, time_vir(:, i),'color', rand(1,3)), hold on
end
title('Virtual Time within Nodes');
xlabel('Iterations');
ylabel('Time (s)');
hold off;
%% varying stochastic matrix
c = 3*10^8; % speed of light
threshold = 1;

alfa = zeros(MAX_TIMESTEPS, NUM_NODES);
alfa(1, :) = ones(1, NUM_NODES);
beta = zeros(MAX_TIMESTEPS, NUM_NODES);


time_loc = zeros(MAX_TIMESTEPS, NUM_NODES); 
time_vir = zeros(MAX_TIMESTEPS, NUM_NODES); 

time_loc(1, :) = offset;
time_vir(1, :) = time_loc(1,:);

Allow2Recieve = -1*ones(1,NUM_NODES);
wait_time = zeros(1,NUM_NODES);


S = zeros(1,NUM_NODES);
A = cell(NUM_NODES,1);

l = 0;

tic
for t = 2:1:MAX_TIMESTEPS
    weights = randfixedsum(10,10,1,0,1)';
    l = l + 1;

    time_loc(t, :) = l*frequency + offset;
    
    % Go through the graph for links to check if 'i' can send
    List = randperm(NUM_NODES);
    for i = List
      if Allow2Recieve(i) ~= -1 %check if i is receiving any messages
          continue
      end
      FLAG = -1;
      for j=1:NUM_NODES
          if i ~= j 
              if GRAPH(i, j) ~= 0 % link is found where i RXs from j
                   if (Allow2Recieve(j) == -1) && wait_time(j) < time_loc(t,j)
                       Allow2Recieve(j) = i;
                       FLAG = -2;
                   end
    
              end
         end
      end
      Allow2Recieve(i) = FLAG;
    end
    % Go through the graph for links send whatever that has been assigned
    for i=1:NUM_NODES
        if Allow2Recieve(i) > 0
            j = Allow2Recieve(i);
            d = dist(j,i);
            delay = d/c;
            contention = randi([0 ,2*15*20])/10^6;

            alfa_i = alfa(t-1,i);
            beta_i = beta(t-1,i);

            alfa_j = alfa(t-1,j);
            beta_j = beta(t-1,j);

            C_j = alfa_j * ((contention + l)*frequency(j) + offset(j)) + beta_j ; %sending time
            C_j = C_j * weights(j,i);
            C_i = alfa_i * ((delay+ contention +l)*frequency(i) + offset(i)) + beta_i ; %receiving time
            wait_time(i) = C_i;
            if abs(C_i - C_j) > threshold
                A_i = A{i};
                Result = conditions_satisfied(A_i,j,S(j));
                if Result
                    [Alfa,Beta] = Complete_Update(A_i,j,C_j,C_i,alfa_i,beta_i);
                    A{i}  = [j,S(j),C_j,C_i];
                else
                    [Alfa,Beta] = Partial_Update(C_j,C_i,alfa_i,beta_i);
                    A{i,end} = [A{i,end};[j,S(j),C_j,C_i]];
                end
                alfa(t,i) = Alfa;
                beta(t,i) = Beta;
                S(i) = S(i) + 1;
            else
                alfa(t,i) = alfa(t-1,i);
                beta(t,i) = beta(t-1,i);
            end
        else
            alfa(t,i) = alfa(t-1,i);
            beta(t,i) = beta(t-1,i);
        end

    end
    time_vir(t, :) = alfa(t, :).*time_loc(t, :) + beta(t,:);
    Allow2Recieve = -1*ones(1,NUM_NODES);
end
time = toc;
fprintf('varying stochastic matrix: %f seconds.\n',time);
%%--Plotting the data
figure;
for i=1:1:NUM_NODES
    plot(1:MAX_TIMESTEPS, time_vir(:, i),'color', rand(1,3)), hold on
end
title('Virtual Time within Nodes');
xlabel('Iterations');
ylabel('Time (s)');
hold off;
%% Doubly Stochastic Matrix
weights = Doubly_Stochastic_Matrix(10);
c = 3*10^8; % speed of light
threshold = 1;

alfa = zeros(MAX_TIMESTEPS, NUM_NODES);
alfa(1, :) = ones(1, NUM_NODES);
beta = zeros(MAX_TIMESTEPS, NUM_NODES);


time_loc = zeros(MAX_TIMESTEPS, NUM_NODES); 
time_vir = zeros(MAX_TIMESTEPS, NUM_NODES); 

time_loc(1, :) = offset;
time_vir(1, :) = time_loc(1,:);

Allow2Recieve = -1*ones(1,NUM_NODES);
wait_time = zeros(1,NUM_NODES);


S = zeros(1,NUM_NODES);
A = cell(NUM_NODES,1);

l = 0;

tic
for t = 2:1:MAX_TIMESTEPS
    l = l + 1;

    time_loc(t, :) = l*frequency + offset;
    
    % Go through the graph for links to check if 'i' can send
    List = randperm(NUM_NODES);
    for i = List
      if Allow2Recieve(i) ~= -1 %check if i is receiving any messages
          continue
      end
      FLAG = -1;
      for j=1:NUM_NODES
          if i ~= j 
              if GRAPH(i, j) ~= 0 % link is found where i RXs from j
                   if (Allow2Recieve(j) == -1) && wait_time(j) < time_loc(t,j)
                       Allow2Recieve(j) = i;
                       FLAG = -2;
                   end
    
              end
         end
      end
      Allow2Recieve(i) = FLAG;
    end
    % Go through the graph for links send whatever that has been assigned
    for i=1:NUM_NODES
        if Allow2Recieve(i) > 0
            j = Allow2Recieve(i);
            d = dist(j,i);
            delay = d/c;
            contention = randi([0 ,2*15*20])/10^6;

            alfa_i = alfa(t-1,i);
            beta_i = beta(t-1,i);

            alfa_j = alfa(t-1,j);
            beta_j = beta(t-1,j);

            C_j = alfa_j * ((contention + l)*frequency(j) + offset(j)) + beta_j ; %sending time
            C_j = C_j * weights(j,i);
            C_i = alfa_i * ((delay+ contention +l)*frequency(i) + offset(i)) + beta_i ; %receiving time
            wait_time(i) = C_i;
            if abs(C_i - C_j) > threshold
                A_i = A{i};
                Result = conditions_satisfied(A_i,j,S(j));
                if Result
                    [Alfa,Beta] = Complete_Update(A_i,j,C_j,C_i,alfa_i,beta_i);
                    A{i}  = [j,S(j),C_j,C_i];
                else
                    [Alfa,Beta] = Partial_Update(C_j,C_i,alfa_i,beta_i);
                    A{i,end} = [A{i,end};[j,S(j),C_j,C_i]];
                end
                alfa(t,i) = Alfa;
                beta(t,i) = Beta;
                S(i) = S(i) + 1;
            else
                alfa(t,i) = alfa(t-1,i);
                beta(t,i) = beta(t-1,i);
            end
        else
            alfa(t,i) = alfa(t-1,i);
            beta(t,i) = beta(t-1,i);
        end

    end
    time_vir(t, :) = alfa(t, :).*time_loc(t, :) + beta(t,:);
    Allow2Recieve = -1*ones(1,NUM_NODES);
end
time = toc;
fprintf(' Doubly Stochastic Matrix %f seconds.\n',time);
%%--Plotting the data
figure;
for i=1:1:NUM_NODES
    plot(1:MAX_TIMESTEPS, time_vir(:, i),'color', rand(1,3)), hold on
end
title('Virtual Time within Nodes');
xlabel('Iterations');
ylabel('Time (s)');
hold off;

%% Functions
function Result = conditions_satisfied(A_i,j,S_j)
    Result = false;
    if isempty(A_i) == 0
        j_s = A_i(:,1);
        where = j_s == j;
        if sum(where) ~= 0
            [r,~] = size(A_i); rows = 1:r; save_rows = rows(where); m = save_rows(end);
            if S_j <= A_i(m,2)
                Result = true;
            end
        end
    end
end
%-------0-------0-------0-------0-------0
function [alfa,beta] = Partial_Update(C_j,C_i,alfa_i,beta_i)
    alfa = alfa_i;
    beta = beta_i + (C_j - C_i)/2;
end
%-------0-------0-------0-------0-------0
function [alfa,beta] = Complete_Update(A_i,j,C_j,C_i,alfa_i,beta_i)
    j_s = A_i(:,1);
    where = j_s == j;
    [r,~] = size(A_i); rows = 1:r; save_rows = rows(where); m = save_rows(end);
    
    Delta = 0;
    for n=m:r
        Delta = Delta + ( A_i(n,3)-A_i(n,4) );
    end
    k = (C_j - A_i(m,3))/(C_i - A_i(m,4) - Delta);
    alfa = alfa_i*(1 + k)/2;
    beta = (C_j - k * C_i)/2 + beta_i*(1 + k)/2;
end
%-------0-------0-------0-------0-------0