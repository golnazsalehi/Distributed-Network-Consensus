function M2 = Doubly_Stochastic_Matrix(n)
    c = randfixedsum(n,n,1,0,1);
    M2 = zeros(n,n);
    tic  % proposed alternative
    for k = 1:n
    %     [~,idx] = sort(F{k});    % Used to compare same results.
        [~,idx] = sort(GRPdur(n)); % Or just idx = GRPdur(n);  ??
        idx = idx + (0:n-1)*n;     
        M2(idx) = M2(idx) + c(k);
    end
end
function p = GRPdur(n)
% -------------------------------------------------------------
% Generate a random permutation p(1:n) using Durstenfeld's 
% O(n) Shuffle Algorithm, CACM, 1964. 
% See Knuth, Section 3.4.2, TAOCP, Vol 2, 3rd Ed.
% Complexity: O(n)
% USE: p = GRPdur(10^7);
% Derek O'Connor, 8 Dec 2010.  derekroconnor@eircom.net
% -------------------------------------------------------------
      p = 1:n;                  % Start with Identity permutation
  for k = n:-1:2    
      r = 1+floor(rand*k);      % random integer between 1 and k
      t    = p(k);
      p(k) = p(r);               % Swap(p(r),p(k)).
      p(r) = t;                  
  end
 end