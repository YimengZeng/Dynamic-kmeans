function W = ComputePearsonCorrelations(Y)

Xo = Y'; clear Y
% X = zeros(size(Xo,1),size(Xo,2)-8);
% global_mean = mean(Xo);
% global_mean = global_mean(9:end);
% global_mean = global_mean - mean(global_mean);
% global_mean = global_mean./norm(global_mean);
% xd = 1:length(global_mean);
% xd = xd'./norm(xd);
% xm = ones(length(global_mean),1);
% xm = xm./norm(xm);
% D = [xm xd global_mean'];
% 
% %Normalization
% for r = 1:size(X,1)
%     x = Xo(r,9:end);
%     beta_hat = D\x';
%     x = x - (D*beta_hat)';
%     X(r,:) = x;
% end
% X = Xo(:,9:end);
% C = X*X'/(size(X,1)-1);
% D = sqrt(diag(C));
% inv_D = sparse(diag(1./D));
% W_tmp = inv_D*C;
% W = W_tmp*inv_D;

X = Xo(:,9:end);
W = corr(X');