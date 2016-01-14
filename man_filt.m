
function [V, lam] = man_filt(eegTp, eegNTp,beta)
D = (mean(eegTp,3) - mean(eegNTp,3))' * (mean(eegTp, 3) - mean(eegNTp, 3)); 
DT = 0; for i = 1:size(eegTp,3), DT = DT + eegTp(:,:,i)'*eegTp(:,:,i); end;
DT = DT / size(eegTp,3);
DT = DT - sum(eegTp,3)'*sum(eegTp,3) / (size(eegTp,3)^2);
DNT = 0; for i = 1:size(eegNTp,3), DNT = DNT + eegNTp(:,:,i)'*eegNTp(:,:,i); end;
DNT = DNT / size(eegNTp,3);
DNT = DNT - sum(eegNTp,3)'*sum(eegNTp,3) / (size(eegNTp,3)^2);
R = (beta*DT + (1-beta)*DNT);
[V, lam] = eig(D, R);
[lam,I] = sort(diag(lam),'descend');
V = V(:, I);
end