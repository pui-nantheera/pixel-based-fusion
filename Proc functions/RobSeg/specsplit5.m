function outmap = specsplit5(W,A,map,level);
%now using new formulation for maximised mean association
%and with propagation of external edges for decision making
level =level+1;
numreg=size(W,1);
thresh1 = 0.8;
thresh2 = 0.65;
%thresh = 1.2;
W(logical(eye(size(W))))=0;

if (numreg>2);  %non-trivial cut computation
    
    d=sum(W,2);
%     [minval,badapple] = min(d);
%     
%     if (minval==0);  %method screws up sometimes if one region isolated, so do it manually
%         inds = [badapple, find([1:numreg]~=badapple)];
%         cutpoint =1;
%         crit = 0;
%         affval = 1;
%         Wsort = W(inds,inds);
%         Asort = A(inds,inds);
%     else
        D = diag(d);
        S = diag(sum(A,2));
        
        Nr = (D-W);
        Dr = (S-A);
        [V,evals]=eig(Nr,Dr); %seems faster than using eigs to get just 2 eigenvectors
        [minval,minind] = min(abs(diag(evals)));
        if( (max(V(:,minind))-min(V(:,minind)))<1e-3 );
            evals(minind,minind) = inf;
            [minval,minind] = min(abs(diag(evals)));;
        end;
        [Y,inds] = sort(V(:,minind));
        Wsort = W(inds,inds);
        Asort = A(inds,inds);
        clear A W D S Nr Dr V evals;
        [crit,affval,cutpoint] = choosecut(Wsort,Asort);
%     end;
    
%     maxassoc = max(affval,bestsofar);

%    finished = ((crit2./maxassoc)>thresh2);
    carryon = (crit<thresh2)|((crit./affval)<thresh1);  %can get divide by zero errors !!!
%    finished = ((crit2)>thresh2);
    finished = ~carryon;
    if(~finished);
        r1 = inds(1:cutpoint);
        r2 = inds(cutpoint+1:end);
        
        aff1 = Wsort(1:cutpoint,1:cutpoint);
        adj1 = Asort(1:cutpoint,1:cutpoint);
        aff2 = Wsort(cutpoint+1:end,cutpoint+1:end);
        adj2 = Asort(cutpoint+1:end,cutpoint+1:end);
        clear Asort Wsort;
        
%         cutmat = Wsort(1:cutpoint,cutpoint+1:end);
%         cutadj = Asort(1:cutpoint,cutpoint+1:end);
%         edgevec1 = sum([cutmat,max(extW(r1),0)],2);  %edge values to propagate with r1
%         denomvec = sum([cutadj,(extW(r1)>=0)],2);
%         edgevec1(denomvec>0) = edgevec1(denomvec>0)./denomvec(denomvec>0);
%         edgevec1(denomvec==0) = -1; %propagate "no external edge" markers
%         
%         edgevec2 = sum([cutmat;max(extW(r2)',0)],1); %same for r2
%         denomvec = sum([cutadj;(extW(r2)>=0)'],1);
%         edgevec2(denomvec>0) = edgevec2(denomvec>0)./denomvec(denomvec>0);
%         edgevec2(denomvec==0) = -1;
%     
%         edgevec2 = edgevec2';
        
        map1 = zeros(size(map));
        map2=map1;
        for k=1:length(r1);
            map1(map==r1(k))=k;
        end;
        for k=1:length(r2);
            map2(map==r2(k))=k;
        end;
        clear map;
        omap1=specsplit5(aff1,adj1,map1,level);
        omap2=specsplit5(aff2,adj2,map2,level);
        omap1(omap1>0) = omap1(omap1>0) + 2^level;
        outmap = omap1+omap2;
    else
        outmap = map;
        outmap(outmap>0)=1;
    end;

elseif (numreg==2);  %best cut is trivial if only two regions
    if (W(2,1)>thresh2);   %slightly arbitrary choice of when to merge two single regions!!
        outmap = map;
        outmap(outmap>0)=1;
    else
        outmap = map;
        outmap(outmap==1) = 1;
        outmap(outmap==2) = 1+2^level;
    end;
    
else  %Have been fed a single region by previous level (could avoid such calls altogether really)
    outmap = map;
    outmap(outmap>0)=1;
        
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function index = choosecut(vec);  
% % This version finds the max jump in the indicator vector
% dvec = diff(vec);
% [dumval,index] = max(dvec);
% return;

% function [breakval,index] = choosecut(affreorder);
% % This jump choose the minimum affinity in the chain as the spot to break
% % NB only works if you have calculated non-adjacent affinities - sometimes
% % the chain goes through non-adjacent links, and these would otherwise be
% % zero
% % This doesn't work on Fiedler-type ordering coz you choose cut on
% % completely different criteria to those which give rise to the ordering
% affchain = diag(affreorder,1);
% [breakval,index] = min(affchain);
% 
% return;

% function [breakval,index] = choosecut(affreorder);
% % choose based on min normalised cut
% N=size(affreorder,1)-1;
% cutvals = zeros(N,1);
% affreorder(logical(eye(size(affreorder))))=0;
% for k=1:N;
%     assocA = sum(sum(affreorder(1:k,1:k)));
%     assocB= sum(sum(affreorder(k+1:end,k+1:end))); 
%     cut = sum(sum(affreorder(1:k,k+1:end))) + sum(sum(affreorder(k+1:end,1:k)));
%     cutvals(k) = (cut./(cut+assocA)) + (cut./(cut+assocB));
% end;
% [breakval,index] = min(cutvals);
% return;

% function [breakval,index] = choosecut(affreorder);
% % choose based on min av cut
% N=size(affreorder,1)-1;
% cutvals = zeros(N,1);
% affreorder(logical(eye(size(affreorder))))=0;
% for k=1:N;
%     crossaff=nonzeros(affreorder(1:k,k+1:end));
%     if isempty(crossaff);
%         whatsupdoc=1;
%         cutvals(k) = 0;
%     else
%         cut = sum(crossaff)./length(crossaff);
%         cutvals(k) = cut;
%     end;
% end;
% [breakval,index] = min(cutvals);
% return;

function [breakval,mval,index] = choosecut(affreorder,adjace);
% choose based on mean between 
N=size(affreorder,1)-1;
mmaffvals = zeros(N,1);

affreorder(logical(eye(size(affreorder))))=0;
for k=1:N;
    crossaff=nonzeros(affreorder(1:k,k+1:end));
   
    A1 = sum(sum(affreorder(1:k,1:k)));
    A2 = sum(sum(affreorder(k+1:end,k+1:end)));
    C = sum(sum(affreorder(1:k,k+1:end)));
    N1 = sum(sum(adjace(1:k,1:k)));
    N2 = sum(sum(adjace(k+1:end,k+1:end)));
    NC = sum(sum(adjace(1:k,k+1:end)));
    if A1>0;
        A1= A1./N1;
    end;
    if A2>0;
        A2 = A2./N2;
    end;
    if C>0;
        C = C./NC;
    end;
    mmaff = max(A1,A2);
    Cvals(k) = C;
    mmaffvals(k) = mmaff;
%     if C==0;
%         cutvals(k) = 0;  % must split sections with no cross affinity
%     elseif mmaff==0;
%         cutvals(k) = inf;  %can't have both regions with no internal affinity!!
%     else
%         cutvals(k) = C./mmaff;
%     end;
end;
%[breakval,index] = min(cutvals);
[breakval,index] = min(Cvals);
mval = mmaffvals(index);
if mval==0;
    nuts=2;
end;
return;

% function [breakval,index] = choosecut(affreorder);
% % choose based on ratio of min between to max internal
% N=size(affreorder,1)-1;
% cutvals = zeros(N,1);
% affreorder(logical(eye(size(affreorder))))=0;
% for k=1:N;
%     crossaff=nonzeros(affreorder(1:k,k+1:end));
%     intaff1 = nonzeros(affreorder(1:k,1:k));
%     intaff2 = nonzeros(affreorder(k+1:end,k+1:end));
%     if isempty(crossaff);
%         whatsupdoc=1;
%         cutvals(k) = 0;
%     else
%         if(isempty([intaff1; intaff2]))
%             cutvals(k)=inf;  %can't have regions with no internal affinity!!
%         elseif isempty(crossaff);
%             cutvals(k)=0;  % must split sections with no cross affinity
%         else
%             cutvals(k) = min(crossaff)./max([intaff1; intaff2]);
%         end;
%     end;
% end;
% [breakval,index] = min(cutvals);
% return;