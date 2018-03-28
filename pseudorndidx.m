function idx = pseudorndidx(ntrl,pseudornd)
% generate indexing to control balancing constraints

idx = ones(1,ntrl);
while abs(diff(hist(idx,1:2)))~=0
    idx = ones(1,ntrl);
    idx(1:pseudornd) = ceil(2*rand(1,pseudornd));
    for i = (pseudornd+1):ntrl
        iidx = ceil(2*rand(1));
        if sum(idx((i-pseudornd):(i-1))==iidx)>=pseudornd
            idx(i) = 3-iidx;
        else
            idx(i) = iidx;
        end
    end
end
end