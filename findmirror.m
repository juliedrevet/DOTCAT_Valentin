function mirror = findmirror(expe)

nprac  = expe.cfg.nprac;
mirror = zeros(1,length(expe.cfg.taskid));
isub   = find(expe.cfg.taskid == 1);
isub   = isub(isub>nprac);
jsub   = find(expe.cfg.taskid == 2);
jsub   = jsub(jsub>nprac);

for i = 1:length(isub)
    iblck = isub(i);
    mirr = jsub(expe(1).cfg.condtn(jsub)==expe(1).cfg.condtn(iblck));
    if (iblck<9)
        jblck = mirr(mirr>8);
    else
        jblck = mirr((mirr<9) & (mirr>nprac));
    end
    mirror(iblck) = i;
    mirror(jblck) = i;
end 

end