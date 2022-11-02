function h = maxhorsign(signrestriction)


for ii = 1 : size(signrestriction,2)
    tmp    = (signrestriction{ii});
    index  = strfind(tmp,','); 
    eval(['tmp1 = ' tmp( index(1)+1 : index(2)-1 ) ';']);
    hh(ii) = max(tmp1);
end

h = max(hh);
