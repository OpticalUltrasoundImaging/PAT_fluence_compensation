% me is nx2 molar extinction matrx, p0 is nx1 spectrum
function output = linear_unmixing(p0,me)
coe = lsqnonneg(me,squeeze(p0(:)));
output = coe(1,1)/sum(coe);

end