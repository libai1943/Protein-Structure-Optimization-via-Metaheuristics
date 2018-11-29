function fFitness = calculateFitness(fObjV)

if (fObjV>=0);
fFitness = 1./(fObjV+1);
else
fFitness = 1+abs(fObjV);
end
