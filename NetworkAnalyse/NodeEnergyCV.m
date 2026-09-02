function value = NodeEnergyCV(En)
En = double(En(:));
value = std(En,0,'omitnan')/(abs(mean(En,'omitnan'))+eps);
end
