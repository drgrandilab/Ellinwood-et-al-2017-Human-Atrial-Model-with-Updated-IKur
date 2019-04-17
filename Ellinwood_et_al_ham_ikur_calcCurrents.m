function currents = Ellinwood_et_al_ham_ikur_calcCurrents(t,y,p)

% After running a simulation, this file calls the model using the time
% vector and state variables to calculate timecourse for ionic currents
% and other intermediates.

currents=[];
for i=1:size(t)
    if ceil(i/1000)==i/1000
        disp(['t = ',num2str(ceil(t(i)))]);
    end
    currents=[currents;Ellinwood_et_al_ham_ikur_drug_model(t(i),y(i,:),p,'currents')];
end
% end calcCurrents