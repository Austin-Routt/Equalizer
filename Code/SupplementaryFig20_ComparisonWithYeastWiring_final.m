% This script generates Supplementary Figure 20. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

copynumber = 1:5:1000;

miRNA = 0.303;
RISC = 1.7e+05;

count = 0;

sbioloadproject ../Models/Equalizer_model m1
csObj = m1.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 1E7 ; 

extracell_inducer = 1;
dox_influx = 0.156.*extracell_inducer;
min_sen = zeros(1,length(miRNA));

POI = zeros(1,length(copynumber));
TetR_Eq = zeros(1,length(copynumber));

for r = 1:length(RISC)
    for m = 1:length(miRNA)
        miRNAdissociation = miRNA(m);

        for i = 1:length(copynumber)
             for j = 1:length(extracell_inducer)
                namevObj1 = strcat('v1_',num2str(count));
                vObj1 = addvariant(m1,namevObj1);
                addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});

                vObj = [vObj1];

                set(m1.Reaction(11).KineticLaw.Parameters,'Value',dox_influx(j));
                set(m1.Reaction(13).KineticLaw.Parameters,'Value',miRNAdissociation);
                set(m1.species(10),'InitialAmount',dox_influx(j)/3.33e-4);
                set(m1.species(13),'InitialAmount',RISC);

                simdata = sbiosimulate(m1,csObj,vObj);
                [~, stateData] = selectbyname(simdata, 'Cell.POI');
                POI(i) = stateData(end);
                Eq_POI = stateData;
                [~, stateData] = selectbyname(simdata, 'Cell.TetR');
                TetR_Eq(i) = stateData(end);
                count = count + 1;
             end
        end
        
        % Numerical differentiation to get local log sensitivity
        dPOIdCN_Eq = zeros(size(POI));

        for j = 1:length(copynumber)
            if j == 1
                dPOIdCN_Eq(j) = (-POI(j+2)+4*POI(j+1)-3*POI(j))/2*(j/POI(j));
            elseif j == length(copynumber)
                dPOIdCN_Eq(j) = (3*POI(j)-4*POI(j-1)+POI(j-2))/2*(j/POI(j));
            else 
                dPOIdCN_Eq(j) = (POI(j+1)-POI(j-1))/2*(j/POI(j));
            end
        end
        min_sen(m) =  mean(dPOIdCN_Eq);
    end
end
figure(1) 
compensation_Eq = 1./dPOIdCN_Eq;
plot(copynumber(2:end),compensation_Eq(2:end),'LineWidth',3)
hold on

sbioloadproject ../Models/Equalizer_with_Yeast_Wiring_model  m1
csObj = m1.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 1E7 ; 

extracell_inducer = [1]; 
dox_influx = 0.156.*extracell_inducer; 

POI = zeros(length(copynumber),length(extracell_inducer));

for i = 1:length(copynumber)
     for j = 1:length(extracell_inducer)
        namevObj1 = strcat('v1_',num2str(count));
        vObj1 = addvariant(m1,namevObj1);
        addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});
        
        vObj = [vObj1];
        
        set(m1.Reaction(11).KineticLaw.Parameters,'Value',dox_influx(j));
        set(m1.Reaction(13).KineticLaw.Parameters,'Value',miRNA);
        set(m1.species(10),'InitialAmount',dox_influx(j)/3.33e-4);
        set(m1.species(13),'InitialAmount',RISC);
        
        simdata = sbiosimulate(m1,csObj,vObj);
        [~, stateData] = selectbyname(simdata, 'Cell.POI');
        POI(i) = stateData(end);
        Eq_POI = stateData;
        [~, stateData] = selectbyname(simdata, 'Cell.TetR');
        TetR_Eq(i) = stateData(end);
        count = count + 1;
     end
end

%Numerical differentiation to get local log sensitivity
dPOIdCN_yeast = zeros(size(POI));

for j = 1:length(copynumber)
   for k = 1:length(extracell_inducer)
       if j == 1
       dPOIdCN_yeast(j,k) = (-POI(j+2,k)+4*POI(j+1,k)-3*POI(j,k))/2*(j/POI(j,k));
       elseif j == length(copynumber)
       dPOIdCN_yeast(j,k) = (3*POI(j,k)-4*POI(j-1,k)+POI(j-2,k))/2*(j/POI(j,k));
       else 
       dPOIdCN_yeast(j,k) = (POI(j+1,k)-POI(j-1,k))/2*(j/POI(j,k));
       end
   end
end
compensation_yeast = 1./dPOIdCN_yeast;
plot(copynumber(2:end),compensation_yeast(2:end),'LineWidth',3)

sbioloadproject ../Models/Equalizer_with_Yeast_Wiring_TetRInhibitor_model  m1
csObj = m1.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 1E7 ; 

inhibition = [1e-9, 1e-8, 1e-7];
dox_influx = 0.156.*extracell_inducer; 

POI = zeros(length(inhibition),length(copynumber),length(extracell_inducer));
for k = 1:length(inhibition)
    for i = 1:length(copynumber)
         for j = 1:length(extracell_inducer)

            namevObj1 = strcat('v1_',num2str(count));
            vObj1 = addvariant(m1,namevObj1);
            addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});

            vObj = [vObj1];

            set(m1.Reaction(11).KineticLaw.Parameters,'Value',dox_influx(j));
            set(m1.Reaction(13).KineticLaw.Parameters,'Value',miRNA);
            set(m1.Reaction(23).KineticLaw.Parameters,'Value',inhibition(k));
            set(m1.species(10),'InitialAmount',dox_influx(j)/3.33e-4);
            set(m1.species(13),'InitialAmount',RISC);

            simdata = sbiosimulate(m1,csObj,vObj);
            [~, stateData] = selectbyname(simdata, 'Cell.POI');
            POI(k,i,j) = stateData(end);
            Eq_POI = stateData;
            count = count + 1;
         end
    end
end

% Numerical differentiation to get local log sensitivity
dPOIdCN_yeast2 = zeros(size(POI));
for i = 1:length(inhibition)
    for j = 1:length(copynumber)
        for k = 1:length(extracell_inducer)
            if j == 1
                dPOIdCN_yeast2(i,j,k) = (-POI(i,j+2,k)+4*POI(i,j+1,k)-3*POI(i,j,k))/2*(j/POI(i,j,k));
            elseif j == length(copynumber)
                dPOIdCN_yeast2(i,j,k) = (3*POI(i,j,k)-4*POI(i,j-1,k)+POI(i,j-2,k))/2*(j/POI(i,j,k));
            else 
                dPOIdCN_yeast2(i,j,k) = (POI(i,j+1,k)-POI(i,j-1,k))/2*(j/POI(i,j,k));
            end
        end
    end
end

compensation_yeast2 = 1./dPOIdCN_yeast2;
for k = 1:length(inhibition)

    plot(copynumber(2:end),compensation_yeast2(k,2:end),'LineWidth',3)
    
end

legend({'Equalizer','Equalizer\_Yeast\_Wiring','Equalizer\_Yeast\_Wiring\_TetRInhibitor','Equalizer\_Yeast\_Wiring\_TetRInhibitor','Equalizer\_Yeast\_Wiring\_TetRInhibitor','Equalizer\_Yeast\_Wiring\_TetRInhibitor'},'Location','northeast','FontSize', 20); 
legend box off 
ylim([0 30])
xt = get(gca, 'XTick');
ax = gca;
ax.YAxis.Color = 'k';
ax.YAxis.LineWidth = 2;
ax.YAxis.FontSize = 18;
ax.XAxis.Color = 'k';
ax.XAxis.LineWidth = 2;
ax.XAxis.FontSize = 18;
xlabel('DNA copy number')
ylabel('Predicted gene dosage compensation')
box off
xticks([1 100 200 300 400 500])
hold off

save('../SimulationData/SupplementaryFig20_Eq_YeastTopology_comparison.mat','compensation_Eq','compensation_yeast','compensation_yeast2','copynumber','extracell_inducer','inhibition')
