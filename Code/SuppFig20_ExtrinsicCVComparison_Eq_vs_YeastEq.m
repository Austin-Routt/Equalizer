% This script generates Supplementary Figure 20. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

miRNAdissociation = 0.303;
RISC = 1.7e5;

celltot = 10000; %Assume a total of 10000 cells
copynumber_cell = round(gamrnd(0.5716,120,celltot,1)); %copynumber in each cell,std of 1.522 fit from PGK A9 well on 20190419
copynumber_cell = copynumber_cell(copynumber_cell>0);
copynumber_cell_sorted = sort(copynumber_cell);
expression_vector=zeros(celltot,1);

sbioloadproject ../Models/Equalizer_model m1
csObj = m1.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 1e7;

extracell_inducer = [0 0.1 0.5 1 5 10 50 100];
dox_influx = 0.156.*extracell_inducer;

copynumber = unique(copynumber_cell_sorted);
POI = zeros(length(copynumber),length(extracell_inducer));

Eq_CV = zeros(length(extracell_inducer),1);
Eq_mean = zeros(length(extracell_inducer),1);
Eq_median = zeros(length(extracell_inducer),1);
count = 0;
for i = 1:length(copynumber)
    for j = 1:length(extracell_inducer)
        namevObj1 = strcat('v1_',num2str(count));
        vObj1 = addvariant(m1,namevObj1);
        addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});
        
        set(m1.Reaction(11).KineticLaw.Parameters,'Value',dox_influx(j));
        set(m1.Reaction(13).KineticLaw.Parameters,'Value',miRNAdissociation);
        set(m1.species(10),'InitialAmount',dox_influx(j)/3.33e-4);
        set(m1.species(13),'InitialAmount',RISC);
        
        simdata = sbiosimulate(m1,csObj,vObj1);
        [~, stateData] = selectbyname(simdata, 'Cell.POI');
        POI(i,j) = stateData(end);
        count = count + 1;
    end
end

for i = 1:length(extracell_inducer)
    explvl = POI(:,i)';
    for j = 1:length(copynumber_cell)
        if copynumber_cell(j)==0
            expression_vector(j)=0;
        else
            expression_vector(j) = explvl(copynumber==copynumber_cell(j));
        end
    end
    Eq_CV(i) = std(expression_vector(:))/mean(expression_vector)*100;
    Eq_mean(i) = mean(expression_vector);
    Eq_median(i) = median(expression_vector);
end

sbioloadproject ../Models/Equalizer_with_Yeast_Wiring_model  m1
csObj = m1.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 1e7;

POI = zeros(length(copynumber),length(extracell_inducer));

Eq_CV_yeast = zeros(length(extracell_inducer),1);
Eq_mean_yeast = zeros(length(extracell_inducer),1);
Eq_median_yeast = zeros(length(extracell_inducer),1);
count = 0;
for i = 1:length(copynumber)
    for j = 1:length(extracell_inducer)
        namevObj1 = strcat('v1_',num2str(count));
        vObj1 = addvariant(m1,namevObj1);
        addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});
        
        set(m1.Reaction(11).KineticLaw.Parameters,'Value',dox_influx(j));
        set(m1.Reaction(13).KineticLaw.Parameters,'Value',miRNAdissociation);
        set(m1.species(10),'InitialAmount',dox_influx(j)/3.33e-4);
        set(m1.species(13),'InitialAmount',RISC);
        
        simdata = sbiosimulate(m1,csObj,vObj1);
        [~, stateData] = selectbyname(simdata, 'Cell.POI');
        POI(i,j) = stateData(end);
        count = count + 1;
    end
end

for i = 1:length(extracell_inducer)
    explvl = POI(:,i)';
    for j = 1:length(copynumber_cell)
        if copynumber_cell(j)==0
            expression_vector(j)=0;
        else
            expression_vector(j) = explvl(copynumber==copynumber_cell(j));
        end
    end
    Eq_CV_yeast(i) = std(expression_vector(:))/mean(expression_vector)*100;
    Eq_mean_yeast(i) = mean(expression_vector);
    Eq_median_yeast(i) = median(expression_vector);
end

sbioloadproject ../Models/Equalizer_with_Yeast_Wiring_TetRInhibitor_model  m1
csObj = m1.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 1e7;
inhibition = [1e-9, 1e-8, 1e-7];
POI = zeros(length(inhibition),length(copynumber),length(extracell_inducer));

Eq_CV_yeast_2 = zeros(length(inhibition),length(extracell_inducer),1);
Eq_mean_yeast_2 = zeros(length(inhibition),length(extracell_inducer),1);
Eq_median_yeast_2 = zeros(length(inhibition),length(extracell_inducer),1);
count = 0;
for k = 1:length(inhibition)
    for i = 1:length(copynumber)
        for j = 1:length(extracell_inducer)
            namevObj1 = strcat('v1_',num2str(count));
            vObj1 = addvariant(m1,namevObj1);
            addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});

            set(m1.Reaction(11).KineticLaw.Parameters,'Value',dox_influx(j));
            set(m1.Reaction(13).KineticLaw.Parameters,'Value',miRNAdissociation);
            set(m1.Reaction(23).KineticLaw.Parameters,'Value',inhibition(k));

            set(m1.species(10),'InitialAmount',dox_influx(j)/3.33e-4);
            set(m1.species(13),'InitialAmount',RISC);

            simdata = sbiosimulate(m1,csObj,vObj1);
            [~, stateData] = selectbyname(simdata, 'Cell.POI');
            POI(k,i,j) = stateData(end);
            count = count + 1;
        end
    end
end

for k = 1:length(inhibition)
    for i = 1:length(extracell_inducer)
        explvl = POI(k,:,i)';
        for j = 1:length(copynumber_cell)
            if copynumber_cell(j)==0
                expression_vector(j)=0;
            else
                expression_vector(j) = explvl(copynumber==copynumber_cell(j));
            end
        end
        Eq_CV_yeast_2(k,i) = std(expression_vector(:))/mean(expression_vector)*100;
        Eq_mean_yeast_2(k,i) = mean(expression_vector);
        Eq_median_yeast_2(k,i) = median(expression_vector);
    end
end

figure(14)

plot(extracell_inducer,Eq_CV,'--bs','LineWidth',4,'MarkerSize',15,'MarkerEdgeColor', 'b')
hold on 
plot(extracell_inducer,Eq_CV_yeast,'--ks','LineWidth',4,'MarkerSize',15,'MarkerEdgeColor', 'k')
plot(extracell_inducer,Eq_CV_yeast_2(1,:),'-gs','LineWidth',4,'MarkerSize',15,'MarkerEdgeColor', 'k')
plot(extracell_inducer,Eq_CV_yeast_2(2,:),'--gs','LineWidth',4,'MarkerSize',15,'MarkerEdgeColor', 'k')
plot(extracell_inducer,Eq_CV_yeast_2(3,:),'-.gs','LineWidth',4,'MarkerSize',15,'MarkerEdgeColor', 'k')

xt = get(gca, 'XTick');
ax = gca;
ax.YAxis.Color = 'k';
ax.YAxis.LineWidth = 2;
ax.YAxis.FontSize = 20;
ax.XAxis.Color = 'k';
ax.XAxis.LineWidth = 2;
ax.XAxis.FontSize = 20;
ax.TickLength = [0.025, 0.025];
ax.XAxis.FontSize = 18;
ax.Title.FontSize = 20;
xlabel('[doxycycline] (ng/mL)')
ylabel('Cell-to-cell variation (%CV)')
title('Total cell-to-cell variation')
legend_list = {'Equalizer','Equalizer\_multipromoter'};
legend(legend_list,'Location','northeast','FontSize', 16); 
legend box off 
pbaspect([1 1 1])
box off
hold off

save('../SimulationData/Eq_vs_Eq_multipromoter_extrinsic_noise.mat','Eq_CV','Eq_CV_yeast','Eq_CV_yeast_2','extracell_inducer','inhibition')
