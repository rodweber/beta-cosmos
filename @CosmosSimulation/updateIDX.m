function updateIDX(this, meanAnomalyFromAN)
%% Update property OrbitSectionNow.
%_____________________________________________________________________
%
% Details here.
%_____________________________________________________________________

%old: this.IDX = find(this.OrbitSections >= meanAnomalyFromAN, 1, 'first');
this.OrbitSectionNow = find(this.OrbitSections >= meanAnomalyFromAN,1,'first');

this.incrementIDX();

end % Function CosmosSimulation.updateIDX.
