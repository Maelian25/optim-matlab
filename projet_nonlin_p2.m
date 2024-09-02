%% *Optimisation ligne de métro*

close all
clear variables
% Chargement des données

% chargement des données du problèmes
trafic = load('data_trafic');
trafic = trafic.trafic;
% taille du problème
capacite = 170;          % capacité d'une rame
nb_station = 8;         % nombre de stations
T = 60;                 % 10 minutes
%% Définition du problème
% Variables de décisions

% variables de décisions
d = optimvar('d',T,1,'Type','integer', 'LowerBound',0,'UpperBound',1);      % rame partant de 1 
quai = optimvar('quai',nb_station,T,'Type','integer', 'LowerBound',0);      % nombre de personne en attente sur le quai s à t
rame = optimvar('rame',nb_station,T,'Type','integer','LowerBound',0);       % nombre de personne dans la rame à la station s, à t
entrant = optimvar('entrant',nb_station,T,'Type','integer','LowerBound',0); % nombre de personne qui sortent à s t
sortant = optimvar('sortant',nb_station,T,'Type','integer','LowerBound',0); % nombre de personne qui rentrent à s t
% *Contraintes de conservation des flux dans les rames*

% Contraintes de conservation des flux dans les rames
contrame1 = rame(1:nb_station,1)==entrant(1:nb_station,1);    % pour s 0 à n-1
contrame2 = rame(1,1:T)==entrant(1,1:T);                      % pour t 0 à T-1
contrame3 = rame(2:nb_station,2:T) == rame(1:nb_station-1,1:T-1)+entrant(2:nb_station,2:T)+sortant(2:nb_station,2:T);
% *Contraintes de conservation des flux sur les quais*

% Contraintes de conservation des fluis sur les quais
contraintes_flux_quais = optimconstr(nb_station,T-1);
for s=1:1:nb_station
    for t=2:1:nb_station
        somme=sum(trafic(s,s+1:nb_station,t));
        contraintes_flux_quais(s,t-1) = quai(s,t) == quai(s,t-1)-entrant(s,t)+somme;
    end
end
% Contraintes de Fréquence

% Contraintes de fréquences pas de plus de rame toutes les 2 minutes
A = zeros(T-2,T);
for i=1:1:T-2
    for j=i:1:i+2
        A(i,j)=1;
    end
end
b = 2*ones(T-2,1);
consfrequence = A*d<=b;
% Contraintes entrée et sortie rame & Contrainte capacité rame

% Contrainte entrée et sortie dans une rame / contrainte capacite de la rame 
contrainte_entree_rame = optimconstr(nb_station,T);
contrainte_sortie_rame = optimconstr(nb_station,T);
contrainte_capacite_rame = optimconstr(nb_station,T);

for t=1:1:T
    for s=1:1:nb_station
        if t-s>=0
            valeur = capacite*d(t-s+1);
        else
            valeur = 0;
        end
        contrainte_entree_rame(s,t) = entrant(s,t)<= valeur;
        contrainte_sortie_rame(s,t) = sortant(s,t)<= valeur;
        contrainte_capacite_rame(s,t)= rame(s,t)<=valeur;
    end
end
% Contraintes de transport

% Contraintes de transport
contrainte_transport_entree = sum(entrant,2)==sum(sum(trafic,3),2);
contrainte_transport_sortie =  sum(sortant,2)==sum(sum(trafic,3),1)';
%% 
% 

prob = optimproblem;
prob.Objective = sum(d);
prob.Constraints.cons1 = contrame1;
prob.Constraints.cons2 = contrame2;
prob.Constraints.cons3 = contrame3; 

prob.Constraints.cons4 = consfrequence;

prob.Constraints.cons5 = contrainte_entree_rame;
prob.Constraints.cons6 = contrainte_sortie_rame;
prob.Constraints.cons7 = contrainte_capacite_rame;

prob.Constraints.cons8 = contrainte_transport_entree;
prob.Constraints.cons9 = contrainte_transport_sortie;
   
prob.Constraints.cons10 = contraintes_flux_quais;
            
solve(prob)
%%