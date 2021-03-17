function [malehist,femalehist,diffhist,R]=huia(alpha,beta,sigmar,sigmag,N,M,mutrate,tmax,delta,p,epsilonR,epsilonGamma,freqSharedLoci)
% alpha: intensity of competition between males
% beta: intensity of competition between females
% sigmar: standard deviation of the resource distribution
% sigmag: standard deviation of the effciency function of resource use
% N: total number of individuals in the population
% M: number of same-sex individuals in each patch
% mutrate: mutation rate
% tmax: the maximum number of generations in simulations
% delta: the distance of mean size between the two different (larger or smaller) resources
% p: the probability of the two different resources coexisting in a patch
% epsilonR: shape parameter of the resource distribution function, when equals 0, the distribtuion function is the standard normal distribution. Details see the manuscript.
% epsilonGamma: shape parameter of the efficiency function of resource use. 
% freqSharedLoci: proportion of loci that is shared between males and females


%%% Choose the number of groups so that the number of individuals (both males and females) is the same in each group.
nofgroups=N/(2*M); 
if round(nofgroups)~=nofgroups disp('N and M - Choose again'); pause; end
if abs(round(p*nofgroups)-p*nofgroups)>0.0001 [p*nofgroups round(p*nofgroups)], disp('N, M and p - Choose again'); pause; end
if abs(round(((1-p)/2)*nofgroups)-((1-p)/2)*nofgroups)>0.0001 [((1-p)/2)*nofgroups round(((1-p)/2)*nofgroups)], disp('N, M and p - Choose again'); pause; end

% create all territories: first the proportion (1-p)/2 that have only the resource with a smaller mean size, then the proportion (1-p)/2 that have only the resource with a larger mean size, then the rest (p) that have both resources
x=linspace(0,1,1001);
replenishment=zeros(nofgroups,length(x));
replenishmentA=generalNormPDF(x,0.5-delta/2,sigmar,epsilonR);
replenishmentB=generalNormPDF(x,0.5+delta/2,sigmar,epsilonR);
for i=1:round(nofgroups*(1-p)/2)
    replenishment(i,:)=2*replenishmentA;
end
for i=round(nofgroups*(1-p)/2)+1:round(nofgroups*(1-p))
    replenishment(i,:)=2*replenishmentB;
end
for i=round(nofgroups*(1-p))+1:nofgroups
    replenishment(i,:)=replenishmentA+replenishmentB;
end

% indexing
loci=100; % total number of loci
malealleles = 1:(1-freqSharedLoci)*loci;
sharedalleles=1+(1-freqSharedLoci)*loci:(1+freqSharedLoci)*loci; 
femalealleles=(1+freqSharedLoci)*loci+1:2*loci; sex=2*loci+1; group=2*loci+2; phenotype=2*loci+3; counterfactual=2*loci+4; cond=2*loci+5;
Pop=unidrnd(2,N,2*loci)-1; % zeroes and ones for the alleles
Pop(:,sex)=[-ones(N/2,1); ones(N/2,1)]; % minus = males, plus = females
Pop(:,group:cond)=NaN; % these are not yet computed

grouping=ones(M,1)*(1:nofgroups); grouping=grouping(:);

malehist=NaN(loci+1,tmax);
femalehist=NaN(loci+1,tmax);

E=NaN(N,1001); % initialize the efficiency matrix
R=NaN(nofgroups,1001);

for t=1:tmax
    
    % create individual phenotypes
    Pop(:,phenotype)=mean(Pop(:,[malealleles,sharedalleles])')'.*(Pop(:,sex)<0)+mean(Pop(:,[femalealleles,sharedalleles])')'.*(Pop(:,sex)>0);
    Pop(:,counterfactual)=mean(Pop(:,[malealleles,sharedalleles])')'.*(Pop(:,sex)>0)+mean(Pop(:,[femalealleles,sharedalleles])')'.*(Pop(:,sex)<0);
            
    % create a table where row is the resource x value and column is the individual; the content is the efficiency of each individual for that resource
	E=generalNormPDF(ones(N,1)*x,Pop(:,phenotype)*ones(size(x)),sigmag,epsilonGamma);
    
	% divide males randomly into groups of size M
    Pop(find(Pop(:,sex)<0),group)=vectperm(grouping);
    % same for females
    Pop(find(Pop(:,sex)>0),group)=vectperm(grouping);
    
    % the equilibrium resource density is replenishment rate/total pressure, but this has to be computed separately for each group
    for i=1:nofgroups
        f=find(Pop(:,group)==i);
        R(i,:)=min(100,replenishment(i,:)./sum(E(f,:)));
        Pop(f,cond)=mean(((ones(length(f),1)*R(i,:)).*E(f,:))')';
    end
    % now the global offspring pool is produced; goal is to make N babies
    % the female propensity to be a mother is proportional to condition^beta
    % offspring structure for now is simply [mother_id sire_id] (where id is the row of the parents)
    females=find(Pop(:,sex)>0);
    offspring(:,1)=females(ddists(Pop(females,cond).^beta,N));
    % each offspring needs a sire - this has to be done within groups
    for i=1:nofgroups
        % who needs a dad in this group?
        f=find(Pop(offspring(:,1),group)==i);
        % who are the potential dads?
        dads=find(Pop(:,group)==i & Pop(:,sex)<0);
        if ~isempty(f)
            offspring(f,2)=dads(ddists(Pop(dads,cond).^alpha,length(f)));
        end
    end
    
    % create the new generation so that the first N/2 are males and the latter ones are females
    % (in our study we would like to keep the sex ratio really 1:1)
    Newgen(:,sex)=[-ones(N/2,1); ones(N/2,1)];
    Newgen(:,group:cond)=NaN;
    for i=1:N
        mom=offspring(i,1);
        dad=offspring(i,2);
        Newgen(i,1:2:(2*loci-1))=Pop(dad,(1:2:(2*loci-1))+unidrnd(2,1,loci)-1); % alleles from dad
        Newgen(i,2:2:(2*loci))=Pop(mom,(1:2:(2*loci-1))+unidrnd(2,1,loci)-1); % alleles from mom
    end
    
    % collect some data of the old population before they die
    malephenotypes=Pop(Pop(:,sex)<0,phenotype);
    malecounterfactuals=Pop(Pop(:,sex)<0,counterfactual);
    femalephenotypes=Pop(Pop(:,sex)>0,phenotype);
	femalecounterfactuals=Pop(Pop(:,sex)>0,counterfactual);
    malehist(:,t)=hist(malephenotypes,linspace(0,1,loci+1))';
    femalehist(:,t)=hist(femalephenotypes,linspace(0,1,loci+1))';
    diffhist(:,t)=hist([malephenotypes-malecounterfactuals; femalecounterfactuals-femalephenotypes],linspace(-1,1,2*loci+1))';
    
    % mutate the offspring somewhat
    Newgen(:,1:2*loci)=xor(Newgen(:,1:2*loci),rand(N,2*loci)<mutrate);
    
    % replace the parental generation with the new one
    Pop=Newgen;
    
    %figure(1); subplot(3,1,1); plot(x,E');
    %subplot(3,1,2); plot(x,replenishment,'k');
    %subplot(3,1,3); plot(x,R');
    %figure(2); hist(Pop(:,cond));
    
end