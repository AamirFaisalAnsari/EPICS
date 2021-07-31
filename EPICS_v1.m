% EPICS=Effective Pairwise Interactions for predicting Community Structure

% EPICS efficiently computes accurate estimates of effective pairwise interactions 
%between species in a microbial community. For unlocking full features of EPICS, see
%here: 

%The following funcion takes n(number of species in the community),
%species abundances in monoculture and leave-one-out subcommunities and
%produces effective pairwise interactions and community structure using EPICS as
%outputs.
clear all;
close all;
clc;
disp('Welcome to EPICS...');

%number of species in the community
number_of_species=input('\nInput number of species in the community=');

%monoculture abundances
ind_abun=input('\nInput monoculture abundances=');
%[0.527000000000000,0.946300000000000,0.705500000000000,0.640100000000000,0.676900000000000,0.849400000000000,0.0899000000000000,0.227800000000000];

%leave-one-out abundances
loo_abun=input('\nInput abundances in leave-one-out subcommunities \n(Rows are subcommunities, columns are species)=');
%[0,0.915277277997681,0.915277277997681,0.915277277997681,0.915277277997681,0.915277277997681,0.915277277997681,0.915277277997681;0.0364102305315166,0,0.0561209565134392,0.0427122755072731,0.0427101879974147,0.0781864267333450,0.0238261628350117,0.0411171070287187;0.470237183352849,0.0622958836158715,0,0.0427122755072731,0.0427101879974147,0.0546752019550849,0.0422339697936855,0.0729313610939620;0.262290406021041,0.0145499683142133,0.0410953357273667,0,0.0181820740132557,0.850792923649617,0.0184078069586738,1.86027023239039e-29;0.0327092179354460,0.0528765566831936,0,0.0446325920280813,0,2.94540337000384e-30,0.0425079844553841,0.0566122942625211;0.261339784414631,0.957577681900169,0.928562758087671,0.947559558298925,0.960303035723643,0,0.954755825395353,0.896179729739669;0.0889459195077008,0.0427323747209261,0.0410953357273667,0.0427122755072731,0.0427101879974147,0.0482751756641382,0,0.0737372002382907;0.483528156105948,0.0913083252363468,0.0726118530437985,0.100006659023204,0.113615556274697,0.0698804209591128,0.114669184814066,0];

%Abundances in the original community-This is required for the comparison
%against EPICS predicted abundances
community_abun=[0.915277277997681,0.0279919384152069,0.0427414353114217,0.0149380406571958,0.0427414353114217,0.955852322733725,0.0427414353114217,0.0947870305092372];



%This code can be used for both GLV and Replicator dynamics model

global n %number of species in the community
global growth_rate_switch %To switch-off/on intrinsic growth rate terms
global total_abundance_constant_switch %To switch-off/on the average term for substracting 

%For GLV set 'growth_rate_switch=1' and 'total_abundance_constant_switch=0'
%For Replicator set 'growth_rate_switch=0' and 'total_abundance_constant_switch=1'

growth_rate_switch=1;total_abundance_constant_switch=0; %GLV

n=number_of_species;
for i=1:n
    A(i,i)=-1/ind_abun(i);
    diag(i)=A(i,i);
end

x_3=zeros(1,n*(n-1));
options = optimset('MaxFunEvals',1000000*n,'MaxIter',1000000*n,'Algorithm','levenberg-marquardt','Display','off');

%This function iteratively computes effective pairwise interactions using
%EPICS
[eff_intn,fval3,exitflag3,output3]  = fsolve(@(x0) obtain_effective_EPICS(x0,loo_abun,diag),x_3,options);

predicted_eff_intn=-1*ones(n,n);
club_effective=[];
for i=1:n
    for j=1:n
        if i~=j
            if j>i
                predicted_eff_intn(i,j)=eff_intn(n*(i-1)+j-i);
            else
                predicted_eff_intn(i,j)=eff_intn(n*(i-1)+j-i+1);
            end
            club_effective=[club_effective,predicted_eff_intn(i,j)];
            
        end
    end
end

for i=1:n
    predicted_eff_intn(i,i)=A(i,i);
end


% to predict abundance in the original community using EPICS-generated effective intns
x_4=0.5*ones(1,n)/n;

options = optimset('MaxFunEvals',1000000*n,'MaxIter',1000000*n,'Algorithm','levenberg-marquardt','Display','off');
[EPICS_orig_abun,fval4,exitflag4,output4]  = fsolve(@(x0) obtain_n_member_EPICS(x0,predicted_eff_intn),x_4,options);

fmt=['\nPredicted Abundances in the original community using EPICS =' repmat(' %1.0f',1,numel(EPICS_orig_abun)) '\n'];
fprintf(fmt,EPICS_orig_abun);

predicted_eff_intn

%To generate the bar plot for the comparison
%bar([community_GR',bb_n_member_new_method']);

%% To obtain effective pairwise interactions using EPICS

function to_return=obtain_effective_EPICS(inputk,N,diag)
global n
global growth_rate_switch
global total_abundance_constant_switch

input=-1*ones(n,n);
for i=1:n
    for j=1:n
        if i~=j
            if j>i
                input(i,j)=inputk(n*(i-1)+j-i);
            else
                input(i,j)=inputk(n*(i-1)+j-i+1);
            end
        end
    end
end

a=input;

for i=1:n
    a(i,i)=diag(i);
end

for tzl=1:n
    for ann=1:n
        if tzl~=ann
            y(tzl,ann)=1*growth_rate_switch+dot(a(tzl,:),N(:,ann));
        end
    end
end

for tzl=1:n
    for ann=1:n
        if tzl~=ann
            y_new(tzl,ann)=y(tzl,ann)-total_abundance_constant_switch*dot(y(:,ann),N(:,ann));%/sum(N(:,ann));
        end
    end
end

to_return=sum(sum(y_new.^2));
end

%% to obtain the abundances in the original community using EPICS

function to_return=obtain_n_member_EPICS(inputk,eff_intn)

global n
global growth_rate_switch
global total_abundance_constant_switch

abund_n_mem=inputk;

for tzl=1:n
    y(tzl)=1*growth_rate_switch+dot(eff_intn(tzl,:),abund_n_mem);
end


for tzl=1:n
    y_new(tzl)=y(tzl)-total_abundance_constant_switch*dot(y,abund_n_mem);%/sum(abund_n_mem);
end

to_return=sum((y_new.^2));
end
