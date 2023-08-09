%% This section creates data.
%
% firingRatesAverage: N x S x D x T
%
% N is the number of neurons
% S is the number of actions (initiate the saccade or wait)
% D is the number of task type i.e. stop signal and choice signal (D=2)
% T is the number of time-points (note that all the trials should have the
% same length in time!)
%
% trialNum -- number of trials for each neuron in each S,D condition (is
% usually different for different conditions and different sessions)
%
%
% firingRatesAverage -- average of firingRates over trials (5th dimension).
% If the firingRates is filled up with nans, then it's simply
%    firingRatesAverage = nanmean(firingRates,5)
% If it's filled up with zeros (as is convenient if it's stored on hard 
% drive as a sparse matrix), then 
%    firingRatesAverage = bsxfun(@times, mean(firingRates,5), size(firingRates,5)./trialNum)

%vars readme
% 1) Choice: 1 when selected go target; 2 when selected stop target
% 2) Correct: Accuracy of Choice; data = 1 when correct, 0 if incorrect; Meaningful in a stop signal task, and is always 1 for economic choice task. Actual reward magnitudes can be found from last two data variables 8 and 9)
% 3) Task type: 1= stop signal task, 2= economic choice task
% 4) Presentation of stop signal / presentation of offer 2: 1= present; 0 = absent
% 5) SSD: Stop signal delay time in secs from the presentation of go offer
% 6) Left / right side presentation of go target (left  = 1 / right = 2)
% 7) reward associated with go target / offer 1
% 8) reward associated with the stop target / offer 2
% 9) Reaction time: reaction time for choosing a target; 

tstart_go = 5000;
mint_go = -250;
maxt_go = 420;
span = 150;

load('/Users/Peeyushi/Documents/MATLAB/T_psth1ms_go.mat')
for neuroni = 1:length(neuroncumT)
    
    if neuroni == 1
        minT = 1;
    else
        minT = neuroncumT(neuroni-1)+1;
    end
    maxT = neuroncumT(neuroni);
    rangeT = [minT:maxT];
    data=horzcat(psth,vars);
    % generating the nested data
    ind1 = data(data(minT:maxT,15003)==1,:);
    %1= stop signal task, 2= choice task
    ind2 = ind1(ind1(:,15004)==0,:);
    ind4 = ind2(ind2(:,15009) <= 0.5000, :);
    raw_data=ind4;
    final = smooth(nanmean(raw_data),span).';
    nesteddata(neuroni,1,1,:)=final(1,mint_go+tstart_go:maxt_go+tstart_go);
    
    
    ind1 = data(data(minT:maxT,15003)==2,:);
    ind2 = ind1(ind1(:,15004)==0,:);
    ind3 = ind2(ind2(:,15007)==0.1000,:); %1= stop cue present, 0= absent (based on stop stim presentation)
    indd = ind3(ind3(:,15009) <= 0.5000,:);
    raw_data=indd(:,mint_go+tstart_go:maxt_go+tstart_go);
    final = smooth(nanmean(raw_data),span).';
    nesteddata(neuroni,1,2,:)=final(1,:);
end
%%
load('/Users/Peeyushi/Documents/MATLAB/T_psth1ms_stop.mat')
data=horzcat(psth,vars);
mint_stop = -250;
maxt_stop = 420;
span = 150;
for neuroni = 1:length(neuroncumT)
    
    if neuroni == 1
        minT = 1;
    else
        minT = neuroncumT(neuroni-1)+1;
    end
    maxT = neuroncumT(neuroni);
    rangeT = [minT:maxT];

    % generating the nested data
    ind1 = data(data(minT:maxT,15003)==1,:); %1= stop signal task, 2= choice task
    ind2 = ind1(ind1(:,15004)==1,:);
    ind4 = ind2(ind2(:,15009)<=0.5000,:);
    trials1 = size(ind4,1);
    for j = 1:trials1
        ssd = ind4(j,15005);
        tstart_stop = 5000-ssd;
        raw_data(j,:) = ind4(j,mint_stop+tstart_stop:maxt_stop+tstart_stop);
        raw_data(any(isnan(raw_data), 2), :) = [];
        raw_data( all(~raw_data,2), : ) = [];
    end
    nesteddata(neuroni,2,1,:) = smooth(mean(raw_data),span).';
end
%%
load('/Users/Peeyushi/Documents/MATLAB/T_psth1ms_stop.mat')
data=horzcat(psth,vars);
mint_go = -250;
maxt_go = 420;
span = 150;
for neuroni = 1:length(neuroncumT)
    
    if neuroni == 1
        minT = 1;
    else
        minT = neuroncumT(neuroni-1)+1;
    end
    maxT = neuroncumT(neuroni);
    rangeT = [minT:maxT];

    % generating the nested data
    ind1 = data(data(minT:maxT,15003)==2,:); %1= stop signal task, 2= choice task
    ind2 = ind1(ind1(:,15004)==1,:);
    ind3 = ind2(ind2(:,15008)==0.2000,:);%1= stop cue present, 0= absent (based on stop stim presentation)
    ind4 = ind3(ind3(:,15009)<=0.5000,:);
    trials1 = size(ind4,1);
    for j = 1:trials1
        ssd = ind4(j,15005);
        tstart_stop = 5000-ssd;
        raw_data(j,:) = ind4(j,mint_stop+tstart_stop:maxt_stop+tstart_stop);
        raw_data(any(isnan(raw_data), 2), :) = [];
        raw_data( all(~raw_data,2), : ) = [];
    end
    nesteddata(neuroni,2,2,:) = smooth(mean(raw_data),span).';
end
% computing PSTHs
firingRatesAverage = nesteddata; %nanmean(firingRates, 4);

%% Define parameter grouping

% firingRates array has [N S D T E] size; herewe ignore the 1st dimension 
% (neurons), i.e. we have the following parameters:
%    1 - action 
%    2 - task
%    3 - time
% There are three pairwise interactions:
%    [1 3] - action/time interaction
%    [2 3] - task/time interaction
%    [1 2] - action/task interaction
% And one three-way interaction:
%    [1 2 3] - rest

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'go-accept/stop-reject', 'stop/choice task', 'Condition-ind', 'T/A Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
time = (1:length(firingRatesAverage)) / 1000;
timeEvents = [0.2500,0.4000,0.5200];

%% dPCA without regularization and ignoring noise covariance

% This is the core function.
% W is the decoder, V is the encoder (ordered by explained variance),
% whichMarg is an array that tells you which component comes from which
% marginalization

tic
[W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams);
toc

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

dpca_plot2(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);



