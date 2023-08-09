function dpca_plot_default(data, time, yspan, explVar, compNum, events, signif, marg)

% Modify this function to adjust how components are plotted.
%
% Parameters are as follows:
%   data      - data matrix, size(data,1)=1 because it's only one component
%   time      - time axis
%   yspan     - y-axis spab
%   explVar   - variance of this component
%   compNum   - component number
%   events    - time events to be marked on the time axis
%   signif    - marks time-point where component is significant
%   marg      - marginalization number


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% displaying legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(data, 'legend')
    
    % if there is only time and no other parameter - do nothing
    if length(time) == 2
        return

    % if there is one parameter
    elseif length(time) == 3
        numOfStimuli = time(2); % time is used to pass size(data) for legend
        colors = lines(numOfStimuli);
        hold on
        
        for f = 1:numOfStimuli
            plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
            text(1.2, f, ['Stimulus ' num2str(f)])
        end
        axis([0 3 -1 1.5+numOfStimuli])
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        set(gca,'Visible','off')
        return

    % two parameters: stimulus and decision (decision can only have two
    % values)
    elseif length(time) == 4 && time(3) == 2
        numOfStimuli = time(2); % time is used to pass size(data) for legend
        colors = lines(numOfStimuli);
        hold on
        
        for f = 1:numOfStimuli
            plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
            text(1.2, 1, ['Go or Accept Trials'])
            text(1.2 , 2, ['Stop or Reject Trials'])
        end
        plot([0.5 1], [-2 -2], 'k', 'LineWidth', 2)
        plot([0.5 1], [-3 -3], 'k--', 'LineWidth', 2)
        text(1.2, -2, 'Stop signal task')
        text(1.2, -3, 'Choice signal task')
        
        axis([0 3 -4.5 1.5+numOfStimuli])
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        set(gca,'Visible','on')
        return
        
    % other cases - do nothing
    else
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up the subplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Yspan = [-260,260];
Xspan = [0,0.6600];
if isempty(time)
    time = 1:size(data, ndims(data));
end
axis([time(1) time(end) Yspan])
hold on

if ~isempty(explVar)
    title(['Component #' num2str(compNum) ' [' num2str(explVar,'%.1f') '%]'])
else
    title(['Component #' num2str(compNum)])
end

if ~isempty(events)
    plot([events; events], Yspan, 'Color', [0.6 0.6 0.6])
end

if ~isempty(signif)
    signif(signif==0) = nan;
    plot(smooth(time, signif + Yspan(1) + (Yspan(2)-Yspan(1))*0.05, 'k', 'LineWidth', 3),100)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting the component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ndims(data) == 2
    % only time - plot it
    just=squeeze(data(1,:));
    plot(time, squeeze(data(1, :)), 'k', 'LineWidth', 2)

elseif ndims(data) == 3
    % different stimuli in different colours
    numOfStimuli = size(data, 2);
    plot(time, squeeze(data(1,:,:)), 'LineWidth', 2)

elseif ndims(data) == 4 && size(data,3)==2
    % different stimuli in different colours and binary condition as
    % solid/dashed
    numOfStimuli = size(data, 2);
    colors = lines(numOfStimuli);
    colors(1,:)=[0.9290 0.6940 0.1250];
    colors(2,:)=[0.6350 0.0780 0.1840];

    transpose1=squeeze(data(1, 1, 1, :)).';
    plot(time, transpose1, 'color', colors(1,:), 'LineWidth', 2)
    transpose2=squeeze(data(1, 2, 1, :)).';
    plot(time, transpose2, 'color', colors(2,:), 'LineWidth', 2)
    transpose3=squeeze(data(1, 1, 2, :));
    plot(time, transpose3, '--', 'color', colors(1,:), 'LineWidth', 2)
    transpose4=squeeze(data(1, 2, 2, :));
    plot(time, transpose4, '--', 'color', colors(2,:), 'LineWidth', 2)

    %for f=1:numOfStimuli 
     %   transpose1=squeeze(data(1, 1, 1, :)).';
      %  transpose2=squeeze(data(1, 2, 1, :)).';
       % plot(time, squeeze(data(1, f, 1, :)), 'color', colors(f,:), 'LineWidth', 2)
        %plot(time, squeeze(data(1, f, 2, :)), '--', 'color', colors(f,:), 'LineWidth', 2)
    %end

else
    % in all other cases pool all conditions and plot them in different
    % colours
    data = squeeze(data);
    dims = size(data);
    data = permute(data, [numel(dims) 1:numel(dims)-1]);
    data = reshape(data, size(data,1), []);
    data = data';
    
    plot(smooth(time, data, 'LineWidth', 2),100)   
end

xlim([0.01 0.77])
ylim([-250 250])
