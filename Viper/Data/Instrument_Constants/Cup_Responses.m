function [out1] = Cup_Responses(Spacecraft)


% Read in a filename and produce the spectral plots
% Hard coding in the dimensions of the arrays at the moment

if Spacecraft == 1;
    n = 2925;
    array = zeros(n,14);
    filename = 'V1_Cup_Responses.txt';
    x = load(filename);
    array(:) = x;
end
if Spacecraft == 2;
    n = 15653;
    array = zeros(n,14);
    filename = 'V2_Cup_Responses.txt';
    x = load(filename);
    array(:) = x;
    %out_of_range = find(array(:,1)<=186);
       
    %cur(zero_val) = 10;
end

Dates = array(:,1);
Radii = array(:,2);
CupA  = array(:,5)./sqrt(array(:,3).^2 + array(:,4).^2 + array(:,5).^2);
CupB  = array(:,8)./sqrt(array(:,6).^2 + array(:,7).^2 + array(:,8).^2);
CupC  = array(:,11)./sqrt(array(:,9).^2 + array(:,10).^2 + array(:,11).^2);
CupD  = array(:,14)./sqrt(array(:,12).^2 + array(:,13).^2 + array(:,14).^2);

% Remove bad points
% Brute force method
if Spacecraft == 1
    indices = [462:508,2268:2276,2318:2328,2591:2606];
    CupA(indices) = -9;
    CupB(indices) = -9;
    CupC(indices) = -9;
    CupD(indices) = -9;
else
    indices = [3605:3627,3633:3635,3983:3985];
    CupA(indices) = -9;
    CupB(indices) = -9;
    CupC(indices) = -9;
    CupD(indices) = -9;
end


if Spacecraft == 1
    Times_ind = find(Dates(:)<=64.85 & Dates>=63);
    Times = Dates(Times_ind);
else
    Times_ind = find(Dates(:)<=191.75 & Dates>=188.5);
    Times = Dates(Times_ind);
end

CupA = CupA(Times_ind);
CupB = CupB(Times_ind);
CupC = CupC(Times_ind);
CupD = CupD(Times_ind);
Radii = Radii(Times_ind);

scatter(Times,CupA,8,'red');
hold on
    scatter(Times,CupB,8,'b');
    scatter(Times,CupC,8,'black');
    scatter(Times,CupD,8,[0 .5 0]);
hold off
ylabel('$R = V_z/|V|$');
xlabel('DOY 1979');

leg = legend('Cup A', 'Cup B', 'Cup C', 'Cup D', 'Location', 'west');
leg_h = findobj(leg);
set(leg_h, 'FontSize', 28);

% Set axis limits
axis([min(Times) max(Times) 0 1]);

ax1 = gca; % current axes
ax1.XColor = 'b';
ax1.YColor = 'b';
ax1.XColor = 'black';
ax1.YColor = 'black';
ylabel('R = V_z/|V|', 'FontSize', 32);
xlabel('DOY 1979', 'FontSize', 32);
if Spacecraft == 1
    title_n = 'Voyager 1 Faraday Cup Response to a Cold Corotating Beam';
    %title(title_n)
else
    title_n = 'Voyager 2 Faraday Cup Response to a Cold Corotating Beam';
    %title(title_n)
end

%title_l = title(title_n, 'FontSize', 20);
%loc_arr = [63.75,189];
%set(title_l, 'Position', [loc_arr(Spacecraft) 1.06]); 
ax1.FontSize = 30;

% Overplot axes
hold on
plot([min(Times) max(Times)],[1 1],'k')
plot(max(Times)*[1 1],[1 0],'k')

% Plot all of the Radial locations
if Spacecraft == 1;
    time_value   = [63.3473 63.6312 63.8857...
                    64.1360 64.4604 64.5428]; % Fill in appropriately
    radial_value = [25 20 15 10 5 5]; % Fill in approriately
else
    time_value   = [188.8540 189.2230 189.5840...
                    189.9440 190.3170 190.8970...
                    190.9810 191.5610]; % Fill in appropriately
    radial_value = [35 30 25 20 15 10 10 15]; % Fill in approriately
end

for kk = 1:numel(time_value);
    text(time_value(kk),1,sprintf('%3.0f',radial_value(kk))...
        ,'VerticalAlignment','Bottom','HorizontalAlignment','Center'...
        ,'FontSize', 26)
    plot(time_value(kk)*[1 1], [1 0.98], 'k')
end


set(gcf, 'Position', [98 182 2547 1260])

% Create text labels for where the S/C rolls and there are data gaps
if Spacecraft == 1;
    %text(62.74,0.6,'S/C Roll','FontSize',15)
    %text(64.75,-2,'S/C Roll','FontSize',15)
    text(64.7,0.35,'S/C Roll','FontSize',18)
    %text(64.92,-3.2,'Data Gap','FontSize',15)
elseif Spacecraft == 2;
    %text(187.7,0.4,'Data Gap','FontSize',15)
    text(189.4,0.92,'S/C Roll','FontSize',18)
end




end