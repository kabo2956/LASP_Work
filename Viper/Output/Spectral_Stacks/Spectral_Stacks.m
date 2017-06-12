function [out1] = Spectral_Stacks(filename,singlecup)

% Read in a filename and produce the spectral plots
% Hard coding in the dimensions of the arrays at the moment



x = load(filename);
sizex = size(x);
ny = sizex(1);
nx = sizex(2);
%fill_val = find(x==12180000);
%x(fill_val) = 1;
array = zeros(ny,nx);
array(:) = x;
nnn = 0;

if singlecup == 1;
    xrang = [2:129];
elseif singlecup == 0;
    xrang = [2:nx];
end

%hold on
    for jjj = 1:ny;
       if nx == 513;
           Date    = array(jjj,1);
           cur     = array(jjj,2:nx);
           if singlecup == 1;
               cur(129) = cur(128); % Remove jump at end
           end
       elseif nx == 516;
           Date    = array(jjj,1);
           Hour    = array(jjj,2);
           Minute  = array(jjj,3);
           Second  = array(jjj,4);
           Date    = Date + Hour/24 + Minute/(24*60);
           cur     = array(jjj,xrang);
       end

       zero_val = find(cur<=10);
       
       cur(zero_val) = 10;
       vec = max(xrang)-min(xrang)+1;
       x1 = linspace(1,vec,vec); %Equivalent of findgen
       j1 = x1; %Get appropriate size
       j1(:) = Date; %Fill with appropriate value

       
       plot3(x1,j1,cur(xrang),'Color','black');
       axis([1 vec min(array(:,1)) max(array(:,1)) 1 10000000]);
       set(gca, 'XScale', 'linear', 'YScale', 'linear', 'ZScale', 'log');
       nnn = nnn+1;
       if nnn == 1, hold on, end
       %if nnn == ny, hold off
       jh = jjj;%
    end
hold off

% Label axes
xlabel('Channel Number', 'FontSize', 18)
set(gca, 'XTick', [1 32 64 96 128])
ylabel('DOY (1979)', 'FontSize', 18)
zlabel('Current (fA)', 'FontSize', 18)

%[az,el] = view Use this command to get the view that is desired and then
%set it with the command below.
if filename == 'Spectral_Stacks_Model_064_06_01_34_064_15_58_28.txt';
    view(58.5,70)
    set(gcf, 'Position', [364 482 1189 864]);
elseif filename == 'Spectral_Stacks_Model_064_10_07_59_064_15_58_28.txt';
    view(35.5,62)
    set(gcf, 'Position', [364 482 1189 864]);
end

end