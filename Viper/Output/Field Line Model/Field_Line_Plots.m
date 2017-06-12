function [out1] = Field_Line_Plots(filename, Paperno)

  x = load(filename);
  sizex = size(x);
  nel   = sizex(1);
  
  if Paperno == 2
    rad = x(:,1);
    z   = x(:,2);
    n1  = x(:,3);
    n2  = x(:,4);
    n3  = x(:,5);
    n4  = x(:,6);
    n5  = x(:,7);
    n6  = x(:,8);
    n7  = x(:,9);
    n8  = x(:,10);
    n9  = x(:,11);
  elseif Paperno == 3
    rad = x(:,1);
    z   = x(:,2);
    n1  = x(:,3);
    n2  = x(:,4);
    n3  = x(:,5);
    n4  = x(:,6);
    n5  = x(:,7);
    n6  = x(:,8);
    n7  = x(:,9);
    n8  = x(:,10);
    n9  = x(:,11);
  elseif Paperno == 2.1
    rad = x(:,1);
    z   = x(:,2);
    n1  = x(:,3);
    n2  = x(:,4);
    n3  = x(:,5);
    n4  = x(:,6);
    n5  = x(:,7);
    n6  = x(:,8);
    n7  = x(:,9);
    n8  = x(:,10);
    n9  = x(:,11);
    n10 = x(:,12);
  elseif Paperno == 3.1
    rad = x(:,1);
    z   = x(:,2);
    n1  = x(:,3);
    n2  = x(:,4);
    n3  = x(:,5);
    n4  = x(:,6);
    n5  = x(:,7);
    n6  = x(:,8);
    n7  = x(:,9);
    n8  = x(:,10);
    n9  = x(:,11);
    n10 = x(:,12);
  end
    
  
  r = [5.9,9.4,15,26.3]; % Moon distances for Io, Europa, Ganymede, and Callisto
  r0 = r(:)*0;
  r0 = transpose(r0);
  
  %Quick reduction to test meshgrid
  %rad = rad(1:20000);
  %z   = z(1:20000);
  %n1  = n1(1:20000);
  
  % Initialize figure;
  %figure;
  if Paperno == 1;
        axestitles = char('n_{e}','n_{e}','n_{O^{++}}','n_{O^{++}}/n_{e}','n_{O^+}','n_{O^+}/n_{e}','n_{s^{++}}','n_{s^{++}}/n_{e}','n_{s^+}','n_{s^+}/n_{e}','n_{s^{+++}}','n_{s^{+++}}/n_{e}');
        axestitles = char('n_{e}','n_{e}','n_{O^{+}}','n_{O^{++}}','n_{S^+}','n_{s^{++}}','n_{s^{+++}}','n_{H^{+}}','n_{Na^{+}}','n_{O^{+}_Hot}');
  end
  %axestitles = char('n_{H^+}','n_{O^{++}}','n_{O^+}','n_{s^{++}}','n_{s^+}','n_{s^{+++}}','n_{e}');
  
  if Paperno == 2;
    % 1 is ne
    % 2 is O+
    % 3 is O++
    % 4 is S+
    % 5 is S++
    % 6 is S+++
    % 7 is O+ Hot
    subarr       = zeros(nel,7);
    subarr(:,1)  = n9(:); % e
    subarr(:,2)  = n1(:)./n9(:); % O+/ne
    subarr(:,3)  = n2(:)./n9(:); % O++/ne
    subarr(:,4)  = n3(:)./n9(:); % S+/ne
    subarr(:,5)  = n4(:)./n9(:); % S++/ne
    subarr(:,6)  = n5(:)./n9(:); % S+++/ne
    subarr(:,7)  = n8(:)./n9(:); % O+ Hot/ne
    subarr(:,2)  = n1(:); % O+
    subarr(:,3)  = n2(:); % O++
    subarr(:,4)  = n3(:); % S+
    subarr(:,5)  = n4(:); % S++
    subarr(:,6)  = n5(:); % S+++
    subarr(:,7)  = n8(:); % O+ Hot
    
    axestitles = char('n_{e}','n_{O^{+}}','n_{O^{++}}','n_{S^+}','n_{s^{++}}','n_{s^{+++}}','n_{O^{+}_{Hot}}');

    
    % Establish a vector of positions
    posit       = zeros(7,4);
    posit(1,:)  = [0.30 0.83 0.36 0.12];
    posit(2,:)  = [0.30 0.70 0.36 0.12];
    posit(3,:)  = [0.30 0.57 0.36 0.12];
    posit(4,:)  = [0.30 0.44 0.36 0.12];
    posit(5,:)  = [0.30 0.31 0.36 0.12];
    posit(6,:)  = [0.30 0.18 0.36 0.12];
    posit(7,:)  = [0.30 0.05 0.36 0.12];



    % Generate a figure
    figure1 = figure(1);
    for j = 1:7;
      hghghghghghg = j %Checks what element
    
    
      % Grid data
      [xq,yq] = meshgrid(0:0.01:40, 0:0.01:10); %0.01 works well for spacing between points  
      vq = griddata(rad,z,subarr(:,j),xq,yq);

      clims = [-2,3.5]; % Have upper and lower limits of 1d4 and 1d-3 respectively for the
      %clims = [-3.5,1]; % Use this one if relative (n/ne)
      %if j == 1;
        axes1 = axes('Parent',figure1,...
        'Position',posit(j,:));
        box(axes1,'on');
        hold(axes1,'on');
      %else
       
      imagesc(xq(1,:),yq(:,1),log10(vq),clims);
      set(gca, 'ydir', 'normal')
      %title(axestitles(j,:));
      xfill = [0,5,5,0,0];
      yfill = [0,0,10,10,0];
      fill(xfill,yfill,[0.5 0.5 0.5],'Edgecolor','none')
      if j ==7;
        text(1,7.5,axestitles(j,:),'Color','White','FontSize',16)
      else
        text(1,8,axestitles(j,:),'Color','White','FontSize',16)
      end
      % Overplot Magnetic Field Lines at 5, 10, 20, 30, and 40 Rj
      % Saved in multiple text files since the arrays are not entirely the
      % same size due to number of points dependent on arc length of field
      % line
      bfield6  = load('CANSheet_6_Rj.txt');
      bfield10 = load('CANSheet_10_Rj.txt');
      bfield20 = load('CANSheet_20_Rj.txt');
      bfield30 = load('CANSheet_30_Rj.txt');
      bfield40 = load('CANSheet_40_Rj.txt');
      
      %xfill = cat(1,bfield6(:,1),0,bfield6(1,1));
      %yfill = cat(1,bfield6(:,2),0,bfield6(1,2));
      %fill(xfill,yfill,[0.2081    0.1663    0.5292],'Edgecolor','none')
      
      plot(bfield6(:,1),bfield6(:,2),'.','Color','White');
      plot(bfield10(:,1),bfield10(:,2),'.','Color','White');
      plot(bfield20(:,1),bfield20(:,2),'.','Color','White');
      plot(bfield30(:,1),bfield30(:,2),'.','Color','White');
      plot(bfield40(:,1),bfield40(:,2),'.','Color','White');

      axis equal
    
      % Remove y axis on the even subplot numbers  
      if j == 7;
          set(gca, 'YTick', [0 5 10])
          ylabel('Z (R_J)');
          set(gca, 'XTick', [5 6 10 20 30 40])
          xlabel('{\rho} (R_J)');
      else
          set(gca, 'YTick', [0 5 10])
          ylabel('Z (R_J)');
          set(gca, 'XTick', [5 6 10 20 30 40])
          set(gca, 'XTickLabel', ['' '' '' '' '' '' ''])
      end
    
      pos = get(gca,'Position');
      h = colorbar;
      ylabel(h,'Log10[Density (/cc)]');
      set(gca, 'FontSize', 14);
      set(h, 'YAxisLocation', 'right')
      if j <= 6;
          % Only do a colorbar at the end
          h.Visible = 'off';
          set(gca, 'Position', pos);
      elseif j == 7;
          h.Location = 'East';
          h.Position = [0.75 0.2 0.037 0.6];
      end
    
      set(gca, 'Layer', 'Top'); %Axis black lines on top
   
      
    end
    % Set final plot position
    plotpos = [72 208 1338 1032];
    set(gcf, 'Position', plotpos);
    set(gca, 'Layer', 'Top'); %Axis black lines on top
  end
  
  
  
  
  
  
  
  
  
  
  if Paperno == 2.1;
    % 1 is ne
    % 2 is O+
    % 3 is O++
    % 4 is S+
    % 5 is S++
    % 6 is S+++
    % 7 is O+ Hot
    subarr       = zeros(nel,8);
    subarr(:,1)  = n10(:); % e
    subarr(:,2)  = n1(:)./n10(:); % O+/ne
    subarr(:,3)  = n2(:)./n10(:); % O++/ne
    subarr(:,4)  = n5(:)./n10(:); % S+/ne
    subarr(:,5)  = n4(:)./n10(:); % S++/ne
    subarr(:,6)  = n3(:)./n10(:); % S+++/ne
    subarr(:,7)  = n9(:)./n10(:); % SO2+/ne
    subarr(:,8)  = n8(:)./n10(:); % O+ Hot/ne
    subarr(:,2)  = n1(:); % O+
    subarr(:,3)  = n2(:); % O++
    subarr(:,4)  = n5(:); % S+
    subarr(:,5)  = n4(:); % S++
    subarr(:,6)  = n3(:); % S+++
    subarr(:,7)  = n9(:); % SO2+
    subarr(:,8)  = n8(:); % O+ Hot
    
    axestitles = char('n_{e}','n_{O^{+}}','n_{O^{++}}','n_{S^+}','n_{s^{++}}','n_{s^{+++}}','n_{SO_2^{+}}','n_{O^{+}_{Hot}}');

    
    % Establish a vector of positions
    posit       = zeros(8,4);
    posit(1,:)  = [0.20 0.70 0.24 0.20];
    posit(2,:)  = [0.65 0.70 0.24 0.20];
    posit(3,:)  = [0.20 0.50 0.24 0.20];
    posit(4,:)  = [0.65 0.50 0.24 0.20];
    posit(5,:)  = [0.20 0.30 0.24 0.20];
    posit(6,:)  = [0.65 0.30 0.24 0.20];
    posit(7,:)  = [0.20 0.10 0.24 0.20];
    posit(8,:)  = [0.65 0.10 0.24 0.20];



    % Generate a figure
    figure1 = figure(1);
    for j = 1:8;
      hghghghghghg = j %Checks what element
    
    
      % Grid data
      [xq,yq] = meshgrid(4.8:0.002:6, 0:0.002:2); %0.01 works well for spacing between points  
      vq = griddata(rad,z,subarr(:,j),xq,yq);

      clims = [0,4.5]; % Have upper and lower limits of 1d4 and 1d-3 respectively for the
      clims = [0,4]; % Use this one if relative (n/ne)
      %if j == 1;
        axes1 = axes('Parent',figure1,...
        'Position',posit(j,:));
        box(axes1,'on');
        hold(axes1,'on');
      %else
       
      imagesc(xq(1,:),yq(:,1),log10(vq),clims);
      set(gca, 'ydir', 'normal')
      %title(axestitles(j,:));
      if j ==7;
        text(1,7.5,axestitles(j,:),'Color','White','FontSize',16)
      else
        text(1,8,axestitles(j,:),'Color','White','FontSize',16)
      end
      % Overplot Magnetic Field Lines at 5, 10, 20, 30, and 40 Rj
      % Saved in multiple text files since the arrays are not entirely the
      % same size due to number of points dependent on arc length of field
      % line
      %bfield6  = load('CANSheet_6_Rj.txt');
      %bfield10 = load('CANSheet_10_Rj.txt');
      %bfield20 = load('CANSheet_20_Rj.txt');
      %bfield30 = load('CANSheet_30_Rj.txt');
      %bfield40 = load('CANSheet_40_Rj.txt');
      
      %xfill = cat(1,bfield6(:,1),0,bfield6(1,1));
      %yfill = cat(1,bfield6(:,2),0,bfield6(1,2));
      %fill(xfill,yfill,[0.2081    0.1663    0.5292],'Edgecolor','none')
      
      %plot(bfield6(:,1),bfield6(:,2),'.','Color','White');
      %plot(bfield10(:,1),bfield10(:,2),'.','Color','White');
      %plot(bfield20(:,1),bfield20(:,2),'.','Color','White');
      %plot(bfield30(:,1),bfield30(:,2),'.','Color','White');
      %plot(bfield40(:,1),bfield40(:,2),'.','Color','White');

      %axis equal
    
      % Remove y axis on the even subplot numbers  
      if j == 8;
          set(gca, 'YTick', [0 2])
          ylabel('Z (R_J)');
          set(gca, 'XTick', [4.8 5 6])
          xlabel('{\rho} (R_J)');
      else
          set(gca, 'YTick', [0 2])
          ylabel('Z (R_J)');
          set(gca, 'XTick', [4.8 5 6])
          set(gca, 'XTickLabel', ['' '' ''])
      end
    
      pos = get(gca,'Position');
      h = colorbar;
      ylabel(h,'Log10[Density (/cc)]');
      set(gca, 'FontSize', 14);
      set(h, 'YAxisLocation', 'right')
      if j <= 7;
          % Only do a colorbar at the end
          h.Visible = 'off';
          set(gca, 'Position', pos);
      elseif j == 8;
          h.Location = 'East';
          h.Position = [0.5 0.2 0.037 0.6];
      end
    
      set(gca, 'Layer', 'Top'); %Axis black lines on top
   
      
    end
    % Set final plot position
    plotpos = [550 135 1995 912];
    set(gcf, 'Position', plotpos);
    set(gca, 'Layer', 'Top'); %Axis black lines on top
  end
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  if Paperno == 3;
    % 1 is ne
    % 2 is H+
    % 3 is H+/e
    subarr       = zeros(nel,7);
    subarr(:,1)  = n9(:); % e
    subarr(:,2)  = n6(:); % H+
    subarr(:,3)  = n2(:)./n9(:); % H+/e
    
    axestitles = char('n_{e}','n_{H^{+}}','n_{H^{+}} / n_{e}');

    
    % Establish a vector of positions
    posit       = zeros(3,4);
    posit(1,:)  = [0.10 0.70 0.75 0.25];
    posit(2,:)  = [0.10 0.40 0.75 0.25];
    posit(3,:)  = [0.10 0.10 0.75 0.25];




    % Generate a figure
    figure1 = figure(1);
    
    % Create ranges for color bars
    climsr      = zeros(3,2);
    climsr(1,:) = [-1,3];
    climsr(2,:) = [-1.5,1.5];
    climsr(3,:) = [-2,-1];
    
    % Loop over the different species
    for j = 1:3;
      hghghghghghg = j %Checks what element
    
    
      % Grid data
      [xq,yq] = meshgrid(0:0.01:40, 0:0.01:10); %0.01 works well for spacing between points  
      vq = griddata(rad,z,subarr(:,j),xq,yq);

      clims = climsr(j,:); % Have upper and lower limits of 1d4 and 1d-3 respectively for the
      %if j == 1;
        axes1 = axes('Parent',figure1,...
        'Position',posit(j,:));
        box(axes1,'on');
        hold(axes1,'on');
      %else
      
      imagesc(xq(1,:),yq(:,1),log10(vq),clims);
      set(gca, 'ydir', 'normal')
      %title(axestitles(j,:));
      xfill = [0,5,5,0,0];
      yfill = [0,0,10,10,0];
      fill(xfill,yfill,[0.5 0.5 0.5],'Edgecolor','none')
      text(1,8,axestitles(j,:),'Color','White','FontSize',16)
      % Overplot Magnetic Field Lines at 5, 10, 20, 30, and 40 Rj
      % Saved in multiple text files since the arrays are not entirely the
      % same size due to number of points dependent on arc length of field
      % line
      bfield6  = load('CANSheet_6_Rj.txt');
      bfield10 = load('CANSheet_10_Rj.txt');
      bfield20 = load('CANSheet_20_Rj.txt');
      bfield30 = load('CANSheet_30_Rj.txt');
      bfield40 = load('CANSheet_40_Rj.txt');
      
      %xfill = cat(1,bfield6(:,1),0,bfield6(1,1));
      %yfill = cat(1,bfield6(:,2),0,bfield6(1,2));

      
      plot(bfield6(:,1),bfield6(:,2),'.','Color','White');
      plot(bfield10(:,1),bfield10(:,2),'.','Color','White');
      plot(bfield20(:,1),bfield20(:,2),'.','Color','White');
      plot(bfield30(:,1),bfield30(:,2),'.','Color','White');
      plot(bfield40(:,1),bfield40(:,2),'.','Color','White');
      
      axis equal
    
      % Remove y axis on the even subplot numbers  
      if j == 3;
        set(gca, 'YTick', [0 5 10])
        ylabel('Z (R_J)');
        set(gca, 'XTick', [5 6 10 20 30 40])
        xlabel('{\rho} (R_J)');
      else
        set(gca, 'YTick', [0 5 10])
        ylabel('Z (R_J)');
        set(gca, 'XTick', [5 6 10 20 30 40])
        set(gca, 'XTickLabel', ['' '' '' '' '' ''])
      end
    
      pos = get(gca,'Position');
      h = colorbar;
      if j == 2;
          ylabel(h,'Log10[Density (/cc)]');
          
      end
      
      set(gca, 'FontSize', 14);
      set(h, 'YAxisLocation', 'right')
      
      h.Visible = 'on';
      h.Location = 'East';
      pos_h      = posit(j,:);
      pos_h      = [0.88 pos_h(2)+0.015 0.037 0.22];
      h.Position = pos_h;
      %if j <= 2;
      %  % Only do a colorbar at the end
      %  h.Visible = 'on';
      %  set(gca, 'Position', pos);
      %elseif j == 3;
      %  h.Location = 'East';
      %  h.Position = [0.88 0.2 0.037 0.6];
      %end
      
      set(gca, 'Layer', 'Top'); %Axis black lines on top
      
    end
    % Set final plot position
    plotpos = [932 350 945 689];
    set(gcf, 'Position', plotpos);
    set(gca, 'Layer', 'Top'); %Axis black lines on top
  end
  

end