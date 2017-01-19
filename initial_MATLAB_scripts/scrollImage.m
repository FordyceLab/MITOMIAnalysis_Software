function [api]=scrollImage(im)

% Create a figure without toolbar and menubar.
hfig = figure('Toolbar','none',...
              'Menubar', 'none',...
              'Name','Graphical Feature Interface',...
              'NumberTitle','off',...
              'IntegerHandle','off');

% Display the image in a figure with imshow.
himage = imshow(im);

% Add the scroll panel.
hpanel = imscrollpanel(hfig,himage);

% Position the scroll panel to accommodate the other tools.
set(hpanel,'Units','normalized','Position',[0 .1 1 .9]);

% Add the Magnification box.
hMagBox = immagbox(hfig,himage);

% Position the Magnification box
pos = get(hMagBox,'Position');
set(hMagBox,'Position',[0 0 pos(3) pos(4)]);

% Add the Overview tool.
hovervw = imoverview(himage);
pos2=get(hovervw,'Position');
set(hovervw,'Position',[0 0 pos2(3)+100 pos2(4)+100]);

%Grab API controls
api = iptgetapi(hpanel);
