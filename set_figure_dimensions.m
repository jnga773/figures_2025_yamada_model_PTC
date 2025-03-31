function set_figure_dimensions(width, height, options)
  % set_figure_dimensions(width, height)
  %
  % Sets the dimensions of the axis borders to width x height. This
  % allows for more precise control over the output of the figure when
  % putting the figure into Inkscape to edit.
  %
  % This is done by setting the width and height of the 'Position' property
  % of the axis to 'width' and 'height', respectively. The actual
  % figure size (fig.Position) is a little bit large to allow for axis
  % labelling and tick labelling.
  %
  % The default units are in centimeters. 
  %
  % Input
  % ----------
  % width : float
  %     Width of the axis plot. Default value is 6cm.
  % height : float
  %     Height of the axis plot. Default value is 4cm.
  % Padding : float
  %     Amount of padding around the axis Position. Defines the area of the
  %     OuterPosition property. Default value is 1cm.
  % Scale : float
  %     Scaling factor for figure dimensions. Default value is 1.
  % fontsize : float
  %     Fontsize of the figure labels and text etc. Default value is 9pt.
  % units : string
  %     Units of the figure. Same options as the 'Units' property of
  %     figures and axes. Default value is 'centimeters'.
  
  % Default value for units
  arguments
    width double   = 6.0;
    height double  = 4.0;

    % Optional arguments
    options.padding double  = 1.0;
    options.scale double    = 1.0;
    options.fontsize double = 9;
    options.units char    = 'centimeters';
  end

  % Get Matlab figure and axis properties
  fig = gcf();
  ax  = gca();

  % Set units
  fig.Units = options.units;
  ax.Units  = options.units;

  % Set fontsize
  ax.FontSize = options.fontsize;

  % Set axis dimensions and position
  ax.Position = [options.padding, options.padding, width, height];

  % Set figure dimensions and position, with extra padding
  fig.Position = [5, 5, width + (2 * options.padding), height + (2 * options.padding)];

  % Scale Positions
  fig.Position = options.scale * fig.Position;
  ax.Position  = options.scale * ax.Position;
  
  % Set paper size for printing
  fig.PaperPosition = [0, 0, fig.Position(3), fig.Position(4)];
  fig.PaperSize     = fig.Position(3:4);

end
