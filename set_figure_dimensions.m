function set_figure_dimensions(width, height, padding, fontsize, units)
  % set_figure_dimensions(width, height)
  %
  % Sets the dimensions of the axis borders to width_in x height_in. This
  % allows for more precise control over the output of the figure when
  % putting the figure into Inkscape to edit.
  %
  % This is done by setting the width and height of the 'Position' property
  % of the axis to 'width_in' and 'height_in', respectively. The actual
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
  % Padding: float
  %     Amount of padding around the axis Position. Defines the area of the
  %     OuterPosition property. Default value is 1cm.
  % fontsize : float
  %     Fontsize of the figure labels and text etc. Default value is 9pt.
  % units : string
  %     Units of the figure. Same options as the 'Units' property of
  %     figures and axes. Default value is 'centimeters'.
  
  % Default value for units
  arguments
    width    = 6.0;
    height   = 4.0;
    padding  = 1.0;
    fontsize = 9;
    units    = 'centimeters';
  end

  % Get Matlab figure and axis properties
  fig = gcf();
  ax  = gca();

  % Set units
  fig.Units = units;
  ax.Units  = units;

  % Set fontsize
  ax.FontSize = fontsize;

  % Set axis dimensions and position
  ax.Position = [padding, padding, width, height];

  % Set figure dimensions and position, with extra padding
  fig.Position = [5, 5, width + (2 * padding), height + (2 * padding)];
  
  % Set paper size for printing
  fig.PaperPosition = [0, 0, fig.Position(3), fig.Position(4)];
  fig.PaperSize     = fig.Position(3:4);


end