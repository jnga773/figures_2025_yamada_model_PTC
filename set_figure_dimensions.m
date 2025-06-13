function set_figure_dimensions(width, height, options)
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
  %     Width of the axis plot. Default value is 6.0cm.
  % height : float
  %     Height of the axis plot. Default value is 4.0cm.
  % Padding : float
  %     Amount of padding around the axis Position. Defines the area of the
  %     OuterPosition property. Default value is 1.0cm.
  % Scale : float
  %     Scaling factor for figure dimensions. Default value is 1.0.
  % fontsize : float
  %     Fontsize of the figure labels and text etc. Default value is 9.0pt.
  % units : string
  %     Units of the figure. Same options as the 'Units' property of
  %     figures and axes. Default value is 'centimeters'.
  % figposition : float
  %     Sets the lower left anchor point of the figure window on the
  %     screen. Default value is 5.0.
  
  % Default value for units
  arguments
    width    = 6.0;
    height   = 4.0;

    % Optional arguments
    options.padding     = 1.0;
    options.scale       = 1.0;
    options.fontsize    = 9.0;
    options.units       = 'centimeters';
    options.figposition = 5.0;
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
  ax.Position = [options.padding, ...
                 options.padding, ...
                 width, ...
                 height];

  % Set figure dimensions and position, with extra padding
  fig.Position = [options.figposition, ...
                  options.figposition, ...
                  width + (2 * options.padding), ...
                  height + (2 * options.padding)];

  % Scale Positions
  % fig.Position = options.scale * fig.Position;
  % ax.Position  = options.scale * ax.Position;
  
  % Set paper size for printing
  % fig.PaperPosition = [0, 0, fig.Position(3), fig.Position(4)];
  % fig.PaperSize     = fig.Position(3:4);

end