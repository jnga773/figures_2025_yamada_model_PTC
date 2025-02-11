function colourmap_out = scale_colour_map(scale_factor)
  % colourmap_out = scale_colour_map(scale_factor)
  % 
  % Rescale the colour map.

  n_colours_in = 2048;
  n_colours_out = 1024;

  % Get colour map
  colour_map = parula(n_colours_in);

  % Create a linear ramp the size of the colormap we actually want
  t = linspace(0,1,n_colours_out);
  % Apply whatever transform you like to the ramp
  t2 = t .^ scale_factor;

  % Use that to scale the big linear colormap into the small stretched one.
  colourmap_out = colour_map(1+floor((n_colours_in-1)*t2'),:);
end
