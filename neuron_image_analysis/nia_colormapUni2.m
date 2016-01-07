function cmap = nia_colormapUni2()
%NIA_COLORMAPUNI2 Create a unidirectional colormap
%   cmap = nia_colormapUni2() returns a colormap which is well-suited
%   for displaying data that is unidirectional (is not centered
%   around zero).  It runs from dark blue colors to red colors, 
%   and saturates at white.  The color scheme has been selected so
%   that the entire colorspace falls within the CYMK gamut, and thus
%   reproduces very well in both screen and print forms.

Chsv = [...
    0/6, 240/360, 0.65, 0.45; ...
    1/6, 215/360, 0.65, 0.70; ...
    2/6, 175/360, 0.40, 0.75; ...
    3/6, 70/360, 0.7, 0.75; ...
    4/6, 55/360, 0.92, 1.0; ...
    5/6, 10/360, 0.85, 0.95; ...
    6/6, -20/360, 0.0, 0.9];

total_pts = 100;
frac = linspace(0, 1, total_pts)';
Chsv_interp = interp1(Chsv(:,1), Chsv(:,2:4), frac);

Chsv_interp_mask = Chsv_interp(:,1) < 0;
Chsv_interp(Chsv_interp_mask, 1) = Chsv_interp(Chsv_interp_mask, 1) + 1;

cmap = hsv2rgb(Chsv_interp);

end