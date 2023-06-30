function color_map = create_color_map3(color1,color2,color3, n_colors)
if nargin<4
    n_colors = 72; 
end

boundaries = [1,round(n_colors/3)+1,round(2*n_colors/3)+1];

color_map(boundaries(1),:) = color1;
color_map(boundaries(2),:) = color2;
color_map(boundaries(3),:) = color3;

for c = 1:3
    color_map(boundaries(1):boundaries(2),c) = linspace(color_map(boundaries(1),c),color_map(boundaries(2),c), boundaries(2)-boundaries(1)+1); 
    color_map(boundaries(2):boundaries(3),c) = linspace(color_map(boundaries(2),c), color_map(boundaries(3),c), boundaries(3)-boundaries(2)+1);
    tmp = linspace(color_map(boundaries(3),c), color_map(boundaries(1),c), n_colors-boundaries(3)+2);
    color_map(boundaries(3):n_colors,c) = tmp(1:(end-1));
end   
    
    
    