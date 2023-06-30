function draw_target(center_angle, angle_width, inner_radius, outer_radius, fill_flag, fill_color)


data = [inner_radius.*cos(linspace(center_angle+angle_width/2,center_angle-angle_width/2,20)) , ...
        cos( center_angle-angle_width/2 )*linspace(inner_radius,outer_radius,20), ...
        outer_radius.*cos(linspace(center_angle-angle_width/2,center_angle+angle_width/2,20)), ...
        cos( center_angle+angle_width/2 )*linspace(outer_radius ,inner_radius,20)];
data(2,:) = [inner_radius.*sin(linspace(center_angle+angle_width/2,center_angle-angle_width/2,20)), ...
             sin( center_angle-angle_width/2 )*linspace(inner_radius,outer_radius,20), ...
             outer_radius.*sin(linspace(center_angle-angle_width/2,center_angle+angle_width/2,20)), ... 
             sin( center_angle+angle_width/2 )*linspace(outer_radius ,inner_radius,20)];
         
 plot(data(1,:), data(2,:), 'Color', 'k')
if nargin>4 && fill_flag
    patch(data(1,:), data(2,:), fill_color)
end




