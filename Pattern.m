classdef Pattern
    %PATTERN class with properties
    %   poly_list   : array of Polygons building the DOT pattern
    %   poly_size   : maximum expected size (circumference) of a Polygon
    %   poly_area   : maximum expected arear of a Polygon
    %   color_prop_true : true color 1 proportion w.r.t. total area
    %   color_prop_expc : expected color 1 proportion w.r.t. toal area
    %   color_balance   : default true, false if expected proportion not reached
    
    properties
        poly_list
        poly_size = 40
        poly_area = 100
        color_prop_true
        color_prop_expc
        color_balance = true
    end   
end

