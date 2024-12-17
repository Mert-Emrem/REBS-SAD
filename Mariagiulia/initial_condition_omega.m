function [om_x0, om_y0, om_z0, om_0] = initial_condition_omega(x,y,z)
om_x0 = x;
om_y0 = y;
om_z0 = z;

om_0 = [om_x0; om_y0; om_z0];