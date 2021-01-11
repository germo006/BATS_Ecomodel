function [lambda, I] = Sunlight(t,z)
%SUNLIGHT Calculates irradiance at a given time (of day, year) and a depth
%in an idealized water column.
%   OUTPUT
%   lambda  =   This is a vector or scalar of wavelength bins. As of 9 Oct,
%               it's just a scalar (length l)
%   I       =   Irradiance (mol photons/s at the corresponding wl). This is
%               a matrix (m x n x l)
%   
%   INPUT
%   t       =   Vector of times in datenum format (length n).
%   z       =   Vector, depth in meters. Used at some future time for 
%               attenuation (length m). 
%

lambda = NaN;
I = zeros(length(z), length(t), length(lambda));
[T, Z] = meshgrid(t, z); 
sunrise = 6; sunset = 20; 

time = (T - floor(T)); % Just get the sub-day data

Im = 10^3;      % Max irradiance (amplitude).

I = Im*(-cos(time*pi*2));
I(I<=0) = 0;
end

