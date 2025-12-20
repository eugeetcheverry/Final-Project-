function LLH = ecef2llh(ECEF)

global REQ_EARTH
global e_EARTH;

x = ECEF(1);
y = ECEF(2);
z = ECEF(3);
r_delta = sqrt(x^2 + y^2);
lon = atan2(y/r_delta, x/r_delta);
gdlat = atan(z/r_delta);
gdlat_old = gdlat-1;
while abs(gdlat-gdlat_old) > 1e-4
    gdlat_old = gdlat;
    C_earth = REQ_EARTH/sqrt(1-e_EARTH^2*sin(gdlat)^2);
    gdlat = atan((z+C_earth*e_EARTH^2*sin(gdlat))/r_delta);
end
if abs(gdlat) > pi/2-1/180*pi
    hellp = z/sin(gdlat) - C_earth*(1-e_EARTH^2);
else
    hellp = r_delta/cos(gdlat) - C_earth;
end
LLH = [gdlat; lon; hellp];