#include<Rcpp.h>
#include<cmath>

using namespace Rcpp;


// convert decimal degrees to radians
double deg2rad(double deg) {
    return (deg * M_PI / 180);
}


// Calculate the area of a slice of the globe from the equator to the parallel
// at latitude f (on WGS84 ellipsoid). Based on:
// https://gis.stackexchange.com/questions/127165/more-accurate-way-to-calculate-area-of-rasters
float slice_area(float f) {
    float a = 6378137; // in meters
    float b = 6356752.3142; // in meters,
    float e = sqrt(1 - pow(b / a, 2));
    float zp = 1 + e * sin(f);
    float zm = 1 - e * sin(f);
    return(M_PI * pow(b, 2) * ((2 * atanh(e * sin(f))) / (2 * e) + sin(f) / (zp * zm)));
}


// Formula to calculate area of a raster cell on WGS84 ellipsoid, following
// https://gis.stackexchange.com/questions/127165/more-accurate-way-to-calculate-area-of-rasters
float calc_cell_area(float ymin, float ymax, float x_width) {
    if (ymin > ymax) {
        float temp;;
        temp = ymax;
        ymax = ymin;
        ymin = temp;
    }
    // ymin: minimum latitude
    // ymax: maximum latitude
    // x_width: width of cell in degrees
    return((slice_area(deg2rad(ymax)) - slice_area(deg2rad(ymin))) * (x_width / 360.));
}


// [[Rcpp::export]]
DataFrame area_ha(DataFrame r, float xres, float yres, bool use_cov_frac) {
    NumericVector cell = r["cell"];
    NumericVector y = r["y"];
    NumericVector cov_frac = r["coverage_fraction"];

    int n_out_rows = cell.length();

    DataFrame out = DataFrame::create(
        Named("cell") = NumericVector(n_out_rows),
        Named("area_ha") = NumericVector(n_out_rows)
    );
                                      

    NumericVector cell_out = out["cell"];
    NumericVector area_out = out["area_ha"];

    // loop over cells in this polygon
    for (int row=0; row < r.nrow(); row++) {
        cell_out[row] = cell[row];
        area_out[row] = calc_cell_area(y[row] - yres/2, y[row] + yres/2, xres) / 10000;
        if (use_cov_frac) {
            area_out[row] = area_out[row] * cov_frac[row];
        }
    }
    return(out);
}
