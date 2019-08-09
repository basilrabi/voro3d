#include <math.h>
#include <stdio.h>
#include <Rcpp.h>
#include <voro++.hh>

#include "dirVector.h"

// If ever x  is too small, use a threshold value for the dimensions of the
// container.
double setThreshold ( double x )
{
  // Hard coded threshold is 2 meters
  const double threshold = 2;

  if ( x < threshold )
    return threshold;
  else
    return x;
}

//' Create Voronoi Diagram
//'
//' Create cell-based voronoi diagram using three-dimensional points. The
//'   polyhedral surface of each cell is defined in well-known text format.
//'
//' @param x numeric vector of the x-coordinates of the points
//' @param y numeric vector of the y-coordinates of the points
//' @param z numeric vector of the z-coordinates of the points
//' @param containerRatio numeric ratio between the length of the container to
//'   be created and the length of the bounding box of the points
//' @return character vector defining the voronoi cells (polyhedral surface)
//'   in well-known text
//' @export
// [[Rcpp::export]]
Rcpp::StringVector voronoi( Rcpp::NumericVector x,
                            Rcpp::NumericVector y,
                            Rcpp::NumericVector z,
                            double containerRatio )
{
  DirVector point;
  double cells, i, j, k;
  double xLength, yLength, zLength;
  double xMin, xMax, yMin, yMax, zMin, zMax;
  double conMarginX, conMarginY, conMarginZ;
  double conXMin, conXMax, conYMin, conYMax, conZMin, conZMax;
  int ii, jj, kk, ll, mm, nn, nx, ny, nz;
  R_xlen_t n, cGCount;
  std::string polygon, polyhedralsurface;
  std::vector< DirVector > points;
  std::vector< double > vertices;
  std::vector< int > faceVertices;
  std::vector< int > vertexOrder;
  voro::voronoicell vc;
  voro::particle_order po;
  voro::wall_list wl;

  n = x.length();

  if ( n != y.length() || n != z.length() )
    Rcpp::stop( "Lengths of coordinate vectors are not equal." );

  if ( n < 2 )
    Rcpp::stop( "Cannot generate cells if points are less than 2." );

  if ( containerRatio < 1 )
    Rcpp::stop( "Invalid containerRatio: Value must not be less than 1." );

  Rcpp::StringVector cellGeometry ( n );
  cGCount = 0;

  // Bounding box vertices
  xMin = Rcpp::min( x );
  yMin = Rcpp::min( y );
  zMin = Rcpp::min( z );
  xMax = Rcpp::max( x );
  yMax = Rcpp::max( y );
  zMax = Rcpp::max( z );

  // Bounding box dimensions
  xLength = setThreshold( xMax - xMin );
  yLength = setThreshold( yMax - yMin );
  zLength = setThreshold( zMax - zMin );

  // Margin of container based on the ratio (multiplying factor)
  conMarginX = xLength * ( containerRatio - 1 ) / 2;
  conMarginY = yLength * ( containerRatio - 1 ) / 2;
  conMarginZ = zLength * ( containerRatio - 1 ) / 2;

  // Container vertices
  conXMin = xMin - conMarginX;
  conYMin = yMin - conMarginY;
  conZMin = zMin - conMarginZ;
  conXMax = xMax + conMarginX;
  conYMax = yMax + conMarginY;
  conZMax = zMax + conMarginZ;

  // Number of divisions per axis
  cells = cbrt( n / ( 5.6 * xLength * yLength * zLength ) );
  nx = int( xLength * cells + 1 );
  ny = int( yLength * cells + 1 );
  nz = int( zLength * cells + 1 );

  // Initialize container
  voro::container con( conXMin, conXMax, conYMin, conYMax, conZMin, conZMax,
                       nx, ny, nz, false, false, false, 8 );
  con.add_wall( wl );

  // Add points to container
  for ( R_xlen_t i = 0; i < n; i++ )
    con.put( po, i, x[i], y[i], z[i] );

  // Compute voronoi cells
  voro::c_loop_order clo( con, po );
  if ( clo.start() )
  {
    do
    {
      if ( con.compute_cell( vc, clo ) )
      {

        // Store the coordinates of the particle
        clo.pos( i, j, k );

        // Store coordinates of each vertex. Each set of vertex coordinates is
        // stored at every 3 elements in `vertices`
        vc.vertices( i, j, k, vertices );

        long unsigned int vCounter = 0;
        for ( ; vCounter < vertices.size(); vCounter += 3 )
        {
          point = DirVector( vertices[vCounter],
                             vertices[vCounter + 1],
                             vertices[vCounter + 2] );
          points.push_back( point );
        }

        polyhedralsurface = "POLYHEDRALSURFACE(";

        // Append unique face data to polyhedral surface. Each face has 3
        // vertices and the data for each face indicates the indices of the 3
        // vertices. The control flow below is copied from voro++ since the API
        // is hard to understand.
        for ( ii = 1; ii < vc.p; ii++ )
        {
          for ( jj = 0; jj < vc.nu[ii]; jj++ )
          {
            kk = vc.ed[ii][jj];
            if ( kk >= 0 )
            {
              vc.ed[ii][jj] = -1 - kk;
              ll = vc.cycle_up( vc.ed[ii][vc.nu[ii] + jj], kk );
              mm = vc.ed[kk][ll];
              vc.ed[kk][ll] = -1 - mm;
              while ( mm != ii )
              {
                nn = vc.cycle_up( vc.ed[kk][vc.nu[kk] + ll], mm );
                polygon = "((" +
                  points[ii].point() + ", " +
                  points[kk].point() + ", " +
                  points[mm].point() + ", " +
                  points[ii].point() + "))";

                if ( polyhedralsurface == "POLYHEDRALSURFACE(" )
                  polyhedralsurface += polygon;
                else
                  polyhedralsurface += ", " + polygon;

                kk = mm;
                ll = nn;
                mm = vc.ed[kk][ll];
                vc.ed[kk][ll] = -1 - mm;
              }
            }
          }
        }

        polyhedralsurface += ")";
        cellGeometry[cGCount] = polyhedralsurface;
        points.clear();
      }

      else
        cellGeometry[cGCount] = NA_STRING;

      cGCount++;

    } while ( clo.inc() );
  }

  return cellGeometry;
}
