#include "dirVector.h"

DirVector::DirVector() { i = 0; j = 0;  k = 0; }

DirVector::DirVector( double x, double y, double z ) { i = x; j = y;  k = z; }

DirVector DirVector::operator-( const DirVector& dir )
{
  return DirVector( i - dir.i, j - dir.j, k - dir.k );
}

DirVector DirVector::operator*( const DirVector& dir )
{
  DirVector temp;
  temp.i = ( j * dir.k ) - ( dir.j * k );
  temp.j = ( k * dir.i ) - ( dir.k * i );
  temp.k = ( i * dir.j ) - ( dir.i * j );
  return temp;
}

DirVector& DirVector::operator=( const DirVector& dir )
{
  i = dir.i; j = dir.j; k = dir.k;
  return *this;
}

double DirVector::magnitude()
{
  return sqrt( pow( i, 2 ) + pow( j, 2 ) + pow( k, 2) );
}

std::string DirVector::point()
{
  return ( std::to_string( i ) + " "
             + std::to_string( j ) + " "
             + std::to_string( k ) );
}

// Dot product of two vectors
double dot( DirVector& a, DirVector& b )
{
  return ( ( a.i * b.i ) + ( a.j + b.j ) + ( a.k + b.k ) );
}

double angle_between( DirVector& a, DirVector& b )
{
  double angle;
  angle = std::acos( dot( a, b ) / ( a.magnitude() * b.magnitude() ) );
  return angle;
}
