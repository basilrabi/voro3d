#ifndef DIRVECTOR_H
#define DIRVECTOR_H

#include <math.h>
#include <string>

class DirVector
{
public:

  double i, j, k;

  // Constructor
  DirVector();
  DirVector( double, double, double );

  // Vector subtraction operator
  DirVector operator- (const DirVector&);

  // Cross-product operator
  DirVector operator* (const DirVector&);

  // Desctructor
  ~DirVector(){}

  // Copy Constructor
  DirVector& operator= ( const DirVector& dir );

  // Magnitude of vector
  double magnitude();

  // Output vector in space delimited string
  std::string point();

};

// Dot product of two vectors
double dot( DirVector& a, DirVector& b );

// Angle between two vectors
double angle_between( DirVector& a, DirVector& b );

#endif
