#include <cassert>
#include <stdexcept>
#include "Polynomial.hpp"

// Create a polynomial of degree 'N - 1' where 'N' is
// the number of objects in the initializer list.
//
// The 'N' coefficients are stored in an array of that
// capacity
//  - The capacity of the array must be at least one
Polynomial::Polynomial(
  std::initializer_list<std::complex<double>> init
):
degree_{ (init.size() == 0) ? 0 : init.size() - 1 },
coeffs_{ new std::complex<double>[degree_ + 1] } {
  // If the user creates a polynomial with no arguments,
  // we will create the zero polynomial.
  if ( init.size() == 0 ) {
    coeffs_[0] = 0.0;

  // Otherwise, the first argument is the constant
  // coefficient, the second is the coefficient of
  // 'x', and so on, until the last is the
  //                 N-1
  // coefficient of x   .
  } else {
    std::size_t k{ 0 };

    for ( auto itr{ init.begin() };
      itr != init.end();
      ++k, ++itr
    ) {
      coeffs_[k] = *itr;
      std::cout << coeffs_[k] << std::endl;
    }

    if ( (degree_ >= 1) && (coeffs_[degree_] == 0.0) ) {
      delete[] coeffs_;
      coeffs_ = nullptr;
      throw std::invalid_argument{ "The leading coefficient must be non-zero" };
    }
  }
}

Polynomial::~Polynomial() {
  delete[] coeffs_;
  coeffs_ = nullptr;
  degree_ = 0;
}

std::complex<double> Polynomial::operator()( std::complex<double> z ) const {
  return evaluate( coeffs_, degree_, z );
}

// This function is writen for you, but you can change it if
// you wish. It is only necessary that it retern a std::vector
// instance with the roots. It must be 'const'.
std::vector<std::complex<double>> Polynomial::roots() const {
  // The vector of roots to be returned by this function
  //   - The argument is the number of entries in the std::vector
  std::vector<std::complex<double>> roots_vec( degree_ );

  // An array of the coefficients
  //  - This will be manipulated by you as you find the roots
  std::complex<double> coeffs[degree_ + 1];

  // The polynomial may have a leading coefficient that is not
  // equal to '1'. It will be easier if the leading coefficient
  // is '1', so the coeffs array will store the coefficients
  // divided by the leadning coefficient, which makes the
  // leading coefficient '1'.

  assert( std::abs( coeffs_[degree_] ) > 1e-16 );
  coeffs[degree_] = 1.0;

  for ( unsigned int k{ 0 }; k < degree_; ++k ) {
    coeffs[k] = coeffs_[k]/coeffs_[degree_];
  }

  // In general, methods like Newton's method and the secant
  // may require a complex initial value. If you implement
  // Muller's method, you do not need this special case and
  // your algorithm can use either real or complex initial
  // values. You could, if you wanted, always use a complex
  // initial value; but this has the nasty side-effect that
  // if you have a polynomial with real coefficients, then the
  // roots, if they are real, will be approximated with values
  // like 1.414213562373095 + 0.3957188605e-53j

  for ( unsigned int k{ 0 }; k < degree_; ++k ) {
    roots_vec[k] = find_root( coeffs, degree_ - k );
    divide( coeffs, degree_ - k, roots_vec[k] );
  }

  return roots_vec;
}

// Evaluate the polynomial at the point 'x'

std::complex<double> Polynomial::evaluate(
  std::complex<double> const *coeffs,
  unsigned int         const  degree,
  std::complex<double> const  z
) {
  // Horner's rule implemented as in lecture 5.1.
  // This is not the most efficient method learned
  // Not worth my time for the more efficient method for a 1% faster run.

    std::complex<double> result {coeffs[degree]};
    for (unsigned int k{degree - 1}; k <= degree; --k){
      result = coeffs[k] + result * z;
    }
    return result;
}

// Given the polynomial of degree 'degree', the coefficients
// of which are given in the array 'coeffs', find an approximation
// of a root.
// If your algorithm requires you to start with a non-real
// complex root, you can generate use the following:
//  std::complex<double> x0;
//
//  if ( use_complex ) {
//    // Generate a pseudo-random complex number where each
//    // component is a randomly generated number between
//    // -10 and 10.
//    x0 = std::complex<double>( rand()*20e-10 - 10.0,
//                               rand()*20e-10 - 10.0 );
//  } else {
//    // Otherwise, generate a pseudo-random real number
//    // from the interval [-10, 10]
//    x0 = rand()*20e-10 - 10.0;
//  }

std::complex<double> Polynomial::find_root(
  std::complex<double>      *coeffs,
  unsigned int         const degree
) {

//////////////////// MULLERS METHOD: DOES NOT WORK ////////////////////
  /*
  // Uses Muller's Method to evaluate roots
  // constants have no roots
  assert(degree > 0);
  // Create initial x approximations
  std::complex<double> X[3];
  for (unsigned int k{0}; k < 3; k++){
   X[k] = rand()*20e-10 - 10.0;
  }
  // Order the initial x approximates
  if ( abs(evaluate(coeffs, degree, X[0])) < abs(evaluate(coeffs, degree, X[1])) )
    std::swap(X[0], X[1]);
  if ( abs(evaluate(coeffs, degree, X[0])) < abs(evaluate(coeffs, degree, X[2])) )
    std::swap(X[0], X[2]);
  if ( abs(evaluate(coeffs, degree, X[1])) < abs(evaluate(coeffs, degree, X[1])) )
    std::swap(X[1], X[2]);
  // Mullers method iteration
  // condition checks for zero to be in range of +/- degree * 10^-12
  while( (degree * 1e-12) < abs(divide(coeffs, degree, X[2])) ){
    // Assigns evaluated polynomial
    std::complex<double> f0 = evaluate(coeffs, degree, X[0]);
    std::complex<double> f1 = evaluate(coeffs, degree, X[1]);
    std::complex<double> f2 = evaluate(coeffs, degree, X[2]);
    // helper functions for interpolated quadratic poly
    std::complex<double> d0 = f0 - f2;
    std::complex<double> d1 = f1 - f2;
    std::complex<double> h0 = X[0] - X[2];
    std::complex<double> h1 = X[1] - X[2];
    // interpolated quadratic poly
    std::complex<double> c = f2;
    std::complex<double> b = ((d1*h0*h0) - (d0*h1*h1)) / ((h0*h1) * (h0 - h1));
    std::complex<double> a = ((d0*h1) - (d1*h0)) / ((h0*h1) * (h0 - h1));
    // quadratic root solutions
    std::complex<double> q_add = (-2.0*c) / (b + sqrt(b*b - 4.0*c*a));
    std::complex<double> q_sub = (-2.0*c) / (b - sqrt(b*b - 4.0*c*a));
    // setting new approximated values to the array
    X[0] = X[1];
    X[1] = X[2];
    // setting quadratic root to the less value abs
    X[2] = (abs(q_add) > abs(q_sub))? q_sub : q_add ;
  }
  return X[2];
  */

  //////////////////// NEWTONS METHOD ////////////////////
  /*
  // constants have no roots unless constant is 0 which has infinite roots
  assert(degree > 0);
  // Initial Approximation
  std::complex<double> x0 {std::complex<double>( rand()*20e-10 - 10.0, rand()*20e-10 - 10.0 )};
  // Finds derivative of poly and stores to new array
  std::complex<double> df_coeffs[degree - 1];
  for (unsigned int k{0}; k < degree; k++){
    df_coeffs[k] = (k + 1.0) * coeffs[k + 1];
  }
  // Evaluates Newtons method until in range of +/- degree * 10^-12
  while( (degree * 1e-12) < std::abs(evaluate(coeffs, degree, x0)) ){
    x0 -= ( evaluate(coeffs, degree, x0) / evaluate(df_coeffs, (degree - 1), x0) );
  }
  // return the found root
  return x0;
  */

  //////////////////// NEWTONS METHOD AS IN LECS, EDITED ////////////////////
  /*
  // constants have no roots unless constant is 0 which has infinite roots
  assert(degree > 0);
  std::complex<double> df_coeffs[degree - 1];
  for (unsigned int k{0}; k < degree; k++){
    df_coeffs[k] = (k + 1.0) * coeffs[k + 1];
  }
  std::size_t max_iterations {10000};
  std::complex<double> x0 {std::complex<double>( rand()*20e-10 - 10.0, rand()*20e-10 - 10.0 )};
  std::complex<double> f0 = {evaluate( coeffs, degree, x0 )};
  if ( f0 == 0.0 ) {return x0;}
  for ( unsigned int k{0}; k < max_iterations; ++k ) {
    std::complex<double> df = evaluate(df_coeffs, degree - 1, x0);
    if ( df == std::complex<double>(0,0) ) {return NAN;} // x0 is at an extrema
    std::complex<double> x1{ x0 - f0/df };
    // Possibly the result of a division-by-zero, or a division by a very small number.
    std::complex<double> f1{ evaluate(coeffs, degree, x1) };
    if ( (std::abs( f1 ) < (degree * 1e-12)) ) {return x1;}
    x0 = x1;
    f0 = f1;
  }
  return NAN;
  */

  //////////////////// SECANT METHOD ////////////////////

  // constants have no roots unless constant is 0 which has infinite roots
  assert(degree > 0);
  // Initial random complex approximations for the starting 2 points
  std::complex<double> x0 {std::complex<double>( rand()*20e-10 - 10.0, rand()*20e-10 - 10.0 )};
  std::complex<double> x1 {std::complex<double>( rand()*20e-10 - 10.0, rand()*20e-10 - 10.0 )};
  // Evaluates secant method until in range of +/- degree * 10^-12
  int count {0};
  while( (degree * 1e-12) < std::abs(evaluate(coeffs, degree, x1))){
    std::complex<double> temp_x {x1};
    x1 -= evaluate(coeffs, degree, x1) *
      ( (x1 - x0)/(evaluate(coeffs, degree, x1) - evaluate(coeffs, degree, x0)) );
    x0 = temp_x;
    count++;
    // Recalls find_root if the function starts to diverge or takes too long to converge
    if (abs(evaluate(coeffs, degree, x1)) > abs(evaluate(coeffs, degree, x1)) || count > 1e20){
    	x1 = find_root(coeffs, degree);
    	break;
    }
  }
  // return the found root when loop breaks
  return x1;
}

// Polynomial division
//
// In place, update the entries in the coefficient array to
// give the quotient when the polynomial defined by
//                                                    degree
//    coeffs[0] + coeffs[1] x + ... + coeffs[degree] x
// is divided by x - r

std::complex<double> Polynomial::divide(
  std::complex<double>      *coeffs,
  unsigned int         const degree,
  std::complex<double> const r
) {
  // Synthetic division is used to divide the polynomial.
  // we keep the highest degree coefficient the same so we start at the next coeff
  for (unsigned int k{degree - 1}; k < degree; k-- ){
    coeffs[k] += r * coeffs[k + 1];
  }
  // shifts array over by one to put coeffs in correct index
  std::complex<double> remainder {coeffs[0]};
  for (unsigned int k{0}; k < degree; k++ ){
    coeffs[k] = coeffs[k + 1];
  }
  //sets the remainder to the end of the array
  coeffs[degree] = remainder;
  return remainder;
}

std::string to_string_term( unsigned int n ) {
  if ( n == 0 ) {
    return "";
  } else if ( n == 1 ) {
    return "z";
  } else {
    return "z^" + std::to_string( n );
  }
}

std::string to_string_complex( std::complex<double> z, bool is_leading, unsigned int n ) {
  if ( z.imag() == 0.0 ) {
    // The coefficient is real (possibly zero)
    if ( z.real() == 0.0 ) {
      if ( is_leading && (n == 0) ) {
        return "0";
      } else {
        return "";
      }
    } else {
      if ( is_leading ) {
        return std::to_string( z.real() ) + to_string_term( n );
      } else {
        if ( z.real() < 0.0 ) {
          return " - " + std::to_string( -z.real() ) + to_string_term( n );
        } else {
          return " + " + std::to_string(  z.real() ) + to_string_term( n );
        }
      }
    }
  } else if ( z.real() == 0.0 ) {
    // The coefficient is a non-zero imaginary number
    if ( is_leading ) {
      return std::to_string( z.imag() ) + "j" + to_string_term( n );
    } else {
      if ( z.imag() < 0.0 ) {
        return " - " + std::to_string( -z.imag() ) + "j" + to_string_term( n );
      } else {
        return " + " + std::to_string(  z.imag() ) + "j" + to_string_term( n );
      }
    }
  } else {
    // Both real and imaginary parts are non-zero
    if ( is_leading && (n == 0) ) {
      if ( z.imag() < 0.0 ) {
        return std::to_string( z.real() ) + " - " + std::to_string( -z.imag() ) + "j";
      } else {
        return std::to_string( z.real() ) + " + " + std::to_string(  z.imag() ) + "j";
      }
    } else {
      bool is_negative{ false };
      std::string out{ "" };

      if ( z.real() < 0.0 ) {
        z = -z;
        is_negative = true;
      }

      if ( z.imag() < 0.0 ) {
        out = out + std::to_string( z.real() ) + " - " + std::to_string( -z.imag() ) + "j";
      } else {
        out = out + std::to_string( z.real() ) + " + " + std::to_string(  z.imag() ) + "j";
      }

      if ( is_negative ) {
        if ( is_leading ) {
          out = "-(" + out + ")";
        } else {
          out = " - (" + out + ")";
        }
      } else if ( n > 0 ) {
        out = "(" + out + ")";

        if ( !is_leading ) {
          out = " + " + out;
        }
      } else {
        out = " + " + out;
      }

      return out + to_string_term( n );
    }
  }
}

std::string Polynomial::to_string() const {
  std::string out{ to_string_complex( coeffs_[degree_], true, degree_ ) };

  for ( unsigned int n{ degree_ - 1 }; n < degree_; --n ) {
    out += to_string_complex( coeffs_[n], false, n );
  }

  return out;
}

std::ostream &operator<<( std::ostream &out, Polynomial const &p ) {
  return out << p.to_string();
}