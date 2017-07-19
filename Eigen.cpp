//
//



#include "Eigen.h"
template<> template<>
Eigen::MatrixBase<Eigen::Block<Matrix3frm, 1, 3, true> >::cross_product_return_type<Eigen::Block<Matrix3frm, 1, 3, true> >::type Eigen::MatrixBase<Eigen::Block<Matrix3frm, 1, 3, true> >::cross<Eigen::Block<Matrix3frm, 1, 3, true> >(Eigen::MatrixBase<Eigen::Block<Matrix3frm, 1, 3, true> > const& vector) const
{
  Vector3f v1(this->x(), this->y(), this->z());
  Vector3f v2(vector.x(), vector.y(), vector.z());

  return Vector3f(v1.y() * v2.z() - v1.z() * v2.y(),
                  v1.z() * v2.x() - v1.x() * v2.z(),
                  v1.x() * v2.y() - v1.y() * v2.x());
}



Vector3f Eigen33::unitOrthogonal(const Vector3f &src)
{
  Vector3f perp;
  /* Let us compute the crossed product of *this with a vector
  * that is not too close to being colinear to *this.
  */

  /* unless the x and y coords are both close to zero, we can
  * simply take ( -y, x, 0 ) and normalize it.
  */
  if(!isMuchSmallerThan(src.x(), src.z()) || !isMuchSmallerThan(src.y(), src.z()))
  {   
    float invnm = 1.0 / sqrtf(src.x() * src.x() + src.y() * src.y());
    perp.x() = -src.y() * invnm;
    perp.y() =  src.x() * invnm;
    perp.z() = 0.0f;
  }   
  /* if both x and y are close to zero, then the vector is close
  * to the z-axis, so it's far from colinear to the x-axis for instance.
  * So we take the crossed product with (1,0,0) and normalize it. 
  */
  else
  {
    float invnm = 1.0 / sqrtf(src.z() * src.z() + src.y() * src.y());
    perp.x() = 0.0f;
    perp.y() = -src.z() * invnm;
    perp.z() =  src.y() * invnm;
  }   

  return perp;
}

void Eigen33::computeRoots2(const float& b, const float& c, Vector3f& roots)
{
 roots.x() = 0.f;
 float d = b * b - 4.f * c;
 if (d < 0.f) // no real roots!!!! THIS SHOULD NOT HAPPEN!
   d = 0.f;

 float sd = sqrtf(d);

 roots.z() = 0.5f * (b + sd);
 roots.y() = 0.5f * (b - sd);
}

void Eigen33::computeRoots3(float c0, float c1, float c2, Vector3f& roots)
{
 if (fabsf(c0) < std::numeric_limits<float>::epsilon())// one root is 0 -> quadratic equation
 {
   computeRoots2 (c2, c1, roots);
 }
 else
 {
   const float s_inv3 = 1.f/3.f;
   const float s_sqrt3 = sqrtf(3.f);
   // Construct the parameters used in classifying the roots of the equation
   // and in solving the equation for the roots in closed form.
   float c2_over_3 = c2 * s_inv3;
   float a_over_3 = (c1 - c2*c2_over_3)*s_inv3;
   if (a_over_3 > 0.f)
     a_over_3 = 0.f;

   float half_b = 0.5f * (c0 + c2_over_3 * (2.f * c2_over_3 * c2_over_3 - c1));

   float q = half_b * half_b + a_over_3 * a_over_3 * a_over_3;
   if (q > 0.f)
     q = 0.f;

   // Compute the eigenvalues by solving for the roots of the polynomial.
   float rho = sqrtf(-a_over_3);
   float theta = atan2f(sqrtf(-q), half_b)*s_inv3;
   float cos_theta = cosf(theta);
   float sin_theta = sinf(theta);
   roots.x() = c2_over_3 + 2.f * rho * cos_theta;
   roots.y() = c2_over_3 - rho * (cos_theta + s_sqrt3 * sin_theta);
   roots.z() = c2_over_3 - rho * (cos_theta - s_sqrt3 * sin_theta);

   // Sort in increasing order.
   if (roots.x() >= roots.y())
     std::swap(roots.x(), roots.y());

   if (roots.y() >= roots.z())
   {
     std::swap(roots.y(), roots.z());

     if (roots.x() >= roots.y())
       std::swap(roots.x(), roots.y());
   }
   if (roots.x() <= 0) // eigenval for symetric positive semi-definite matrix can not be negative! Set it to 0
     computeRoots2 (c2, c1, roots);
 }
}

void Eigen33::compute(Matrix3frm& tmp, Matrix3frm& vec_tmp, Matrix3frm& evecs, Vector3f& evals)
{
  // Scale the matrix so its entries are in [-1,1].  The scaling is applied
  // only when at least one matrix entry has magnitude larger than 1.

  float max01 = fmaxf( fabsf(mat_pkg[0]), fabsf(mat_pkg[1]) );
  float max23 = fmaxf( fabsf(mat_pkg[2]), fabsf(mat_pkg[3]) );
  float max45 = fmaxf( fabsf(mat_pkg[4]), fabsf(mat_pkg[5]) );
  float m0123 = fmaxf( max01, max23);
  float scale = fmaxf( max45, m0123);

  if (scale <= std::numeric_limits<float>::min())
   scale = 1.f;

  mat_pkg[0] /= scale;
  mat_pkg[1] /= scale;
  mat_pkg[2] /= scale;
  mat_pkg[3] /= scale;
  mat_pkg[4] /= scale;
  mat_pkg[5] /= scale;

  // The characteristic equation is x^3 - c2*x^2 + c1*x - c0 = 0.  The
  // eigenvalues are the roots to this equation, all guaranteed to be
  // real-valued, because the matrix is symmetric.
  float c0 = m00() * m11() * m22() 
           + m01() * m02() * m12() * 2.0
           - m00() * m12() * m12()
           - m11() * m02() * m02()
           - m22() * m01() * m01();

  float c1 =  m00() * m11() -
              m01() * m01() +
              m00() * m22() -
              m02() * m02() +
              m11() * m22() -
              m12() * m12();

  float c2 =  m00() + m11() + m22();

  computeRoots3(c0, c1, c2, evals);

  if(evals.z() - evals.x() <= std::numeric_limits<float>::epsilon()) {
    evecs.row(0) = Vector3f(1.f, 0.f, 0.f);
    evecs.row(1) = Vector3f(0.f, 1.f, 0.f);
    evecs.row(2) = Vector3f(0.f, 0.f, 1.f);
  }
  else if(evals.y() - evals.x() <= std::numeric_limits<float>::epsilon()) {
    // first and second equal
    tmp.row(0) = row0();  tmp.row(1) = row1();  tmp.row(2) = row2();
    tmp.row(0).x() -= evals.z(); tmp.row(1).y() -= evals.z(); tmp.row(0).z() -= evals.z();

    vec_tmp.row(0) = tmp.row(0).cross(tmp.row(1));
    vec_tmp.row(1) = tmp.row(0).cross(tmp.row(2));
    vec_tmp.row(2) = tmp.row(1).cross(tmp.row(2));

    float len1 = vec_tmp.row(0).dot(vec_tmp.row(0));
    float len2 = vec_tmp.row(1).dot(vec_tmp.row(1));
    float len3 = vec_tmp.row(2).dot(vec_tmp.row(2));

    if (len1 >= len2 && len1 >= len3)
    {
      evecs.row(2) = vec_tmp.row(0) / sqrtf(len1);
    }
    else if (len2 >= len1 && len2 >= len3)
    {
      evecs.row(2) = vec_tmp.row(1) / sqrtf(len2);
    }
    else
    {
      evecs.row(2) = vec_tmp.row(2) / sqrtf(len3);
    }

    evecs.row(1) = unitOrthogonal(evecs.row(2));
    evecs.row(0) = evecs.row(1).cross(evecs.row(2));
  }
  else if (evals.z() - evals.y() <= std::numeric_limits<float>::epsilon())
  {
    // second and third equal
    tmp.row(0) = row0();  tmp.row(1) = row1();  tmp.row(2) = row2();
    tmp.row(0).x() -= evals.x(); tmp.row(1).y() -= evals.x(); tmp.row(0).z() -= evals.x();

    vec_tmp.row(0) = tmp.row(0).cross(tmp.row(1));
    vec_tmp.row(1) = tmp.row(0).cross(tmp.row(2));
    vec_tmp.row(2) = tmp.row(1).cross(tmp.row(2));

    float len1 = vec_tmp.row(0).dot(vec_tmp.row(0));
    float len2 = vec_tmp.row(1).dot(vec_tmp.row(1));
    float len3 = vec_tmp.row(2).dot(vec_tmp.row(2));

    if (len1 >= len2 && len1 >= len3)
      evecs.row(0) = vec_tmp.row(0) / sqrtf(len1);

    else if (len2 >= len1 && len2 >= len3)
      evecs.row(0) = vec_tmp.row(1) / sqrtf(len2);

    else
      evecs.row(0) = vec_tmp.row(2) / sqrtf(len3);

    evecs.row(1) = unitOrthogonal(evecs.row(2));
    evecs.row(2) = evecs.row(0).cross(evecs.row(1));
  }
  else
  {
    tmp.row(0) = row0();  tmp.row(1) = row1();  tmp.row(2) = row2();
    tmp.row(0).x() -= evals.z(); tmp.row(1).y() -= evals.z(); tmp.row(2).z() -= evals.z();

    vec_tmp.row(0) = tmp.row(0).cross(tmp.row(1));
    vec_tmp.row(1) = tmp.row(0).cross(tmp.row(2));
    vec_tmp.row(2) = tmp.row(1).cross(tmp.row(2));

    float len1 = vec_tmp.row(0).dot(vec_tmp.row(0));
    float len2 = vec_tmp.row(1).dot(vec_tmp.row(1));
    float len3 = vec_tmp.row(2).dot(vec_tmp.row(2));

    float mmax[3];

    unsigned int min_el = 2;
    unsigned int max_el = 2;
    if (len1 >= len2 && len1 >= len3)
    {
      mmax[2] = len1;
      evecs.row(2) = vec_tmp.row(0) / sqrtf(len1);
    }
    else if (len2 >= len1 && len2 >= len3)
    {
      mmax[2] = len2;
      evecs.row(2) = vec_tmp.row(1) / sqrtf(len2);
    }
    else
    {
      mmax[2] = len3;
      evecs.row(2) = vec_tmp.row(2) / sqrtf(len3);
    }

    tmp.row(0) = row0();  tmp.row(1) = row1();  tmp.row(2) = row2();
    tmp.row(0).x() -= evals.y(); tmp.row(1).y() -= evals.y(); tmp.row(2).z() -= evals.y();

//      vec_tmp[0] = cross(tmp[0], tmp[1]);
//      vec_tmp[1] = cross(tmp[0], tmp[2]);
//      vec_tmp[2] = cross(tmp[1], tmp[2]);
//
//      len1 = dot(vec_tmp[0], vec_tmp[0]);
//      len2 = dot(vec_tmp[1], vec_tmp[1]);
//      len3 = dot(vec_tmp[2], vec_tmp[2]);

    if (len1 >= len2 && len1 >= len3)
    {
      mmax[1] = len1;
      evecs.row(1) = vec_tmp.row(0) / sqrtf(len1);
      min_el = len1 <= mmax[min_el] ? 1 : min_el;
      max_el = len1  > mmax[max_el] ? 1 : max_el;
    }
    else if (len2 >= len1 && len2 >= len3)
    {
      mmax[1] = len2;
      evecs.row(1) = vec_tmp.row(1) / sqrtf(len2);
      min_el = len2 <= mmax[min_el] ? 1 : min_el;
      max_el = len2  > mmax[max_el] ? 1 : max_el;
    }
    else
    {
      mmax[1] = len3;
      evecs.row(1) = vec_tmp.row(2) / sqrtf(len3);
      min_el = len3 <= mmax[min_el] ? 1 : min_el;
      max_el = len3 >  mmax[max_el] ? 1 : max_el;
    }

//      tmp[0] = row0();  tmp[1] = row1();  tmp[2] = row2();
//      tmp[0].x -= evals.x; tmp[1].y -= evals.x; tmp[2].z -= evals.x;
//
//      vec_tmp[0] = cross(tmp[0], tmp[1]);
//      vec_tmp[1] = cross(tmp[0], tmp[2]);
//      vec_tmp[2] = cross(tmp[1], tmp[2]);
//
//      len1 = dot (vec_tmp[0], vec_tmp[0]);
//      len2 = dot (vec_tmp[1], vec_tmp[1]);
//      len3 = dot (vec_tmp[2], vec_tmp[2]);


    if (len1 >= len2 && len1 >= len3)
    {
      mmax[0] = len1;
      evecs.row(0) = vec_tmp.row(0) / sqrtf(len1);
      min_el = len3 <= mmax[min_el] ? 0 : min_el;
      max_el = len3  > mmax[max_el] ? 0 : max_el;
    }
    else if (len2 >= len1 && len2 >= len3)
    {
      mmax[0] = len2;
      evecs.row(0) = vec_tmp.row(1) / sqrtf(len2);
      min_el = len3 <= mmax[min_el] ? 0 : min_el;
      max_el = len3  > mmax[max_el] ? 0 : max_el;
    }
    else
    {
      mmax[0] = len3;
      evecs.row(0) = vec_tmp.row(2) / sqrtf(len3);
      min_el = len3 <= mmax[min_el] ? 0 : min_el;
      max_el = len3  > mmax[max_el] ? 0 : max_el;
    }

    unsigned mid_el = 3 - min_el - max_el;
    evecs.row(min_el) = (evecs.row((min_el+1) % 3).cross(evecs.row((min_el+2) % 3)));
    evecs.row(min_el).normalize();
    evecs.row(mid_el) = (evecs.row((mid_el+1) % 3).cross(evecs.row((mid_el+2) % 3)));
    evecs.row(mid_el).normalize();
  }
  // Rescale back to the original size.
  evals *= scale;
}



bool Eigen33::isMuchSmallerThan(float x, float y)
{
  // copied from <eigen>/include/Eigen/src/Core/NumTraits.h
  const float prec_sqr = std::numeric_limits<float>::epsilon() * std::numeric_limits<float>::epsilon();
  return x * x <= prec_sqr * y * y;
}
