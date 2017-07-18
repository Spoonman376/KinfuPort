//
//

#ifndef Utils_h
#define Utils_h


struct float3
{
  float x, y, z;

  float3()
  {
    x = y = z = 0.0;
  }

  float3(float a, float b, float c)
  {
    x = a;
    y = b;
    z = c;
  }
};

float dot(const float3& v1, const float3& v2)
{
  return v1.x * v2.x + v1.y*v2.y + v1.z*v2.z;
}

float3& operator +=(float3& vec, const float& v)
{
  vec.x += v;  vec.y += v;  vec.z += v; return vec;
}

float3 operator +(const float3& v1, const float3& v2)
{
  return float3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

float3& operator *=(float3& vec, const float& v)
{
  vec.x *= v;  vec.y *= v;  vec.z *= v; return vec;
}

float3 operator -(const float3& v1, const float3& v2)
{
  return float3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

float3 operator *(const float3& v1, const float& v)
{
  return float3(v1.x * v, v1.y * v, v1.z * v);
}

float3 operator /(const float3& v1, const float& v)
{
  return float3(v1.x / v, v1.y / v, v1.z / v);
}

float norm(const float3& v)
{
  return sqrt(dot(v, v));
}

float3 normalized(const float3& v)
{
  return v / sqrt(dot(v, v));
}

float3 cross(const float3& v1, const float3& v2)
{
  return float3(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}


void computeRoots2(const float& b, const float& c, float3& roots)
{
 roots.x = 0.f;
 float d = b * b - 4.f * c;
 if (d < 0.f) // no real roots!!!! THIS SHOULD NOT HAPPEN!
   d = 0.f;

 float sd = sqrtf(d);

 roots.z = 0.5f * (b + sd);
 roots.y = 0.5f * (b - sd);
}

void computeRoots3(float c0, float c1, float c2, float3& roots)
{
 if ( fabsf(c0) < std::numeric_limits<float>::epsilon())// one root is 0 -> quadratic equation
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
   roots.x = c2_over_3 + 2.f * rho * cos_theta;
   roots.y = c2_over_3 - rho * (cos_theta + s_sqrt3 * sin_theta);
   roots.z = c2_over_3 - rho * (cos_theta - s_sqrt3 * sin_theta);

   // Sort in increasing order.
   if (roots.x >= roots.y)
     std::swap(roots.x, roots.y);

   if (roots.y >= roots.z)
   {
     std::swap(roots.y, roots.z);

     if (roots.x >= roots.y)
       std::swap(roots.x, roots.y);
   }
   if (roots.x <= 0) // eigenval for symetric positive semi-definite matrix can not be negative! Set it to 0
     computeRoots2 (c2, c1, roots);
 }
}

struct Eigen33
{
public:
  template<int Rows>
  struct MiniMat
  {
    float3 data[Rows];
    float3& operator[](int i) { return data[i]; }
    const float3& operator[](int i) const { return data[i]; }
  };
  typedef MiniMat<3> Mat33;
  typedef MiniMat<4> Mat43;




 
  static float3 unitOrthogonal (const float3& src)
  {
    float3 perp;
    /* Let us compute the crossed product of *this with a vector
    * that is not too close to being colinear to *this.
    */

    /* unless the x and y coords are both close to zero, we can
    * simply take ( -y, x, 0 ) and normalize it.
    */
    if(!isMuchSmallerThan(src.x, src.z) || !isMuchSmallerThan(src.y, src.z))
    {   
      float invnm = 1.0 / sqrtf(src.x*src.x + src.y*src.y);
      perp.x = -src.y * invnm;
      perp.y =  src.x * invnm;
      perp.z = 0.0f;
    }   
    /* if both x and y are close to zero, then the vector is close
    * to the z-axis, so it's far from colinear to the x-axis for instance.
    * So we take the crossed product with (1,0,0) and normalize it. 
    */
    else
    {
      float invnm = 1.0 / sqrtf(src.z * src.z + src.y * src.y);
      perp.x = 0.0f;
      perp.y = -src.z * invnm;
      perp.z =  src.y * invnm;
    }   

    return perp;
  }

  Eigen33(volatile float* mat_pkg_arg) : mat_pkg(mat_pkg_arg) {}


  void compute(Mat33& tmp, Mat33& vec_tmp, Mat33& evecs, float3& evals)
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

    if(evals.z - evals.x <= std::numeric_limits<float>::epsilon()) {
      evecs[0] = float3(1.f, 0.f, 0.f);
      evecs[1] = float3(0.f, 1.f, 0.f);
      evecs[2] = float3(0.f, 0.f, 1.f);
    }
    else if(evals.y - evals.x <= std::numeric_limits<float>::epsilon()) {
      // first and second equal
      tmp[0] = row0();  tmp[1] = row1();  tmp[2] = row2();
      tmp[0].x -= evals.z; tmp[1].y -= evals.z; tmp[2].z -= evals.z;

      vec_tmp[0] = cross(tmp[0], tmp[1]);
      vec_tmp[1] = cross(tmp[0], tmp[2]);
      vec_tmp[2] = cross(tmp[1], tmp[2]);

      float len1 = dot (vec_tmp[0], vec_tmp[0]);
      float len2 = dot (vec_tmp[1], vec_tmp[1]);
      float len3 = dot (vec_tmp[2], vec_tmp[2]);

      if (len1 >= len2 && len1 >= len3)
      {
        evecs[2] = vec_tmp[0] / sqrtf(len1);
      }
      else if (len2 >= len1 && len2 >= len3)
      {
        evecs[2] = vec_tmp[1] / sqrtf(len2);
      }
      else
      {
        evecs[2] = vec_tmp[2] / sqrtf(len3);
      }

      evecs[1] = unitOrthogonal(evecs[2]);
      evecs[0] = cross(evecs[1], evecs[2]);
    }
    else if (evals.z - evals.y <= std::numeric_limits<float>::epsilon() )
    {
      // second and third equal
      tmp[0] = row0();  tmp[1] = row1();  tmp[2] = row2();
      tmp[0].x -= evals.x; tmp[1].y -= evals.x; tmp[2].z -= evals.x;

      vec_tmp[0] = cross(tmp[0], tmp[1]);
      vec_tmp[1] = cross(tmp[0], tmp[2]);
      vec_tmp[2] = cross(tmp[1], tmp[2]);

      float len1 = dot(vec_tmp[0], vec_tmp[0]);
      float len2 = dot(vec_tmp[1], vec_tmp[1]);
      float len3 = dot(vec_tmp[2], vec_tmp[2]);

      if (len1 >= len2 && len1 >= len3)
        evecs[0] = vec_tmp[0] / sqrtf(len1);

      else if (len2 >= len1 && len2 >= len3)
        evecs[0] = vec_tmp[1] / sqrtf(len2);

      else
        evecs[0] = vec_tmp[2] / sqrtf(len3);

      evecs[1] = unitOrthogonal( evecs[0] );
      evecs[2] = cross(evecs[0], evecs[1]);
    }
    else
    {

      tmp[0] = row0();  tmp[1] = row1();  tmp[2] = row2();
      tmp[0].x -= evals.z; tmp[1].y -= evals.z; tmp[2].z -= evals.z;

      vec_tmp[0] = cross(tmp[0], tmp[1]);
      vec_tmp[1] = cross(tmp[0], tmp[2]);
      vec_tmp[2] = cross(tmp[1], tmp[2]);

      float len1 = dot(vec_tmp[0], vec_tmp[0]);
      float len2 = dot(vec_tmp[1], vec_tmp[1]);
      float len3 = dot(vec_tmp[2], vec_tmp[2]);

      float mmax[3];

      unsigned int min_el = 2;
      unsigned int max_el = 2;
      if (len1 >= len2 && len1 >= len3)
      {
        mmax[2] = len1;
        evecs[2] = vec_tmp[0] / sqrtf(len1);
      }
      else if (len2 >= len1 && len2 >= len3)
      {
        mmax[2] = len2;
        evecs[2] = vec_tmp[1] / sqrtf(len2);
      }
      else
      {
        mmax[2] = len3;
        evecs[2] = vec_tmp[2] / sqrtf(len3);
      }

      tmp[0] = row0();  tmp[1] = row1();  tmp[2] = row2();
      tmp[0].x -= evals.y; tmp[1].y -= evals.y; tmp[2].z -= evals.y;

      vec_tmp[0] = cross(tmp[0], tmp[1]);
      vec_tmp[1] = cross(tmp[0], tmp[2]);
      vec_tmp[2] = cross(tmp[1], tmp[2]);

      len1 = dot(vec_tmp[0], vec_tmp[0]);
      len2 = dot(vec_tmp[1], vec_tmp[1]);
      len3 = dot(vec_tmp[2], vec_tmp[2]);

      if (len1 >= len2 && len1 >= len3)
      {
        mmax[1] = len1;
        evecs[1] = vec_tmp[0] / sqrtf(len1);
        min_el = len1 <= mmax[min_el] ? 1 : min_el;
        max_el = len1  > mmax[max_el] ? 1 : max_el;
      }
      else if (len2 >= len1 && len2 >= len3)
      {
        mmax[1] = len2;
        evecs[1] = vec_tmp[1] / sqrtf(len2);
        min_el = len2 <= mmax[min_el] ? 1 : min_el;
        max_el = len2  > mmax[max_el] ? 1 : max_el;
      }
      else
      {
        mmax[1] = len3;
        evecs[1] = vec_tmp[2] / sqrtf(len3);
        min_el = len3 <= mmax[min_el] ? 1 : min_el;
        max_el = len3 >  mmax[max_el] ? 1 : max_el;
      }

      tmp[0] = row0();  tmp[1] = row1();  tmp[2] = row2();
      tmp[0].x -= evals.x; tmp[1].y -= evals.x; tmp[2].z -= evals.x;

      vec_tmp[0] = cross(tmp[0], tmp[1]);
      vec_tmp[1] = cross(tmp[0], tmp[2]);
      vec_tmp[2] = cross(tmp[1], tmp[2]);

      len1 = dot (vec_tmp[0], vec_tmp[0]);
      len2 = dot (vec_tmp[1], vec_tmp[1]);
      len3 = dot (vec_tmp[2], vec_tmp[2]);


      if (len1 >= len2 && len1 >= len3)
      {
        mmax[0] = len1;
        evecs[0] = vec_tmp[0] / sqrtf(len1);
        min_el = len3 <= mmax[min_el] ? 0 : min_el;
        max_el = len3  > mmax[max_el] ? 0 : max_el;
      }
      else if (len2 >= len1 && len2 >= len3)
      {
        mmax[0] = len2;
        evecs[0] = vec_tmp[1] / sqrtf(len2);
        min_el = len3 <= mmax[min_el] ? 0 : min_el;
        max_el = len3  > mmax[max_el] ? 0 : max_el;
      }
      else
      {
        mmax[0] = len3;
        evecs[0] = vec_tmp[2] / sqrtf(len3);
        min_el = len3 <= mmax[min_el] ? 0 : min_el;
        max_el = len3  > mmax[max_el] ? 0 : max_el;
      }

      unsigned mid_el = 3 - min_el - max_el;
      evecs[min_el] = normalized( cross( evecs[(min_el+1) % 3], evecs[(min_el+2) % 3] ) );
      evecs[mid_el] = normalized( cross( evecs[(mid_el+1) % 3], evecs[(mid_el+2) % 3] ) );
    }
    // Rescale back to the original size.
    evals *= scale;
  }
private:
  volatile float* mat_pkg;

  float m00() const { return mat_pkg[0]; }
  float m01() const { return mat_pkg[1]; }
  float m02() const { return mat_pkg[2]; }
  float m10() const { return mat_pkg[1]; }
  float m11() const { return mat_pkg[3]; }
  float m12() const { return mat_pkg[4]; }
  float m20() const { return mat_pkg[2]; }
  float m21() const { return mat_pkg[4]; }
  float m22() const { return mat_pkg[5]; }

  float3 row0() const { return float3( m00(), m01(), m02() ); }
  float3 row1() const { return float3( m10(), m11(), m12() ); }
  float3 row2() const { return float3( m20(), m21(), m22() ); }

  static bool isMuchSmallerThan (float x, float y)
  {
    // copied from <eigen>/include/Eigen/src/Core/NumTraits.h
    const float prec_sqr = std::numeric_limits<float>::epsilon() * std::numeric_limits<float>::epsilon();
    return x * x <= prec_sqr * y * y;
  }
};


float3 operator *(const Eigen33::Mat33& m, const float3& v)
{
  float3 r;
  r.x = dot(m.data[0], v);
  r.y = dot(m.data[1], v);
  r.z = dot(m.data[2], v);
  return r;
}








#endif /* Utils_h */
