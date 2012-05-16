#ifndef MAT_HPP
#define MAT_HPP
#include "mesh.hpp"

class Mat3
{
 public:
  /**
   * The identity matrix.
   */
  static const Mat3 Identity;

  /**
   * The zero matrix.
   */
  static const Mat3 Zero;

  static const int DIM = 3;
  static const int SIZE = DIM*DIM;
  /**@brief COLUMN Major format !!!!!*/
  explicit Mat3(real_t r[SIZE]);
  Mat3() {}
  Mat3(const Mat3 & in){for(int ii=0;ii<9;ii++){m[ii]=in.m[ii];}}
  /**@brief ROW Major format!!!!*/
  Mat3(real_t m00, real_t m10, real_t m20,
       real_t m01, real_t m11, real_t m21,
       real_t m02, real_t m12, real_t m22);
  /**@brief each ROW is passed in as a vector*/
  Mat3(const Vec3 & v1, const Vec3&v2,const Vec3 & v3);
  void inverse(Mat3& rv) const;
  //original matrix is destroyed
  void gauss_pivot(Mat3 & rv);
  union{
    real_t m[SIZE];
    real_t _m[DIM][DIM];
  };
  void transpose();
  void mult(const Mat3 & rv, Mat3& result);
  /**@TODO: not working yet, need to figure out
  how to use permutation matrix p and q*/
  void invert_full(Mat3 & rv);
private:
  void h_solve(real_t *x, int *p, int *q);
  void h_pivot_decomp(int *p, int *q);
};

class Vec4
{
public:
    /**
     * The zero vector.
     */
    static const Vec4 Zero;

    /**
     * The vector (1,1,1,1).
     */
    static const Vec4 Ones;

    /**
     * The vector (1,0,0,0)
     */
    static const Vec4 UnitX;

    /**
     * The vector (0,1,0,0)
     */
    static const Vec4 UnitY;

    /**
     * The vector (0,0,1,0)
     */
    static const Vec4 UnitZ;

    /**
     * The vector (0,0,0,1)
     */
    static const Vec4 UnitW;

    /**
     * Components of this vector.
     */
    real_t x, y, z, w;

    /**
     * Default constructor. Leaves values unitialized.
     */
    Vec4() {}

    /**
     * Create a vector with the given values.
     */
    Vec4(real_t x, real_t y, real_t z, real_t w)
        :x(x), y(y), z(z), w(w) {}

    /**
     * Create the vector (v.x, v.y, v.z, w).
     */
    Vec4(const Vec3& v, real_t w)
        :x(v.x[0]), y(v.x[1]), z(v.x[2]), w(w) {}

    /**
     * Create a vector from the elements of the given array.
     */
    explicit Vec4(const float p[4])
        : x(p[0]), y(p[1]), z(p[2]), w(p[3]) {}

    /**
     * Create a vector from the elements of the given array.
     */
    explicit Vec4(const double p[4])
        : x(p[0]), y(p[1]), z(p[2]), w(p[3]) {}

    // also uses default copy and assignment

    Vec4 operator+(const Vec4& rhs) const
    {
        return Vec4( x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.w );
    }

    Vec4& operator+=(const Vec4& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        w += rhs.w;
        return *this;
    }

    Vec4 operator-(const Vec4& rhs) const
    {
        return Vec4( x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w );
    }

    Vec4& operator-=(const Vec4& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        w -= rhs.w;
        return *this;
    }

    // not really mathematically sound, but makes sense for colors
    Vec4 operator*(const Vec4& rhs) const
    {
        return Vec4( x * rhs.x, y * rhs.y, z * rhs.z, w * rhs.w );
    }

    // not really mathematically sound, but makes sense for colors
    Vec4& operator*=(const Vec4& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        z *= rhs.z;
        w *= rhs.w;
        return *this;
    }

    Vec4 operator*(real_t s) const
    {
        return Vec4( x * s, y * s, z * s, w * s );
    }

    Vec4& operator*=(real_t s)
    {
        x *= s;
        y *= s;
        z *= s;
        w *= s;
        return *this;
    }

    // not really mathematically sound, but makes sense for colors
    Vec4 operator/(const Vec4& rhs) const
    {
        return Vec4( x / rhs.x, y / rhs.y, z / rhs.z, w / rhs.w );
    }

    // not really mathematically sound, but makes sense for colors
    Vec4& operator/=(const Vec4& rhs)
    {
        x /= rhs.x;
        y /= rhs.y;
        z /= rhs.z;
        w /= rhs.w;
        return *this;
    }

    Vec4 operator/(real_t s) const
    {
        return Vec4( x / s, y / s, z / s, w / s );
    }

    Vec4& operator/=(real_t s)
    {
        x /= s;
        y /= s;
        z /= s;
        w /= s;
        return *this;
    }

    Vec4 operator-() const
    {
        return Vec4( -x, -y, -z, -w);
    }

    /**
     * @remark No bounds checking.
     */
    real_t operator[](size_t i) const
    {
        return (&x)[i];
    }

    /**
     * @remark No bounds checking.
     */
    real_t& operator[](size_t i)
    {
        return (&x)[i];
    }

    bool operator==(const Vec4& rhs) const
    {
        return x == rhs.x && y == rhs.y && z == rhs.z && w == rhs.w;
    }

    bool operator!=(const Vec4& rhs) const
    {
        return !operator==(rhs);
    }

    /**
     * Returns the 3d vector corresponding to this 4d vector.
     * @remark If w==0, returns (x,y,z).
     */
    Vec3 projection() const
    {
        real_t tmpw = w == 0 ? 1 : w;
        return Vec3( x/tmpw, y/tmpw, z/tmpw );
    }

    /**
     * Returns the first three components, ignoring the fourth
     */
    Vec3 xyz() const
    {
        return Vec3(x, y, z);
    }

    /**
     * Returns the dot product of two vectors
     */
    real_t dot(const Vec4& rhs) const
    {
        return x*rhs.x + y*rhs.y + z*rhs.z + w*rhs.w;
    }

    /**
     * Returns the magnitude of a vector
     */
    real_t magnitude() const
    {
        return sqrt(x*x + y*y + z*z + w*w);
    }

    /**
     * Efficiency function: does not require square root operation.
     */
    real_t squared_magnitude() const
    {
        return x*x + y*y + z*z + w*w;
    }

    /**
     * Calculate the positive distance between two vectors.
     */
    real_t distance(const Vec4& rhs) const
    {
        return (*this-rhs).magnitude();
    }

    /**
     * Efficiency function: does not require square root operation.
     */
    real_t squared_distance(const Vec4& rhs) const
    {
        return (*this-rhs).squared_magnitude();
    }

    /**
     * Returns the unit vector pointing in the same direction as this vector.
     */
    Vec4 unit() const
    {
        return *this / magnitude();
    }

    /**
     * Noramlizes this vector; that is, sets its magnitude to 1.
     */
    Vec4& normalize()
    {
        return *this /= magnitude();
    }

    /**
     * Returns a vector whose elements are the absolute values of all the
     * elements of this vector.
     */
    Vec4 abs() const
    {
        return Vec4(fabs(x), fabs(y), fabs(z), fabs(w));
    }

    /**
     * Returns a vector which is the point exactly between this and the given
     * vector.
     */
    Vec4 midpoint(const Vec4& rhs) const
    {
        return (*this + rhs) * .5;
    }

    /**
     * Clamps the lower bound of this vector; that is, sets this vector's values
     * to the max of the current values and the given vector's values.
     * @param rhs The vector representing the desired lower bound.
     */
    Vec4& clamp_min(const Vec4& rhs)
    {
        x = std::max(x, rhs.x);
        y = std::max(y, rhs.y);
        z = std::max(z, rhs.z);
        w = std::max(w, rhs.w);
        return *this;
    }

    /**
     * Clamps the upper bound of this vector; that is, sets this vector's values
     * to the min of the current values and the given vector's values.
     * @param rhs The vector representing the desired upper bound.
     */
    Vec4& clamp_max(const Vec4& rhs)
    {
        x = std::min(x, rhs.x);
        y = std::min(y, rhs.y);
        z = std::min(z, rhs.z);
        w = std::min(w, rhs.w);
        return *this;
    }

    /**
     * Returns a vector whose values are the maximum of this vector's values and
     * the given vector's values.
     */
    Vec4 maximum(const Vec4& rhs) const
    {
        return Vec4(
            std::max(x, rhs.x),
            std::max(y, rhs.y),
            std::max(z, rhs.z),
            std::max(w, rhs.w)
        );
    }

    /**
     * Returns a vector whose values are the minimum of this vector's values and
     * the given vector's values.
     */
    Vec4 minimum(const Vec4& rhs) const
    {
        return Vec4(
            std::min(x, rhs.x),
            std::min(y, rhs.y),
            std::min(z, rhs.z),
            std::min(w, rhs.w)
        );
    }

    void to_array(float arr[4]) const
    {
        arr[0] = float(x);
        arr[1] = float(y);
        arr[2] = float(z);
        arr[3] = float(w);
    }
};

inline Vec4 operator*(real_t s, const Vec4& rhs)
{
    return rhs * s;
}

/**
 * Outputs a vector text formatted as "(x,y,z)".
 */
std::ostream& operator<<(std::ostream& os, const Vec4& rhs);

/**
 * Reads in a vector text formatted as "(x,y,z)". Whitespaces are fine, as is
 * using '<' and '>' instead of '(' and ')'.
 */
std::istream& operator>>(std::istream& is, Vec4& rhs);

class Quat;
class Mat4
{
public:
    /**
     * The identity matrix.
     */
    static const Mat4 Identity;

    /**
     * The zero matrix.
     */
    static const Mat4 Zero;

    static const int DIM = 4;
    static const int SIZE = DIM*DIM;

    /**
     * The values of this matrix. Named as both a 1d array and a 2d array.
     */
    union{
        real_t m[SIZE];
        real_t _m[DIM][DIM]; // _m[column][row]
    };

    // constructors

    /**
     * Leaves values unitialized.
     */
    Mat4() {}

    /**
     * Construct a matrix from the given array, in COLUMN MAJOR format
     */
    explicit Mat4(real_t r[SIZE]);

    /**
     * Construct a matrix from the given values in ROW MAJOR format (A,B,C,D,E,F,G,H,I).
     */
    Mat4(real_t m00, real_t m10, real_t m20, real_t m30,
         real_t m01, real_t m11, real_t m21, real_t m31,
         real_t m02, real_t m12, real_t m22, real_t m32,
         real_t m03, real_t m13, real_t m23, real_t m33);

    // basic operations
    Mat4 operator+(const Mat4& rhs) const;
    Mat4& operator+=(const Mat4& rhs);
    Mat4 operator-(const Mat4& rhs) const;
    Mat4& operator-=(const Mat4& rhs);
    Mat4 operator*(const Mat4& rhs) const;
    Vec4 operator*(const Vec4& v) const;
    Mat4& operator*=(const Mat4& rhs);
    Mat4 operator*(real_t r) const;
    Mat4& operator*=(real_t r);
    Mat4 operator/(real_t r) const;
    Mat4& operator/=(real_t r);
    Mat4 operator-() const;

    // comparaisons
    bool operator==(const Mat4& rhs) const;
    bool operator!=(const Mat4& rhs) const;

    // accessors

    /**
     * Mat4(i,j) gives the element at the ith column and jth row.
     */
    real_t operator()(int col, int row) const
    {
        return _m[col][row];
    }

    /**
     * Mat4(i,j) gives the element at the ith column and jth row.
     */
    real_t& operator()(int col, int row)
    {
        return _m[col][row];
    }

    /**
     * Transform the given vector using this matrix
     */
    Vec4 transform(const Vec4& v) const
    {
        return operator*(v);
    }

    Vec3 transform_point(const Vec3& v) const
    {
        return transform(Vec4(v, 1)).projection();
    }

    Vec3 transform_vector(const Vec3& v) const
    {
        return transform(Vec4(v, 0)).projection();
    }

    /**
     * Combines two transformations into one, with this Matrix being the first
     * to be applied (rhs) and the given matrix the second (lhs). If this matrix
     * is A and the given lhs is B, this function is equivalent to A=B*A.
     */
    Mat4& concatenate(const Mat4& lhs)
    {
        return *this = lhs*(*this);
    }

    /**
     * Returns the transpose matrix.
     */
    void transpose(Mat4& rv) const;

    /**
     * Returns the inverse matrix. Is undefined if determinant is 0.
     */
    void inverse(Mat4& rv) const;
};

inline Mat4 operator*(real_t r, const Mat4& m)
{
    return m * r;
}

/**
 * Constructs a matrix that traslates vectors by the given position.
 */
void make_translation_matrix(Mat4& mat, const Vec3& pos);

/**
 * Constructs a matrix that rotates vectors by the given orientation.
 */
void make_rotation_matrix(Mat4& mat, const Quat& ori);

/**
 * Constructs a matrix that scales vectors by the given scale.
 */
void make_scaling_matrix(Mat4& mat, const Vec3& scl);

/**
 * Constructs a matrix that scales, rotates, then translates vectors.
 */
void make_transformation_matrix(Mat4& mat,
                                const Vec3& pos,
                                const Quat& ori,
                                const Vec3& scl);

/**
 * Constructs the matrix to transform normal vectors corresponding to the
 * given transformation matrix.
 */
void make_normal_matrix(Mat3& rv, const Mat4& tmat);

/**
 * Constructs the matrix to transform normal vectors corresponding to the
 * given transformation matrix.
 */
void make_normal_matrix(Mat3& rv, const Mat3& tmat);

#endif
