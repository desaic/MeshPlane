/**
 * @file quat.h
 * @brief A quaternion class.
 *
 * @author Eric Butler (edbutler)
 * @author Zeyang Li (zeyangl)
 */

#ifndef _VEC_QUAT_H_
#define _VEC_QUAT_H_


#include "vec.h"
#include "math.hpp"
class Mat4;

/*
   This code is loosely based on quaternion code from Ogre3d (www.ogre3d.org).
 */

class Quat
{
public:
    static const Quat Zero;
    static const Quat Identity;

    real_t w, x, y, z;

    Quat() {}

    /**
     * Construct a quaternion with the given values.
     */
    Quat(real_t w, real_t x, real_t y, real_t z):
    w(w),
    x(x),
    y(y),
    z(z)
    {

    }

    /**
     * Constructs a quaternion representing a rotation about the given axis
     * by the given angle.
     */
    Quat(const Vec3& axis, real_t radians);

    /**
     *  Constructs a quaternion from a rotation matrix.
     */
    Quat(const Mat4& mat);

    real_t operator[]( const size_t i ) const
    {
        return *(&w+i);
    }

    real_t& operator[]( const size_t i )
    {
        return *(&w+i);
    }

    Quat operator*(const Quat& rhs) const;

    Vec3 operator* (const Vec3& rhs) const;

    bool operator== (const Quat& rhs) const
    {
        return rhs.x == x && rhs.y == y &&
               rhs.z == z && rhs.w == w;
    }

    bool operator!= (const Quat& rhs) const
    {
        return !operator==(rhs);
    }

    real_t norm() const
    {
        return x*x + y*y + z*z + w*w;
    }

    real_t magnitude() const
    {
        return sqrt(norm());
    }

    Quat& normalize();

    /**
     * Returns the inverse of this quaternion.
     */
    Quat inverse() const;

    /**
     * Convert this quaternion into an angle and axis.
     * Returns the rotation in radians about an axis.
     */
    void to_angle_axis(Vec3& axis, real_t* angle) const;

    /**
     * Converts this quaternion to a matrix.
     */
    void to_matrix(Mat4& mat) const;

    /**
     * Returns the X,Y,Z axes rotated by this quaternion.
     */
    void to_axes(Vec3* axes) const;
};

std::ostream& operator <<( std::ostream& o, const Quat& q );

inline Vec3 operator*(const Vec3& v, const Quat& q)
{
    return q * v;
}

inline Vec3& operator*=(Vec3& v, const Quat& q)
{
    return v = q * v;
}

#endif /* _VEC_QUAT_H_ */

