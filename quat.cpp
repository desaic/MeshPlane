/**
 * @file quat.cpp
 * @brief A quaternion class.
 *
 * @author Eric Butler (edbutler)
 * @author Zeyang Li (zeyangl)
 */

#include "quat.h"
#include "mat.hpp"

const Quat Quat::Zero(0.0,0.0,0.0,0.0);

const Quat Quat::Identity(1.0,0.0,0.0,0.0);

Quat::Quat(const Vec3& axis, real_t radians)
{
    Vec3 naxis = axis.unit();
    radians *= 0.5;
    real_t sine = sin(radians);
    w = cos(radians);
    x = sine * naxis[0];
    y = sine * naxis[1];
    z = sine * naxis[2];
    normalize();
}

Quat::Quat(const Mat4& mat)
{
    // Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
    // article "Quaternion Calculus and Fast Animation".

    real_t trace = mat._m[0][0] + mat._m[1][1] + mat._m[2][2];
    real_t root;

    if ( trace > 0.0 )
    {
        // |w| > 1/2, may as well choose w > 1/2
        root = sqrt(trace + 1.0);  // 2w
        w = 0.5 * root;
        root = 0.5 / root;  // 1/(4w)
        x = (mat._m[2][1] - mat._m[1][2]) * root;
        y = (mat._m[0][2] - mat._m[2][0]) * root;
        z = (mat._m[1][0] - mat._m[0][1]) * root;
    }
    else
    {
        // |w| <= 1/2
        static size_t next[3] = { 1, 2, 0 };
        size_t i = 0;
        if ( mat._m[1][1] > mat._m[0][0] )
            i = 1;
        if ( mat._m[2][2] > mat._m[i][i] )
            i = 2;
        size_t j = next[i];
        size_t k = next[j];

        root = sqrt(mat._m[i][i] - mat._m[j][j] - mat._m[k][k] + 1.0);
        x = 0.5 * root;
        root = 0.5 / root;
        w = (mat._m[k][j] - mat._m[j][k]) * root;
        y = (mat._m[j][i] + mat._m[i][j]) * root;
        z = (mat._m[k][i] + mat._m[i][k]) * root;
    }
}

Quat Quat::operator*(const Quat& rhs) const
{
    return Quat(
        w * rhs.w - x * rhs.x - y * rhs.y - z * rhs.z,
        w * rhs.x + x * rhs.w + y * rhs.z - z * rhs.y,
        w * rhs.y + y * rhs.w + z * rhs.x - x * rhs.z,
        w * rhs.z + z * rhs.w + x * rhs.y - y * rhs.x
    );
}

Vec3 Quat::operator*(const Vec3& v) const
{
    // nVidia SDK implementation
    Vec3 qvec(x, y, z);
    Vec3 uv = qvec.cross(v);
    Vec3 uuv = qvec.cross(uv);
    uv *= (2.0 * w);
    uuv *= 2.0;

    return v + uv + uuv;
}

Quat& Quat::normalize()
{
    real_t mag = magnitude();
    x /= mag;
    y /= mag;
    z /= mag;
    w /= mag;
    return *this;
}

Quat Quat::inverse() const
{
    real_t norm = this->norm();
    if (norm > 0) {
        real_t inv = 1.0 / norm;
        return Quat( w*inv, -x*inv, -y*inv, -z*inv );
    } else {
        return Zero;
    }
}

void Quat::to_angle_axis(Vec3& axis, real_t* angle) const
{
    // The quaternion representing the rotation is
    // q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)
    real_t norm = x*x + y*y + z*z;
    if ( norm > 0.0 ) {
        *angle = 2.0 * acos(w);
        real_t inverse_length = 1 / sqrt(norm);
        axis[0] = x * inverse_length;
        axis[1] = y * inverse_length;
        axis[2] = z * inverse_length;
    } else {
        // angle is 0 (mod 2*pi), so any axis will do
        *angle = 0.0;
        axis = Vec3(1,0,0);
    }
}

void Quat::to_matrix(Mat4& mat) const
{
    real_t x2  = 2.0 * x;
    real_t y2  = 2.0 * y;
    real_t z2  = 2.0 * z;
    real_t xw2 = x2 * w;
    real_t yw2 = y2 * w;
    real_t zw2 = z2 * w;
    real_t xx2 = x2 * x;
    real_t xy2 = y2 * x;
    real_t xz2 = z2 * x;
    real_t yy2 = y2 * y;
    real_t yz2 = z2 * y;
    real_t zz2 = z2 * z;

    mat(0,0) = 1.0 - (yy2 + zz2);
    mat(0,1) = xy2 + zw2;
    mat(0,2) = xz2 - yw2;
    mat(0,3) = 0;

    mat(1,0) = xy2 - zw2;
    mat(1,1) = 1.0 - (xx2 + zz2);
    mat(1,2) = yz2 + xw2;
    mat(1,3) = 0;

    mat(2,0) = xz2 + yw2;
    mat(2,1) = yz2 - xw2;
    mat(2,2) = 1.0 - (xx2 + yy2);
    mat(2,3) = 0;

    mat(3,0) = 0;
    mat(3,1) = 0;
    mat(3,2) = 0;
    mat(3,3) = 1;
}

void Quat::to_axes(Vec3* axes) const
{
    real_t x2  = 2.0 * x;
    real_t y2  = 2.0 * y;
    real_t z2  = 2.0 * z;
    real_t xw2 = x2 * w;
    real_t yw2 = y2 * w;
    real_t zw2 = z2 * w;
    real_t xx2 = x2 * x;
    real_t xy2 = y2 * x;
    real_t xz2 = z2 * x;
    real_t yy2 = y2 * y;
    real_t yz2 = z2 * y;
    real_t zz2 = z2 * z;

    axes[0] = Vec3( 1.0 - (yy2 + zz2), xy2 + zw2, xz2 - yw2 );
    axes[1] = Vec3( xy2 - zw2, 1.0 - (xx2 + zz2), yz2 + xw2 );
    axes[2] = Vec3( xz2 + yw2, yz2 - xw2, 1.0 - (xx2 + yy2) );
}

std::ostream& operator <<( std::ostream& o, const Quat& q )
{
    o << "Quat(" << q.w << ", " << q.x << ", " << q.y << ", " << q.z << ")";
    return o;
}
