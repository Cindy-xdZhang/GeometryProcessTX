#ifndef UTILS_QUATERNIONS_H
#define UTILS_QUATERNIONS_H

#include <cmath>
#include <cassert>

#include "utils/constants.h"


namespace quaternions
{
    ////////////////////////////////////////////////////////////////////////////
    /// Canonical Quaternions
    ////////////////////////////////////////////////////////////////////////////

    // Define some canonical quaternions for floats and doubles
    // A Quaternion, q, is defined here as an arrays of four scalars (x,y,z,w),
    // such that q = x*i + y*j + z*k + w

    // Float versions
#define SQRT_2_OVER_2 0.707106781f

    // Identity
    const float IDENTITY_QUAT_F[4] = {0, 0, 0, 1};

    // The following match the Matlab canonical views

    // X point right, Y pointing up and Z point out
    const float XY_PLANE_QUAT_F[4] = {0, 0, 0, 1};

    // X points right, Y points *in* and Z points up
    const float XZ_PLANE_QUAT_F[4] = {-SQRT_2_OVER_2, 0, 0, SQRT_2_OVER_2};

    // X points out, Y points right, and Z points up
    const float YZ_PLANE_QUAT_F[4] = {-0.5, -0.5, -0.5, 0.5};

    const float CANONICAL_VIEW_QUAT_F[][4] =
            {
                    {0,              0,              0,              1}, // 0
                    {0,              0,              SQRT_2_OVER_2, SQRT_2_OVER_2}, // 1
                    {0,              0,              1,              0}, // 2
                    {0,              0,              SQRT_2_OVER_2,  -SQRT_2_OVER_2}, // 3

                    {0,              -1,             0,              0}, // 4
                    {-SQRT_2_OVER_2, SQRT_2_OVER_2,  0,              0}, // 5
                    {-1,             0,              0,              0}, // 6
                    {-SQRT_2_OVER_2, -SQRT_2_OVER_2, 0,              0}, // 7

                    {-0.5,           -0.5,           -0.5,           0.5}, // 8
                    {0,              -SQRT_2_OVER_2, 0,             SQRT_2_OVER_2}, // 9
                    {0.5,            -0.5,           0.5,            0.5}, // 10
                    {SQRT_2_OVER_2,  0,              SQRT_2_OVER_2,  0}, // 11

                    {SQRT_2_OVER_2,  0,              -SQRT_2_OVER_2, 0}, // 12
                    {0.5,            0.5,            -0.5,           0.5}, // 13
                    {0,              SQRT_2_OVER_2,  0,             SQRT_2_OVER_2}, // 14
                    {-0.5,           0.5,            0.5,            0.5}, // 15

                    {0,              SQRT_2_OVER_2,  SQRT_2_OVER_2,  0}, // 16
                    {-0.5,           0.5,            0.5,            -0.5}, // 17
                    {-SQRT_2_OVER_2, 0,              0,              -SQRT_2_OVER_2}, // 18
                    {-0.5,           -0.5,           -0.5,           -0.5}, // 19

                    {-SQRT_2_OVER_2, 0,              0,             SQRT_2_OVER_2}, // 20
                    {-0.5,           -0.5,           0.5,            0.5}, // 21
                    {0,              -SQRT_2_OVER_2, SQRT_2_OVER_2,  0}, // 22
                    {0.5,            -0.5,           0.5,            -0.5}  // 23
            };

#undef SQRT_2_OVER_2

    // Double versions
#define SQRT_2_OVER_2 0.70710678118654757

    // Identity
    const double IDENTITY_QUAT_D[4] = {0, 0, 0, 1};

    // The following match the Matlab canonical views

    // X point right, Y pointing up and Z point out
    const double XY_PLANE_QUAT_D[4] = {0, 0, 0, 1};

    // X points right, Y points *in* and Z points up
    const double XZ_PLANE_QUAT_D[4] = {-SQRT_2_OVER_2, 0, 0, SQRT_2_OVER_2};

    // X points out, Y points right, and Z points up
    const double YZ_PLANE_QUAT_D[4] = {-0.5, -0.5, -0.5, 0.5};

    const double CANONICAL_VIEW_QUAT_D[][4] =
            {
                    {0,              0,              0,              1},
                    {0,              0,              SQRT_2_OVER_2, SQRT_2_OVER_2},
                    {0,              0,              1,              0},
                    {0,              0,              SQRT_2_OVER_2,  -SQRT_2_OVER_2},

                    {0,              -1,             0,              0},
                    {-SQRT_2_OVER_2, SQRT_2_OVER_2,  0,              0},
                    {-1,             0,              0,              0},
                    {-SQRT_2_OVER_2, -SQRT_2_OVER_2, 0,              0},

                    {-0.5,           -0.5,           -0.5,           0.5},
                    {0,              -SQRT_2_OVER_2, 0,             SQRT_2_OVER_2},
                    {0.5,            -0.5,           0.5,            0.5},
                    {SQRT_2_OVER_2,  0,              SQRT_2_OVER_2,  0},

                    {SQRT_2_OVER_2,  0,              -SQRT_2_OVER_2, 0},
                    {0.5,            0.5,            -0.5,           0.5},
                    {0,              SQRT_2_OVER_2,  0,             SQRT_2_OVER_2},
                    {-0.5,           0.5,            0.5,            0.5},

                    {0,              SQRT_2_OVER_2,  SQRT_2_OVER_2,  0},
                    {-0.5,           0.5,            0.5,            -0.5},
                    {-SQRT_2_OVER_2, 0,              0,              -SQRT_2_OVER_2},
                    {-0.5,           -0.5,           -0.5,           -0.5},

                    {-SQRT_2_OVER_2, 0,              0,             SQRT_2_OVER_2},
                    {-0.5,           -0.5,           0.5,            0.5},
                    {0,              -SQRT_2_OVER_2, SQRT_2_OVER_2,  0},
                    {0.5,            -0.5,           0.5,            -0.5}
            };

#undef SQRT_2_OVER_2

#define NUM_CANONICAL_VIEW_QUAT 24

    // NOTE: I want to rather be able to return a Q_type[][] but C++ is not
    // making it easy. So instead I've written a per-element accessor

    // Return element [i][j] of the corresponding CANONICAL_VIEW_QUAT_* of the
    // given templated type
    // Inputs:
    //   i  index of quaternion
    //   j  index of coordinate in quaternion i
    // Returns values of CANONICAL_VIEW_QUAT_*[i][j]
    template<typename Q_type>
    inline Q_type CANONICAL_VIEW_QUAT(int i, int j);

    // Template specializations for float and double
    template<>
    inline float CANONICAL_VIEW_QUAT<float>(int i, int j)
    {
        return (float) CANONICAL_VIEW_QUAT_F[i][j];
    }

    template<>
    inline double CANONICAL_VIEW_QUAT<double>(int i, int j)
    {
        return (double) CANONICAL_VIEW_QUAT_D[i][j];
    }

    ////////////////////////////////////////////////////////////////////////////
    /// Operators
    ////////////////////////////////////////////////////////////////////////////

    template<typename Q_type>
    inline void quat_conjugate(const Q_type* q1, Q_type* out)
    {
        out[0] = -q1[0];
        out[1] = -q1[1];
        out[2] = -q1[2];
        out[3] = q1[3];
    }

    template<typename Q_type>
    inline bool normalize_quat(const Q_type* q, Q_type* out)
    {
        // Get length
        Q_type len = sqrt(
                q[0] * q[0] +
                q[1] * q[1] +
                q[2] * q[2] +
                q[3] * q[3]);

        // Normalize each coordinate
        out[0] = q[0] / len;
        out[1] = q[1] / len;
        out[2] = q[2] / len;
        out[3] = q[3] / len;

        // Test whether length was below Epsilon
        return (len > utils::EPS<Q_type>());
    }

    template<typename Q_type>
    inline void quat_mult(const Q_type* q1, const Q_type* q2, Q_type* out)
    {
        // output can't be either of the inputs
        assert(q1 != out);
        assert(q2 != out);

        out[0] = q1[3] * q2[0] + q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1];
        out[1] = q1[3] * q2[1] + q1[1] * q2[3] + q1[2] * q2[0] - q1[0] * q2[2];
        out[2] = q1[3] * q2[2] + q1[2] * q2[3] + q1[0] * q2[1] - q1[1] * q2[0];
        out[3] = q1[3] * q2[3] - (q1[0] * q2[0] + q1[1] * q2[1] + q1[2] * q2[2]);
    }

    template<typename Q_type>
    inline void rotate_by_quat(const Q_type* v, const Q_type* q, Q_type* out)
    {
        // Quaternion form of v, copy data in v, (as a result out can be same pointer
        // as v)
        Q_type quat_v[4] = {v[0], v[1], v[2], 0};

        // normalize input
        Q_type normalized_q[4];

#ifndef NDEBUG
        bool normalized =
#endif
                normalize_quat<Q_type>(q, normalized_q);
#ifndef NDEBUG
        assert(normalized);
#endif

        // Conjugate of q
        Q_type q_conj[4];
        quat_conjugate<Q_type>(normalized_q, q_conj);

        // Rotate of vector v by quaternion q is:
        // q*v*conj(q)
        // Compute q*v
        Q_type q_mult_quat_v[4];
        quat_mult<Q_type>(normalized_q, quat_v, q_mult_quat_v);
        // Compute (q*v) * conj(q)
        quat_mult<Q_type>(q_mult_quat_v, q_conj, out);
    }

    ////////////////////////////////////////////////////////////////////////////
    /// Conversions
    ////////////////////////////////////////////////////////////////////////////

    // http://www.antisphere.com/Wiki/tools:anttweakbar

    template<typename Q_type>
    inline void axis_angle_to_quat(const Q_type* axis, const Q_type angle, Q_type* out)
    {
        Q_type n = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2];
        if (fabs(n) > utils::EPS<Q_type>())
        {
            Q_type f = 0.5 * angle;
            out[3] = cos(f);
            f = sin(f) / sqrt(n);
            out[0] = axis[0] * f;
            out[1] = axis[1] * f;
            out[2] = axis[2] * f;
        }
        else
        {
            out[3] = 1.0;
            out[0] = out[1] = out[2] = 0.0;
        }
    }

    template<typename Q_type>
    inline void quat_to_axis_angle(const Q_type* q, Q_type* axis, Q_type& angle)
    {
        if (fabs(q[3]) > (1.0 + utils::EPS<Q_type>()))
        {
            //axis[0] = axis[1] = axis[2] = 0; // no, keep the previous value
            angle = 0;
        }
        else
        {
            double a;
            if (q[3] >= 1.0f)
            {
                a = 0; // and keep V
            }
            else if (q[3] <= -1.0f)
            {
                a = utils::PI; // and keep V
            }
            else if (fabs(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]) < utils::EPS_SQ<Q_type>())
            {
                a = 0;
            }
            else
            {
                a = acos(q[3]);
                if (a * angle < 0)
                { // Preserve the sign of angle
                    a = -a;
                }
                double f = 1.0f / sin(a);
                axis[0] = q[0] * f;
                axis[1] = q[1] * f;
                axis[2] = q[2] * f;
            }
            angle = 2.0 * a;
        }

        //  if( angle>FLOAT_PI )
        //      angle -= 2.0f*FLOAT_PI;
        //  else if( angle<-FLOAT_PI )
        //      angle += 2.0f*FLOAT_PI;
        //angle = RadToDeg(angle);

        if (fabs(angle) < utils::EPS<Q_type>() &&
            fabs(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]) < utils::EPS_SQ<Q_type>())
        {
            axis[0] = 1.0e-7;    // all components cannot be null
        }
    }

    template<typename Q_type>
    inline void quat_to_axis_angle_deg(const Q_type* q, Q_type* axis, Q_type& angle)
    {
        quat_to_axis_angle(q, axis, angle);
        angle = angle * (180.0 / utils::PI);
    }

    // https://bmgame.googlecode.com/svn/idlib/math/Simd_AltiVec.cpp

    template<typename Q_type>
    inline void mat4_to_quat(const Q_type* mat, Q_type* q)
    {
        Q_type trace;
        Q_type s;
        Q_type t;
        int i;
        int j;
        int k;

        static int next[3] = {1, 2, 0};

        trace = mat[0 * 4 + 0] + mat[1 * 4 + 1] + mat[2 * 4 + 2];

        if (trace > 0.0f)
        {

            t = trace + 1.0f;
            s = ReciprocalSqrt(t) * 0.5f;

            q[3] = s * t;
            q[0] = (mat[1 * 4 + 2] - mat[2 * 4 + 1]) * s;
            q[1] = (mat[2 * 4 + 0] - mat[0 * 4 + 2]) * s;
            q[2] = (mat[0 * 4 + 1] - mat[1 * 4 + 0]) * s;

        }
        else
        {

            i = 0;
            if (mat[1 * 4 + 1] > mat[0 * 4 + 0])
            {
                i = 1;
            }
            if (mat[2 * 4 + 2] > mat[i * 4 + i])
            {
                i = 2;
            }
            j = next[i];
            k = next[j];

            t = (mat[i * 4 + i] - (mat[j * 4 + j] + mat[k * 4 + k])) + 1.0f;
            s = ReciprocalSqrt(t) * 0.5f;

            q[i] = s * t;
            q[3] = (mat[j * 4 + k] - mat[k * 4 + j]) * s;
            q[j] = (mat[i * 4 + j] + mat[j * 4 + i]) * s;
            q[k] = (mat[i * 4 + k] + mat[k * 4 + i]) * s;
        }

        //// Unused translation
        //jq.t[0] = mat[0 * 4 + 3];
        //jq.t[1] = mat[1 * 4 + 3];
        //jq.t[2] = mat[2 * 4 + 3];
    }

    template<typename Q_type>
    inline void mat3_to_quat(const Q_type* mat, Q_type* q)
    {
        Q_type trace;
        Q_type s;
        Q_type t;
        int i;
        int j;
        int k;

        static int next[3] = {1, 2, 0};

        trace = mat[0 * 3 + 0] + mat[1 * 3 + 1] + mat[2 * 3 + 2];

        if (trace > 0.0f)
        {

            t = trace + 1.0f;
            s = ReciprocalSqrt(t) * 0.5f;

            q[3] = s * t;
            q[0] = (mat[1 * 3 + 2] - mat[2 * 3 + 1]) * s;
            q[1] = (mat[2 * 3 + 0] - mat[0 * 3 + 2]) * s;
            q[2] = (mat[0 * 3 + 1] - mat[1 * 3 + 0]) * s;

        }
        else
        {

            i = 0;
            if (mat[1 * 3 + 1] > mat[0 * 3 + 0])
            {
                i = 1;
            }
            if (mat[2 * 3 + 2] > mat[i * 3 + i])
            {
                i = 2;
            }
            j = next[i];
            k = next[j];

            t = (mat[i * 3 + i] - (mat[j * 3 + j] + mat[k * 3 + k])) + 1.0f;
            s = ReciprocalSqrt(t) * 0.5f;

            q[i] = s * t;
            q[3] = (mat[j * 3 + k] - mat[k * 3 + j]) * s;
            q[j] = (mat[i * 3 + j] + mat[j * 3 + i]) * s;
            q[k] = (mat[i * 3 + k] + mat[k * 3 + i]) * s;
        }

        //// Unused translation
        //jq.t[0] = mat[0 * 4 + 3];
        //jq.t[1] = mat[1 * 4 + 3];
        //jq.t[2] = mat[2 * 4 + 3];
    }

    template<typename Q_type>
    inline void quat_to_mat(const Q_type* quat, Q_type* mat)
    {
        Q_type yy2 = 2.0f * quat[1] * quat[1];
        Q_type xy2 = 2.0f * quat[0] * quat[1];
        Q_type xz2 = 2.0f * quat[0] * quat[2];
        Q_type yz2 = 2.0f * quat[1] * quat[2];
        Q_type zz2 = 2.0f * quat[2] * quat[2];
        Q_type wz2 = 2.0f * quat[3] * quat[2];
        Q_type wy2 = 2.0f * quat[3] * quat[1];
        Q_type wx2 = 2.0f * quat[3] * quat[0];
        Q_type xx2 = 2.0f * quat[0] * quat[0];
        mat[0 * 4 + 0] = -yy2 - zz2 + 1.0f;
        mat[0 * 4 + 1] = xy2 + wz2;
        mat[0 * 4 + 2] = xz2 - wy2;
        mat[0 * 4 + 3] = 0;
        mat[1 * 4 + 0] = xy2 - wz2;
        mat[1 * 4 + 1] = -xx2 - zz2 + 1.0f;
        mat[1 * 4 + 2] = yz2 + wx2;
        mat[1 * 4 + 3] = 0;
        mat[2 * 4 + 0] = xz2 + wy2;
        mat[2 * 4 + 1] = yz2 - wx2;
        mat[2 * 4 + 2] = -xx2 - yy2 + 1.0f;
        mat[2 * 4 + 3] = 0;
        mat[3 * 4 + 0] = mat[3 * 4 + 1] = mat[3 * 4 + 2] = 0;
        mat[3 * 4 + 3] = 1;
    }
}

#endif // UTILS_QUATERNIONS_H
