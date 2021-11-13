#ifndef UTILS_CAMERA_H
#define UTILS_CAMERA_H

#include <string>
#include <cstdio>
#include <cstring>
#include <Eigen/Eigen>

#include "utils/constants.h"
#include "utils/maths.h"
#include "utils/quaternions.h"


namespace utils
{
    template<
            typename DerivedEye,
            typename DerivedCenter,
            typename DerivedUp,
            typename DerivedR>
    inline void look_at(
            const Eigen::PlainObjectBase<DerivedEye>& eye,
            const Eigen::PlainObjectBase<DerivedCenter>& center,
            const Eigen::PlainObjectBase<DerivedUp>& up,
            Eigen::PlainObjectBase<DerivedR>& R)
    {
        typedef Eigen::Matrix<typename DerivedR::Scalar, 3, 1> Vector3S;
        Vector3S f = (center - eye).normalized();
        Vector3S s = f.cross(up).normalized();
        Vector3S u = s.cross(f);
        R = Eigen::Matrix<typename DerivedR::Scalar, 4, 4>::Identity();
        R(0, 0) = s(0);
        R(0, 1) = s(1);
        R(0, 2) = s(2);
        R(1, 0) = u(0);
        R(1, 1) = u(1);
        R(1, 2) = u(2);
        R(2, 0) = -f(0);
        R(2, 1) = -f(1);
        R(2, 2) = -f(2);
        R(0, 3) = -s.transpose() * eye;
        R(1, 3) = -u.transpose() * eye;
        R(2, 3) = f.transpose() * eye;
    }

    template<typename DerivedP>
    inline void ortho(
            const typename DerivedP::Scalar left,
            const typename DerivedP::Scalar right,
            const typename DerivedP::Scalar bottom,
            const typename DerivedP::Scalar top,
            const typename DerivedP::Scalar nearVal,
            const typename DerivedP::Scalar farVal,
            Eigen::PlainObjectBase<DerivedP>& P)
    {
        P.setIdentity();
        P(0, 0) = 2. / (right - left);
        P(1, 1) = 2. / (top - bottom);
        P(2, 2) = -2. / (farVal - nearVal);
        P(0, 3) = -(right + left) / (right - left);
        P(1, 3) = -(top + bottom) / (top - bottom);
        P(2, 3) = -(farVal + nearVal) / (farVal - nearVal);
    }

    template<typename DerivedP>
    inline void frustum(
            const typename DerivedP::Scalar left,
            const typename DerivedP::Scalar right,
            const typename DerivedP::Scalar bottom,
            const typename DerivedP::Scalar top,
            const typename DerivedP::Scalar nearVal,
            const typename DerivedP::Scalar farVal,
            Eigen::PlainObjectBase<DerivedP>& P)
    {
        P.setConstant(4, 4, 0.);
        P(0, 0) = (2.0 * nearVal) / (right - left);
        P(1, 1) = (2.0 * nearVal) / (top - bottom);
        P(0, 2) = (right + left) / (right - left);
        P(1, 2) = (top + bottom) / (top - bottom);
        P(2, 2) = -(farVal + nearVal) / (farVal - nearVal);
        P(3, 2) = -1.0;
        P(2, 3) = -(2.0 * farVal * nearVal) / (farVal - nearVal);
    }

    template<typename Scalar>
    inline void snap_to_fixed_up(
            const Eigen::Quaternion<Scalar>& q,
            Eigen::Quaternion<Scalar>& s)
    {
        typedef Eigen::Matrix<Scalar, 3, 1> Vector3Q;
        const Vector3Q up = q.matrix() * Vector3Q(0, 1, 0);
        Vector3Q proj_up(0, up(1), up(2));
        if (proj_up.norm() == 0)
        {
            proj_up = Vector3Q(0, 1, 0);
        }
        proj_up.normalize();
        Eigen::Quaternion<Scalar> dq;
        dq = Eigen::Quaternion<Scalar>::FromTwoVectors(up, proj_up);
        s = dq * q;
    }

    template<typename Scalar>
    inline Eigen::Matrix<Scalar, 3, 1> project(
            const Eigen::Matrix<Scalar, 3, 1>& obj,
            const Eigen::Matrix<Scalar, 4, 4>& model,
            const Eigen::Matrix<Scalar, 4, 4>& proj,
            const Eigen::Matrix<Scalar, 4, 1>& viewport)
    {
        Eigen::Matrix<Scalar, 4, 1> tmp;
        tmp << obj, 1;

        tmp = model * tmp;

        tmp = proj * tmp;

        tmp = tmp.array() / tmp(3);
        tmp = tmp.array() * 0.5f + 0.5f;
        tmp(0) = tmp(0) * viewport(2) + viewport(0);
        tmp(1) = tmp(1) * viewport(3) + viewport(1);

        return tmp.head(3);
    }

    template<typename DerivedV, typename DerivedM, typename DerivedN, typename DerivedO, typename DerivedP>
    inline void project(
            const Eigen::MatrixBase<DerivedV>& V,
            const Eigen::MatrixBase<DerivedM>& model,
            const Eigen::MatrixBase<DerivedN>& proj,
            const Eigen::MatrixBase<DerivedO>& viewport,
            Eigen::PlainObjectBase<DerivedP>& P)
    {
        typedef typename DerivedP::Scalar PScalar;
        Eigen::Matrix<PScalar, DerivedV::RowsAtCompileTime, 4> HV(V.rows(), 4);
        HV.leftCols(3) = V.template cast<PScalar>();
        HV.col(3).setConstant(1);
        HV = (HV * model.template cast<PScalar>().transpose() *
              proj.template cast<PScalar>().transpose()).eval();
        HV = (HV.array().colwise() / HV.col(3).array()).eval();
        HV = (HV.array() * 0.5 + 0.5).eval();
        HV.col(0) = (HV.array().col(0) * viewport(2) + viewport(0)).eval();
        HV.col(1) = (HV.array().col(1) * viewport(3) + viewport(1)).eval();
        P = HV.leftCols(3);
    }

    template<
            typename DerivedWin,
            typename DerivedModel,
            typename DerivedProj,
            typename DerivedViewport,
            typename DerivedScene>
    inline void unproject(
            const Eigen::MatrixBase<DerivedWin>& win,
            const Eigen::MatrixBase<DerivedModel>& model,
            const Eigen::MatrixBase<DerivedProj>& proj,
            const Eigen::MatrixBase<DerivedViewport>& viewport,
            Eigen::PlainObjectBase<DerivedScene>& scene)
    {
        if (win.cols() != 3)
        {
            assert(win.rows() == 3);
            // needless transposes
            Eigen::Matrix<typename DerivedScene::Scalar, 1, 3> sceneT;
            unproject(win.transpose().eval(), model, proj, viewport, sceneT);
            scene = sceneT.head(3);
            return;
        }
        assert(win.cols() == 3);
        const int n = win.rows();
        scene.resize(n, 3);
        for (int i = 0; i < n; i++)
        {
            typedef typename DerivedScene::Scalar Scalar;
            Eigen::Matrix<Scalar, 4, 4> Inverse =
                    (proj.template cast<Scalar>() * model.template cast<Scalar>()).inverse();

            Eigen::Matrix<Scalar, 4, 1> tmp;
            tmp << win.row(i).head(3).transpose(), 1;
            tmp(0) = (tmp(0) - viewport(0, 0)) / viewport(2, 0);
            tmp(1) = (tmp(1) - viewport(1, 0)) / viewport(3, 0);
            tmp = tmp.array() * 2.0f - 1.0f;

            Eigen::Matrix<Scalar, 4, 1> obj = Inverse * tmp;
            obj /= obj(3);

            scene.row(i).head(3) = obj.head(3);
        }
    }

    template<typename Scalar>
    inline Eigen::Matrix<Scalar, 3, 1> unproject(
            const Eigen::Matrix<Scalar, 3, 1>& win,
            const Eigen::Matrix<Scalar, 4, 4>& model,
            const Eigen::Matrix<Scalar, 4, 4>& proj,
            const Eigen::Matrix<Scalar, 4, 1>& viewport)
    {
        Eigen::Matrix<Scalar, 3, 1> scene;
        unproject(win, model, proj, viewport, scene);
        return scene;
    }

    template<typename ScalarQuatDown, typename ScalarQuat>
    inline void two_axis_evaluator_fixed_up(
            const double w,
            const double h,
            const double speed,
            const Eigen::Quaternion<ScalarQuatDown>& down_quat,
            const int down_x,
            const int down_y,
            const int mouse_x,
            const int mouse_y,
            Eigen::Quaternion<ScalarQuat>& quat)
    {
        {
            Eigen::Matrix<ScalarQuat, 3, 1> axis(0, 1, 0);
            quat = down_quat *
                   Eigen::Quaternion<ScalarQuat>(
                           Eigen::AngleAxis<ScalarQuat>(
                                   PI * ((ScalarQuat) (mouse_x - down_x)) / (ScalarQuat) w * speed / 2.0,
                                   axis.normalized()));
            quat.normalize();
        }
        {
            Eigen::Matrix<ScalarQuat, 3, 1> axis(1, 0, 0);
            if (axis.norm() != 0)
            {
                quat =
                        Eigen::Quaternion<ScalarQuat>(
                                Eigen::AngleAxis<ScalarQuat>(
                                        PI * (mouse_y - down_y) / (ScalarQuat) h * speed / 2.0,
                                        axis.normalized())) * quat;
                quat.normalize();
            }
        }
    }

    // Note: For the canonical view quaternions it should be completely possible to
    // determine this analytically. That is the max_distance should be a
    // theoretical known value
    // Also: I'm not sure it matters in this case, but. We are dealing with
    // quaternions on the 4d unit sphere, but measuring distance in general 4d
    // space (i.e. not geodesics on the sphere). Probably something with angles
    // would be better.
    template<typename Q_type>
    inline bool snap_to_canonical_view_quat(
            const Q_type* q,
            const Q_type threshold,
            Q_type* s)
    {
        // Copy input into output
        // CANNOT use std::copy here according to:
        // http://www.cplusplus.com/reference/algorithm/copy/
        s[0] = q[0];
        s[1] = q[1];
        s[2] = q[2];
        s[3] = q[3];

        // Normalize input quaternion
        Q_type qn[4];
        bool valid_len = quaternions::normalize_quat(q, qn);
        // If normalizing valid then don't bother
        if (!valid_len)
        {
            return false;
        }

        // 0.290019
        const Q_type MAX_DISTANCE = 0.4;
        Q_type min_distance = 2 * MAX_DISTANCE;
        int min_index = -1;
        double min_sign = 0;
        // loop over canonical view quaternions
        for (int sign = -1; sign <= 1; sign += 2)
        {
            for (int i = 0; i < NUM_CANONICAL_VIEW_QUAT; i++)
            {
                Q_type distance = 0.0;
                // loop over coordinates
                for (int j = 0; j < 4; j++)
                {
                    // Double cast because of bug in llvm version 4.2 with -O3
                    distance +=
                            (qn[j] - (double) sign * quaternions::CANONICAL_VIEW_QUAT<Q_type>(i, j)) *
                            (qn[j] - (double) sign * quaternions::CANONICAL_VIEW_QUAT<Q_type>(i, j));
                }
                if (min_distance > distance)
                {
                    min_distance = distance;
                    min_index = i;
                    min_sign = (double) sign;
                }
            }
        }

        if (MAX_DISTANCE < min_distance)
        {
            fprintf(
                    stderr,
                    "ERROR: found new max MIN_DISTANCE: %g\n"
                    "PLEASE update snap_to_canonical_quat()",
                    min_distance);
        }

        assert(min_distance < MAX_DISTANCE);
        assert(min_index >= 0);

        if (min_distance / MAX_DISTANCE <= threshold)
        {
            // loop over coordinates
            for (int j = 0; j < 4; j++)
            {
                s[j] = min_sign * quaternions::CANONICAL_VIEW_QUAT<Q_type>(min_index, j);
            }
            return true;
        }
        return false;
    }

    template<typename ScalarQ, typename ScalarS>
    inline bool snap_to_canonical_view_quat(
            const Eigen::Quaternion<ScalarQ>& q,
            const double threshold,
            Eigen::Quaternion<ScalarS>& s)
    {
        return snap_to_canonical_view_quat<ScalarS>(
                q.coeffs().data(), threshold, s.coeffs().data());
    }

    // Utility functions
    template<typename Q_type>
    static inline Q_type QuatD(double w, double h)
    {
        using namespace std;
        return (Q_type) (std::abs(w) < std::abs(h) ? std::abs(w) : std::abs(h)) - 4;
    }

    template<typename Q_type>
    static inline Q_type QuatIX(double x, double w, double h)
    {
        return (2.0f * (Q_type)
                x - (Q_type)
                        w - 1.0f) / QuatD<Q_type>(w, h);
    }

    template<typename Q_type>
    static inline Q_type QuatIY(double y, double w, double h)
    {
        return (-2.0f * (Q_type)
                y + (Q_type)
                        h - 1.0f) / QuatD<Q_type>(w, h);
    }

    // This is largely the trackball as implemented in AntTweakbar. Much of the
    // code is straight from its source in TwMgr.cpp
    // http://www.antisphere.com/Wiki/tools:anttweakbar
    template<typename Q_type>
    inline void trackball(
            const double w,
            const double h,
            const Q_type speed_factor,
            const double down_mouse_x,
            const double down_mouse_y,
            const double mouse_x,
            const double mouse_y,
            Q_type* quat)
    {
        assert(speed_factor > 0);

        double original_x =
                QuatIX<Q_type>(speed_factor * (down_mouse_x - w / 2) + w / 2, w, h);
        double original_y =
                QuatIY<Q_type>(speed_factor * (down_mouse_y - h / 2) + h / 2, w, h);

        double x = QuatIX<Q_type>(speed_factor * (mouse_x - w / 2) + w / 2, w, h);
        double y = QuatIY<Q_type>(speed_factor * (mouse_y - h / 2) + h / 2, w, h);

        double z = 1;
        double n0 = sqrt(original_x * original_x + original_y * original_y + z * z);
        double n1 = sqrt(x * x + y * y + z * z);
        if (n0 > utils::DOUBLE_EPS && n1 > utils::DOUBLE_EPS)
        {
            double v0[] = {original_x / n0, original_y / n0, z / n0};
            double v1[] = {x / n1, y / n1, z / n1};
            double axis[3];
            maths::cross(v0, v1, axis);
            double sa = sqrt(maths::dot(axis, axis));
            double ca = maths::dot(v0, v1);
            double angle = atan2(sa, ca);
            if (x * x + y * y > 1.0)
            {
                angle *= 1.0 + 0.2f * (sqrt(x * x + y * y) - 1.0);
            }
            double q_rot[4];
            quaternions::axis_angle_to_quat(axis, angle, q_rot);
            quat[0] = q_rot[0];
            quat[1] = q_rot[1];
            quat[2] = q_rot[2];
            quat[3] = q_rot[3];
        }
    }

    template<typename Q_type>
    inline void trackball(
            const double w,
            const double h,
            const Q_type speed_factor,
            const Q_type* down_quat,
            const double down_mouse_x,
            const double down_mouse_y,
            const double mouse_x,
            const double mouse_y,
            Q_type* quat)
    {
        double q_rot[4], q_res[4], q_orig[4];
        trackball<double>(
                w, h,
                speed_factor,
                down_mouse_x, down_mouse_y,
                mouse_x, mouse_y,
                q_rot);
        double n_q_orig =
                sqrt(down_quat[0] * down_quat[0] +
                     down_quat[1] * down_quat[1] +
                     down_quat[2] * down_quat[2] +
                     down_quat[3] * down_quat[3]);

        if (fabs(n_q_orig) > utils::DOUBLE_EPS_SQ)
        {
            q_orig[0] = down_quat[0] / n_q_orig;
            q_orig[1] = down_quat[1] / n_q_orig;
            q_orig[2] = down_quat[2] / n_q_orig;
            q_orig[3] = down_quat[3] / n_q_orig;
            quaternions::quat_mult<double>(q_rot, q_orig, q_res);
            quat[0] = q_res[0];
            quat[1] = q_res[1];
            quat[2] = q_res[2];
            quat[3] = q_res[3];
        }
        else
        {
            quat[0] = q_rot[0];
            quat[1] = q_rot[1];
            quat[2] = q_rot[2];
            quat[3] = q_rot[3];
        }
    }

    template<typename ScalarQuatDown, typename ScalarQuat>
    inline void trackball(
            const double w,
            const double h,
            const double speed_factor,
            const Eigen::Quaternion<ScalarQuatDown>& down_quat,
            const double down_mouse_x,
            const double down_mouse_y,
            const double mouse_x,
            const double mouse_y,
            Eigen::Quaternion<ScalarQuat>& quat)
    {
        return trackball(
                w,
                h,
                (ScalarQuat) speed_factor,
                down_quat.coeffs().data(),
                down_mouse_x,
                down_mouse_y,
                mouse_x,
                mouse_y,
                quat.coeffs().data());
    }
}

#endif // UTILS_CAMERA_H
