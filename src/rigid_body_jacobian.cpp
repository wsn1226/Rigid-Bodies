#include <rigid_body_jacobian.h>

// Input:
//   R - rotation matrix for rigid body
//   p - world space position of center-of-mass
//   X -  undeformed position at which to compute the Jacobian.
// Output:
//   J - the rigid body jacobian acting on the undeformed space point X.
void rigid_body_jacobian(Eigen::Matrix36d &J,
                         Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p,
                         Eigen::Ref<const Eigen::Vector3d> x)
{
    Eigen::Matrix3d x_cross_mat;
    Eigen::Matrix3d I;
    I.setIdentity();
    x_cross_mat.setZero();
    x_cross_mat(0, 1) = -1.0 * x.z();
    x_cross_mat(0, 2) = x.y();
    x_cross_mat(1, 0) = x.z();
    x_cross_mat(1, 2) = -1.0 * x.x();
    x_cross_mat(2, 0) = -1.0 * x.y();
    x_cross_mat(2, 1) = x.x();

    J.block(0, 0, 3, 3) = R * x_cross_mat.transpose() * R.transpose();
    J.block(3, 0, 3, 3) = I;
}
