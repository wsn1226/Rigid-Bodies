#include <rodrigues.h>
#include <cmath>

//  R - rotation matrix
//  omega - angular velocity vector
void rodrigues(Eigen::Matrix3d &R, Eigen::Ref<const Eigen::Vector3d> omega)
{
    Eigen::Matrix3d omega_cross_mat;
    Eigen::Matrix3d I;
    I.setIdentity();
    omega_cross_mat.setZero();
    omega_cross_mat(0, 1) = -1.0 * omega.z();
    omega_cross_mat(0, 2) = omega.y();
    omega_cross_mat(1, 0) = omega.z();
    omega_cross_mat(1, 2) = -1.0 * omega.x();
    omega_cross_mat(2, 0) = -1.0 * omega.y();
    omega_cross_mat(2, 1) = omega.x();

    R = I + sin(1.0) * omega_cross_mat + (1.0 - cos(1.0)) * omega_cross_mat * omega_cross_mat;
}