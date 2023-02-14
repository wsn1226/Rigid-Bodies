#include <Eigen/Dense>
#include <EigenTypes.h>
#include <rodrigues.h>
#include <iostream>

// Input:
//   q - 12n vector where n is the number of rigid bodies. Each rigid body is stored as 12 doubles.
//       The first 9 doubles are the columns of a 3x3 rotation matrix and the final 3 doubles are the world space position of the object's center of mass.
//   qdot - 6n vector of generalied velocities. The first 3 doubles of each body are the world space angular velocity and
//          the second 3 are the world space linear velocity.
//   dt - the integration time step
//   masses - a vector to mass matrices for each rigid body
//   forces - a 6n vector of generalized forces for n rigid bodies. The first 3 doubles of each rigid body are the torques acting on the object
//            while the second 3 doubles are the linear forces.
// Output:
//   q - updated generalized coordinates
//   qdot - updated generalized velocities
inline void exponential_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt,
                              std::vector<Eigen::Matrix66d> &masses, Eigen::Ref<const Eigen::VectorXd> forces)
{
    Eigen::Vector3d p, omega, center, torq, fext, omega_next, pdot, pdot_next, p_next;
    Eigen::Matrix3d R, I, RIR_t, R_next, matrix_exp;
    double mass;
    for (int i = 0; i < q.size() / 12; i++)
    {
        Eigen::Matrix3d R = Eigen::Map<const Eigen::Matrix3d>(q.segment<9>(12 * i).data());
        // Eigen::Vector3d p = Eigen::Map<const Eigen::Vector3d>(q.segment<3>(12 * i).data());
        p = q.segment(12 * i + 9, 3);

        omega = qdot.segment(6 * i, 3);
        pdot = qdot.segment(6 * i + 3, 3);

        torq = forces.segment(6 * i, 3);
        fext = forces.segment(6 * i + 3, 3);
        I = masses[i].block(0, 0, 3, 3);
        mass = masses[i](3, 3);
        RIR_t = R * I * R.transpose();
        omega_next = omega + dt * RIR_t.inverse() * ((omega.cross(RIR_t * omega)) + torq);
        pdot_next = pdot + dt * (1.0 / mass) * fext;
        p_next = p + dt * pdot;

        rodrigues(matrix_exp, omega * dt);
        R_next = matrix_exp * R;

        // update
        q.segment(12 * i, 3) = R_next.block(0, 0, 3, 1);
        q.segment(12 * i + 3, 3) = R_next.block(0, 1, 3, 1);
        q.segment(12 * i + 6, 3) = R_next.block(0, 2, 3, 1);
        q.segment(12 * i + 9, 3) = p_next;
        std::cout << qdot << std::endl;
        qdot.segment(6 * i, 3) = omega_next;
        qdot.segment(6 * i + 3, 3) = pdot_next;
    }
}