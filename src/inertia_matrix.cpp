#include <inertia_matrix.h>
#include <cassert>
#include <iostream>
// Input:
//   V - the nx3 matrix of vertices.
//   F - the mx3 matrix of triangle vertex indices.
//   density - the material density.
// Output:
//   I - the 3x3 angular inertia matrix
//   center - the center of mass of the object
//   mass - the total mass of the object
// compute inertia matrix and volume by integrating on surfaces
void inertia_matrix(Eigen::Matrix3d &I, Eigen::Vector3d &center, double &mass, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double density)
{
    mass = 0;
    center.setZero();
    I.setZero();
    double area, x0, x1, x2, y0, y1, y2, z0, z1, z2;
    int ind0, ind1, ind2;
    Eigen::Vector3d AB, AC, n;
    for (int i = 0; i < F.rows(); i++)
    {
        ind0 = F(i, 0);
        ind1 = F(i, 1);
        ind2 = F(i, 2);
        x0 = V(ind0, 0);
        x1 = V(ind1, 0);
        x2 = V(ind2, 0);

        y0 = V(ind0, 1);
        y1 = V(ind1, 1);
        y2 = V(ind2, 1);

        z0 = V(ind0, 2);
        z1 = V(ind1, 2);
        z2 = V(ind2, 2);

        AB = V.row(ind1) - V.row(ind0);
        AC = V.row(ind2) - V.row(ind0);
        //(yAB⋅zAC−zAB⋅yAC)2+(zAB⋅xAC−xAB⋅zAC)2+(xAB⋅yAC−yAB⋅xAC)2
        n = AB.cross(AC);
        area = n.norm() / 2.0;
        n.normalize();
        mass += density * area * n.x() * 1 / 3.0 * (x0 + y0 + z0);
        center(0) += density * area * n.x() * (x0 * x0 / 12 + (x0 * x1) / 12 + (x0 * x2) / 12 + x1 * x1 / 12 + (x1 * x2) / 12 + x2 * x2 / 12);
        center(1) += density * area * n.y() * (y0 * y0 / 12 + (y0 * y1) / 12 + (y0 * y2) / 12 + y1 * y1 / 12 + (y1 * y2) / 12 + y2 * y2 / 12);
        center(2) += density * area * n.z() * (z0 * z0 / 12 + (z0 * z1) / 12 + (z0 * z2) / 12 + z1 * z1 / 12 + (z1 * z2) / 12 + z2 * z2 / 12);
    }
    center /= mass;
    Eigen::MatrixXd V_C;
    V_C.resize(V.rows(), V.cols());
    for (int i = 0; i < V.rows(); i++)
    {
        V_C.row(i) = V.row(i) - center.transpose();
    }

    double xx, yy, zz, xy, xz, yz;
    for (int i = 0; i < F.rows(); i++)
    {
        ind0 = F(i, 0);
        ind1 = F(i, 1);
        ind2 = F(i, 2);
        x0 = V_C(ind0, 0);
        x1 = V_C(ind1, 0);
        x2 = V_C(ind2, 0);

        y0 = V_C(ind0, 1);
        y1 = V_C(ind1, 1);
        y2 = V_C(ind2, 1);

        z0 = V_C(ind0, 2);
        z1 = V_C(ind1, 2);
        z2 = V_C(ind2, 2);

        AB = (V_C.row(ind1) - V_C.row(ind0)).transpose();
        AC = (V_C.row(ind2) - V_C.row(ind0)).transpose();
        n = AB.cross(AC);
        area = n.norm() / 2.0;

        n.normalize();

        xx = (2.0 / 3.0) * area * density * n.x() * (x0 * x0 * x0 / 20 + (x0 * x0 * x1) / 20 + (x0 * x0 * x2) / 20 + (x0 * x1 * x1) / 20 + (x0 * x1 * x2) / 20 + (x0 * x2 * x2) / 20 + x1 * x1 * x1 / 20 + (x1 * x1 * x2) / 20 + (x1 * x2 * x2) / 20 + x2 * x2 * x2 / 20);
        yy = (2.0 / 3.0) * area * density * n.y() * (y0 * y0 * y0 / 20 + (y0 * y0 * y1) / 20 + (y0 * y0 * y2) / 20 + (y0 * y1 * y1) / 20 + (y0 * y1 * y2) / 20 + (y0 * y2 * y2) / 20 + y1 * y1 * y1 / 20 + (y1 * y1 * y2) / 20 + (y1 * y2 * y2) / 20 + y2 * y2 * y2 / 20);
        zz = (2.0 / 3.0) * area * density * n.z() * (z0 * z0 * z0 / 20 + (z0 * z0 * z1) / 20 + (z0 * z0 * z2) / 20 + (z0 * z1 * z1) / 20 + (z0 * z1 * z2) / 20 + (z0 * z2 * z2) / 20 + z1 * z1 * z1 / 20 + (z1 * z1 * z2) / 20 + (z1 * z2 * z2) / 20 + z2 * z2 * z2 / 20);
        xy = area * density * n.x() * (x0 * x0 * y0) / 20 + (x0 * x0 * y1) / 60 + (x1 * x1 * y0) / 60 + (x0 * x0 * y2) / 60 + (x1 * x1 * y1) / 20 + (x2 * x2 * y0) / 60 + (x1 * x1 * y2) / 60 + (x2 * x2 * y1) / 60 + (x2 * x2 * y2) / 20 + (x0 * x1 * y0) / 30 + (x0 * x1 * y1) / 30 + (x0 * x2 * y0) / 30 + (x0 * x1 * y2) / 60 + (x0 * x2 * y1) / 60 + (x1 * x2 * y0) / 60 + (x0 * x2 * y2) / 30 + (x1 * x2 * y1) / 30 + (x1 * x2 * y2) / 30;
        xz = area * density * n.x() * (x0 * x0 * z0) / 20 + (x0 * x0 * z1) / 60 + (x1 * x1 * z0) / 60 + (x0 * x0 * z2) / 60 + (x1 * x1 * z1) / 20 + (x2 * x2 * z0) / 60 + (x1 * x1 * z2) / 60 + (x2 * x2 * z1) / 60 + (x2 * x2 * z2) / 20 + (x0 * x1 * z0) / 30 + (x0 * x1 * z1) / 30 + (x0 * x2 * z0) / 30 + (x0 * x1 * z2) / 60 + (x0 * x2 * z1) / 60 + (x1 * x2 * z0) / 60 + (x0 * x2 * z2) / 30 + (x1 * x2 * z1) / 30 + (x1 * x2 * z2) / 30;
        yz = area * density * n.y() * (y0 * y0 * z0) / 20 + (y0 * y0 * z1) / 60 + (y1 * y1 * z0) / 60 + (y0 * y0 * z2) / 60 + (y1 * y1 * z1) / 20 + (y2 * y2 * z0) / 60 + (y1 * y1 * z2) / 60 + (y2 * y2 * z1) / 60 + (y2 * y2 * z2) / 20 + (y0 * y1 * z0) / 30 + (y0 * y1 * z1) / 30 + (y0 * y2 * z0) / 30 + (y0 * y1 * z2) / 60 + (y0 * y2 * z1) / 60 + (y1 * y2 * z0) / 60 + (y0 * y2 * z2) / 30 + (y1 * y2 * z1) / 30 + (y1 * y2 * z2) / 30;

        I(0, 0) += (yy + zz);
        I(0, 1) += (-1.0 * xy);
        I(0, 2) += (-1.0 * xz);
        I(1, 0) += (-1.0 * xy);
        I(1, 1) += (xx + zz);
        I(1, 2) += (-1.0 * yz);
        I(2, 0) += (-1.0 * xz);
        I(2, 1) += (-1.0 * yz);
        I(2, 2) += (xx + yy);
    }
}