#include "toolfunc.h"
#include <cmath>
#include <QVector2D>
#include <numeric>  

HS_Math::MatrixXd HS_Math::BSplineFitting(const std::vector<QPointF> &points, int degree)
{
    // 型值点数量
    int n = points.size();
    // 控制点数量
    int num_control_points = n + degree - 1;
    // 节点矢量数量
    int num_knot = n + degree + 3;
    // 节点矢量u，并将所有元素初始化为0。
    std::vector<double> u(num_knot,0.);

    // 积累弦长法进行型值点的参数化 n-1个弦长
    std::vector<double> t(n-1, 0);

    for (int i = 0; i < n-1; i++) {
        t[i] = std::sqrt(std::pow(points[i+1].x() - points[i].x(), 2) + std::pow(points[i+1].y() - points[i].y(), 2));
    }

    // 计算弦长和
    double sum = std::accumulate(t.begin(), t.end(), 0.0);

    // 计算节点矢量u
    // 前degree + 1个节点矢量u为0
    for (int i = 0; i < degree + 1; i++) {
        u[i] = 0;
    }
    // 后degree + 1个节点矢量u为1
    for (int i = num_knot - degree - 1; i < num_knot; i++) {
        u[i] = 1;
    }
    // n个数据点，内节点为n-2个,u[degree+1]作为初始值
    for (int i = degree + 1; i < num_knot - degree - 1; i++) {
        u[i] = u[i-1] + t[i-degree-1] / sum;
    }

    // 通过points构造pt
    Eigen::MatrixXd pt(n, 2);
    for (int i = 0; i < n; i++) {
        pt(i, 0) = points[i].x();
        pt(i, 1) = points[i].y();
    }

    // 反算n + degree - 1个控制点 (切矢边界)
    // 首数据点切矢
    Eigen::RowVector2d dpt1(0, 1);

    // 末数据点切矢
    Eigen::RowVector2d dptn(-1, 0);
    std::vector<double> dU(num_knot, 0.);
    for (int i = degree; i < n + degree -1; i++) {
        dU[i] = u[i+1] - u[i];
    }
    //求解线性方程组获得控制顶点向量，A*D=E,A为系数矩阵，元素为B样条基函数的值；
    //D是控制顶点列向量；E是列向量

    Eigen::MatrixXd A(n, n);
    Eigen::MatrixXd E(n, 2);
    // 切矢条件 a1 = 1, b1 = c1 = 0
    A(0,0) = 1;

    // 切矢条件 an = bn = 0, cn = 1
    A(n - 1, n - 1) = 1;

    // 首端点条件
    E.row(0) = pt.row(0) + (dU[3] / 3) * dpt1;

    // 末端点条件
    E.row(n - 1) = pt.row(n - 1) - (dU[n + 1] / 3) * dptn;

    for (int i = 1; i < n - 1; i++) {
        double a = dU[i + 3] * dU[i + 3] / (dU[i + 1] + dU[i + 2] + dU[i + 3]);
        double b = dU[i + 3] * (dU[i + 1] + dU[i + 2]) / (dU[i + 1] + dU[i + 2] + dU[i + 3])
                 + dU[i + 2] * (dU[i + 3] + dU[i + 4]) / (dU[i + 2] + dU[i + 3] + dU[i + 4]);
        double c = dU[i + 2] * dU[i + 2] / (dU[i + 2] + dU[i + 3] + dU[i + 4]);

        A(i, i - 1) = a;
        A(i, i) = b;
        A(i, i + 1) = c;

        E.row(i) = (dU[i + 2] + dU[i + 3]) * pt.row(i);
    }
    // 去除首末端点的控制顶点向量
    Eigen::VectorXd D = A.fullPivLu().solve(E.col(0)); // 解线性方程组
    Eigen::VectorXd M = A.fullPivLu().solve(E.col(1)); // 解线性方程组

    // 由于数据点比控制顶点少两个，添加首末端点
    Eigen::MatrixXd controlPoints(num_control_points , 2);
    controlPoints.row(0) = pt.row(0);
    for (int i = 1; i < num_control_points -1; i++) {
        controlPoints(i, 0) = D(i-1);  // 设置 x 坐标
        controlPoints(i, 1) = M(i-1);  // 设置 y 坐标
    }
    controlPoints.row(num_control_points - 1) = pt.row(n - 1);

    // 控制点数量
    int rrr = controlPoints.rows();
    int ccc = controlPoints.cols();

    return controlPoints;
}

std::vector<QPointF> HS_Math::EigenMatrixToVector(const Eigen::MatrixXd &matrix)
{

    int rows = matrix.rows();
    int cols = matrix.cols();
    std::vector<QPointF> points;

    if (cols == 2) { // 确保矩阵的列数是2，以确保可以转换为Point
        for (int i = 0; i < rows; i++) {
            QPointF point(matrix(i, 0), matrix(i, 1));
            points.push_back(point);
        }
    }
    return points;

}
