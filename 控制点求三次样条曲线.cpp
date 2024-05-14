#include "hs_spline.h"
#include "hs_line.h"
#include <numeric>
#include <QDebug>

HS_SplineData::HS_SplineData(int degree, bool closed)
    :degree(degree), closed(closed){}

void HS_Spline::addControlPoint(const std::vector<QPointF> &points)
{
    for (auto const& point: points) {
        data.controlPoints.push_back(HS_Point(point));
    }
}

void HS_Spline::update()
{
    clearChilds();
    if(data.degree < 1 || data.degree > 3)
    {
        qDebug() << "invalid degree";
        return;
    }
    if(data.controlPoints.size() < size_t(data.degree) + 1)
    {
        qDebug() << "not enough control points";
        return;
    }
    resetBox();

    std::vector<HS_Point> tControlPoints = data.controlPoints;
    if (data.closed && (data.degree == 2 || !hasWrappedControlPoints())) {
        std::vector<HS_Point> wrappedPoints{data.controlPoints.cbegin(), data.controlPoints.cbegin() + data.degree};
        tControlPoints.insert(tControlPoints.end(), wrappedPoints.cbegin(), wrappedPoints.cend());
    }
    // 控制点数量
    const size_t npts = tControlPoints.size();
    // order:
    const size_t  k = data.degree + 1;
    // resolution: number of points on the curve
	const size_t  p1 = 8 * npts;
    std::vector<double> h(npts+1, 1.);
	std::vector<HS_Point> p(p1, {0., 0.});
    if(data.closed)
    {
        rbsplinu(npts,k,p1,tControlPoints,h,p);
    }
    else
    {
        rbspline(npts,k,p1,tControlPoints,h,p);
    }

    // child line
    HS_Point prev{false};
	for (auto const& vp: p) {
		if (prev.valid) {
            auto t = std::make_shared<HS_Line>(prev, vp);
            t->setEntityAttributes(this->attribute);
            childs.push_back(t);
		}
		prev = vp;
		minV = HS_Point::minimum(prev, minV);
		maxV = HS_Point::maximum(prev, maxV);
	}
}

/**
 * 计算样条曲线的非有理基函数
 * c: 阶数 + 1
 * t: 要计算的点位置
 * npts: 控制点数量
 * x: 控制点的位置 (节点矢量)
 * h: 控制点的权重
 */
std::vector<double> HS_Spline::rbasis(int c, double t, int npts, const std::vector<double> &x, const std::vector<double> &h) const
{
    // 控制点数量 + 阶数
    int const nplusc = npts + c;

    // 临时向量temp，并将所有元素初始化为0。 
	std::vector<double> temp(nplusc,0.);

    // 一阶非有理基函数n[i],如果t在x[i]和x[i+1]之间，temp[i]被设置为1
	for (int i = 0; i< nplusc-1; i++)
		if ((t >= x[i]) && (t < x[i+1])) temp[i] = 1;

    // 计算高阶的非有理基函数
	for (int k = 2; k <= c; k++) {
		for (int i = 0; i < nplusc-k; i++) {
			//  temp[i] = 0，也就是低阶基函数为0，跳过计算
            //  现在CAD一直样条曲线重推，不然需要判断(x[i+k-1]-x[i])分母为0
            if (temp[i] != 0)
				temp[i] = ((t-x[i])*temp[i])/(x[i+k-1]-x[i]);
            // temp[i+1] = 0，也就是低阶基函数为0，跳过计算
            // 现在CAD一直样条曲线重推，不然需要判断(x[i+k]-x[i+1])分母为0
            if (temp[i+1] != 0)
				temp[i] += ((x[i+k]-t)*temp[i+1])/(x[i+k]-x[i+1]);
        }
    }

    // t是否大于等于最后一个点x[nplusc-1]，如果是，则将temp[npts-1]设置为1。 
	if (t >= x[nplusc-1]) temp[npts-1] = 1;

    // 计算有理基函数的分母的和sum
	double sum = 0.;
	for (int i = 0; i < npts; i++) {
		sum += temp[i]*h[i];
    }

	std::vector<double> r(npts, 0);
    // 对于每个i从0到npts的值，计算r[i]的值。 
	if (sum != 0) {
		for (int i = 0; i < npts; i++)
			r[i] = (temp[i]*h[i])/sum;
	}
	return r;
}

std::vector<double> HS_Spline::knot(size_t num, size_t order) const
{
   if (data.knotslist.size() == num + order) {
		//use custom knot vector
		return data.knotslist;
	}

	std::vector<double> knotVector(num + order, 0.);
	//use uniform knots
	std::iota(knotVector.begin() + order, knotVector.begin() + num + 1, 1);
	std::fill(knotVector.begin() + num + 1, knotVector.end(), knotVector[num]);
	return knotVector;
}

void HS_Spline::rbspline(size_t npts, size_t k, size_t p1, const std::vector<HS_Point> &b, const std::vector<double> &h, std::vector<HS_Point> &p) const
{
    size_t const nplusc = npts + k;

    //  knot vector 控制点数量 + 阶数 + 1
	auto const x = knot(npts, k);

    // calculate the points on the rational B-spline curve
    double t {x[0]};
    double const step {(x[nplusc-1] - t) / (p1-1)};

	for (auto& vp: p) {
		if (x[nplusc-1] - t < 5e-6) t = x[nplusc-1];

        // 生成t的基函数
		auto const nbasis = rbasis(k, t, npts, x, h);

        // generate a point on the curve
        auto vp_p = vp.getPointF();
		for (size_t i = 0; i < npts; i++)
        {
            vp_p += b[i].getPointF() * nbasis[i];
        }
        vp.setPointF(vp_p);

		t += step;
    }
}
std::vector<double> HS_Spline::knotu(size_t num, size_t order) const
{
    if (data.knotslist.size() == num + order) {
		//use custom knot vector
		return data.knotslist;
	}
	std::vector<double> knotVector(num + order, 0.);
	std::iota(knotVector.begin(), knotVector.end(), 0);
	return knotVector;
}
void HS_Spline::rbsplinu(size_t npts, size_t k, size_t p1, const std::vector<HS_Point> &b, const std::vector<double> &h, std::vector<HS_Point> &p) const
{
    size_t const nplusc = npts + k;

	/* generate the periodic knot vector */
	std::vector<double> const x = knotu(npts, k);

    /*    calculate the points on the rational B-spline curve */
	double t = k-1;
	double const step = double(npts - k + 1)/(p1 - 1);

	for (auto& vp: p) {
		if (x[nplusc-1] - t < 5e-6) t = x[nplusc-1];

		/* generate the basis function for this value of t */
		auto const nbasis = rbasis(k, t, npts, x, h);
		/* generate a point on the curve, for x, y, z */
        auto vp_p = vp.getPointF();
		for (size_t i = 0; i < npts; i++)
        {
            vp_p += b[i].getPointF() * nbasis[i];
        }
        vp.setPointF(vp_p);
		t += step;
    }
}
bool HS_Spline::hasWrappedControlPoints() const
{
    std::vector<QPointF> controlPoints;
    for (const auto& point: data.controlPoints) {
        controlPoints.push_back(point.getPointF()); 
    }
    if (!data.closed || data.degree < 3 || controlPoints.size() < size_t(2 * data.degree) + 1)
        return false;

    return std::equal(controlPoints.cbegin(), controlPoints.cbegin() + data.degree,
               controlPoints.cbegin() + controlPoints.size() - data.degree);
}

void HS_Spline::resetBox()
{
    double maxd = HS_MAXDOUBLE;
    double mind = HS_MINDOUBLE;

    minV = {maxd, maxd};
    maxV = {mind, mind};
}

QPainterPath HS_Spline::getPath() 
{
   QPainterPath path;
   for (auto & child: childs) {
       // child转换为HS_Line
        auto line = static_cast<HS_Line*>(child.get());
        if(line == nullptr)
            continue;
        if (path.isEmpty()) {
            path.moveTo(line->getStart().getPointF());
        }
        path.lineTo(line->getEnd().getPointF());
    }
    return path;

}
