/**
 * @file WCSTNX.h 声明文件, 基于非标准WCS格式TNX, 计算图像坐标与WCS坐标之间的对应关系
 * @version 0.1
 * @date 2017年11月9日
 * - 从FITS文件头加载TNX参数项
 * - 从文本文件加载TNX参数项
 * - 将TNX参数项写入FITS文件
 * - 计算(x,y)对应的WCS坐标(ra, dec)
 * - 计算(ra,dec)对应的图像坐标(x,y)
 *
 * @version 0.2
 * @date 2020年6月30日
 * - 增加: 拟合功能. 使用匹配的XY和RA/DEC拟合模型
 */
#include <stdlib.h>
#include <string.h>
#include "ADefine.h"
#include "WCSTNX.h"

namespace AstroUtil {
//////////////////////////////////////////////////////////////////////////////
PrjTNXRes::PrjTNXRes() {
	func = xterm = 0;
	xorder = yorder = 0;
	xmin = ymin = 0.0;
	xmax = ymax = 0.0;
	basefunc = NULL;
	nitem = 0;
	coef = xv = yv = NULL;
}

PrjTNXRes::~PrjTNXRes() {
	free_array(&coef);
	free_array(&xv);
	free_array(&yv);
}

void PrjTNXRes::SetParam(int func, int xterm, int xorder, int yorder) {
	this->func   = func;
	this->xterm  = xterm;
	this->xorder = xorder;
	this->yorder = yorder;
}

void PrjTNXRes::SetRange(double xmin, double ymin, double xmax, double ymax) {
	this->xmin = xmin;
	this->ymin = ymin;
	this->xmax = xmax;
	this->ymax = ymax;
}

bool PrjTNXRes::Initialize() {
	if      (func == FUNC_LINEAR)    basefunc = &power_array;
	else if (func == FUNC_LEGENDRE)  basefunc = &legendre_array;
	else if (func == FUNC_CHEBYSHEV) basefunc = &chebyshev_array;
	else return false;
	if (xterm < X_NONE || xterm > X_FULL) return false;
	if (xorder < 2 || xorder > 10) return false;
	if (yorder < 2 || yorder > 10) return false;

	int n;
	int order = xorder < yorder ? xorder : yorder;
	if      (xterm == X_NONE) n = xorder + yorder - 1;
	else if (xterm == X_FULL) n = xorder * yorder;
	else                      n = xorder * yorder - order * (order - 1) / 2;
	if (n != nitem) UnInitialize();
	nitem = n;
	if (coef == NULL) coef = (double *) calloc(nitem, sizeof(double));
	if (xv == NULL)   xv = (double *) calloc(nitem, sizeof(double));
	if (yv == NULL)   yv = (double *) calloc(nitem, sizeof(double));

	return true;
}

void PrjTNXRes::UnInitialize() {
	free_array(&coef);
	free_array(&xv);
	free_array(&yv);
}

bool PrjTNXRes::ItemVector(double x, double y, double *vptr) {
	if (!xv) return false;
	int maxorder = xorder > yorder ? xorder : yorder;
	int i, j, k, imax(xorder);
	double t;

	(*basefunc)(x, xmin, xmax, xorder, xv);
	(*basefunc)(y, ymin, ymax, yorder, yv);

	for (j = k = 0; j < yorder; ++j) {
		if (j) {
			if (xterm == X_NONE && imax != 1) imax = 1;
			else if (xterm == X_HALF && (j + xorder) > maxorder) --imax;
		}
		for (i = 0, t = yv[j]; i < imax; ++i, ++k) vptr[k] = xv[i] * t;
	}
	return true;
}

double PrjTNXRes::PolyVal(double x, double y) {
	if (!coef) return 0.0;

	int maxorder = xorder > yorder ? xorder : yorder;
	int i, j, k, imax(xorder);
	double t, sum;

	(*basefunc)(x, xmin, xmax, xorder, xv);
	(*basefunc)(y, ymin, ymax, yorder, yv);

	for (j = k = 0, sum = 0.0; j < yorder; ++j) {
		if (j) {
			if      (xterm == X_NONE && imax != 1) imax = 1;
			else if (xterm == X_HALF && (j + xorder) > maxorder) --imax;
		}

		for (i = 0, t = 0.0; i < imax; ++i, ++k) t += (coef[k] * xv[i]);
		sum += (t * yv[j]);
	}
	return sum;
}

void PrjTNXRes::free_array(double **ptr) {
	if (*ptr) {
		free(*ptr);
		*ptr = NULL;
	}
}

void PrjTNXRes::power_array(double value, double min, double max, int order, double *ptr) {
	ptr[0] = 1.0;
	for (int i = 1; i < order; ++i) ptr[i] = value * ptr[i - 1];
}

void PrjTNXRes::legendre_array(double value, double min, double max, int order, double *ptr) {
	double norm = ((max + min) - 2.0 * value) / (max - min);

	ptr[0] = 1.0;
	ptr[1] = norm;
	for (int i = 2; i < order; ++i) {
		ptr[i] = ((2 * i - 1) * norm * ptr[i - 1] - (i - 1) * ptr[i - 2]) / i;
	}
}

void PrjTNXRes::chebyshev_array(double value, double min, double max, int order, double *ptr) {
	double norm = ((max + min) - 2.0 * value) / (max - min);

	ptr[0] = 1.0;
	ptr[1] = norm;
	for (int i = 2; i < order; ++i) ptr[i] = 2 * norm * ptr[i - 1] - ptr[i - 2];
}
//////////////////////////////////////////////////////////////////////////////
void PrjTNX::Image2WCS(double x, double y, double &ra, double &dc) {

}

void PrjTNX::WCS2Image(double ra, double dc, double &x, double &y) {

}

void PrjTNX::image_to_plane(double x, double y, double &xi, double &eta) {
	double dx(x - ref_pixx), dy(y - ref_pixy);
	xi  = cd[0][0] * dx + cd[0][1] * dy;
	eta = cd[1][0] * dx + cd[1][1] * dy;
}

//////////////////////////////////////////////////////////////////////////////
WCSTNX::WCSTNX() {

}

WCSTNX::~WCSTNX() {

}
//////////////////////////////////////////////////////////////////////////////
} /* namespace AstroUtil */
