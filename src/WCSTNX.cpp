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
#include <stdio.h>
#include <longnam.h>
#include <fitsio.h>
#include <string>
#include "ADefine.h"
#include "AMath.h"
#include "WCSTNX.h"

using namespace std;

namespace AstroUtil {
//////////////////////////////////////////////////////////////////////////////
/* TNX投影的畸变改正 */
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
	int n;
	int order = xorder < yorder ? xorder : yorder;

	if      (func == FUNC_LEGENDRE)  basefunc = &legendre_array;
	else if (func == FUNC_LINEAR)    basefunc = &power_array;
	else if (func == FUNC_CHEBYSHEV) basefunc = &chebyshev_array;
	else return false;
	if (xorder < 2 || xorder > 10) return false;
	if (yorder < 2 || yorder > 10) return false;
	if      (xterm == X_NONE) n = xorder + yorder - 1;
	else if (xterm == X_FULL) n = xorder * yorder;
	else if (xterm == X_HALF) n = xorder * yorder - order * (order - 1) / 2;
	else return false;
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
/* TNX投影+模型 */
void PrjTNX::Image2WCS(double x, double y, double &ra, double &dc) {
	double xi, eta;
	Image2Plane(x, y, xi, eta);
	xi  += res[0].PolyVal(x, y) * AS2R;
	eta += res[1].PolyVal(x, y) * AS2R;
	Plane2WCS(xi, eta, ra, dc);
	ra *= R2D;
	dc *= R2D;
}

void PrjTNX::WCS2Image(double ra, double dc, double &x, double &y) {
	double xi, eta, dxi, deta, xt, yt;
	int cnt(0);

	WCS2Plane(ra * D2R, dc * D2R, xi, eta);
	Plane2Image(xi, eta, x, y);

	do {
		xt = x;
		yt = y;
		dxi  = res[0].PolyVal(x, y) * AS2R;
		deta = res[1].PolyVal(x, y) * AS2R;
		Plane2Image(xi - dxi, eta - deta, x, y);
	} while ((fabs(xt - x) > 1E-3 || fabs(yt - y) > 1E-3) && ++cnt < 10);
}

void PrjTNX::Image2Plane(double x, double y, double &xi, double &eta) {
	double dx(x - ref_pixx), dy(y - ref_pixy);
	xi  = (cd[0][0] * dx + cd[0][1] * dy) * D2R;
	eta = (cd[1][0] * dx + cd[1][1] * dy) * D2R;
}

void PrjTNX::Plane2WCS(double xi, double eta, double &ra, double &dc) {
	double fract = cos(ref_wcsy) - eta * sin(ref_wcsy);
	ra = cyclemod(ref_wcsx + atan2(xi, fract), A2PI);
	dc = atan(((eta * cos(ref_wcsy) + sin(ref_wcsy)) * cos(ra - ref_wcsx)) / fract);
}

void PrjTNX::WCS2Plane(double ra, double dc, double &xi, double &eta) {
	double fract = sin(ref_wcsy) * sin(dc) + cos(ref_wcsy) * cos(dc) * cos(ra - ref_wcsx);
	xi  = cos(dc) * sin(ra - ref_wcsx) / fract;
	eta = (cos(ref_wcsy) * sin(dc) - sin(ref_wcsy) * cos(dc) * cos(ra - ref_wcsx)) / fract;
}

void PrjTNX::Plane2Image(double xi, double eta, double &x, double &y) {
	xi  *= R2D;
	eta *= R2D;
	x = (ccd[0][0] * xi + ccd[0][1] * eta) + ref_pixx;
	y = (ccd[1][0] * xi + ccd[1][1] * eta) + ref_pixy;
}

//////////////////////////////////////////////////////////////////////////////
/* TNX WCS接口 */
WCSTNX::WCSTNX() {

}

WCSTNX::~WCSTNX() {

}

const PrjTNX *WCSTNX::GetModel() {
	return &model_;
}

void WCSTNX::SetNormalRange(double xmin, double ymin, double xmax, double ymax) {
	model_.res[0].SetRange(xmin, ymin, xmax, ymax);
	model_.res[1].SetRange(xmin, ymin, xmax, ymax);
}

void WCSTNX::SetParamRes(int func, int xterm, int xorder, int yorder) {
	model_.res[0].SetParam(func, xterm, xorder, yorder);
	model_.res[1].SetParam(func, xterm, xorder, yorder);
}

int WCSTNX::LoadImage(const char* filepath) {
	fitsfile *fitsptr;
	int status(0);
	fits_open_file(&fitsptr, filepath, 0, &status);
	if (status) return -1; // FITS文件打开失败

	char keyword[10], value[70], CTYPE1[10], CTYPE2[10];
	char strcor[2][2048], (*ptr)[2048];
	int i, j, n;

	// 检查TNX关键字
	fits_read_key(fitsptr, TSTRING, "CTYPE1", CTYPE1, NULL, &status);
	fits_read_key(fitsptr, TSTRING, "CTYPE2", CTYPE2, NULL, &status);
	if (status || strcasecmp(CTYPE1, "RA---TNX") || strcasecmp(CTYPE2, "DEC--TNX"))
		return -2;

	// 检查参考点及转换矩阵
	fits_read_key(fitsptr, TDOUBLE, "CRVAL1", &model_.ref_wcsx, NULL, &status);
	fits_read_key(fitsptr, TDOUBLE, "CRVAL2", &model_.ref_wcsx, NULL, &status);
	fits_read_key(fitsptr, TDOUBLE, "CRPIX1", &model_.ref_pixx, NULL, &status);
	fits_read_key(fitsptr, TDOUBLE, "CRPIX2", &model_.ref_pixy, NULL, &status);
	fits_read_key(fitsptr, TDOUBLE, "CD1_1",  &model_.cd[0][0], NULL, &status);
	fits_read_key(fitsptr, TDOUBLE, "CD1_2",  &model_.cd[0][1], NULL, &status);
	fits_read_key(fitsptr, TDOUBLE, "CD2_1",  &model_.cd[1][0], NULL, &status);
	fits_read_key(fitsptr, TDOUBLE, "CD2_2",  &model_.cd[1][1], NULL, &status);
	if (status) return -3;
	// 生成逆转换矩阵
	AMath math;
	memcpy(&model_.ccd[0][0], &model_.cd[0][0], sizeof(model_.cd));
	math.MatrixInvert(2, &model_.ccd[0][0]);
	// 读取畸变改正项
	fits_read_key(fitsptr, TSTRING, "WAT1_001", value, NULL, &status);
	fits_read_key(fitsptr, TSTRING, "WAT2_001", value, NULL, &status);
	fits_close_file(fitsptr, &status);
	if (status) return 0;

	// 解析残差修正模型
	bool rslt;
	for (j = 1; j <= 2; ++j) {
		ptr = &strcor[j - 1];
		i = n = 0;
		while(!status) {
			sprintf(keyword, "WAT%d_%03d", j, ++i);
			fits_read_key(fitsptr, TSTRING, keyword, value, NULL, &status);
			if (!status) {
				/* 作为字符串读出时, 结尾空格被删除, 需要手动填充, 避免数据被截断 */
				if (strlen(value) == 68) n += sprintf((*ptr) + n, "%s", value);
				else n += sprintf((*ptr) + n, "%s ", value);
			}
		}
		status = 0;
	}
	rslt = !(resolve_residual(&strcor[0][0], &model_.res[0])
			|| resolve_residual(&strcor[1][0], &model_.res[1]));

	return rslt ? 0 : -5;
}

bool WCSTNX::LoadText(const char* filepath) {
	FILE *fp = fopen(filepath, "r");
	if (!fp) return false;
	const int size(100);
	char line[size];
	char seps[] = " \t\r\n";
	char *token, *token1;

	while(!feof(fp)) {
		if (fgets(line, size, fp) == NULL) continue;
		token = strtok(line, seps);

		if      (!strcmp(token, "xpixref")) model_.ref_pixx = atof(strtok(NULL, seps));
		else if (!strcmp(token, "ypixref")) model_.ref_pixy = atof(strtok(NULL, seps));
		else if (!strcmp(token, "lngref"))  model_.ref_wcsx = atof(strtok(NULL, seps)) * D2R;
		else if (!strcmp(token, "latref"))  model_.ref_wcsy = atof(strtok(NULL, seps)) * D2R;
		else if (token[0] == 's') {
			int srfc = !strcmp(token, "surface1") ? 1 : (!strcmp(token, "surface2") ? 2 : 0);
			if (srfc) {
				wcs_tnx *tnx = srfc == 1 ? &param_.tnx1[0] : &param_.tnx2[0];
				int i, j, n;

				n = atoi(strtok(NULL, seps));
				for (i = j = 0; i < n && !feof(fp); ++i) {
					fgets(line, size, fp);
					token = strtok(line, seps);
					token1= strtok(NULL, seps);

					if (i > 7) {
						tnx[0].coef[j] = atof(token);
						tnx[1].coef[j] = atof(token1);
						++j;
					}
					else if (i == 0) {
						tnx[0].function = int(atof(token)  + 0.5);
						tnx[1].function = int(atof(token1) + 0.5);
					}
					else if (i == 1) {
						tnx[0].set_orderx(int(atof(token)  + 0.5));
						tnx[1].set_orderx(int(atof(token1) + 0.5));
					}
					else if (i == 2) {
						tnx[0].set_ordery(int(atof(token)  + 0.5));
						tnx[1].set_ordery(int(atof(token1) + 0.5));
					}
					else if (i == 3) {
						tnx[0].set_xterm(int(atof(token)  + 0.5));
						tnx[1].set_xterm(int(atof(token1) + 0.5));
					}
					else if (i == 4) { tnx[0].xmin = atof(token); tnx[1].xmin = atof(token1); }
					else if (i == 5) { tnx[0].xmax = atof(token); tnx[1].xmax = atof(token1); }
					else if (i == 6) { tnx[0].ymin = atof(token); tnx[1].ymin = atof(token1); }
					else if (i == 7) { tnx[0].ymax = atof(token); tnx[1].ymax = atof(token1); }
				}
				param_.valid[srfc - 1] = n && i == n
						&& TNX_CHEBYSHEV <= tnx[0].function && tnx[0].function <= TNX_LINEAR
						&& TNX_CHEBYSHEV <= tnx[1].function && tnx[1].function <= TNX_LINEAR
						&& tnx[0].xorder > 0 && tnx[0].yorder > 0
						&& tnx[1].xorder > 0 && tnx[1].yorder > 0
						&& TNX_XNONE <= tnx[0].xterm && tnx[0].xterm <= TNX_XHALF
						&& TNX_XNONE <= tnx[1].xterm && tnx[1].xterm <= TNX_XHALF
						&& tnx[0].xmax > tnx[0].xmin && tnx[0].ymax > tnx[0].ymin
						&& tnx[1].xmax > tnx[1].xmin && tnx[1].ymax > tnx[1].ymin;
			}
		}
	}
	fclose(fp);

	if (param_.valid[0]) {// 计算旋转矩阵
		wcs_tnx *tnx = &param_.tnx1[0];
		param_.cd[0] = 2 * tnx[0].coef[1] * AS2D / (tnx[0].xmax - tnx[0].xmin);
		param_.cd[1] = 2 * tnx[0].coef[2] * AS2D / (tnx[0].ymax - tnx[0].ymin);
		param_.cd[2] = 2 * tnx[1].coef[1] * AS2D / (tnx[1].xmax - tnx[1].xmin);
		param_.cd[3] = 2 * tnx[1].coef[2] * AS2D / (tnx[1].ymax - tnx[1].ymin);
		if (invert_matrix()) param_.valid[0] = false;
	}

	return (param_.valid[0]);
}

int WCSTNX::WriteImage(const char* filepath) {
	if (!param_.valid[0]) return -1; // 至少需要线性项
	fitsfile *fitsptr;
	if (!hfits(filepath, 1)) return -2; // FITS文件打开失败

	string WCSASTRM = "ct4m.19990714T012701 (USNO-K V) by F. Valdes 1999-08-02";
	string CTYPE1 = "RA---TNX";
	string CTYPE2 = "DEC--TNX";
	string WAT0_001 = "system=image";
	string WAT1 = "wtype=tnx axtype=ra lngcor = ";
	string WAT2 = "wtype=tnx axtype=dec latcor = ";
	int WCSDIM = 2;
	double CRVAL1 = model_.ref_wcsx;
	double CRVAL2 = model_.ref_wcsy;
	double CRPIX1 = model_.ref_pixx;
	double CRPIX2 = model_.ref_pixy;
	double CD1_1  = model_.cd[0][0];
	double CD1_2  = model_.cd[0][1];
	double CD2_1  = model_.cd[1][0];
	double CD2_2  = model_.cd[1][1];
	int status(0);

	if (param_.valid[1]) {
		wcs_tnx *tnx = &param_.tnx2[0];
		int n, nc, i, j;
		char strcor[2][2048];
		char txtdbl[20];
		char (*ptr)[2048];
		for (j = 0, tnx = &param_.tnx2[0]; j < 2; ++j, ++tnx) {
			ptr = &strcor[j];
			n = sprintf(*ptr, "\"%d %d %d %d ", tnx->function, tnx->xorder, tnx->yorder, tnx->xterm);
			output_precision_double(txtdbl, tnx->xmin);
			n += sprintf((*ptr) + n, "%s ", txtdbl);
			output_precision_double(txtdbl, tnx->xmax);
			n += sprintf((*ptr) + n, "%s ", txtdbl);
			output_precision_double(txtdbl, tnx->ymin);
			n += sprintf((*ptr) + n, "%s ", txtdbl);
			output_precision_double(txtdbl, tnx->ymax);
			n += sprintf((*ptr) + n, "%s ", txtdbl);

			nc = tnx->ncoef;
			for (i = 0; i < nc; ++i) {
				output_precision_double(txtdbl, tnx->coef[i]);
				n += sprintf((*ptr) + n, "%s ", txtdbl);
			}
			n += sprintf((*ptr) + n, "\"");
			(*ptr)[n] = 0;
		}

		// 合并修正项
		WAT1 += strcor[0];
		WAT2 += strcor[1];
	}

	// 写入FITS头
	fits_update_key(fitsptr, TSTRING, "WCSASTRM", (void*) WCSASTRM.c_str(), "WCS Source",         &status);
	fits_update_key(fitsptr, TINT,    "WCSDIM",   &WCSDIM,                  "WCS dimensionality", &status);
	fits_update_key(fitsptr, TSTRING, "CTYPE1",   (void*) CTYPE1.c_str(),   "Coordinate type",    &status);
	fits_update_key(fitsptr, TSTRING, "CTYPE2",   (void*) CTYPE2.c_str(),   "Coordinate type",    &status);
	fits_update_key(fitsptr, TDOUBLE, "CRVAL1",   &CRVAL1, "Coordinate reference value", &status);
	fits_update_key(fitsptr, TDOUBLE, "CRVAL2",   &CRVAL2, "Coordinate reference value", &status);
	fits_update_key(fitsptr, TDOUBLE, "CRPIX1",   &CRPIX1, "Coordinate reference pixel", &status);
	fits_update_key(fitsptr, TDOUBLE, "CRPIX2",   &CRPIX2, "Coordinate reference pixel", &status);
	fits_update_key(fitsptr, TDOUBLE, "CD1_1",    &CD1_1,  "Coordinate matrix",          &status);
	fits_update_key(fitsptr, TDOUBLE, "CD1_2",    &CD1_2,  "Coordinate matrix",          &status);
	fits_update_key(fitsptr, TDOUBLE, "CD2_1",    &CD2_1,  "Coordinate matrix",          &status);
	fits_update_key(fitsptr, TDOUBLE, "CD2_2",    &CD2_2,  "Coordinate matrix",          &status);
	if (param_.valid[1]) {// 畸变改正项
		char item[70]; // 每行实际可存储数据68字节
		char keyword[10]; // 关键字
		int i, j, byteleft, byteitem(68);
		int len[] = { WAT1.size(), WAT2.size() };
		const char *head[] = { WAT1.c_str(), WAT2.c_str() };
		const char *ptr;

		fits_update_key(fitsptr, TSTRING, "WAT0_001", (void*) WAT0_001.c_str(), "Coordinate system", &status);
		for (j = 0; j < 2; ++j) {// 0: R.A.; 1: DEC.
			for (byteleft = len[j], i = 0, ptr = head[j]; byteleft; ptr += byteitem) {
				sprintf(keyword, "WAT%d_%03d", j + 1, ++i);
				if (byteleft > byteitem) {
					strncpy(item, ptr, byteitem);
					item[byteitem] = 0;
					fits_update_key(fitsptr, TSTRING, keyword, (void*) item, "", &status);
				}
				else {
					fits_update_key(fitsptr, TSTRING, keyword, (void*) ptr, "", &status);
				}
				byteleft = byteleft > byteitem ? byteleft - byteitem : 0;
			}
		}
	}
	fits_close_file(fitsptr, &status);

	return status;
}

int WCSTNX::resolve_residual(char *strcor, PrjTNXRes *res)  {
	char *pstr, *ptr;
	const char seps[] = " ";
	const char flags[] = "1234567890-+.";
	int i, nitem;
	double *coef;

	if (NULL == (pstr = strpbrk(strcor, flags))) return 1;
	res->func    = int(atof(strtok_r(pstr, seps, &ptr)) + 0.5);
	res->xorder  = (pstr = strtok_r(NULL, seps, &ptr)) != NULL ? int(atof(pstr) + 0.5) : 0;
	res->yorder  = (pstr = strtok_r(NULL, seps, &ptr)) != NULL ? int(atof(pstr) + 0.5) : 0;
	res->xterm   = (pstr = strtok_r(NULL, seps, &ptr)) != NULL ? int(atof(pstr) + 0.5) : 0;
	res->xmin    = (pstr = strtok_r(NULL, seps, &ptr)) != NULL ? atof(pstr) : 0.0;
	res->xmax    = (pstr = strtok_r(NULL, seps, &ptr)) != NULL ? atof(pstr) : 0.0;
	res->ymin    = (pstr = strtok_r(NULL, seps, &ptr)) != NULL ? atof(pstr) : 0.0;
	res->ymax    = (pstr = strtok_r(NULL, seps, &ptr)) != NULL ? atof(pstr) : 0.0;
	if (!res->Initialize()) return 2;

	nitem = res->nitem;
	coef  = res->coef;
	for (i = 0; i < nitem && (pstr = strtok_r(NULL, seps, &ptr)) != NULL; ++i, ++coef)
		*coef = atof(pstr);

	return (i == nitem) ? 0 : 3;
}

int WCSTNX::output_precision_double(char *output, double value) {
	int n(0), pos(0), valid(0), intv;
	if (value < 0) {
		output[0] = '-';
		++pos;
		value = -value;
	}

	// 整数部分
	intv = int(value);
	n = sprintf(output + pos, "%d.", intv);
	pos += n;
	valid = n - 1;
	// 小数部分
	while(valid < 17 && (value = (value - intv) * 10) > AEPS) {
		intv = int(value);
		sprintf(output + pos, "%d", intv);
		++pos;
		++valid;
	}
	output[pos] = 0;

	return n;
}

bool WCSTNX::PrepareFit(double refx, double refy) {
	if (!(model_.res[0].Initialize() && model_.res[1].Initialize())) return false;

	stars_.clear();
	model_.ref_pixx = refx;
	model_.ref_pixy = refy;

	return true;
}

void WCSTNX::AddSample(const MatchedStar &matstar) {
	stars_.push_back(matstar);
}

int WCSTNX::ProcessFit() {
	if (stars_.size() <= model_.res[0].nitem) return 1;	// 样本不足

	PrjTNXRes *res = model_.res;
	double refx(model_.ref_pixx), refy(model_.ref_pixy);
	bool xauto = refx < res->xmin || refy > res->xmax;
	bool yauto = refy < res->ymin || refy > res->ymax;
	// 选择中心点作为投影中心
	if (xauto) model_.ref_pixx = (res->xmin + res->xmax) * 0.5;
	if (yauto) model_.ref_pixy = (res->ymin + res->ymax) * 0.5;
	// 选择最接近参考点的参考星作为投影中心. (refx, refy)保存输入值
	find_nearest(model_.ref_pixx, model_.ref_pixy, model_.ref_wcsx, model_.ref_wcsy);
	if (!try_fit()) return 2;
	if (!(xauto && yauto)) {// 使用新的参考点重新拟合模型
		double x, y, refr, refd;
		x = xauto ? model_.ref_pixx : refx;
		y = yauto ? model_.ref_pixy : refy;
		model_.Image2WCS(x, y, refr, refd);
		model_.ref_pixx = x;
		model_.ref_pixy = y;
		model_.ref_wcsx = refr;
		model_.ref_wcsy = refd;
		try_fit();
	}
	// 计算拟合结果
	AMath math;
	// 计算逆旋转矩阵
	memcpy(&model_.ccd[0][0], &model_.cd[0][0], sizeof(model_.ccd));
	math.MatrixInvert(2, &model_.ccd[0][0]);
	// 计算倾角
	model_.rotation = atan2(model_.cd[0][1], model_.cd[0][0]) * R2D;
	// 计算像元比例尺
	double tmp[2][2];
	memcpy(&tmp[0][0], &model_.cd[0][0], sizeof(model_.ccd));
	model_.scale = R2AS * sqrt(math.LUDet(2, &tmp[0][0]));
	// 计算拟合残差
	double esum(0.0), esq(0.0);
	int i, n(stars_.size());
	double ra, dc, t;

	for (i = 0; i < n; ++i) {
		model_.Image2WCS(stars_[i].x, stars_[i].y, ra, dc);
		t = SphereRange(ra * D2R, dc * D2R, stars_[i].ra * D2R, stars_[i].dc * D2R);
		esum += t;
		esq  += (t * t);
	}
	model_.errfit = sqrt((esq - esum * esum / n) / n) * R2AS;

	return 0;
}

void WCSTNX::find_nearest(double &refx, double &refy, double &refr, double &refd) {
	double dxy, dxymin(1E30);
	int i, j, n(stars_.size());

	for (i = 0; i < n; ++i) {
		dxy = fabs(stars_[i].x - refx) + fabs(stars_[i].y - refy);
		if (dxy < dxymin) {
			dxymin = dxy;
			j = i;
		}
	}
	refx = stars_[j].x;
	refy = stars_[j].y;
	refr = stars_[j].ra;
	refd = stars_[j].dc;
}

bool WCSTNX::try_fit() {
	bool rslt;
	int n(stars_.size()), i;
	double *X, *Y, xi, eta;
	double x0(model_.ref_pixx), y0(model_.ref_pixy);
	AMath math;

	X = (double*) calloc(2*n, sizeof(double));
	Y = (double*) calloc(2*n, sizeof(double));
	for (i = 0; i < n; ++i) {
		model_.WCS2Plane(stars_[i].ra * D2R, stars_[i].dc * D2R, xi, eta);
		X[i]     = stars_[i].x - x0;
		X[n + i] = stars_[i].y - y0;
		Y[i]     = xi * R2D;
		Y[n + i] = eta * R2D;
	}
	// 拟合转换矩阵
	rslt = math.LSFitLinear(n, 2, X, Y, &model_.cd[0][0])
			&& math.LSFitLinear(n, 2, X, Y + n, &model_.cd[1][0]);

	if (rslt) {// 拟合畸变残差
		PrjTNXRes *res = model_.res;
		int nitem = model_.res[0].nitem;
		double *X1 = (double *) calloc(nitem * n, sizeof(double));
		double *X2 = (double *) calloc(nitem * n, sizeof(double));
		double *ptr;

		for (i = 0, ptr = X2; i < n; ++i, ptr += nitem) {
			model_.Image2Plane(stars_[i].x, stars_[i].y, xi, eta);
			Y[i]     = (Y[i] * D2R - xi) * R2AS;
			Y[n + i] = (Y[n + i] * D2R - eta) * R2AS;
			res->ItemVector(stars_[i].x, stars_[i].y, ptr);
		}
		math.MatrixTranspose(n, nitem, X2, X1);
		rslt = math.LSFitLinear(n, nitem, X1, Y, res[0].coef)
				&& math.LSFitLinear(n, nitem, X1, Y + n, res[1].coef);

		free(X1);
		free(X2);
	}

	free(X);
	free(Y);

	return rslt;
}
//////////////////////////////////////////////////////////////////////////////
} /* namespace AstroUtil */
