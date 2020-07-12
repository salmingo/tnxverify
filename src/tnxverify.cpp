/**
 Name        : tnxverify.cpp
 Author      : Xiaomeng Lu
 Description : 验证TNX投影精度
 Version     : 0.1
 - 从SEx结果读取XY
 - 从Ast结果读取WCS
 - 由XY和WCS计算(ra1, dc1)
 - 从UCAC4星表读取坐标(ra2, dc2), 匹配精度?角秒且具有唯一性
 - ProjectTNX拟合(x,y)和(ra2,dc2), 建立TNX模型
 - 由XY和TNX模型计算(ra3, dc3)
 - 验证模型可行性和精度
 - 由(ra3, dc3)匹配UCAC4星表, 验证与WCS模型的提升程度

 Version     : 0.2
 - 从Ast生成的.xyls和.rdls中提取(x1,y1)和(ra1,dc1)
 - 使用(x1,y1)和(ra1,dc1)拟合生成TNX模型, 评估模型精度
 - 从SEx读取(x2,y2)
 - 从Ast生成的.wcs中提取WCS模型
 - 由WCS模型计算(ra2,dc2)
 - 由TNX模型计算(x2,y2)对应的(ra3,dc3)
 - 从UCAC4中查找与(ra3, dc3)匹配的参考星
 - 评估定位精度；输出匹配恒星信息
 **/

#include <cstdio>
#include <boost/smart_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <longnam.h>
#include <fitsio.h>
#include "ADefine.h"
#include "ADIData.h"
#include "ACatUCAC4.h"
#include "FITSHandler.hpp"
#include "WCSTNX.h"

using namespace std;
using namespace boost::filesystem;
using namespace boost::posix_time;
using namespace AstroUtil;

/*--------------------------------------------------------------------------*/
struct wcsinfo {
	double x0, y0;	//< XY参考点
	double r0, d0;	//< RA/DEC参考点, 量纲: 弧度
	double cd[2][2];	//< 转换矩阵
	int orderA, orderB;	//< SIP改正阶数
	int ncoefA, ncoefB;	//< SIP系数数量
	double *A, *B;	//< 线性改正系数

public:
	wcsinfo() {
		x0 = y0 = 0.0;
		r0 = d0 = 0.0;
		orderA = orderB = 0;
		ncoefA = ncoefB = 0;
		A = B = NULL;
	}

	virtual ~wcsinfo() {
		if (A) {
			delete[] A;
			A = NULL;
		}
		if (B) {
			delete[] B;
			B = NULL;
		}
	}

	wcsinfo &operator=(const wcsinfo &other) {
		if (this != &other) {
			x0   = other.x0;
			y0   = other.y0;
			r0   = other.r0;
			d0   = other.d0;
			orderA = other.orderA;
			orderB = other.orderB;
			ncoefA = other.ncoefA;
			ncoefB = other.ncoefB;
			memcpy(&cd[0][0], &other.cd[0][0], sizeof(double) * 4);
			if (ncoefA) {
				A = new double[ncoefA];
				memcpy(A, other.A, sizeof(double) * ncoefA);
			}
			if (ncoefB) {
				B = new double[ncoefB];
				memcpy(B, other.B, sizeof(double) * ncoefB);
			}
		}

		return *this;
	}

protected:
	/*!
	 * @brief 计算SIP改正模型中与阶数对应的系数数量
	 * @return
	 * 系数数量
	 */
	int term_count(int order) {
		return (order + 1) * (order + 2) / 2;
	}

	/*!
	 * @brief 为SIP改正系数分配存储空间
	 * @param n   系数数量
	 * @param ptr 系数存储地址
	 */
	void alloc_coef(int n, double **ptr) {
		if ((ptr == &A && n != ncoefA) || (ptr == &B && n != ncoefA)) {
			if (ptr == &A)
				ncoefA = n;
			else
				ncoefB = n;
			if (*ptr) {
				delete[] (*ptr);
				(*ptr) = NULL;
			}
		}
		if (*ptr == NULL)
			(*ptr) = new double[n];
	}

	/*!
	 * @brief 参考平面坐标转换为世界坐标
	 * @param xi    参考平面X坐标, 量纲: 弧度
	 * @param eta   参考平面Y坐标, 量纲: 弧度
	 * @param newr  参考点(x0, y0)对应的新的赤经, 量纲: 弧度
	 * @param newd  参考点(x0, y0)对应的新的赤纬, 量纲: 弧度
	 * @param ra    世界坐标赤经, 量纲: 弧度
	 * @param dec   世界坐标赤纬, 量纲: 弧度
	 */
	void plane_to_wcs(double xi, double eta, double newr, double newd, double &ra, double &dec) {
		double fract = cos(newd) - eta * sin(newd);
		ra = cyclemod(newr + atan2(xi, fract), A2PI);
		dec = atan2(((eta * cos(newd) + sin(newd)) * cos(ra - newr)), fract);
	}

	double poly_val(double x, double y, double *coef, int order) {
		int i, j, k, m;
		double val(0.0), t, px(1.0), py;

		for (i = 0, k = 0; i <= order; ++i) {
			for (j = 0, py = 1.0, t = 0.0, m = order - i; j <= m; ++j, ++k) {
				t += coef[k] * py;
				py *= y;
			}

			val += t * px;
			px *= x;
		}

		return val;
	}

	void project_correct(double &x, double &y) {
		double dx(0.0), dy(0.0);
		dx = poly_val(x, y, A, orderA);
		dy = poly_val(x, y, B, orderB);
		x += dx;
		y += dy;
	}

public:
	bool load_wcs(const string &filepath) {
		fitsfile *fitsptr;	//< 基于cfitsio接口的文件操作接口
		char key[10];
		int status(0), ncoef, i, j, k, m;

		fits_open_file(&fitsptr, filepath.c_str(), 0, &status);
		fits_read_key(fitsptr, TDOUBLE, "CRPIX1", &x0, NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CRPIX2", &y0, NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CRVAL1", &r0, NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CRVAL2", &d0, NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CD1_1", &cd[0][0], NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CD1_2", &cd[0][1], NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CD2_1", &cd[1][0], NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CD2_2", &cd[1][1], NULL, &status);

		fits_read_key(fitsptr, TINT, "A_ORDER", &orderA, NULL, &status);
		if (status)
			return false;
		ncoef = term_count(orderA);
		alloc_coef(ncoef, &A);
		for (i = 0, k = 0; i <= orderA; ++i) {
			for (j = 0, m = orderA - i; j <= m; ++j, ++k) {
				sprintf(key, "A_%d_%d", i, j);
				fits_read_key(fitsptr, TDOUBLE, key, A + k, NULL, &status);
			}
		}

		fits_read_key(fitsptr, TINT, "B_ORDER", &orderB, NULL, &status);
		if (status)
			return false;
		ncoef = term_count(orderB);
		alloc_coef(ncoef, &B);
		for (i = 0, k = 0; i <= orderB; ++i) {
			for (j = 0, m = orderB - i; j <= m; ++j, ++k) {
				sprintf(key, "B_%d_%d", i, j);
				fits_read_key(fitsptr, TDOUBLE, key, B + k, NULL, &status);
			}
		}

		fits_close_file(fitsptr, &status);

		r0 *= D2R;
		d0 *= D2R;
		return !status;
	}

	/*!
	 * @brief 图像坐标转换为世界坐标
	 * @param x     图像X坐标, 量纲: 像素
	 * @param y     图像Y坐标, 量纲: 像素
	 * @param newr  参考点(x0, y0)对应的新的赤经, 量纲: 弧度
	 * @param newd  参考点(x0, y0)对应的新的赤纬, 量纲: 弧度
	 * @param ra    世界坐标赤经, 量纲: 角度
	 * @param dec   世界坐标赤纬, 量纲: 角度
	 */
	void image_to_wcs(double x, double y, double newr, double newd, double &ra, double &dec) {
		double xi, eta;
		x -= x0;
		y -= y0;
		project_correct(x, y);
		xi = (cd[0][0] * x + cd[0][1] * y) * D2R;
		eta = (cd[1][0] * x + cd[1][1] * y) * D2R;
		plane_to_wcs(xi, eta, newr, newd, ra, dec);
		ra *= R2D;
		dec *= R2D;
	}

	void image_to_wcs(double x, double y, double &ra, double &dec) {
		image_to_wcs(x, y, r0, d0, ra, dec);
	}
};

/*--------------------------------------------------------------------------*/
bool all_exists(const string &filepath) {
	path pathname = filepath;
	pathname.replace_extension(".cat");
	if (!exists(pathname)) return false;
	pathname.replace_extension(".wcs");
	if (!exists(pathname)) return false;

	return true;
}

/*!
 * @brief 从FITS头提取曝光起始时间
 * @param filepath  文件路径
 * @param frame     数据结构
 */
void resolve_dateobs(const string &filepath, ImgFrmPtr frame) {
	path pathname = filepath;
	fitsfile *fitsptr;	//< 基于cfitsio接口的文件操作接口
	int status(0);
	char dateobs[30], timeobs[30], tmfull[40];
	bool datefull;

	fits_open_file(&fitsptr, filepath.c_str(), 0, &status);
	fits_read_key(fitsptr, TINT, "NAXIS1", &frame->wdim, NULL, &status);
	fits_read_key(fitsptr, TINT, "NAXIS2", &frame->hdim, NULL, &status);
	fits_read_key(fitsptr, TSTRING, "DATE-OBS", dateobs,  NULL, &status);
	if (!(datefull = NULL != strstr(dateobs, "T")))
		fits_read_key(fitsptr, TSTRING, "TIME-OBS", timeobs,  NULL, &status);
	fits_close_file(fitsptr, &status);

	if (!datefull) sprintf(tmfull, "%sT%s", dateobs, timeobs);
	frame->filename = pathname.filename().string();
	frame->dateobs  = datefull ? dateobs : tmfull;

	ptime tmobs  = from_iso_extended_string(frame->dateobs);
	frame->mjdobs = tmobs.date().modjulian_day(); // 精确到天. 精度满足自行
}

/*!
 * @brief 从Ast生成的.xyls和.rdls文件中提取XY和赤道坐标
 * @param filepath  文件路径
 * @param frame     数据结构
 */
void refstar_from_list(const string &filepath, WCSTNX& wcstnx) {
	FITSHandler hfit;
	int status(0);
	int nhdu;
	const int* hdutype;
	long nrows;
	double *ra, *dc, *x, *y;

	// 从.rdls读取赤经赤纬
	path rdls = filepath;
	rdls.replace_extension(".rdls");
	hfit(rdls.c_str());
	fits_get_num_rows(hfit(), &nrows, &status);
	ra = new double[nrows];
	dc = new double[nrows];
	x  = new double[nrows];
	y  = new double[nrows];
	fits_read_col(hfit(), TDOUBLE, 1, 1, 1, nrows, NULL, ra, NULL, &status);
	fits_read_col(hfit(), TDOUBLE, 2, 1, 1, nrows, NULL, dc, NULL, &status);

	// 从.xyls读取XY
	path filename = path(filepath).stem();
	path xyls = path(filepath).parent_path();

	filename += "-indx.xyls";
	xyls /= filename;
	hfit(xyls.c_str());
	fits_read_col(hfit(), TDOUBLE, 1, 1, 1, nrows, NULL, x, NULL, &status);
	fits_read_col(hfit(), TDOUBLE, 2, 1, 1, nrows, NULL, y, NULL, &status);

	// 填充参考星
	wcstnx.PrepareFit();
	for (long i = 0; i < nrows; ++i) {
		MatchedStar star;
		star.x   = x[i];
		star.y   = y[i];
		star.ra  = ra[i];
		star.dc  = dc[i];
		wcstnx.AddSample(star);
	}

	delete []ra;
	delete []dc;
	delete []x;
	delete []y;
}

void refstar_from_frame(ImgFrmPtr frame, WCSTNX& wcstnx) {
	NFObjVector &nfobj = frame->nfobj;
	wcstnx.PrepareFit();

	for (NFObjVector::iterator x = nfobj.begin(); x != nfobj.end(); ++x) {
		if (x->matched != 1) continue;
		MatchedStar star;
		star.x   = x->ptbc.x;
		star.y   = x->ptbc.y;
		star.ra  = x->ra_cat;
		star.dc  = x->dec_cat;
		wcstnx.AddSample(star);
	}
}

/*!
 * @brief 提取XY坐标
 * @param filepath  文件路径
 * @param frame     数据结构
 */
void take_xy(const string &filepath, ImgFrmPtr frame) {
	path pathname = filepath;
	FILE *fp;
	char line[200];
	double x, y, flux;
	NFObjVector &nfobj = frame->nfobj;
	nfobj.clear();

	pathname.replace_extension(".cat");
	fp = fopen(pathname.c_str(), "r");
	while (!feof(fp)) {
		if (NULL == fgets(line, 200, fp) || line[0] == '#') continue;

		ObjectInfo obj;
		sscanf(line, "%lf %lf %lf", &x, &y, &flux);
		obj.ptbc.x = x;
		obj.ptbc.y = y;
		obj.flux   = flux;
		nfobj.push_back(obj);
	}
	fclose(fp);

	frame->nobjs = frame->nfobj.size();
	printf (">> %d objects are found\n", frame->nobjs);
}

/*!
 * @brief 由WCS计算(ra1, dc1)
 * @param frame  数据结构
 */
void rd_from_wcs(const string &filepath, ImgFrmPtr frame) {
	path pathname = filepath;
	wcsinfo wcs;
	NFObjVector &nfobj = frame->nfobj;

	pathname.replace_extension(".wcs");
	wcs.load_wcs(pathname.string());
	for (NFObjVector::iterator x = nfobj.begin(); x != nfobj.end(); ++x) {
		wcs.image_to_wcs(x->ptbc.x, x->ptbc.y, x->ra_inst, x->dec_inst);
	}
}

/*!
 * @brief 由TNX模型拟合(ra3, dc3)
 * @param frame  数据结构
 * @param model  模型
 */
void rd_from_tnx(ImgFrmPtr frame, PrjTNX& model) {
	NFObjVector &nfobj = frame->nfobj;

	for (NFObjVector::iterator x = nfobj.begin(); x != nfobj.end(); ++x) {
		model.Image2WCS(x->ptbc.x, x->ptbc.y, x->ra_fit, x->dec_fit);
		x->ra_fit  *= R2D;
		x->dec_fit *= R2D;
	}
}

/*!
 * @brief 与UCAC4匹配, 提取星表坐标
 * @param frame  数据结构
 */
void match_ucac4(ImgFrmPtr frame) {
	boost::shared_ptr<ACatUCAC4> ucac4 = boost::make_shared<ACatUCAC4>();
	NFObjVector &nfobj = frame->nfobj;
	int n(0), nfound, i, j;
	double r(0.8);	//1角分半径
	double t; // 儒略年 * 0.1毫角秒转角秒系数
	ucac4item_ptr starptr;

	t = (frame->mjdobs - 51544.5) * 1E-4 / 365.25;
	ucac4->SetPathRoot("/Users/lxm/Catalogue/UCAC4");
	for (NFObjVector::iterator x = nfobj.begin(); x != nfobj.end(); ++x) {
		if ((nfound = ucac4->FindStar(x->ra_fit, x->dec_fit, r))) {
			starptr = ucac4->GetResult();
			if (nfound == 1) j = 0;
			else {
				double er, ed, emin(1E30), t, t1, cosd;
				cosd = cos(x->dec_fit * D2R);
				for (i = 0; i < nfound; ++i) {
					ed = starptr[i].spd / MILLISEC - 90.0 - x->dec_fit;
					er = starptr[i].ra / MILLISEC - x->ra_fit;
					if (er > 180.0) er -= 360.0;
					else if (er < -180.0) er += 360.0;
					t1 = er * cosd;
					t = t1 * t1 + ed * ed;
					if (t < emin) {
						j = i;
						emin = t;
					}
				}
			}

			x->dec_pm   = starptr[j].pmdc * t;
			x->dec_cat  = (double) starptr[j].spd / MILLISEC - 90.0 + x->dec_pm * AS2D;
			x->ra_pm    = starptr[j].pmrac * t / cos(x->dec_cat * D2R);
			x->ra_cat   = (double) starptr[j].ra / MILLISEC + x->ra_pm * AS2D;
			x->mag_cat  = starptr[j].apasm[1] * 0.001;
			x->matched = 1;
			++n;
		}
	}
	printf (">> %d stars are matched\n", n);
}

/*!
 * @brief 评估ProjectTNX拟合结构
 * @param frame  数据结构
 */
void final_stat(const string &filepath, ImgFrmPtr frame) {
	path pathname = filepath;
	NFObjVector &nfobj = frame->nfobj;
	FILE *fprslt;
	double er1, er2, ed1, ed2;
	double esq_ra(0.0), esum_ra(0.0), esq_dc(0.0), esum_dc(0.0);
	double mean_ra, mean_dc;
	int n(0);

	pathname.replace_extension(".txt");
	fprslt = fopen(pathname.c_str(), "w");

	fprintf (fprslt, "%6s %6s %8s %8s %7s %8s %8s %8s %8s\n",
			"X", "Y", "RA_CAT", "DEC_CAT", "MAG_CAT", "RA_WCS", "DEC_WCS", "RA_TNX", "DEC_TNX");

	for (NFObjVector::iterator x = nfobj.begin(); x != nfobj.end(); ++x) {
		if (x->matched != 1) continue;
		er1 = x->ra_fit - x->ra_cat;
		ed1 = (x->dec_fit - x->dec_cat) * 3600.0;
		er2 = x->ra_inst - x->ra_fit;
		ed2 = (x->dec_inst - x->dec_fit) * 3600.0;
		if (er1 > 180.0) er1 -= 360.0;
		else if (er1 < -180.0) er1 += 360.0;
		if (er2 > 180.0) er2 -= 360.0;
		else if (er2 < -180.0) er2 += 360.0;
		er1 *= 3600.0;
		er2 *= 3600.0;

		esum_ra += er1;
		esq_ra  += er1 * er1;
		esum_dc += ed1;
		esq_dc  += ed1 * ed1;
		++n;

		fprintf (fprslt, "%6.1f %6.1f %8.4f %8.4f %7.3f %8.4f %8.4f %8.4f %8.4f %5.1f %5.1f | %5.1f %5.1f | %5.1f %5.1f\n",
				x->ptbc.x, x->ptbc.y,
				x->ra_cat, x->dec_cat, x->mag_cat,
				x->ra_inst, x->dec_inst,
				x->ra_fit, x->dec_fit,
				x->ra_pm, x->dec_pm,
				er1, ed1,
				er2, ed2);
	}
	if (n > 2) {
		mean_ra = esum_ra / n;
		mean_dc = esum_dc / n;
		printf ("RA\n\tMean  = %.2f\n\tStdev = %.2f\n",
				mean_ra, sqrt((esq_ra - mean_ra * esum_ra) / n));
		printf ("DEC\n\tMean  = %.2f\n\tStdev = %.2f\n",
				mean_dc, sqrt((esq_dc - mean_dc * esum_dc) / n));
	}

	fclose(fprslt);
}

int main(int argc, char **argv) {
	ImgFrmPtr frame = boost::make_shared<ImageFrame>();
	PrjTNX model;
	WCSTNX wcstnx;
	string extfit(".fit");	// 扩展名
	string filepath;
//	double *ptr;
//	int i, j;

	wcstnx.SetModel(&model);
	model.SetParamRes(FUNC_LEGENDRE, X_FULL, 6, 6);

	for (int i = 1; i < argc; ++i) {
		path pathname = argv[i];
		if (is_directory(pathname)) {
			printf ("--------------------------------------\n");
			printf ("scan directory: %s\n", pathname.c_str());
			for (directory_iterator x = directory_iterator(pathname); x != directory_iterator(); ++x) {
				if (x->path().extension().string() == extfit && all_exists(filepath = x->path().string())) {
					printf ("#### %s ####\n", filepath.c_str());
					resolve_dateobs(filepath, frame);
					model.SetNormalRange(1, 1, frame->wdim, frame->hdim);
					refstar_from_list(filepath, wcstnx);
					if (!wcstnx.ProcessFit()) {
						printf ("FIT 1 >> scale: %.3f, orientation: %.2f, error: %.3f\n",
								model.scale, model.rotation, model.errfit);
//						printf ("Rotation Matrix:\n");
//						for (j = 0, ptr = &model.cd[0][0]; j < 2; ++j) {
//							for (i = 0; i < 2; ++i, ++ptr) {
//								printf ("%G  ", *ptr);
//							}
//							printf ("\n");
//						}
//						printf ("Residual 0:\n");
//						for (i = 0, ptr = model.res[0].coef; i < model.res[0].nitem; ++i, ++ptr)
//							printf ("\t %f\n", *ptr);
//						printf ("Residual 1:\n");
//						for (i = 0, ptr = model.res[1].coef; i < model.res[1].nitem; ++i, ++ptr)
//							printf ("\t %f\n", *ptr);
						take_xy(filepath, frame);
						rd_from_wcs(filepath, frame);
						rd_from_tnx(frame, model);
						match_ucac4(frame);

						refstar_from_frame(frame, wcstnx);
						if (!wcstnx.ProcessFit()) {
							printf ("FIT 2 >> scale: %.3f, orientation: %.2f, error: %.3f\n",
									model.scale, model.rotation, model.errfit);
							rd_from_tnx(frame, model);
							final_stat(filepath, frame);
						}
					}
					printf ("\n");
				}
			}
			printf ("--------------------------------------\n");
		}
		else if (is_regular_file(pathname) && pathname.extension().string() == extfit) {
			filepath = pathname.string();
			if (all_exists(filepath)) {
				printf ("#### %s ####\n", filepath.c_str());
				resolve_dateobs(filepath, frame);
				model.SetNormalRange(1, 1, frame->wdim, frame->hdim);
				refstar_from_list(filepath, wcstnx);
				if (!wcstnx.ProcessFit()) {
					printf ("FIT 1 >> scale: %.3f, orientation: %.2f, error: %.3f\n",
							model.scale, model.rotation, model.errfit);
//					printf ("Rotation Matrix:\n");
//					for (j = 0, ptr = &model.cd[0][0]; j < 2; ++j) {
//						for (i = 0; i < 2; ++i, ++ptr) {
//							printf ("%G  ", *ptr);
//						}
//						printf ("\n");
//					}
//					printf ("Residual 0:\n");
//					for (i = 0, ptr = model.res[0].coef; i < model.res[0].nitem; ++i, ++ptr)
//						printf ("\t %f\n", *ptr);
//					printf ("Residual 1:\n");
//					for (i = 0, ptr = model.res[1].coef; i < model.res[1].nitem; ++i, ++ptr)
//						printf ("\t %f\n", *ptr);

					take_xy(filepath, frame);
					rd_from_wcs(filepath, frame);
					rd_from_tnx(frame, model);
					match_ucac4(frame);

					refstar_from_frame(frame, wcstnx);
					if (!wcstnx.ProcessFit()) {
						printf ("FIT 2 >> scale: %.3f, orientation: %.2f, error: %.3f\n",
								model.scale, model.rotation, model.errfit);
						rd_from_tnx(frame, model);
						final_stat(filepath, frame);
					}
				}
				printf ("\n");
			}
		}
	}

	return 0;
}
