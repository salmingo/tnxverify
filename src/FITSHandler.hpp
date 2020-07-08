/**
 * @file FITSHandler 基于cfitsio封装FITS文件基本操作
 */

#ifndef FITS_HANDLER_H_
#define FITS_HANDLER_H_

#include <stdlib.h>
#include <longnam.h>
#include <fitsio.h>

struct FITSHandler {
	fitsfile *fitsptr;	//< 基于cfitsio接口的文件操作接口
	int errcode;		//< 错误代码
	char errmsg[100];	//< 错误提示

protected:
	int hdunum;		//< HDU数量
	int *hdutype;	//< HDU类型

protected:
	void update_hdunum(int n) {
		if (n != hdunum) {
			hdunum = n;
			if (hdutype) {
				free(hdutype);
				hdutype = NULL;
			}
		}
		if (hdutype == NULL && n) hdutype = (int*) calloc(n, sizeof(int));
	}

public:
	FITSHandler() {
		fitsptr = NULL;
		errcode = 0;
		hdunum  = 0;
		hdutype = NULL;
	}

	virtual ~FITSHandler() {
		Close();
		if (hdutype) free(hdutype);
	}

	bool Close() {
		errcode = 0;
		if (fitsptr)  fits_close_file(fitsptr, &errcode);
		if (!errcode) fitsptr = NULL;
		return !errcode;
	}

	const char *GetError() {
		fits_get_errstatus(errcode, errmsg);
		return &errmsg[0];
	}

	const int* HDUType(int &n) {
		n = hdunum;
		return n == 0 ? NULL : &hdutype[0];
	}

	/*!
	 * @brief 获得基于cfitsio的FITS文件访问接口
	 * @return
	 * FITS文件访问接口
	 */
	fitsfile* operator()() {
		return fitsptr;
	}

	/*!
	 * @brief 打开或创建FITS文件并获得基于cfitsio的FITS文件访问接口
	 * @param pathname    文件路径
	 * @param mode        0: 只读; 1: 读写; 2: 新文件
	 * @return
	 * 操作结果
	 */
	bool operator()(const char *filepath, const int mode = 0) {
		if (!Close()) return false;
		errcode = 0;
		if (mode == 0 || mode == 1) {
			int n, i;

			fits_open_file(&fitsptr, filepath, mode, &errcode);
			fits_get_num_hdus(fitsptr, &n, &errcode);
			if (!errcode) {
				update_hdunum(n);
				for (i = 0; i < n; ++i) {
					fits_movabs_hdu(fitsptr, i + 1, hdutype + i, &errcode);
				}
			}
		}
		else {
			fits_create_file(&fitsptr, filepath, &errcode);
			update_hdunum(0);
		}
		return !errcode;
	}

	/*!
	 * @brief 改变当前HDU
	 * @param idx 索引, [1, hdunum]
	 * @return
	 * 改变后HDU的头区中关键字数量
	 */
	int MovetoHDU(int idx) {
		if (idx < 1 || idx > hdunum) return 0;
		int n0, n1;
		errcode = 0;
		fits_movabs_hdu(fitsptr, idx, hdutype + idx - 1, &errcode);
		if (!errcode) fits_get_hdrspace(fitsptr, &n0, &n1, &errcode);
		return (errcode ? 0 : n0);
	}
};
typedef FITSHandler HFITS;
typedef FITSHandler* HFITSPtr;

#endif
