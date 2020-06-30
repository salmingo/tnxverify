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

#ifndef WCSTNX_H_
#define WCSTNX_H_

namespace AstroUtil {
//////////////////////////////////////////////////////////////////////////////
/*---------- 声明相关数据结构 ----------*/
/*!
 * @struct MatchedStar 建立匹配关系的参考星坐标
 */
struct MatchedStar {
	/* 参考星坐标 */
	double x, y;	/// XY坐标
	double ra, dc;	/// RA/DEC坐标, 量纲: 角度
	/* 拟合结果 */
	double ra_fit, dc_fit;	/// RA/DEC坐标, 量纲: 角度
};

enum {// 畸变修正函数类型
	FUNC_LINEAR = 1,	//< 线性
	FUNC_LEGENDRE,		//< 勒让德
	FUNC_CHEBYSHEV		//< 契比雪夫
};

enum {// 多项式交叉系数类型
	X_NONE,	//< 无交叉项
	X_HALF,	//< 半交叉
	X_FULL	//< 全交叉
};

/*!
 * @struct PrjTNXRes  TNX投影残差模型参数
 */
struct PrjTNXRes {
	/*!
	 * @function WCSTNXFunc WCSTNX中模型拟合的基函数
	 */
	typedef void (*WCSTNXFunc)(double, double, double, int, double *);
	/* 用户 */
	int func;	/// 基函数类型
	int xterm;	/// 交叉项类型
	int xorder, yorder;	/// XY自变量的阶次
	double xmin, ymin;	/// 归一化范围: 最小值
	double xmax, ymax;	/// 归一化范围: 最大值
	/* 内部 */
	int nitem;			/// 多项式单元数量
	double *coef;		/// 系数
	double *xv, *yv;	/// 多项式单项存储区

private:
	WCSTNXFunc basefunc;	/// 基函数

/* 构造函数 */
public:
	PrjTNXRes();
	virtual ~PrjTNXRes();

/* 接口 */
public:
	/*!
	 * @brief 设置多项式参数
	 */
	void SetParam(int func, int xterm, int xorder, int yorder);
	/*!
	 * @brief 设置归一化范围
	 */
	void SetRange(double xmin, double ymin, double xmax, double ymax);
	/*!
	 * @brief 从外部设置参数后, 初始化数据结构的存储区
	 * @return
	 * 初始化结果
	 * - 初始化过程中检查用户输入的有效性
	 */
	bool Initialize();
	/*!
	 * @brief 释放存储区
	 */
	void UnInitialize();
	/*!
	 * @brief 由函数基生成多项式各项数值
	 * @param x      自变量X
	 * @param y      自变量Y
	 * @param vptr   由函数基生成的
	 * @return
	 * 生成结果
	 * @note
	 * - vptr空间由用户分配, 其单元数量等于nitem
	 */
	bool ItemVector(double x, double y, double *vptr);
	/*!
	 * @brief 由拟合系数生成多项式数值
	 * @param x  自变量X
	 * @param y  自变量Y
	 * @return
	 * 多项式数值. 当系数无效时返回0.0
	 */
	double PolyVal(double x, double y);

/* 功能 */
protected:
	/*!
	 * @brief 释放缓存区
	 * @param ptr  缓存区地址
	 */
	void free_array(double **ptr);
	/*!
	 * @brief 基函数
	 * @param value   自变量
	 * @param min     最小值
	 * @param max     最大值
	 * @param order   阶次
	 * @param ptr     存储区
	 * @note
	 * - power_array:     幂函数
	 * - legendre_array:  勒让德函数
	 * - chebyshev_array: 切比雪夫函数
	 */
	void power_array(double value, double min, double max, int order, double *ptr);
	void legendre_array(double value, double min, double max, int order, double *ptr);
	void chebyshev_array(double value, double min, double max, int order, double *ptr);
};

/*!
 * @struct PrjTNX TNX投影模型参数
 */
struct PrjTNX {
	double ref_pixx, ref_pixy;	/// 参考点XY坐标
	double ref_wcsx, ref_wcsy;	/// 参考点世界坐标, 量纲: 角度
	double cd[2][2];			/// XY->WCS的转换矩阵, 量纲: 角度/像素
	double ccd[2][2];			/// WCS->XY的转换矩阵, 量纲: 像素/角度
	PrjTNXRes res[2];			/// 残差/畸变参数及模型

/* 接口 */
public:
	/*!
	 * @brief 使用模型, 由XY计算世界坐标
	 * @param x    坐标X
	 * @param y    坐标Y
	 * @param ra   世界坐标X, 量纲: 角度
	 * @param dc   世界坐标Y, 量纲: 角度
	 */
	void Image2WCS(double x, double y, double &ra, double &dc);
	/*!
	 * @brief 使用模型, 由世界坐标计算XY
	 * @param ra   世界坐标X, 量纲: 角度
	 * @param dc   世界坐标Y, 量纲: 角度
	 * @param x    坐标X
	 * @param y    坐标Y
	 */
	void WCS2Image(double ra, double dc, double &x, double &y);

/* 功能 */
protected:
	/*!
	 * @brief 使用转换矩阵, 由XY计算投影平面坐标
	 * @param x    坐标X
	 * @param y    坐标Y
	 * @param xi   投影平面X坐标, 量纲: 角度
	 * @param eta  投影平面Y坐标, 量纲: 角度
	 */
	void image_to_plane(double x, double y, double &xi, double &eta);
};
/*---------- 声明相关数据结构 ----------*/
//////////////////////////////////////////////////////////////////////////////
class WCSTNX {
public:
	WCSTNX();
	virtual ~WCSTNX();
};

} /* namespace AstroUtil */
//////////////////////////////////////////////////////////////////////////////
#endif /* WCSTNX_H_ */
