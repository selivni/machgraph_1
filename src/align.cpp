#include "align.h"
#include <string>

using std::string;
using std::cout;
using std::endl;
using std::tie;
using std::tuple;
using std::make_tuple;

enum
{
	w_crcor,
	w_avdev
};

enum
{
	c_r=0,
	c_g=1,
	c_b=2
};

Image equalize_hist(Image src)
{
	int hist[3][256];
	int cdf[3][256],cdfmin[3];
	int sum,i,j,k;
	int c[3];
	int rows = src.n_rows, cols = src.n_cols;
	for (i=0; i<3; i++)
		for (j=0; j<256; j++)
			hist[i][j] = 0;
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
		{
			tie(c[0],c[1],c[2]) = src(i,j);
			for (k=0; k<3; k++)
				hist[k][c[k]]++;
		}
	for (k=0; k<3; k++)
		for (i=0; i<256; i++)
		{
			sum = 0;
			for (j=0; j<=i; j++)
				sum += hist[k][j];
			cdf[k][i] = sum;
		}
	for (k=0; k<3; k++)
	{
		i=0;
		while (cdf[k][i] == 0 && i < 256)
			i++;
		cdfmin[k] = cdf[k][i];
	}
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
		{
			tie(c[0],c[1],c[2]) = src(i,j);
			for (k=0; k<3; k++)
			{
				c[k] = round(255*(cdf[k][c[k]]-cdfmin[k]+0.0)/(rows*cols-cdfmin[k]));
				if (c[k] > 255)
					c[k] = 255;
			}
			src(i,j) = tie(c[0],c[1],c[2]);
		}
	return src;
}

Image convert(Image src)
{
	Image result(src.n_rows, src.n_cols);
	int i,j,c[3],r[3];
	int rows = src.n_rows, cols = src.n_cols;
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
		{
			tie(c[0],c[1],c[2]) = src(i,j);
			r[0] = round(0.299*c[0]+0.587*c[1]+0.114*c[2]);
			r[1] = round(0.564*(c[2]-r[0]));
			r[2] = round(0.713*(c[0]-r[0]));
			result(i,j) = tie(r[0],r[1],r[2]);
		}
	return result;
}

Image white_balance(Image source)
{
	float Yavg=0,Cravg=0,Cbavg=0,Ravg=0,Gavg=0,Bavg=0;
	float Yl,Crl,Cbl,Yu,Cru,Cbu;
	int Ybgt=-1,Cbbgt=0,Crbgt=0;
	int rows = source.n_rows, cols = source.n_cols;
	Image res(rows,cols);
	Image c_res(rows,cols);
	int i,j,k,c[3],r[3],prob_count=0,ref_count=0;
	bool probable_exist=false,ref_exist=false;
	float Yw=0,Rw=0,Gw=0,Bw=0;
	float Rscale, Gscale, Bscale, Rgwa, Ggwa, Bgwa,R,B,G;
	int colorcast;
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			res(i,j) = source(i,j);
	res = equalize_hist(res);
	c_res = convert(res);
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
		{
			tie (c[0],c[1],c[2]) = c_res(i,j);
			if (c[0] >= 210 && abs(c[1]) <= 3 && abs(c[2]) <= 3)
			{
				tie(r[0],r[1],r[2]) = res(i,j);
				probable_exist = true;
				prob_count++;
				Yavg += c[0];
				Cravg += c[2];
				Cbavg += c[1];
				Ravg += r[0];
				Gavg += r[1];
				Bavg += r[2];
				if (c[0] > Ybgt)
				{
					Ybgt = c[0];
					Cbbgt = c[1];
					Crbgt = c[2];
				} else
				if (c[0] == Ybgt)
					if (c[1]*c[2] < Cbbgt*Crbgt)
					{
						Cbbgt = c[1];
						Crbgt = c[2];
					}
			}
		}
	if (!probable_exist)
		return source;
	Yavg /= prob_count;
	Cravg /= prob_count;
	Cbavg /= prob_count;
	Ravg /= prob_count;
	Gavg /= prob_count;
	Bavg /= prob_count;
	if (Ybgt <= Yavg)
	{
		Yl = Ybgt;
		Yu = Yavg;
	} else
	{
		Yl = Yavg;
		Yu = Ybgt;
	}
	if (Crbgt <= Cravg)
	{
		Crl = Crbgt;
		Cru = Cravg;
	} else
	{
		Crl = Cravg;
		Cru = Crbgt;
	}
	if (Cbbgt <= Cbavg)
	{
		Cbl = Cbbgt;
		Cbu = Cbavg;
	} else
	{
		Cbl = Cbavg;
		Cbu = Cbbgt;
	}
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
		{
			tie(c[0],c[1],c[2]) = c_res(i,j);
			if (Yl <= c[0] && Yu >= c[0] && Crl <= c[2] && Cru >= c[2] && Cbl <= c[1] && Cbu >= c[1])
			{
				ref_count++;
				ref_exist = true;
				tie(r[0],r[1],r[2]) = res(i,j);
				Rw += r[0];
				Gw += r[1];
				Bw += r[2];
			}
		}
	if (!ref_exist)
		return source;
	Rw /= ref_count;
	Gw /= ref_count;
	Bw /= ref_count;
	Yw = 0.299*Rw+0.587*Gw+0.114*Bw;
	Rscale = Yw/Rw;
	Gscale = Yw/Gw;
	Bscale = Yw/Bw;
	Rgwa = Yavg/Ravg;
	Ggwa = Yavg/Gavg;
	Bgwa = Yavg/Bavg;
	if (Bavg+3 >= Gavg && Bavg >= Ravg)
		colorcast = 1;
	else if (Gavg+3 > Ravg && Ravg > Bavg)
		colorcast = 2;
	else if (Ravg > Gavg && Gavg > Bavg)
		colorcast = 3;
	else colorcast = 0;
	if (!colorcast)
		return source;
	if (colorcast == 1)
	{
		R = Rscale;
		G = Gscale;
		B = Bgwa;
	} else
	if (colorcast == 2)
	{
		R = Rscale;
		G = Ggwa;
		B = Bscale;
	} else
	if (colorcast == 3)
	{
		R = Rgwa;
		G = Gscale;
		B = Bscale;
	}
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
		{
			tie(c[0],c[1],c[2]) = res(i,j);
			c[0] = round(c[0]*R);
			c[1] = round(c[1]*G);
			c[2] = round(c[2]*B);
			for (k=0; k<3; k++)
				if (c[k] > 255)
					c[k] = 255;
			res(i,j) = tie(c[0],c[1],c[2]);
		}
	return res;
}

Image mirrorize(Image src_image, int radius)
{
	int rows = src_image.n_rows, cols = src_image.n_cols;
	Image result(rows+radius*2, cols+radius*2);
	int i,j,k;
	for (i=radius; i<rows+radius; i++)
		for (j=radius; j<cols+radius; j++)
			result(i,j) = src_image(i-radius,j-radius);
	for (i=radius; i<rows+radius; i++)
		for (k=0; k<radius; k++)
		{
			result(i,radius-1-k) = result(i,radius+k);
			result(i,cols+radius+k) = result(i,cols+radius-1-k);
		}
	for (j=radius; j<cols+radius; j++)
		for (k=0; k<radius; k++)
		{
			result(radius-1-k,j) = result(radius+k,j);
			result(rows+radius+k,j) = result(rows+radius-1-k,j);
		}
// Now for the corners
	for (i=0; i<radius; i++)
		for (j=0; j<radius; j++)
		{
//top left
			result(radius-1-i,radius-1-j) = result(radius+i,radius+j);
//bottom right
			result(rows+radius+i,cols+radius+j) = result(rows+radius-1-i,cols+radius-1-j);
//bottom left
			result(rows+radius+i,radius-1-j) = result(rows+radius-1-i,radius+j);
//top right
			result(radius-1-i, cols+radius+j) = result(radius+i,cols+radius-1-j);
		}
	return result;
}

Image demirrorize(Image src, int radius)
{
	Image result(src.n_rows-radius*2,src.n_cols-radius*2);
	int i,j,rows = src.n_rows, cols = src.n_cols;
	for (i=radius; i<rows-radius; i++)
		for (j=radius; j<cols-radius; j++)
			result(i-radius,j-radius) = src(i,j);
	return result;
}

class KernelSharpen
{
public:
	tuple<uint, uint, uint> operator () (const Image &m) const
	{
    	int size = 2 * radius + 1;
		int r,g,b, sum_r=0, sum_g=0, sum_b=0;
    	Matrix<int> kernel = {{-1, -4, -1},
    	                         {-4, 26, -4},
    	                         {-1, -4, -1}};
		for (int i=0; i<size; i++)
			for (int j=0; j<size; j++)
			{
				tie(r,g,b) = m(i,j);
				sum_r += r * kernel(i,j);
				sum_g += g * kernel(i,j);
				sum_b += b * kernel(i,j);
			}
		if (sum_r < 0)
			sum_r = 0;
		if (sum_g < 0)
			sum_g = 0;
		if (sum_b < 0)
			sum_b = 0;
		if (sum_r/6 >= 255)
			sum_r = 255;
		else sum_r /= 6;
		if (sum_g/6 >= 255)
			sum_g = 255;
		else sum_g /= 6;
		if (sum_b/6 >= 255)
			sum_b = 255;
		else sum_b /= 6;
		return make_tuple(sum_r,sum_g,sum_b);
/*       uint r, g, b, sum_r = 0, sum_g = 0, sum_b = 0;
        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                // Tie is useful for taking elements from tuple
                tie(r, g, b) = m(i, j);
                sum_r += r;
                sum_g += g;
                sum_b += b;
            }
        }
        auto norm = size * size;
        sum_r /= norm;
        sum_g /= norm;
        sum_b /= norm;
        return make_tuple(sum_r, sum_g, sum_b);*/
    }
    // Radius of neighbourhoud, which is passed to that operator
    static const int radius = 1;
};

class Imgbox
{
	int sqr(int i);
	bool r_ins;
	bool g_ins;
	bool b_ins;
	Image r;
	Image g;
	Image b;
	Image result;
	int only_pos(int i);
	int only_neg(int i);
	void put_color(Image src, uint color);
	float mse_count(int i, int j, Image src);
	float ccor_count(int i, int j, Image src);
public:
	static const uint way;
	Image get_r();
	Image get_g();
	Image get_b();
	Image get_res();
	Imgbox(Image source);
	void put(uint color);
};

int Imgbox::sqr(int i)
{
	return i*i;
}

int Imgbox::only_pos(int i)
{
	if (i > 0)
		return i;
	else return 0;
}

int Imgbox::only_neg(int i)
{
	if (i < 0)
		return i;
	else return 0;
}

const uint Imgbox::way = w_avdev;

Image Imgbox::get_r() {return r;}
Image Imgbox::get_g() {return g;}
Image Imgbox::get_b() {return b;}
Image Imgbox::get_res() {return result;}

Imgbox::Imgbox(Image source): r_ins(false), g_ins(false), b_ins(false), r(source.n_rows/3,source.n_cols), g(source.n_rows/3,source.n_cols), b(source.n_rows/3,source.n_cols), result(source.n_rows/3,source.n_cols)
{
	uint rows = source.n_rows/3, cols = source.n_cols;
	for (uint i=0; i<rows; i++)
		for (uint j=0; j<cols; j++)
		r(i,j) = source(rows*2+i,j);
	
	for (uint i=0; i<rows; i++)
		for (uint j=0; j<cols; j++)
		g(i,j) = source(rows+i,j);

	for (uint i=0; i<rows; i++)
		for (uint j=0; j<cols; j++)
		b(i,j) = source(i,j);
}

void Imgbox::put(uint color)
{
	Image& red=r;
	Image& green=g;
	Image& blue=b;
	if (color == c_r)
	{
		put_color(red,color);
		r_ins = true;
	}
	else if (color == c_g)
	{
		put_color(green,color);
		g_ins = true;
	}
	else if (color == c_b)
	{
		put_color(blue,color);
		b_ins = true;
	}
}

void Imgbox::put_color(Image src, uint color)
{
	uint y[3],z[3];
	int rows=src.n_rows, cols=src.n_cols;
	if (!(r_ins || g_ins || b_ins))
	{
		for (int i=0; i<rows; i++)
			for (int j=0; j<cols; j++)
			{
				tie(y[0],y[1],y[2]) = result(i,j);
				tie(z[0],z[1],z[2]) = src(i,j);
//				result(i,j) = src(i,j);
				y[color] = z[color];
				result(i,j) = tie(y[0],y[1],y[2]);
			}
		return;
	}
	if (way == w_avdev)
	{
		Image& source = src;
		int mini=0, minj=0;
		float min=1048576, res;
		for (int i=-15; i<=15; i++)
			for (int j=-15; j<=15; j++)
			{
				res = mse_count(i,j,src);
//				printf("%f\n",res);
//				if (i>=-4 && i<=0 && j >= -1 && j <= 3)
//					printf("%f,%d,%d\n",res,i,j);
				if (res < min)
				{
					min = res;
					mini = i;
					minj = j;
				}
			}
		rows = src.n_rows-abs(mini); cols = src.n_cols-abs(minj);
		int negi = only_neg(mini), negj = only_neg(minj), posi = only_pos(mini), posj = only_pos(minj);
		for (int k=0; k<rows; k++)
			for (int t=0; t<cols; t++)
			{
				tie(y[0],y[1],y[2]) = result(k+posi,t+posj);
				tie(z[0],z[1],z[2]) = src(k-negi,t-negj);
				y[color] = z[0];
				result(k+posi,t+posj) = tie(y[0],y[1],y[2]);
			}
		return;
	} else
	{
		Image& source = src;
		int maxi=0,maxj=0;
		float max=0, res;
		for (int i=-15; i<=15; i++)
			for (int j=-15; j<=15; j++)
			{
				res=ccor_count(i,j,source);
				if (res >= max)
				{
					max = res;
					maxi = i;
					maxj = j;
				}
			}
		rows = src.n_rows-abs(maxi); cols = src.n_cols-abs(maxj);
		int negi = only_neg(maxi), negj = only_neg(maxj), posi = only_pos(maxi), posj = only_pos(maxj);
		for (int k=0; k<rows; k++)
			for (int t=0; t<cols; t++)
			{
				tie(y[0],y[1],y[2]) = result(k+posi,t+posj);
				tie(z[0],z[1],z[2]) = src(k-negi,t-negj);
				y[color] = z[color];
				result(k+posi,t+posj) = tie(y[0],y[1],y[2]);
			}
	}
}

float Imgbox::mse_count(int i, int j, Image src)
{
	int y[3],z[3];
	int rows = src.n_rows-abs(i), cols = src.n_cols-abs(j);
	long long int sum=0;
	int negi = only_neg(i), negj = only_neg(j), posi = only_pos(i), posj = only_pos(j);
	int color;
	if (r_ins)
		color = c_r;
	else if (g_ins)
		color = c_g;
	else 
		color = c_b;
//	if (g_ins)
//		for (int k=0; k<rows; k++)
//			for (int t=0; t<cols; t++)
//			{
//				tie(y[0],y[1],y[2]) = result(k+posi,t+posj);
//				tie(z[0],z[1],z[2]) = src(k-negi,t-negj);
//				sum = sum + sqr(z[0]+z[1]+z[2]-y[0]-y[1]-y[2]);
//				sum = sum + sqr(z[color]/*+z[color-1])/2*/-y[color]);
//			}
//	else
	for (int k=30; k<rows-30; k++)
		for (int t=30; t<cols-30; t++)
		{
			tie(y[0],y[1],y[2]) = result(k+posi,t+posj);
			tie(z[0],z[1],z[2]) = src(k-negi,t-negj);
//			if (k>=0&&k<10&&t>=0&&t<10)
//				printf("Comparing result[%d,%d]=%d and src[%d,%d]=%d...(%d,%d)\n",k+posi,t+posj,y[color],k-negi,t-negj,z[color],i,j);
			sum += sqr(y[color]-z[color]);
		}
//	return (sum + 0.0)/(src.n_rows*src.n_cols);
	return (sum + 0.0)/(rows*cols);
}

float Imgbox::ccor_count(int i, int j, Image src)
{
	uint y[3],z[3];
	int rows = src.n_rows-abs(i), cols = src.n_cols-abs(j);
	long long int sum=0;
	int negi = only_neg(i), negj = only_neg(j), posi = only_pos(i), posj = only_pos(j);
	for (int k=0; k<rows; k++)
		for (int t=0; t<cols; t++)
		{
			tie(y[0],y[1],y[2]) = result(k+posi,t+posj);
			tie(z[0],z[1],z[2]) = src(k-negi,t-negj);
			sum = sum + z[0]*y[0];
		}
	return sum;
}

Image align(Image srcImage, bool isPostprocessing, std::string postprocessingType, double fraction, bool isMirror, 
            bool isInterp, bool isSubpixel, double subScale)
{
	Imgbox box(srcImage);
	box.put(c_r);
	box.put(c_g);
	box.put(c_b);
	srcImage = box.get_res();
	if (isPostprocessing)
	{
		if (postprocessingType == "--gray-world")
			srcImage = gray_world(srcImage);
		else if (postprocessingType == "--unsharp")
		{
			if (isMirror)
				srcImage = mirrorize(srcImage, 1);
			srcImage = unsharp(srcImage);
			if (isMirror)
				srcImage = demirrorize(srcImage, 1);
		}
		else if (postprocessingType == "--autocontrast")
			srcImage = autocontrast(srcImage, fraction);
		else if (postprocessingType == "--white-balance")
			srcImage = white_balance(srcImage);
	}
    return srcImage;
}



Image sobel_x(Image src_image) {
    Matrix<double> kernel = {{-1, 0, 1},
                             {-2, 0, 2},
                             {-1, 0, 1}};
    return custom(src_image, kernel);
}

Image sobel_y(Image src_image) {
    Matrix<double> kernel = {{ 1,  2,  1},
                             { 0,  0,  0},
                             {-1, -2, -1}};
    return custom(src_image, kernel);
}

Image unsharp(Image src_image) {
	src_image = src_image.unary_map(KernelSharpen());
    return src_image;
}

Image gray_world(Image src_image) {
	long int sum_r=0,sum_g=0,sum_b=0,sum;
	int r,g,b;
	int rows=src_image.n_rows, cols=src_image.n_cols;
	double c_r,c_g,c_b;
	for (int i=0; i<rows; i++)
		for (int j=0; j<cols; j++)
		{
			tie(r,g,b) = src_image(i,j);
			sum_r += r;
			sum_g += g;
			sum_b += b;
		}
	sum = sum_r + sum_g + sum_b;
	c_r = (sum+0.0)/(3*sum_r);
	c_g = (sum+0.0)/(3*sum_g);
	c_b = (sum+0.0)/(3*sum_b);
	for (int i=0; i<rows; i++)
		for (int j=0; j<cols; j++)
		{
			tie(r,g,b) = src_image(i,j);
			if (r*c_r >= 255)
				r = 255;
			else r = round(r*c_r);
			if (g*c_g >= 255)
				g = 255;
			else g = round(g*c_g);
			if (b*c_b >= 255)
				b = 255;
			else b = round(b*c_b);
			src_image(i,j) = tie(r,g,b);
		}
    return src_image;
}

Image resize(Image src_image, double scale) {
    return src_image;
}

Image custom(Image src_image, Matrix<double> kernel) {
    // Function custom is useful for making concrete linear filtrations
    // like gaussian or sobel. So, we assume that you implement custom
    // and then implement other filtrations using this function.
    // sobel_x and sobel_y are given as an example.
    return src_image;
}

Image autocontrast(Image src_image, double fraction) {
	long int br[256];
	int r,g,b,s,i,j;
	int rows = src_image.n_rows, cols = src_image.n_cols,mini,maxi;
	float y, max=0, min=1024;
	for (int t=0; t<256; t++)
		br[t] = 0;
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
		{
			tie (r,g,b) = src_image(i,j);
			y = r*0.2125+g*0.7154+b*0.0721;
			(br[static_cast<int>(round(y))])++;
		}
	s = rows*cols*fraction;
	i = 0;
	while (i < 256 && s > 0)
	{
		s -= br[i];
		if (s > 0)
			i++;
	}
	i = 0;
	while (i < 256)
		if (br[i++] > 0)
			break;
		else i++;
	mini = i;
	i = 255;
	s = rows*cols*fraction;
	while (i > 0 && s > 0)
	{
		s -= br[i];
		if (s > 0)
			i++;
	}
	i = 255;
	while (i > 0)
		if (br[i] > 0)
			break;
		else i--;
	maxi = i;
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
		{
			tie(r,g,b) = src_image(i,j);
			y = r*0.2125+g*0.7154+b*0.0721;
			if (round(y) > mini && round(y) < maxi)
			{
				if (y > max)
					max = y;
				else if (y < min)
					min = y;
			}
		}
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
		{
			tie (r,g,b) = src_image(i,j);
			y = r*0.2125+g*0.7154+b*0.0721;
			r = round(255*(r-min)+0.0)/(max-min);
			if (r >= 255)
				r = 255;
			else if (r < 0)
				r = 0;
			g = round(255*(g-min)+0.0)/(max-min);
			if (g >= 255)
				g = 255;
			else if (g < 0)
				g = 0;
			b = round(255*(b-min)+0.0)/(max-min);
			if (b >= 255)
				b = 255;
			else if (b < 0)
				b = 0;
			src_image(i,j) = tie(r,g,b);
		}
    return src_image;
}

Image gaussian(Image src_image, double sigma, int radius)  {
    return src_image;
}

Image gaussian_separable(Image src_image, double sigma, int radius) {
    return src_image;
}

Image median(Image src_image, int radius) {
	int hist[3][256];
	int rows=src_image.n_rows, cols=src_image.n_cols;
	int i,j,a,b,k,c[3],sum_c,sum;
	Image out_image(src_image.n_rows,src_image.n_cols);
	sum_c = (2*radius+1)*(2*radius+1)/2/*+(2*radius+1)*(2*radius+1)%2*/;
	printf("%d",sum_c);
	for (i=radius; i<rows-radius; i++)
		for (j=radius; j<cols-radius; j++)
		{
			for (a=0; a<3; a++)
				for (b=0; b < 256; b++)
				{
					hist[a][b] = 0;
				}
			for (a = i-radius; a <= i+radius; a++)
				for (b = j-radius; b <= j+radius; b++)
				{
					tie(c[0],c[1],c[2]) = src_image(a,b);
					for (k=0; k<3; k++)
						hist[k][c[k]]++;
				}
			for (a=0; a<3; a++)
			{
				b = 0;
				sum = sum_c;
				while (1)
				{
					sum -= hist[a][b];
					if (sum > 0)
						b++;
					else break;
				}
				c[a] = b;
			}
			out_image(i,j) = tie(c[0],c[1],c[2]);
		}
    return out_image;
}

Image median_linear(Image src_image, int radius) {
	int hist[3][256];
	int rows=src_image.n_rows, cols=src_image.n_cols;
	int i,j,k,c[3],a,b,sum_c,sum;
	Image out_image(src_image.n_rows,src_image.n_cols);
	sum_c = (2*radius+1)*(2*radius+1)/2/*+(2*radius+1)*(2*radius+1)%2*/;
	for (i=radius; i<rows-radius; i++)
	{
		for (k=0; k<3; k++)
			for (j=0; j<256; j++)
				hist[k][j] = 0;
		j=radius;
		for (a = i-radius; a <= i+radius; a++)
			for (b=j-radius; b <= j+radius; b++)
			{
				tie(c[0],c[1],c[2]) = src_image(a,b);
				for (k=0; k<3; k++)
					hist[k][c[k]]++;
			}
		for (a=0; a<3; a++)
		{
			b = 0;
			sum = sum_c;
			while (1)
			{
				sum -= hist[a][b];
				if (sum > 0)
					b++;
				else break;
			}
			c[a] = b;
		}
		out_image(i,j) = tie(c[0],c[1],c[2]);
		for (j=radius+1; j<cols-radius; j++)
		{
			for (k=-radius; k<=radius; k++)
			{
				tie(c[0],c[1],c[2]) = src_image(i+k,j-radius-1);
				for (a=0; a<3; a++)
					hist[a][c[a]]--;
				tie(c[0],c[1],c[2]) = src_image(i+k,j+radius);
				for (a=0; a<3; a++)
					hist[a][c[a]]++;
			}
			for (a=0; a<3; a++)
			{
				b = 0;
				sum = sum_c;
				while (1)
				{
					sum -= hist[a][b];
					if (sum > 0)
						b++;
					else break;
				}
				c[a] = b;
			}
			out_image(i,j) = tie(c[0],c[1],c[2]);
		}
	}
    return out_image;
}

Image median_const(Image src_image, int radius) {
    return src_image;
}

Image canny(Image src_image, int threshold1, int threshold2) {
    return src_image;
}
