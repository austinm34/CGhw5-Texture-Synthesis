#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "lodepng.h"


struct samp_wind {
	int x;
	int y;
	double dist;
};

std::vector<std::vector<double>> getwp(std::vector<std::vector<double>> &image, int x, int y, int u, int v);	// all of this will eventually have another parameter for window size
std::vector<std::vector<double>> getkp(std::vector<std::vector<double>> &wp);
std::vector<std::vector<double>> getgauss();
std::vector<std::vector<double>> getws(std::vector<std::vector<double>> &sample, int x, int y);
std::vector<std::vector<double>> mult(std::vector<std::vector<double>> a1, std::vector<std::vector<double>> a2, int s);
std::vector<std::vector<double>> sub(std::vector<std::vector<double>> a1, std::vector<std::vector<double>> a2, int s);
double vec_sum(std::vector<std::vector<double>> a1, int s);
double dist(std::vector<std::vector<double>> &wp, std::vector<std::vector<double>> &kp, std::vector<std::vector<double>> &ws, int s);
std::vector<samp_wind> min_3(std::vector<samp_wind>a1, int s);

void synthesize_pixel(std::vector<std::vector<double>> &syn_image, std::vector<std::vector<double>> &samp_image, int x, int y, int u, int v, int wsize);



int main() {

	std::vector<unsigned char> image; // raw pixels
	unsigned width, height;

	// ================================= decode image ==================================
	unsigned error = lodepng::decode(image, width, height, "small_sand.png");
	int length = width * height * 4;

	if (error) {
		std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
		exit(1);
	}
	
	// ======================================= gray scale image ==============================
	for (int i = 0; i < length; i++) {
		if (i % 4 == 0) {
			int gray_scale = (0.3 * +image[i]) + (0.59 * +image[i + 1]) + (0.11 * +image[i + 2]);
			image[i] = image[i + 1] = image[i + 2] = (unsigned char)gray_scale;
		}
	}
	// ====================================== prepare for texture synthesis ==============================
	// put image into an array
	int u = (int)width;
	int v = (int)height;
	std::cout << "width " << u << " height " << v << std::endl;
	std::vector< std::vector<double> > sample_image(u, std::vector<double>(v, 0));
	int c1 = 0;	// a counter variable to match to image vector;
	for (int i = 0; i < u; i++) {
		for (int j = 0; j < v; j++) {
			sample_image[i][j] = (double)image[c1 * 4];
			c1++;
		}
	}

	// ===================================== synthesis start ===============================================

	// this is our array that will hold the synthesized texture. it has been initiallized to 0s.

	std::vector<std::vector<double>> synthesized_image(u, std::vector<double>(v, 0));
	// lets start with placing a seed in the center. the seed will be 3x3
	//*
	int middle_x = u / 2;
	int middle_y = v / 2;
	synthesized_image[middle_x - 1][middle_y - 1] = sample_image[middle_x - 1][middle_y - 1];	//x-1 y-1
	synthesized_image[middle_x][middle_y - 1] = sample_image[middle_x][middle_y - 1]; //x y-1
	synthesized_image[middle_x + 1][middle_y - 1] = sample_image[middle_x + 1][middle_y - 1]; // x+1 y-1

	synthesized_image[middle_x - 1][middle_y] = sample_image[middle_x - 1][middle_y]; //x-1 y
	synthesized_image[middle_x][middle_y] = sample_image[middle_x][middle_y]; // x y
	synthesized_image[middle_x + 1][middle_y] = sample_image[middle_x + 1][middle_y]; //x+1 y

	synthesized_image[middle_x - 1][middle_y + 1] = sample_image[middle_x - 1][middle_y + 1]; // x-1 y+1
	synthesized_image[middle_x][middle_y + 1] = sample_image[middle_x][middle_y + 1]; // x y+1
	synthesized_image[middle_x + 1][middle_y + 1] = sample_image[middle_x + 1][middle_y + 1]; // x+1 y+1

	// =============================== iterate around seed ===================
	
	// to do this we save four values corelating to the bounds of a square
	int x1 = middle_x - 2;	// lower bound on x
	int x2 = middle_x + 2;	// upper bound on x
	int y1 = middle_y - 2;	// lower bound on y
	int y2 = middle_y + 2;	// upper bound on y
	//*
	bool x1_good, x2_good, y1_good, y2_good;
	x1_good = x2_good = y1_good = y2_good = true;

	// x1_good && x2_good && y1_good && y2_good
	// bool test = true;
	while (x1_good || x2_good || y1_good || y2_good) {			// possible issue, this calculates the corners twice each time the loop goes around.

		if (x1 < 10) {
			std::cout << "x1: " << x1 << " x2: " << x2 << " y1: " << y1 << " y2: " << y2 << std::endl;
		}

		for (int j = y1; j <= y2; j++) {			// top
			//synthesized_image[x1][j] = 255;
			synthesize_pixel(synthesized_image, sample_image, x1, j, u, v, 3);
		}	
		std::cout << "top done" << std::endl;
		for (int i = x1; i <= x2; i++) {			// right
			//synthesized_image[i][y2] = 255;
			synthesize_pixel(synthesized_image, sample_image, i, y2, u, v, 3);
		}
		std::cout << "right done" << std::endl;
		for (int j = y2; j > y1; j--) {			// bottom	I don't know why only a > works, it shouldn't but it does
			//synthesized_image[x2][j] = 255;
			synthesize_pixel(synthesized_image, sample_image, x2, j, u, v, 3);
		}
		std::cout << "bottom done" << std::endl;
		for (int i = x2; i > x1; i--) {			// left
			//synthesized_image[i][y1] = 255;
			synthesize_pixel(synthesized_image, sample_image, i, y1, u, v, 3);
		}

		// check that we are in bounds
		if (x1 == 0)
			x1_good = false;
		if (x2 == u-1)
			x2_good = false;
		if (y1 == 0)
			y1_good = false;
		if (y2 == v-1)
			y2_good = false;

		// if we're in bounds, iterate the coords of the rectangle
		if (x1_good)
			x1--;
		if (x2_good)
			x2++;
		if (y1_good)
			y1--;
		if (y2_good)
			y2++;
	}
	//*/
	// ========================================== test, lets calculate one pixel =====================================
	/*
	int wsize = 3;
	int wind_wid = 3 / 2;
	// were going to synthesize at point x2, middle_y

	// step one: get wp
	std::vector<std::vector<double>> wp = getwp(synthesized_image, x2, middle_y, u, v);
	// step two: get kp
	std::vector<std::vector<double>> kp = getkp(wp);
	// step three: iterate through sample and get the top 3 matches
	std::vector<samp_wind> sample_windows;
	std::vector<std::vector<double>> ws;
	samp_wind temp_wind;
	for (int i = wind_wid; i < u - wind_wid; i++) {
		for (int j = wind_wid; j < v - wind_wid; j++) {
		
			ws = getws(sample_image, i, j);
			temp_wind.x = i;
			temp_wind.y = j;
			temp_wind.dist = dist(wp, kp, ws, wsize);
			sample_windows.push_back(temp_wind);

		}
	}
	std::vector<samp_wind> top3 = min_3(sample_windows, wsize);
	int win_ind = rand() % 3;
	synthesized_image[x2][middle_y] = sample_image[top3[win_ind].x][top3[win_ind].y];

	// we've synthesized 1 pixel !
	*/
	// ========================================= synthesis end, prepare to encode ===================================

	//put image array back into vector to be encoded. we'll use the original one we decoded so we don't
	//have to make a new one. we'll assume we're synthesizing an image of the same size as a sample, using only a 3x3 seed to start.
	int c2 = 0;	// a counter variable
	int ind = 0; // stores the index of interest in the output image vector
	for (int i = 0; i < u; i++) {
		for (int j = 0; j < v; j++) {
			ind = c2 * 4;
			image[ind] = image[ind + 1] = image[ind + 2] = (unsigned char)synthesized_image[i][j];
			c2++;
		}
	}

	//================================= encode image =====================================
	error = lodepng::encode("sand_out.png", image, width, height);

	if (error) {
		std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
		exit(1);
	}

	std::cout << "Press any key to continue...";
	std::getchar();
	return 0;
}


std::vector<std::vector<double>> getwp(std::vector<std::vector<double>> &image, int x, int y, int u, int v) {

	// image is an array storing the image to be synthesized
	// x and y are the center pixel we want to generate a window around
	// u and v are the width and height of the image so we can check bounds

	std::vector<std::vector<double>> wp(3, std::vector<double>(3, -1));	// both "3s" will eventually be the window size

	int inval_l, inval_r, inval_t, inval_b;
	inval_l = inval_r = inval_t = inval_b = 0;

	if (x - 1 < 0) {	// eventually 1 will be floor(size/2)
		inval_t = 1;	// this will eventually be how many pixels we need to not update in the mask, aka abs(x-1)
	}					// for now we're using the simple case of a 3x3 window
	if (x + 1 > u) {
		inval_b = 1;
	}
	if (y - 1 < 0) {
		inval_l = 1;
	}
	if (y + 1 > u) {
		inval_r = 1;
	}

	for (int i = inval_t; i < 3 - inval_b; i++) {			// "3" will eventually be the window size
		for (int j = inval_l; j < 3 - inval_r; j++) {
			wp[i][j] = image[i][j];
		}
	}

	return wp;

}

std::vector<std::vector<double>> getkp(std::vector<std::vector<double>> &wp) {	// gets "known" vector for the pixel, or kp.
	std::vector<std::vector<double>> kp(3, std::vector<double>(3, 0));
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (wp[i][j] != -1) {
				kp[i][j] = 1;		// if not -1, then we know it. therefore set known = 1 for that position.
			}
		}
	}

	return kp;
}

std::vector<std::vector<double>> getgauss() {	// this will eventually take size and sigma
	std::vector<std::vector<double>> gauss(3, std::vector<double>(3, 0));

	double sig = 0;
	double s = 2.0 * sig * sig;

	int width = 1; // this will eventually be size/2

	for (int i = 0 - width; i < width; i++) {
		for (int j = 0 - width; j < width; j++) {
			gauss[i + width][j + width] = (exp(-(i*i + j*j) / s)) / (3.14159265 * s);
		}
	}

	//we don't need to normalize, it will be done later.

	return gauss;
}

std::vector<std::vector<double>> getws(std::vector<std::vector<double>> &sample, int x, int y) {
	std::vector<std::vector<double>> ws(3, std::vector<double>(3, 0));

	// this function assumes that x and y are the center pixel, and are chosen in such away that edge
	// cases will not arise for the given window size.
	for (int i = 0; i < 3; i++) {		// 3 will be replaced by window size
		for (int j = 0; j < 3; j++) {
			ws[i][j] = sample[x - 1 + i][y - 1 + j];	// the ones will be replaced by a width = floor(size/2)
		}
	}

	return ws;
}

std::vector<std::vector<double>> mult(std::vector<std::vector<double>> a1, std::vector<std::vector<double>> a2, int s) {	// element by element multiplication of square matrices
	std::vector<std::vector<double>> prod(s, std::vector<double>(s, 0));

	for (int i = 0; i < s; i++) {
		for (int j = 0; j < s; j++) {
			prod[i][j] = a1[i][j] * a2[i][j];
		}
	}

	return prod;
}

std::vector<std::vector<double>> sub(std::vector<std::vector<double>> a1, std::vector<std::vector<double>> a2, int s) {	// element by element subtraction of square matrices
	std::vector<std::vector<double>> diff(s, std::vector<double>(s, 0));

	for (int i = 0; i < s; i++) {
		for (int j = 0; j < s; j++) {
			diff[i][j] = a1[i][j] - a2[i][j];
		}
	}

	return diff;
}

double vec_sum(std::vector<std::vector<double>> a1, int s) {
	double sum = 0;

	for (int i = 0; i < s; i++) {
		for (int j = 0; j < s; j++) {
			sum += a1[i][j];
		}
	}

	return sum;
}

double dist(std::vector<std::vector<double>> &wp, std::vector<std::vector<double>> &kp, std::vector<std::vector<double>> &ws, int s) {
	std::vector<std::vector<double>> dif1 = sub(ws, wp, s);			// difference between wp and sampled window
	std::vector<std::vector<double>> dif_known = mult(dif1, kp, s);	// difference only where wp is known
	std::vector<std::vector<double>> sq_dif = mult(dif_known, dif_known, s);	// squared differences

	std::vector<std::vector<double>> gaus = getgauss();
	std::vector<std::vector<double>> gaus_known = mult(gaus, kp, s);	// gaussian kernel, only where known
	double gaus_sum = vec_sum(gaus_known, s);	// sum of known gaussian, used for normalization

	std::vector<std::vector<double>> weighted_sq_dif = mult(gaus_known, sq_dif, s);	// weighted squared differences
	double sum_sq_dif = vec_sum(weighted_sq_dif, s);	// sum of weighted squared differences

	double norm_ssd = sum_sq_dif / gaus_sum;	// normalized sum of weighted squared differences

	return norm_ssd;

}

std::vector<samp_wind> min_3(std::vector<samp_wind> a1, int s) {	// returns the sample windows with the minimum distances
	std::vector<samp_wind> top3;
	
	double min = a1[0].dist;
	int min_ind = 0;
	for (int n = 0; n < 3; n++) {
		for (int i = 0; i < s; i++) {
			if (a1[i].dist < min) {
				min = a1[i].dist;
				min_ind = i;
			}
		}
		top3.push_back(a1[min_ind]);
		a1.erase(a1.begin() + min_ind);
	}

	return top3;
}

void synthesize_pixel(std::vector<std::vector<double>> &syn_image, std::vector<std::vector<double>> &samp_image, int x, int y, int u, int v, int wsize) {
	int wind_wid = wsize / 2;
	// were going to synthesize at point x2, middle_y

	// step one: get wp
	std::vector<std::vector<double>> wp = getwp(syn_image, x, y, u, v);
	// step two: get kp
	std::vector<std::vector<double>> kp = getkp(wp);
	// step three: iterate through sample and get the top 3 matches
	std::vector<samp_wind> sample_windows;
	std::vector<std::vector<double>> ws;
	samp_wind temp_wind;
	for (int i = wind_wid; i < u - wind_wid; i++) {
		for (int j = wind_wid; j < v - wind_wid; j++) {

			ws = getws(samp_image, i, j);
			temp_wind.x = i;
			temp_wind.y = j;
			temp_wind.dist = dist(wp, kp, ws, wsize);
			sample_windows.push_back(temp_wind);

		}
	}
	std::vector<samp_wind> top3 = min_3(sample_windows, wsize);
	int win_ind = rand() % 3;
	syn_image[x][y] = samp_image[top3[win_ind].x][top3[win_ind].y];
}