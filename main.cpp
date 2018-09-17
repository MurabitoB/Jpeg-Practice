#include<iostream>
#include<fstream>
#include<string>
#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdio.h>
#include<iomanip>
#include<map>
#include"jpgtable.h"
#define SQH 0.707106781186547  /* square root of 2 */
#define SWAP(a,b)  tempr=(a); (a) = (b); (b) = tempr
#define MSB(a) a < 0 ? 1:0
#define PI  3.1415927
#define PRINT_BLOCK 63
#define PRINT_BLOCK_COL 63
bool dc_is_decoding = true;
int decode_row = 0, decode_col = 0;
using namespace std;
unsigned char output_counter = 0;
unsigned char output = 0;
int diff_pre_value = 0;
int decode_codeword = 0;
bool dc_codeword_decoding = true, dc_diff_decoding = false, ac_codeword_decoding = false, ac_diff_decoding = false, complete_decoding = false;
int decode_ssss_index = 0;
int decode_counter = 0;
int decode_run = 0;
int decode_size = 0;
int decode_ac_counter = 1;
map < int, map <int, ac_table_de>> dec_ac_t;//AC table for decode
FILE *file_out;
void output_method(char output_value)
{
	output <<= 1;
	output += output_value;
	output_counter++;
	if (output_counter == 8)
	{
		fwrite(&output, sizeof(char), 1, file_out);
		output_counter = 0;
		output = 0;
	}
}
struct block
{
	int value[8][8];
};
static block input_block[64][64];


static void fft1(float *data, int nn, int isign)
{
	int n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	float tempr, tempi;
	n = nn << 1;
	j = 1;
	for (i = 1; i<n; i += 2) {
		if (j>i) {
			SWAP(data[j], data[i]);
			SWAP(data[j + 1], data[i + 1]);
		}
		m = n >> 1;
		while (m >= 2 && j>m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax = 2;
	while (n>mmax) {
		istep = 2 * mmax;
		theta = 6.28318530717959 / (-isign*mmax);
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m<mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempr = wr*data[j] - wi*data[j + 1];
				tempi = wr*data[j + 1] + wi*data[j];
				data[j] = data[i] - tempr;
				data[j + 1] = data[i + 1] - tempi;
				data[i] += tempr;
				data[i + 1] += tempi;
			}
			wr = (wtemp = wr)*wpr - wi*wpi + wr;
			wi = wi*wpr + wtemp*wpi + wi;
		}
		mmax = istep;
	}

	if (isign == -1) {
		for (i = 1; i <= n; ++i)
			data[i] = data[i] / nn;
	}
}


void dct1(float *x, int n)
{
	int i, ii, nn, mm;
	float tc, ts, sqn, temp2;
	double temp1;
	float *v;


	nn = n >> 1;
	mm = n << 1;
	sqn = (float)sqrt((double)n);

	v = (float *)calloc(mm, sizeof(float));
	if (v == NULL) {
		printf("allocation failure\n");
		exit(1);
	}

	for (i = 0; i<nn; i++) {
		ii = i << 1;
		v[ii] = x[ii];
		v[ii + 1] = 0.0;
	}
	for (i = nn; i<n; i++) {
		ii = i << 1;
		v[ii] = x[mm - ii - 1];
		v[ii + 1] = 0.0;
	}

	fft1(v - 1, n, 1);

	temp2 = SQH / sqn;
	x[0] = v[0] / sqn;
	for (i = 1; i <= nn; i++) {
		ii = i << 1;
		temp1 = (double)(PI*i / mm);
		tc = (float)cos(temp1);
		ts = (float)sin(temp1);
		x[i] = 2.0*(tc*v[ii] + ts*v[ii + 1])*temp2;
		x[n - i] = 2.0*(ts*v[ii] - tc*v[ii + 1])*temp2;
	}

	free(v);
}
/* -------------------------------------------------- */

void idct1(float *x, int n)
{
	int i, ii, mm, nn;
	float *v;
	float temp2, tc, ts, sqn;
	double temp1;
	nn = n >> 1;
	mm = n << 1;
	sqn = (float)sqrt((double)n);

	v = (float *)calloc(mm, sizeof(float));
	if (v == NULL) {
		printf("allocation failure\n");
		exit(1);
	}

	temp2 = sqn / SQH;
	v[0] = x[0] * sqn;
	v[1] = 0.0;
	for (i = 1; i<n; i++) {
		ii = i << 1;
		temp1 = (double)(PI*i / mm);
		tc = (float)cos(temp1);
		ts = (float)sin(temp1);
		v[ii] = 0.5*(tc*x[i] + ts*x[n - i])*temp2;
		v[ii + 1] = 0.5*(ts*x[i] - tc*x[n - i])*temp2;
	}

	fft1(v - 1, n, -1);

	for (i = 0; i<nn; i++) {
		ii = i << 1;
		x[ii] = v[ii];
	}
	for (i = nn; i<n; i++) {
		ii = i << 1;
		x[mm - ii - 1] = v[ii];
	}
	free(v);
}
void readblock(FILE *file_i)
{
	unsigned char readbyte;
	for (int i = 0; i < 64; i++)
	{
		for (int l = 0; l < 8; l++)
		{
			for (int j = 0; j < 64; j++)
			{
				for (int k = 0; k < 8; k++)
				{
					fread(&readbyte, sizeof(char), 1, file_i);
					input_block[i][j].value[l][k] = readbyte;
				}
			}
		}
	}
}
void writeblock()
{
	unsigned char writebyte;
	for (int i = 0; i < 64; i++)
	{
		for (int l = 0; l < 8; l++)
		{
			for (int j = 0; j < 64; j++)
			{
				for (int k = 0; k < 8; k++)
				{
					writebyte = input_block[i][j].value[l][k];
					fwrite(&writebyte, sizeof(char), 1, file_out);
				}
			}
		}
	}
}
void dct2(float **x, int n)
{
	int i, j;
	float *y;
	y = (float *)calloc(n, sizeof(float));
	if (y == NULL) {
		printf("allocation failure\n");
		exit(1);
	}
	for (i = 0; i<n; i++) {
		for (j = 0; j<n; j++)
			y[j] = x[j][i];
		dct1(y, n);
		for (j = 0; j<n; j++)
			x[j][i] = y[j];
	}   /* end of loop i */

	for (i = 0; i<n; i++) {
		for (j = 0; j<n; j++)
			y[j] = x[i][j];
		dct1(y, n);
		for (j = 0; j<n; j++)
			x[i][j] = y[j];
	}   /* end of loop i */

	free(y);
}

/* ----------------------------------------------- */

void idct2(float **x, int n)
{
	int i, j;
	float *y;
	y = (float *)calloc(n, sizeof(float));
	if (y == NULL) {
		printf("allocation failure\n");
		exit(1);
	}
	for (i = 0; i<n; i++) {
		for (j = 0; j<n; j++)
			y[j] = x[j][i];
		idct1(y, n);
		for (j = 0; j<n; j++)
			x[j][i] = y[j];
	}   /* end of loop i */

	for (i = 0; i<n; i++) {
		for (j = 0; j<n; j++)
			y[j] = x[i][j];
		idct1(y, n);
		for (j = 0; j<n; j++)
			x[i][j] = y[j];
	}   /* end of loop i */

	free(y);
}

void DCT()
{
	float **dct_block = new float *[8];
	for (int i = 0; i < 8; i++){ dct_block[i] = new float[8]; }
	auto get_ptr = [](char a, char b, float **ptr)
	{
		for (int i = 0; i < 8; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				ptr[i][j] = input_block[a][b].value[i][j] - 128;
			}
		}
		return ptr;
	};
	auto return_ptr = [](char a, char b, float **ptr){
		for (int i = 0; i < 8; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				input_block[a][b].value[i][j] = round(ptr[i][j]);
			}
		}
		return ptr; };
	for (int i = 0; i < 64; i++)
	{
		for (int j = 0; j < 64; j++)
		{
			dct2(get_ptr(i, j, dct_block), 8);
			return_ptr(i, j, dct_block);
		}
	}
	for (int i = 0; i < 8; i++){ delete dct_block[i]; }
	delete dct_block;
}
void I_DCT()
{
	float **dct_block = new float *[8];
	for (int i = 0; i < 8; i++){ dct_block[i] = new float[8]; }
	auto get_ptr = [](char a, char b, float **ptr)
	{
		for (int i = 0; i < 8; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				ptr[i][j] = input_block[a][b].value[i][j];
			}
		}
		return ptr;
	};
	auto return_ptr = [](char a, char b, float **ptr){
		for (int i = 0; i < 8; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				input_block[a][b].value[i][j] = round(ptr[i][j]) + 128;
			}
		}
		return ptr; };
	for (int i = 0; i < 64; i++)
	{
		for (int j = 0; j < 64; j++)
		{
			idct2(get_ptr(i, j, dct_block), 8);
			return_ptr(i, j, dct_block);
		}
	}
	for (int i = 0; i < 8; i++){ delete dct_block[i]; }
	delete dct_block;
}
void quantization(int QF)
{
	float factor;
	float q_block[8][8];
	if (QF >= 50)
	{
		factor = 200 - 2 * QF;
	}
	else
	{
		factor = 5000 / QF;
	}
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			q_block[i][j] = quantization_table[i][j] * factor / 100;
		}
	}
	for (int i = 0; i < 64; i++)
	{
		for (int j = 0; j < 64; j++)
		{
			for (int k = 0; k < 8; k++)
			{
				for (int l = 0; l < 8; l++)
				{
					input_block[i][j].value[k][l] = round((float)input_block[i][j].value[k][l] / q_block[k][l]);
				}
			}
		}
	}
}
void de_quantization(int QF)
{
	float factor;
	float q_block[8][8];
	if (QF >= 50)
	{
		factor = 200 - 2 * QF;
	}
	else
	{
		factor = 5000 / QF;
	}
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			q_block[i][j] = quantization_table[i][j] * factor / 100;
		}
	}
	for (int i = 0; i < 64; i++)
	{
		for (int j = 0; j < 64; j++)
		{
			for (int k = 0; k < 8; k++)
			{
				for (int l = 0; l < 8; l++)
				{
					input_block[i][j].value[k][l] = round((float)input_block[i][j].value[k][l] * q_block[k][l]);
				}
			}
		}
	}
}
void dc_encode(unsigned char row, unsigned char col)
{
	unsigned char ssss_index = 0;
	int output_codeword = 0;
	int lower_bound = 0;
	int diff_val_codeword;
	int diff_now_value = input_block[row][col].value[0][0] - diff_pre_value;
	auto get_diff_codeword = [](unsigned char ssss, int diff_value, int lb)
	{
		int diff_size = powl(2, ssss);
		int ret_val;
		if (diff_value > 0)
		{
			ret_val = (diff_value - lb) + diff_size / 2;
		}
		else
		{
			ret_val = lb + diff_value + diff_size / 2 - 1;
		}
		ret_val <<= (32 - ssss);
		return ret_val;
	};
	for (int i = 0; i < 12; i++)
	{
		if (pow(2, i) > abs(diff_now_value))
		{
			ssss_index = i;
			if (i > 0)
			{
				lower_bound = powl(2, i - 1);
			}
			break;
		}
	}
	/*FIX CODEWORD*/
	output_codeword = dc_t[ssss_index].codeword << (32 - dc_t[ssss_index].code_length);
	/*OUT PUT CODEWORD*/
	for (int i = 0; i < dc_t[ssss_index].code_length; i++)
	{
		output_method(MSB(output_codeword));
		output_codeword <<= 1;
	}
	/*GET DIFF VALUE CODEWORD*/
	diff_val_codeword = get_diff_codeword(ssss_index, diff_now_value, lower_bound);
	/*OUTPUT DIFF VALUE CODEWORD*/
	if (ssss_index == 0)
	{
		output_method(0);
	}
	for (int i = 0; i < ssss_index; i++)
	{
		output_method(MSB(diff_val_codeword));
		diff_val_codeword <<= 1;
	}
	diff_pre_value = input_block[row][col].value[0][0];
}
void ac_encode(unsigned char row, unsigned char col)
{
	unsigned char run = 0;
	unsigned char size = 0;
	int output_codeword = 0;
	int lower_bound;
	int diff_val_codeword;
	auto get_diff_codeword = [](unsigned char ssss, int diff_value, int lb)
	{
		int diff_size = pow(2, ssss);
		int ret_val;
		if (diff_value > 0)
		{
			ret_val = (diff_value - lb) + diff_size / 2;
		}
		else
		{
			ret_val = lb + diff_value + diff_size / 2 - 1;
		}
		ret_val <<= (32 - ssss);
		return ret_val;
	};
	for (int i = 1; i < 64; i++)
	{
		int block_inner_row_index, block_inner_col_index;
		/*set Z-ZAG */
		block_inner_row_index = zigZagOrder[i] / 8;
		block_inner_col_index = zigZagOrder[i] % 8;

		if (input_block[row][col].value[block_inner_row_index][block_inner_col_index] == 0)
		{
			run++;
		}
		/*if encount nonzero*/
		else
		{
			while (run > 15)
			{
				/*output ZRL*/
				int zrl_codeword = ac_t[15][0].codeword << (32 - ac_t[15][0].length);
				for (int i = 0; i < ac_t[15][0].length; i++)
				{
					output_method(MSB(zrl_codeword));
					zrl_codeword <<= 1;
				}
				/*run -15*/
				run -= 15;
			}
			/*compute size*/
			for (int i = 0; i < 12; i++)
			{
				if (pow(2, i) > abs(input_block[row][col].value[block_inner_row_index][block_inner_col_index]))
				{
					size = i;
					if (i > 0)
					{
						lower_bound = pow(2, i - 1);
					}
					break;
				}
			}
			/*output codeword*/
			output_codeword = ac_t[run][size].codeword << (32 - ac_t[run][size].length);
			for (int i = 0; i < ac_t[run][size].length; i++)
			{
				output_method(MSB(output_codeword));
				output_codeword <<= 1;
			}
			/*get diffvalue codeword*/
			diff_val_codeword = get_diff_codeword(size, input_block[row][col].value[block_inner_row_index][block_inner_col_index], lower_bound);
			/*output diffvalue codeword*/
			for (int i = 0; i < size; i++)
			{
				output_method(MSB(diff_val_codeword));
				diff_val_codeword <<= 1;
			}
			/*run reset*/
			run = 0;
		}
	}
	/*EOB*/
	int eob_codeword = ac_t[0][0].codeword << (32 - ac_t[0][0].length);
	for (int i = 0; i < ac_t[0][0].length; i++)
	{
		output_method(MSB(eob_codeword));
		eob_codeword <<= 1;
	}

}
void print_block(int row, int col)
{
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			cout << setw(4) << input_block[row][col].value[i][j];
		}
		cout << endl;
	}
}
void set_decode_actable()
{
	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			if (ac_t[i][j].length != 0)
			{
				dec_ac_t[(ac_t[i][j].codeword)][ac_t[i][j].length] = { i, j };
			}
		}
	}
}
void set_decode_flag(bool dc1, bool dc2, bool ac1, bool ac2)
{
	dc_codeword_decoding = dc1;
	dc_diff_decoding = dc2;
	ac_codeword_decoding = ac1;
	ac_diff_decoding = ac2;
	decode_counter = 0;
	decode_codeword = 0;
}
void dc_decode()
{
	if (dc_codeword_decoding)
	{
		if (decode_counter >= 2)
		{
			for (int i = 0; i < 12; i++)
			{
				if (decode_counter == dc_t[i].code_length)
				{
					if ((decode_codeword << (32 - decode_counter)) == (dc_t[i].codeword << (32 - dc_t[i].code_length)))
					{
						decode_ssss_index = i;
						set_decode_flag(false, true, false, false);
					}
				}
			}
		}
	}
	else if (dc_diff_decoding)
	{
		int diff_now_value = 0;
		if (decode_counter == decode_ssss_index)
		{
			/*output previous value + now value(dpcm)*/
			if (decode_codeword < powl(2, decode_ssss_index) / 2)// if decodeword is negative
			{
				diff_now_value = -(powl(2, decode_ssss_index) - 1 - decode_codeword) + diff_pre_value;
				input_block[decode_row][decode_col].value[0][0] = diff_now_value;
				diff_pre_value = diff_now_value;
				set_decode_flag(false, false, true, false);
			}
			else
			{
				diff_now_value = decode_codeword + diff_pre_value;
				input_block[decode_row][decode_col].value[0][0] = diff_now_value;
				diff_pre_value = diff_now_value;
				set_decode_flag(false, false, true, false);
			}

		}
		else if (decode_counter == 1 && decode_ssss_index == 0)
		{
			/*output previous dc value*/
			input_block[decode_row][decode_col].value[0][0] = diff_pre_value;
			set_decode_flag(false, false, true, false);
		}
	}
}	

void ac_decode()
{
	if (ac_codeword_decoding)
	{
		if (decode_counter >= 2)
		{
			map < int, map <int, ac_table_de>> ::iterator ac_codeword = dec_ac_t.find(decode_codeword);
			/*if find codeword*/
			if (ac_codeword != dec_ac_t.end())
			{
				map <int, ac_table_de> ::iterator run_and_size = ac_codeword->second.find(decode_counter);
				decode_run = run_and_size->second.run;
				decode_size = run_and_size->second.size;
				set_decode_flag(false, false, false, true);
			}
		}
	}
	else if (ac_diff_decoding)
	{
		int block_inner_row_index, block_inner_col_index;
		int diff_now_value = 0;
		/*set Z-ZAG */
		block_inner_row_index = zigZagOrder[decode_ac_counter] / 8;
		block_inner_col_index = zigZagOrder[decode_ac_counter] % 8;

		if (decode_size != 0)
		{
			if (decode_counter == decode_size)
			{
				for (int i = 0; i < decode_run; i++)
				{
					input_block[decode_row][decode_col].value[block_inner_row_index][block_inner_col_index] = 0;
					/*refresh block index*/
					decode_ac_counter++;
					block_inner_row_index = zigZagOrder[decode_ac_counter] / 8;
					block_inner_col_index = zigZagOrder[decode_ac_counter] % 8;
				}
				/*output  now value*/
				if (decode_codeword < powl(2, decode_size) / 2)// if decodeword is negative
				{
					diff_now_value = -(powl(2, decode_size) - 1 - decode_codeword);
					input_block[decode_row][decode_col].value[block_inner_row_index][block_inner_col_index] = diff_now_value;
				}
				else
				{
					diff_now_value = decode_codeword;
					input_block[decode_row][decode_col].value[block_inner_row_index][block_inner_col_index] = diff_now_value;
				}
				decode_ac_counter++;
				set_decode_flag(false, false, true, false);
			}
		}
		else
		{
			if (decode_run == 0)//EOB
			{
				while (decode_ac_counter < 64)
				{
					block_inner_row_index = zigZagOrder[decode_ac_counter] / 8;
					block_inner_col_index = zigZagOrder[decode_ac_counter] % 8;
					input_block[decode_row][decode_col].value[block_inner_row_index][block_inner_col_index] = 0;
					decode_ac_counter++;
				}
				int temp = decode_codeword;
				set_decode_flag(true, false, false, false);
				decode_codeword = temp;
				decode_counter = 1;
				decode_ac_counter = 1;
				if (decode_col < 63)
				{
						decode_col++;
				}
				else
				{
		
					decode_col = 0;
					decode_row++;
					if (decode_row == 64)
					{
						complete_decoding = true;
					}
				}

			}
			else if (decode_run == 15)
			{
				
				for (int i = 0; i < decode_run; i++)
				{
					input_block[decode_row][decode_col].value[block_inner_row_index][block_inner_col_index] = 0;
					/*refresh block index*/
					decode_ac_counter++;
					block_inner_row_index = zigZagOrder[decode_ac_counter] / 8;
					block_inner_col_index = zigZagOrder[decode_ac_counter] % 8;
				}
				int temp = decode_codeword;
				set_decode_flag(false, false, true, false);
				decode_codeword = temp;
				decode_counter = 1;
			}
		}
	}
}
void decoder(int most_symbol_bit)
{
	decode_counter++;
	decode_codeword <<= 1;
	decode_codeword += most_symbol_bit;
	if (decode_counter > 0)
	{
		if (dc_codeword_decoding || dc_diff_decoding)
		{
			dc_decode();
		}
		else
		{
			ac_decode();
		}
	}
}
int main()
{

	string filename;
	int mode, QF;
	cout << "Please select your mode " << endl;
	cout << "1) Encode " << endl << "2) Decode " <<endl<<"3) Count PSNR"<< endl;
	cin >> mode;
	if (mode == 1)
	{
		cout << "Please input your filename (The filename must be xxx.raw)	" << endl;
		cin >> filename;
		cout << "Please input your QF" << endl;
		//input file 
		FILE *file = fopen(filename.c_str(), "rb");
		filename[filename.length() - 3] = 'o', filename[filename.length() - 2] = 'u', filename[filename.length() - 1] = 't';
		file_out = fopen(filename.c_str(), "wb");
		int counter = 0;
		//read file
		readblock(file);
		DCT();
		/*Quantization*/
		cin >> QF;
		quantization(QF);
		for (int i = 0; i < 64; i++)
		{
			for (int j = 0; j < 64; j++)
			{
				dc_encode(i, j);
				ac_encode(i, j);
			}
		}
		if (output_counter > 0)
		{
			output <<= (8 - output_counter);
			fwrite(&output, sizeof(char), 1, file_out);
		}
		fclose(file);
		fclose(file_out);
	}
	else if (mode == 2)
	{
		cout << "Please input your filename (The filename must be xxx.out)" << endl;
		cin >> filename;
		cout << "Please input your QF" << endl;
		cin >> QF;
		//input file 
		FILE *file = fopen(filename.c_str(), "rb");
		filename[filename.length() - 3] = 'r', filename[filename.length() - 2] = 'e', filename[filename.length() - 1] = 't';
		file_out = fopen(filename.c_str(), "wb");
		set_decode_actable();
		char temp;
		while (!feof(file))
		{
			fread(&temp, sizeof(char), 1, file);
			for (int i = 0; i < 8; i++)
			{

				decoder(MSB(temp));
				if (complete_decoding)
				{
					break;
				}
				temp <<= 1;
			}
			if (complete_decoding)
			{
				break;
			}
		}
		de_quantization(QF);
		I_DCT();
		writeblock();
		fclose(file);
		fclose(file_out);
	/*	while (1)
		{
			int a = 0;
			int b = 0;
			cin >> a >> b;
			print_block(a, b);
		}*/
	}
	else if (mode == 3)
	{
		string filename_jpg;
		cout << "Please input your filename (The filename must be xxx.raw)" << endl;
		cin >> filename;
		cout << "Please input your 'JPG' filename (The filename must be xxx.ret)" << endl;
		cin >> filename_jpg;
		FILE *file = fopen(filename.c_str(), "rb");
		FILE *file_jpg = fopen(filename_jpg.c_str(), "rb");
		unsigned char raw_temp;
		unsigned char jpg_temp;
		double MSE = 0;
		double PSNR;
		while (!feof(file))
		{
			fread(&raw_temp, sizeof(char), 1, file);
			fread(&jpg_temp, sizeof(char), 1, file_jpg);
			MSE += pow(raw_temp-jpg_temp,2	) / 262144.0;
		}

		PSNR = 10 * log10( 65025.0 / MSE);
		cout << PSNR << endl;
	}
	system("PAUSE");
	return 0;
}