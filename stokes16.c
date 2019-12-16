#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include <time.h>


#define DATA_SYMBOL	(6250)
#define DATA_BIT	(DATA_SYMBOL *4)
#define TR_SYMBOL	(1000)
#define TR_BIT	(TR_SYMBOL *4)
#define DEL_SYMBOL	(942)
#define DEL_BIT	(DEL_SYMBOL * 4)
#define N_SYMBOL (TR_SYMBOL + DATA_SYMBOL + DEL_SYMBOL)
#define N_BIT (N_SYMBOL*4)
#define	N_BIT_PAM8	(4)
#define P_PAD (64)
#define P_ANALOG (P_PAD*N_SYMBOL)
#define RATE_BAUD (25e9) //ボーレート：25GBaud



#define ROLLOFF_FACTOR_RO (2)
#define BAND_RATE_ROLLOFF (0.05)

#define DIST_Z	(5.0E+3)
#define DIST_H	(1.0E+2)
#define SSFM_COUNT	((int)DIST_Z) / ((int)DIST_H)
#define SENSITIVITY	(0.6)
#define flag_SHOT (1) 	//ショット雑音on off
#define flag_THERMAL (1) //熱雑音フラグ

#define STEP_COEF   (0.1)
#define STEP_MIN	(1)
#define STEP_MAX	(3)
#define TAP_MIN		(1)
#define TAP_MAX		(3)
#define SFT_MIN		(0)
#define SFT_MAX		(2)				/////////ここ変えた

#define STEP_NUM	(STEP_MAX-STEP_MIN+1)		// int型、EVM変数格納配列数
#define TAP_NUM		(TAP_MAX-TAP_MIN+1)
#define SFT_NUM		(SFT_MAX-SFT_MIN+1)

#define PI (3.141592653589793238)
#define PMD_PARAM	(0.02)
#define C_SPEED		(2.99792458E+8)
#define BETA1       (PMD_PARAM / (sqrt((double)DIST_H*1.0E-3)))		    // double型、定偏波ファイバ毎の群遅延の平均値[ps/km]
#define ALPHA     	(2.0E-4)					//ファイバ損失のパラメータ
#define BETA2       (-20.0E-27)				//波長分散のパラメータ
#define GAMMA       (1.27E-3)				//非線形のパラメータ
#define CN_WAVE	    (1550.0E-9)						// double型、中心波長[m]
#define RESISTANCE	(1000.0)							//負荷抵抗
#define BOLTZMAN	(1.38064852e-23)				//ボルツマン定数
#define	TEMPRATURE	(300.0)							//気温[K]
#define ST_ASE_SEED	(10)							//ASE雑音シード開始値
#define H_PLANCK 	(6.62607004e-34)				//プランク定数
#define CENTER_FREQ	((double)C_SPEED / CN_WAVE)		//中心周波数
#define PMD_SEED (2)

#define CONST_D (0.631475730333305)				//ポアンカレ球に内接する、正十二面体の内接する、正二十面体の、原点から面の距離ｄ（各面の方程式の定数項）
#define CONST_A (0.577350269189626)				//面の方程式の定数A
#define CONST_B (0.577350269189626)				//面の方程式の定数B
#define CONST_C (0.577350269189626)				//面の方程式の定数C = 1/sqrt(3)

double dBm;							//入射パワー[dBm]
double POWER_LAUNCH;	//入射パワー[mW]
double S0;	//全光パワー[W]

int LMS_TAP;
int SFT_SYMBOL;
double LMS_STEP;
// double* test_e = (double*)malloc(sizeof(double) * TR_SYMBOL *3);

void CSV_freq(char *name, double *deta1);
void CSV_5(char *name,int num, double *deta_four);
void map(int *bitStream, double *S_digit);
void StoE(double *S_data, double *E_data);
void EtoS(double *E_data, double *S_data);
void RaisedCosine(double *E_digit, double *E_ana);
void SSFM(double *O_data, double *save_deltabeta);
double nrand();
void EtoO(double *E_ana,double *O_ana);
void Stokes_nomal(double *S_data, double *S_nomal, int num);
void Sample(double *analogData,double *digitData, int n);
void coupler(double *analogData);
void SHOTnoise(double *analogData, int noise_seed, int flag);
double gaussian(double mu, double sgm);
double genRand();
void res_sen(double *O_data, double *E_data);
void THERMALnoise(double *S_data, int noise_seed, int flag);
int FIR_filter(double *S_aft, double *_pre);
int DEL_TSYMBOL(double *S_data);
int Stokes_mod(int *sequence, double *S_data);
void catch_Symbol(double *S_data, int Symbol, int num);


int main()
{
	printf(" main ok\n DIST: %.2e\n flag_SHOT: %d\n flag_THERMAL: %d\n",DIST_Z,flag_SHOT,flag_THERMAL);
	int end,k;
	end = 3;				//////変えた
	double* ber_stock = (double*)malloc(sizeof(double) * N_SYMBOL);
	int a=0;
	int loop;
	int loop_end = 3;					//変えた
	double progress;
	double evm_sum = 0.0;
	double evm_ave[100];
	double distance[100];
	double SN_ratio[100];
	double thermal_sigma;
	double sigma_c;
	//int PMD_SEED;
	char* filename = (char*)malloc(sizeof(char)*1024);

	//サンプル点ずらしループ(aループ)
	for(a = 0 ; a < 1 ; a++){					//変えた
		printf("a loop ok\n");

		//伝送距離ずらしループ(kループ)
		for(k = 0;k<end;k++){
			printf("k loop ok\n");


			dBm = (5.0 * (double)k) - 20.0; //入射パワー[dBm]
			POWER_LAUNCH = ((1.0e-3) * pow(10.0,dBm/10.0));	//[W]
			S0 = POWER_LAUNCH * 1000; //[mW]

			int i,j;
			FILE* fp;
			int bitStream[N_BIT];
			int seed = 1;
			double* S_parameter_tx = (double*)malloc(sizeof(double) * N_SYMBOL * 4);
			double* S_parameter_rx = (double*)malloc(sizeof(double) * N_SYMBOL * 4);
			double* S_normaled_preFIR_rx = (double*)malloc(sizeof(double) * N_SYMBOL * 4);
			double* fil_S = (double*)malloc(sizeof(double)* N_SYMBOL *4);
			double* fil_S_nomaled = (double*)malloc(sizeof(double)* N_SYMBOL *4);
			double* S_FIRed_rx = (double*)malloc(sizeof(double) * N_SYMBOL * 4);
			double* S_FIRed_normaled_rx = (double*)malloc(sizeof(double) * N_SYMBOL * 4);
			double* E_sym = (double*)malloc(sizeof(double) * N_SYMBOL *4);								//xy偏波パラメータ　0~1:Exre,Exim,Eyre,Eyim
			double* E_analog = (double*)malloc(sizeof(double) * P_ANALOG * 4);
			double* O_analog = (double*)malloc(sizeof(double) * P_ANALOG * 4);
			double* E_analog_rx = (double*)malloc(sizeof(double) * P_ANALOG * 4);
			double* S_analog_rx = (double*)malloc(sizeof(double) * P_ANALOG * 4);
			double* E_sym_rx = (double*)malloc(sizeof(double) * N_SYMBOL * 4);
			int bit;
			double frequency_x_re[P_ANALOG],frequency_x_im[P_ANALOG],frequency_y_re[P_ANALOG],frequency_y_im[P_ANALOG];
			double rxdigit_x_re[N_SYMBOL],rxdigit_x_im[N_SYMBOL],rxdigit_y_re[N_SYMBOL],rxdigit_y_im[N_SYMBOL];								//受信シンボル
			double s0[N_SYMBOL],s1[N_SYMBOL],s2[N_SYMBOL],s3[N_SYMBOL];														//受信コンスタ
			double y[N_SYMBOL];
			double rxdigit_y, rxdigit_x;																															//受信x,yパワー
			double power_rx=0.0;																																		//受信パワー


			double* filter = (double*)malloc(sizeof(double) * P_ANALOG); //フィルター配列

			double optical_re[P_ANALOG], optical_im[P_ANALOG];
			double power_dsp_x,power_dsp_y;			//電圧信号の平均パワー
			double power_adjust_x,power_adjust_y;		//パワー変換係数α
			double power_adjust;

			double *save_deltabeta = (double*)malloc(sizeof(double)*SSFM_COUNT);

			//FIRフィルタ用の配列

			double fil_s1[N_SYMBOL];
			double fil_s2[N_SYMBOL];
			double fil_s3[N_SYMBOL];
			double AFTER_EVM_1[STEP_NUM][TAP_NUM][SFT_NUM];			// 重ね合わせ後のEVM格納
			double AFTER_EVM_2[STEP_NUM][TAP_NUM][SFT_NUM];
			double AFTER_EVM_3[STEP_NUM][TAP_NUM][SFT_NUM];

			int step_lp,tap_lp,sft_lp;
			double evm;
			int rece_bit[N_BIT];									// 受信ビット列
			int error_bit;											//　エラービットの数
			double BER;

			//雑音定義用関数
			static int noise_seed = ST_ASE_SEED;

			//BER loop
			for(loop = 0; loop < loop_end; loop++){

				printf("loop ok %d\n",loop);

				map(bitStream, S_parameter_tx);								//送信波形作成
					strcpy(filename,"01_tx_Stokes.csv");
					CSV_5(filename,N_SYMBOL,S_parameter_tx);
				StoE(S_parameter_tx, E_sym);
					strcpy(filename,"02_tx_symbol.csv");
					CSV_5(filename,N_SYMBOL,E_sym);
				RaisedCosine(E_sym, E_analog);							//RC関数
					strcpy(filename,"03_tx_analog.csv");
					CSV_5(filename,P_ANALOG,E_analog);
				EtoO(E_analog,O_analog);
					strcpy(filename,"04tx_optical.csv");
					CSV_5(filename,P_ANALOG,O_analog);

				//伝送路
				SSFM(O_analog, save_deltabeta);
					strcpy(filename,"05_txSSFMed_optical.csv");
					CSV_5(filename,P_ANALOG,O_analog);

				//受信部
				//光カプラ
				coupler(O_analog);
				//ショット雑音
				SHOTnoise(O_analog,noise_seed,flag_SHOT);

				//受光感度の計算[W]to[A]
				res_sen(O_analog,E_analog_rx);
					strcpy(filename,"06_rx_analog.csv");
					CSV_5(filename,P_ANALOG,E_analog_rx);

				//偏波パラからストークスパラメータ
				EtoS(E_analog_rx,S_analog_rx);
					strcpy(filename,"07_rx_StokesAnalog.csv");
					CSV_5(filename,P_ANALOG,S_analog_rx);

				//熱雑音
				THERMALnoise(S_analog_rx,noise_seed,flag_THERMAL);
				if(flag_THERMAL == 1){
					strcpy(filename,"08_rxTHERMALed_StokesAnalog.csv");
					CSV_5(filename,P_ANALOG,S_analog_rx);
				}

				//サンプリング(a点だけずらす)
				Sample(S_analog_rx, S_parameter_rx,a);
					strcpy(filename,"09_rx_StokesSymbol.csv");
					CSV_5(filename,N_SYMBOL,S_parameter_rx);

				//規格化
				Stokes_nomal(S_parameter_rx, S_normaled_preFIR_rx, N_SYMBOL);
					strcpy(filename,"10_rx_normaledStokes_preFIR.csv");
					CSV_5(filename,N_SYMBOL,S_normaled_preFIR_rx);

				//FIR filter
				for(step_lp = STEP_MIN; step_lp <= STEP_MAX; step_lp++){
					for(tap_lp = TAP_MIN; tap_lp <= TAP_MAX; tap_lp++){
						for(sft_lp = SFT_MIN; sft_lp <= SFT_MAX; sft_lp++){
							LMS_STEP = (double)step_lp * STEP_COEF;
							LMS_TAP = tap_lp;
							SFT_SYMBOL = sft_lp;

							printf("STEP: %.2f\nTAP: %d\nSFT: %d\n",LMS_STEP,tap_lp,sft_lp);

							memcpy(fil_S, S_parameter_rx, sizeof(double) * N_SYMBOL *4);

							FIR_filter(fil_S, S_parameter_tx);
						}
					}
					printf("%d\n",step_lp);
				}

				memcpy(S_FIRed_rx, fil_S, sizeof(double) * N_SYMBOL *4);
					Stokes_nomal(fil_S,fil_S_nomaled,N_SYMBOL);
					strcpy(filename,"11_rx_FIRedStokes.csv");
					CSV_5(filename,N_SYMBOL,fil_S_nomaled);

				//データシンボルのみ取り出す
				DEL_TSYMBOL(S_FIRed_rx);
					memcpy(S_FIRed_rx, fil_S, sizeof(double) * N_SYMBOL *4);
					Stokes_nomal(fil_S,fil_S_nomaled,N_SYMBOL);
					strcpy(filename,"12_rx_DataSymbol.csv");
					CSV_5(filename,DATA_SYMBOL,fil_S_nomaled);

				if(loop == loop_end - 1){
					Stokes_nomal(S_FIRed_rx,S_FIRed_normaled_rx, DATA_SYMBOL);
					snprintf(filename,1024,"13_Sift_%d,TR_%d,step_%.2f_%.1fkm_%.1fdBm.csv",a,TR_SYMBOL,LMS_STEP,DIST_Z/1000.0,dBm);
					CSV_5(filename,DATA_SYMBOL,S_FIRed_normaled_rx);
				}


				//ストークスから便と列に直し、エラービットを数える
				// Stokes_mod(rece_bit,S_FIRed_rx);
				//
				// error_bit += calc_error(bitStream, rece_bit);
				// printf("error bit : %d\n",error_bit);
				// }
				//

			}
			// BER = log10((double)error_bit / (double)(DATA_BIT * (loop_end-1)));
			//("logBER : %.2f\n",BER);


			//SNRの導出
			thermal_sigma = sqrt(4.0 * BOLTZMAN * TEMPRATURE * RATE_BAUD / RESISTANCE);
			sigma_c = thermal_sigma/3;
			SN_ratio[k] = dBm;
			printf("SN_ratio[%d] = %lf\n",k,SN_ratio[k]);


			//BERの格納
			//ber_stock[k] = BER;
			//evm_ave[k] = evm_sum / ((double)loop_end - 1.0);
			//evm = 0.0;
			//evm_sum = 0.0;



			//error_bit = 0.0;


			free(S_parameter_tx);
			free(S_parameter_rx);
			free(S_normaled_preFIR_rx);
			free(S_FIRed_rx);
			free(S_FIRed_normaled_rx);
			free(fil_S);
			free(fil_S_nomaled);
			free(E_sym);
			free(E_analog);
			free(O_analog);
			free(E_analog_rx);
			free(E_sym_rx);
			free(filter);
			free(save_deltabeta);


			//エラービットの初期化
			//error_bit = 0.0;
		}//kループの閉じ
		//	snprintf(filename,1024,"Sift_%d,BER_result.csv",a);
		//	CSV_6(k,SN_ratio,ber_stock,evm_ave,filename);

	}//aループの閉じ


	free(filename);
	free(ber_stock);
	return 0;
}






//============================================
//ストークスパラメータ to x,yシンボル
//============================================

void StoE(double *S_data, double *E_data){
	int i;
  for(i=0;i<N_SYMBOL;i++){
	  E_data[4*i +0] = sqrt((S_data[4*i + 0] + S_data[4*i +1])/2);
		E_data[4*i +1] = 0;
	  E_data[4*i +2] = S_data[4*i +2] / sqrt(2 * (S_data[4*i +0]+ S_data[4*i +1]));
	  E_data[4*i +3] = S_data[4*i +3] / sqrt(2 * (S_data[4*i +0]+ S_data[4*i +1]));
  }

	return;
}


//==========
//周波数軸csv作成
//==========

void CSV_freq(char *name, double *deta1){
	int i;

	FILE* fp;
	fp = fopen(name,"w");
	fprintf(fp,"点,周波数,power\n");
	for(i=0;i<P_ANALOG;i++){
		if(i<P_ANALOG/2){
			fprintf(fp,"%d,%f,%f\n",i,(double)RATE_BAUD/(double)N_SYMBOL * (double)i, deta1[i]);
		}
		else{
			fprintf(fp,"%d,%f,%f\n",i - (int)P_ANALOG,(double)RATE_BAUD/(double)N_SYMBOL * (double)(i-P_ANALOG), deta1[i]);
		}
	}
}

//===========
//CSV_5:点数、データ4種
//===========

void CSV_5(char *name,int num,double *data_four){
	int i;
	FILE* fp;

	fp = fopen(name,"w");
	for(i=0;i<num;i++){
		fprintf(fp,"%d,%f,%f,%f,%f\n",i,data_four[4 * i +0],data_four[4 * i +1],data_four[4 * i +2],data_four[4 * i +3]);
	}
	fclose(fp);
	return;
}

void map(int *bitStream, double *S_digit){
	int i;
	int bit;
//=====================
//ビット列生成
//=====================

	srand(1);
	for(i = 0;i < N_BIT; i++)
	{
		bitStream[i] = rand() % 2;
	}

//=====================
//マッピング
//=====================

	for(i=0;i<N_SYMBOL;i++)
	{
		bit = 8*bitStream[i*4+0] + 4*bitStream[i*4+1] + 2*bitStream[i*4+2] + bitStream[i*4+3];
		S_digit[4*i +0] = 1;
		if(bit==0)
		{
			S_digit[4*i +1] = 1/sqrt(3);
			S_digit[4*i +2] = 1/sqrt(3);
			S_digit[4*i +3] = 1/sqrt(3);
		}
		else if(bit==1)
		{
			S_digit[4*i +1] = 1/sqrt(3);
			S_digit[4*i +2] = -1/sqrt(3);
			S_digit[4*i +3] = 1/sqrt(3);
		}
		else if(bit==2)
		{
			S_digit[4*i +1] = -1/sqrt(3);
			S_digit[4*i +2] = 1/sqrt(3);
			S_digit[4*i +3] = 1/sqrt(3);
		}
		else if(bit==3)
		{
			S_digit[4*i +1] = -1/sqrt(3);
			S_digit[4*i +2] = -1/sqrt(3);
			S_digit[4*i +3] = 1/sqrt(3);
		}
		else if(bit==4)
		{
			S_digit[4*i +1] = 1/sqrt(3);
			S_digit[4*i +2] = 1/sqrt(3);
			S_digit[4*i +3] = -1/sqrt(3);
		}
		else if(bit==5)
		{
			S_digit[4*i +1] = 1/sqrt(3);
			S_digit[4*i +2] = -1/sqrt(3);
			S_digit[4*i +3] = -1/sqrt(3);
		}
		else if(bit==6)
		{
			S_digit[4*i +1] = -1/sqrt(3);
			S_digit[4*i +2] = 1/sqrt(3);
			S_digit[4*i +3] = -1/sqrt(3);
		}
		else if(bit==7)
		{
			S_digit[4*i +1] = -1/sqrt(3);
			S_digit[4*i +2] = -1/sqrt(3);
			S_digit[4*i +3] = -1/sqrt(3);
		}
		else if(bit==8)
		{
			S_digit[4*i +1] = 0;
			S_digit[4*i +2] = 2/(1+sqrt(5)) /sqrt(3);
			S_digit[4*i +3] = (1+sqrt(5))/2 /sqrt(3);
		}
		else if(bit==9)
		{
			S_digit[4*i +1] = 0;
			S_digit[4*i +2] = 2/(1+sqrt(5)) /sqrt(3);
			S_digit[4*i +3] = -1*(1+sqrt(5))/2 /sqrt(3);
		}
		else if(bit==10)
		{
			S_digit[4*i +1] = 0;
			S_digit[4*i +2] = -2/(1+sqrt(5)) /sqrt(3);
			S_digit[4*i +3] = (1+sqrt(5))/2 /sqrt(3);
		}
		else if(bit==11)
		{
			S_digit[4*i +1] = 0;
			S_digit[4*i +2] = -2/(1+sqrt(5)) /sqrt(3);
			S_digit[4*i +3] = -1*(1+sqrt(5))/2 /sqrt(3);
		}
		else if(bit==12)
		{
			S_digit[4*i +1] = (1+sqrt(5))/2 /sqrt(3);
			S_digit[4*i +2] = 0;
			S_digit[4*i +3] = 2/(1+sqrt(5)) /sqrt(3);
		}
		else if(bit==13)
		{
			S_digit[4*i +1] = (1+sqrt(5))/2 /sqrt(3);
			S_digit[4*i +2] = 0;
			S_digit[4*i +3] = -2/(1+sqrt(5)) /sqrt(3);
		}
		else if(bit==14)
		{
			S_digit[4*i +1] = -1*(1+sqrt(5))/2 /sqrt(3);
			S_digit[4*i +2] = 0;
			S_digit[4*i +3] = 2/(1+sqrt(5)) /sqrt(3);
		}
		else if(bit==15)
		{
			S_digit[4*i +1] = -1*(1+sqrt(5))/2 /sqrt(3);
			S_digit[4*i +2] = 0;
			S_digit[4*i +3] = -2/(1+sqrt(5)) /sqrt(3);
		}
	}

	return;
}


//===========
//RaisedCosine生成 シンボルn_symbolからアナログ波形p_analog
//===========

void RaisedCosine(double *E_digit, double *E_ana){
	int i;
	int i_digit, i_internal, i_analog;
	int stable, raise;
	double x_re_vs, x_im_vs, x_re_vg, x_im_vg;
	double y_re_vs, y_im_vs, y_re_vg, y_im_vg;
	double a, b;
	stable = P_PAD / 2.0;
	raise = P_PAD - stable;

	for(i_digit = 0; i_digit < N_SYMBOL; i_digit++){
		x_re_vs = E_digit[4*i_digit + 0];
		x_im_vs = E_digit[4*i_digit + 1];
		y_re_vs = E_digit[4*i_digit + 2];
		y_im_vs = E_digit[4*i_digit + 3];

		if(i_digit != N_SYMBOL -1 )
		{
			x_re_vg = E_digit[4*i_digit + 4];
			x_im_vg = E_digit[4*i_digit + 5];
			y_re_vg = E_digit[4*i_digit + 6];
			y_im_vg = E_digit[4*i_digit + 7];
		}
		else
		{
			x_re_vg = E_digit[0 +0];
			x_im_vg = E_digit[0 +1];
			y_re_vg = E_digit[0 +2];
			y_im_vg = E_digit[0 +3];
		}

		for(i_internal = 0; i_internal < P_PAD; i_internal++){

			i_analog = i_digit * P_PAD + i_internal;

			a = (double)raise / 2.0;
			b = (2.0*(double)stable + (double)raise ) /4.0;

			if(i_internal < stable){
				E_ana[4*i_analog+0] = x_re_vs;
				E_ana[4*i_analog+1] = x_im_vs;
				E_ana[4*i_analog+2] = y_re_vs;
				E_ana[4*i_analog+3] = y_im_vs;
			}

			else{
				E_ana[4*i_analog+0] = (x_re_vs - x_re_vg)*(0.5 + 0.5*cos(PI*((double)i_internal - 2.0*b + a) / (2.0*a))) + x_re_vg;
				E_ana[4*i_analog+1] = (x_im_vs - x_im_vg)*(0.5 + 0.5*cos(PI*((double)i_internal - 2.0*b + a) / (2.0*a))) + x_im_vg;
				E_ana[4*i_analog+2] = (y_re_vs - y_re_vg)*(0.5 + 0.5*cos(PI*((double)i_internal - 2.0*b + a) / (2.0*a))) + y_re_vg;
				E_ana[4*i_analog+3] = (y_im_vs - y_im_vg)*(0.5 + 0.5*cos(PI*((double)i_internal - 2.0*b + a) / (2.0*a))) + y_im_vg;
			}
		}
	}
	return ;
}

void EtoS(double *E_data, double *S_data){	//Ex,EyからStokes
	double power_rx = 0.0;
	double rxdigit_y = 0.0 , rxdigit_x = 0.0;
	double power_adjust;
	int i;
	for(i=0;i<P_ANALOG;i++){
		rxdigit_y = sqrt(E_data[4 * i + 2] * E_data[4 * i + 2] + E_data[4 * i + 3] * E_data[4 * i + 3]);
		rxdigit_x = sqrt(E_data[4 * i + 0] * E_data[4 * i + 0] + E_data[4 * i + 1] * E_data[4 * i + 1]);
		S_data[4 * i +0] = rxdigit_x * rxdigit_x + rxdigit_y * rxdigit_y;
		S_data[4 * i +1] = rxdigit_x * rxdigit_x - rxdigit_y * rxdigit_y;
		S_data[4 * i +2] = 2 * rxdigit_x * rxdigit_y * (E_data[4 * i + 2] / rxdigit_y);
		S_data[4 * i +3] = 2 * rxdigit_x * rxdigit_y * (E_data[4 * i + 3] / rxdigit_y);
	}

	return ;
}

void SSFM(double *O_data, double *save_deltabeta){
	int i, j;                                       //ループ用
	int index;										//線形効果用変数
	double now_distance;					      	// 計算後の現在位置を代入
	double power_x, power_y; 						//振幅の二乗値
	double xpower_re,xpower_im,ypower_re,ypower_im; //各偏波のパワー
	double x_re, x_im, y_re, y_im;
	double coef,coef0;
	double df,alphaC;								//周波数分解能, ファイバ損失係数
	double omega, omega2;							//線形効果ω
	static int noise_seed = ST_ASE_SEED;			// ASE雑音シード
	double amp_gain;							// EDFAの利得を代入

	double beta1x, beta1y, deltabeta2;			// 各過程における分散パラメータを代入
	double x_diff, y_diff;						// x、yのパワーの差を代入
	double Lx, Ly, Mx, My, Nx, Ny;				// 計算に使用する変数
	double Ex1, Ex2, Ey1, Ey2, DIx, DIy;		// 計算に使用する変数
	double *deltabeta ;	// SSFM変数Δβを格納
	double *phi;	    // SSFM変数φを格納
	double *theta;	    // SSFM変数θを格納
	deltabeta = (double*)malloc(sizeof(double) * SSFM_COUNT);
	phi  = (double*)malloc(sizeof(double) * SSFM_COUNT);
	theta = (double*)malloc(sizeof(double) * SSFM_COUNT);

	double *xfreq_re;
	double *xfreq_im;
	double *yfreq_re;
	double *yfreq_im;
	xfreq_re = (double*)malloc(sizeof(double)*P_ANALOG);
	xfreq_im = (double*)malloc(sizeof(double)*P_ANALOG);
	yfreq_re = (double*)malloc(sizeof(double)*P_ANALOG);
	yfreq_im = (double*)malloc(sizeof(double)*P_ANALOG);

	df = (RATE_BAUD / (double)N_SYMBOL); //周波数分解能の計算
	alphaC = ALPHA * 0.10 * log(10.0); //ファイバ損失係数の計算

	fftw_complex *in, *out; //FFTの変数宣言
	fftw_plan fft, ifft;
	in = fftw_malloc(sizeof(fftw_complex) * P_ANALOG);
	out = fftw_malloc(sizeof(fftw_complex) * P_ANALOG);
	fft = fftw_plan_dft_1d(P_ANALOG, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	ifft = fftw_plan_dft_1d(P_ANALOG, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

// PMD分布乱数生成
	srand48(PMD_SEED);

	/* ランダムパラメータ計算 */
	for (i = 0; i < (int)SSFM_COUNT; i++){
		deltabeta[i] = (nrand() * BETA1 * 0.2 + BETA1) * 1.0E-15;
		save_deltabeta[i] = deltabeta[i];
		phi[i] = deltabeta[i] * DIST_H * 2.0 * PI * C_SPEED / CN_WAVE;
		//phi[i] = 0.0;//群速度遅延に関係
		//theta[i] = 0.0;//PI/4.0;//軸の回転に関係
		theta[i] = ((PI / 50.0) * ((double)rand() / ((double)RAND_MAX + 1) * 101.0) - 50.0);
	}

	for (i = 0; i < (int)SSFM_COUNT; i++){
		/* パラメータの設定 */
		now_distance = (double)(i + 1) * DIST_H;//相関長に関するパラメータ
		beta1x = deltabeta[i] / 2.0;
		beta1y = -deltabeta[i] / 2.0;
		deltabeta2 = deltabeta[i] * DIST_H * 2.0 * PI * C_SPEED / CN_WAVE;

		//非線形効果
		for (j = 0; j < P_ANALOG; j++)
		{
			x_re = O_data[4*j+0];
			x_im = O_data[4*j+1];
			y_re = O_data[4*j+2];
			y_im = O_data[4*j+3];

			xpower_re = x_re * x_re;
			xpower_im = x_im * x_im;
			ypower_re = y_re * y_re;
			ypower_im = y_im * y_im;

			power_x = xpower_re + xpower_im; //|A|^2
			power_y = ypower_re + ypower_im;
			x_diff = xpower_re - xpower_im;
			y_diff = ypower_re - ypower_im;

			Lx = x_diff * x_diff + 4.0 * xpower_re * xpower_im;
			Ly = y_diff * y_diff + 4.0 * ypower_re * ypower_im;
			Mx = x_diff * y_diff + 4.0 * x_re * x_im * y_re * y_im;
			My = x_diff * y_diff + 4.0 * x_re * x_im * y_re * y_im;
			Nx = x_re * x_im * y_diff - y_re * y_im * x_diff;
			Ny = y_re * y_im * x_diff - x_re * x_im * y_diff;

			coef = DIST_H * GAMMA; //h*γ
			coef0 = (double)DIST_H * (double)GAMMA * power_x;

			// パラメトリック効果有りの場合（未整理）

			// エラー回避用分岐
			/* if (Lx == 0)
			{
			Ex1 = 0;
			Ex2 = 0;
			}
			else
			{
			Ex1 = exp(coef * (2.0 / 3.0 * power_x * cos(2.0*deltabeta2*now_distance) * Nx / Lx + 1.0 / 3.0 * power_x * sin(2.0*deltabeta2*now_distance) * Mx / Lx));
			Ex2 = coef * ((power_x + 2.0 / 3.0 * power_y) - 2.0 / 3.0 * power_x * sin(2.0*deltabeta2*now_distance) * Nx / Lx + 1.0 / 3.0 * power_x * cos(2.0*deltabeta2*now_distance) * Mx / Lx);
			}

			if (Ly == 0)
			{
			Ey1 = 0;
			Ey2 = 0;
			}
			else
			{
			Ey1 = exp(coef * (2.0 / 3.0 * power_y * cos(2.0*deltabeta2*now_distance) * Ny / Ly + 1.0 / 3.0 * power_y * sin(2.0*deltabeta2*now_distance) * My / Ly));
			Ey2 = coef * ((power_y + 2.0 / 3.0 * power_x) + 2.0 / 3.0 * power_y * sin(2.0*deltabeta2*now_distance) * Ny / Ly + 1.0 / 3.0 * power_y * cos(2.0*deltabeta2*now_distance) * My / Ly);
			}

			// 非線形効果の代入
			analog_re_x[j] = Ex1 * (x_re * cos(Ex2) - x_im * sin(Ex2));
			analog_im_x[j] = Ex1 * (x_re * sin(Ex2) + x_im * cos(Ex2));
			analog_re_y[j] = Ey1 * (y_re * cos(Ey2) - y_im * sin(Ey2));
			analog_im_y[j] = Ey1 * (y_re * sin(Ey2) + y_im * cos(Ey2)); */

			/* パラメトリック効果無しの場合 */
		 	Ex1 = 0.0;
			Ex2 = coef * (power_x + 2.0 / 3.0 * power_y);
			Ey1 = 0.0;
			Ey2 = coef * (power_y + 2.0 / 3.0 * power_x);

			O_data[4*j+0] = (x_re * cos(Ex2)) - (x_im * sin(Ex2));
			O_data[4*j+1] = (x_re * sin(Ex2)) + (x_im * cos(Ex2));
			O_data[4*j+2] = (y_re * cos(Ey2)) - (y_im * sin(Ey2));
			O_data[4*j+3] = (y_re * sin(Ey2)) + (y_im * cos(Ey2));
		}

		//SSFMcheck
		//sprintf(filename2,"SSFM_check_NLprosess_x_NL_step dist %d m_.csv", i);
		//make_CSV2(P_ANALOG, analog_re_x, analog_im_x, filename2);

		//フーリエ変換シミュレーションで計算できるようにフーリエ変換で周波数領域にする
		for (j = 0; j < P_ANALOG; j++)
		{
			in[j][0] = O_data[4*j+0];
			in[j][1] = O_data[4*j+1];
		}
		fftw_execute(fft);
		for (j = 0; j < P_ANALOG; j++)
		{
			xfreq_re[j] = out[j][0] / (double)P_ANALOG;
			xfreq_im[j] = out[j][1] / (double)P_ANALOG;
		}

		for (j = 0; j < P_ANALOG; j++)
		{
 			in[j][0] = O_data[4*j+2];
			in[j][1] = O_data[4*j+3];
		}
		fftw_execute(fft);
		for (j = 0; j < P_ANALOG; j++)
		{
			yfreq_re[j] = out[j][0] / (double)P_ANALOG;
			yfreq_im[j] = out[j][1] / (double)P_ANALOG;
		}

		//SSFMcheck
		//sprintf(filename2,"SSFM_check_NLprosess_x_FFT_step dist %d m_.csv", i);
		//make_CSV2(P_ANALOG, xfreq_re, xfreq_im, filename2);

		//線形効果
		for (j = 0; j < P_ANALOG; j++)
		{
			//正しい周波数処理が行われるためのインデックス計算
			if (j >= P_ANALOG / 2)
			{
				index = j - P_ANALOG;
			}
			else
			{
				index = j;
			}
		    //線形効果処理をかける角周波数の計算
			omega = 2.0 * PI * df *(double)index;   // 角周波数ωの計算
			omega2 = omega*omega;                   // ω二乗を計算
			DIx = (BETA2 * omega2 / 2.0 - beta1x * omega) * DIST_H;
			DIy = (BETA2 * omega2 / 2.0 - beta1y * omega) * DIST_H;

			x_re = xfreq_re[j];
			x_im = xfreq_im[j];
			y_re = yfreq_re[j];
			y_im = yfreq_im[j];

			xfreq_re[j] = (x_re * cos(DIx) - x_im * sin(DIx)) * exp(-DIST_H*alphaC / 2.0);
			xfreq_im[j] = (x_re * sin(DIx) + x_im * cos(DIx)) * exp(-DIST_H*alphaC / 2.0);
			yfreq_re[j] = (y_re * cos(DIy) - y_im * sin(DIy)) * exp(-DIST_H*alphaC / 2.0);
			yfreq_im[j] = (y_re * sin(DIy) + y_im * cos(DIy)) * exp(-DIST_H*alphaC / 2.0);
		}
		//SSFMcheck
		//sprintf(filename2,"SSFM_check_NLprosess_x_L_step dist %d m_.csv", i);
		//make_CSV2(P_ANALOG, xfreq_re, xfreq_im, filename2);

		//逆フーリエ変換
		for (j = 0; j < P_ANALOG; j++) //逆フーリエ変換する周波数領域信号を代入
		{
			in[j][0] = xfreq_re[j];
			in[j][1] = xfreq_im[j];
		}
		fftw_execute(ifft); //逆フーリエ変換
		for (j = 0; j < P_ANALOG; j++)
		{
			O_data[4*j+0] = out[j][0];
			O_data[4*j+1] = out[j][1];
		}

		for (j = 0; j < P_ANALOG; j++) //逆フーリエ変換する周波数領域信号を代入
		{
			in[j][0] = yfreq_re[j];
			in[j][1] = yfreq_im[j];
		}
		fftw_execute(ifft); //逆フーリエ変換
		for (j = 0; j < P_ANALOG; j++)
		{
			O_data[4*j+2] = out[j][0];
			O_data[4*j+3] = out[j][1];
		}

 		// 偏波回転
		// if((int)now_distance % 8 == 0)//相関長0.1kmに設定(0.1kmごとに偏波回転する)
		// {
  	// 	for(j = 0; j < P_ANALOG; j++)
    //  	{
		//
		// 	x_re = O_data[4*j+0];
		// 	x_im = O_data[4*j+1];
		// 	y_re = O_data[4*j+2];
		// 	y_im = O_data[4*j+3];
		//
		//
		// 	O_data[4*j+0] = x_re * cos(theta[i]) + y_re * sin(theta[i]) * cos(phi[i]) - y_im * sin(theta[i]) * sin(phi[i]);
		// 	O_data[4*j+1] = x_im * cos(theta[i]) + y_im * sin(theta[i]) * cos(phi[i]) + y_re * sin(theta[i]) * sin(phi[i]);
		// 	O_data[4*j+2] = y_re * cos(theta[i]) - x_re * sin(theta[i]) * cos(phi[i]) - x_im * sin(theta[i]) * sin(phi[i]);
		// 	O_data[4*j+3] = y_im * cos(theta[i]) - x_im * sin(theta[i]) * cos(phi[i]) + x_re * sin(theta[i]) * sin(phi[i]);
		//
		//
		// 	/* //スライドのほうの式(Exからどれくらいの位相差があるかをφに設定),どっちも(θとφ)ランダムにしたときにバグる
		// 	analog_re_x[j] = x_re * cos(theta[i]) + y_re * sin(theta[i]) * cos(phi[i]) + y_im * sin(theta[i]) * sin(phi[i]);
		// 	analog_im_x[j] = x_im * cos(theta[i]) + y_re * sin(theta[i]) * sin(phi[i]) + y_im * sin(theta[i]) * cos(phi[i]);
		// 	analog_re_y[j] = -x_re * sin(theta[i]) + y_re * cos(theta[i]) * cos(phi[i]) + y_im * cos(theta[i]) * sin(phi[i]);
		// 	analog_im_y[j] = -x_im * sin(theta[i]) - y_re * cos(theta[i]) * sin(phi[i]) + y_im * cos(theta[i]) * cos(phi[i]); */
		// 	}
		// }

		/* 		//SSFMcheck
		//sprintf(filename2,"SSFM_check_NLprosess_x_IFFT_step dist %d m_.csv", i);
		//make_CSV2(P_ANALOG, analog_re_x, analog_im_x, filename2);

		if(i % 10 == 0)
		{
			printf("SSFM Progress:%.2f m\n", now_distance);
		}    */
	}

	// FFTW3ライブラリの終了処理
	fftw_destroy_plan(fft);
	fftw_destroy_plan(ifft);
	fftw_free(in);
	fftw_free(out);
	free(xfreq_re);free(xfreq_im);
	free(yfreq_re);free(yfreq_im);
	free(deltabeta);
	free(phi);
	free(theta);

}

double nrand()
{
	static int sw = 0;
	static double r1, r2, s;

	if (sw == 0)
	{
		sw = 1;
		do
		{
			r1 = 2.0 * drand48() - 1.0;
			r2 = 2.0 * drand48() - 1.0;
			s = r1 * r1 + r2 * r2;
		} while (s > 1.0 || s == 0.0);
		s = sqrt(-2.0 * log(s) / s);
		return(r1 * s);
	}
	else
	{
		sw = 0;
		return(r2 * s);
	}
}

void EtoO(double *E_ana, double *O_ana){			//パワー調整
	int i;
	double power_dsp_x = 0.0;
	double power_dsp_y = 0.0;
	double power_adjust_x = 0.0;
	double power_adjust_y =0.0;
	for(i=0;i<P_ANALOG;i++){
		power_dsp_x += E_ana[4*i +0] * E_ana[4*i +0] + E_ana[4*i +1] * E_ana[4*i +1];
		power_dsp_y += E_ana[4*i +2] * E_ana[4*i +2] + E_ana[4*i +3] * E_ana[4*i +3];
	}

	power_dsp_x /= (double)P_ANALOG;
	power_dsp_y /= (double)P_ANALOG;
	power_adjust_x = POWER_LAUNCH / (2 * power_dsp_x);
	power_adjust_y = POWER_LAUNCH / (2 * power_dsp_y); //偏波多重のため2で割る


	for(i=0;i<P_ANALOG;i++){
		O_ana[4*i +0] = E_ana[4*i +0] * sqrt(power_adjust_x);
		O_ana[4*i +1] = E_ana[4*i +1] * sqrt(power_adjust_x);
		O_ana[4*i +2] = E_ana[4*i +2] * sqrt(power_adjust_y);
		O_ana[4*i +3] = E_ana[4*i +3] * sqrt(power_adjust_y);
	}

	return ;
}

void Stokes_nomal(double *S_data, double *S_nomal, int num){
	int i;
	double power_rx = 0.0;
	double power_adjust = 0.0;

	// 	for(i=0;i<N_SYMBOL;i++){
	// 		power_rx += S_data[4 * i +0];	//S0合計
	// 	}
	// 	power_rx /= N_SYMBOL;						//S0平均
	//
	// 	power_adjust = 1/power_rx;			//1で規格化
	//
	// 	for(i=0;i<N_SYMBOL;i++){
	// 		S_nomal[4 * i + 0] *= power_adjust;
	// 		S_nomal[4 * i + 1] *= power_adjust;
	// 		S_nomal[4 * i + 2] *= power_adjust;
	// 		S_nomal[4 * i + 3] *= power_adjust;
	// 	}
	for(i=0;i<num;i++){
		//s0
		power_rx = sqrt((S_data[4*i+1]*S_data[4*i+1]) + (S_data[4*i+2]*S_data[4*i+2]) + (S_data[4*i+3]*S_data[4*i+3]));

		S_nomal[4*i+1] = S_data[4*i+1] / power_rx;
		S_nomal[4*i+2] = S_data[4*i+2] / power_rx;
		S_nomal[4*i+3] = S_data[4*i+3] / power_rx;
		S_nomal[4*i+0] = 1;
	}

	return ;
}

void Sample(double *analogData,double *digitData, int n){
	int i = 0;
	int j = 0;

	printf("a: %d\n",n);
	for(i=0;i<P_ANALOG;i++)		//アナログ信号から64点ピックアップ
	{
		if(i%P_PAD==0){
			j = i + n;
			digitData[4 * i/P_PAD +0] = analogData[4 * j + 0];
			digitData[4 * i/P_PAD +1] = analogData[4 * j + 1];
			digitData[4 * i/P_PAD +2] = analogData[4 * j + 2];
			digitData[4 * i/P_PAD +3] = analogData[4 * j + 3];
		}
	}
	return ;
}

void coupler(double *analogData){
	int i;
	for(i=0;i<P_ANALOG;i++){
		analogData[i * 4 + 0] /= sqrt(3.0);
		analogData[i * 4 + 1] /= sqrt(3.0);
		analogData[i * 4 + 2] /= sqrt(3.0);
		analogData[i * 4 + 3] /= sqrt(3.0);
	}
	return;
}

void SHOTnoise(double *analogData, int noise_seed, int flag){

	if(flag == 0){
		return;
	}
	int i;							// ループ用変数
	static int gaussian_seed = 0;	// ガウス分布のシード
	double* gaussian_rand = (double*)malloc(sizeof(double) * P_ANALOG *4);			// ガウス分布の乱数格納
	double ase_value[4];	// ASE雑音の値を格納
	double ase_sigma;				// ASE雑音の標準偏差を格納
	double ase_average;				// ASE雑音の平均値を格納
	double gain;					// リニアスケールの利得
	double myu;						// Box-Muller法で発生させる分布の平均
	double sigma;					// Box-Muller法で発生させる分布の標準偏差

	/* 初期値を代入 */
	ase_average = 0.0;
	myu = 0.0;
	sigma = 1.0;

	/* ガウス分布の乱数生成 */
	srand(gaussian_seed + noise_seed);
	for (i = 0;i < P_ANALOG; i++) gaussian_rand[4*i +0] = gaussian(myu, sigma);
	srand(gaussian_seed + 8192 + noise_seed);
	for (i = 0;i < P_ANALOG; i++) gaussian_rand[4*i +1] = gaussian(myu, sigma);
	srand(gaussian_seed + 16384 + noise_seed);
	for (i = 0;i < P_ANALOG; i++) gaussian_rand[4*i +2] = gaussian(myu, sigma);
	srand(gaussian_seed + 24576 + noise_seed);
	for (i = 0;i < P_ANALOG; i++) gaussian_rand[4*i +3] = gaussian(myu, sigma);


	gaussian_seed++;

	/* ASE雑音の標準偏差を計算 */
	// 全帯域をかける
	ase_sigma = 0.5 * sqrt(1.0 * H_PLANCK * CENTER_FREQ * RATE_BAUD);//単位[√W]
	//両側電力スペクトル(光領域)に雑音を付加するのでボーレート分だけ付加する

	//全帯域をかけるとき、(double)P_ANALOG / T_ALL

	printf("shotnoise = %e\n",ase_sigma);
	printf("gaussian_rand = %f\n",gaussian_rand[4*i+0]);
	printf("ase_average = %f\n",ase_average);

	/* 増幅・ASE雑音の付加 */
	for (i = 0; i < P_ANALOG; i++)
	{
	// ASE雑音の生成
		ase_value[0] = ase_sigma * gaussian_rand[4*i+0] + ase_average;
		ase_value[1] = ase_sigma * gaussian_rand[4*i+1] + ase_average;
		ase_value[2] = ase_sigma * gaussian_rand[4*i+2] + ase_average;
		ase_value[3] = ase_sigma * gaussian_rand[4*i+3] + ase_average;

	// 信号の増幅＆雑音の付加
		analogData[4*i+0] += ase_value[0];//元の信号+雑音
		analogData[4*i+1] += ase_value[1];
		analogData[4*i+2] += ase_value[2];
		analogData[4*i+3] += ase_value[3];
	}

	/* メモリ解放 */
	free(gaussian_rand);
	return ;
}

double gaussian(double mu,double sgm){
	static int flag = 0;
	static double save = 0.0;
	double u_0, u_1, v_0, v_1;

	if(flag == 0){
		u_0 = genRand();
		u_1 = genRand();
		v_0 = mu + sgm * sqrt(-2.0 * log(u_0)) * cos(2.0 * PI * u_1);
		v_1 = mu + sgm * sqrt(-2.0 * log(u_0)) * sin(2.0 * PI * u_1);
		save = v_1;
		flag = 1;
		return v_0;
	}
	else{
		flag = 0;
		return save;
	}
}

double genRand(){
	return ((double)rand()) / (RAND_MAX);
}

void res_sen(double *O_data,double *E_data){
	int i;
	for(i=0;i<P_ANALOG;i++){
		E_data[4 * i + 0] = O_data[4 * i + 0] * sqrt(SENSITIVITY);
		E_data[4 * i + 1] = O_data[4 * i + 1] * sqrt(SENSITIVITY);
		E_data[4 * i + 2] = O_data[4 * i + 2] * sqrt(SENSITIVITY);
		E_data[4 * i + 3] = O_data[4 * i + 3] * sqrt(SENSITIVITY);
	}
	return;
}

void THERMALnoise(double *S_data, int noise_seed, int flag){

	if(flag == 0){
		return;
	}

	int i;
	static int gaussian_seed=0;			// ガウス分布のシード
	double* gaussian_rand = (double*)malloc(sizeof(double)*P_ANALOG *3);	// Sパラメータガウス分布の乱数格納
	double thermal_value[3];				//ASE雑音の値
	double thermal_sigma;						//ASE雑音の標準偏差
	double thermal_ave;							//ASE雑音の平均値を格納
	double gain;			//リニアスケールの利用
	double mu;				//Box-Muller法で発生させる分布の平均
	double sigma;			//Box-Muller法で発生させる分布の標準偏差

	thermal_ave = 0.0;
	mu = 0.0;
	sigma = 1.0;

	srand(gaussian_seed + noise_seed);
	for(i=0;i<P_ANALOG;i++){
		gaussian_rand[3 * i +0] = gaussian(mu,sigma);
	}
	srand(gaussian_seed * 8192 + noise_seed);
	for(i=0;i<P_ANALOG;i++){
		gaussian_rand[3 * i +1] = gaussian(mu,sigma);
	}
	srand(gaussian_seed * 16384 + noise_seed);
	for(i=0;i<P_ANALOG;i++){
		gaussian_rand[3 * i +2] = gaussian(mu,sigma);
	}

	thermal_sigma = sqrt(4.0 + BOLTZMAN * TEMPRATURE * RATE_BAUD / (2.0 * RESISTANCE));
	printf("thermal noise = %e[A]\n",thermal_sigma);

	//ASE雑音
	for(i=0;i<P_ANALOG;i++){
		thermal_value[0] = thermal_sigma * gaussian_rand[3 * i + 0] + thermal_ave;
		thermal_value[1] = thermal_sigma * gaussian_rand[3 * i + 1] + thermal_ave;
		thermal_value[2] = thermal_sigma * gaussian_rand[3 * i + 2] + thermal_ave;

		S_data[3 * i + 1] += thermal_value[0];
		S_data[3 * i + 2] += thermal_value[1];
		S_data[3 * i + 3] += thermal_value[2];

		S_data[3 * i + 0] = sqrt(S_data[3 * i + 1] * S_data[3 * i + 1] + S_data[3 * i + 2] *S_data[3 * i + 2] + S_data[3 * i + 3] * S_data[3 * i + 3]);
	}

	free(gaussian_rand);

	return;
}

int FIR_filter(double *S_aft, double *S_pre){
	int i,j,k;
	double dpwr_x, dpwr_y;

	double h1, h2, h3;

	double ueS11,ueS12,ueS13,ueS21,ueS22,ueS23,ueS31,ueS32,ueS33;

	double* FIR_S = (double*)malloc(sizeof(double) * TR_SYMBOL *3);  //各入力信号

	double *H_11 = (double*)malloc(sizeof(double)*LMS_TAP);		// 各フィルタ係数
	double *H_12 = (double*)malloc(sizeof(double)*LMS_TAP);
	double *H_13 = (double*)malloc(sizeof(double)*LMS_TAP);
	double *H_21 = (double*)malloc(sizeof(double)*LMS_TAP);
	double *H_22 = (double*)malloc(sizeof(double)*LMS_TAP);
	double *H_23 = (double*)malloc(sizeof(double)*LMS_TAP);
	double *H_31 = (double*)malloc(sizeof(double)*LMS_TAP);
	double *H_32 = (double*)malloc(sizeof(double)*LMS_TAP);
	double *H_33 = (double*)malloc(sizeof(double)*LMS_TAP);

	double e_S[3];	//誤差(ストークス)
	double e_E[4];	//誤差
	double ue_S[3] = {0,0,0};
	double ue_E[4];

	for(i = 0; i<LMS_TAP; i++){
		H_11[i] = 0.0;	// 各フィルタ係数
		H_12[i] = 0.0;
		H_13[i] = 0.0;
		H_21[i] = 0.0;
		H_22[i] = 0.0;
		H_23[i] = 0.0;
		H_31[i] = 0.0;
		H_32[i] = 0.0;
		H_33[i] = 0.0;
	}

	for(i = 0; i<TR_SYMBOL; i++){
		FIR_S[3*i +0] = 0.0;	// 各入力信号
		FIR_S[3*i +1] = 0.0;
		FIR_S[3*i +2] = 0.0;
	}

	//FIR応答処理

	for(i=0; i<TR_SYMBOL; i++){
		FIR_S[0] = S_aft[4*i +1];	// 各入力信号
		FIR_S[1] = S_aft[4*i +2];
		FIR_S[2] = S_aft[4*i +3];

		h1 = 0.0;
		h2 = 0.0;
		h3 = 0.0;

		for(j = 0; j<LMS_TAP; j++){
			h1 += (H_11[j] * FIR_S[4*j +1]) + (H_21[j] * FIR_S[4*j +2]) + (H_31[j] * FIR_S[4*j +3]);//最初は0
			h2 += (H_12[j] * FIR_S[4*j +1]) + (H_22[j] * FIR_S[4*j +2]) + (H_32[j] * FIR_S[4*j +3]);
			h3 += (H_13[j] * FIR_S[4*j +1]) + (H_23[j] * FIR_S[4*j +2]) + (H_33[j] * FIR_S[4*j +3]);
		}

		if(i>= SFT_SYMBOL){
			e_S[0] = S_pre[4*(i - SFT_SYMBOL)] -h1;
			e_S[1] = S_pre[4*(i - SFT_SYMBOL)] -h2;
			e_S[2] = S_pre[4*(i - SFT_SYMBOL)] -h3;

			ue_S[0] = e_S[0] * LMS_STEP;
			ue_S[1] = e_S[1] * LMS_STEP;
			ue_S[2] = e_S[2] * LMS_STEP;

			for(j=0; j<LMS_TAP; j++){
				ueS11 = FIR_S[3*j +0] * ue_S[0];//μe(n)×S(n)
				ueS21 = FIR_S[3*j +1] * ue_S[0];
				ueS31 = FIR_S[3*j +2] * ue_S[0];
				ueS12 = FIR_S[3*j +0] * ue_S[1];
				ueS22 = FIR_S[3*j +1] * ue_S[1];
				ueS32 = FIR_S[3*j +2] * ue_S[1];
				ueS13 = FIR_S[3*j +0] * ue_S[2];
				ueS23 = FIR_S[3*j +1] * ue_S[2];
				ueS33 = FIR_S[3*j +2] * ue_S[2];

				ueS11 /= S0;//規格化
				ueS12 /= S0;
				ueS13 /= S0;
				ueS21 /= S0;
				ueS22 /= S0;
				ueS23 /= S0;
				ueS31 /= S0;
				ueS32 /= S0;
				ueS33 /= S0;


				H_11[j] += ueS11;//フィルタ係数に新たなものを加算代入、H_11に数値が入る
				H_12[j] += ueS12;
				H_13[j] += ueS13;
				H_21[j] += ueS21;
				H_22[j] += ueS22;
				H_23[j] += ueS23;
				H_31[j] += ueS31;
				H_32[j] += ueS32;
				H_33[j] += ueS33;
			}

		}
		for(j=LMS_TAP - 1; j>0; j--){
			FIR_S[3*j+0] = FIR_S[3*(j-1)+0];
			FIR_S[3*j+1] = FIR_S[3*(j-1)+1];
			FIR_S[3*j+2] = FIR_S[3*(j-1)+2];
		}
		}

		//最終的にTR_SYMBOLをすべて回した上でこちらに来る

		// 係数を用いて、ペイロードのフィルタリングを適用
		// 本信号のみに適応
		for (i = TR_SYMBOL; i < N_SYMBOL; i++)
		{
			FIR_S[0] = S_aft[4*i +1];
			FIR_S[1] = S_aft[4*i +2];
			FIR_S[2] = S_aft[4*i +3];
			h1 = 0.0;
			h2 = 0.0;
			h3 = 0.0;

			for (j = 0; j < LMS_TAP; j++)
			{
				h1 += (H_11[j] * FIR_S[3*j+0]) + (H_21[j] * FIR_S[3*j+1]) + (H_31[j] * FIR_S[3*j+2]);
				h2 += (H_12[j] * FIR_S[3*j+0]) + (H_22[j] * FIR_S[3*j+1]) + (H_32[j] * FIR_S[3*j+2]);
				h3 += (H_13[j] * FIR_S[3*j+0]) + (H_23[j] * FIR_S[3*j+1]) + (H_33[j] * FIR_S[3*j+2]);

			}
			S_aft[4*i +1] = h1;
			S_aft[4*i +2] = h2;
			S_aft[4*i +3] = h3;

			//データシフト

			for(j=LMS_TAP - 1; j>0; j--){
				FIR_S[3*j+0] = FIR_S[3*(j-1)+0];
				FIR_S[3*j+1] = FIR_S[3*(j-1)+1];
				FIR_S[3*j+2] = FIR_S[3*(j-1)+2];
			}
	}

	// メモリ解放
	free(H_11);
	free(H_12);
	free(H_13);
	free(H_21);
	free(H_22);
	free(H_23);
	free(H_31);
	free(H_32);
	free(H_33);

	return 0;

}

int DEL_TSYMBOL(double *S_data){
	int i;

	for(i=0; i<DATA_SYMBOL; i++){
		S_data[4*i +0] = S_data[4*(i+TR_SYMBOL) +0];
		S_data[4*i +1] = S_data[4*(i+TR_SYMBOL) +1];
		S_data[4*i +2] = S_data[4*(i+TR_SYMBOL) +2];
		S_data[4*i +3] = S_data[4*(i+TR_SYMBOL) +3];
	}

	for(i=DATA_SYMBOL;i<N_SYMBOL;i++){
		S_data[4*i +0] = 0.0;
		S_data[4*i +1] = 0.0;
		S_data[4*i +2] = 0.0;
		S_data[4*i +3] = 0.0;
	}
	return 0;
}

int Stokes_mod(int *sequence, double *S_data){
	int i,j;
	int *number = (int*)malloc(sizeof(int)* DATA_SYMBOL);

	for(i=0; i<DATA_SYMBOL; i++) number[i] = 0;

	for(i=0; i< DATA_SYMBOL; i++){

	}

	return 0;
}

void catch_Symbol(double *S_data, int bitStream, int num){
	int i;
	int SymbolStream[num];			//0~8のシンボル列
	int binary;

	for(i=0 ; i<num ; i++){
			if(S_data[4*i +1]/sqrt(3) +S_data[4*i +2]/sqrt(3) + S_data[4*i +3]/sqrt(3) > CONST_D*(-1)) SymbolStream[i] = 0;

			else if(S_data[4*i +1]/sqrt(3) +S_data[4*i +2]*-1/sqrt(3) + S_data[4*i +3]/sqrt(3) > CONST_D*(-1)) SymbolStream[i] = 1;

			else if(S_data[4*i +1]*-1/sqrt(3) +S_data[4*i +2]/sqrt(3) + S_data[4*i +3]/sqrt(3) > CONST_D*(-1)) SymbolStream[i] = 2;

			else if(S_data[4*i +1]*-1/sqrt(3) +S_data[4*i +2]*-1/sqrt(3) + S_data[4*i +3]/sqrt(3) > CONST_D*(-1)) SymbolStream[i] = 3;

			else if(S_data[4*i +1]/sqrt(3) +S_data[4*i +2]/sqrt(3) + S_data[4*i +3]*-1/sqrt(3) > CONST_D*(-1)) SymbolStream[i] = 4;

			else if(S_data[4*i +1]/sqrt(3) +S_data[4*i +2]*-1/sqrt(3) + S_data[4*i +3]-1/sqrt(3) > CONST_D*(-1)) SymbolStream[i] = 5;

			else if(S_data[4*i +1]*-1/sqrt(3) +S_data[4*i +2]/sqrt(3) + S_data[4*i +3]*-1/sqrt(3) > CONST_D*(-1)) SymbolStream[i] = 6;

			else if(S_data[4*i +1]*-1/sqrt(3) +S_data[4*i +2]*-1/sqrt(3) + S_data[4*i +3]*-1/sqrt(3) > CONST_D*(-1)) SymbolStream[i] = 7;

			else if(S_data[4*i +1]*0 +S_data[4*i +2]*2/(1+sqrt(5))/sqrt(3) + S_data[4*i +3]*(1+sqrt(5))/2 /sqrt(3) > CONST_D*(-1)) SymbolStream[i] = 8;

			else if(S_data[4*i +1]*0 +S_data[4*i +2]*2/(1+sqrt(5))/sqrt(3) + S_data[4*i +3]*-1*(1+sqrt(5))/2 /sqrt(3) > CONST_D*(-1)) SymbolStream[i] = 9;

			else if(S_data[4*i +1]*0 +S_data[4*i +2]*-2/(1+sqrt(5))/sqrt(3) + S_data[4*i +3]*(1+sqrt(5))/2 /sqrt(3) > CONST_D*(-1)) SymbolStream[i] = 10;

			else if(S_data[4*i +1]*0 +S_data[4*i +2]*-2/(1+sqrt(5))/sqrt(3) + S_data[4*i +3]*-1*(1+sqrt(5))/2 /sqrt(3) > CONST_D*(-1)) SymbolStream[i] = 11;

			else if(S_data[4*i +1]*(1+sqrt(5))/2 /sqrt(3) + S_data[4*i +2]*0 + S_data[4*i +3]*2/(1+sqrt(5)) /sqrt(3) > CONST_D*(-1)) SymbolStream[i] = 12;

			else if(S_data[4*i +1]*(1+sqrt(5))/2 /sqrt(3) + S_data[4*i +2]*0 + S_data[4*i +3]*-2/(1+sqrt(5)) /sqrt(3) > CONST_D*(-1)) SymbolStream[i] = 13;

			else if(S_data[4*i +1]*-1*(1+sqrt(5))/2 /sqrt(3) + S_data[4*i +2]*0 + S_data[4*i +3]*2/(1+sqrt(5)) /sqrt(3) > CONST_D*(-1)) SymbolStream[i] = 14;

			else if(S_data[4*i +1]*-1*(1+sqrt(5))/2 /sqrt(3) + S_data[4*i +2]*0 + S_data[4*i +3]*-2/(1+sqrt(5)) /sqrt(3) > CONST_D*(-1)) SymbolStream[i] = 15;

			else {
				SymbolStream[i] = 0;
				printf("catch Symbol error!! \n S0: %lf, S1: %lf, S2:% lf, S3: %lf\n",S_data[4*i+0],S_data[4*i+1],S_data[4*i+2],S_data[4*i+3]);
			}

		}

		//シンボル列からビット列
	for(i=0;i<num;i++){
		for(j=0;j<4;j++){
  	 	bitStream[4*i+(3-j)] = ( SymbolStream[i] % 2 ) ;
  		SymbolStream[i] /= 2;
		}
	}

	return;
}
