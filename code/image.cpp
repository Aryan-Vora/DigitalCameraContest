#include <algorithm>
#include <cmath>
#include <iostream>
#include <math.h>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <utility>
#include <vector>

int pic[2460][3360];
void OB();
void WB();
void color_interpolation();
void colormatrix();
void apply_gamma(double);
void color();
void edge_enhance();
void rgb2bmp(char *);

void changeVibrance(double change, double sat);
void adaptiveVibrance(double change, int window);
void getColor();
void NearestNeighborInterpolation();
void sharpen(double change);
void changeSaturation(double change);
void RGB_2HSV_2RGB(double change);

void CFA_to_krkb();
int truncate(double value);

void adjustContrast(double contrast);
int r[2460][3360];
int g[2460][3360];
int b[2460][3360];
int e[2460][3360];
int kr[2460][3360];
int kb[2460][3360];

unsigned short int gamma_value[1024] = {
    0,   1,   3,   5,   8,   9,   11,  12,  13,  14,  16,  17,  18,  20,  21,
    22,  24,  25,  26,  28,  29,  30,  32,  33,  34,  36,  37,  38,  39,  41,
    42,  43,  45,  46,  47,  49,  50,  51,  53,  54,  55,  57,  58,  59,  61,
    62,  63,  64,  66,  67,  68,  70,  71,  72,  74,  75,  76,  78,  79,  80,
    82,  83,  84,  85,  87,  87,  88,  89,  89,  90,  91,  91,  92,  93,  93,
    94,  95,  95,  96,  97,  97,  98,  99,  99,  100, 101, 101, 102, 103, 103,
    104, 105, 105, 106, 107, 107, 108, 108, 109, 110, 110, 111, 112, 112, 113,
    114, 114, 115, 116, 116, 117, 118, 118, 119, 120, 120, 121, 122, 122, 123,
    124, 124, 125, 126, 126, 127, 128, 128, 129, 129, 130, 130, 130, 130, 131,
    131, 131, 132, 132, 132, 133, 133, 133, 134, 134, 134, 135, 135, 135, 136,
    136, 136, 137, 137, 137, 138, 138, 138, 139, 139, 139, 140, 140, 140, 141,
    141, 141, 142, 142, 142, 143, 143, 143, 144, 144, 144, 145, 145, 145, 146,
    146, 146, 147, 147, 147, 148, 148, 148, 149, 149, 149, 150, 150, 150, 151,
    151, 151, 152, 152, 152, 153, 153, 153, 154, 154, 154, 154, 155, 155, 155,
    156, 156, 156, 157, 157, 157, 158, 158, 158, 159, 159, 159, 160, 160, 160,
    161, 161, 161, 162, 162, 162, 163, 163, 163, 164, 164, 164, 165, 165, 165,
    166, 166, 166, 167, 167, 167, 168, 168, 168, 169, 169, 169, 170, 170, 170,
    171, 171, 171, 171, 171, 172, 172, 172, 172, 172, 172, 173, 173, 173, 173,
    173, 173, 174, 174, 174, 174, 174, 174, 175, 175, 175, 175, 175, 175, 176,
    176, 176, 176, 176, 176, 176, 177, 177, 177, 177, 177, 177, 178, 178, 178,
    178, 178, 178, 179, 179, 179, 179, 179, 179, 180, 180, 180, 180, 180, 180,
    181, 181, 181, 181, 181, 181, 182, 182, 182, 182, 182, 182, 183, 183, 183,
    183, 183, 183, 184, 184, 184, 184, 184, 184, 185, 185, 185, 185, 185, 185,
    186, 186, 186, 186, 186, 186, 187, 187, 187, 187, 187, 187, 188, 188, 188,
    188, 188, 188, 189, 189, 189, 189, 189, 189, 189, 190, 190, 190, 190, 190,
    190, 191, 191, 191, 191, 191, 191, 192, 192, 192, 192, 192, 192, 193, 193,
    193, 193, 193, 193, 194, 194, 194, 194, 194, 194, 195, 195, 195, 195, 195,
    195, 196, 196, 196, 196, 196, 196, 197, 197, 197, 197, 197, 197, 198, 198,
    198, 198, 198, 198, 199, 199, 199, 199, 199, 199, 200, 200, 200, 200, 200,
    200, 201, 201, 201, 201, 201, 201, 201, 202, 202, 202, 202, 202, 202, 203,
    203, 203, 203, 203, 203, 204, 204, 204, 204, 204, 204, 205, 205, 205, 205,
    205, 205, 206, 206, 206, 206, 206, 206, 207, 207, 207, 207, 207, 207, 208,
    208, 208, 208, 208, 208, 209, 209, 209, 209, 209, 209, 210, 210, 210, 210,
    210, 210, 211, 211, 211, 211, 211, 211, 212, 212, 212, 212, 212, 212, 213,
    213, 213, 213, 213, 213, 213, 213, 213, 213, 214, 214, 214, 214, 214, 214,
    214, 214, 214, 214, 214, 214, 215, 215, 215, 215, 215, 215, 215, 215, 215,
    215, 215, 215, 216, 216, 216, 216, 216, 216, 216, 216, 216, 216, 216, 216,
    217, 217, 217, 217, 217, 217, 217, 217, 217, 217, 217, 217, 217, 218, 218,
    218, 218, 218, 218, 218, 218, 218, 218, 218, 218, 219, 219, 219, 219, 219,
    219, 219, 219, 219, 219, 219, 219, 220, 220, 220, 220, 220, 220, 220, 220,
    220, 220, 220, 220, 221, 221, 221, 221, 221, 221, 221, 221, 221, 221, 221,
    221, 222, 222, 222, 222, 222, 222, 222, 222, 222, 222, 222, 222, 223, 223,
    223, 223, 223, 223, 223, 223, 223, 223, 223, 223, 223, 224, 224, 224, 224,
    224, 224, 224, 224, 224, 224, 224, 224, 225, 225, 225, 225, 225, 225, 225,
    225, 225, 225, 225, 225, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226,
    226, 226, 227, 227, 227, 227, 227, 227, 227, 227, 227, 227, 227, 227, 228,
    228, 228, 228, 228, 228, 228, 228, 228, 228, 228, 228, 229, 229, 229, 229,
    229, 229, 229, 229, 229, 229, 229, 229, 229, 230, 230, 230, 230, 230, 230,
    230, 230, 230, 230, 230, 230, 231, 231, 231, 231, 231, 231, 231, 231, 231,
    231, 231, 231, 232, 232, 232, 232, 232, 232, 232, 232, 232, 232, 232, 232,
    233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 234, 234, 234,
    234, 234, 234, 234, 234, 234, 234, 234, 234, 235, 235, 235, 235, 235, 235,
    235, 235, 235, 235, 235, 235, 235, 236, 236, 236, 236, 236, 236, 236, 236,
    236, 236, 236, 236, 237, 237, 237, 237, 237, 237, 237, 237, 237, 237, 237,
    237, 238, 238, 238, 238, 238, 238, 238, 238, 238, 238, 238, 238, 239, 239,
    239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 240, 240, 240, 240, 240,
    240, 240, 240, 240, 240, 240, 240, 241, 241, 241, 241, 241, 241, 241, 241,
    241, 241, 241, 241, 241, 242, 242, 242, 242, 242, 242, 242, 242, 242, 242,
    242, 242, 243, 243, 243, 243, 243, 243, 243, 243, 243, 243, 243, 243, 244,
    244, 244, 244, 244, 244, 244, 244, 244, 244, 244, 244, 245, 245, 245, 245,
    245, 245, 245, 245, 245, 245, 245, 245, 246, 246, 246, 246, 246, 246, 246,
    246, 246, 246, 246, 246, 247, 247, 247, 247, 247, 247, 247, 247, 247, 247,
    247, 247, 247, 248, 248, 248, 248, 248, 248, 248, 248, 248, 248, 248, 248,
    249, 249, 249, 249, 249, 249, 249, 249, 249, 249, 249, 249, 250, 250, 250,
    250, 250, 250, 250, 250, 250, 250, 250, 250, 251, 251, 251, 251, 251, 251,
    251, 251, 251, 251, 251, 251, 252, 252, 252, 252, 252, 252, 252, 252, 252,
    252, 252, 252, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253,
    254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 255, 255, 255,
    255, 255, 255, 255};

void adjustContrast(double contrast) {
  int i, j;
  double factor = (259 * (contrast + 255)) / (255 * (259 - contrast));

  for (i = 0; i < 2460; i++) {
    for (j = 0; j < 3360; j++) {
      r[i][j] = truncate(factor * (r[i][j] - 128) + 128);
      g[i][j] = truncate(factor * (g[i][j] - 128) + 128);
      b[i][j] = truncate(factor * (b[i][j] - 128) + 128);
    }
  }
  printf("Contrast adjustment ---OK!\n");
}

int truncate(double value) {
  if (value < 0)
    return 0;
  if (value > 255)
    return 255;
  return (int)value;
}

void OB() {
  int ob_sum = 0;
  int ob_avg = 0;
  int count = 0, i, j;

  // OB for E region
  ob_sum = 0;
  count = 0;

  for (i = 0; i < 2460; i = i + 2) {
    for (j = 0; j < 5; j = j + 2) {
      ob_sum = ob_sum + pic[i][j];
      count++;
    }
    for (j = 3294; j < 3343; j = j + 2) {
      ob_sum = ob_sum + pic[i][j];
      count++;
    }
  }

  ob_avg = ob_sum / count;

  for (i = 0; i < 2460; i = i + 2) {
    for (j = 6; j < 3293; j = j + 2) {
      if (pic[i][j] > ob_avg)
        pic[i][j] = pic[i][j] - ob_avg;
      else
        pic[i][j] = 0;
    }
  }

  ob_sum = 0;
  count = 0;

  // OB for R region
  for (i = 0; i < 2460; i = i + 2) {
    for (j = 1; j < 5; j = j + 2) {
      ob_sum = ob_sum + pic[i][j];
      count++;
    }
    for (j = 3293; j < 3343; j = j + 2) {
      ob_sum = ob_sum + pic[i][j];
      count++;
    }
  }
  ob_avg = ob_sum / count;

  for (i = 0; i < 2460; i = i + 2) {
    for (j = 5; j < 3293; j = j + 2) {
      if (pic[i][j] > ob_avg)
        pic[i][j] = pic[i][j] - ob_avg;
      else
        pic[i][j] = 0;
    }
  }

  // OB for B region
  ob_sum = 0;
  count = 0;

  for (i = 1; i < 2460; i = i + 2) {
    for (j = 0; j < 5; j = j + 2) {
      ob_sum = ob_sum + pic[i][j];
      count++;
    }
    for (j = 3294; j < 3343; j = j + 2) {
      ob_sum = ob_sum + pic[i][j];
      count++;
    }
  }

  ob_avg = ob_sum / count;

  for (i = 1; i < 2460; i = i + 2) {
    for (j = 6; j < 3293; j = j + 2) {
      if (pic[i][j] > ob_avg)
        pic[i][j] = pic[i][j] - ob_avg;
      else
        pic[i][j] = 0;
    }
  }

  // OB for G region
  ob_sum = 0;
  count = 0;

  for (i = 1; i < 2460; i = i + 2) {
    for (j = 1; j < 5; j = j + 2) {
      ob_sum = ob_sum + pic[i][j];
      count++;
    }
    for (j = 3293; j < 3343; j = j + 2) {
      ob_sum = ob_sum + pic[i][j];
      count++;
    }
  }

  ob_avg = ob_sum / count;

  for (i = 1; i < 2460; i = i + 2) {
    for (j = 5; j < 3293; j = j + 2) {
      if (pic[i][j] > ob_avg)
        pic[i][j] = pic[i][j] - ob_avg;
      else
        pic[i][j] = 0;
    }
  }

  printf("Optical Black ---OK!\n");

} // OB() end

void WB() {
  int i, j;

  double r_hist[16384], g_hist[16384], b_hist[16384]; // 16384=2^14
  double r_accu[16384], g_accu[16384], b_accu[16384];

  for (i = 0; i < 16384; i++) {
    // initialization histogram
    r_hist[i] = 0;
    g_hist[i] = 0;
    b_hist[i] = 0;
  }

  for (i = 0; i < 2460; i += 2) {
    for (j = 5; j < 3293; j += 2) {
      r_hist[pic[i][j]]++;           // R
      b_hist[pic[i + 1][j + 1]]++;   // B
      g_hist[pic[i + 1][j]]++;       // G
      g_hist[pic[i][j + 1] / 2]++;   // E
      r_hist[(pic[i][j + 1] / 3)]++; // E
      b_hist[(pic[i][j + 1] / 2)]++; // E
    }
  }

  // histogram accumulate �ֿn
  r_accu[0] = r_hist[0];
  g_accu[0] = g_hist[0];
  b_accu[0] = b_hist[0];

  for (i = 1; i < 16384; i++) {
    r_accu[i] = r_hist[i] + r_accu[i - 1];
    g_accu[i] = g_hist[i] + g_accu[i - 1];
    b_accu[i] = b_hist[i] + b_accu[i - 1];
  }

  // normalization
  for (i = 0; i < 16384; i++) {
    r_hist[i] = ((float)r_accu[i] / 1997568) *
                16384; // 2066400 = 3360*2460/4, 16384 = 2^14
    g_hist[i] = ((float)g_accu[i] / 1997568) *
                16384; // 2066400 = 3360*2460/4, 16384 = 2^14
    b_hist[i] = ((float)b_accu[i] / 1997568) *
                16384; // 2066400 = 3360*2460/4, 16384 = 2^14
  }

  /* White Point Detection */
  float r_equal, g_equal, b_equal;
  float Y_hist, Cr_hist, Cb_hist;
  float Y_bright = 0.0, Cr_bright = 192.0,
        Cb_bright = 192.0; // probably white point 192
  double Y_avg_hist = 0.0, Cr_avg_hist1 = 0.0, Cb_avg_hist1 = 0.0,
         Cr_avg_hist2 = 0.0, Cb_avg_hist2 = 0.0;
  unsigned int counterCb1 = 0, counterCb2 = 0, counterCr1 = 0, counterCr2 = 0;

  for (i = 0; i < 2460; i += 2) {
    for (j = 5; j < 3293; j += 2) {
      r_equal = r_hist[pic[i][j]];         // r histogram equalization
      b_equal = b_hist[pic[i + 1][j + 1]]; // b histogram equalization
      g_equal = g_hist[pic[i + 1][j]];     // g histogram equalization

      Y_hist = (0.299 * r_equal) + (0.587 * g_equal) + (0.114 * b_equal);
      Cr_hist = (0.500 * r_equal) - (0.4187 * g_equal) - (0.0813 * b_equal);
      Cb_hist = (-0.1687 * r_equal) - (0.3313 * g_equal) + (0.500 * b_equal);

      if (Y_hist >= 12000 && (Cr_hist <= 192 && Cr_hist >= -192 &&
                              (Cb_hist <= 192 && Cb_hist >= -192))) {
        if (Y_hist > Y_bright) // the bigger the better
          Y_bright = Y_hist;
        if (fabs(Cr_hist) < fabs(Cr_bright)) // the smaller the better
          Cr_bright = Cr_hist;
        if (fabs(Cb_hist) < fabs(Cb_bright)) // the smaller the better
          Cb_bright = Cb_hist;

        Y_avg_hist += Y_hist;

        if (Cr_hist > 0) {
          Cr_avg_hist1 += Cr_hist;
          counterCr1++;
        } else {
          Cr_avg_hist2 += Cr_hist;
          counterCr2++;
        }

        if (Cb_hist > 0) {
          Cb_avg_hist1 += Cb_hist;
          counterCb1++;
        } else {
          Cb_avg_hist2 += Cb_hist;
          counterCb2++;
        }
      }
    }
  }

  if (counterCr1 + counterCr2 + counterCb1 + counterCb2 == 0) {
    printf("================error message============\n");
    printf("error in varsha awb count awb white point\n");
    printf("=========================================\n");
    // Y_bright=11000;
    // Cr_bright=-400;
    // Cb_bright=-400;
  } else {
    Y_avg_hist /= (counterCr1 + counterCr2 + counterCb1 + counterCb2);
    Cr_avg_hist1 /= counterCr1;
    Cr_avg_hist2 /= counterCr2;
    Cb_avg_hist1 /= counterCb1;
    Cb_avg_hist2 /= counterCb2;
  }

  // select the minimum and maximum between
  // brightest of histtogram in Y Cr Cb
  // average of histogram in Y Cr Cb
  // as Y_low,Cr_low,Cb_low Y_up,Cr_up,Cb_up
  float Y_low, Cr_low, Cb_low, Y_up, Cr_up, Cb_up;

  if (Y_bright > Y_avg_hist) {
    Y_up = Y_bright;
    Y_low = Y_avg_hist;
  } else {
    Y_low = Y_bright;
    Y_up = Y_avg_hist;
  }

  if (Cr_bright > Cr_avg_hist1) {
    Cr_up = Cr_bright;
    Cr_low = Cr_avg_hist2;
  } else if (Cr_bright < Cr_avg_hist2) {
    Cr_low = Cr_bright;
    Cr_up = Cr_avg_hist1;
  } else {
    Cr_up = Cr_avg_hist1;
    Cr_low = Cr_avg_hist2;
  }

  if (Cb_bright > Cb_avg_hist1) {
    Cb_up = Cb_bright;
    Cb_low = Cb_avg_hist2;
  } else if (Cb_bright < Cb_avg_hist2) {
    Cb_low = Cb_bright;
    Cb_up = Cb_avg_hist1;
  } else {
    Cb_up = Cb_avg_hist1;
    Cb_low = Cb_avg_hist2;
  }

  /* find the reference white point */
  long int R_white = 0, G_white = 0, B_white = 0;
  double R_avg = 0, G_avg = 0, B_avg = 0;
  float R_hist_avg = 0, G_hist_avg = 0, B_hist_avg = 0;
  unsigned int counter = 0;

  for (i = 0; i < 2460; i += 2) {
    for (j = 5; j < 3293; j += 2) {
      // calculate the original R G B
      R_avg += pic[i][j];         // r
      B_avg += pic[i + 1][j + 1]; // b
      G_avg += pic[i + 1][j];     // g

      // calculate the R G B after histogram equalization
      r_equal = r_hist[pic[i][j]];         // r histogram equalization
      b_equal = b_hist[pic[i + 1][j + 1]]; // b histogram equalization
      g_equal = g_hist[pic[i + 1][j]];     // g histogram equalization
      R_hist_avg += r_equal;
      G_hist_avg += g_equal;
      B_hist_avg += b_equal;

      Y_hist = (0.299 * r_equal) + (0.587 * g_equal) + (0.114 * b_equal);
      Cr_hist = (0.500 * r_equal) - (0.4187 * g_equal) - (0.0813 * b_equal);
      Cb_hist = (-0.1687 * r_equal) - (0.3313 * g_equal) + (0.500 * b_equal);

      // find the pixel which Y_up<=Y_hist>=Y_low and Cr_up<=Cr_hist>=Cr_low and
      // Cb_up<=Cb_hist>=Cb_low
      if ((Y_hist >= Y_low && Y_hist <= Y_up) &&
          (Cr_hist >= Cr_low && Cr_hist <= Cr_up) &&
          (Cb_hist >= Cb_low && Cb_hist <= Cb_up)) {
        R_white += pic[i][j];         // r
        B_white += pic[i + 1][j + 1]; // b
        G_white += pic[i + 1][j];     // g
        counter++;
      }
    }
  }

  if (counter != 0) {
    // find out the reference white point
    R_avg /= 2066400;
    G_avg /= 2066400;
    B_avg /= 2066400;

    R_hist_avg /= 2066400;
    G_hist_avg /= 2066400;
    B_hist_avg /= 2066400;

    R_white /= counter;
    G_white /= counter;
    B_white /= counter;

    /* White Balance Adjustment */
    // calculate the scale factor from reference white point
    float R_scale, G_scale, B_scale;
    R_scale = ((0.299 * R_white) + (0.587 * G_white) + (0.114 * B_white)) /
              R_white; // Y_white/R_white
    G_scale = ((0.299 * R_white) + (0.587 * G_white) + (0.114 * B_white)) /
              G_white; // Y_white/G_white
    B_scale = ((0.299 * R_white) + (0.587 * G_white) + (0.114 * B_white)) /
              B_white; // Y_white/B_white

    // calculate the scale factor from gray world assumption
    float R_GWA, G_GWA, B_GWA;
    R_GWA = ((0.299 * R_avg) + (0.587 * G_avg) + (0.114 * B_avg)) /
            R_avg; // Y_avg/R_avg
    G_GWA = ((0.299 * R_avg) + (0.587 * G_avg) + (0.114 * B_avg)) /
            G_avg; // Y_avg/G_avg
    B_GWA = ((0.299 * R_avg) + (0.587 * G_avg) + (0.114 * B_avg)) /
            B_avg; // Y_avg/B_avg

    // calculate the color cast
    float R_factor, G_factor, B_factor;
    if (B_hist_avg + 192 >= G_hist_avg && B_hist_avg > R_hist_avg) {
      // blusih
      R_factor = R_scale;
      G_factor = G_scale;
      B_factor = B_GWA;
      printf("bluish\n");
    } else if (G_hist_avg + 192 > R_hist_avg && R_hist_avg > B_hist_avg) {
      // greenish
      R_factor = R_scale;
      G_factor = G_GWA;
      B_factor = B_scale;
      printf("greenish\n");
    } else if (R_hist_avg > G_hist_avg && G_hist_avg > B_hist_avg) {
      // reddish
      R_factor = R_GWA;
      G_factor = R_scale;
      B_factor = B_scale;
      printf("redish\n");
    } else {
      R_factor = R_scale;
      G_factor = G_scale;
      B_factor = B_scale;
      printf("normal\n");
    }

    // white balance gain factor
    //		cout << R_factor << endl;
    //		cout << G_factor << endl;
    //		cout << B_factor << endl;

    if (!(R_factor == 1 && G_factor == 1 && B_factor == 1)) {
      for (i = 0; i < 2460; i += 2) {
        for (j = 0; j < 3360; j += 2) {
          // calculate the original R G B
          if (pic[i][j + 1] * R_factor < 0)
            pic[i][j + 1] = 0;
          else
            pic[i][j + 1] *= R_factor; // r

          if (pic[i + 1][j] * B_factor < 0)
            pic[i + 1][j] = 0;
          else
            pic[i + 1][j] *= B_factor; // b

          if (pic[i + 1][j + 1] * G_factor < 0)
            pic[i + 1][j + 1] = 0;
          else
            pic[i + 1][j + 1] *= G_factor; // g
        }
      }
    }
  } else {
    // not find out the reference white point
    printf("@_@ not do the varsha white balance... @_@\n");
  }

  printf("White balance ---OK!\n");

} // WB() end

void color_interpolation() {
  int i, j;

  // R region
  for (i = 0; i < 2460; i = i + 2) {
    for (j = 5; j < 3293; j = j + 2) {
      r[i][j] = pic[i][j];
      g[i][j] = (pic[i - 1][j] + pic[i + 1][j]) / 2;
      b[i][j] = (pic[i - 1][j - 1] + pic[i - 1][j + 1] + pic[i + 1][j - 1] +
                 pic[i + 1][j + 1]) /
                4;
    }
  }

  // G region
  for (i = 1; i < 2460; i = i + 2) {
    for (j = 5; j < 3293; j = j + 2) {
      g[i][j] = pic[i][j];
      r[i][j] = (pic[i - 1][j] + pic[i + 1][j]) / 2;
      b[i][j] = (pic[i][j - 1] + pic[i][j + 1]) / 2;
    }
  }

  // E region
  for (i = 0; i < 2460; i = i + 2) {
    for (j = 6; j < 3293; j = j + 2) {
      b[i][j] = (pic[i - 1][j] + pic[i + 1][j]) / 2;
      g[i][j] = (pic[i - 1][j - 1] + pic[i - 1][j + 1] + pic[i + 1][j - 1] +
                 pic[i + 1][j + 1]) /
                4;
      r[i][j] = (pic[i][j - 1] + pic[i][j + 1]) / 2;
    }
  }

  // B region
  for (i = 1; i < 2460; i = i + 2) {
    for (j = 6; j < 3293; j = j + 2) {
      b[i][j] = pic[i][j];
      r[i][j] = (pic[i - 1][j - 1] + pic[i - 1][j + 1] + pic[i + 1][j - 1] +
                 pic[i + 1][j + 1]) /
                4;
      g[i][j] = (pic[i][j - 1] + pic[i][j + 1]) / 2;
    }
  }

  printf("Color interpolation ---OK!\n");

} // color_interpolation() end

void colormatrix() {
  //	[R']=[ 1+sa+sb   -sa       -sb     ][R]
  //	[G']=[ -sc       1+sc+sd   -sd     ][G]
  //	[B']=[ -se       -sf       1+se+sf ][B]

  unsigned short int i, j;
  double temp_r, temp_g, temp_b;

  double sa = 0.265625, sb = 0.17;
  double sc = 0.21, sd = 0.21;
  double se = 0.19, sf = 0.19;

  for (i = 0; i < 2460; i++) {
    for (j = 0; j < 3360; j++) {
      temp_r = r[i][j];
      temp_g = g[i][j];
      temp_b = b[i][j];

      r[i][j] =
          (double)((6.3329) * temp_r + (-2.7900) * temp_g + (1.5002) * temp_b);

      if (r[i][j] > 16360) // 1022.5*16
        r[i][j] = 16360;
      else if (r[i][j] < 0)
        r[i][j] = 0;

      g[i][j] =
          (double)((0.1362) * temp_r + (2.6597) * temp_g + (0.0751) * temp_b);

      if (g[i][j] > 16360) // 1022.5*16
        g[i][j] = 16360;
      else if (g[i][j] < 0)
        g[i][j] = 0;

      b[i][j] =
          (double)((1.2656) * temp_r + (-1.5919) * temp_g + (4.7861) * temp_b);

      if (b[i][j] > 16360) // 1022.5*16
        b[i][j] = 16360;
      else if (b[i][j] < 0)
        b[i][j] = 0;
    }
  }
  printf("colormatrix ---OK!\n");

} // colormatrix() end

void apply_gamma(double gamma_cur) {
  int i, j;
  for (i = 0; i < 2460; i++) {
    for (j = 0; j < 3360; j++) {
      if (pow(r[i][j], gamma_cur) > 16383) {
        r[i][j] = 16383;
      } else {
        r[i][j] = pow(r[i][j], gamma_cur);
      }

      if (pow(g[i][j], gamma_cur) > 16383) {
        g[i][j] = 16383;
      } else {
        g[i][j] = pow(g[i][j], gamma_cur);
      }

      if (pow(b[i][j], gamma_cur) > 16383) {
        b[i][j] = 16383;
      } else {
        b[i][j] = pow(b[i][j], gamma_cur);
      }

      r[i][j] = gamma_value[r[i][j] / 16]; // 16=2^14/2^10
      g[i][j] = gamma_value[g[i][j] / 16];
      b[i][j] = gamma_value[b[i][j] / 16];
    }
  }

  printf("gamma ---OK!\n");

} // apply_gamma() end

void edge_enhance() {
  for (int i = 0; i < 2460; i++) {
    for (int j = 5; j < 3293; j++) {
      if ((i == 0) || (i == 2459) || (j == 5) || (j == 3292)) {
        pic[i][j] = r[i][j];
      } else {
        pic[i][j] = 2 * r[i][j] -
                    (r[i - 1][j] + r[i][j - 1] + r[i][j + 1] + r[i + 1][j]) / 4;
      }
    }
  }

  for (int i = 1; i < 2459; i++) {
    for (int j = 6; j < 3292; j++) {
      if (pic[i][j] > 255) {
        pic[i][j] = 255;
      } else if (pic[i][j] < 0) {
        pic[i][j] = 0;
      }

      r[i][j] = pic[i][j];
    }
  }

  for (int i = 0; i < 2460; i++) {
    for (int j = 5; j < 3293; j++) {
      if ((i == 0) || (i == 2459) || (j == 5) || (j == 3292)) {
        pic[i][j] = g[i][j];
      } else {
        pic[i][j] = 2 * g[i][j] -
                    (g[i - 1][j] + g[i][j - 1] + g[i][j + 1] + g[i + 1][j]) / 4;
      }
    }
  }

  for (int i = 1; i < 2459; i++) {
    for (int j = 6; j < 3292; j++) {
      if (pic[i][j] > 255) {
        pic[i][j] = 255;
      } else if (pic[i][j] < 0) {
        pic[i][j] = 0;
      }

      g[i][j] = pic[i][j];
    }
  }

  for (int i = 0; i < 2460; i++) {
    for (int j = 5; j < 3293; j++) {
      if ((i == 0) || (i == 2459) || (j == 5) || (j == 3292)) {
        pic[i][j] = b[i][j];
      } else {
        pic[i][j] = 2 * b[i][j] -
                    (b[i - 1][j] + b[i][j - 1] + b[i][j + 1] + b[i + 1][j]) / 4;
      }
    }
  }

  for (int i = 1; i < 2459; i++) {
    for (int j = 6; j < 3292; j++) {
      if (pic[i][j] > 255) {
        pic[i][j] = 255;
      } else if (pic[i][j] < 0) {
        pic[i][j] = 0;
      }

      b[i][j] = pic[i][j];
    }
  }

  printf("Edge Enhancement ---OK!\n");

} // edge_enhance() end

void rgb2bmp(char *name) {
  FILE *fp;
  int zero = 0;
  char filename[100];
  strcpy(filename, name);
  strcat(filename, ".bmp");

  fp = fopen(filename, "wb");

  // HEADER
  fwrite(&"BM", 1, 2, fp);
  int temp = 54 + 3264 * 2448 * 3;
  fwrite(&temp, 4, 1, fp);

  for (int bmpcount = 0; bmpcount < 4; bmpcount++) // reserved
    fwrite(&zero, 1, 1, fp);

  temp = 54; // start address of datas
  fwrite(&temp, 4, 1, fp);

  temp = 40; // size of [0E-35]
  fwrite(&temp, 1, 1, fp);
  fwrite(&zero, 1, 1, fp);
  fwrite(&zero, 1, 1, fp);
  fwrite(&zero, 1, 1, fp);

  int width = 3264; // width
  fwrite(&width, 4, 1, fp);

  int height = 2448; // height
  fwrite(&height, 4, 1, fp);

  temp = 1; // color planes
  fwrite(&temp, 2, 1, fp);

  temp = 24; // bits/pixel
  fwrite(&temp, 2, 1, fp);
  fwrite(&zero, 4, 1, fp); // Compression format

  temp = 23970816; // data size->3264*2448*3
  fwrite(&temp, 4, 1, fp);
  fwrite(&zero, 4, 1, fp); // width pixels/meter
  fwrite(&zero, 4, 1, fp); // height pixels/meter
  fwrite(&zero, 2, 1, fp); // color tables
  fwrite(&zero, 1, 1, fp);
  fwrite(&zero, 1, 1, fp);
  fwrite(&zero, 4, 1, fp); // color indexes

  // body

  for (int bmpi = 2453; bmpi >= 6; bmpi--) {
    for (int bmpj = 17; bmpj < 3281; bmpj++) {
      fwrite(&(b[bmpi][bmpj]), 1, 1, fp);
      fwrite(&(g[bmpi][bmpj]), 1, 1, fp);
      fwrite(&(r[bmpi][bmpj]), 1, 1, fp);
    }
  }

  fclose(fp);
  printf("%s output finished.\n", name);

} // rgb2bmp() end

void CFA_to_krkb() {
  printf("start get color interpolation\n");

  // init color => bilinear
  // get three channel RGB

  int i, j;
  // bilinear
  for (i = 1; i < 2459; i = i + 2) // Height
  {
    for (j = 2; j < 3359; j = j + 2) // Width
    {
      r[i][j] = (pic[i - 1][j - 1] + pic[i - 1][j + 1] + pic[i + 1][j - 1] +
                 pic[i + 1][j + 1]) /
                4;
      g[i][j] =
          (pic[i - 1][j] + pic[i][j - 1] + pic[i][j + 1] + pic[i + 1][j]) / 4;
      b[i][j] = pic[i][j];
    }
  }

  for (i = 1; i < 2459; i = i + 2) // Height
  {
    for (j = 1; j < 3359; j = j + 2) // Width
    {
      r[i][j] = (pic[i - 1][j] + pic[i + 1][j]) / 2;
      g[i][j] = pic[i][j];
      b[i][j] = (pic[i][j - 1] + pic[i][j + 1]) / 2;
    }
  }

  for (i = 2; i < 2459; i = i + 2) // Height
  {
    for (j = 2; j < 3359; j = j + 2) // Width
    {
      r[i][j] = (pic[i][j - 1] + pic[i][j + 1]) / 2;
      g[i][j] = pic[i][j];
      b[i][j] = (pic[i - 1][j] + pic[i + 1][j]) / 2;
    }
  }

  for (i = 2; i < 2459; i = i + 2) // Height
  {
    for (j = 1; j < 3359; j = j + 2) // Width
    {
      r[i][j] = pic[i][j];
      g[i][j] =
          (pic[i - 1][j] + pic[i][j - 1] + pic[i][j + 1] + pic[i + 1][j]) / 4;
      b[i][j] = (pic[i - 1][j - 1] + pic[i - 1][j + 1] + pic[i + 1][j - 1] +
                 pic[i + 1][j + 1]) /
                4;
    }
  }

  for (i = 1; i < 2459; i = i + 1) {   // height
    for (j = 1; j < 3359; j = j + 1) { // width
      r[i][j] = 0;
      g[i][j] = 0;
      b[i][j] = 0;
      if (i % 2 == 0) {   // even height
        if (j % 2 == 0) { // g=> r b

          r[i][j] = (pic[i][j - 1] + pic[i][j + 1]) / 2;
          g[i][j] = pic[i][j];
          b[i][j] = (pic[i - 1][j] + pic[i + 1][j]) / 2;
        } else { // r=>g b

          r[i][j] = pic[i][j];
          g[i][j] =
              (pic[i][j - 1] + pic[i][j + 1] + pic[i - 1][j] + pic[i + 1][j]) /
              4;
          b[i][j] = (pic[i - 1][j - 1] + pic[i + 1][j - 1] + pic[i - 1][j + 1] +
                     pic[i + 1][j + 1]) /
                    4;
        }
      } else {            // odd height
        if (j % 2 == 0) { // b=>g r

          r[i][j] = (pic[i - 1][j - 1] + pic[i + 1][j - 1] + pic[i - 1][j + 1] +
                     pic[i + 1][j + 1]) /
                    4;
          g[i][j] =
              (pic[i][j - 1] + pic[i][j + 1] + pic[i - 1][j] + pic[i + 1][j]) /
              4;
          b[i][j] = pic[i][j];
        } else { // g=> r b

          r[i][j] = (pic[i - 1][j] + pic[i + 1][j]) / 2;
          g[i][j] = pic[i][j];
          b[i][j] = (pic[i][j - 1] + pic[i][j + 1]) / 2;
        }
      }
    }
  }

  // change to kr kb
  // kr = g - r
  // kb = g - b

  for (i = 0; i < 2459; i = i + 1) {   // width
    for (j = 0; j < 3359; j = j + 1) { // height
      kr[i][j] = g[i][j] - r[i][j];
      kb[i][j] = g[i][j] - b[i][j];
    }
  }
  // horizontal and vertical edge-filters
  // threshold don't know so use another method
  int dhg, dvg, dhb, dvb, dhr, dvr;
  for (i = 3; i < 2458; i = i + 1) {   // height
    for (j = 3; j < 3358; j = j + 1) { // width
      if (i % 2 == 0) {                // even height
        if (j % 2 == 0) {              // g=> kr kb
          kb[i][j] = (kb[i + 1][j] + kb[i - 1][j]) / 2;
          kr[i][j] = (kr[i][j + 1] + kr[i][j - 1]) / 2;
        } else { // r=>g kb
          dhg = abs(g[i][j + 1] - g[i][j - 1]) +
                abs(2 * (kr[i][j]) - kr[i][j - 2] - kr[i][j + 2]);
          dvg = abs(g[i + 1][j] - g[i - 1][j]) +
                abs(2 * (kr[i][j]) - kr[i - 2][j] - kr[i + 2][j]);
          dhb = abs(kb[i - 1][j + 1] - kb[i + 1][j - 1]) +
                abs(2 * (kr[i][j]) - kr[i + 2][j - 2] - kr[i - 2][j + 2]);
          dvb = abs(kb[i + 1][j + 1] - kb[i - 1][j - 1]) +
                abs(2 * (kr[i][j]) - kr[i - 2][j - 2] - kr[i + 2][j + 2]);
          if (dhg < dvg) {
            g[i][j] = (g[i][j + 1] + g[i][j - 1]) / 2;
          } else {
            if (dhg > dvg) {
              g[i][j] = (g[i + 1][j] + g[i - 1][j]) / 2;
            } else { //==
              g[i][j] =
                  (g[i][j + 1] + g[i][j - 1] + g[i][j + 1] + g[i][j - 1]) / 4;
            }
          }

          if (dhb < dvb) {
            kb[i][j] = (kb[i - 1][j + 1] + kb[i + 1][j - 1]) / 2;
          } else {
            if (dhb > dvb) {
              kb[i][j] = (kb[i + 1][j + 1] + kb[i - 1][j - 1]) / 2;
            } else { //==
              kb[i][j] = (kb[i - 1][j + 1] + kb[i + 1][j - 1] +
                          kb[i + 1][j + 1] + kb[i - 1][j - 1]) /
                         4;
            }
          }
        }
      } else {            // odd height
        if (j % 2 == 0) { // b=>g kr
          dhg = abs(g[i][j + 1] - g[i][j - 1]) +
                abs(2 * (kb[i][j]) - kb[i][j - 2] - kb[i][j + 2]);
          dvg = abs(g[i + 1][j] - g[i - 1][j]) +
                abs(2 * (kb[i][j]) - kb[i - 2][j] - kb[i + 2][j]);
          dhr = abs(kr[i - 1][j + 1] - kr[i + 1][j - 1]) +
                abs(2 * (kb[i][j]) - kb[i + 2][j - 2] - kb[i - 2][j + 2]);
          dvr = abs(kr[i + 1][j + 1] - kr[i - 1][j - 1]) +
                abs(2 * (kb[i][j]) - kb[i - 2][j - 2] - kb[i + 2][j + 2]);
          if (dhg < dvg) {
            g[i][j] = (g[i][j + 1] + g[i][j - 1]) / 2;
          } else {
            if (dhg > dvg) {
              g[i][j] = (g[i + 1][j] + g[i - 1][j]) / 2;
            } else { //==
              g[i][j] =
                  (g[i][j + 1] + g[i][j - 1] + g[i][j + 1] + g[i][j - 1]) / 4;
            }
          }

          if (dhr < dvr) {
            kr[i][j] = (kr[i - 1][j + 1] + kr[i + 1][j - 1]) / 2;
          } else {
            if (dhr > dvr) {
              kr[i][j] = (kr[i + 1][j + 1] + kr[i - 1][j - 1]) / 2;
            } else { //==
              kr[i][j] = (kr[i - 1][j + 1] + kr[i + 1][j - 1] +
                          kr[i + 1][j + 1] + kr[i - 1][j - 1]) /
                         4;
            }
          }
        } else { // g=> kr kb
          kb[i][j] = (kb[i][j + 1] + kb[i][j - 1]) / 2;
          kr[i][j] = (kr[i + 1][j] + kr[i - 1][j]) / 2;
        }
      }
    }
  }
  // back rgb domain
  for (i = 0; i < 2459; i = i + 1) {   // width
    for (j = 0; j < 3359; j = j + 1) { // height
      r[i][j] = g[i][j] - kr[i][j];
      b[i][j] = g[i][j] - kb[i][j];
    }
  }
}

void getColor() {
  int i, j;
  // 1, 1
  for (i = 1; i < 2459; i = i + 2) // Height
  {
    for (j = 2; j < 3359; j = j + 2) // Width
    {
      r[i][j] = (pic[i - 1][j - 1] + pic[i - 1][j + 1] + pic[i + 1][j - 1] +
                 pic[i + 1][j + 1]) /
                4;
      g[i][j] =
          (pic[i - 1][j] + pic[i][j - 1] + pic[i][j + 1] + pic[i + 1][j]) / 4;
      b[i][j] = pic[i][j];
    }
  }
  printf("1, 1 ---OK!\n");

  // 2, 1
  for (i = 1; i < 2459; i = i + 2) // Height
  {
    for (j = 1; j < 3359; j = j + 2) // Width
    {
      r[i][j] = (pic[i - 1][j] + pic[i + 1][j]) / 2;
      g[i][j] = pic[i][j];
      b[i][j] = (pic[i][j - 1] + pic[i][j + 1]) / 2;
    }
  }
  printf("2, 1 ---OK!\n");

  // 1, 2
  for (i = 2; i < 2459; i = i + 2) // Height
  {
    for (j = 2; j < 3359; j = j + 2) // Width
    {
      r[i][j] = (pic[i][j - 1] + pic[i][j + 1]) / 2;
      g[i][j] = pic[i][j];
      b[i][j] = (pic[i - 1][j] + pic[i + 1][j]) / 2;
    }
  }
  printf("1, 2 ---OK!\n");

  // 2, 2
  for (i = 2; i < 2459; i = i + 2) // Height
  {
    for (j = 1; j < 3359; j = j + 2) // Width
    {
      r[i][j] = pic[i][j];
      g[i][j] =
          (pic[i - 1][j] + pic[i][j - 1] + pic[i][j + 1] + pic[i + 1][j]) / 4;
      b[i][j] = (pic[i - 1][j - 1] + pic[i - 1][j + 1] + pic[i + 1][j - 1] +
                 pic[i + 1][j + 1]) /
                4;
    }
  }
  printf("2, 2 ---OK!\n");
}

void NearestNeighborInterpolation() {
  int i, j;
  for (i = 1; i < 2459; i = i + 2) // Height
  {
    for (j = 2; j < 3359; j = j + 2) // Width
    {
      r[i][j] = pic[i - 1][j + 1];
      g[i][j] = pic[i][j - 1];
      b[i][j] = pic[i][j];
    }
  }
  printf("1, 1 ---OK!\n");

  for (i = 1; i < 2459; i = i + 2) // Height
  {
    for (j = 1; j < 3359; j = j + 2) // Width
    {
      r[i][j] = pic[i - 1][j];
      g[i][j] = pic[i][j];
      b[i][j] = pic[i][j - 1];
    }
  }
  printf("2, 1 ---OK!\n");

  for (i = 2; i < 2459; i = i + 2) // Height
  {
    for (j = 2; j < 3359; j = j + 2) // Width
    {
      r[i][j] = pic[i][j + 1];
      g[i][j] = pic[i][j];
      b[i][j] = pic[i + 1][j];
    }
  }
  printf("1, 2 ---OK!\n");

  for (i = 2; i < 2459; i = i + 2) // Height
  {
    for (j = 1; j < 3359; j = j + 2) // Width
    {
      r[i][j] = pic[i][j];
      g[i][j] = pic[i][j - 1];
      b[i][j] = pic[i + 1][j - 1];
    }
  }
  printf("2, 2 ---OK!\n");
}

void changeSaturation(double change) {
  for (int i = 1; i < 2459; i++) // Height
  {
    for (int j = 1; j < 3359; j++) // Width
    {
      double P = sqrt(r[i][j] * r[i][j] * 0.299     // 0.299
                      + g[i][j] * g[i][j] * 0.587   // 0.587
                      + b[i][j] * b[i][j] * 0.114); // 0.114

      r[i][j] * 0.299        // 0.299
          + g[i][j] * 0.587  // 0.587
          + b[i][j] * 0.114; // 0.114

      r[i][j] = P + (r[i][j] - P) * change;
      g[i][j] = P + (g[i][j] - P) * change;
      b[i][j] = P + (b[i][j] - P) * change;
    }
  }
}

void RGB_2HSV_2RGB(double change) {
  int themin, themax, delta;
  int h;
  double s, v;

  for (int i = 1; i < 2459; i++) // Height
  {

    for (int j = 1; j < 3359; j++) // Width
    {
      //
      //   Calculate HSV from RGB
      //   Hue is in degrees
      //   Lightness is betweeen 0 and 1
      //   Saturation is between 0 and 1
      //

      int r1 = r[i][j];
      int g1 = g[i][j];
      int b1 = b[i][j];

      themin = (r1 < g1 ? r1 : g1) < b1 ? (r1 < g1 ? r1 : g1) : b1;
      themax = (r1 > g1 ? r1 : g1) > b1 ? (r1 > g1 ? r1 : g1) : b1;
      delta = themax - themin;

      h = 0;
      s = 0;
      v = themax;

      if (themax > 0)
        s = delta / themax;

      if (delta > 0) {
        if (themax == r[i][j] && themax != g[i][j])
          h += (g[i][j] - b[i][j]) / delta;

        if (themax == g[i][j] && themax != b[i][j])
          h += (2 + (b[i][j] - r[i][j]) / delta);

        if (themax == b[i][j] && themax != r[i][j])
          h += (4 + (r[i][j] - g[i][j]) / delta);

        h *= 60;
      }

      s = s * change;

      //   Calculate RGB from HSV, reverse of RGB2HSV()
      //   Hue is in degrees
      //   Lightness is between 0 and 1
      //   Saturation is between 0 and 1

      float m, n, f;

      if ((j > 100) && (j < 120)) {
        // printf("%d ", h);
      }

      int i;

      i = floor(h);
      i = h % 6;

      if ((j > 100) && (j < 120)) {
        // printf("=>  %d \n", h);
      }

      f = h - i;

      if (!(i & 1))
        f = 1 - f; // if i is even

      m = v * (1 - s);
      n = v * (1 - s * f);

      switch (i) {
      case 6:
      case 0:
        r[i][j] = v;
        g[i][j] = n;
        b[i][j] = m;
        break;
      case 1:
        r[i][j] = n;
        g[i][j] = v;
        b[i][j] = m;
        break;
      case 2:
        r[i][j] = m;
        g[i][j] = v;
        b[i][j] = n;
        break;
      case 3:
        r[i][j] = m;
        g[i][j] = n;
        b[i][j] = v;
        break;
      case 4:
        r[i][j] = n;
        g[i][j] = m;
        b[i][j] = v;
        break;
      case 5:
        r[i][j] = v;
        g[i][j] = m;
        b[i][j] = n;
        break;
      }

      while (h < 0)
        h += 360;

      while (h > 360)
        h -= 360;

      if (h < 120) {
        r[i][j] = (120 - h) / 60;
        g[i][j] = h / 60;
        b[i][j] = 0;
      } else if (h < 240) {
        r[i][j] = 0;
        g[i][j] = (240 - h) / 60;
        b[i][j] = (h - 120) / 60;
      } else {
        r[i][j] = (h - 240) / 60;
        g[i][j] = 0;
        b[i][j] = (360 - h) / 60;
      }

      r[i][j] = (r[i][j] < 1) ? r[i][j] : 1;
      g[i][j] = (g[i][j] < 1) ? g[i][j] : 1;
      b[i][j] = (b[i][j] < 1) ? b[i][j] : 1;

      r[i][j] = (1 - s + s * r[i][j]) * v;
      g[i][j] = (1 - s + s * g[i][j]) * v;
      b[i][j] = (1 - s + s * b[i][j]) * v;
    }
  }
}

void RGBtoHSV(int r, int g, int b, double &h, double &s, double &v) {
  int themin, themax, delta;
  themin = std::min({r, g, b});
  themax = std::max({r, g, b});
  delta = themax - themin;

  v = themax / 255.0;
  s = (themax == 0) ? 0 : (double)delta / themax;

  if (s == 0)
    h = 0;
  else {
    if (r == themax)
      h = (g - b) / (double)delta;
    else if (g == themax)
      h = 2 + (b - r) / (double)delta;
    else
      h = 4 + (r - g) / (double)delta;

    h *= 60;
    if (h < 0)
      h += 360;
  }
}

void HSVtoRGB(double h, double s, double v, int &r, int &g, int &b) {
  int i;
  double f, p, q, t;

  if (s == 0) {
    r = g = b = round(v * 255);
    return;
  }

  h /= 60;
  i = floor(h);
  f = h - i;
  p = v * (1 - s);
  q = v * (1 - s * f);
  t = v * (1 - s * (1 - f));

  switch (i) {
  case 0:
    r = round(v * 255);
    g = round(t * 255);
    b = round(p * 255);
    break;
  case 1:
    r = round(q * 255);
    g = round(v * 255);
    b = round(p * 255);
    break;
  case 2:
    r = round(p * 255);
    g = round(v * 255);
    b = round(t * 255);
    break;
  case 3:
    r = round(p * 255);
    g = round(q * 255);
    b = round(v * 255);
    break;
  case 4:
    r = round(t * 255);
    g = round(p * 255);
    b = round(v * 255);
    break;
  default:
    r = round(v * 255);
    g = round(p * 255);
    b = round(q * 255);
    break;
  }
}

void changeVibrance(double vibrance_factor, double saturation_threshold = 0.2) {
  vibrance_factor *= -1; // Invert vibrance factor

  for (int i = 1; i < 2459; i++) {   // Height
    for (int j = 1; j < 3359; j++) { // Width
      double h, s, v;
      RGBtoHSV(r[i][j], g[i][j], b[i][j], h, s, v);

      // Selective saturation adjustment
      if (s > saturation_threshold) {
        s = std::min(1.0, s + (vibrance_factor * s));
      }

      HSVtoRGB(h, s, v, r[i][j], g[i][j], b[i][j]);
    }
  }
}

void adaptiveVibrance(double vibrance_factor, int window_size = 3) {
  vibrance_factor *= -1; // Invert vibrance factor

  for (int i = 1; i < 2459; i++) {   // Height
    for (int j = 1; j < 3359; j++) { // Width
      double h, s, v;
      RGBtoHSV(r[i][j], g[i][j], b[i][j], h, s, v);

      // Calculate local average saturation
      double local_avg_saturation = 0.0;
      int count = 0;
      for (int k = -window_size / 2; k <= window_size / 2; k++) {
        for (int l = -window_size / 2; l <= window_size / 2; l++) {
          if (i + k >= 0 && i + k < 2459 && j + l >= 0 && j + l < 3359) {
            double h_neighbor, s_neighbor, v_neighbor;
            RGBtoHSV(r[i + k][j + l], g[i + k][j + l], b[i + k][j + l],
                     h_neighbor, s_neighbor, v_neighbor);
            local_avg_saturation += s_neighbor;
            count++;
          }
        }
      }
      local_avg_saturation /= count;

      // Adjust saturation based on local average
      if (s > local_avg_saturation) {
        s = std::max(0.0, s + (vibrance_factor * (s - local_avg_saturation)));
      } else {
        s = std::min(1.0, s + (vibrance_factor * (s - local_avg_saturation)));
      }

      HSVtoRGB(h, s, v, r[i][j], g[i][j], b[i][j]);
    }
  }
}

void sharpen(double factor) {
  int width = 3360;
  int height = 2460;
  int kernel[3][3] = {{0, -1, 0}, {-1, 5, -1}, {0, -1, 0}};

  int **temp_r = new int *[height];
  int **temp_g = new int *[height];
  int **temp_b = new int *[height];
  for (int i = 0; i < height; ++i) {
    temp_r[i] = new int[width];
    temp_g[i] = new int[width];
    temp_b[i] = new int[width];
  }
  for (int y = 1; y < height - 1; ++y) {
    for (int x = 1; x < width - 1; ++x) {
      int sumR = 0, sumG = 0, sumB = 0;
      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          int pixelR = r[y + ky][x + kx];
          int pixelG = g[y + ky][x + kx];
          int pixelB = b[y + ky][x + kx];
          sumR += pixelR * kernel[ky + 1][kx + 1];
          sumG += pixelG * kernel[ky + 1][kx + 1];
          sumB += pixelB * kernel[ky + 1][kx + 1];
        }
      }
      temp_r[y][x] =
          std::min(std::max(int(r[y][x] + factor * (sumR - r[y][x])), 0), 255);
      temp_g[y][x] =
          std::min(std::max(int(g[y][x] + factor * (sumG - g[y][x])), 0), 255);
      temp_b[y][x] =
          std::min(std::max(int(b[y][x] + factor * (sumB - b[y][x])), 0), 255);
    }
  }
  for (int y = 1; y < height - 1; ++y) {
    for (int x = 1; x < width - 1; ++x) {
      r[y][x] = temp_r[y][x];
      g[y][x] = temp_g[y][x];
      b[y][x] = temp_b[y][x];
    }
  }

  for (int i = 0; i < height; ++i) {
    delete[] temp_r[i];
    delete[] temp_g[i];
    delete[] temp_b[i];
  }
  delete[] temp_r;
  delete[] temp_g;
  delete[] temp_b;

  printf("sharpen ---OK!\n");
}

int main(int argc, char *argv[]) {
  for (int k = 1; k < argc; k++) {
    FILE *fp;
    fp = fopen(argv[k], "rb"); // Open file
    fseek(fp, 862144, SEEK_SET);

    for (int i = 0; i < 2460; i++) {
      for (int j = 0; j < 3360; j++) {
        pic[i][j] = 0;
        fread(pic[i] + j, 1, 1, fp);
        pic[i][j] <<= 8;
        fread(pic[i] + j, 1, 1, fp);
      }
    }
    fclose(fp); // Close file
    OB();
    WB();

    getColor();
    // color_interpolation();
    NearestNeighborInterpolation();

    // colormatrix();
    apply_gamma(1);
    edge_enhance();

    changeSaturation(1.1);
    RGB_2HSV_2RGB(1.3);
    adjustContrast(1.1);
    changeVibrance(-0.7, 0.2);
    // adaptiveVibrance(-1, 3);
    // sharpen(1.2);

    rgb2bmp(argv[k]);
  }

  return 0;
}