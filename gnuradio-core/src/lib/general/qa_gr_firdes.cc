/* -*- c++ -*- */
/*
 * Copyright 2002 Free Software Foundation, Inc.
 * 
 * This file is part of GNU Radio
 * 
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#include <qa_gr_firdes.h>
#include <gr_firdes.h>
#include <cppunit/TestAssert.h>
#include <gr_complex.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>

#define	NELEM(x) (sizeof (x) / sizeof (x[0]))

using std::vector;

#if 0
static void
print_taps (std::ostream &s, vector<float> &v)
{
  
  for (unsigned int i = 0; i < v.size (); i++){
    printf ("tap[%2d] = %16.7e\n", i, v[i]);
  }
}
#endif

static void
check_symmetry (vector<float> &v)
{
  int	n = v.size ();
  int	m = n / 2;

  for (int i = 0; i < m; i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL (v[i], v[n - i - 1], 1e-9);
}

const static float t1_exp[53] = {
 -9.0525491e-04,
  2.0713841e-04,
  1.2388536e-03,
  2.9683491e-04,
 -1.7744775e-03,
 -1.3599906e-03,
  2.2031884e-03,
  3.2744040e-03,
 -1.8868084e-03,
 -5.9935520e-03,
  6.4301129e-18,
  8.9516686e-03,
  4.2178580e-03,
 -1.0998557e-02,
 -1.1173409e-02,
  1.0455756e-02,
  2.0686293e-02,
 -5.2032238e-03,
 -3.1896964e-02,
 -7.4998410e-03,
  4.3362070e-02,
  3.2502845e-02,
 -5.3328082e-02,
 -8.5621715e-02,
  6.0117975e-02,
  3.1128189e-01,
  4.3769023e-01,
  3.1128189e-01,
  6.0117975e-02,
 -8.5621715e-02,
 -5.3328082e-02,
  3.2502845e-02,
  4.3362070e-02,
 -7.4998410e-03,
 -3.1896964e-02,
 -5.2032238e-03,
  2.0686293e-02,
  1.0455756e-02,
 -1.1173409e-02,
 -1.0998557e-02,
  4.2178580e-03,
  8.9516686e-03,
  6.4301129e-18,
 -5.9935520e-03,
 -1.8868084e-03,
  3.2744040e-03,
  2.2031884e-03,
 -1.3599906e-03,
 -1.7744775e-03,
  2.9683491e-04,
  1.2388536e-03,
  2.0713841e-04,
 -9.0525491e-04
};

const static float t2_exp[53] = {
  9.0380036e-04,
 -2.0680559e-04,
 -1.2368630e-03,
 -2.9635796e-04,
  1.7716263e-03,
  1.3578053e-03,
 -2.1996482e-03,
 -3.2691427e-03,
  1.8837767e-03,
  5.9839217e-03,
 -6.4197810e-18,
 -8.9372853e-03,
 -4.2110807e-03,
  1.0980885e-02,
  1.1155456e-02,
 -1.0438956e-02,
 -2.0653054e-02,
  5.1948633e-03,
  3.1845711e-02,
  7.4877902e-03,
 -4.3292396e-02,
 -3.2450620e-02,
  5.3242393e-02,
  8.5484132e-02,
 -6.0021374e-02,
 -3.1078172e-01,
  5.6184036e-01,
 -3.1078172e-01,
 -6.0021374e-02,
  8.5484132e-02,
  5.3242393e-02,
 -3.2450620e-02,
 -4.3292396e-02,
  7.4877902e-03,
  3.1845711e-02,
  5.1948633e-03,
 -2.0653054e-02,
 -1.0438956e-02,
  1.1155456e-02,
  1.0980885e-02,
 -4.2110807e-03,
 -8.9372853e-03,
 -6.4197810e-18,
  5.9839217e-03,
  1.8837767e-03,
 -3.2691427e-03,
 -2.1996482e-03,
  1.3578053e-03,
  1.7716263e-03,
 -2.9635796e-04,
 -1.2368630e-03,
 -2.0680559e-04,
  9.0380036e-04
};

const static float t3_exp[107] = {
  -1.8970841e-06,
  -7.1057165e-04,
   5.4005696e-04,
   4.6233178e-04,
   2.0572044e-04,
   3.5209916e-04,
  -1.4098573e-03,
   1.1279077e-04,
  -6.2994129e-04,
   1.1450432e-03,
   1.3637283e-03,
  -6.4360141e-04,
   3.6509900e-04,
  -3.2864159e-03,
   7.0192874e-04,
   3.7524730e-04,
   2.0256115e-03,
   3.0641893e-03,
  -3.6618244e-03,
   7.5592739e-05,
  -5.5586505e-03,
   2.3849572e-03,
   4.0114378e-03,
   1.6636450e-03,
   4.7835698e-03,
  -1.0191196e-02,
  -3.8158931e-04,
  -5.5551580e-03,
   5.3901658e-03,
   1.1366769e-02,
  -3.0000482e-03,
   4.9341680e-03,
  -2.0093076e-02,
   5.5752542e-17,
   1.2093617e-03,
   8.6089745e-03,
   2.2382140e-02,
  -1.6854567e-02,
   1.6913920e-03,
  -3.1222520e-02,
   3.2711059e-03,
   2.2604836e-02,
   8.1451107e-03,
   3.7583180e-02,
  -5.2293688e-02,
  -8.0551542e-03,
  -4.0092729e-02,
   1.5582236e-02,
   9.7452506e-02,
  -1.6183170e-02,
   8.3281815e-02,
  -2.8196752e-01,
  -1.0965768e-01,
   5.2867508e-01,
  -1.0965768e-01,
  -2.8196752e-01,
   8.3281815e-02,
  -1.6183170e-02,
   9.7452506e-02,
   1.5582236e-02,
  -4.0092729e-02,
  -8.0551542e-03,
  -5.2293688e-02,
   3.7583180e-02,
   8.1451107e-03,
   2.2604836e-02,
   3.2711059e-03,
  -3.1222520e-02,
   1.6913920e-03,
  -1.6854567e-02,
   2.2382140e-02,
   8.6089745e-03,
   1.2093617e-03,
   5.5752542e-17,
  -2.0093076e-02,
   4.9341680e-03,
  -3.0000482e-03,
   1.1366769e-02,
   5.3901658e-03,
  -5.5551580e-03,
  -3.8158931e-04,
  -1.0191196e-02,
   4.7835698e-03,
   1.6636450e-03,
   4.0114378e-03,
   2.3849572e-03,
  -5.5586505e-03,
   7.5592739e-05,
  -3.6618244e-03,
   3.0641893e-03,
   2.0256115e-03,
   3.7524730e-04,
   7.0192874e-04,
  -3.2864159e-03,
   3.6509900e-04,
  -6.4360141e-04,
   1.3637283e-03,
   1.1450432e-03,
  -6.2994129e-04,
   1.1279077e-04,
  -1.4098573e-03,
   3.5209916e-04,
   2.0572044e-04,
   4.6233178e-04,
   5.4005696e-04,
  -7.1057165e-04,
  -1.8970841e-06
};


const static float t4_exp[] = { // low pass
 0.001059958362,
0.0002263929928,
-0.001277606934,
-0.0009675776237,
 0.001592264394,
  0.00243603508,
-0.001451682881,
-0.004769335967,
5.281541594e-18,
 0.007567512803,
 0.003658855334,
-0.009761494584,
 -0.01011830103,
 0.009636915289,
   0.0193619132,
-0.004935568199,
 -0.03060629964,
-0.007267376408,
  0.04236677289,
  0.03197422624,
 -0.05274848267,
  -0.0850463286,
  0.05989059806,
     0.31065014,
   0.4370569289,
     0.31065014,
  0.05989059806,
  -0.0850463286,
 -0.05274848267,
  0.03197422624,
  0.04236677289,
-0.007267376408,
 -0.03060629964,
-0.004935568199,
   0.0193619132,
 0.009636915289,
 -0.01011830103,
-0.009761494584,
 0.003658855334,
 0.007567512803,
5.281541594e-18,
-0.004769335967,
-0.001451682881,
  0.00243603508,
 0.001592264394,
-0.0009675776237,
-0.001277606934,
0.0002263929928,
 0.001059958362,
};


const static float t5_exp[] = { //high pass
-0.001062123571,
-0.0002268554381,
 0.001280216733,
 0.000969554123,
-0.001595516922,
-0.002441011136,
 0.001454648213,
 0.004779078532,
-5.292330097e-18,
-0.007582970895,
 -0.00366632943,
 0.009781434201,
  0.01013896987,
-0.009656600654,
 -0.01940146461,
 0.004945650231,
  0.03066881932,
  0.00728222169,
 -0.04245331511,
 -0.03203954175,
  0.05285623297,
  0.08522006124,
 -0.06001294032,
  -0.3112847209,
   0.5630782247,
  -0.3112847209,
 -0.06001294032,
  0.08522006124,
  0.05285623297,
 -0.03203954175,
 -0.04245331511,
  0.00728222169,
  0.03066881932,
 0.004945650231,
 -0.01940146461,
-0.009656600654,
  0.01013896987,
 0.009781434201,
 -0.00366632943,
-0.007582970895,
-5.292330097e-18,
 0.004779078532,
 0.001454648213,
-0.002441011136,
-0.001595516922,
 0.000969554123,
 0.001280216733,
-0.0002268554381,
-0.001062123571,
};

const static float t6_exp[] = { // bandpass
0.0002809273137,
-0.001047327649,
7.936541806e-05,
-0.0004270860809,
0.0007595835486,
0.0008966081077,
-0.0004236323002,
0.0002423936094,
-0.002212299034,
0.0004807534278,
0.0002620361629,
 0.001443728455,
 0.002229931997,
-0.002720607212,
5.731141573e-05,
-0.004297634587,
 0.001878833398,
 0.003217151389,
 0.001357055153,
 0.003965090029,
-0.008576190099,
-0.0003257228818,
-0.004805727862,
 0.004721920472,
  0.01007549558,
-0.002688719891,
 0.004467967432,
 -0.01837076992,
5.119658377e-17,
 0.001125075156,
 0.008071650751,
  0.02113764361,
 -0.01602453552,
 0.001618095324,
 -0.03004053794,
 0.003163811285,
   0.0219683405,
 0.007950295694,
  0.03682873398,
 -0.05142467469,
 -0.00794606097,
 -0.03965795785,
  0.01544955093,
  0.09681399167,
 -0.01610304788,
  0.08297294378,
  -0.2811714709,
  -0.1094062924,
   0.5275565982,
  -0.1094062924,
  -0.2811714709,
  0.08297294378,
 -0.01610304788,
  0.09681399167,
  0.01544955093,
 -0.03965795785,
 -0.00794606097,
 -0.05142467469,
  0.03682873398,
 0.007950295694,
   0.0219683405,
 0.003163811285,
 -0.03004053794,
 0.001618095324,
 -0.01602453552,
  0.02113764361,
 0.008071650751,
 0.001125075156,
5.119658377e-17,
 -0.01837076992,
 0.004467967432,
-0.002688719891,
  0.01007549558,
 0.004721920472,
-0.004805727862,
-0.0003257228818,
-0.008576190099,
 0.003965090029,
 0.001357055153,
 0.003217151389,
 0.001878833398,
-0.004297634587,
5.731141573e-05,
-0.002720607212,
 0.002229931997,
 0.001443728455,
0.0002620361629,
0.0004807534278,
-0.002212299034,
0.0002423936094,
-0.0004236323002,
0.0008966081077,
0.0007595835486,
-0.0004270860809,
7.936541806e-05,
-0.001047327649,
0.0002809273137,
};

void
qa_gr_firdes::t1 ()
{
  vector<float> taps =
    gr_firdes::low_pass ( 1.0,
			  8000,
			  1750,
			  500,
			  gr_firdes::WIN_HAMMING);

  // cout << "ntaps: " << taps.size () << endl;
  // print_taps (cout, taps);

  CPPUNIT_ASSERT_EQUAL (NELEM (t1_exp), taps.size ());
  for (unsigned int i = 0; i < taps.size (); i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL (t1_exp[i], taps[i], 1e-9);

  check_symmetry (taps);
}

void
qa_gr_firdes::t2 ()
{
  vector<float> taps =
    gr_firdes::high_pass ( 1.0,
			   8000,
			   1750,
			   500,
			   gr_firdes::WIN_HAMMING);

  // cout << "ntaps: " << taps.size () << endl;
  // print_taps (cout, taps);

  CPPUNIT_ASSERT_EQUAL (NELEM (t2_exp), taps.size ());

  for (unsigned int i = 0; i < taps.size (); i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL (t2_exp[i], taps[i], 1e-9);

  check_symmetry (taps);
}

void
qa_gr_firdes::t3 ()
{
  vector<float> taps =
    gr_firdes::band_pass ( 1.0,
			   20e6,
			   5.75e6 - (5.28e6/2),
			   5.75e6 + (5.28e6/2),
			   0.62e6,
			   gr_firdes::WIN_HAMMING);

  // cout << "ntaps: " << taps.size () << endl;
  // print_taps (cout, taps);

  CPPUNIT_ASSERT_EQUAL (NELEM (t3_exp), taps.size ());

  for (unsigned int i = 0; i < taps.size (); i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL (t3_exp[i], taps[i], 1e-7);

  check_symmetry (taps);
}

void
qa_gr_firdes::t4 ()
{
  vector<float> taps =
    gr_firdes::low_pass_2 ( 1.0,
			  8000,
			  1750,
			  500,
			  66,
			  gr_firdes::WIN_HAMMING);

  //  std::cout << "ntaps: " << taps.size () << std::endl;
  //  print_taps (std::cout, taps);

  CPPUNIT_ASSERT_EQUAL (NELEM (t4_exp), taps.size ());
  for (unsigned int i = 0; i < taps.size (); i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL (t4_exp[i], taps[i], 1e-9);


  check_symmetry (taps);
}

void
qa_gr_firdes::t5 ()
{
  vector<float> taps =
    gr_firdes::high_pass_2 ( 1.0,
			   8000,
			   1750,
			   500,
			   66,
			   gr_firdes::WIN_HAMMING);

  //  std::cout << "ntaps: " << taps.size () << std::endl;
  //  print_taps (std::cout, taps);

  CPPUNIT_ASSERT_EQUAL (NELEM (t5_exp), taps.size ());

  for (unsigned int i = 0; i < taps.size (); i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL (t5_exp[i], taps[i], 1e-9);

  check_symmetry (taps);
}

void
qa_gr_firdes::t6 ()
{
  vector<float> taps =
    gr_firdes::band_pass_2 ( 1.0,
			   20e6,
			   5.75e6 - (5.28e6/2),
			   5.75e6 + (5.28e6/2),
			   0.62e6,
			   66,
			   gr_firdes::WIN_HAMMING);

  //  std::cout << "ntaps: " << taps.size () << std::endl;
  //  print_taps (std::cout, taps);

  CPPUNIT_ASSERT_EQUAL (NELEM (t6_exp), taps.size ());

  for (unsigned int i = 0; i < taps.size (); i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL (t6_exp[i], taps[i], 1e-7);

  check_symmetry (taps);
}

void
qa_gr_firdes::t7 ()
{
}
