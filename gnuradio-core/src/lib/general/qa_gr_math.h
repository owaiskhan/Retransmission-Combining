/* -*- c++ -*- */
/*
 * Copyright 2008 Free Software Foundation, Inc.
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
#ifndef _QA_GR_MATH_H_
#define _QA_GR_MATH_H_

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

class qa_gr_math : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE(qa_gr_math);
  CPPUNIT_TEST(test_binary_slicer1);
  CPPUNIT_TEST(test_quad_0deg_slicer1);
  CPPUNIT_TEST(test_quad_45deg_slicer1);
  CPPUNIT_TEST_SUITE_END();

 private:
  void test_binary_slicer1();
  void test_quad_0deg_slicer1();
  void test_quad_45deg_slicer1();
};

#endif /* _QA_GR_MATH_H_ */
