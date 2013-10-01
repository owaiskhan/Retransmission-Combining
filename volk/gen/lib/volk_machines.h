
// This file is automatically generated by make_machines_h.py.
// Do not edit this file.

#ifndef INCLUDED_LIBVOLK_MACHINES_H
#define INCLUDED_LIBVOLK_MACHINES_H

#include <volk/volk_common.h>
#include <volk/volk_typedefs.h>

__VOLK_DECL_BEGIN

struct volk_machine {
   const unsigned int caps; //capabilities (i.e., archs compiled into this machine, in the volk_get_lvarch format)
   const char *name;
   const unsigned int alignment; //the maximum byte alignment required for functions in this library
    const char *volk_16ic_s32f_magnitude_32f_a_name;
    const char *volk_16ic_s32f_magnitude_32f_a_indices[18];
    const int volk_16ic_s32f_magnitude_32f_a_arch_defs[18];
    const p_16ic_s32f_magnitude_32f_a volk_16ic_s32f_magnitude_32f_a_archs[18];
    const int volk_16ic_s32f_magnitude_32f_a_n_archs;
    const char *volk_32f_s32f_multiply_32f_a_name;
    const char *volk_32f_s32f_multiply_32f_a_indices[18];
    const int volk_32f_s32f_multiply_32f_a_arch_defs[18];
    const p_32f_s32f_multiply_32f_a volk_32f_s32f_multiply_32f_a_archs[18];
    const int volk_32f_s32f_multiply_32f_a_n_archs;
    const char *volk_64u_byteswap_a_name;
    const char *volk_64u_byteswap_a_indices[18];
    const int volk_64u_byteswap_a_arch_defs[18];
    const p_64u_byteswap_a volk_64u_byteswap_a_archs[18];
    const int volk_64u_byteswap_a_n_archs;
    const char *volk_32fc_x2_conjugate_dot_prod_32fc_u_name;
    const char *volk_32fc_x2_conjugate_dot_prod_32fc_u_indices[18];
    const int volk_32fc_x2_conjugate_dot_prod_32fc_u_arch_defs[18];
    const p_32fc_x2_conjugate_dot_prod_32fc_u volk_32fc_x2_conjugate_dot_prod_32fc_u_archs[18];
    const int volk_32fc_x2_conjugate_dot_prod_32fc_u_n_archs;
    const char *volk_8ic_x2_s32f_multiply_conjugate_32fc_a_name;
    const char *volk_8ic_x2_s32f_multiply_conjugate_32fc_a_indices[18];
    const int volk_8ic_x2_s32f_multiply_conjugate_32fc_a_arch_defs[18];
    const p_8ic_x2_s32f_multiply_conjugate_32fc_a volk_8ic_x2_s32f_multiply_conjugate_32fc_a_archs[18];
    const int volk_8ic_x2_s32f_multiply_conjugate_32fc_a_n_archs;
    const char *volk_32f_stddev_and_mean_32f_x2_a_name;
    const char *volk_32f_stddev_and_mean_32f_x2_a_indices[18];
    const int volk_32f_stddev_and_mean_32f_x2_a_arch_defs[18];
    const p_32f_stddev_and_mean_32f_x2_a volk_32f_stddev_and_mean_32f_x2_a_archs[18];
    const int volk_32f_stddev_and_mean_32f_x2_a_n_archs;
    const char *volk_32fc_x2_dot_prod_32fc_u_name;
    const char *volk_32fc_x2_dot_prod_32fc_u_indices[18];
    const int volk_32fc_x2_dot_prod_32fc_u_arch_defs[18];
    const p_32fc_x2_dot_prod_32fc_u volk_32fc_x2_dot_prod_32fc_u_archs[18];
    const int volk_32fc_x2_dot_prod_32fc_u_n_archs;
    const char *volk_8i_s32f_convert_32f_u_name;
    const char *volk_8i_s32f_convert_32f_u_indices[18];
    const int volk_8i_s32f_convert_32f_u_arch_defs[18];
    const p_8i_s32f_convert_32f_u volk_8i_s32f_convert_32f_u_archs[18];
    const int volk_8i_s32f_convert_32f_u_n_archs;
    const char *volk_32f_x2_add_32f_u_name;
    const char *volk_32f_x2_add_32f_u_indices[18];
    const int volk_32f_x2_add_32f_u_arch_defs[18];
    const p_32f_x2_add_32f_u volk_32f_x2_add_32f_u_archs[18];
    const int volk_32f_x2_add_32f_u_n_archs;
    const char *volk_32f_s32f_convert_32i_a_name;
    const char *volk_32f_s32f_convert_32i_a_indices[18];
    const int volk_32f_s32f_convert_32i_a_arch_defs[18];
    const p_32f_s32f_convert_32i_a volk_32f_s32f_convert_32i_a_archs[18];
    const int volk_32f_s32f_convert_32i_a_n_archs;
    const char *volk_8ic_deinterleave_real_8i_a_name;
    const char *volk_8ic_deinterleave_real_8i_a_indices[18];
    const int volk_8ic_deinterleave_real_8i_a_arch_defs[18];
    const p_8ic_deinterleave_real_8i_a volk_8ic_deinterleave_real_8i_a_archs[18];
    const int volk_8ic_deinterleave_real_8i_a_n_archs;
    const char *volk_32f_x2_dot_prod_32f_u_name;
    const char *volk_32f_x2_dot_prod_32f_u_indices[18];
    const int volk_32f_x2_dot_prod_32f_u_arch_defs[18];
    const p_32f_x2_dot_prod_32f_u volk_32f_x2_dot_prod_32f_u_archs[18];
    const int volk_32f_x2_dot_prod_32f_u_n_archs;
    const char *volk_32fc_s32fc_multiply_32fc_u_name;
    const char *volk_32fc_s32fc_multiply_32fc_u_indices[18];
    const int volk_32fc_s32fc_multiply_32fc_u_arch_defs[18];
    const p_32fc_s32fc_multiply_32fc_u volk_32fc_s32fc_multiply_32fc_u_archs[18];
    const int volk_32fc_s32fc_multiply_32fc_u_n_archs;
    const char *volk_16ic_s32f_deinterleave_32f_x2_a_name;
    const char *volk_16ic_s32f_deinterleave_32f_x2_a_indices[18];
    const int volk_16ic_s32f_deinterleave_32f_x2_a_arch_defs[18];
    const p_16ic_s32f_deinterleave_32f_x2_a volk_16ic_s32f_deinterleave_32f_x2_a_archs[18];
    const int volk_16ic_s32f_deinterleave_32f_x2_a_n_archs;
    const char *volk_16i_convert_8i_a_name;
    const char *volk_16i_convert_8i_a_indices[18];
    const int volk_16i_convert_8i_a_arch_defs[18];
    const p_16i_convert_8i_a volk_16i_convert_8i_a_archs[18];
    const int volk_16i_convert_8i_a_n_archs;
    const char *volk_8i_convert_16i_a_name;
    const char *volk_8i_convert_16i_a_indices[18];
    const int volk_8i_convert_16i_a_arch_defs[18];
    const p_8i_convert_16i_a volk_8i_convert_16i_a_archs[18];
    const int volk_8i_convert_16i_a_n_archs;
    const char *volk_16ic_s32f_deinterleave_real_32f_a_name;
    const char *volk_16ic_s32f_deinterleave_real_32f_a_indices[18];
    const int volk_16ic_s32f_deinterleave_real_32f_a_arch_defs[18];
    const p_16ic_s32f_deinterleave_real_32f_a volk_16ic_s32f_deinterleave_real_32f_a_archs[18];
    const int volk_16ic_s32f_deinterleave_real_32f_a_n_archs;
    const char *volk_32u_byteswap_a_name;
    const char *volk_32u_byteswap_a_indices[18];
    const int volk_32u_byteswap_a_arch_defs[18];
    const p_32u_byteswap_a volk_32u_byteswap_a_archs[18];
    const int volk_32u_byteswap_a_n_archs;
    const char *volk_16u_byteswap_a_name;
    const char *volk_16u_byteswap_a_indices[18];
    const int volk_16u_byteswap_a_arch_defs[18];
    const p_16u_byteswap_a volk_16u_byteswap_a_archs[18];
    const int volk_16u_byteswap_a_n_archs;
    const char *volk_16ic_deinterleave_16i_x2_a_name;
    const char *volk_16ic_deinterleave_16i_x2_a_indices[18];
    const int volk_16ic_deinterleave_16i_x2_a_arch_defs[18];
    const p_16ic_deinterleave_16i_x2_a volk_16ic_deinterleave_16i_x2_a_archs[18];
    const int volk_16ic_deinterleave_16i_x2_a_n_archs;
    const char *volk_32u_popcnt_a_name;
    const char *volk_32u_popcnt_a_indices[18];
    const int volk_32u_popcnt_a_arch_defs[18];
    const p_32u_popcnt_a volk_32u_popcnt_a_archs[18];
    const int volk_32u_popcnt_a_n_archs;
    const char *volk_8ic_s32f_deinterleave_32f_x2_a_name;
    const char *volk_8ic_s32f_deinterleave_32f_x2_a_indices[18];
    const int volk_8ic_s32f_deinterleave_32f_x2_a_arch_defs[18];
    const p_8ic_s32f_deinterleave_32f_x2_a volk_8ic_s32f_deinterleave_32f_x2_a_archs[18];
    const int volk_8ic_s32f_deinterleave_32f_x2_a_n_archs;
    const char *volk_8ic_x2_multiply_conjugate_16ic_a_name;
    const char *volk_8ic_x2_multiply_conjugate_16ic_a_indices[18];
    const int volk_8ic_x2_multiply_conjugate_16ic_a_arch_defs[18];
    const p_8ic_x2_multiply_conjugate_16ic_a volk_8ic_x2_multiply_conjugate_16ic_a_archs[18];
    const int volk_8ic_x2_multiply_conjugate_16ic_a_n_archs;
    const char *volk_64f_convert_32f_a_name;
    const char *volk_64f_convert_32f_a_indices[18];
    const int volk_64f_convert_32f_a_arch_defs[18];
    const p_64f_convert_32f_a volk_64f_convert_32f_a_archs[18];
    const int volk_64f_convert_32f_a_n_archs;
    const char *volk_32fc_magnitude_32f_u_name;
    const char *volk_32fc_magnitude_32f_u_indices[18];
    const int volk_32fc_magnitude_32f_u_arch_defs[18];
    const p_32fc_magnitude_32f_u volk_32fc_magnitude_32f_u_archs[18];
    const int volk_32fc_magnitude_32f_u_n_archs;
    const char *volk_32fc_magnitude_squared_32f_a_name;
    const char *volk_32fc_magnitude_squared_32f_a_indices[18];
    const int volk_32fc_magnitude_squared_32f_a_arch_defs[18];
    const p_32fc_magnitude_squared_32f_a volk_32fc_magnitude_squared_32f_a_archs[18];
    const int volk_32fc_magnitude_squared_32f_a_n_archs;
    const char *volk_32f_s32f_convert_8i_u_name;
    const char *volk_32f_s32f_convert_8i_u_indices[18];
    const int volk_32f_s32f_convert_8i_u_arch_defs[18];
    const p_32f_s32f_convert_8i_u volk_32f_s32f_convert_8i_u_archs[18];
    const int volk_32f_s32f_convert_8i_u_n_archs;
    const char *volk_8ic_deinterleave_16i_x2_a_name;
    const char *volk_8ic_deinterleave_16i_x2_a_indices[18];
    const int volk_8ic_deinterleave_16i_x2_a_arch_defs[18];
    const p_8ic_deinterleave_16i_x2_a volk_8ic_deinterleave_16i_x2_a_archs[18];
    const int volk_8ic_deinterleave_16i_x2_a_n_archs;
    const char *volk_8i_convert_16i_u_name;
    const char *volk_8i_convert_16i_u_indices[18];
    const int volk_8i_convert_16i_u_arch_defs[18];
    const p_8i_convert_16i_u volk_8i_convert_16i_u_archs[18];
    const int volk_8i_convert_16i_u_n_archs;
    const char *volk_32f_x2_multiply_32f_a_name;
    const char *volk_32f_x2_multiply_32f_a_indices[18];
    const int volk_32f_x2_multiply_32f_a_arch_defs[18];
    const p_32f_x2_multiply_32f_a volk_32f_x2_multiply_32f_a_archs[18];
    const int volk_32f_x2_multiply_32f_a_n_archs;
    const char *volk_32f_x3_sum_of_poly_32f_a_name;
    const char *volk_32f_x3_sum_of_poly_32f_a_indices[18];
    const int volk_32f_x3_sum_of_poly_32f_a_arch_defs[18];
    const p_32f_x3_sum_of_poly_32f_a volk_32f_x3_sum_of_poly_32f_a_archs[18];
    const int volk_32f_x3_sum_of_poly_32f_a_n_archs;
    const char *volk_32f_x2_dot_prod_32f_a_name;
    const char *volk_32f_x2_dot_prod_32f_a_indices[18];
    const int volk_32f_x2_dot_prod_32f_a_arch_defs[18];
    const p_32f_x2_dot_prod_32f_a volk_32f_x2_dot_prod_32f_a_archs[18];
    const int volk_32f_x2_dot_prod_32f_a_n_archs;
    const char *volk_32f_s32f_stddev_32f_a_name;
    const char *volk_32f_s32f_stddev_32f_a_indices[18];
    const int volk_32f_s32f_stddev_32f_a_arch_defs[18];
    const p_32f_s32f_stddev_32f_a volk_32f_s32f_stddev_32f_a_archs[18];
    const int volk_32f_s32f_stddev_32f_a_n_archs;
    const char *volk_32f_sqrt_32f_a_name;
    const char *volk_32f_sqrt_32f_a_indices[18];
    const int volk_32f_sqrt_32f_a_arch_defs[18];
    const p_32f_sqrt_32f_a volk_32f_sqrt_32f_a_archs[18];
    const int volk_32f_sqrt_32f_a_n_archs;
    const char *volk_16i_permute_and_scalar_add_a_name;
    const char *volk_16i_permute_and_scalar_add_a_indices[18];
    const int volk_16i_permute_and_scalar_add_a_arch_defs[18];
    const p_16i_permute_and_scalar_add_a volk_16i_permute_and_scalar_add_a_archs[18];
    const int volk_16i_permute_and_scalar_add_a_n_archs;
    const char *volk_16ic_deinterleave_real_8i_a_name;
    const char *volk_16ic_deinterleave_real_8i_a_indices[18];
    const int volk_16ic_deinterleave_real_8i_a_arch_defs[18];
    const p_16ic_deinterleave_real_8i_a volk_16ic_deinterleave_real_8i_a_archs[18];
    const int volk_16ic_deinterleave_real_8i_a_n_archs;
    const char *volk_32f_s32f_convert_8i_a_name;
    const char *volk_32f_s32f_convert_8i_a_indices[18];
    const int volk_32f_s32f_convert_8i_a_arch_defs[18];
    const p_32f_s32f_convert_8i_a volk_32f_s32f_convert_8i_a_archs[18];
    const int volk_32f_s32f_convert_8i_a_n_archs;
    const char *volk_32f_x2_add_32f_a_name;
    const char *volk_32f_x2_add_32f_a_indices[18];
    const int volk_32f_x2_add_32f_a_arch_defs[18];
    const p_32f_x2_add_32f_a volk_32f_x2_add_32f_a_archs[18];
    const int volk_32f_x2_add_32f_a_n_archs;
    const char *volk_16i_max_star_horizontal_16i_a_name;
    const char *volk_16i_max_star_horizontal_16i_a_indices[18];
    const int volk_16i_max_star_horizontal_16i_a_arch_defs[18];
    const p_16i_max_star_horizontal_16i_a volk_16i_max_star_horizontal_16i_a_archs[18];
    const int volk_16i_max_star_horizontal_16i_a_n_archs;
    const char *volk_32f_index_max_16u_a_name;
    const char *volk_32f_index_max_16u_a_indices[18];
    const int volk_32f_index_max_16u_a_arch_defs[18];
    const p_32f_index_max_16u_a volk_32f_index_max_16u_a_archs[18];
    const int volk_32f_index_max_16u_a_n_archs;
    const char *volk_32fc_x2_square_dist_32f_a_name;
    const char *volk_32fc_x2_square_dist_32f_a_indices[18];
    const int volk_32fc_x2_square_dist_32f_a_arch_defs[18];
    const p_32fc_x2_square_dist_32f_a volk_32fc_x2_square_dist_32f_a_archs[18];
    const int volk_32fc_x2_square_dist_32f_a_n_archs;
    const char *volk_32fc_s32f_magnitude_16i_a_name;
    const char *volk_32fc_s32f_magnitude_16i_a_indices[18];
    const int volk_32fc_s32f_magnitude_16i_a_arch_defs[18];
    const p_32fc_s32f_magnitude_16i_a volk_32fc_s32f_magnitude_16i_a_archs[18];
    const int volk_32fc_s32f_magnitude_16i_a_n_archs;
    const char *volk_32fc_conjugate_32fc_a_name;
    const char *volk_32fc_conjugate_32fc_a_indices[18];
    const int volk_32fc_conjugate_32fc_a_arch_defs[18];
    const p_32fc_conjugate_32fc_a volk_32fc_conjugate_32fc_a_archs[18];
    const int volk_32fc_conjugate_32fc_a_n_archs;
    const char *volk_64f_x2_min_64f_a_name;
    const char *volk_64f_x2_min_64f_a_indices[18];
    const int volk_64f_x2_min_64f_a_arch_defs[18];
    const p_64f_x2_min_64f_a volk_64f_x2_min_64f_a_archs[18];
    const int volk_64f_x2_min_64f_a_n_archs;
    const char *volk_32fc_magnitude_squared_32f_u_name;
    const char *volk_32fc_magnitude_squared_32f_u_indices[18];
    const int volk_32fc_magnitude_squared_32f_u_arch_defs[18];
    const p_32fc_magnitude_squared_32f_u volk_32fc_magnitude_squared_32f_u_archs[18];
    const int volk_32fc_magnitude_squared_32f_u_n_archs;
    const char *volk_32fc_x2_multiply_32fc_a_name;
    const char *volk_32fc_x2_multiply_32fc_a_indices[18];
    const int volk_32fc_x2_multiply_32fc_a_arch_defs[18];
    const p_32fc_x2_multiply_32fc_a volk_32fc_x2_multiply_32fc_a_archs[18];
    const int volk_32fc_x2_multiply_32fc_a_n_archs;
    const char *volk_32i_x2_or_32i_a_name;
    const char *volk_32i_x2_or_32i_a_indices[18];
    const int volk_32i_x2_or_32i_a_arch_defs[18];
    const p_32i_x2_or_32i_a volk_32i_x2_or_32i_a_archs[18];
    const int volk_32i_x2_or_32i_a_n_archs;
    const char *volk_32fc_deinterleave_real_32f_a_name;
    const char *volk_32fc_deinterleave_real_32f_a_indices[18];
    const int volk_32fc_deinterleave_real_32f_a_arch_defs[18];
    const p_32fc_deinterleave_real_32f_a volk_32fc_deinterleave_real_32f_a_archs[18];
    const int volk_32fc_deinterleave_real_32f_a_n_archs;
    const char *volk_64f_convert_32f_u_name;
    const char *volk_64f_convert_32f_u_indices[18];
    const int volk_64f_convert_32f_u_arch_defs[18];
    const p_64f_convert_32f_u volk_64f_convert_32f_u_archs[18];
    const int volk_64f_convert_32f_u_n_archs;
    const char *volk_32f_x2_min_32f_a_name;
    const char *volk_32f_x2_min_32f_a_indices[18];
    const int volk_32f_x2_min_32f_a_arch_defs[18];
    const p_32f_x2_min_32f_a volk_32f_x2_min_32f_a_archs[18];
    const int volk_32f_x2_min_32f_a_n_archs;
    const char *volk_32f_x2_s32f_interleave_16ic_a_name;
    const char *volk_32f_x2_s32f_interleave_16ic_a_indices[18];
    const int volk_32f_x2_s32f_interleave_16ic_a_arch_defs[18];
    const p_32f_x2_s32f_interleave_16ic_a volk_32f_x2_s32f_interleave_16ic_a_archs[18];
    const int volk_32f_x2_s32f_interleave_16ic_a_n_archs;
    const char *volk_8ic_s32f_deinterleave_real_32f_a_name;
    const char *volk_8ic_s32f_deinterleave_real_32f_a_indices[18];
    const int volk_8ic_s32f_deinterleave_real_32f_a_arch_defs[18];
    const p_8ic_s32f_deinterleave_real_32f_a volk_8ic_s32f_deinterleave_real_32f_a_archs[18];
    const int volk_8ic_s32f_deinterleave_real_32f_a_n_archs;
    const char *volk_16ic_magnitude_16i_a_name;
    const char *volk_16ic_magnitude_16i_a_indices[18];
    const int volk_16ic_magnitude_16i_a_arch_defs[18];
    const p_16ic_magnitude_16i_a volk_16ic_magnitude_16i_a_archs[18];
    const int volk_16ic_magnitude_16i_a_n_archs;
    const char *volk_32fc_deinterleave_imag_32f_a_name;
    const char *volk_32fc_deinterleave_imag_32f_a_indices[18];
    const int volk_32fc_deinterleave_imag_32f_a_arch_defs[18];
    const p_32fc_deinterleave_imag_32f_a volk_32fc_deinterleave_imag_32f_a_archs[18];
    const int volk_32fc_deinterleave_imag_32f_a_n_archs;
    const char *volk_16i_convert_8i_u_name;
    const char *volk_16i_convert_8i_u_indices[18];
    const int volk_16i_convert_8i_u_arch_defs[18];
    const p_16i_convert_8i_u volk_16i_convert_8i_u_archs[18];
    const int volk_16i_convert_8i_u_n_archs;
    const char *volk_32fc_32f_multiply_32fc_a_name;
    const char *volk_32fc_32f_multiply_32fc_a_indices[18];
    const int volk_32fc_32f_multiply_32fc_a_arch_defs[18];
    const p_32fc_32f_multiply_32fc_a volk_32fc_32f_multiply_32fc_a_archs[18];
    const int volk_32fc_32f_multiply_32fc_a_n_archs;
    const char *volk_32i_s32f_convert_32f_a_name;
    const char *volk_32i_s32f_convert_32f_a_indices[18];
    const int volk_32i_s32f_convert_32f_a_arch_defs[18];
    const p_32i_s32f_convert_32f_a volk_32i_s32f_convert_32f_a_archs[18];
    const int volk_32i_s32f_convert_32f_a_n_archs;
    const char *volk_32f_s32f_multiply_32f_u_name;
    const char *volk_32f_s32f_multiply_32f_u_indices[18];
    const int volk_32f_s32f_multiply_32f_u_arch_defs[18];
    const p_32f_s32f_multiply_32f_u volk_32f_s32f_multiply_32f_u_archs[18];
    const int volk_32f_s32f_multiply_32f_u_n_archs;
    const char *volk_32i_x2_and_32i_a_name;
    const char *volk_32i_x2_and_32i_a_indices[18];
    const int volk_32i_x2_and_32i_a_arch_defs[18];
    const p_32i_x2_and_32i_a volk_32i_x2_and_32i_a_archs[18];
    const int volk_32i_x2_and_32i_a_n_archs;
    const char *volk_16i_x4_quad_max_star_16i_a_name;
    const char *volk_16i_x4_quad_max_star_16i_a_indices[18];
    const int volk_16i_x4_quad_max_star_16i_a_arch_defs[18];
    const p_16i_x4_quad_max_star_16i_a volk_16i_x4_quad_max_star_16i_a_archs[18];
    const int volk_16i_x4_quad_max_star_16i_a_n_archs;
    const char *volk_32fc_s32f_power_32fc_a_name;
    const char *volk_32fc_s32f_power_32fc_a_indices[18];
    const int volk_32fc_s32f_power_32fc_a_arch_defs[18];
    const p_32fc_s32f_power_32fc_a volk_32fc_s32f_power_32fc_a_archs[18];
    const int volk_32fc_s32f_power_32fc_a_n_archs;
    const char *volk_32fc_index_max_16u_a_name;
    const char *volk_32fc_index_max_16u_a_indices[18];
    const int volk_32fc_index_max_16u_a_arch_defs[18];
    const p_32fc_index_max_16u_a volk_32fc_index_max_16u_a_archs[18];
    const int volk_32fc_index_max_16u_a_n_archs;
    const char *volk_32fc_x2_multiply_conjugate_32fc_u_name;
    const char *volk_32fc_x2_multiply_conjugate_32fc_u_indices[18];
    const int volk_32fc_x2_multiply_conjugate_32fc_u_arch_defs[18];
    const p_32fc_x2_multiply_conjugate_32fc_u volk_32fc_x2_multiply_conjugate_32fc_u_archs[18];
    const int volk_32fc_x2_multiply_conjugate_32fc_u_n_archs;
    const char *volk_16i_s32f_convert_32f_u_name;
    const char *volk_16i_s32f_convert_32f_u_indices[18];
    const int volk_16i_s32f_convert_32f_u_arch_defs[18];
    const p_16i_s32f_convert_32f_u volk_16i_s32f_convert_32f_u_archs[18];
    const int volk_16i_s32f_convert_32f_u_n_archs;
    const char *volk_32fc_magnitude_32f_a_name;
    const char *volk_32fc_magnitude_32f_a_indices[18];
    const int volk_32fc_magnitude_32f_a_arch_defs[18];
    const p_32fc_magnitude_32f_a volk_32fc_magnitude_32f_a_archs[18];
    const int volk_32fc_magnitude_32f_a_n_archs;
    const char *volk_32fc_s32fc_multiply_32fc_a_name;
    const char *volk_32fc_s32fc_multiply_32fc_a_indices[18];
    const int volk_32fc_s32fc_multiply_32fc_a_arch_defs[18];
    const p_32fc_s32fc_multiply_32fc_a volk_32fc_s32fc_multiply_32fc_a_archs[18];
    const int volk_32fc_s32fc_multiply_32fc_a_n_archs;
    const char *volk_32f_convert_64f_u_name;
    const char *volk_32f_convert_64f_u_indices[18];
    const int volk_32f_convert_64f_u_arch_defs[18];
    const p_32f_convert_64f_u volk_32f_convert_64f_u_archs[18];
    const int volk_32f_convert_64f_u_n_archs;
    const char *volk_32fc_s32f_atan2_32f_a_name;
    const char *volk_32fc_s32f_atan2_32f_a_indices[18];
    const int volk_32fc_s32f_atan2_32f_a_arch_defs[18];
    const p_32fc_s32f_atan2_32f_a volk_32fc_s32f_atan2_32f_a_archs[18];
    const int volk_32fc_s32f_atan2_32f_a_n_archs;
    const char *volk_32i_s32f_convert_32f_u_name;
    const char *volk_32i_s32f_convert_32f_u_indices[18];
    const int volk_32i_s32f_convert_32f_u_arch_defs[18];
    const p_32i_s32f_convert_32f_u volk_32i_s32f_convert_32f_u_archs[18];
    const int volk_32i_s32f_convert_32f_u_n_archs;
    const char *volk_32f_s32f_normalize_a_name;
    const char *volk_32f_s32f_normalize_a_indices[18];
    const int volk_32f_s32f_normalize_a_arch_defs[18];
    const p_32f_s32f_normalize_a volk_32f_s32f_normalize_a_archs[18];
    const int volk_32f_s32f_normalize_a_n_archs;
    const char *volk_32fc_x2_s32f_square_dist_scalar_mult_32f_a_name;
    const char *volk_32fc_x2_s32f_square_dist_scalar_mult_32f_a_indices[18];
    const int volk_32fc_x2_s32f_square_dist_scalar_mult_32f_a_arch_defs[18];
    const p_32fc_x2_s32f_square_dist_scalar_mult_32f_a volk_32fc_x2_s32f_square_dist_scalar_mult_32f_a_archs[18];
    const int volk_32fc_x2_s32f_square_dist_scalar_mult_32f_a_n_archs;
    const char *volk_32f_x2_interleave_32fc_a_name;
    const char *volk_32f_x2_interleave_32fc_a_indices[18];
    const int volk_32f_x2_interleave_32fc_a_arch_defs[18];
    const p_32f_x2_interleave_32fc_a volk_32f_x2_interleave_32fc_a_archs[18];
    const int volk_32f_x2_interleave_32fc_a_n_archs;
    const char *volk_32fc_deinterleave_real_64f_a_name;
    const char *volk_32fc_deinterleave_real_64f_a_indices[18];
    const int volk_32fc_deinterleave_real_64f_a_arch_defs[18];
    const p_32fc_deinterleave_real_64f_a volk_32fc_deinterleave_real_64f_a_archs[18];
    const int volk_32fc_deinterleave_real_64f_a_n_archs;
    const char *volk_8i_s32f_convert_32f_a_name;
    const char *volk_8i_s32f_convert_32f_a_indices[18];
    const int volk_8i_s32f_convert_32f_a_arch_defs[18];
    const p_8i_s32f_convert_32f_a volk_8i_s32f_convert_32f_a_archs[18];
    const int volk_8i_s32f_convert_32f_a_n_archs;
    const char *volk_32fc_x2_multiply_conjugate_32fc_a_name;
    const char *volk_32fc_x2_multiply_conjugate_32fc_a_indices[18];
    const int volk_32fc_x2_multiply_conjugate_32fc_a_arch_defs[18];
    const p_32fc_x2_multiply_conjugate_32fc_a volk_32fc_x2_multiply_conjugate_32fc_a_archs[18];
    const int volk_32fc_x2_multiply_conjugate_32fc_a_n_archs;
    const char *volk_32fc_x2_multiply_32fc_u_name;
    const char *volk_32fc_x2_multiply_32fc_u_indices[18];
    const int volk_32fc_x2_multiply_32fc_u_arch_defs[18];
    const p_32fc_x2_multiply_32fc_u volk_32fc_x2_multiply_32fc_u_archs[18];
    const int volk_32fc_x2_multiply_32fc_u_n_archs;
    const char *volk_32f_s32f_calc_spectral_noise_floor_32f_a_name;
    const char *volk_32f_s32f_calc_spectral_noise_floor_32f_a_indices[18];
    const int volk_32f_s32f_calc_spectral_noise_floor_32f_a_arch_defs[18];
    const p_32f_s32f_calc_spectral_noise_floor_32f_a volk_32f_s32f_calc_spectral_noise_floor_32f_a_archs[18];
    const int volk_32f_s32f_calc_spectral_noise_floor_32f_a_n_archs;
    const char *volk_32fc_x2_dot_prod_32fc_a_name;
    const char *volk_32fc_x2_dot_prod_32fc_a_indices[18];
    const int volk_32fc_x2_dot_prod_32fc_a_arch_defs[18];
    const p_32fc_x2_dot_prod_32fc_a volk_32fc_x2_dot_prod_32fc_a_archs[18];
    const int volk_32fc_x2_dot_prod_32fc_a_n_archs;
    const char *volk_32f_accumulator_s32f_a_name;
    const char *volk_32f_accumulator_s32f_a_indices[18];
    const int volk_32f_accumulator_s32f_a_arch_defs[18];
    const p_32f_accumulator_s32f_a volk_32f_accumulator_s32f_a_archs[18];
    const int volk_32f_accumulator_s32f_a_n_archs;
    const char *volk_64f_x2_max_64f_a_name;
    const char *volk_64f_x2_max_64f_a_indices[18];
    const int volk_64f_x2_max_64f_a_arch_defs[18];
    const p_64f_x2_max_64f_a volk_64f_x2_max_64f_a_archs[18];
    const int volk_64f_x2_max_64f_a_n_archs;
    const char *volk_32f_s32f_convert_32i_u_name;
    const char *volk_32f_s32f_convert_32i_u_indices[18];
    const int volk_32f_s32f_convert_32i_u_arch_defs[18];
    const p_32f_s32f_convert_32i_u volk_32f_s32f_convert_32i_u_archs[18];
    const int volk_32f_s32f_convert_32i_u_n_archs;
    const char *volk_32fc_deinterleave_32f_x2_a_name;
    const char *volk_32fc_deinterleave_32f_x2_a_indices[18];
    const int volk_32fc_deinterleave_32f_x2_a_arch_defs[18];
    const p_32fc_deinterleave_32f_x2_a volk_32fc_deinterleave_32f_x2_a_archs[18];
    const int volk_32fc_deinterleave_32f_x2_a_n_archs;
    const char *volk_32f_s32f_power_32f_a_name;
    const char *volk_32f_s32f_power_32f_a_indices[18];
    const int volk_32f_s32f_power_32f_a_arch_defs[18];
    const p_32f_s32f_power_32f_a volk_32f_s32f_power_32f_a_archs[18];
    const int volk_32f_s32f_power_32f_a_n_archs;
    const char *volk_16i_s32f_convert_32f_a_name;
    const char *volk_16i_s32f_convert_32f_a_indices[18];
    const int volk_16i_s32f_convert_32f_a_arch_defs[18];
    const p_16i_s32f_convert_32f_a volk_16i_s32f_convert_32f_a_archs[18];
    const int volk_16i_s32f_convert_32f_a_n_archs;
    const char *volk_32f_x2_multiply_32f_u_name;
    const char *volk_32f_x2_multiply_32f_u_indices[18];
    const int volk_32f_x2_multiply_32f_u_arch_defs[18];
    const p_32f_x2_multiply_32f_u volk_32f_x2_multiply_32f_u_archs[18];
    const int volk_32f_x2_multiply_32f_u_n_archs;
    const char *volk_32fc_x2_conjugate_dot_prod_32fc_a_name;
    const char *volk_32fc_x2_conjugate_dot_prod_32fc_a_indices[18];
    const int volk_32fc_x2_conjugate_dot_prod_32fc_a_arch_defs[18];
    const p_32fc_x2_conjugate_dot_prod_32fc_a volk_32fc_x2_conjugate_dot_prod_32fc_a_archs[18];
    const int volk_32fc_x2_conjugate_dot_prod_32fc_a_n_archs;
    const char *volk_8ic_deinterleave_real_16i_a_name;
    const char *volk_8ic_deinterleave_real_16i_a_indices[18];
    const int volk_8ic_deinterleave_real_16i_a_arch_defs[18];
    const p_8ic_deinterleave_real_16i_a volk_8ic_deinterleave_real_16i_a_archs[18];
    const int volk_8ic_deinterleave_real_16i_a_n_archs;
    const char *volk_16i_branch_4_state_8_a_name;
    const char *volk_16i_branch_4_state_8_a_indices[18];
    const int volk_16i_branch_4_state_8_a_arch_defs[18];
    const p_16i_branch_4_state_8_a volk_16i_branch_4_state_8_a_archs[18];
    const int volk_16i_branch_4_state_8_a_n_archs;
    const char *volk_16i_x5_add_quad_16i_x4_a_name;
    const char *volk_16i_x5_add_quad_16i_x4_a_indices[18];
    const int volk_16i_x5_add_quad_16i_x4_a_arch_defs[18];
    const p_16i_x5_add_quad_16i_x4_a volk_16i_x5_add_quad_16i_x4_a_archs[18];
    const int volk_16i_x5_add_quad_16i_x4_a_n_archs;
    const char *volk_16ic_deinterleave_real_16i_a_name;
    const char *volk_16ic_deinterleave_real_16i_a_indices[18];
    const int volk_16ic_deinterleave_real_16i_a_arch_defs[18];
    const p_16ic_deinterleave_real_16i_a volk_16ic_deinterleave_real_16i_a_archs[18];
    const int volk_16ic_deinterleave_real_16i_a_n_archs;
    const char *volk_32f_x2_max_32f_a_name;
    const char *volk_32f_x2_max_32f_a_indices[18];
    const int volk_32f_x2_max_32f_a_arch_defs[18];
    const p_32f_x2_max_32f_a volk_32f_x2_max_32f_a_archs[18];
    const int volk_32f_x2_max_32f_a_n_archs;
    const char *volk_32f_s32f_convert_16i_a_name;
    const char *volk_32f_s32f_convert_16i_a_indices[18];
    const int volk_32f_s32f_convert_16i_a_arch_defs[18];
    const p_32f_s32f_convert_16i_a volk_32f_s32f_convert_16i_a_archs[18];
    const int volk_32f_s32f_convert_16i_a_n_archs;
    const char *volk_32fc_s32f_x2_power_spectral_density_32f_a_name;
    const char *volk_32fc_s32f_x2_power_spectral_density_32f_a_indices[18];
    const int volk_32fc_s32f_x2_power_spectral_density_32f_a_arch_defs[18];
    const p_32fc_s32f_x2_power_spectral_density_32f_a volk_32fc_s32f_x2_power_spectral_density_32f_a_archs[18];
    const int volk_32fc_s32f_x2_power_spectral_density_32f_a_n_archs;
    const char *volk_32fc_s32f_deinterleave_real_16i_a_name;
    const char *volk_32fc_s32f_deinterleave_real_16i_a_indices[18];
    const int volk_32fc_s32f_deinterleave_real_16i_a_arch_defs[18];
    const p_32fc_s32f_deinterleave_real_16i_a volk_32fc_s32f_deinterleave_real_16i_a_archs[18];
    const int volk_32fc_s32f_deinterleave_real_16i_a_n_archs;
    const char *volk_32f_x2_divide_32f_a_name;
    const char *volk_32f_x2_divide_32f_a_indices[18];
    const int volk_32f_x2_divide_32f_a_arch_defs[18];
    const p_32f_x2_divide_32f_a volk_32f_x2_divide_32f_a_archs[18];
    const int volk_32f_x2_divide_32f_a_n_archs;
    const char *volk_16i_max_star_16i_a_name;
    const char *volk_16i_max_star_16i_a_indices[18];
    const int volk_16i_max_star_16i_a_arch_defs[18];
    const p_16i_max_star_16i_a volk_16i_max_star_16i_a_archs[18];
    const int volk_16i_max_star_16i_a_n_archs;
    const char *volk_32fc_deinterleave_64f_x2_a_name;
    const char *volk_32fc_deinterleave_64f_x2_a_indices[18];
    const int volk_32fc_deinterleave_64f_x2_a_arch_defs[18];
    const p_32fc_deinterleave_64f_x2_a volk_32fc_deinterleave_64f_x2_a_archs[18];
    const int volk_32fc_deinterleave_64f_x2_a_n_archs;
    const char *volk_32fc_conjugate_32fc_u_name;
    const char *volk_32fc_conjugate_32fc_u_indices[18];
    const int volk_32fc_conjugate_32fc_u_arch_defs[18];
    const p_32fc_conjugate_32fc_u volk_32fc_conjugate_32fc_u_archs[18];
    const int volk_32fc_conjugate_32fc_u_n_archs;
    const char *volk_32fc_s32f_power_spectrum_32f_a_name;
    const char *volk_32fc_s32f_power_spectrum_32f_a_indices[18];
    const int volk_32fc_s32f_power_spectrum_32f_a_arch_defs[18];
    const p_32fc_s32f_power_spectrum_32f_a volk_32fc_s32f_power_spectrum_32f_a_archs[18];
    const int volk_32fc_s32f_power_spectrum_32f_a_n_archs;
    const char *volk_32f_s32f_convert_16i_u_name;
    const char *volk_32f_s32f_convert_16i_u_indices[18];
    const int volk_32f_s32f_convert_16i_u_arch_defs[18];
    const p_32f_s32f_convert_16i_u volk_32f_s32f_convert_16i_u_archs[18];
    const int volk_32f_s32f_convert_16i_u_n_archs;
    const char *volk_32f_x2_subtract_32f_a_name;
    const char *volk_32f_x2_subtract_32f_a_indices[18];
    const int volk_32f_x2_subtract_32f_a_arch_defs[18];
    const p_32f_x2_subtract_32f_a volk_32f_x2_subtract_32f_a_archs[18];
    const int volk_32f_x2_subtract_32f_a_n_archs;
    const char *volk_32f_convert_64f_a_name;
    const char *volk_32f_convert_64f_a_indices[18];
    const int volk_32f_convert_64f_a_arch_defs[18];
    const p_32f_convert_64f_a volk_32f_convert_64f_a_archs[18];
    const int volk_32f_convert_64f_a_n_archs;
    const char *volk_64u_popcnt_a_name;
    const char *volk_64u_popcnt_a_indices[18];
    const int volk_64u_popcnt_a_arch_defs[18];
    const p_64u_popcnt_a volk_64u_popcnt_a_archs[18];
    const int volk_64u_popcnt_a_n_archs;
    const char *volk_32f_s32f_32f_fm_detect_32f_a_name;
    const char *volk_32f_s32f_32f_fm_detect_32f_a_indices[18];
    const int volk_32f_s32f_32f_fm_detect_32f_a_arch_defs[18];
    const p_32f_s32f_32f_fm_detect_32f_a volk_32f_s32f_32f_fm_detect_32f_a_archs[18];
    const int volk_32f_s32f_32f_fm_detect_32f_a_n_archs;
};
    
#if LV_MACHINE_AVX_ONLY
extern struct volk_machine volk_machine_avx_only;
#endif
#if LV_MACHINE_SSSE3_32
extern struct volk_machine volk_machine_ssse3_32;
#endif
#if LV_MACHINE_SSE3_64
extern struct volk_machine volk_machine_sse3_64;
#endif
#if LV_MACHINE_SSE2_32
extern struct volk_machine volk_machine_sse2_32;
#endif
#if LV_MACHINE_GENERIC
extern struct volk_machine volk_machine_generic;
#endif
#if LV_MACHINE_SSE4_2_64
extern struct volk_machine volk_machine_sse4_2_64;
#endif
#if LV_MACHINE_SSE4_A_64
extern struct volk_machine volk_machine_sse4_a_64;
#endif
#if LV_MACHINE_NEON
extern struct volk_machine volk_machine_neon;
#endif
#if LV_MACHINE_AVX_64
extern struct volk_machine volk_machine_avx_64;
#endif
#if LV_MACHINE_SSE4_1_32
extern struct volk_machine volk_machine_sse4_1_32;
#endif
#if LV_MACHINE_SSE2_64
extern struct volk_machine volk_machine_sse2_64;
#endif
#if LV_MACHINE_SSE4_A_32
extern struct volk_machine volk_machine_sse4_a_32;
#endif
#if LV_MACHINE_ALTIVEC
extern struct volk_machine volk_machine_altivec;
#endif
#if LV_MACHINE_SSE4_2_32
extern struct volk_machine volk_machine_sse4_2_32;
#endif
#if LV_MACHINE_AVX_32
extern struct volk_machine volk_machine_avx_32;
#endif
#if LV_MACHINE_SSE2_ONLY
extern struct volk_machine volk_machine_sse2_only;
#endif
#if LV_MACHINE_SSE4_1_64
extern struct volk_machine volk_machine_sse4_1_64;
#endif
#if LV_MACHINE_SSE3_32
extern struct volk_machine volk_machine_sse3_32;
#endif
#if LV_MACHINE_SSSE3_64
extern struct volk_machine volk_machine_ssse3_64;
#endif


__VOLK_DECL_END

#endif //INCLUDED_LIBVOLK_MACHINES_H