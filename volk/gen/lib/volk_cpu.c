/*this file is auto_generated by volk_register.py*/

#include <volk/volk_cpu.h>
#include <volk/volk_config_fixed.h>

struct VOLK_CPU volk_cpu;

#if defined(__i386__) || defined(__x86_64__) || defined(_M_IX86) || defined(_M_X64)
#  define VOLK_CPU_x86
#endif

#if defined(VOLK_CPU_x86)

//implement get cpuid for gcc compilers using a copy of cpuid.h
#if defined(__GNUC__)
#include <gcc_x86_cpuid.h>
#define cpuid_x86(op, r) __get_cpuid(op, (unsigned int *)r+0, (unsigned int *)r+1, (unsigned int *)r+2, (unsigned int *)r+3)

//implement get cpuid for MSVC compilers using __cpuid intrinsic
#elif defined(_MSC_VER)
#include <intrin.h>
#define cpuid_x86(op, r) __cpuid(r, op)

#else
#error "A get cpuid for volk is not available on this compiler..."
#endif

static inline unsigned int cpuid_eax(unsigned int op) {
    int regs[4];
    cpuid_x86 (op, regs);
    return regs[0];
}

static inline unsigned int cpuid_ebx(unsigned int op) {
    int regs[4];
    cpuid_x86 (op, regs);
    return regs[1];
}

static inline unsigned int cpuid_ecx(unsigned int op) {
    int regs[4];
    cpuid_x86 (op, regs);
    return regs[2];
}

static inline unsigned int cpuid_edx(unsigned int op) {
    int regs[4];
    cpuid_x86 (op, regs);
    return regs[3];
}
#endif

int i_can_has_generic () {
    return 1;
}

int i_can_has_altivec () {
#ifdef __PPC__
    return 1;
#else
    return 0;
#endif
}

#if defined(__arm__) && defined(__linux__)
#include <asm/hwcap.h>
#include <linux/auxvec.h>
#include <stdio.h>
#define LOOK_FOR_NEON
#endif

int i_can_has_neon () {
//it's linux-specific, but if you're compiling libvolk for NEON
//on Windows you have other problems

#ifdef LOOK_FOR_NEON
    FILE *auxvec_f;
    unsigned long auxvec[2];
    unsigned int found_neon = 0;
    auxvec_f = fopen("/proc/self/auxv", "rb");
    if(!auxvec_f) return 0;
    
    //so auxv is basically 32b of ID and 32b of value
    //so it goes like this
    while(!found_neon && auxvec_f) {
      fread(auxvec, sizeof(unsigned long), 2, auxvec_f);
      if((auxvec[0] == AT_HWCAP) && (auxvec[1] & HWCAP_NEON))
        found_neon = 1;
    }
    
    fclose(auxvec_f);
    return found_neon;

#else
    return 0;
#endif
}

int i_can_has_32 () {
#if defined(VOLK_CPU_x86)
    return 1;
#else
    return 0;
#endif
}
                
int i_can_has_64 () {
#if defined(VOLK_CPU_x86)
    unsigned int extended_fct_count = cpuid_eax(0x80000000);
    if (extended_fct_count < 0x80000001)
        return 1^1;
    unsigned int extended_features = cpuid_edx (0x80000001);
    return ((extended_features >> 29) & 1) == 1;
#else
    return 0;
#endif
}

int i_can_has_3dnow () {
#if defined(VOLK_CPU_x86)
    unsigned int extended_fct_count = cpuid_eax(0x80000000);
    if (extended_fct_count < 0x80000001)
        return 1^1;
    unsigned int extended_features = cpuid_edx (0x80000001);
    return ((extended_features >> 31) & 1) == 1;
#else
    return 0;
#endif
}

int i_can_has_abm () {
#if defined(VOLK_CPU_x86)
    unsigned int extended_fct_count = cpuid_eax(0x80000000);
    if (extended_fct_count < 0x80000001)
        return 1^1;
    unsigned int extended_features = cpuid_edx (0x80000001);
    return ((extended_features >> 5) & 1) == 1;
#else
    return 0;
#endif
}

int i_can_has_popcount () {
#if defined(VOLK_CPU_x86)
    unsigned int ecx = cpuid_ecx (1);
    return ((ecx >> 23) & 1) == 1;
#else
    return 0;
#endif
}

int i_can_has_mmx () {
#if defined(VOLK_CPU_x86)
    unsigned int edx = cpuid_edx (1);
    return ((edx >> 23) & 1) == 1;
#else
    return 0;
#endif
}

int i_can_has_sse () {
#if defined(VOLK_CPU_x86)
    unsigned int edx = cpuid_edx (1);
    return ((edx >> 25) & 1) == 1;
#else
    return 0;
#endif
}

int i_can_has_sse2 () {
#if defined(VOLK_CPU_x86)
    unsigned int edx = cpuid_edx (1);
    return ((edx >> 26) & 1) == 1;
#else
    return 0;
#endif
}

int i_can_has_orc () {
    return 1;
}

int i_can_has_sse3 () {
#if defined(VOLK_CPU_x86)
    unsigned int ecx = cpuid_ecx (1);
    return ((ecx >> 0) & 1) == 1;
#else
    return 0;
#endif
}

int i_can_has_ssse3 () {
#if defined(VOLK_CPU_x86)
    unsigned int ecx = cpuid_ecx (1);
    return ((ecx >> 9) & 1) == 1;
#else
    return 0;
#endif
}

int i_can_has_sse4_a () {
#if defined(VOLK_CPU_x86)
    unsigned int extended_fct_count = cpuid_eax(0x80000000);
    if (extended_fct_count < 0x80000001)
        return 1^1;
    unsigned int extended_features = cpuid_ecx (0x80000001);
    return ((extended_features >> 6) & 1) == 1;
#else
    return 0;
#endif
}

int i_can_has_sse4_1 () {
#if defined(VOLK_CPU_x86)
    unsigned int ecx = cpuid_ecx (1);
    return ((ecx >> 19) & 1) == 1;
#else
    return 0;
#endif
}

int i_can_has_sse4_2 () {
#if defined(VOLK_CPU_x86)
    unsigned int ecx = cpuid_ecx (1);
    return ((ecx >> 20) & 1) == 1;
#else
    return 0;
#endif
}

int i_can_has_avx () {
#if defined(VOLK_CPU_x86)
    unsigned int ecx = cpuid_ecx (1);
    return ((ecx >> 28) & 1) == 1;
#else
    return 0;
#endif
}

void volk_cpu_init() {
    volk_cpu.has_generic = &i_can_has_generic;
    volk_cpu.has_altivec = &i_can_has_altivec;
    volk_cpu.has_neon = &i_can_has_neon;
    volk_cpu.has_32 = &i_can_has_32;
    volk_cpu.has_64 = &i_can_has_64;
    volk_cpu.has_3dnow = &i_can_has_3dnow;
    volk_cpu.has_abm = &i_can_has_abm;
    volk_cpu.has_popcount = &i_can_has_popcount;
    volk_cpu.has_mmx = &i_can_has_mmx;
    volk_cpu.has_sse = &i_can_has_sse;
    volk_cpu.has_sse2 = &i_can_has_sse2;
    volk_cpu.has_orc = &i_can_has_orc;
    volk_cpu.has_sse3 = &i_can_has_sse3;
    volk_cpu.has_ssse3 = &i_can_has_ssse3;
    volk_cpu.has_sse4_a = &i_can_has_sse4_a;
    volk_cpu.has_sse4_1 = &i_can_has_sse4_1;
    volk_cpu.has_sse4_2 = &i_can_has_sse4_2;
    volk_cpu.has_avx = &i_can_has_avx;
}

unsigned int volk_get_lvarch() {
    unsigned int retval = 0;
    volk_cpu_init();
    retval += volk_cpu.has_generic() << LV_GENERIC;
    retval += volk_cpu.has_altivec() << LV_ALTIVEC;
    retval += volk_cpu.has_neon() << LV_NEON;
    retval += volk_cpu.has_32() << LV_32;
    retval += volk_cpu.has_64() << LV_64;
    retval += volk_cpu.has_3dnow() << LV_3DNOW;
    retval += volk_cpu.has_abm() << LV_ABM;
    retval += volk_cpu.has_popcount() << LV_POPCOUNT;
    retval += volk_cpu.has_mmx() << LV_MMX;
    retval += volk_cpu.has_sse() << LV_SSE;
    retval += volk_cpu.has_sse2() << LV_SSE2;
    retval += volk_cpu.has_orc() << LV_ORC;
    retval += volk_cpu.has_sse3() << LV_SSE3;
    retval += volk_cpu.has_ssse3() << LV_SSSE3;
    retval += volk_cpu.has_sse4_a() << LV_SSE4_A;
    retval += volk_cpu.has_sse4_1() << LV_SSE4_1;
    retval += volk_cpu.has_sse4_2() << LV_SSE4_2;
    retval += volk_cpu.has_avx() << LV_AVX;
    return retval;
}

