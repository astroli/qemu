/*
 * SPDX-License-Identifier: GPL-2.0-or-later
 *
 * QEMU model of the Virtual MXU Control.
 */


#include "qemu/osdep.h"
#include "qapi/error.h"
#include "qemu/log.h"
#include "migration/vmstate.h"
#include "hw/qdev-properties.h"
#include "hw/sysbus.h"
#include "hw/irq.h"
#include "hw/register.h"

#include "qemu/bitops.h"
#include "qapi/qmp/qerror.h"

#include "hw/misc/virt_mxu.h"

#define CMD_A_ADDR    0x00
#define CMD_A_H_ADDR  0x04
#define CMD_B_ADDR    0x08
#define CMD_B_H_ADDR  0x0C
#define CMD_C_ADDR    0x10
#define CMD_C_H_ADDR  0x14
#define CMD_M         0x18
#define CMD_N         0x1C
#define CMD_K         0x20
#define CMD_LDA       0x24
#define CMD_LDB       0x28
#define CMD_LDC       0x2C
#define CMD_START     0x30


#define HPL_rone  1.0
#define HPL_rzero 0.0

static void HPL_dscal
(
   const int                        N,
   const double                     ALPHA,
   double *                         X,
   const int                        INCX
)
{
   register double           x0, x1, x2, x3, x4, x5, x6, x7;
   register const double     alpha = ALPHA;
   const double              * StX;
   register int              i;
   int                       nu;
   const int                 incX2 = 2 * INCX, incX3 = 3 * INCX,
                             incX4 = 4 * INCX, incX5 = 5 * INCX,
                             incX6 = 6 * INCX, incX7 = 7 * INCX,
                             incX8 = 8 * INCX;

   if( ( N > 0 ) && ( alpha != HPL_rone ) )
   {
      if( alpha == HPL_rzero )
      {
         if( ( nu = ( N >> 3 ) << 3 ) != 0 )
         {
            StX = (double *)X + nu * INCX;

            do
            {
               (*X)     = HPL_rzero; X[incX4] = HPL_rzero;
               X[INCX ] = HPL_rzero; X[incX5] = HPL_rzero;
               X[incX2] = HPL_rzero; X[incX6] = HPL_rzero;
               X[incX3] = HPL_rzero; X[incX7] = HPL_rzero; X += incX8;

            } while( X != StX );
         }

         for( i = N - nu; i != 0; i-- ) { *X = HPL_rzero; X += INCX; }
      }
      else
      {
         if( ( nu = ( N >> 3 ) << 3 ) != 0 )
         {
            StX = X + nu * INCX;

            do
            {
               x0 = (*X);     x4 = X[incX4]; x1 = X[INCX ]; x5 = X[incX5];
               x2 = X[incX2]; x6 = X[incX6]; x3 = X[incX3]; x7 = X[incX7];

               x0 *= alpha;   x4 *= alpha;   x1 *= alpha;   x5 *= alpha;
               x2 *= alpha;   x6 *= alpha;   x3 *= alpha;   x7 *= alpha;

               (*X)     = x0; X[incX4] = x4; X[INCX ] = x1; X[incX5] = x5;
               X[incX2] = x2; X[incX6] = x6; X[incX3] = x3; X[incX7] = x7;

               X  += incX8;

            } while( X != StX );
         }

         for( i = N - nu; i != 0; i-- )
         { x0 = (*X); x0 *= alpha; *X = x0; X += INCX; }
      }
   }
}

static void mxu_dgemm_mkn
(
   const int                  M,
   const int                  N,
   const int                  K,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   const double               * B,
   const int                  LDB,
   const double               BETA,
   double                     * C,
   const int                  LDC
)
{
   register double            t0;
   int                        i, iail, iblj, icij, j, jal, jbj, jcj, l;

   for( j = 0, jbj = 0, jcj  = 0; j < N; j++, jbj += LDB, jcj += LDC )
   {
      HPL_dscal( M, BETA, C+jcj, 1 );
      for( l = 0, jal = 0, iblj = jbj; l < K; l++, jal += LDA, iblj += 1 )
      {
         t0 = ALPHA * B[iblj];
         for( i = 0, iail = jal, icij = jcj; i < M; i++, iail += 1, icij += 1 )
         { C[icij] += A[iail] * t0; }
      }
   }
}

static void mxu_dump_cfg(VirtMXUCtrl *s)
{
     qemu_log("    MXU A addr      0x%x\n", s->a_addr);
     qemu_log("    MXU A high addr 0x%x\n", s->a_h_addr);
     qemu_log("    MXU B addr      0x%x\n", s->b_addr);
     qemu_log("    MXU B high addr 0x%x\n", s->b_h_addr);
     qemu_log("    MXU C addr      0x%x\n", s->c_addr);
     qemu_log("    MXU C high addr 0x%x\n", s->c_h_addr);
     qemu_log("    MXU M           %d\n", s->m);
     qemu_log("    MXU N           %d\n", s->n);
     qemu_log("    MXU K           %d\n", s->k);
     qemu_log("    MXU LDA         %d\n", s->lda);
     qemu_log("    MXU LDB         %d\n", s->ldb);
     qemu_log("    MXU LDC         %d\n", s->ldc);
}

static uint64_t mxu_read_memory(void *opaque, hwaddr addr, unsigned int size)
{
    VirtMXUCtrl *s = VIRT_MXU_CTRL(opaque);

    switch (addr) {
    case CMD_A_ADDR:
        qemu_log("MXU A addr 0x%x\n", s->a_addr);
        break;
    case CMD_A_H_ADDR:
        qemu_log("MXU A high addr 0x%x\n", s->a_h_addr);
        break;
    case CMD_B_ADDR:
        qemu_log("MXU B addr 0x%x\n", s->b_addr);
        break;
    case CMD_B_H_ADDR:
        qemu_log("MXU B high addr 0x%x\n", s->b_h_addr);
        break;
    case CMD_C_ADDR:
        qemu_log("MXU C addr 0x%x\n", s->c_addr);
        break;
    case CMD_C_H_ADDR:
        qemu_log("MXU C high addr 0x%x\n", s->c_h_addr);
        break;
    case CMD_M:
        qemu_log("MXU M %d\n", s->m);
        break;
    case CMD_N:
        qemu_log("MXU N %d\n", s->n);
        break;
    case CMD_K:
        qemu_log("MXU K %d\n", s->k);
        break;
    case CMD_LDA:
        qemu_log("MXU LDA %d\n", s->lda);
        break;
    case CMD_LDB:
        qemu_log("MXU LDB %d\n", s->ldb);
        break;
    case CMD_LDC:
        qemu_log("MXU LDC %d\n", s->ldc);
        break;
    case CMD_START:
        qemu_log("MXU start %d\n", 0);
        break;
    default:
        qemu_log("MXU unknown rd cmd %ld\n", addr);
        break;
    }
   return 0;
}

static void mxu_handle_start_cmd(VirtMXUCtrl *s)
{
    int M = s->m;
    int N = s->n;
    int K = s->k;
    double ALPHA = 1.0;
    double BETA = 1.0;
    double LDA = s->lda;
    double LDB = s->ldb;
    double LDC = s->ldc;
    hwaddr a_addr = ((hwaddr)s->a_h_addr << 32) | s->a_addr;
    hwaddr b_addr = ((hwaddr)s->b_h_addr << 32) | s->b_addr;
    hwaddr c_addr = ((hwaddr)s->c_h_addr << 32) | s->c_addr;
    double *A = g_malloc(M * K * sizeof(double));
    double *B = g_malloc(K * N * sizeof(double));
    double *C = g_malloc(M * N * sizeof(double));

    if (M < 4 || N < 4 || K < 4 || !a_addr || !b_addr || !c_addr) {
        qemu_log("invalid MXU parameters\n");
        return;
    }
    if (!A || !B || !C) {
        qemu_log("MXU A, B, C allocation fail\n");
        return;
    }
    
    cpu_physical_memory_read(a_addr, A, M * K * sizeof(double));
    cpu_physical_memory_read(b_addr, B, K * N * sizeof(double));
    cpu_physical_memory_read(c_addr, C, M * N * sizeof(double));

    mxu_dgemm_mkn(M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);

    qemu_log("MXU done\n");

    g_free(A);
    g_free(B);
    g_free(C);
}

static void mxu_write_memory(void *opaque, hwaddr addr,
                             uint64_t val64, unsigned int size)
{
    VirtMXUCtrl *s = VIRT_MXU_CTRL(opaque);

    switch (addr) {
    case CMD_A_ADDR:
        qemu_log("MXU A addr 0x%lx\n", val64);
        s->a_addr = (uint32_t)val64;
        break;
    case CMD_A_H_ADDR:
        qemu_log("MXU A high addr 0x%lx\n", val64);
        s->a_h_addr = (uint32_t)val64;
        break;
    case CMD_B_ADDR:
        qemu_log("MXU B addr 0x%lx\n", val64);
        s->b_addr = (uint32_t)val64;
        break;
    case CMD_B_H_ADDR:
        qemu_log("MXU B high addr 0x%lx\n", val64);
        s->b_h_addr = (uint32_t)val64;
        break;
    case CMD_C_ADDR:
        qemu_log("MXU C addr 0x%lx\n", val64);
        s->c_addr = (uint32_t)val64;
        break;
    case CMD_C_H_ADDR:
        qemu_log("MXU C high addr 0x%lx\n", val64);
        s->c_h_addr = (uint32_t)val64;
        break;
    case CMD_M:
        qemu_log("MXU M %ld\n", val64);
        s->m = (uint32_t)val64;
        break;
    case CMD_N:
        qemu_log("MXU N %ld\n", val64);
        s->n = (uint32_t)val64;
        break;
    case CMD_K:
        qemu_log("MXU K %ld\n", val64);
        s->k = (uint32_t)val64;
        break;
    case CMD_LDA:
        qemu_log("MXU LDA %ld\n", val64);
        s->lda = (uint32_t)val64;
        break;
    case CMD_LDB:
        qemu_log("MXU LDB %ld\n", val64);
        s->ldb = (uint32_t)val64;
        break;
    case CMD_LDC:
        qemu_log("MXU LDC %ld\n", val64);
        s->ldc = (uint32_t)val64;
        break;
    case CMD_START:
        qemu_log("MXU start %ld\n", val64);
        mxu_dump_cfg(s);
	mxu_handle_start_cmd(s);
        break;
    default:
        qemu_log("MXU unknown wr cmd %ld\n", addr);
        break;
    }
}

static const MemoryRegionOps virt_mxu_ops = {
    .read = mxu_read_memory,
    .write = mxu_write_memory,
    .endianness = DEVICE_LITTLE_ENDIAN,
    .valid = {
        .min_access_size = 4,
        .max_access_size = 4,
    }
};

static void virt_mxu_init(Object *obj)
{
    VirtMXUCtrl *s = VIRT_MXU_CTRL(obj);

    memory_region_init_io(&s->mmio, obj, &virt_mxu_ops, s,
                          TYPE_VIRT_MXU_CTRL, MXU_REG_SIZE);
    sysbus_init_mmio(SYS_BUS_DEVICE(obj), &s->mmio);
    sysbus_init_irq(SYS_BUS_DEVICE(obj), &s->irq);
}

static const VMStateDescription vmstate_virt_mxu = {
    .name = TYPE_VIRT_MXU_CTRL,
    .version_id = 1,
    .minimum_version_id = 1,
    .fields = (VMStateField[]) {
        VMSTATE_UINT32(a_addr, VirtMXUCtrl),
        VMSTATE_UINT32(a_h_addr, VirtMXUCtrl),
        VMSTATE_UINT32(b_addr, VirtMXUCtrl),
        VMSTATE_UINT32(b_h_addr, VirtMXUCtrl),
        VMSTATE_UINT32(c_addr, VirtMXUCtrl),
        VMSTATE_UINT32(c_h_addr, VirtMXUCtrl),
        VMSTATE_UINT32(m, VirtMXUCtrl),
        VMSTATE_UINT32(n, VirtMXUCtrl),
        VMSTATE_UINT32(k, VirtMXUCtrl),
        VMSTATE_UINT32(lda, VirtMXUCtrl),
        VMSTATE_UINT32(ldb, VirtMXUCtrl),
        VMSTATE_UINT32(ldc, VirtMXUCtrl),
        VMSTATE_END_OF_LIST(),
    }
};

static void virt_mxu_class_init(ObjectClass *klass, void *data)
{
    DeviceClass *dc = DEVICE_CLASS(klass);

    dc->vmsd = &vmstate_virt_mxu;
}

static const TypeInfo virt_mxu_info = {
    .name              = TYPE_VIRT_MXU_CTRL,
    .parent            = TYPE_SYS_BUS_DEVICE,
    .instance_size     = sizeof(VirtMXUCtrl),
    .class_init        = virt_mxu_class_init,
    .instance_init     = virt_mxu_init,
};


static void virt_mxu_register_types(void)
{
    type_register_static(&virt_mxu_info);
}

type_init(virt_mxu_register_types)

/*
 * Create MXU device.
 */
void virt_mxu_create(hwaddr addr)
{
    DeviceState *dev = qdev_new(TYPE_VIRT_MXU_CTRL);
    sysbus_realize_and_unref(SYS_BUS_DEVICE(dev), &error_fatal);
    sysbus_mmio_map(SYS_BUS_DEVICE(dev), 0, addr);
}
