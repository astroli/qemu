/*
 * SPDX-License-Identifier: GPL-2.0-or-later
 *
 * QEMU model of the Virtual MXU Control.
 */
#ifndef HW_MISC_VIRT_MXU_H
#define HW_MISC_VIRT_MXU_H

#include "hw/sysbus.h"
#include "hw/register.h"

#define TYPE_VIRT_MXU_CTRL "virt.mxu"
OBJECT_DECLARE_SIMPLE_TYPE(VirtMXUCtrl, VIRT_MXU_CTRL)

#define MXU_REG_SIZE 0x1000

struct VirtMXUCtrl {
    /*< private >*/
    SysBusDevice parent_obj;

    /*< public >*/
    MemoryRegion mmio;
    qemu_irq irq;

    uint32_t a_addr;
    uint32_t a_h_addr;
    uint32_t b_addr;
    uint32_t b_h_addr;
    uint32_t c_addr;
    uint32_t c_h_addr;
    uint32_t m;
    uint32_t k;
    uint32_t n;
    uint32_t lda;
    uint32_t ldb;
    uint32_t ldc;
};

void virt_mxu_create(hwaddr addr);

#endif
