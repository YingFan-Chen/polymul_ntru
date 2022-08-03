LDSCRIPT   = ../stm32f4discovery.ld
LIBNAME    = opencm3_stm32f4
ARCH_FLAGS = -mthumb -mcpu=cortex-m4 -mfloat-abi=hard -mfpu=fpv4-sp-d16
DEFINES    = -DSTM32F4 -DSTM32F407VG
OBJS       = ../hal-opencm3.o

HEADER = params.h
SRC = poly_rq_mul.c
ASM_I =
ASM =

all: test.bin speed.bin fast_speed.bin

test.elf: $(HEADER) $(SRC) $(ASM_I) $(ASM) test.o $(OBJS) $(LDSCRIPT)
		$(LD) -o $@ test.o $(SRC) $(ASM) $(OBJS) $(LDFLAGS) -l$(LIBNAME)

speed.elf: $(HEADER) $(SRC) $(ASM_I) $(ASM) speed.o $(OBJS) $(LDSCRIPT)
		$(LD) -o $@ speed.o $(SRC) $(ASM) $(OBJS) $(LDFLAGS) -l$(LIBNAME)

fast_speed.elf: $(HEADER) $(SRC) $(ASM_I) $(ASM) fast_speed.o $(OBJS) $(LDSCRIPT)
		$(LD) -o $@ fast_speed.o $(SRC) $(ASM) $(OBJS) $(LDFLAGS) -l$(LIBNAME)


include ../make.mk