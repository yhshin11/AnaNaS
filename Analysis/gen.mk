
# module.mk appends to FILES and DICTFILES
FILES :=
DICTFILES :=

# including the module.mk file
-include $(patsubst %, %/module.mk,core)

# appending the values found in FILES to the variable SRC
SRC += $(patsubst %,core/%,$(FILES))

# appending the values found in DICTFILES to DICTH_modulename
DICTH_core := $(foreach i, $(patsubst %,core/%,$(DICTFILES)), $(wildcard $(i).h) $(wildcard $(i).hh))
# if dict header files exist, append to variable SRC
ifneq ($(DICTH_core),)
SRC += core/$(PKGNAME)_dict_core
endif

PROG += $(patsubst %,$(BINDIR)/%,$(PROGRAMS))

# appending the values found in FILES to the variable SRC
PROGSRC += $(patsubst %,core/%,$(PROGRAMS))

# VPATH += :core

# make sure the output directories are there
__dummy := $(shell mkdir -p $(DEPDIR)/core $(OBJDIR)/core)

# a couple of rules to copy executable files correctly
$(BINDIR)/%: core/%
	cp $^ $@

#$(BINDIR)/%: core/%.bin
#	cp $^ $@

$(BINDIR)/%: core/bin/%
	cp $^ $@

$(BINDIR)/%: ${OBJDIR}/core/%.bin
	cp $^ $@


# module.mk appends to FILES and DICTFILES
FILES :=
DICTFILES :=

# including the module.mk file
-include $(patsubst %, %/module.mk,exe)

# appending the values found in FILES to the variable SRC
SRC += $(patsubst %,exe/%,$(FILES))

# appending the values found in DICTFILES to DICTH_modulename
DICTH_exe := $(foreach i, $(patsubst %,exe/%,$(DICTFILES)), $(wildcard $(i).h) $(wildcard $(i).hh))
# if dict header files exist, append to variable SRC
ifneq ($(DICTH_exe),)
SRC += exe/$(PKGNAME)_dict_exe
endif

PROG += $(patsubst %,$(BINDIR)/%,$(PROGRAMS))

# appending the values found in FILES to the variable SRC
PROGSRC += $(patsubst %,exe/%,$(PROGRAMS))

# VPATH += :exe

# make sure the output directories are there
__dummy := $(shell mkdir -p $(DEPDIR)/exe $(OBJDIR)/exe)

# a couple of rules to copy executable files correctly
$(BINDIR)/%: exe/%
	cp $^ $@

#$(BINDIR)/%: exe/%.bin
#	cp $^ $@

$(BINDIR)/%: exe/bin/%
	cp $^ $@

$(BINDIR)/%: ${OBJDIR}/exe/%.bin
	cp $^ $@


# module.mk appends to FILES and DICTFILES
FILES :=
DICTFILES :=

# including the module.mk file
-include $(patsubst %, %/module.mk,fit)

# appending the values found in FILES to the variable SRC
SRC += $(patsubst %,fit/%,$(FILES))

# appending the values found in DICTFILES to DICTH_modulename
DICTH_fit := $(foreach i, $(patsubst %,fit/%,$(DICTFILES)), $(wildcard $(i).h) $(wildcard $(i).hh))
# if dict header files exist, append to variable SRC
ifneq ($(DICTH_fit),)
SRC += fit/$(PKGNAME)_dict_fit
endif

PROG += $(patsubst %,$(BINDIR)/%,$(PROGRAMS))

# appending the values found in FILES to the variable SRC
PROGSRC += $(patsubst %,fit/%,$(PROGRAMS))

# VPATH += :fit

# make sure the output directories are there
__dummy := $(shell mkdir -p $(DEPDIR)/fit $(OBJDIR)/fit)

# a couple of rules to copy executable files correctly
$(BINDIR)/%: fit/%
	cp $^ $@

#$(BINDIR)/%: fit/%.bin
#	cp $^ $@

$(BINDIR)/%: fit/bin/%
	cp $^ $@

$(BINDIR)/%: ${OBJDIR}/fit/%.bin
	cp $^ $@


# module.mk appends to FILES and DICTFILES
FILES :=
DICTFILES :=

# including the module.mk file
-include $(patsubst %, %/module.mk,ghmzz)

# appending the values found in FILES to the variable SRC
SRC += $(patsubst %,ghmzz/%,$(FILES))

# appending the values found in DICTFILES to DICTH_modulename
DICTH_ghmzz := $(foreach i, $(patsubst %,ghmzz/%,$(DICTFILES)), $(wildcard $(i).h) $(wildcard $(i).hh))
# if dict header files exist, append to variable SRC
ifneq ($(DICTH_ghmzz),)
SRC += ghmzz/$(PKGNAME)_dict_ghmzz
endif

PROG += $(patsubst %,$(BINDIR)/%,$(PROGRAMS))

# appending the values found in FILES to the variable SRC
PROGSRC += $(patsubst %,ghmzz/%,$(PROGRAMS))

# VPATH += :ghmzz

# make sure the output directories are there
__dummy := $(shell mkdir -p $(DEPDIR)/ghmzz $(OBJDIR)/ghmzz)

# a couple of rules to copy executable files correctly
$(BINDIR)/%: ghmzz/%
	cp $^ $@

#$(BINDIR)/%: ghmzz/%.bin
#	cp $^ $@

$(BINDIR)/%: ghmzz/bin/%
	cp $^ $@

$(BINDIR)/%: ${OBJDIR}/ghmzz/%.bin
	cp $^ $@


# module.mk appends to FILES and DICTFILES
FILES :=
DICTFILES :=

# including the module.mk file
-include $(patsubst %, %/module.mk,selectors)

# appending the values found in FILES to the variable SRC
SRC += $(patsubst %,selectors/%,$(FILES))

# appending the values found in DICTFILES to DICTH_modulename
DICTH_selectors := $(foreach i, $(patsubst %,selectors/%,$(DICTFILES)), $(wildcard $(i).h) $(wildcard $(i).hh))
# if dict header files exist, append to variable SRC
ifneq ($(DICTH_selectors),)
SRC += selectors/$(PKGNAME)_dict_selectors
endif

PROG += $(patsubst %,$(BINDIR)/%,$(PROGRAMS))

# appending the values found in FILES to the variable SRC
PROGSRC += $(patsubst %,selectors/%,$(PROGRAMS))

# VPATH += :selectors

# make sure the output directories are there
__dummy := $(shell mkdir -p $(DEPDIR)/selectors $(OBJDIR)/selectors)

# a couple of rules to copy executable files correctly
$(BINDIR)/%: selectors/%
	cp $^ $@

#$(BINDIR)/%: selectors/%.bin
#	cp $^ $@

$(BINDIR)/%: selectors/bin/%
	cp $^ $@

$(BINDIR)/%: ${OBJDIR}/selectors/%.bin
	cp $^ $@


# module.mk appends to FILES and DICTFILES
FILES :=
DICTFILES :=

# including the module.mk file
-include $(patsubst %, %/module.mk,src)

# appending the values found in FILES to the variable SRC
SRC += $(patsubst %,src/%,$(FILES))

# appending the values found in DICTFILES to DICTH_modulename
DICTH_src := $(foreach i, $(patsubst %,src/%,$(DICTFILES)), $(wildcard $(i).h) $(wildcard $(i).hh))
# if dict header files exist, append to variable SRC
ifneq ($(DICTH_src),)
SRC += src/$(PKGNAME)_dict_src
endif

PROG += $(patsubst %,$(BINDIR)/%,$(PROGRAMS))

# appending the values found in FILES to the variable SRC
PROGSRC += $(patsubst %,src/%,$(PROGRAMS))

# VPATH += :src

# make sure the output directories are there
__dummy := $(shell mkdir -p $(DEPDIR)/src $(OBJDIR)/src)

# a couple of rules to copy executable files correctly
$(BINDIR)/%: src/%
	cp $^ $@

#$(BINDIR)/%: src/%.bin
#	cp $^ $@

$(BINDIR)/%: src/bin/%
	cp $^ $@

$(BINDIR)/%: ${OBJDIR}/src/%.bin
	cp $^ $@


# module.mk appends to FILES and DICTFILES
FILES :=
DICTFILES :=

# including the module.mk file
-include $(patsubst %, %/module.mk,tools)

# appending the values found in FILES to the variable SRC
SRC += $(patsubst %,tools/%,$(FILES))

# appending the values found in DICTFILES to DICTH_modulename
DICTH_tools := $(foreach i, $(patsubst %,tools/%,$(DICTFILES)), $(wildcard $(i).h) $(wildcard $(i).hh))
# if dict header files exist, append to variable SRC
ifneq ($(DICTH_tools),)
SRC += tools/$(PKGNAME)_dict_tools
endif

PROG += $(patsubst %,$(BINDIR)/%,$(PROGRAMS))

# appending the values found in FILES to the variable SRC
PROGSRC += $(patsubst %,tools/%,$(PROGRAMS))

# VPATH += :tools

# make sure the output directories are there
__dummy := $(shell mkdir -p $(DEPDIR)/tools $(OBJDIR)/tools)

# a couple of rules to copy executable files correctly
$(BINDIR)/%: tools/%
	cp $^ $@

#$(BINDIR)/%: tools/%.bin
#	cp $^ $@

$(BINDIR)/%: tools/bin/%
	cp $^ $@

$(BINDIR)/%: ${OBJDIR}/tools/%.bin
	cp $^ $@


# module.mk appends to FILES and DICTFILES
FILES :=
DICTFILES :=

# including the module.mk file
-include $(patsubst %, %/module.mk,utils)

# appending the values found in FILES to the variable SRC
SRC += $(patsubst %,utils/%,$(FILES))

# appending the values found in DICTFILES to DICTH_modulename
DICTH_utils := $(foreach i, $(patsubst %,utils/%,$(DICTFILES)), $(wildcard $(i).h) $(wildcard $(i).hh))
# if dict header files exist, append to variable SRC
ifneq ($(DICTH_utils),)
SRC += utils/$(PKGNAME)_dict_utils
endif

PROG += $(patsubst %,$(BINDIR)/%,$(PROGRAMS))

# appending the values found in FILES to the variable SRC
PROGSRC += $(patsubst %,utils/%,$(PROGRAMS))

# VPATH += :utils

# make sure the output directories are there
__dummy := $(shell mkdir -p $(DEPDIR)/utils $(OBJDIR)/utils)

# a couple of rules to copy executable files correctly
$(BINDIR)/%: utils/%
	cp $^ $@

#$(BINDIR)/%: utils/%.bin
#	cp $^ $@

$(BINDIR)/%: utils/bin/%
	cp $^ $@

$(BINDIR)/%: ${OBJDIR}/utils/%.bin
	cp $^ $@

