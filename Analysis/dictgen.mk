core/$(PKGNAME)_dict_core.cc: $(DICTH_core)
	@echo "Generating dictionary $@..."
	@$(ROOTCINT) -f $@ -c $(filter -I% -D%,$(CXXFLAGS)) $^
exe/$(PKGNAME)_dict_exe.cc: $(DICTH_exe)
	@echo "Generating dictionary $@..."
	@$(ROOTCINT) -f $@ -c $(filter -I% -D%,$(CXXFLAGS)) $^
fit/$(PKGNAME)_dict_fit.cc: $(DICTH_fit)
	@echo "Generating dictionary $@..."
	@$(ROOTCINT) -f $@ -c $(filter -I% -D%,$(CXXFLAGS)) $^
ghmzz/$(PKGNAME)_dict_ghmzz.cc: $(DICTH_ghmzz)
	@echo "Generating dictionary $@..."
	@$(ROOTCINT) -f $@ -c $(filter -I% -D%,$(CXXFLAGS)) $^
selectors/$(PKGNAME)_dict_selectors.cc: $(DICTH_selectors)
	@echo "Generating dictionary $@..."
	@$(ROOTCINT) -f $@ -c $(filter -I% -D%,$(CXXFLAGS)) $^
src/$(PKGNAME)_dict_src.cc: $(DICTH_src)
	@echo "Generating dictionary $@..."
	@$(ROOTCINT) -f $@ -c $(filter -I% -D%,$(CXXFLAGS)) $^
tools/$(PKGNAME)_dict_tools.cc: $(DICTH_tools)
	@echo "Generating dictionary $@..."
	@$(ROOTCINT) -f $@ -c $(filter -I% -D%,$(CXXFLAGS)) $^
utils/$(PKGNAME)_dict_utils.cc: $(DICTH_utils)
	@echo "Generating dictionary $@..."
	@$(ROOTCINT) -f $@ -c $(filter -I% -D%,$(CXXFLAGS)) $^
