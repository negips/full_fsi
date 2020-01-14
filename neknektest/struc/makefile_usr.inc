SRC_INC := ../src
SRC_USR := ../src
$(OBJDIR)/solid.o     :$(SRC_USR)/solid.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/plan_struc3.o     :$(SRC_USR)/plan_struc3.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/debug_struc.o     :$(SRC_USR)/debug_struc.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/fsi_plan3.o       :$(SRC_USR)/fsi_plan3.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/fsi_neknek.o       :$(SRC_USR)/fsi_neknek.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@

