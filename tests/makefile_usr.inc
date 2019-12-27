SRC_INC := ./
SRC_USR := ./
$(OBJDIR)/solid.o     :$(SRC_USR)/solid.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@

