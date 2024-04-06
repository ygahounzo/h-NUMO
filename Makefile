RM = /bin/rm -rf
MV = /bin/mv -f
CP = /bin/cp -f

export PLATFORMS:= macbook-p4est-intel borah borah-debug macbook-p4est-brew falcon

export NUMO_DIR    := $(CURDIR)
export DEPEND_FILE := $(NUMO_DIR)/depend.mk
export CONFIG_USER := $(NUMO_DIR)/config.user
export P4EST_LOC   := $(NUMO_DIR)/p4est

export TARGET=$(NUMO_DIR)/bin/numo3d

include $(NUMO_DIR)/config.p4est

SUBDIR := $(NUMO_DIR)/bin
SUBDIR += $(NUMO_DIR)/include
SUBDIR += $(NUMO_DIR)/libs

.NOTPARALLEL:

.PHONY: p4est

$(PLATFORMS): %: %-config p4est/local/lib/libp4est.a depend $(SUBDIR)
	$(MAKE) -C src $@ TARGET=$(TARGET)

$(SUBDIR):
	if [ ! -d $@ ]; then mkdir -p $@;fi

p4est:
	$(NUMO_DIR)/p4est/get_p4est_install_dir.sh $(P4EST_LOC)

p4est/local/lib/libp4est.a:
	cd $(NUMO_DIR)/p4est && \
	tar -xjf p4est_feature.tar.bj && \
	./configure $(P4EST_CONF) && \
	make -j && \
	make check && \
	make install

p4est_clean:
	cd $(NUMO_DIR)/p4est && \
	mv p4est_feature.tar.bj ../unused/. && \
	mv get_* ../unused/. && \
	rm -rf * && \
	mv ../unused/p4est_feature.tar.bj . && \
	mv ../unused/get_p4est* .

depend:
	$(NUMO_DIR)/make_depend.pl -s -l -o $(DEPEND_FILE) $(NUMO_DIR)/src

emptyrule:

clean : depend
	$(MAKE) -C src $@
	for dir in $(SUBDIRS); do (echo $(RM) $$dir ; $(RM) $$dir) ; done

deepclean : clean depend p4est_clean occa_clean
	$(MAKE) -C src $@
	(if [ -e $(DEPEND_FILE) ] ; then $(RM) $(DEPEND_FILE) ; fi )
